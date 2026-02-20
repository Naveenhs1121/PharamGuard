"""
PharmaGuard Risk Engine
=======================
CPIC-aligned pharmacogenomic risk prediction engine.

Pipeline:
  VCF variants → rsID annotation → diplotype inference
  → phenotype classification → drug-phenotype rule lookup
  → confidence scoring → structured clinical output

Supported Genes : CYP2D6, CYP2C19, CYP2C9, SLCO1B1, TPMT, DPYD
Supported Drugs : Codeine, Warfarin, Clopidogrel, Simvastatin,
                  Azathioprine, Fluorouracil
"""

import logging
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Any

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("PharmaGuard.RiskEngine")

# ---------------------------------------------------------------------------
# Constants — Phenotype codes
# ---------------------------------------------------------------------------
PM  = "Poor Metabolizer"
IM  = "Intermediate Metabolizer"
NM  = "Normal Metabolizer"
RM  = "Rapid Metabolizer"
URM = "Ultrarapid Metabolizer"
INDETERMINATE = "Indeterminate"

# Metabolizer activity scores (used for diplotype → phenotype mapping)
ACTIVITY_SCORE = {
    "LOF": 0.0,   # Loss-of-function allele
    "DEF": 0.5,   # Decreased-function allele
    "NF_CYP2C19": 0.5, # Normal-function on CYP2C19 specific scale
    "NF":  1.0,   # Normal-function allele (standard)
    "INF_CYP2C19": 1.0, # Increased-function on CYP2C19 specific scale
    "INF": 2.0,   # Increased-function allele (gene duplication / gain)
}

# Default Activity Scores (if no variants found, i.e., *1/*1)
DEFAULT_ACTIVITY_SCORES = {
    "CYP2D6": 2.0,
    "CYP2C19": 1.0,  # Per user prompt: *1/*1 = 1.0
    "CYP2C9": 2.0,
    "SLCO1B1": 2.0,  # *1/*1 (NF/NF)
    "TPMT": 2.0,
    "DPYD": 2.0,
}

# ---------------------------------------------------------------------------
# rsID → Allele Annotation Database
# Maps known pharmacogenomic rsIDs to:
#   gene, star_allele, function, activity_score, evidence_strength
# Evidence strength: "high" | "moderate" | "low"
# ---------------------------------------------------------------------------
RSID_DATABASE: Dict[str, Dict] = {
    # ── CYP2D6 ──────────────────────────────────────────────────────────────
    "rs3892097":  {"gene": "CYP2D6",  "star": "*4",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs35742686": {"gene": "CYP2D6",  "star": "*3",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs5030655":  {"gene": "CYP2D6",  "star": "*6",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs28371725": {"gene": "CYP2D6",  "star": "*41", "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs1065852":  {"gene": "CYP2D6",  "star": "*10", "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs16947":    {"gene": "CYP2D6",  "star": "*2",  "function": "NF",  "activity": 1.0, "evidence": "high"},
    "rs1135840":  {"gene": "CYP2D6",  "star": "*2",  "function": "NF",  "activity": 1.0, "evidence": "moderate"},
    "rs5030865":  {"gene": "CYP2D6",  "star": "*8",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    # Add *5 (Gene deletion) placeholder if RSID known, but for now rely on existing.

    # ── CYP2C19 ─────────────────────────────────────────────────────────────
    # Note: Using *1=0.5 scale per prompt.
    "rs4244285":  {"gene": "CYP2C19", "star": "*2",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs4986893":  {"gene": "CYP2C19", "star": "*3",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs28399504": {"gene": "CYP2C19", "star": "*4",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs56337013": {"gene": "CYP2C19", "star": "*5",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs12248560": {"gene": "CYP2C19", "star": "*17", "function": "INF", "activity": 1.0, "evidence": "high"},
    "rs41291556": {"gene": "CYP2C19", "star": "*17", "function": "INF", "activity": 1.0, "evidence": "moderate"},

    # ── CYP2C9 ──────────────────────────────────────────────────────────────
    "rs1799853":  {"gene": "CYP2C9",  "star": "*2",  "function": "DEF", "activity": 0.5, "evidence": "high"},
    "rs1057910":  {"gene": "CYP2C9",  "star": "*3",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs28371686": {"gene": "CYP2C9",  "star": "*5",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs9332131":  {"gene": "CYP2C9",  "star": "*6",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs7900194":  {"gene": "CYP2C9",  "star": "*8",  "function": "DEF", "activity": 0.5, "evidence": "moderate"},

    # ── VKORC1 ──────────────────────────────────────────────────────────────
    "rs9923231":  {"gene": "VKORC1",  "star": "-1639G>A", "function": "BS", "activity": 0.0, "evidence": "high"},

    # ── SLCO1B1 ─────────────────────────────────────────────────────────────
    # Mapping *5 to 0.0 effectively makes it "Poor Function" contribution relative to *1 (Normal=1.0)
    "rs4149056":  {"gene": "SLCO1B1", "star": "*5",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs2306283":  {"gene": "SLCO1B1", "star": "*1B", "function": "NF",  "activity": 1.0, "evidence": "moderate"},
    "rs11045819": {"gene": "SLCO1B1", "star": "*15", "function": "LOF", "activity": 0.0, "evidence": "moderate"},
    "rs4363657":  {"gene": "SLCO1B1", "star": "*5",  "function": "LOF", "activity": 0.0, "evidence": "high"},

    # ── TPMT ────────────────────────────────────────────────────────────────
    "rs1142345":  {"gene": "TPMT",    "star": "*3C", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs1800460":  {"gene": "TPMT",    "star": "*3B", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs1800462":  {"gene": "TPMT",    "star": "*2",  "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs1800584":  {"gene": "TPMT",    "star": "*4",  "function": "LOF", "activity": 0.0, "evidence": "moderate"},

    # ── DPYD ────────────────────────────────────────────────────────────────
    "rs3918290":  {"gene": "DPYD",    "star": "*2A", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs55886062": {"gene": "DPYD",    "star": "*13", "function": "LOF", "activity": 0.0, "evidence": "high"},
    "rs67376798": {"gene": "DPYD",    "star": "HapB3","function": "DEF","activity": 0.5, "evidence": "high"},
    "rs75017182": {"gene": "DPYD",    "star": "HapB3","function": "DEF","activity": 0.5, "evidence": "moderate"},
}

# ---------------------------------------------------------------------------
# Gene → Phenotype rules (activity-score based, CPIC aligned)
# Each rule: (lower_bound_inclusive, upper_bound_exclusive, phenotype)
# ---------------------------------------------------------------------------
GENE_PHENOTYPE_RULES: Dict[str, List] = {
    "CYP2D6": [
        (0.0, 0.01, PM),
        (0.01, 1.25, IM),
        (1.25, 2.25, NM),
        (2.25, 99.0, URM),
    ],
    "CYP2C19": [
        (0.0, 0.01, PM),
        (0.01, 0.9, IM),   # Target 0.5
        (0.9, 1.25, NM),   # Target 1.0
        (1.25, 1.75, RM),  # Target 1.5
        (1.75, 99.0, URM), # Target 2.0+
    ],
    "CYP2C9": [
        (0.0, 0.01, PM),
        (0.01, 1.5, IM),   # Target 0.5 - 1.0 (Strictly < 2.0)
        (1.5, 99.0, NM),   # Target 2.0
    ],
    "SLCO1B1": [
        (0.0, 0.5, "Poor Function"),     # 0 alleles (*5/*5 = 0)
        (0.5, 1.5, "Decreased Function"), # 1 allele (*1/*5 = 1)
        (1.5, 99.0, "Normal Function"),   # 2 alleles (*1/*1 = 2)
    ],
    "TPMT": [
        (0.0, 0.5, PM),
        (0.5, 1.5, IM),
        (1.5, 99.0, NM),
    ],
    "DPYD": [
        (0.0, 0.5, PM),
        (0.5, 1.75, IM),
        (1.75, 99.0, NM),
    ],
}

# ---------------------------------------------------------------------------
# CPIC Drug–Phenotype Rules
# ---------------------------------------------------------------------------
DRUG_RULES: Dict[str, Dict[str, Dict[str, Dict]]] = {
    # ── CYP2D6 Drugs ────────────────────────────────────────────────────────
    "codeine": {
        "CYP2D6": {
            PM: {"label": "Ineffective", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Codeine (2021)", "clinical_action": "Avoid codeine. Risk of no analgesia. Use non-opioid or non-CYP2D6 opioid."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Codeine (2021)", "clinical_action": "Use lowest effective dose with caution. Monitor for reduced efficacy."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Codeine (2021)", "clinical_action": "Standard dosing."},
            RM: {"label": "Safe", "severity": "low", "confidence_base": 0.85, "cpic_guideline": "CPIC Codeine (2021)", "clinical_action": "Standard dosing. Monitor side effects."},
            URM: {"label": "Toxic", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Codeine (2021)", "clinical_action": "Avoid. Risk of life-threatening toxicity (excess morphine)."},
        }
    },
    "tramadol": {
        "CYP2D6": {
            PM: {"label": "Ineffective", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Tramadol (2021)", "clinical_action": "Avoid. Use alternative analgesic."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.80, "cpic_guideline": "CPIC Tramadol (2021)", "clinical_action": "Use with caution; reduced efficacy possible."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Tramadol (2021)", "clinical_action": "Standard dosing."},
            URM: {"label": "Toxic", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Tramadol (2021)", "clinical_action": "Avoid. Risk of serious adverse events."},
        }
    },
    "amitriptyline": {
        "CYP2D6": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Avoid TCA or reduce dose by 50%. Consider alternative."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.80, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Consider 25% dose reduction. Monitor drug levels."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Standard dosing."},
            URM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Avoid TCA or increase dose with monitoring."},
        },
        "CYP2C19": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Avoid or reduce dose by 50%. Monitor levels."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.80, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Reduce starting dose by 25%. Monitor."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Standard dosing."},
            RM: {"label": "Ineffective", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Consider dose titration. Monitor for reduced efficacy."},
            URM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Consider dose titration. Monitor for reduced efficacy."},
        }
    },
    "nortriptyline": {
        "CYP2D6": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Avoid TCA or reduce dose by 50%. Consider alternative."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.80, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Consider 25% dose reduction. Monitor drug levels."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Standard dosing."},
            URM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC TCA Guideline", "clinical_action": "Avoid TCA or increase dose with monitoring."},
        }
    },
    "ondansetron": {
        "CYP2D6": {
            PM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Antiemetics", "clinical_action": "Use alternative antiemetic (e.g., granisetron). Standard dose may be ineffective."},
            IM: {"label": "Adjust Dosage", "severity": "low", "confidence_base": 0.75, "cpic_guideline": "CPIC Antiemetics", "clinical_action": "Use with caution; may have reduced efficacy."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Antiemetics", "clinical_action": "Standard dosing."},
            URM: {"label": "Safe", "severity": "none", "confidence_base": 0.80, "cpic_guideline": "CPIC Antiemetics", "clinical_action": "Standard dosing."},
        }
    },
    "tamoxifen": {
        "CYP2D6": {
            PM: {"label": "Ineffective", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Tamoxifen", "clinical_action": "Avoid if possible. Inadequate active metabolite. Recommend alternative (e.g. aromatase inhibitor)."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Tamoxifen", "clinical_action": "Consider higher dose (40mg/day) or alternative. Monitor endoxifen."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Tamoxifen", "clinical_action": "Standard 20 mg/day dosing."},
            URM: {"label": "Safe", "severity": "none", "confidence_base": 0.85, "cpic_guideline": "CPIC Tamoxifen", "clinical_action": "Standard 20 mg/day dosing."},
        }
    },
    "atomoxetine": {
        "CYP2D6": {
            PM: {"label": "Toxic", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC Atomoxetine", "clinical_action": "Initiate at 50% of normal dose; titrate slowly. Increased exposure risk."},
            IM: {"label": "Adjust Dosage", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC Atomoxetine", "clinical_action": "Standard starting dose with close monitoring."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Atomoxetine", "clinical_action": "Standard dosing."},
            URM: {"label": "Ineffective", "severity": "low", "confidence_base": 0.75, "cpic_guideline": "CPIC Atomoxetine", "clinical_action": "Standard dosing (limited data)."},
        }
    },
    "haloperidol": {
        "CYP2D6": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.85, "cpic_guideline": "CPIC Antipsychotics (Level B)", "clinical_action": "Use lowest effective dose; high risk of ADRs. Consider alternatives."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.75, "cpic_guideline": "CPIC Antipsychotics", "clinical_action": "Monitor carefully; consider dose reduction."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Antipsychotics", "clinical_action": "Standard dosing."},
            URM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.80, "cpic_guideline": "CPIC Antipsychotics", "clinical_action": "May need higher doses; monitor for reduced efficacy."},
        }
    },

    # ── CYP2C19 Drugs ───────────────────────────────────────────────────────
    "clopidogrel": {
        "CYP2C19": {
            PM: {"label": "Ineffective", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Clopidogrel", "clinical_action": "Avoid. Minimal antiplatelet effect. Use prasugrel or ticagrelor."},
            IM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Clopidogrel", "clinical_action": "Use with caution; consider alternative. If used, consider higher dose."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Clopidogrel", "clinical_action": "Standard 75 mg/day."},
            RM: {"label": "Safe", "severity": "none", "confidence_base": 0.85, "cpic_guideline": "CPIC Clopidogrel", "clinical_action": "Standard dosing."},
            URM: {"label": "Safe", "severity": "none", "confidence_base": 0.85, "cpic_guideline": "CPIC Clopidogrel", "clinical_action": "Standard dosing."},
        }
    },
    "voriconazole": {
        "CYP2C19": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC Voriconazole", "clinical_action": "High exposure; reduce dose and monitor levels. Risk of ADRs."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.80, "cpic_guideline": "CPIC Voriconazole", "clinical_action": "Monitor trough levels; may need dose reduction."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Voriconazole", "clinical_action": "Standard dosing with TDM."},
            RM: {"label": "Ineffective", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC Voriconazole", "clinical_action": "Significantly reduced exposure; ineffective. Use alternative antifungal."},
            URM: {"label": "Ineffective", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Voriconazole", "clinical_action": "Significantly reduced exposure; ineffective. Use alternative antifungal."},
        }
    },
    "citalopram": {
        "CYP2C19": {
            PM: {"label": "Toxic", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Reduce dose by 50% (max 20mg/day). Increased QT risk."},
            IM: {"label": "Adjust Dosage", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Use lowest effective dose; monitor QT."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Standard dosing."},
            RM: {"label": "Ineffective", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Standard dosing (consider alternative if no response)."},
            URM: {"label": "Ineffective", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Standard dosing (consider alternative if no response)."},
        }
    },
    "escitalopram": {
        "CYP2C19": {
            PM: {"label": "Toxic", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Reduce dose by 50% (max 10mg/day). Increased QT risk."},
            IM: {"label": "Adjust Dosage", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Use lowest effective dose; monitor QT."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Standard dosing."},
            RM: {"label": "Ineffective", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Standard dosing (consider alternative if no response)."},
            URM: {"label": "Ineffective", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC SSRIs", "clinical_action": "Standard dosing (consider alternative if no response)."},
        }
    },
    "omeprazole": {
        "CYP2C19": {
            PM: {"label": "Toxic", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC PPIs", "clinical_action": "Initiate at 50% of standard dose. Monitor for ADRs (increased exposure)."},
            IM: {"label": "Adjust Dosage", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC PPIs", "clinical_action": "Consider dose reduction; start at lower end."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC PPIs", "clinical_action": "Standard dosing."},
            RM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC PPIs", "clinical_action": "Consider doubling dose for H. pylori or erosive esophagitis."},
            URM: {"label": "Ineffective", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC PPIs", "clinical_action": "Consider doubling dose for H. pylori or erosive esophagitis."},
        }
    },

    # ── CYP2C9 Drugs ────────────────────────────────────────────────────────
    "warfarin": {
        "CYP2C9": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Warfarin", "clinical_action": "Start at significantly reduced dose (~5-6 mg/week). Very slow titration. High bleeding risk. Check VKORC1."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.88, "cpic_guideline": "CPIC Warfarin", "clinical_action": "Reduce starting dose by 25-50%. Close INR monitoring. Check VKORC1 status."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Warfarin", "clinical_action": "Standard dosing. Consider VKORC1 genotype for fine-tuning."},
        }
    },
    "celecoxib": {
        "CYP2C9": {
            PM: {"label": "Toxic", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC NSAIDs", "clinical_action": "Reduce starting dose by 25-50%. Use lowest effective dose."},
            IM: {"label": "Adjust Dosage", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC NSAIDs", "clinical_action": "Consider 25% dose reduction. Monitor."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC NSAIDs", "clinical_action": "Standard dosing."},
        }
    },
    "phenytoin": {
        "CYP2C9": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Phenytoin", "clinical_action": "Reduce dose by 25-50%. Monitor levels. Risk of toxicity."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Phenytoin", "clinical_action": "Consider 25% reduction. Monitor."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Phenytoin", "clinical_action": "Standard dosing."},
        }
    },
    "siponimod": {
        "CYP2C9": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Siponimod", "clinical_action": "Contraindicated; avoid use."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC Siponimod", "clinical_action": "Use 1 mg/day maintenance (vs 2 mg)."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.95, "cpic_guideline": "CPIC Siponimod", "clinical_action": "Standard 2 mg/day."},
        }
    },

    # ── SLCO1B1 Drugs ───────────────────────────────────────────────────────
    "simvastatin": {
        "SLCO1B1": {
            "Poor Function": {"label": "Toxic", "severity": "high", "confidence_base": 0.95, "cpic_guideline": "CPIC Statins", "clinical_action": "Avoid. High myopathy risk. Use alternative (rosuvastatin, pravastatin)."},
            "Decreased Function": {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Statins", "clinical_action": "Use max 20 mg/day. If higher needed, use alternative."},
            "Normal Function": {"label": "Safe", "severity": "none", "confidence_base": 0.95, "cpic_guideline": "CPIC Statins", "clinical_action": "Standard dosing (max 40 mg/day)."},
        }
    },
    "atorvastatin": {
        "SLCO1B1": {
            "Poor Function": {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Statins", "clinical_action": "Avoid or use lowest possible dose. Use alternative."},
            "Decreased Function": {"label": "Safe", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC Statins", "clinical_action": "Use max 40 mg/day. Monitor for myopathy."},
            "Normal Function": {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Statins", "clinical_action": "Standard dosing."},
        }
    },
    "rosuvastatin": {
        "SLCO1B1": {
            "Poor Function": {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.85, "cpic_guideline": "CPIC Statins", "clinical_action": "Avoid or use max 20 mg/day. Use alternative."},
            "Decreased Function": {"label": "Safe", "severity": "low", "confidence_base": 0.80, "cpic_guideline": "CPIC Statins", "clinical_action": "Use max 20 mg/day."},
            "Normal Function": {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Statins", "clinical_action": "Standard dosing."},
        }
    },

    # ── TPMT Drugs ──────────────────────────────────────────────────────────
    "azathioprine": {
        "TPMT": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.98, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Reduce dose 10-fold or avoid. Risk of life-threatening myelosuppression."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Start at 30-70% full dose. Titrate based on toxicity/efficacy. Monitor CBC."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.95, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Standard dosing (2-3 mg/kg/day)."},
        }
    },
    "mercaptopurine": {
        "TPMT": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.98, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Reduce dose to 10% of normal. Titrate cautiously."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Start at 30-80% standard dose. Monitor CBC."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.95, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Standard dosing."},
        }
    },
    "thioguanine": {
        "TPMT": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.98, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Reduce dose to 10% of standard. Risk of fatal myelosuppression."},
            IM: {"label": "Adjust Dosage", "severity": "moderate", "confidence_base": 0.90, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Start at 30-50% of normal dose."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.95, "cpic_guideline": "CPIC Thiopurines", "clinical_action": "Standard dosing."},
        }
    },

    # ── DPYD Drugs ──────────────────────────────────────────────────────────
    "fluorouracil": {
        "DPYD": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.98, "cpic_guideline": "CPIC Fluoropyrimidines", "clinical_action": "Avoid. Fatal toxicity risk. If no alt, use extreme caution."},
            IM: {"label": "Adjust Dosage", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC Fluoropyrimidines", "clinical_action": "Reduce starting dose by 50%. Monitor closely."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Fluoropyrimidines", "clinical_action": "Standard dosing."},
        }
    },
    "capecitabine": {
        "DPYD": {
            PM: {"label": "Toxic", "severity": "high", "confidence_base": 0.98, "cpic_guideline": "CPIC Fluoropyrimidines", "clinical_action": "Avoid. Fatal toxicity risk."},
            IM: {"label": "Adjust Dosage", "severity": "high", "confidence_base": 0.90, "cpic_guideline": "CPIC Fluoropyrimidines", "clinical_action": "Reduce starting dose by 50%. Monitor closely."},
            NM: {"label": "Safe", "severity": "none", "confidence_base": 0.90, "cpic_guideline": "CPIC Fluoropyrimidines", "clinical_action": "Standard dosing."},
        }
    },
}

# ---------------------------------------------------------------------------
# Data classes for structured output
# ---------------------------------------------------------------------------
@dataclass
class VariantAnnotation:
    rsid: str
    gene: str
    star_allele: str
    function: str
    activity_score: float
    evidence_strength: str
    chrom: str = ""
    pos: int = 0
    ref: str = ""
    alt: List[str] = field(default_factory=list)


@dataclass
class GeneResult:
    gene: str
    detected_rsids: List[str]
    annotated_variants: List[VariantAnnotation]
    total_activity_score: float
    diplotype: str
    phenotype: str
    phenotype_reasoning: str


@dataclass
class DrugRiskResult:
    drug: str
    gene_used: str
    phenotype: str
    risk_label: str
    severity: str
    confidence_score: float
    clinical_action: str
    cpic_guideline: str
    supporting_variants: List[str]
    reasoning: str
    evidence_strength: str


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _normalise_drug(drug: str) -> str:
    """Lowercase and strip the drug name for rule lookup."""
    return drug.strip().lower().replace("-", "").replace(" ", "")


def _annotate_variant(variant: Dict) -> Optional[VariantAnnotation]:
    """
    Match a VCF variant dict to the rsID database.
    Returns VariantAnnotation or None if not in database.
    """
    rsid_raw = variant.get("rsid")
    if not rsid_raw:
        return None

    rsid_list = rsid_raw if isinstance(rsid_raw, list) else [rsid_raw]

    for rsid in rsid_list:
        if rsid and rsid in RSID_DATABASE:
            db = RSID_DATABASE[rsid]
            return VariantAnnotation(
                rsid=rsid,
                gene=db["gene"],
                star_allele=db["star"],
                function=db["function"],
                activity_score=db["activity"],
                evidence_strength=db["evidence"],
                chrom=str(variant.get("chrom", "")),
                pos=int(variant.get("pos", 0)),
                ref=str(variant.get("ref", "")),
                alt=[str(a) for a in variant.get("alt", [])],
            )
    return None


def _classify_phenotype(gene: str, activity_score: float) -> str:
    """Map a gene's total activity score to its metabolizer phenotype."""
    rules = GENE_PHENOTYPE_RULES.get(gene)
    if rules is None:
        return INDETERMINATE

    for lo, hi, phenotype in rules:
        if lo <= activity_score < hi:
            return phenotype

    # Edge case: score exactly at the max end
    return rules[-1][2]


def _build_diplotype_label(annotated: List[VariantAnnotation], gene: str) -> str:
    """
    Build a human-readable diplotype string from detected star alleles.
    Assumes biallelic diploid (two allele slots).
    """
    stars = sorted(set(a.star_allele for a in annotated))
    if not stars:
        return f"{gene}:*1/*1 (reference assumed)"
    if len(stars) == 1:
        return f"{gene}:{stars[0]}/{stars[0]}"
    return f"{gene}:{stars[0]}/{stars[1]}"


def _evidence_strength_factor(annotated: List[VariantAnnotation]) -> float:
    """
    Returns a multiplier [0.6, 1.0] based on the weakest evidence level
    in the supporting variants, penalising low-evidence calls.
    """
    if not annotated:
        return 0.6
    weights = {"high": 1.0, "moderate": 0.85, "low": 0.65}
    return min(weights.get(a.evidence_strength, 0.65) for a in annotated)


def _phenotype_reasoning(gene: str, phenotype: str, score: float,
                          annotated: List[VariantAnnotation]) -> str:
    """Generate a plain-English phenotype reasoning string."""
    variant_list = ", ".join(
        f"{a.rsid} ({a.star_allele}, {a.function})" for a in annotated
    ) or "no database-annotated variants"
    return (
        f"{gene} activity score = {score:.2f}. "
        f"Detected variants: {variant_list}. "
        f"Phenotype classified as {phenotype} per CPIC activity-score model."
    )


def _risk_reasoning(drug: str, gene: str, phenotype: str,
                    rule: Dict, score: float) -> str:
    """Generate a plain-English risk reasoning string."""
    return (
        f"{drug.capitalize()} risk assessment based on {gene} phenotype "
        f"({phenotype}, activity score {score:.2f}). "
        f"CPIC classification: '{rule['label']}' "
        f"(severity: {rule['severity']})."
    )


def _compute_confidence(
    base: float,
    activity_score: float,
    n_variants: int,
    evidence_factor: float,
    has_variants: bool,
) -> float:
    """
    Compute final confidence score on [0.0, 1.0].

    Factors:
    - base: rule-level CPIC confidence
    - variant_support: max 0.05 bonus for ≥ 2 supporting variants
    - evidence_factor: penalty for low-evidence rsIDs
    - no_variant_penalty: small deduction if phenotype inferred from absence
    """
    if not has_variants:
        # Inferred as NM from absence of known LOF variants — less certain
        return round(min(base * evidence_factor * 0.85, 1.0), 3)

    variant_bonus = min(n_variants * 0.02, 0.05)
    confidence = (base + variant_bonus) * evidence_factor
    return round(min(confidence, 1.0), 3)


# ---------------------------------------------------------------------------
# Core prediction functions
# ---------------------------------------------------------------------------

def analyse_gene(gene: str, variants: List[Dict]) -> GeneResult:
    """
    Analyse all variants for a single gene.
    Handles zygosity (Het/Hom) and infers missing alleles as Reference (*1).
    """
    logger.info(f"Analysing gene: {gene} with {len(variants)} candidate variants")

    gene_variants = [v for v in variants if v.get("gene") == gene]
    detected_alleles: List[VariantAnnotation] = []

    for v in gene_variants:
        ann = _annotate_variant(v)
        if ann and ann.gene == gene:
            # Parse Genotype
            gt = v.get("gt", "0/1") # Default to Het if missing
            is_hom = "1/1" in gt or "1|1" in gt or "2/2" in gt # Basic check
            
            # Add allele instances
            detected_alleles.append(ann) # First allele
            if is_hom:
                detected_alleles.append(ann) # Second allele (same)
            
            logger.debug(f"  {gene}: matched {ann.rsid} ({gt}) -> {ann.star_allele} (x{2 if is_hom else 1})")
        else:
            rsid = v.get("rsid", "unknown")
            logger.debug(f"  {gene}: rsID {rsid!r} not in database — skipped")

    # Calculate Activity Score
    # 1. Get default diplotype score (*1/*1)
    default_diplotype_score = DEFAULT_ACTIVITY_SCORES.get(gene, 2.0)
    ref_allele_score = default_diplotype_score / 2.0
    
    # 2. Determine final alleles (max 2 for diploidy)
    # Heuristic: If >2 alleles found, take top 2 (usually severe ones). 
    # But for now, just take first 2 or pad with Ref.
    
    final_alleles = []
    
    # Copy detected
    for a in detected_alleles:
        final_alleles.append(a)
    
    # Pad with Reference if needed (target 2 alleles)
    # Note: We represent Reference as None or a dummy Annotation?
    # Let's count score directly.
    
    activity_score = 0.0
    
    # Sum up detected variants (up to 2)
    # Issue: If 3 variants? e.g. *2 and *4. 
    # If unphased, simpler to just sum them? No, *1/*2 + *4 is impossible (3 alleles).
    # Assume 2 slots.
    
    # Logic: Start with score = 0. Fill 2 slots.
    # Slot 1: detected[0] or Ref
    # Slot 2: detected[1] or Ref
    
    slots_filled = 0
    annotated_variants_list = [] # For reporting
    
    # Sort detected alleles? Maybe LOF first?
    # If we have *2 (NF) and *4 (PM) and *1 (Ref).
    # If VCF has *2 (Het) and *4 (Het).
    # Diplotype *2/*4.
    # If VCF has *2 (Hom). *2/*2.
    
    # We use the 'detected_alleles' list which has duplicated entries for Hom.
    # We take the first 2.
    
    for allele in detected_alleles[:2]:
        activity_score += allele.activity_score
        slots_filled += 1
        annotated_variants_list.append(allele)
        
    while slots_filled < 2:
        activity_score += ref_allele_score
        slots_filled += 1
        # Implicit *1
        
    phenotype = _classify_phenotype(gene, activity_score)
    diplotype = _build_diplotype_label(annotated_variants_list, gene)
    
    # Special handling for Diplotype string to show *1 if implicit
    if len(annotated_variants_list) == 0:
        diplotype = f"{gene}:*1/*1" # Fully Reference
    elif len(annotated_variants_list) == 1:
        diplotype = f"{gene}:{annotated_variants_list[0].star_allele}/*1"
    
    reasoning = _phenotype_reasoning(gene, phenotype, activity_score, annotated_variants_list)

    logger.info(f"  {gene}: activity={activity_score:.2f}, phenotype={phenotype}, diplotype={diplotype}")

    return GeneResult(
        gene=gene,
        detected_rsids=[a.rsid for a in annotated_variants_list],
        annotated_variants=annotated_variants_list,
        total_activity_score=activity_score,
        diplotype=diplotype,
        phenotype=phenotype,
        phenotype_reasoning=reasoning,
    )


def predict_drug_risk(drug_name: str, gene_results: Dict[str, GeneResult]) -> DrugRiskResult:
    """
    Predict risk for a single drug given pre-computed gene results.

    Returns a DrugRiskResult with risk label, severity, confidence, and
    clinical action.
    """
    drug_key = _normalise_drug(drug_name)
    logger.info(f"Predicting risk for drug: {drug_name!r} (key={drug_key!r})")

    drug_rules = DRUG_RULES.get(drug_key)
    if not drug_rules:
        logger.warning(f"Drug '{drug_name}' not in rule database.")
        return DrugRiskResult(
            drug=drug_name,
            gene_used="N/A",
            phenotype=INDETERMINATE,
            risk_label="Unknown",
            severity="unknown",
            confidence_score=0.0,
            clinical_action=(
                f"No pharmacogenomic data available for '{drug_name}'. "
                "Proceed per standard clinical guidelines."
            ),
            cpic_guideline="N/A",
            supporting_variants=[],
            reasoning=f"Drug '{drug_name}' is not in the current CPIC rule set.",
            evidence_strength="none",
        )

    # Try each gene defined for this drug (priority order matters)
    for gene, phenotype_rules in drug_rules.items():
        gene_result = gene_results.get(gene)
        if gene_result is None:
            logger.debug(f"  No gene result available for {gene} — skipping rule")
            continue

        phenotype = gene_result.phenotype
        if phenotype == INDETERMINATE:
            logger.debug(f"  {gene} phenotype is Indeterminate — skipping")
            continue

        rule = phenotype_rules.get(phenotype)
        if rule is None:
            logger.debug(f"  No rule for {gene}/{phenotype} combo — defaulting to Safe")
            rule = {
                "label": "Safe",
                "severity": "none",
                "confidence_base": 0.70,
                "clinical_action": (
                    f"No specific {drug_name} recommendation for {phenotype} "
                    f"{gene} phenotype. Use standard dosing with caution."
                ),
                "cpic_guideline": "No specific CPIC guidance",
            }

        has_variants = len(gene_result.annotated_variants) > 0
        ev_factor = _evidence_strength_factor(gene_result.annotated_variants)
        confidence = _compute_confidence(
            base=rule["confidence_base"],
            activity_score=gene_result.total_activity_score,
            n_variants=len(gene_result.annotated_variants),
            evidence_factor=ev_factor,
            has_variants=has_variants,
        )

        ev_strength = (
            min((a.evidence_strength for a in gene_result.annotated_variants),
                key=lambda x: {"high": 0, "moderate": 1, "low": 2}.get(x, 3))
            if has_variants else "inferred"
        )

        reasoning = _risk_reasoning(
            drug_name, gene, phenotype, rule,
            gene_result.total_activity_score,
        )

        logger.info(
            f"  {drug_name}: gene={gene}, phenotype={phenotype}, "
            f"label={rule['label']}, confidence={confidence:.3f}"
        )

        return DrugRiskResult(
            drug=drug_name,
            gene_used=gene,
            phenotype=phenotype,
            risk_label=rule["label"],
            severity=rule["severity"],
            confidence_score=confidence,
            clinical_action=rule["clinical_action"],
            cpic_guideline=rule["cpic_guideline"],
            supporting_variants=gene_result.detected_rsids,
            reasoning=reasoning,
            evidence_strength=ev_strength,
        )

    # Reached here: no applicable gene result found
    logger.warning(f"No applicable gene/phenotype found for {drug_name}")
    return DrugRiskResult(
        drug=drug_name,
        gene_used="N/A",
        phenotype=INDETERMINATE,
        risk_label="Unknown",
        severity="unknown",
        confidence_score=0.40,
        clinical_action=(
            f"Insufficient genomic data to assess {drug_name} risk. "
            "No variants detected in relevant genes. Proceed with standard care."
        ),
        cpic_guideline="N/A",
        supporting_variants=[],
        reasoning=(
            f"Relevant genes for {drug_name} were not detected in the VCF. "
            "Risk set to Unknown."
        ),
        evidence_strength="none",
    )


# ---------------------------------------------------------------------------
# Public API — main entry points
# ---------------------------------------------------------------------------

TARGET_GENES = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]
SUPPORTED_DRUGS = [
    # CYP2D6
    "codeine", "tramadol", "amitriptyline", "nortriptyline", "ondansetron",
    "tamoxifen", "atomoxetine", "haloperidol",
    # CYP2C19
    "clopidogrel", "voriconazole", "citalopram", "escitalopram", "omeprazole",
    # CYP2C9
    "warfarin", "celecoxib", "phenytoin", "siponimod",
    # SLCO1B1
    "simvastatin", "atorvastatin", "rosuvastatin",
    # TPMT
    "azathioprine", "mercaptopurine", "thioguanine",
    # DPYD
    "fluorouracil", "capecitabine",
]


def predict_risk(variants: List[Dict], drug: str) -> Dict[str, Any]:
    """
    Legacy single-drug API (backwards-compatible with main.py).

    Args:
        variants : list of variant dicts from vcf_parser.extract_variants()
        drug     : drug name string

    Returns:
        dict with keys: label, severity, confidence (+ full result nested)
    """
    drug_key = _normalise_drug(drug)
    if drug_key not in SUPPORTED_DRUGS:
        logger.warning(f"Drug '{drug}' not supported — returning Unknown")
        return {
            "label": "Unknown",
            "severity": "unknown",
            "confidence": 0.0,
            "full_result": None,
        }

    gene_results = {
        gene: analyse_gene(gene, variants) for gene in TARGET_GENES
    }
    result = predict_drug_risk(drug, gene_results)
    return {
        "label": result.risk_label,
        "severity": result.severity,
        "confidence": result.confidence_score,
        "full_result": asdict(result),
    }


def predict_multi_drug(variants: List[Dict], drugs: List[str]) -> Dict[str, Any]:
    """
    Multi-drug prediction API.

    Args:
        variants : list of variant dicts from vcf_parser.extract_variants()
        drugs    : list of drug name strings

    Returns:
        dict with:
          gene_profiles  — per-gene analysis (diplotype, phenotype, reasoning)
          drug_results   — per-drug risk prediction
    """
    if not variants and not drugs:
        logger.warning("predict_multi_drug called with empty variants and drugs")

    # Validate inputs
    validated_drugs: List[str] = []
    skipped_drugs: List[str] = []
    for d in drugs:
        key = _normalise_drug(d)
        if key in SUPPORTED_DRUGS:
            validated_drugs.append(d)
        else:
            logger.warning(f"Drug '{d}' is not supported — skipped")
            skipped_drugs.append(d)

    # Step 1: Analyse all genes once
    gene_results: Dict[str, GeneResult] = {}
    for gene in TARGET_GENES:
        gene_results[gene] = analyse_gene(gene, variants)

    # Step 2: Predict risk for each validated drug
    drug_results: List[Dict] = []
    for drug in validated_drugs:
        drug_risk = predict_drug_risk(drug, gene_results)
        drug_results.append(asdict(drug_risk))

    # Step 3: Add Unknown entries for unsupported drugs
    for drug in skipped_drugs:
        drug_results.append({
            "drug": drug,
            "gene_used": "N/A",
            "phenotype": INDETERMINATE,
            "risk_label": "Unknown",
            "severity": "unknown",
            "confidence_score": 0.0,
            "clinical_action": f"'{drug}' is not in the supported drug list.",
            "cpic_guideline": "N/A",
            "supporting_variants": [],
            "reasoning": f"Drug '{drug}' is not supported by this version of PharmaGuard.",
            "evidence_strength": "none",
        })

    return {
        "gene_profiles": {
            gene: asdict(gr) for gene, gr in gene_results.items()
        },
        "drug_results": drug_results,
        "skipped_drugs": skipped_drugs,
    }
