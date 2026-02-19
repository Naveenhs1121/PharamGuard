"""
PharmaGuard — LLM Explanation Engine
Generates clinical explanations using Google Gemini (primary)
with a structured template fallback if LLM is unavailable.
"""

import os
import logging
from typing import Dict, List, Any, Optional

logger = logging.getLogger("PharmaGuard.ExplanationEngine")

# ---------------------------------------------------------------------------
# LLM setup — Google Gemini
# ---------------------------------------------------------------------------
_gemini_configured = False

def _configure_gemini():
    """Lazy-initialise Gemini client (only once)."""
    global _gemini_configured
    if _gemini_configured:
        return True

    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        logger.warning("GEMINI_API_KEY not set — LLM explanations disabled, using fallback.")
        return False

    try:
        import google.generativeai as genai
        genai.configure(api_key=api_key)
        _gemini_configured = True
        logger.info("Google Gemini client configured successfully.")
        return True
    except Exception as e:
        logger.error(f"Failed to configure Gemini client: {e}")
        return False


# ---------------------------------------------------------------------------
# Prompt template
# ---------------------------------------------------------------------------
CLINICAL_PROMPT = """You are a pharmacogenomics clinical assistant.
Given the detected gene variants, inferred diplotype, phenotype, and predicted drug risk, generate a concise clinical explanation describing how the genetic variation affects drug metabolism, efficacy, or toxicity.

Explain the biological mechanism, mention the gene involved, and connect phenotype to the predicted drug response. Keep the explanation medically accurate, structured, and understandable for clinicians.

Input:
Gene: {gene}
Diplotype: {diplotype}
Phenotype: {phenotype}
Drug: {drug}
Risk: {risk_label}
Detected variants: {variant_list}"""


def _build_prompt(gene: str, diplotype: str, phenotype: str,
                  drug: str, risk_label: str, variant_list: List[str]) -> str:
    """Fill the clinical prompt template with provided values."""
    variants_str = ", ".join(variant_list) if variant_list else "None detected"
    return CLINICAL_PROMPT.format(
        gene=gene,
        diplotype=diplotype,
        phenotype=phenotype,
        drug=drug,
        risk_label=risk_label,
        variant_list=variants_str,
    )


# ---------------------------------------------------------------------------
# LLM call
# ---------------------------------------------------------------------------
def _call_gemini(prompt: str) -> Optional[str]:
    """Call Gemini and return the response text, or None on failure."""
    if not _configure_gemini():
        return None
    try:
        import google.generativeai as genai
        model = genai.GenerativeModel("gemini-1.5-flash")
        response = model.generate_content(prompt)
        text = response.text.strip()
        logger.info("Gemini explanation generated successfully.")
        return text
    except Exception as e:
        logger.error(f"Gemini API call failed: {e}")
        return None


# ---------------------------------------------------------------------------
# Fallback — structured template explanation
# ---------------------------------------------------------------------------
PHENOTYPE_DESC = {
    "Poor Metabolizer":          "cannot efficiently metabolize this drug due to severely reduced enzyme activity",
    "Intermediate Metabolizer":  "has partially reduced enzyme activity, leading to slower drug metabolism than normal",
    "Normal Metabolizer":        "metabolizes this drug at a standard rate with no expected pharmacogenomic interaction",
    "Rapid Metabolizer":         "metabolizes this drug faster than average, which may affect drug plasma levels",
    "Ultrarapid Metabolizer":    "metabolizes this drug at an exceptionally high rate, risking sub-therapeutic levels or toxicity",
    "Indeterminate":             "has an uncertain metabolizer status based on the available genomic data",
}

SEVERITY_CONTEXT = {
    "none":     "No clinically significant pharmacogenomic interaction is expected.",
    "low":      "A minor pharmacogenomic interaction is noted; routine monitoring is advised.",
    "moderate": "A clinically significant interaction exists; dose adjustment is recommended.",
    "high":     "A serious pharmacogenomic interaction is identified; immediate clinical action is required.",
    "unknown":  "The interaction risk is uncertain due to insufficient genomic data.",
}


def _fallback_explanation(gene: str, diplotype: str, phenotype: str,
                           drug: str, risk_label: str,
                           variant_list: List[str],
                           severity: str, clinical_action: str) -> str:
    """Rule-based clinical explanation when LLM is unavailable."""
    phenotype_desc = PHENOTYPE_DESC.get(phenotype, "has an atypical metabolizer status")
    severity_ctx   = SEVERITY_CONTEXT.get(severity, SEVERITY_CONTEXT["unknown"])

    variant_text = (
        f"Detected variants ({', '.join(variant_list)}) mapped to diplotype {diplotype}."
        if variant_list
        else f"No annotated variants found in {gene}; reference diplotype {diplotype} assumed."
    )

    parts = [
        f"Pharmacogenomic analysis of {drug.capitalize()} based on {gene} genotyping:",
        f"The patient {phenotype_desc} ({phenotype}). {variant_text}",
        f"Risk classification: '{risk_label}'. {severity_ctx}",
    ]
    if clinical_action:
        parts.append(f"Clinical recommendation: {clinical_action}")
    parts.append(
        "This assessment follows CPIC (Clinical Pharmacogenomics Implementation Consortium) "
        "guidelines and should be interpreted alongside the patient's full clinical context."
    )
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def generate_explanation(
    variants: List[Dict],
    drug: str,
    risk: Dict[str, Any],
    gene_profile: Optional[Dict] = None,
    clinical_action: str = "",
) -> str:
    """
    Generate a clinical explanation for a drug-gene interaction.

    Tries Gemini LLM first; falls back to structured template if unavailable.

    Args:
        variants        : list of VCF variant dicts
        drug            : drug name
        risk            : dict with keys label, severity, confidence
        gene_profile    : gene analysis result dict (diplotype, phenotype, rsIDs...)
        clinical_action : CPIC-recommended action string

    Returns:
        Clinical explanation string.
    """
    risk_label = risk.get("label", "Unknown")
    severity   = risk.get("severity", "unknown")

    # Extract gene profile fields
    gene       = (gene_profile or {}).get("gene", "Unknown gene")
    diplotype  = (gene_profile or {}).get("diplotype", "Unknown")
    phenotype  = (gene_profile or {}).get("phenotype", "Indeterminate")
    rsids      = (gene_profile or {}).get("detected_rsids", [])

    # Handle no-variant case early
    if not variants:
        return (
            f"No pharmacogenomic variants relevant to {drug.capitalize()} were detected "
            "in the uploaded VCF file. In the absence of known risk variants, standard "
            f"{drug.capitalize()} dosing is generally appropriate. Clinical judgment "
            "should guide prescribing decisions."
        )

    # Build prompt and try LLM
    prompt = _build_prompt(
        gene=gene,
        diplotype=diplotype,
        phenotype=phenotype,
        drug=drug,
        risk_label=risk_label,
        variant_list=rsids,
    )

    llm_response = _call_gemini(prompt)
    if llm_response:
        return llm_response

    # Fallback
    logger.info("Using template-based fallback explanation.")
    return _fallback_explanation(
        gene=gene,
        diplotype=diplotype,
        phenotype=phenotype,
        drug=drug,
        risk_label=risk_label,
        variant_list=rsids,
        severity=severity,
        clinical_action=clinical_action,
    )
