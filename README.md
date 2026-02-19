# PharmaGuard — Pharmacogenomics (PGx) Risk Intelligence

![PharmaGuard Banner](https://via.placeholder.com/1200x400.png?text=PharmaGuard+Banner)

> **Know Your Genes. Know Your Risks.**  
> A Next-Gen Pharmacogenomic Risk Predictor aligned with CIPC Guidelines.

---

## 🚀 Live Demo & Video
- **🌐 Live Application**: [https://pharmaguard-demo.vercel.app](https://pharmaguard-demo.vercel.app) *(Replace with your deployed URL)*
- **🎥 LinkedIn Demo**: [Watch the Video](https://linkedin.com) *(Replace with your video link)*

---

## 🧬 Problem Statement
Adverse drug reactions (ADRs) are a leading cause of hospitalizations globally. Many of these are preventable through **Pharmacogenomics (PGx)** — tailoring drug prescriptions to an individual's genetic makeup.

**PharmaGuard** solves the complexity of PGx by:
1.  **Parsing** raw genetic data (VCF files).
2.  **Predicting** drug risks using CPIC (Clinical Pharmacogenomics Implementation Consortium) guidelines.
3.  **Explaining** the "Why" using AI, bridging the gap between complex genetics and clinical decision-making.

---

## 🏗️ Architecture

PharmaGuard employs a **Hybrid Intelligence** approach:

1.  **Deterministic Risk Engine (Safety Layer)**:
    - **100% Rule-Based**: Uses hardcoded CPIC allele tables and phenotype rules.
    - **Zero Hallucination**: Risk labels (e.g., "Toxic", "Safe") are never guessed by AI.
2.  **Generative AI Layer (Explanation Layer)**:
    - Uses **Google Gemini 1.5 Flash** to translate technical diplotypes into patient-friendly explanations.

### Tech Stack
- **Backend**: Python 3.11, FastAPI
- **Frontend**: Vanilla HTML5/CSS3 (Obsidian Theme), JavaScript (ES6+)
- **Genomics**: `vcfpy` for parsing
- **AI**: Google Gemini API

---

## 🛠️ Installation & Setup

### Prerequisites
- Python 3.9+
- Google Gemini API Key (Optional, for AI explanations)

### Steps
1.  **Clone the Repository**
    ```bash
    git clone https://github.com/YourUsername/PharmaGuard.git
    cd PharmaGuard
    ```

2.  **Install Dependencies**
    ```bash
    pip install -r backend/requirements.txt
    ```

3.  **Configure Environment**
    Create a `.env` file in the root directory:
    ```env
    GEMINI_API_KEY=your-gemini-api-key-here
    ```

4.  **Run the Application**
    ```bash
    uvicorn backend.main:app --reload
    ```

5.  **Access the App**
    Open your browser and navigate to: `http://127.0.0.1:8000`

---

## 📖 Usage Guide

1.  **Upload VCF**: Drag & drop your VCF file (e.g., `sample_vcf/test.vcf`).
2.  **Select Drugs**: Choose from the supported drug list (Codeine, Warfarin, etc.).
3.  **Analyze**: Click "Analyse Pharmacogenomic Risk".
4.  **View Results**:
    - **Risk Cards**: Color-coded cards showing risk levels.
    - **AI Explanation**: Expansion panel with clinical reasoning.
    - **Genomic Summary**: Technical details of diplotypes and phenotypes.
5.  **Export**: Click "⬇️ Download Report JSON" to save the analysis.

---

## 🔌 API Documentation

### `POST /analyze/multi`
Analyzes a VCF file against a list of drugs.

**Request**: `multipart/form-data`
- `file`: .vcf file
- `drugs`: Comma-separated string (e.g., "codeine,warfarin")

**Response**:
```json
{
  "patient_id": "PATIENT_DEMO",
  "timestamp": "2026-02-20T10:00:00Z",
  "drug_reports": [
    {
      "drug": "codeine",
      "risk_assessment": {
        "risk_label": "Toxic",
        "severity": "high",
        "confidence_score": 0.95
      },
      "pharmacogenomic_profile": {
        "primary_gene": "CYP2D6",
        "diplotype": "*1/*1",
        "phenotype": "Ultrarapid Metabolizer",
        "detected_variants": []
      },
      "clinical_recommendation": {
        "action": "Avoid codeine...",
        "cpic_guideline": "CPIC Codeine Guideline..."
      },
      "llm_generated_explanation": {
        "summary": "Patient is an Ultrarapid Metabolizer..."
      },
      "quality_metrics": { ... }
    }
  ],
  ...
}
```

---

## 👥 Team
- **Your Name** - Full Stack Developer & AI Engineer

---
*Built for RIFT 2026 Hackathon*
