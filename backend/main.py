from fastapi import FastAPI, UploadFile, Form, File
import shutil
from vcf_parser import extract_variants
from risk_engine import predict_risk
from explanation_engine import generate_explanation
from datetime import datetime

app = FastAPI()

@app.get("/")
def home():
    return {"message": "PharmaGuard backend running"}

@app.post("/analyze")
async def analyze(file: UploadFile = File(...), drug: str = Form(...)):
    file_location = f"temp_{file.filename}"

    with open(file_location, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)

    variants = []

    if file.filename.endswith(".vcf"):
        variants = extract_variants(file_location)

    risk = predict_risk(variants, drug)
    explanation = generate_explanation(variants, drug, risk)

    return {
        "patient_id": "PATIENT_DEMO",
        "drug": drug,
        "timestamp": datetime.utcnow().isoformat(),

        "risk_assessment": {
            "risk_label": risk["label"],
            "severity": risk["severity"],
            "confidence_score": risk["confidence"]
        },

        "pharmacogenomic_profile": {
            "detected_variants": variants
        },

        "explanation": explanation,

        "quality_metrics": {
            "vcf_parsing_success": True
        }
    }
