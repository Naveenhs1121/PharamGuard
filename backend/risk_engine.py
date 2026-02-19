def predict_risk(variants, drug):
    drug = drug.lower()

    for v in variants:
        gene = v["gene"]

        if drug == "warfarin" and gene == "CYP2C9":
            return {"label": "Adjust Dosage", "severity": "moderate", "confidence": 0.8}

        if drug == "clopidogrel" and gene == "CYP2C19":
            return {"label": "Ineffective", "severity": "high", "confidence": 0.9}

        if drug == "simvastatin" and gene == "SLCO1B1":
            return {"label": "Toxic", "severity": "high", "confidence": 0.85}

        if drug == "codeine" and gene == "CYP2D6":
            return {"label": "Adjust Dosage", "severity": "moderate", "confidence": 0.8}

    return {"label": "Safe", "severity": "none", "confidence": 0.7}
