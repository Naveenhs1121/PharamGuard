def generate_explanation(variants, drug, risk):
    if not variants:
        return "No pharmacogenomic variants affecting this drug."

    genes = [v["gene"] for v in variants]

    return f"Detected variants in {', '.join(genes)} influence response to {drug}, leading to {risk['label']} risk."
