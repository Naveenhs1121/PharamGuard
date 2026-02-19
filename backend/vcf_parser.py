import vcfpy

TARGET_GENES = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]

def extract_variants(vcf_path):
    reader = vcfpy.Reader.from_path(vcf_path)

    variants = []

    for record in reader:
        gene_info = record.INFO.get("GENE")

        if gene_info:
            gene = gene_info[0] if isinstance(gene_info, list) else gene_info

            if gene in TARGET_GENES:
                variant_info = {
                    "gene": gene,
                    "chrom": record.CHROM,
                    "pos": record.POS,
                    "rsid": record.ID,
                    "ref": record.REF,
                    "alt": [str(a) for a in record.ALT]
                }
                variants.append(variant_info)

    return variants
    