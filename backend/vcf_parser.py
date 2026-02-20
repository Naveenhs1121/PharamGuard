import vcfpy

TARGET_GENES = ["CYP2D6", "CYP2C19", "CYP2C9", "SLCO1B1", "TPMT", "DPYD"]

def extract_variants(vcf_path):
    reader = vcfpy.Reader.from_path(vcf_path)

    variants = []

    for record in reader:
        # basic filtering
        if not record.ID:
            continue
            
        # Check genotype for the first sample (assuming single-sample VCF)
        if not record.calls:
            continue
            
        call = record.calls[0]
        gt = call.data.get("GT")
        
        # Skip if 0/0 (Hom-Ref), ./., or None
        # VCFPy usually returns GT as a string '0/0' or integer 0
        # If it's a Call object, we need to check is_variant or is_het/is_hom_alt
        
        # Using vcfpy's helper properties on the Call object if available, 
        # or manual check. Safest is manual check on the GT string/integers.
        
        # '0/0' or '0|0' -> Reference
        # './.' -> Missing
        
        # Let's filter:
        if not call.called:
            continue
        
        # If it is reference (all alleles are 0), skip
        # call.gt_type specific to vcfpy? Let's check call.data['GT']
        # Typically vcfpy returns raw GT value.
        
        # More robust check:
        # If no alleles are alternate (all 0), skip.
        has_alt = False
        if isinstance(gt, str):
            if "1" in gt or "2" in gt: # naive check for alt allele index
                has_alt = True
        elif isinstance(gt, (list, tuple)): # specific integers
            if any(x > 0 for x in gt if x is not None):
                has_alt = True
        
        if not has_alt:
             continue

        # Get gene info from INFO field
        gene_info = record.INFO.get("GENE")
        if not gene_info:
            continue
            
        gene = gene_info[0] if isinstance(gene_info, list) else gene_info

        if gene in TARGET_GENES:
            variant_info = {
                "gene": gene,
                "chrom": record.CHROM,
                "pos": record.POS,
                "rsid": record.ID,
                "ref": record.REF,
                "alt": [str(a) for a in record.ALT],
                "gt": str(gt)
            }
            variants.append(variant_info)

    return variants