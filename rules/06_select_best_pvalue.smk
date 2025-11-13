# Step 6: Keep only the best (lowest) p-value for each SNP
# For each SNP, select the record with the lowest p-value across all sources

rule select_best_pvalue:
    input:
        gwas_info = "gwas_snp_info.tsv"
    output:
        best_pvalue = "gwas_snp_best_pvalue.tsv"
    log:
        "logs/06_select_best_pvalue.log"
    shell:
        """
        python3 << 'PYEOF'
import sys

# Read the harmonized GWAS data
snp_records = {{}}  # {{snpid: [record with best p-value]}}

print("Reading harmonized GWAS data...")

with open("{input.gwas_info}", 'r') as f:
    header = f.readline().strip()
    header_fields = header.split('\\t')
    
    # Find column indices
    snpid_idx = header_fields.index('SNPID')
    p_idx = header_fields.index('P')
    
    for line in f:
        fields = line.strip().split('\\t')
        snpid = fields[snpid_idx]
        
        try:
            p_value = float(fields[p_idx])
        except:
            # Skip records with invalid p-values
            continue
        
        # Check if this is the best (lowest) p-value for this SNP
        if snpid not in snp_records:
            snp_records[snpid] = (p_value, fields)
        else:
            current_best_p = snp_records[snpid][0]
            if p_value < current_best_p:
                snp_records[snpid] = (p_value, fields)
                print(f"  Updated best p-value for {{snpid}}: {{p_value}} < {{current_best_p}}")

print(f"\\nProcessed {{len(snp_records)}} unique SNPs")

# Write output - keep only the best record for each SNP
with open("{output.best_pvalue}", 'w') as out:
    out.write(header + '\\n')
    
    # Sort by SNP ID
    for snpid in sorted(snp_records.keys()):
        best_p, best_record = snp_records[snpid]
        out.write('\\t'.join(best_record) + '\\n')

print(f"\\nWrote {{len(snp_records)}} records with best p-values to output")

# Write log
with open("{log}", 'w') as log_file:
    log_file.write(f"Best p-value selection completed\\n")
    log_file.write(f"Total unique SNPs: {{len(snp_records)}}\\n")
    log_file.write(f"\\nP-value range:\\n")
    
    all_pvalues = [p for p, _ in snp_records.values()]
    if all_pvalues:
        log_file.write(f"  Minimum: {{min(all_pvalues)}}\\n")
        log_file.write(f"  Maximum: {{max(all_pvalues)}}\\n")
    
    # Count sources
    log_file.write(f"\\nBest records by source:\\n")
    source_counts = {{}}
    for _, record in snp_records.values():
        source = record[-1]  # SOURCE is last column
        source_counts[source] = source_counts.get(source, 0) + 1
    
    for source in sorted(source_counts.keys()):
        log_file.write(f"  {{source}}: {{source_counts[source]}}\\n")

PYEOF
        """
