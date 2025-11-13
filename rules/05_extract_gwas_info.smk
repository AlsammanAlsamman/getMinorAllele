# Step 5: Extract GWAS information for each SNP from GWAS tables
# Extracts SNPID, CHR, POS, EA, NEA, BETA, SE, P, OR (calculated from BETA if missing)

rule extract_gwas_info:
    input:
        snplist = "snpList_unix.txt"
    output:
        gwas_info = "gwas_snp_info.tsv"
    params:
        gwas_dir = "gwas_tables"
    log:
        "logs/05_extract_gwas_info.log"
    shell:
        """
        python3 << 'PYEOF'
import sys
import math
import os
import glob

# Read SNP list
with open("{input.snplist}", 'r') as f:
    snps_to_extract = set(line.strip() for line in f if line.strip())

print(f"Looking for {{len(snps_to_extract)}} SNPs in GWAS tables...")

# Find all TSV files in gwas_tables directory
gwas_dir = "{params.gwas_dir}"
gwas_files = glob.glob(os.path.join(gwas_dir, "*.tsv"))

print(f"Found {{len(gwas_files)}} TSV files in {{gwas_dir}}/")

all_results = []

for table_file in gwas_files:
    # Get table name from filename (remove .tsv extension)
    table_name = os.path.basename(table_file).replace('.tsv', '')
    print(f"Processing {{table_name}}...")
    
    try:
        with open(table_file, 'r') as f:
            header = f.readline().strip().split('\\t')
            
            # Find column indices
            snpid_idx = header.index('SNPID')
            chr_idx = header.index('CHR')
            pos_idx = header.index('POS')
            ea_idx = header.index('EA')
            nea_idx = header.index('NEA')
            beta_idx = header.index('BETA')
            se_idx = header.index('SE')
            p_idx = header.index('P')
            or_idx = header.index('OR') if 'OR' in header else -1
            
            count = 0
            for line in f:
                fields = line.strip().split('\\t')
                snpid = fields[snpid_idx]
                
                if snpid in snps_to_extract:
                    # Get OR or calculate from BETA
                    if or_idx >= 0 and or_idx < len(fields) and fields[or_idx] and fields[or_idx] != 'NA':
                        or_val = fields[or_idx]
                    else:
                        # Calculate OR = exp(BETA)
                        try:
                            beta = float(fields[beta_idx])
                            or_val = str(math.exp(beta))
                        except:
                            or_val = 'NA'
                    
                    # Extract fields
                    result = [
                        fields[snpid_idx],
                        fields[chr_idx],
                        fields[pos_idx],
                        fields[ea_idx],
                        fields[nea_idx],
                        fields[beta_idx],
                        fields[se_idx],
                        fields[p_idx],
                        or_val,
                        table_name
                    ]
                    all_results.append(result)
                    count += 1
            
            print(f"  Found {{count}} SNPs in {{table_name}}")
    
    except Exception as e:
        print(f"Warning: Error processing {{table_file}}: {{e}}")

# Harmonize alleles - use first occurrence as reference
# Store reference alleles for each SNP
snp_reference = {{}}  # {{snpid: (ea, nea)}}
harmonized_results = []

print("\\nHarmonizing alleles across datasets...")

for result in sorted(all_results):
    snpid = result[0]
    ea = result[3]
    nea = result[4]
    beta = result[5]
    or_val = result[8]
    
    if snpid not in snp_reference:
        # First occurrence - set as reference
        snp_reference[snpid] = (ea, nea)
        harmonized_results.append(result)
    else:
        # Check if alleles need flipping
        ref_ea, ref_nea = snp_reference[snpid]
        
        if ea == ref_ea and nea == ref_nea:
            # Alleles match - no flipping needed
            harmonized_results.append(result)
        elif ea == ref_nea and nea == ref_ea:
            # Alleles are flipped - need to flip BETA and OR
            print(f"  Flipping alleles for {{snpid}} in {{result[9]}}")
            
            # Flip BETA (multiply by -1)
            try:
                flipped_beta = str(-float(beta))
            except:
                flipped_beta = beta
            
            # Flip OR (1/OR)
            try:
                flipped_or = str(1.0 / float(or_val))
            except:
                flipped_or = or_val
            
            # Create flipped result
            flipped_result = [
                result[0],  # SNPID
                result[1],  # CHR
                result[2],  # POS
                ref_ea,     # EA (use reference)
                ref_nea,    # NEA (use reference)
                flipped_beta,  # BETA (flipped)
                result[6],  # SE (unchanged)
                result[7],  # P (unchanged)
                flipped_or, # OR (flipped)
                result[9]   # SOURCE
            ]
            harmonized_results.append(flipped_result)
        else:
            # Alleles don't match at all - keep as is but warn
            print(f"  Warning: {{snpid}} in {{result[9]}} has incompatible alleles ({{ea}}/{{nea}} vs {{ref_ea}}/{{ref_nea}})")
            harmonized_results.append(result)

# Write output
with open("{output.gwas_info}", 'w') as out:
    out.write('SNPID\\tCHR\\tPOS\\tEA\\tNEA\\tBETA\\tSE\\tP\\tOR\\tSOURCE\\n')
    for result in harmonized_results:
        out.write('\\t'.join(result) + '\\n')

print(f"\\nTotal records extracted: {{len(harmonized_results)}}")
unique_snps = len(set(r[0] for r in harmonized_results))
print(f"Unique SNPs found: {{unique_snps}}")

# Write log
with open("{log}", 'w') as log_file:
    log_file.write(f"GWAS information extraction completed\\n")
    log_file.write(f"SNPs requested: {{len(snps_to_extract)}}\\n")
    log_file.write(f"Total records extracted: {{len(harmonized_results)}}\\n")
    log_file.write(f"Unique SNPs found: {{unique_snps}}\\n")
    log_file.write(f"Alleles harmonized based on first occurrence\\n")
    log_file.write(f"\\nRecords per source:\\n")
    sources = set(r[9] for r in harmonized_results)
    for source in sorted(sources):
        count = sum(1 for r in harmonized_results if r[9] == source)
        log_file.write(f"  {{source}}: {{count}}\\n")

PYEOF
        """
