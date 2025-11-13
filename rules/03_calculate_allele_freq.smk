# Step 3: Calculate allele frequencies for both alleles
# This calculates the frequency of both alleles (major and minor) for each SNP

rule calculate_allele_freq:
    input:
        bed = f"{OUTPUT_DIR}/{{panel}}.bed",
        bim = f"{OUTPUT_DIR}/{{panel}}.bim",
        fam = f"{OUTPUT_DIR}/{{panel}}.fam"
    output:
        freq = f"{OUTPUT_DIR}/{{panel}}.frq"
    params:
        bfile = f"{OUTPUT_DIR}/{{panel}}"
    log:
        "logs/03_calculate_allele_freq_{panel}.log"
    resources:
        mem_mb = 8000
    shell:
        """
        ml plink2/1.90b3w 
        plink --bfile {params.bfile} \
              --freq \
              --out {params.bfile} 2> {log}
        
        echo "Allele frequency calculation completed for {wildcards.panel}" >> {log}
        echo "Output file: {output.freq}" >> {log}
        """
