# Step 2: Extract SNPs from reference panels using PLINK
# This extracts the specified SNPs from each reference panel (EUR, CHI, Hisp)

rule extract_snps:
    input:
        snplist = "snpList_unix.txt"
    output:
        bed = f"{OUTPUT_DIR}/{{panel}}.bed",
        bim = f"{OUTPUT_DIR}/{{panel}}.bim",
        fam = f"{OUTPUT_DIR}/{{panel}}.fam",
        log = f"{OUTPUT_DIR}/{{panel}}.log"
    params:
        bfile = lambda wildcards: PANELS[wildcards.panel],
        out_prefix = f"{OUTPUT_DIR}/{{panel}}"
    log:
        "logs/02_extract_snps_{panel}.log"
    resources:
        mem_mb = 64000
    shell:
        """
        ml plink2/1.90b3w 
        plink --bfile {params.bfile} \
              --extract {input.snplist} \
              --make-bed \
              --out {params.out_prefix} \
              --chr 1-22 \
              --allow-extra-chr 2> {log}
        
        echo "Extraction completed for {wildcards.panel}" >> {log}
        """
