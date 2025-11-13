# Main Snakemake pipeline for minor allele frequency calculation
import yaml

# Load reference panel paths from YAML
with open("repanels.yml", "r") as f:
    PANELS = yaml.safe_load(f)

# Define output directory
OUTPUT_DIR = "extracted_data"

# Get panel names
PANEL_NAMES = list(PANELS.keys())

# Include modular snakemake files
include: "rules/01_prepare_snplist.smk"
include: "rules/02_extract_snps.smk"
include: "rules/03_calculate_allele_freq.smk"
include: "rules/04_missing_snps_report.smk"

# Main rule to run the entire pipeline
rule all:
    input:
        # Step 1: Unix-formatted SNP list
        "snpList_unix.txt",
        # Step 2: Extracted PLINK files for each panel
        expand(f"{OUTPUT_DIR}/{{panel}}.bed", panel=PANEL_NAMES),
        expand(f"{OUTPUT_DIR}/{{panel}}.bim", panel=PANEL_NAMES),
        expand(f"{OUTPUT_DIR}/{{panel}}.fam", panel=PANEL_NAMES),
        expand(f"{OUTPUT_DIR}/{{panel}}.log", panel=PANEL_NAMES),
        # Step 3: Allele frequency files for each panel
        expand(f"{OUTPUT_DIR}/{{panel}}.frq", panel=PANEL_NAMES),
        # Step 4: Missing SNPs report
        "missing_snps_report.txt"
