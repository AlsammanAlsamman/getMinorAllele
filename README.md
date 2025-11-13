# Minor Allele Frequency Pipeline

This Snakemake pipeline calculates minor allele frequencies for SNPs across multiple reference panels.

## Pipeline Structure

The pipeline is organized into modular steps:

### Step 1: Prepare SNP List (`rules/01_prepare_snplist.smk`)
- Converts the SNP list to Unix format (removes Windows line endings)
- Input: `snpList.txt`
- Output: `snpList_unix.txt`

### Step 2: Extract SNPs (`rules/02_extract_snps.smk`)
- Extracts specified SNPs from reference panels using PLINK
- Input: `snpList_unix.txt`, reference panels from `repanels.yml`
- Output: PLINK binary files (`.bed`, `.bim`, `.fam`) for each panel in `extracted_data/`

## Usage

Run the entire pipeline:
```bash
snakemake --cores 1
```

Run specific steps:
```bash
# Only prepare SNP list
snakemake snpList_unix.txt --cores 1

# Only extract SNPs for a specific panel
snakemake extracted_data/EUR.bed --cores 1
```

Dry run to see what will be executed:
```bash
snakemake -n
```

## Reference Panels

Configured in `repanels.yml`:
- EUR: European
- CHI: East Asian (EAS)
- Hisp: Hispanic/Latino (AMR)

## Directory Structure

```
.
├── Snakefile                 # Main pipeline file
├── repanels.yml             # Reference panel paths
├── snpList.txt              # Input SNP list
├── rules/                   # Modular step definitions
│   ├── 01_prepare_snplist.smk
│   └── 02_extract_snps.smk
├── logs/                    # Log files
├── extracted_data/          # Output directory (created by pipeline)
└── README.md               # This file
```
