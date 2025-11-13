# Step 7: Create final Excel report with annotations and MAF
# Combines GWAS data, SNP annotations, and MAF from all reference panels

rule create_final_report:
    input:
        gwas_best = "gwas_snp_best_pvalue.tsv",
        annotations = "snp_annotations.txt",
        freq_eur = f"{OUTPUT_DIR}/EUR.frq",
        freq_chi = f"{OUTPUT_DIR}/CHI.frq",
        freq_hisp = f"{OUTPUT_DIR}/Hisp.frq"
    output:
        excel = "final_snp_report.xlsx"
    log:
        "logs/07_create_final_report.log"
    shell:
        """
        Rscript - << 'REOF'
# Load required libraries
suppressPackageStartupMessages({{
    library(openxlsx)
    library(dplyr)
}})

cat("Reading input files...\\n")

# Read GWAS best p-value data
gwas <- read.table("{input.gwas_best}", header=TRUE, sep="\\t", stringsAsFactors=FALSE, quote="")

# Read SNP annotations
annotations <- read.table("{input.annotations}", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

# Read frequency files for each panel
freq_eur <- read.table("{input.freq_eur}", header=TRUE, stringsAsFactors=FALSE)
freq_chi <- read.table("{input.freq_chi}", header=TRUE, stringsAsFactors=FALSE)
freq_hisp <- read.table("{input.freq_hisp}", header=TRUE, stringsAsFactors=FALSE)

cat("Processing MAF data...\\n")

# Function to get MAF aligned with effect allele
get_aligned_maf <- function(snp, ea, freq_data) {{
    row <- freq_data[freq_data$SNP == snp, ]
    if (nrow(row) == 0) return(NA)
    
    # Check if A1 matches EA
    if (row$A1 == ea) {{
        return(row$MAF)
    }} else {{
        # A1 doesn't match EA, so return 1-MAF
        return(1 - row$MAF)
    }}
}}

# Add MAF columns for each panel
gwas$MAF_EUR <- sapply(1:nrow(gwas), function(i) {{
    get_aligned_maf(gwas$SNPID[i], gwas$EA[i], freq_eur)
}})

gwas$MAF_CHI <- sapply(1:nrow(gwas), function(i) {{
    get_aligned_maf(gwas$SNPID[i], gwas$EA[i], freq_chi)
}})

gwas$MAF_Hisp <- sapply(1:nrow(gwas), function(i) {{
    get_aligned_maf(gwas$SNPID[i], gwas$EA[i], freq_hisp)
}})

# Calculate average MAF (ignoring NAs)
gwas$MAF_Average <- rowMeans(gwas[, c("MAF_EUR", "MAF_CHI", "MAF_Hisp")], na.rm=TRUE)

# Create SNP_ID column in chr:pos format
gwas$SNP_ID <- paste0(gwas$CHR, ":", gwas$POS)

cat("Merging with annotations...\\n")

# Merge with annotations (left join - keep all GWAS SNPs)
final_data <- merge(gwas, annotations, by.x="SNPID", by.y="SNP", all.x=TRUE)

# Reorder columns (add SNP_ID after SNPID)
final_data <- final_data[, c("SNPID", "SNP_ID", "gene", "Locus", "CHR", "POS", "Start", "End",
                              "EA", "NEA", "BETA", "SE", "P", "OR", 
                              "MAF_EUR", "MAF_CHI", "MAF_Hisp", "MAF_Average", "SOURCE")]

# Sort by chromosome and position
final_data <- final_data[order(as.numeric(final_data$CHR), as.numeric(final_data$POS)), ]

cat("Creating Excel workbook...\\n")

# Create workbook
wb <- createWorkbook()
addWorksheet(wb, "SNP_Report")

# Write data
writeData(wb, "SNP_Report", final_data, startRow=1, startCol=1)

# Format header
headerStyle <- createStyle(
    fontSize = 12,
    fontColour = "#FFFFFF",
    halign = "center",
    valign = "center",
    textDecoration = "bold",
    fgFill = "#4F81BD",
    border = "TopBottomLeftRight",
    borderColour = "#FFFFFF"
)

addStyle(wb, "SNP_Report", headerStyle, rows=1, cols=1:ncol(final_data), gridExpand=TRUE)

# Format P-value column as scientific notation
pvalueStyle <- createStyle(numFmt = "0.00E+00")
addStyle(wb, "SNP_Report", pvalueStyle, rows=2:(nrow(final_data)+1), cols=which(names(final_data)=="P"), gridExpand=TRUE)

# Format numeric columns
numericStyle <- createStyle(numFmt = "0.0000")
numeric_cols <- c("BETA", "SE", "OR", "MAF_EUR", "MAF_CHI", "MAF_Hisp", "MAF_Average")
for (col_name in numeric_cols) {{
    col_idx <- which(names(final_data) == col_name)
    if (length(col_idx) > 0) {{
        addStyle(wb, "SNP_Report", numericStyle, rows=2:(nrow(final_data)+1), cols=col_idx, gridExpand=TRUE)
    }}
}}

# Auto-size columns
setColWidths(wb, "SNP_Report", cols=1:ncol(final_data), widths="auto")

# Freeze first row
freezePane(wb, "SNP_Report", firstRow=TRUE)

# Add filters
addFilter(wb, "SNP_Report", rows=1, cols=1:ncol(final_data))

# Save workbook
saveWorkbook(wb, "{output.excel}", overwrite=TRUE)

cat(paste0("\\nFinal report created: {output.excel}\\n"))
cat(paste0("Total SNPs in report: ", nrow(final_data), "\\n"))

# Write log
sink("{log}")
cat("Final report generation completed\\n")
cat(paste0("Total SNPs: ", nrow(final_data), "\\n"))
cat(paste0("SNPs with annotations: ", sum(!is.na(final_data$gene)), "\\n"))
cat(paste0("\\nMAF availability:\\n"))
cat(paste0("  EUR: ", sum(!is.na(final_data$MAF_EUR)), "\\n"))
cat(paste0("  CHI: ", sum(!is.na(final_data$MAF_CHI)), "\\n"))
cat(paste0("  Hisp: ", sum(!is.na(final_data$MAF_Hisp)), "\\n"))
cat(paste0("\\nP-value range:\\n"))
cat(paste0("  Min: ", min(final_data$P, na.rm=TRUE), "\\n"))
cat(paste0("  Max: ", max(final_data$P, na.rm=TRUE), "\\n"))
sink()

REOF
        """
