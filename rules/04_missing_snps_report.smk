# Step 4: Generate report on missing SNPs
# This creates a report showing how many SNPs are missing in each panel and overall

rule missing_snps_report:
    input:
        original_snplist = "snpList_unix.txt",
        freq_files = expand(f"{OUTPUT_DIR}/{{panel}}.frq", panel=PANEL_NAMES)
    output:
        report = "missing_snps_report.txt"
    params:
        panels = PANEL_NAMES,
        output_dir = OUTPUT_DIR
    log:
        "logs/04_missing_snps_report.log"
    run:
        import os
        
        # Read original SNP list
        with open(input.original_snplist, 'r') as f:
            original_snps = set(line.strip() for line in f if line.strip())
        
        total_snps = len(original_snps)
        
        # Dictionary to store found SNPs per panel
        panel_snps = {}
        
        # Read SNPs found in each panel
        for panel in params.panels:
            frq_file = f"{params.output_dir}/{panel}.frq"
            found_snps = set()
            
            with open(frq_file, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        snp_id = parts[1]  # SNP column
                        found_snps.add(snp_id)
            
            panel_snps[panel] = found_snps
        
        # Find SNPs missing in each panel
        panel_missing = {}
        for panel, found in panel_snps.items():
            panel_missing[panel] = original_snps - found
        
        # Find SNPs missing in ALL panels
        all_found = set.union(*panel_snps.values()) if panel_snps else set()
        missing_in_all = original_snps - all_found
        
        # Write report
        with open(output.report, 'w') as out:
            out.write("=" * 80 + "\n")
            out.write("MISSING SNPs REPORT\n")
            out.write("=" * 80 + "\n\n")
            
            out.write(f"Total SNPs requested: {total_snps}\n\n")
            
            # Per panel statistics
            out.write("-" * 80 + "\n")
            out.write("MISSING SNPs PER REFERENCE PANEL\n")
            out.write("-" * 80 + "\n\n")
            
            for panel in params.panels:
                found = len(panel_snps[panel])
                missing = len(panel_missing[panel])
                missing_pct = (missing / total_snps * 100) if total_snps > 0 else 0
                
                out.write(f"Panel: {panel}\n")
                out.write(f"  Found: {found} SNPs\n")
                out.write(f"  Missing: {missing} SNPs ({missing_pct:.2f}%)\n")
                
                if missing > 0 and missing <= 20:  # List missing SNPs if <= 20
                    out.write(f"  Missing SNPs: {', '.join(sorted(panel_missing[panel]))}\n")
                elif missing > 20:
                    out.write(f"  Missing SNPs: {', '.join(sorted(list(panel_missing[panel])[:10]))}...\n")
                    out.write(f"  (showing first 10 of {missing} missing SNPs)\n")
                out.write("\n")
            
            # Overall missing
            out.write("-" * 80 + "\n")
            out.write("SNPs MISSING IN ALL PANELS\n")
            out.write("-" * 80 + "\n\n")
            
            missing_all_count = len(missing_in_all)
            missing_all_pct = (missing_all_count / total_snps * 100) if total_snps > 0 else 0
            
            out.write(f"Total SNPs missing in ALL panels: {missing_all_count} ({missing_all_pct:.2f}%)\n")
            
            if missing_all_count > 0:
                out.write(f"\nMissing SNPs:\n")
                for snp in sorted(missing_in_all):
                    out.write(f"  - {snp}\n")
            else:
                out.write("\nAll SNPs were found in at least one panel!\n")
            
            out.write("\n" + "=" * 80 + "\n")
            out.write("END OF REPORT\n")
            out.write("=" * 80 + "\n")
        
        # Log summary
        with open(log[0], 'w') as log_file:
            log_file.write(f"Missing SNPs report generated successfully\n")
            log_file.write(f"Total SNPs: {total_snps}\n")
            log_file.write(f"SNPs missing in all panels: {missing_all_count}\n")
            for panel in params.panels:
                log_file.write(f"{panel}: {len(panel_missing[panel])} missing\n")
