# Step 1: Convert SNP list to Unix format (dos2unix functionality)
# This removes Windows line endings (CRLF) and converts to Unix format (LF)

rule prepare_snplist:
    input:
        snplist = "snpList.txt"
    output:
        unix_snplist = "snpList_unix.txt"
    log:
        "logs/01_prepare_snplist.log"
    shell:
        """
        # Convert DOS/Windows line endings to Unix format
        # Using tr to remove carriage returns
        tr -d '\\r' < {input.snplist} > {output.unix_snplist} 2> {log}
        
        echo "SNP list converted to Unix format" >> {log}
        echo "Total SNPs: $(wc -l < {output.unix_snplist})" >> {log}
        """
