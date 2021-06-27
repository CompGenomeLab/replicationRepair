
rule organize:
    input:
        muts="results/mutation/{samples}/{samples}_subs.tsv",
        genome="resources/ref_genomes/hg19/genome_hg19.fa",
    output:
        org="results/mutation/{samples}/{samples}_subs_org.tsv",
        fa="results/mutation/{samples}/{samples}_subs_org.fa",
        bed="results/mutation/{samples}/{samples}.bed",     
    log:
        "logs/{samples}/{samples}_organize.log",
    benchmark:
        "logs/{samples}/{samples}_organize.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Organize and remove chrY, chrMT..." &&
        grep -v -e "chrY" -e "chrMT" {input.muts} | \
        awk '{{print $2"\\t"($3-2)"\\t"($4+1)"\\t"$6"_"$7"::"$2":"($3-2)"-"($4+1)"\\t""0""\\t""+"}}' > \
        {output.org} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        
        (echo "`date -R`: Creating fasta to retrieve mutated nucleotides..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.org} \
        -fo {output.fa} \
        -name &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1    

        (echo "`date -R`: Creating final bed..." &&
        awk -F'[>_:-]' \
        'NR%2==1{{a=$5"\\t"$6"\\t"$7"\\t"$2"_"$3}} NR%2==0{{print a"_"$0}}' \
        {output.fa} | \
        sort -k1,1 -k2,2n -k3,3n > \
        {output.bed} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1    
        """
