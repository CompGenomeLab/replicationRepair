rule get_sbs_muts:
    input:
        "resources/samples/mutation/{samples}.tsv.gz",
    output:
        unzip=temp("results/mutation/{samples}/{samples}.tsv"),
        subs="results/mutation/{samples}/{samples}_subs.tsv",     
    log:
        "logs/{samples}/{samples}_get_sbs_muts.log",
    benchmark:
        "logs/{samples}/{samples}_get_sbs_muts.benchmark.txt",
    shell:
        """
        (echo "`date -R`: Unzipping..." &&
        zcat {input} | \
        sed 's/single base sub/single_base_sub/g' > \
        {output.unzip} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1
        
        (echo "`date -R`: Extracting the subs..." &&
        awk -v t="\\t" '{{print $2t"chr"$9t$10t$11t$14t$16t$17}}' \
        {output.unzip} | \
        sort -k1,1 -k2,2 -k3,3n -k4,4n | \
        uniq | \
        grep single_base_sub > \
        {output.subs} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1    
        """
