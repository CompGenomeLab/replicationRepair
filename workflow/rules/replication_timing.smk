rule replication_timing:
    input:
        er1="results/edu/R21061297-EdUrep-2hrls_combined/R21061297-EdUrep-2hrls_combined_hg19.bdg",
        lr1="results/edu/R21061297-EdUrep-4hrls_combined/R21061297-EdUrep-4hrls_combined_hg19.bdg",
        euvr1="results/edu/R21061297-EdUrep-UV1-5hrls_combined/R21061297-EdUrep-UV1-5hrls_combined_hg19.bdg",
        luvr1="results/edu/R21061297-EdUrep-UV3-5hrls_combined/R21061297-EdUrep-UV3-5hrls_combined_hg19.bdg",
        er2="results/edu/R21071354-EdUrep2-2hrls2_combined/R21071354-EdUrep2-2hrls2_combined_hg19.bdg",
        lr2="results/edu/R21071354-EdUrep2-4hrls2_combined/R21071354-EdUrep2-4hrls2_combined_hg19.bdg",
        euvr2="results/edu/R21071354-EdUrep2-UV1-5hrls2_combined/R21071354-EdUrep2-UV1-5hrls2_combined_hg19.bdg",
        luvr2="results/edu/R21071354-EdUrep2-UV3-5hrls2_combined/R21071354-EdUrep2-UV3-5hrls2_combined_hg19.bdg",
        #euvr3="results/edu/R21102739-EdUrep3-UV1-5hrls3_combined/R21102739-EdUrep3-UV1-5hrls3_combined_hg19.bdg",
        #luvr3="results/edu/R21102739-EdUrep3-UV3-5hrls3_combined/R21102739-EdUrep3-UV3-5hrls3_combined_hg19.bdg",
    output:
        r1="results/edu/rep1_T.bg",
        r2="results/edu/rep2_T.bg",
        uvr1="results/edu/rep1_uv_T.bg",
        uvr2="results/edu/rep2_uv_T.bg",
        #uvr3="results/edu/rep3_uv_T.bg",
        untreated="results/edu/merge_RT.txt",
        uv="results/edu/merge_uv_RT.txt",
    log:
        "logs/rule/analysis/replication_timing.log",
    benchmark:
        "logs/rule/analysis/replication_timing.benchmark.txt",
    conda:
        "../envs/sambedtools.yaml"
    shell:  
        """
        (echo "`date -R`: Replication timing replicate 1..." &&
        paste {input.er1} {input.lr1} | awk '{{if($8 != 0 && $4 != 0){{print $1,$2,$3,log($4/$8)/log(2)}}}}' OFS='\\t' > {output.r1} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Replication timing replicate 2..." &&
        paste {input.er2} {input.lr2} | awk '{{if($8 != 0 && $4 != 0){{print $1,$2,$3,log($4/$8)/log(2)}}}}' OFS='\\t' > {output.r2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Replication timing combine replicates..." &&
        bedtools unionbedg -filler “NA” -i {output.r1} {output.r2} > {output.untreated} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Replication timing uv replicate 1..." &&
        paste {input.euvr1} {input.luvr1} | awk '{{if($8 != 0 && $4 != 0){{print $1,$2,$3,log($4/$8)/log(2)}}}}' OFS='\\t' > {output.uvr1} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Replication timing uv replicate 2..." &&
        paste {input.euvr2} {input.luvr2} | awk '{{if($8 != 0 && $4 != 0){{print $1,$2,$3,log($4/$8)/log(2)}}}}' OFS='\\t' > {output.uvr2} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        
        (echo "`date -R`: Replication timing combine replicates..." &&
        bedtools unionbedg -filler “NA” -i {output.uvr1} {output.uvr2} > {output.uv} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """