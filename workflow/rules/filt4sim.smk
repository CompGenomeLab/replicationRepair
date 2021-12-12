rule filt4sim_xr_IZ:
    input:
        xr_plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        xr_minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
        region="results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100.bed",
    output:
        xr_plus_IZ="resources/samples/XR/{samples}_{build}_sorted_plus_IZ.bed",
        xr_minus_IZ="resources/samples/XR/{samples}_{build}_sorted_minus_IZ.bed",
        #xr_plus_rmTTTT="resources/samples/XR/{samples}_{build}_sorted_plus_rmTTTT.bed",
        #xr_minus_rmTTTT="resources/samples/XR/{samples}_{build}_sorted_minus_rmTTTT.bed",
    log:
        "logs/{samples}/{samples}_{build}_filt4sim_xr_IZ.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filt4sim_xr_IZ.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Getting reads from IZ region (xr plus strand)..." &&
        bedtools intersect \
        -a {input.xr_plus} \
        -b {input.region} \
        -f 0.5 > {output.xr_plus_IZ} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Getting reads from IZ region (minus strand)..." &&
        bedtools intersect \
        -a {input.xr_minus} \
        -b {input.region} \
        -v -f 0.5 > {output.xr_minus_IZ} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule filt4sim_ds_IZ:
    input:
        ds_plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        ds_minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        region="results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100.bed",
    output:
        ds_plus_IZ="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_IZ.bed",
        ds_minus_IZ="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_IZ.bed",
        #ds_plus_rmTTTT="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_rmTTTT.bed",
        #ds_minus_rmTTTT="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_rmTTTT.bed",
    log:
        "logs/{samples}/{samples}_{build}_filt4sim_ds_IZ.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filt4sim_ds_IZ.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Getting reads from IZ region (plus strand)..." &&
        bedtools intersect \
        -a {input.ds_plus} \
        -b {input.region} \
        -v -f 0.5 > {output.ds_plus_IZ} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Getting reads from IZ region (minus strand)..." &&
        bedtools intersect \
        -a {input.ds_minus} \
        -b {input.region} \
        -v -f 0.5 > {output.ds_minus_IZ} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule filt4sim_xr_rmTTTT:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        fa_plus=temp("resources/samples/XR/{samples}_{build}_sorted_plus_rmTTTT.fa"),
        fa_minus=temp("resources/samples/XR/{samples}_{build}_sorted_minus_rmTTTT.fa"),
        plus_rmTTTT="resources/samples/XR/{samples}_{build}_sorted_plus_rmTTTT.bed",
        minus_rmTTTT="resources/samples/XR/{samples}_{build}_sorted_minus_rmTTTT.bed",
    log:
        "logs/{samples}/{samples}_{build}_filt4sim_xr_rmTTTT.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filt4sim_xr_rmTTTT.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.fa_plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.fa_minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (plus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_plus} \
        -o {output.plus_rmTTTT} \
        -r 'TTTT'\
        -m "XR" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (minus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_minus} \
        -o {output.minus_rmTTTT} \
        -r 'TTTT' \
        -m "XR" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule filt4sim_ds_rmTTTT:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        fa_plus=temp("resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_rmTTTT.fa"),
        fa_minus=temp("resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_rmTTTT.fa"),
        plus_rmTTTT="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_rmTTTT.bed",
        minus_rmTTTT="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_rmTTTT.bed",
    log:
        "logs/{samples}/{samples}_{build}_filt4sim_ds_rmTTTT.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filt4sim_ds_rmTTTT.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.fa_plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.fa_minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (plus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_plus} \
        -o {output.plus_rmTTTT} \
        -r 'TTTT'\
        -m "DS" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (minus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_minus} \
        -o {output.minus_rmTTTT} \
        -r 'TTTT' \
        -m "DS" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule filt4sim_xr_IZ_rmTTTT:
    input:
        plus="resources/samples/XR/{samples}_{build}_sorted_plus_IZ.bed",
        minus="resources/samples/XR/{samples}_{build}_sorted_minus_IZ.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        fa_plus=temp("resources/samples/XR/{samples}_{build}_sorted_plus_IZ_rmTTTT.fa"),
        fa_minus=temp("resources/samples/XR/{samples}_{build}_sorted_minus_IZ_rmTTTT.fa"),
        plus_rmTTTT="resources/samples/XR/{samples}_{build}_sorted_plus_IZ_rmTTTT.bed",
        minus_rmTTTT="resources/samples/XR/{samples}_{build}_sorted_minus_IZ_rmTTTT.bed",
    log:
        "logs/{samples}/{samples}_{build}_filt4sim_xr_IZ_rmTTTT.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filt4sim_xr_IZ_rmTTTT.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.fa_plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.fa_minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (plus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_plus} \
        -o {output.plus_rmTTTT} \
        -r 'TTTT'\
        -m "XR" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (minus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_minus} \
        -o {output.minus_rmTTTT} \
        -r 'TTTT' \
        -m "XR" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule filt4sim_ds_IZ_rmTTTT:
    input:
        plus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_IZ.bed",
        minus="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_IZ.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        fa_plus=temp("resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_IZ_rmTTTT.fa"),
        fa_minus=temp("resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_IZ_rmTTTT.fa"),
        plus_rmTTTT="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_plus_IZ_rmTTTT.bed",
        minus_rmTTTT="resources/samples/DS/{samples}_{build}_sorted_ds_dipyrimidines_minus_IZ_rmTTTT.bed",
    log:
        "logs/{samples}/{samples}_{build}_filt4sim_ds_IZ_rmTTTT.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_filt4sim_ds_IZ_rmTTTT.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """
        (echo "`date -R`: Converting {input.plus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.plus} \
        -fo {output.fa_plus} \
        -s &&
        echo "`date -R`: Success! {input.plus} is converted." || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Converting {input.minus} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {input.minus} \
        -fo {output.fa_minus} \
        -s &&
        echo "`date -R`: Success! {input.minus} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (plus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_plus} \
        -o {output.plus_rmTTTT} \
        -r 'TTTT'\
        -m "DS" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Filtering by the given motif (minus)..." &&
        python3 workflow/scripts/rmMotifFromFasta.py \
        -i {output.fa_minus} \
        -o {output.minus_rmTTTT} \
        -r 'TTTT' \
        -m "DS" &&
        echo "`date -R`: Success! Filtering is done." ||
        {{ echo "`date -R`: Process failed..."; rm {output}; exit 1; }}  ) >> {log} 2>&1
        """

rule filt4sim_input_IZ:
    input:
        inpfile="resources/samples/input/{inputs}_{build}.fasta",
        region="results/regions/iz_hela_repdomains_uv_mean0.5_windows_201_100.bed",
        genome="resources/ref_genomes/{build}/genome_{build}.fa",
    output:
        inpfile_bed="resources/samples/input/{inputs}_{build}.bed",
        inpfile_IZ_bed="resources/samples/input/{inputs}_{build}_IZ.bed",
        inpfile_IZ="resources/samples/input/{inputs}_{build}_IZ.fasta",
    log:
        "logs/{inputs}/{inputs}_{build}_filt4sim_input_IZ.log",
    benchmark:
        "logs/{inputs}/{inputs}_{build}_filt4sim_input_IZ.benchmark.txt",
    conda:
        "../envs/bed2fasta.yaml"
    shell:
        """  
        (echo "`date -R`: Convert fasta to bed..." &&
        python3 workflow/scripts/fa2bed.py \
        -i {input.inpfile} \
        -o {output.inpfile_bed} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1       

        (echo "`date -R`: Getting input reads from IZ region..." &&
        bedtools intersect \
        -a {output.inpfile_bed} \
        -b {input.region} \
        -f 0.5 > {output.inpfile_IZ_bed} &&
        echo "`date -R`: Success!" || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1

        (echo "`date -R`: Converting {output.inpfile_IZ_bed} to fasta format..." &&
        bedtools getfasta \
        -fi {input.genome} \
        -bed {output.inpfile_IZ_bed} \
        -fo {output.inpfile_IZ} \
        -s &&
        echo "`date -R`: Success! {output.inpfile_IZ_bed} is converted." || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """