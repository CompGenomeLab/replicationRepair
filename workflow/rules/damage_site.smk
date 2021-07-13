rule damage_site_xr:
    input:
        plus="results/XR/{samples}/{samples}_{build}_xr_plus_sorted.txt",
        minus="results/XR/{samples}/{samples}_{build}_xr_minus_sorted.txt",
    output:
        plus="results/XR/{samples}/{samples}_{build}_sorted_xr_plus_damSite.txt",
        minus="results/XR/{samples}/{samples}_{build}_sorted_xr_minus_damSite.txt",
    log:
        "logs/{samples}/{samples}_{build}_damage_site_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_damage_site_xr.benchmark.txt",
    shell:
        """   
        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-8"\\t"$3-7"\\t""reads""\\t"".""\\t"$6}}' {input.plus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+7"\t"$2+8"\\t""reads""\\t"".""\\t"$6}}' {input.minus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule damage_site_ds:
    input:
        plus="results/DS/{samples}/{samples}_{build}_ds_dipyrimidines_plus_sorted.txt",
        minus="results/DS/{samples}/{samples}_{build}_ds_dipyrimidines_minus_sorted.txt", 
    output:
        plus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_plus_damSite.txt",
        minus="results/DS/{samples}/{samples}_{build}_sorted_ds_dipyrimidines_minus_damSite.txt",
    log:
        "logs/{samples}/{samples}_{build}_damage_site_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_damage_site_ds.benchmark.txt",
    shell:
        """   
        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-5"\t"$3-4"\\t""reads""\\t"".""\\t"$6}}' {input.plus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+4"\t"$2+5"\\t""reads""\\t"".""\\t"$6}}' {input.minus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule damage_site_xr_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_xr_sim_plus_sorted.txt",
        minus="results/sim/{samples}/{samples}_{build}_xr_sim_minus_sorted.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_xr_sim_plus_damSite.txt",
        minus="results/sim/{samples}/{samples}_{build}_xr_sim_minus_damSite.txt",
    log:
        "logs/{samples}/{samples}_{build}_damage_site_sim_xr.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_damage_site_sim_xr.benchmark.txt",
    shell:
        """   
        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-8"\\t"$3-7"\\t""reads""\\t"".""\\t"$6}}' {input.plus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+7"\t"$2+8"\\t""reads""\\t"".""\\t"$6}}' {input.minus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """

rule damage_site_ds_sim:
    input:
        plus="results/sim/{samples}/{samples}_{build}_ds_sim_plus_sorted.txt",
        minus="results/sim/{samples}/{samples}_{build}_ds_sim_minus_sorted.txt",
    output:
        plus="results/sim/{samples}/{samples}_{build}_ds_sim_plus_damSite.txt",
        minus="results/sim/{samples}/{samples}_{build}_ds_sim_minus_damSite.txt",
    log:
        "logs/{samples}/{samples}_{build}_damage_site_sim_ds.log",
    benchmark:
        "logs/{samples}/{samples}_{build}_damage_site_sim_ds.benchmark.txt",
    shell:
        """   
        (echo "`date -R`: Extracting exact damage sites (plus strand)..." &&
        awk '{{print $1"\\t"$3-5"\t"$3-4"\\t""reads""\\t"".""\\t"$6}}' {input.plus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.plus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") > {log} 2>&1

        (echo "`date -R`: Extracting exact damage sites (minus strand)..." &&
        awk '{{print $1"\\t"$2+4"\t"$2+5"\\t""reads""\\t"".""\\t"$6}}' {input.minus} |
        sort -u -k1,1 -k2,2n -k3,3n > {output.minus} &&
        echo "`date -R`: Success!" || 
        echo "`date -R`: Process failed...") >> {log} 2>&1
        """