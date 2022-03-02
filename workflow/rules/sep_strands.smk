rule sep_strands_mutation:
    input:
        "results/mutation/{samples}/{samples}_target_mut.tsv",
    output:
        plus="results/mutation/{samples}/{samples}_target_mut_plus.tsv",
        minus="results/mutation/{samples}/{samples}_target_mut_minus.tsv",
    log:
        "logs/rule/analysis/{samples}/{samples}_sep_strands_mutation.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_sep_strands_mutation.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded mutations..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Mutations are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Separating minus stranded mutations..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Mutations are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """

rule sep_strands_edu:
    input:
        "results/edu/{samples}/{samples}_hg19_sorted_chr.bed",
    output:
        plus="results/edu/{samples}/{samples}_hg19_sorted_plus.bed",
        minus="results/edu/{samples}/{samples}_hg19_sorted_minus.bed",
    log:
        "logs/rule/analysis/{samples}/{samples}_hg19_sep_strands_edu.log",
    benchmark:
        "logs/rule/analysis/{samples}/{samples}_hg19_sep_strands_edu.benchmark.txt",
    shell:  
        """
        (echo "`date -R`: Separating plus stranded reads..." &&
        awk '{{if($6=="+"){{print}}}}' {input} > {output.plus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) > {log} 2>&1

        (echo "`date -R`: Separating minus stranded reads..." &&
        awk '{{if($6=="-"){{print}}}}' {input} > {output.minus} &&
        echo "`date -R`: Success! Reads are separated." || 
        {{ echo "`date -R`: Process failed..."; exit 1; }}  ) >> {log} 2>&1
        """