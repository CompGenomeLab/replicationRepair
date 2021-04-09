
source ../1_code/source_dir.sh

sample="" # name of the okseq data without the extension (must be .fastq.gz). 

okseqRaw="${mainPath}/0_data/OKseq"
okseqPath="${mainPath}/3_output/OKseq"
mkdir -p ${okseqPath} # directory for OKseq data

cutadapt -a ACACTCTTTCCCTACACGACGCTCTTCC -o ${okseqPath}/${sample}_cutadapt.fastq.gz ${okseqRaw}/${sample}.fastq.gz 

bowtie2 -p 4 -x ${genomePath}/Bowtie2/genome -U ${okseqPath}/${sample}_cutadapt.fastq.gz -S ${okseqPath}/${sample}_cutadapt.sam

samtools view -Sbq 20 -o ${okseqPath}/${sample}_cutadapt.bam ${okseqPath}/${sample}_cutadapt.sam 

samtools rmdup -s ${okseqPath}/${sample}_cutadapt.bam ${okseqPath}/${sample}_cutadapt_filt.bam

samtools sort ${okseqPath}/${sample}_cutadapt_filt.bam -o ${okseqPath}/${sample}_cutadapt_filt_sorted.bam 

samtools index ${okseqPath}/${sample}_cutadapt_filt_sorted.bam 