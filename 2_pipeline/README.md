## XRseq and Damage-seq Pipeline Steps  

  

* XRseq and Damage-seq pipelines are executed in two different shell scripts. In order to use this pipeline, the name of the sample name must be provided together with the executed script and some other information unique to the sample that will be important for the analysis. If sample name is sample_name.fastq the script can be executed as  


```  
SAMPLE=sample_name  

~/exonRepair_DS.sh $SAMPLE  

~/exonRepair_XR.sh $SAMPLE
``` 


* Paths of the codes, samples, and genome are stored in a separate file which are loaded into script with `source` built-in command when the script is executed:

```  
source ../1_code/source_dir.sh 
``` 

* The first step will be the adaptor trimming using Cutadapt. For XRseq, 3' adaptors will be trimmed as

```
cutadapt -a ${ADAPTOR} -o ${SAMPLE}_cutadapt.fastq ${SAMPLE}.fastq
```
* If the layout is paired:

```
cutadapt -a ${ADAPTOR} -A ${ADAPTOR} -o ${SAMPLE}_1_cutadapt.fastq -p ${SAMPLE}_2_cutadapt.fastq ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq 
```

* For Damage-seq, the reads that contains adaptors at 5' end will be discarded as

```
cutadapt --discard-trimmed -g ${ADAPTOR} -o ${SAMPLE}_cutadapt.fastq ${SAMPLE}.fastq 
```

* Again, if the layout is paired:

```
cutadapt --discard-trimmed -g ${ADAPTOR} -G ${ADAPTOR} -o ${SAMPLE}_1_cutadapt.fastq -p ${SAMPLE}_2_cutadapt.fastq ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq
```

* For more information about the trimming process, you can check the manual of [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html). After the adaptor trimming we will align our reads to the genome by Bowtie2. The command used for the alignment is the same for both XRseq and Damage-seq pipelines. If the layout is single:

```
bowtie2 -p 4 -x ${GENOME} -U ${SAMPLE}_cutadapt.fastq -S ${SAMPLE}_cutadapt.sam
```
* And for the paired end reads, the command is 

```
bowtie2 -p 4 -X 1000 -x ${GENOME} -1 ${SAMPLE}_1_cutadapt.fastq -2 ${SAMPLE}_2_cutadapt.fastq -S ${SAMPLE}_cutadapt.sam
```
* For more about Bowtie2 and the alignment process, you can check the manual of [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

* Sam files can be converted into bam files via SAMtools while filtering out quality scores lower than 20. SAMtools includes many features to visualize and manipulating your reads including sorting, merging, indexing and many more. You can check its manual to find out more about [SAMtools](http://www.htslib.org/doc/samtools-1.2.html). Later on we are using BEDTools to convert bam files to bed. For single-end layout:

```
samtools view -q 20 -b -o ${SAMPLE}_cutadapt.bam ${SAMPLE}_cutadapt.sam

bedtools bamtobed -i ${SAMPLE}_cutadapt.bam > ${SAMPLE}_cutadapt.bed
```

* For paired-end:

```
samtools view -Sb -o ${SAMPLE}_cutadapt.bam ${SAMPLE}_cutadapt.sam # samtools: sam to bam

samtools sort -n ${SAMPLE}_cutadapt.bam -o ${SAMPLE}_cutadapt_sorted.bam

samtools view -q 20 -bf 0x2 ${SAMPLE}_cutadapt_sorted.bam | bedtools bamtobed -bedpe -mate1 > ${SAMPLE}_cutadapt.bedpe

awk '{if($9=="+"){print $1"\t"$2"\t"$6"\t"$7"\t"$8"\t"$9}}' ${SAMPLE}_cutadapt.bedpe > ${SAMPLE}_cutadapt_plus.bed

awk '{if($9=="-"){print $1"\t"$5"\t"$3"\t"$7"\t"$8"\t"$9}}' ${SAMPLE}_cutadapt.bedpe > ${SAMPLE}_cutadapt_minus.bed

cat ${SAMPLE}_cutadapt_plus.bed ${SAMPLE}_cutadapt_minus.bed > ${SAMPLE}_cutadapt.bed
```

* Important Note: After obtaining the bed files, all the process coming on will be same for both layouts.

* We are sorting the coordinates and removing the duplicates. The removal is needed to avoid having the same read multiple times. Then we are checking the number of mapped reads:

``` 
sort -u -k1,1 -k2,2n -k3,3n ${SAMPLE}_cutadapt.bed > ${SAMPLE}_cutadapt_sorted.bed

grep -c "^" ${SAMPLE}_cutadapt_sorted.bed > ${SAMPLE}_cutadapt_sorted_readCount.txt
```

* After this point, the commands of XRseq and Damage-seq are differ a little bit. However, the reasoning behind processes are still similar to each other. For XRseq, we are checking the length distribution of the reads. The experimental design of XRseq includes the immunoprecipitation of excised oligomer, which is cut off from the genome by the Nucleotide Excision Repair Mechanism (NER). Therefore, the length of the oligomer can vary. The damage that responsible for the NER activation, usually positioned at 6 to 7 nucleotide before the last base of the oligomer, which means the damage position vary between the oligomers of different length. 

```
awk '{print $3-$2}' ${SAMPLE}_cutadapt_sorted.bed | sort -k1,1n | uniq -c | sed 's/\s\s*/ /g' | awk '{print $2"\t"$1}' > ${SAMPLE}_cutadapt_length_distribution.txt
```       

* If the most occurred read length is 26, our next step will be checking the dinucleotide content of the reads with 26 base pairs, to detect the dinucleotide distribution at the site of damage. 	

```
awk '{ if ($3-$2 == 26) { print } }' ${SAMPLE}_cutadapt_sorted.bed > ${SAMPLE}_cutadapt_sorted_26.bed

bedtools getfasta -fi ${GENOME} -bed ${SAMPLE}_cutadapt_sorted_26.bed -fo ${SAMPLE}_cutadapt_sorted_26.fa

fa2kmerAbundanceTable.py -i ${SAMPLE}_cutadapt_sorted_26.fa -k 2 -o ${SAMPLE}_cutadapt_sorted_26_dinucleotideTable.txt
```

* `fa2kmerAbundanceTable.py` is a python script that created by Ogun Adebali to obtain the counts of dinucleotides (or a motif in any length) at each position. You can click [here](https://github.com/adebali/NGStoolkit/blob/master/bin/README.md) for more information. 

* Furthermore, we are filtering out the chromosomes with invalid names and separating our reads to plus and minus strand. We also create a variable which will store the count of mapped reads for later use. 

```
grep "chr" ${SAMPLE}_cutadapt_sorted.bed | grep -v -e "chrY" -e "chrM" > ${SAMPLE}_cutadapt_sorted_chr.bed

awk '{if($6=="+"){print}}' ${SAMPLE}_cutadapt_sorted_chr.bed > ${SAMPLE}_cutadapt_sorted_plus.bed

awk '{if($6=="-"){print}}' ${SAMPLE}_cutadapt_sorted_chr.bed > ${SAMPLE}_cutadapt_sorted_minus.bed

mappedReads="$(grep -c '^' ${SAMPLE}_cutadapt_sorted_chr.bed)"
```

* In Damage-seq, we already know that the length of the reads are 50 base pairs because of the experimental design. However, the methodology of Damage-seq based on not to obtaining the damage itself but the site that just before that damage. So, in order to get the damage site we use BEDTools features called `flank` and `slop`. Basically with the `flank` option, you can get the flanking sites of your read and with the `slop` option, you can increase or decrease the boundaries of your reads. For more information about many features of BEDTools, you can check the [manual](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html). The critical part here is that the position of damage in the reads at the opposite strand will be assymetric. Therefore we must separate our reads to plus and minus strand and then continue with these steps.

```
bedtools flank -i ${SAMPLE}_cutadapt_sorted_plus.bed -g ${GENOME} -l 6 -r 0 > ${SAMPLE}_cutadapt_flanked_plus.bed

bedtools flank -i ${SAMPLE}_cutadapt_sorted_minus.bed -g ${GENOME} -l 0 -r 6 > ${SAMPLE}_cutadapt_flanked_minus.bed

bedtools slop -i ${SAMPLE}_cutadapt_flanked_plus.bed -g ${GENOME} -l 0 -r 4 > ${SAMPLE}_cutadapt_slopped_plus.bed

bedtools slop -i ${SAMPLE}_cutadapt_flanked_minus.bed -g ${GENOME} -l 4 -r 0 > ${SAMPLE}_cutadapt_slopped_minus.bed
```

* After executing these commands, our reads become 10 base pair long containing damage site at the middle positions (5-6). Then we can check the dinucleotide distribution of plus and minus strand files separately. In Damage-seq method, because we know the exact site of the damage, we can further filter our reads based on the dinucleotides on the damage site. For example, we prepared a Damage-seq experiment to obtain cyclobutane pyrimidine dimers (CPDs) which can occur on TT, TC, CT, CC dinucleotides. All the other dinucleotides located at the damage site of our reads can be considered as noise and can be discarded. 

```
fa2bedByChoosingReadMotifs.py -i ${SAMPLE}_cutadapt_sorted_plus_10.fa -o ${SAMPLE}_cutadapt_sorted_plus_dipyrimidines.bed -r ".{4}(c|t|C|T){2}.{4}"

fa2bedByChoosingReadMotifs.py -i ${SAMPLE}_cutadapt_sorted_minus_10.fa -o ${SAMPLE}_cutadapt_sorted_minus_dipyrimidines.bed -r ".{4}(c|t|C|T){2}.{4}"
```

* `fa2bedByChoosingReadMotifs.py` is a python script for obtaining the the reads with the chosen motifs which is written by Ogun Adebali. The source code of the script is available at the [link](https://github.com/adebali/NGStoolkit/blob/master/bin/fa2bedByChoosingReadMotifs.py). Then, as we did at XRseq, we store the mapped read count in a variable for later use.

```
minus_line="$(grep -c '^' ${SAMPLE}_cutadapt_sorted_minus_dipyrimidines.bed)"
plus_line="$(grep -c '^' ${SAMPLE}_cutadapt_sorted_plus_dipyrimidines.bed)"
mappedReads=`echo "$minus_line + $plus_line" | bc`  
```

* Important Note: After this step, all the process coming on will be same for both XRseq and Damage-seq.  

* We are intersecting our reads with bed files containing target regions and counting the number of overlap.

```
bedtools intersect -a ${DATA} -b ${SAMPLE}_cutadapt_sorted_plus.bed -wa -c -F 0.5 >  ${SAMPLE}_cutadapt_sorted_plus_intersected.bed

bedtools intersect -a ${DATA} -b ${SAMPLE}_cutadapt_sorted_minus.bed -wa -c -F 0.5 > ${SAMPLE}_cutadapt_sorted_minus_intersected.bed
```

* Some bed files with the target regions are separated into multiple windows and each of these windows of a region tagged with a window number. If we want to aggregate all the windows of the same number, `combineWindows.py` command can used for the task:

```
combinewindows.py -i ${SAMPLE}_cutadapt_sorted_plus_intersected.txt ${combine_options[i]} -o ${SAMPLE}_cutadapt_sorted_plus_intersected_combined.txt

combinewindows.py -i ${SAMPLE}_cutadapt_sorted_minus_intersected.txt ${combine_options[i]} -o ${SAMPLE}_cutadapt_sorted_minus_intersected_combined.txt
```

* Also we add some information about the sample by retrieving it from a csv file we prepared before, called `samples.csv`. 

```
moreinfo="$(grep $1 samples.csv | sed 's/,/\t/g' | awk '{print $1}')"

addColumns.py -i ${SAMPLE}_cutadapt_sorted_plus_intersected.bed -o ${SAMPLE}_cutadapt_sorted_plus_intersected_info.txt -c + "${moreinfo}" "${mappedReads}"

addColumns.py -i ${SAMPLE}_cutadapt_sorted_minus_intersected.bed -o ${SAMPLE}_cutadapt_sorted_minus_intersected_info.txt -c - "${moreinfo}" "${mappedReads}"
```

* `addColumns.py` is a python script which is written by Ogun Adebali. The source code of the script is available at the [link](https://github.com/adebali/NGStoolkit/blob/master/bin/addColumns.py).

* Then, we are calculating the RPKM values of these regions:

```
python RPKM.py -i ${SAMPLE}_cutadapt_sorted_plus_intersected_info.txt -chse 2 3 -c 7 -mr 0 -o ${SAMPLE}_cutadapt_sorted_plus_intersected_rpkm.txt

python RPKM.py -i ${SAMPLE}_cutadapt_sorted_minus_intersected_info.txt -chse 2 3 -c 7 -mr 0 -o ${SAMPLE}_cutadapt_sorted_minus_intersected_rpkm.txt
```
* Lastly, we combine the rpkm files of all samples in a master report file using `cat` command for further processing.
  
