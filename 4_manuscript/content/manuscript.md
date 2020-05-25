## Abstract 

Maintaining genome integrity is crucial for healthy cells to avoid cancer. 
Considering that DNA damages occur approximately 70,000 times per cell per day, 
repair of these damages is vital for the maintenance of genome stability [@tag:dnaDamage]. 
On the other hand, replication is the mechanism that causes unrepaired DNA damages to turn into mutations that might lead to cancer. 
The role of replication on DNA repair in general is yet to be clarified. 
Recently developed methods Damage-seq and XR-seq map damage formation and nucleotide excision repair events respectively, in various conditions [@tag:UVdamageRepair]. 
In this study, we analyze the Damage-seq and XR-seq results of cyclobutane pyrimidine dimers (CPDs) and pyrimidine-pyrimidone (6-4) photoproducts) ((6-4)PPs) 
from UV-irradiated HeLa cells synchronized at two stages of the cell cycle: early S phase, and late S phase. 
We compare these datasets with the localized replication origins and replication domains of HeLa cells [@tag:repDomain]. 
With the motivation to reveal how replication itself is influencing the repair preferences, 
we found out that in both early and late S phased cells, early replication domains are more efficiently repaired relative to late replication domains and 
repair efficiency favors leading strand around replication origins. 


## Introduction

DNA replication is a highly conserved and regulated temporal process that is essential 
to the inheritance and the maintenance of the genome. 
In fact, stochastic effects of DNA replication are a major cause of mutagenesis, 
that contributes to cancer [12]. 
Therefore, an accurate and properly coordinated DNA replication is needed 
to prevent and to correct any errors for preserving the DNA fidelity
which is constantly threatened by both endogenous and exogenous sources during DNA replication.
Considering 70,000 lesions are occurred at a single cell per day [16],
it is essential to repair those lesions before the next cell division, 
to avoid their permanent conversion into mutations. 

In about 8-10 hours, a mammalian DNA can replicates itself during the S phase of the cell cycle, 
thanks to many potential replication origins (ORIs) that are distributed throughout the genome.
As not each of these origins are fired at a cell division, 
they are not fired at the same time point as well.
However ORIs that are close to each other tend to fire in a synchronized manner,
resulting in simultaneous replication of Mb-sized regions called "replication domains" [28].
The regions that are replicated at the early S phase are called early replication domains (ERDs), 
which are usually structurally located at the interior of the nucleus,
while the nuclear periphery regions are replicated at late S phase, 
and therefore called late replication domains (LRDs) [29-35]. 
Genome-wide analysis of mutation rate suggest that 
mutation frequency increases at LRDs [17,18].
Moreover, LRDs at majority of cancers contain uplifted frequency of base substitution mutations 
than ERDs' hold [19].
One possible explanation can be the time limitation at LRDs 
for recognizing and correcting the lesions before the cell division ends.  

In eukaryotes, replication fires from predefined ORIs [8,9]
and once the firing occurs, 
replication fork forms to coordinate DNA replication [10]. 
Replication fork proceeds in both directions,
while polymerase epsilon continuesly synthesizes leading strand, 
polymerase delta discontinuesly synthesizes lagging strand until the fork collides to another 
that is fired from an adjacent ORI.
Due to the bidirectional movement of replication forks,
polymerases works asymmetrically on the strands according to the direction of the fork.
In other words, polymerase epsilon will be using plus strand as template at the left replicating fork,
minus strand at right replicating fork, and vice versa for polymerase delta. 
Expectedly, this asymmetric labor of polymerases during replication reflects the DNA fidelity as well.
A study revealed that 
there is a strong correlation between the replication timing and mutational signatures,
so that, a significant replication strand asymmetry is prominent 
for the majority of the mutational signatures [20].
In addition, it is demonstrated that 
there are strand specific mutation rate differences 
caused by the asymmetric synthesis of DNA [21,22,23,24].
More recent studies suggest that lagging strand is prone to mutations,
caused by the error-prone mechanism which bypassses the damage site [46].
Moreover, leading strand is more responsive to damage 
with the sensitive damage recognition of the attached helicase [47,48]. 
All these effects caused by the strand asymmetric labor of polymerases around ORIs
contributes to a mutational footprints at those regions [45]. 

DNA excision repair mechanisms are known to relentlessly coup with lesions 
that are potential sites of mutations. 
While oxidation, deamination and alkylation damages are repaired by Base Excision Repair (BER) [13], 
mismatches that escape proofreading are identified and corrected by Mismatch Repair (MMR) [14],
and bulky adducts caused by UV irradiation, environmental mutagens, and chemotherapeutic agents are removed by Nucleotide Excision Repair (NER) [15].
NER contains two subpathways called Global Repair (GR) and Transcription-Coupled Repair (TCR)
that differ from each other at the recognition of damage.
TCR is specialized into recognizing adducts at transcribed regions,
while GR is able to recognize bulky adducts at any site. 
All three of these mechanisms are essential for the repair of ssDNA breaks 
that account for the 75% of daily occurred lesions.
In deficiencies of both MMR and NER mechanisms, 
there are specific mutational signatures associated 
which contribute to different cancer types [25]. 
NER associated signature 7 that represents UV induced mutations, 
display replication timing differences and replication related strand asymmetry,
suggesting that NER itself can be biased as well.
MMR is suggested to create a mutation difference between ERDs and LRDs 
by reducing mutations at ERDs more effectively 
due to the reachability of ERDs [26].
Similarly, TCR creates a transcriptional strand asymmetry 
by repairing adducts only at transcribed strand
and leaving the opposite strand untouched [27]. 
Even though signature 7 is linked with DNA replication timing and strand asymmetry [20], 
the contribution of NER mechanism to mutation differences during replication
is still unclear.  

In this paper, to study the influence of NER mechanism on mutation distribution across the genome,
we analyzed the Damage-seq and eXcision Repair sequencing (XR-seq) data of 
cyclobutane pyrimidine dimers (CPDs) and pyrimidine-pyrimidone (6-4) pyrimidine-pyrimidone photoproducts ((6-4)PPs) 
from UV-irradiated HeLa cells synchronized at two stages of the cell cycle: 
early S phase, and late S phase. 
Damage-seq locates and quantifies the regions of UV induced CPD and (64)PP damages,
while XR-seq captures excised oligomers of the damage site 
that are removed by the NER mechanism.
With this approach, we are able to obtain the genome-wide distribution of UV-induced damages
and the differential repair frequency of these damage sites.
Our observations suggest that NER contributes to the mutational differences 
between ERDs and LRDs due to its elevated activity at ERDs. 
In addition, repair efficiency of NER displays a replication related strand asymmetry near the origin of replication sites, 
favoring the leading template strands, 
thus reducing the mutation rate of leading strands.

## Results

### ERDs are repaired more efficiently than LRDs, however the repair rate of LRDs elevate during the course of replication timing.

To determine how NER preferences influenced 
by replication domains during replication,
we compared repair efficiency of early replication domains and late replication domains
at two time points: early and late S phases of the cell cycle.
We obtained replication domains of HeLa cells from a study 
that developed a supervised method called DNN-HMM (Deep Neural Network-Hidden Markov Model) 
to define domains from repli-seq data [36].
Repair and damage samples are mapped to corresponding replication domains 
that split into 10 kb long bins covering total of 2 Mbp (200 bins) regions 
where replication domains are centered. 
Because damages might not be uniformly distributed throughout the genome,
we divided RPKM values of repair to the damage (XR-seq/Damage-seq) at each bin
to eliminate any such increase in repair caused by high damage frequency. 
By the approach, we are able to visualize the efficiency of repair at a given region, 
which we refer to as repair rate. 

Consistent across the samples, repair rate of ERDs are elevated at the center 
and gradually reduced when we move away from the center,
while LRDs exhibited an opposite pattern.
Considering that ERDs and LRDs strongly correlate with A/B compartments of Hi-C data respectively [37,38],
LRDs can be less reachable for NER complexes to identify damage sites.
Moreover, LRDs are known to contain higher mutation frequency than other regions,
hence; low repair rate of UV damages located at LRDs can be a key factor of mutagenesis
in the case of melanoma cancers.
On the other hand, the difference between early and late S phases suggest that 
repair rate head towards in favor of LRDs 
when replication timing continues early to late. 
Although we are observing a reduction of repair rate at ERDs and increase at LRDs,
XR-seq is a relative method. 
In other words, an increase at LRDs will affect genome-wide skew,  
cause aa attenuated view at the valley of ERDs 
even though the actual repair level might be the same as before. 
Therefore, the difference in repair rate might be indicating the reduced levels at ERDs,
increased levels at LRDs or both. 

Despite the findings being valid for all of the samples, the repair rate differs between them.
(6-4)PP damages at 12 minutes show very minor differences between phases, 
because of its almost instant repair after the damage occurrence [39].
In fact, (6-4)PP damages at 12 minutes are the only samples 
that do not show a significant increase of repair rate from early to late S phase (p=0.091).
CPD damages at 12 minutes and 2 hours demonstrate significant increase for LRDs 
and decrease for ERDs with p-values smaller than 0.00001.
Comparing 12 minutes to 2 hours, the accumulation of repair mechanism at LRDs lifted 
while at ERDs it is declined.
For global repair (GG-NER), CPD damages are much harder to dealt with. 
Therefore, the repair efficiency of CPD damages are not high enough 
to make an influence on damage levels at 1 hour after exposure [39].
However CPD damages at 2 hours have a closer pattern to (6-4)PP damages at 12 minutes,
suggesting that repair efficiency of CPD damages increases over time.      

### Variety of chromatin states are associated with dominant repair. 

Active chromatin states known to be repaired better, 
basically because those regions are more accessible for NER [40].
However, we wonder how repair machinery at ERDs and LRDs is influenced 
by the chromatin states during the replication.
We retrieved chromatin states of HeLa cells segmented by ChromHMM from UCSC website.
Then, intersected the chromatin states data to replication domains 
and mapped our samples to these regions.
Later on, we aggregated the data so that for each chromosome 
there are 25 chromatin states at ERDs and 25 chromatin states at LRDs.
After calculating the repair rates, we further divided early/late phases
to observe at which time point those regions repaired with higher efficiency.

!!!Will be added!!!

### Asymmetric damage around initiation zones causes asymmetric repair profiles.

Even so replication domains can give a general idea 
about the replication affect on NER,
these domains are 1-2 Mb-sized chunks of DNA,
which represents the accumulated involvement of many small ORIs they harbor. 
To better understand the effect of replication at a smaller scale,
we decided to check the damage and repair profiles 
at the sites of ORIs where the replication process starts.
To achieve our goal, we retrieved two publicly available, independent dataset 
that holds the information of ORIs.
One is derived from OK-seq data which sequence highly purified Okazaki fragments
to quantify initiation zones 
that are the sets of closely positioned replication origin sites [5]. 
The other one is derived from short nascent strand sequencing (SNS-seq) [6], 
which precisely identifies the ORIs [7].
Remarkably, we observed explicit strand asymmetry at both directions around the initiation zones
that is consistent across the samples. 
The asymmetry suggests that lagging template strand (minus strand at left direction, plus strand at right direction) harbors more damage sites and attracts more repair accordingly. 
Due to suspecting that the asymmetry might be caused by a strand bias,
we simulated our damage and repair signals using their nucleotide frequencies, 
and compared the simulated reads with real ones to see if they will be correlate.
Our simulated reads indeed have the same asymmetry pattern, 
with lower RPKM values in general, but with an elevated asymmetry.
It is not surprising that simulated reads for all samples are lower, 
because they are not the real damage or repair signals, 
though they are similar in nucleotide composition.    
On the other hand, the elevated strand asymmetry compared to real reads is somewhat interesting,
suggesting that NER might not create, 
oppositely, might be reducing the strand asymmetry at initiation zones.
MMR have shown to have a balancing effect on mutational asymmetries during replication [41].
In addition, for both simulated and real reads, the flip between strands is not centered,
sliding slightly to the right direction.
The reason is, initiation zones contains multiple potential ORIs
that is not known which one will be fired in a given replication process. 
We centered the initiation zones for our analysis, however as it seems, 
initiated ORIs mostly accumulated at 1-2 kb rightward of the center in our experiments. 


### Mutation profile on replication domains

Any asymmetry at NER should directly impact the mutation frequency at melanoma cancers.
With this reasoning, we obtained melanoma mutation dataset 
from International Cancer Genome Consortium (ICGC) data portal 
to see if mutation frequency on initiation zones correlate with the repair and damage signals.
After filtering (see Methods), we mapped the mutation data to 20 kb regions centered at initiation zones 
that are separated into their corresponding replication domains. 
By the approach, we are able to compare both the effect of replication domains 
and initiation zones within them. 
We also extended the analyzed region to 200 kb 
to observe the mutation profile around initiation zones on a wider scale.

In the context of replication domains, the mutation count of ERDs and LRDs are exact opposite of 
what we have seen at NER rate of these domains. 
This makes perfect sense, 
simply because high repair rates at ERDs will lead to less mutation in return, vice versa for LRDs.
Moreover, LRDs are replicated at the end of the replication process 
and until then it is mostly condensed, 
it is expected to be less repaired considering 
that NER does not have the time to coup with the damage levels.
Thus, our results regarding ERDs and LRDs very much correlates with the mutation data 
and in agreement with previous studies [17-19].
Mutation counts of transition zones are in between ERDs and LRDs.
Furthermore, in the wider scale mutation counts of transition zones increases 
at the regions that are distant to ERDs. 
In addition, mutation counts reveal a strand asymmetry at all the initiation zones 
no matter which replication domain they belong to. 
While LRDs show the sharpest difference between strands, 
the asymmetry at ERDs is wider than it is at LRDs.
One possible reason can be the amount of ORIs they harbor.
Earlier studies suggest that ERDs contain significantly higher replication origin [42], 
and the cumulative effect of these origins might create a wider asymmetry.
Additionally, replication fork movement at LRDs (1.5–2.3 kb/min) are known to be faster than
ERDs (1.1–1.2 kb/min) [43],
which might cause more mutation and increased asymmetry between strands.

The asymmetry of strands suggest that lagging strand (minus strand at left direction, plus strand at right direction) have more mutations than leading.
It might seem contradictory with the previous figure which we indicated that
lagging templates have more repair.
But more repair does not necessarily mean that it is repaired better.
In fact, when we check the repair rate (XR-seq/Damage-seq) of our samples,
we observe that although lagging strand template repaired more, 
the leading strand template is repaired more efficiently than lagging.   


### Repair rate of samples

To see whether the repair rate of our samples are correlating with the mutation data,
as we did with mutation data, we mapped our samples to 20 kb regions centered at initiation zones 
that are separated into their corresponding replication domains.
Expectedly, our samples are inversely correlated with the mutation counts.
Strand asymmetry appears at all samples with CPD damages
in the favor of leading template strands, 
which supports less mutation count on the leading strand.
Samples with (6-4)PP damages are not showing any strand difference, 
because of NER mechanism's fast repair ability at these damage site,
further suggesting that NER might be a key factor 
on balancing the strand asymmetry.

!!!!Difference between CPD 12 and 120?? 

In addition, repair rate differences of replication domains are exactly opposite of mutation counts at those regions which implies that efficient excision repair by the NER machinery is a major effector on preventing melanoma cancer.


## Methods 

### Cell culture and treatments

HeLa S3 cells (purchased from ATCC) were cultured in DMEM medium supplemented with 10%
FBS and 1% penicillin/streptomycin at 37 °C in a 5% atmosphere CO 2 humidified chamber.
Cells were synchronized at late G1 phase by double-thymidine treatment, then released into S
phase by removing thymidine. In brief, for first thymidine treatment, thymidine was added to
cells at 50% confluence to a final concentration of 2mM. After 18 hours, cells were released by
washing with PBS and cultured in fresh medium for 9 hours. Then cells were treated with 2mM
thymidine for 15 hours and released into S phase for indicated time before UV irradiation. Cells
were irradiated with 20J/m 2 of UVC as described, then collected immediately or after incubation
at 37 °C for indicated time for the following assays.

### Flow cytometry analysis

Cells were trypsinized, PBS washed and fixed in 70 % (v/v) ethanol at −20°C for at least 2 hours,
then stained in the staining solution (0.1% (v/v) Triton X-100, 0.2 mg/ml RNase A and 20 g/ml
propidium iodide in PBS) for 30 min at room temperature. Then S phase progression of the cells
was analyzed by a flow cytometer (Beckman Coulter CyAn or Becton Dickinson FACSCanto II).

### Damage-seq and XR-seq libraries preparation and sequencing

Cells were collected in ice-cold PBS at indicated time and subjected to Damage-seq (HS-
Damage-seq) and XR-seq as previously described.
For Damage-seq, briefly, genomic DNA were extracted with PureLink Genomic DNA Mini Kit
(Thermo) and sheared by sonication with a Q800 Sonicator (Qsonica). DNA fragments (1μg)
were used for end repair, dA-tailing and ligation with Ad1 by NEBNext UltraII DNA Library
Prep kit (New England Biolabs). Samples were then denatured and subjected to
immunoprecipitation with anti-(6-4)PP or anti-CPD antibody (Cosmo Bio), respectively. In order
to detect the precise position of lesions, primer Bio3U was attached to purified DNA and
extended by NEBNext UltraII Q5 DNA polymerase (New England Biolabs). After purification,
the extension products were annealed to oligo SH for subtractive hybridization. Oilgo SH was
then removed by incubating with Dynabeads MyOne Streptavidin C1 (Thermo) and the rest of
the sample was ligated to Ad2, followed by PCR amplification.
For XR-seq, in brief, fresh cells were lysed by Douce Homogenizer in Buffer A (ref). After
centrifuging to remove chromatin DNA, co-immunoprecipitation with anti-XPG (Santa Cruz)
antibody was performed to pull down primary excision products generated by nucleotide excision
repair. Purified excised fragments were ligated to both 5’ and 3’ adaptors. Ligation products were
further purified by immunoprecipitation with either anti-(6-4)PP or anti-CPD antibody and
repaired by corresponding photolyase, followed by PCR amplification and gel purification.
Libraries with different indexes were pooled and sequenced in SE50 form on Hiseq 2000/2500
platform by the University of North Carolina High-Throughput Sequencing Facility, or in PE150
form on Hiseq X platform by the WuXiNextCODE Company.

### Damage-seq sequence pre-analysis

The sequenced reads with adapter sequence GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT 
at 5' end, 
were discarded via cutadapt with default parameters 
for both single-end and paired-end reads [@tag:cutadapt]. 
The remaining reads were aligned to the hg19 human genome using bowtie2 with 4 threads (`-p`) [@tag:bowtie2].
For paired-end reads, maximum fragment length (`-X`), which means the maximum accepted total length of mated reads and the gap between them, was chosen as 1000. 
Using samtools, aligned paired-end reads were converted to bam format, 
sorted using `samtools sort -n` command, 
and properly mapped reads with a mapping quality greater than 20 were filtered 
using the command `samtools view -q 20 -bf 0x2` in the respective order [referans]. 
Then, resulting bam files were converted into bed format using 
`bedtools bamtobed -bedpe -mate1` command [@tag:bedtools].
The aligned single-end reads were directly converted into bam format 
after the removal of low quality reads (mapping quality smaller than 20) 
and further converted into bed format with `bedtools bamtobed` command [@tag:bedtools].
Because the exact damage sites should to be positioned at two nucleotides upstream of the reads [referans], bedtools flank and slop command were used to obtain 10 nucleotide long positions 
bearing damage sites at the center (5. and 6. positions) [@tag:bedtools].
The reads that have the same starting and ending positions, were reduced to a single read 
for deduplication and remaining reads were sorted with the command 
`sort -u -k1,1 -k2,2n -k3,3n`. 
Then, reads that did not contain dipyrimidines (TT, TC, CT, CC) at their damage site (5. and 6. positions) were filtered out to eliminate all the reads that do not harbor a UV damage.
Lastly, only the reads that were aligned to common chromosomes (chromosome 1-22 + X)
were held for further analysis.


### XR-seq sequence pre-analysis

The adaptor sequence TGGAATTCTCGGGTGCCAAGGAACTCCAGTNNNNNNACGATCTCGTATGCCGTCTTCTGCTTG at the 3' of the reads were trimmed and sequences without the adaptor sequences were discarded using cutadapt with default parameters [@tag:cutadapt].
Bowtie2 was used with 4 threads (`-p`) to align the reads to the hg19 human genome [@tag:bowtie2]. 
Then reads with mapping quality smaller than 20 were removed by samtools [referans]. 
Bam files obtained from samtools were converted into bed format by bedtools [@tag:bedtools].
Multiple reads that were aligned to the same position, were reduced to a single read to 
prevent duplication affect and remaining reads were sorted with the command 
`sort -u -k1,1 -k2,2n -k3,3n`.  
Lastly, only the reads that were aligned to common chromosomes 
were held for further analysis.


### Dna-seq sequence pre-analysis

Paired-end reads were aligned to hg19 human genome via bowtie2 with 4 threads (`-p`) 
and maximum fragment length (`-X`) chosen as 1000 [@tag:bowtie2]. 
Sam files were converted into bed format 
as it was performed at damage-seq paired-end reads.
Duplicates were removed 
and reads were sorted with `sort -u -k1,1 -k2,2n -k3,3n` command.
Lastly, the reads that did not align to the common chromosomes were discarded.


### XR-seq and Damage-seq simulation

Art simulator was used to produce synthetic reads with the parameters
`-l 26 -f 2`, `-l 10 -f 2` for XR-seq and Damage-seq, respectively [44].
To better represent our filtered real reads,
read length (`-l`) parameter was chosen as the most frequent read length after pre-analysis done.
The `.fastq` file that Art produced, were filtered according to our reads by 
calculating a score using nucleotide frequency of the real reads and
obtaining most similar 10 million simulated reads. 
The filtering was done by `filter_syn_fasta.go` script, 
which is available at the repository: https://github.com/compGenomeLab/lemurRepair.
Filtered files were preceded by pre-analysis again for further analysis.


### Data sets

Replication domains obtained from Gene Expression Omnibus (GEO), accession no: GSE53984. 
High replication fork direction and initiation zones datasets 
which kindly given us can be obtained from accession no: SRP065949 as raw OK-seq data. 
Lastly SNS-seq data can be obtained from Gene Expression Omnibus (GEO), accession no: GSE37757. 

Melanoma somatic mutations of 183 tumor samples were obtained 
from the data portal of International Cancer Genome
Consortium (ICGC) as compressed `.tsv` files which is publicly available at 
https://dcc.icgc.org/releases/release_28/Projects/MELA-AU.
Single base substitution mutations were extracted, 
and only the mutations of common chromosomes were used.
To obtain the mutations that could be caused by UV-induced photoproducts, 
C -> T mutations that have a pyrimidine at the upstream nucleotide 
was further extracted.  
Later on, mutations were quantified on 20 kb long initiation zones 
that were separated into their corresponding replication domains 
using `bedtools intersect` command with the `-wa -c -F 0.5` options.


### Further Analysis

In order to separate a region data (replication domains, initiation zones, or replication origins)
into chosen number of (201) bins,
the start and end positions of all the regions set to a desired range with the unix command:
`awk -v a="$intervalLen" -v b="$windowNum" -v c="$name" '{print $1"\t"int(($2+$3)/2-a/2-a*(b-1)/2)"\t"int(($2+$3)/2+a/2+a*(b-1)/2)"\t"$4"\t"".""\t"$6}'`
Then, any intersecting regions or regions crossing the borders of its chromosomes were filtered  
to eliminate the possibility of signal's canceling out effect.
After that, `bedtools makewindows` command was used with the `-n 201 -i srcwinnum` options 
to create a `.bed` file contaiing the bins.

To quantify the XR-seq and Damage-seq profiles on the prepared `.bed` file,
`bedtools intersect` command was used to intersect  
as it was performed for mutation data.
Then all bins were aggregated given their bin numbers,
and the mean of the total value of each bin were calculated.
Lastly RPKM normalization was performed and the plots were produced  
using ggplot2 in R programming language.



# References
1. https://www.nature.com/articles/ng.363.pdf
2. https://www.sciencedirect.com/science/article/pii/S0002929712005770
3. https://dcc.icgc.org/releases/release_28/Projects/MELA-AU
4. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34399
5. https://www.ncbi.nlm.nih.gov/pubmed/26751768
6. https://www.nature.com/articles/nsmb.2339#accessions
7. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5137433/
8. https://www.ncbi.nlm.nih.gov/pubmed/27147572
9. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5342515/
10. https://www.ncbi.nlm.nih.gov/pubmed/19652539
11. https://www.ncbi.nlm.nih.gov/pubmed/23603017
12. https://www.ncbi.nlm.nih.gov/pubmed/25554788
13. https://cshperspectives.cshlp.org/content/5/4/a012583.short
14. https://www.annualreviews.org/doi/abs/10.1146/annurev.biochem.74.082803.133243
15. https://cshperspectives.cshlp.org/content/5/10/a012609.short
16. https://www.ncbi.nlm.nih.gov/pubmed/12760027
17. https://www.nature.com/articles/nature12213
18. https://www.nature.com/articles/ng.363
19. https://www.nature.com/articles/nature11273
20. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1509-y
21. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3469411/ 
22. https://www.sciencedirect.com/science/article/pii/S0092867415017146?via%3Dihub 
23. https://www.nature.com/articles/nature14183 
24. https://www.ncbi.nlm.nih.gov/pubmed/25228659
25. https://www.ncbi.nlm.nih.gov/pubmed/24981601
26. https://www.nature.com/articles/nature14173
27. https://www.ncbi.nlm.nih.gov/pubmed/25456125
28. https://www.ncbi.nlm.nih.gov/pubmed/9508763
29. https://www.ncbi.nlm.nih.gov/pubmed/12356909
30. https://www.ncbi.nlm.nih.gov/pubmed/2910875
31. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2289349/
32. https://www.ncbi.nlm.nih.gov/pubmed/25416942
33. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2556267/
34. https://www.ncbi.nlm.nih.gov/pubmed/18842067
35. https://www.pnas.org/content/107/1/139
36. https://www.ncbi.nlm.nih.gov/pubmed/26545821
37. https://www.ncbi.nlm.nih.gov/pubmed/20430782
38. https://www.nature.com/articles/nature13986
39. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5495279/
40. https://www.ncbi.nlm.nih.gov/pubmed/27036006
41. https://www.cell.com/abstract/S0092-8674(15)01714-6
42. https://www.ncbi.nlm.nih.gov/pubmed/22751019
43. https://www.sciencedirect.com/science/article/pii/S001448270400655X?casa_token=oZ-XePmnVusAAAAA:BAgEC0C-UqTgqVJGIglRWOhhIv1qKXF-OKQezUr07Q-sSAbQbToWt-1ipgqjsFDhdz9bFeUUkpY
44. https://www.ncbi.nlm.nih.gov/pubmed/22199392
45. https://www.sciencedirect.com/science/article/pii/S0962892416300381?casa_token=9wTh7tOmy1QAAAAA:ixkZ8YuO8pMl-NlA-EgxH2iDmK3ZNL-HYSOxPkvvGSYAQMMuk2qi8dh_6m9Thn_lJ6J8-mXIMdI
46. https://www.nature.com/articles/s41588-018-0285-7
47. https://cshperspectives.cshlp.org/content/5/5/a012815.short
48. https://pubs.acs.org/doi/abs/10.1021/acs.chemrev.7b00046?casa_token=GGyTMY1W4N0AAAAA:T8-MbcOg8Coem5IPrmjt84dch6xXi5sM39Jrdsl4LjPlMiz7KWuZub489ciO_WbgiIZkuNgITgt_jCJZ#
