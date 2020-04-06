
# Introduction

DNA replication is a highly conserved and regulated temporal process that is essential 
to the inheritance and the maintenance of the genome. 
In eukaryotes, replication fires from predefined locations called origin of replication (ORI) [8,9]
and once the firing occurs, 
replication fork forms to coordinate DNA replication [10]. 
Replication fork proceeds in both directions 
until it collides to another fork 
that is fired from an adjacent ORI.
Due to the fact that DNA polymerases can only elongate in a 5' to 3' manner,
the strand that is synthesized at the direction of the replication fork elongates continuously (leading strand),
while the opposite strand elongates in small pieces (lagging strand).
Moreover, replication of the strands are asymmetrical at opposite directions, 
so that forward strand acts as leading template (synthesizing lagging strand) at right replicating direction, 
and lagging template (synthesizing leading strand) at left replicating direction.
Variety of these replication forks, 
usually associated with high transcription, 
fire earlier than other forks 
which makes the process of replication a temporal process.

Stochastic effects of DNA replication are a major cause of mutagenesis, 
that contributes to cancer [12]. 
Therefore, an accurate and properly coordinated DNA replication is needed 
to prevent and to correct any errors for preserving the DNA fidelity
which is constantly threatened by both endogenous and exogenous sources during DNA replication.
Considering 70,000 lesions are occurred at a single cell per day [16],
it is assential to repair those lesions before the next cell division 
to avoid their permanent convertion into mutations. 
Genome-wide analysis of mutation rate suggest that 
mutation frequency increases at late replicating regions [17,18].
Moreover, late replicating regions at majority of cancers contain uplifted frequency of base substitution mutations 
than early replicating regions hold [19].
One possible explanation can be the time limitation of the lesions at late replicating regions 
for being recognized and corrected before the cell division ends. 
Furthermore, studies indicate that 
there is a strong correlation between the replication timing and mutational signatures.
In addition, a study revealed a significant replication strand asymmetry at most of the mutational signatures [20]. 
Previous studies demonstrated that there are strand specific mutation rate differences 
caused by the asymmetric synthesis of DNA [21,22,23,24].  

DNA excision repair mechanisms are known to relentlessly coup with different types of lesions. 
While oxidation, deamination and alkylation damages are repaired by Base Excision Repair (BER) [13], 
mismatches that escape proofreading are identified and corrected by Mismatch Repair (MMR) [14],
and bulky adducts caused by UV irradiation, environmental mutagens, and chemotherapeutic agents are removed by Nucleotide excision repair (NER) [15].
NER contains two subpathways called global genomic repair (GG-NER) and transcription coupled repair (TC-NER)
that differ from each other at the recognition of damage.
TC-NER is specialized into recognizing adducts at transcribed regions,
while GG-NER is able to recognize bulky adducts at any site. 
All three of these mechanisms are essential for the repair of ssDNA breaks 
that account for the 75% of daily occurred lesions.
In deficiencies of both MMR and NER mechanisms, 
there are specific mutational signatures associated 
which contribute to different cancer types [25]. 
NER associated signature 7 that represents UV induced mutations, 
display replication timing differences and replication related strand asymmetry,
suggesting that NER itself can be biased as well.
MMR suggested to create a mutation difference between early and late replicating regions 
by reducing mutations at early replicating regions more effectively 
due to the reachability of early replicating regions [26].
Similarly, TC-NER creates a transcriptional strand asymmetry 
by repairing adducts only at transcribed strand
and leaving the opposite strand untouched [27]. 
Even though signature 7 is linked with DNA replication timing and strand asymmetry [20], 
the contribution of GG-NER mechanism to mutation differences during replication
is yet to be known.  


In this paper, to study the influence of GG-NER mechanism on mutation distribution across the genome,
we analyzed the Damage-seq and eXcision Repair sequencing (XR-seq) data of 
cyclobutane pyrimidine dimers (CPDs) and pyrimidine-pyrimidone (6-4) photoproducts) [(6-4)PPs] 
from UV-irradiated HeLa cells synchronized at two stages of the cell cycle: 
early S phase, and late S phase. 
Damage-seq locates and quantifies the regions of UV induced CPD and (64)PP damages,
while XR-seq captures excised oligomers of the damage site 
that are removed by the NER mechanism.
With this approach, we are able to obtain the genome-wide distribution of UV-induced damages
and the differential repair frequency of these damage sites.
Our observations suggest that NER contributes to the mutational differences 
between early and late replicating regions due to its elevated activity at early replicating regions. 
In addition, repair efficiency of NER displays a replication related strand asymmetry near the origin of replication sites, 
favoring the lagging template strands, 
thus reducing the mutation rate of leading strands.

# Results

It is already known that the stratification of DNA replication and mutation rate is highly correlated, with increased mutation rate at late replicating domains [1,2]. 

We obtained replication origin sites from two independent datasets,
that are derived from Ok-seq data [5] which is sequenced locations of Okazaki fragments 
and SNS-seq data [6] which is the sequenced locations of small nascent DNA strands [7]. 
Then intersected these replication origin sites (Ok-seq and SNS-seq data separately) with replication domain data [4],
to differentiate early and late firing replication forks.  

To begin with, we wanted to validate the mutation asymmetry with our methods. 
Therefore, we mapped ICGC melanoma mutation dataset [3] to replication fork regions  
that are separated into their replication domains.   

Mutation counts of replication forks suggest that 
mutation burden of late replicating domains are much higher than it is in early replicating domains,
which is proven by multiple studies (referans ekle).
Moreover, mutation counts of transition zones are in between of ERDs and LRDs, 
and as we widen the region, it becomes more clear that
from the ERDs near side to LRDs near side, the mutation burden increases.  
The forkes that are derived from OK-seq demonstrates a much higher strand asymmetry,
while SNS-seq no or very little asymmetry at all replication domains. 


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