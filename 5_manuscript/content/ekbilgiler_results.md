## Figures

Briefly: Explain replication domains. Talk about method, data obtained. fig1


[46] - Our experimental results show that the number of damage-induced mutations reduces with the increasing time span between introduction of DNA damage and cell division. 



### images_ek/TS_NTS_repair_rate_log2norm

RPKM olarak da bak

64PP damages, and CPD 12 early do not a sign of trasncription bias.

CPD 12 late shows a little.

CPD async 12 and CPD 120 shows considerable amount of transcription bias. Why might be? Is TCR active?

### images_ek/fig2_new

Mutation and replication domain connection. Mutation strand asymmetry. Fig2

Bu plota update gelecek!! (Sadece dyprimidenler ile plotlanacak, ve trinuc ayrı şekilde de bakılacak.)

Belki initiation zoneların n sayısı eklenebilir

Generally we are seeing differences of mutation frequency on replication domains. ERDs are lower than LRDs and transition zones are in between. Because LRDs are replicated at the end of the replication process and until then it is mostly condensed, it is expected to be less repaired considering that nucleotide excision repair does not have the time to reach every single damage site. This is consistent with the known literature and visible at both initiation zone data (ok-seq and sns-seq).

There is a strand asymmetry on ok-seq data. For the ORIs corresponds to ERD and LRD the asymmetry is wider (ERD>LRD>UTZ=DTZ). Why exactly? Median length of domains might be affecting this? 

At right, plus strand is higher (have more mutation).

Also at the ok-seq data, the asymmetry is slightly slipping into right direction. Is it about ok-seq data? Both this fig, fig 3 and fig 4 have the same pattern. Check ok-seq for that.

SNS-seq is not showing an asymmetry. What might be the reason?

Close to the center, LRDs have higher strand asymmetry, but it is limited to a ~ 80 kb region. At a longer scale, ERD's strand asymmetry exceeds LRD. 

Very middle of the ERD has higher mutation. Why?

### images_ek/fig3

Strand asymmetry on damage and xr, possible reasons, maybe refer to supplementary. Fig3

first line of this figure matters. Others garbage!

hem simulated, hem de damage ve xr için replication domainler ayrıyken bak.


The asymmetry on the other hand is higher at the simulated ones. Why? What is causing the asymmetry, which is seen even higher at the fake reads?

At right, plus strand is higher (have more damage and repair).



### images_ek/fig4

Repair rate differences, domination of one strand over the other at the context of repair rate. Fig 4

200 kb alanda da  bi bakalım aynı plota

paired t-test gibi bişi uygulanabilir mi her windowdaki asymmetry için?

Generally, plot gives almost opposite of the mutation plot which makes sense. Meaning that NER is in fact a dominant factor to prevent mutations occurs on melanoma. 

ERDs are highest, LRDs are lowest. Because ERDs are repaired better, they end up having less mutations.

64PP damages do not show any strand asymmetry. Most of the damages are repaired probably.

CPD at 12 minutes has a small strand asymmetry, it gets higher at 120 minutes. So, when you wait longer asymmetry increases?

At right, minus strand is higher. Even though plus strand has more repair, because of the amount of damage it contains, minus strand is repaired more efficiently at plus strand.

Strand asymmetry increases from early to late.




### images_ek/chrom

Early ve late çizgilerini neye göre belirledin, bi mantığa oturt.

Berkayın yaptığı işten yürüyebilirsin

At faireW, low and transcription associated chromatin states, ERDs are repaired better at early, LRDs repaired better at late phase.

### 