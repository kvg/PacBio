A first look at the PacBio 3D7 data
===================================







We prepared roughly 15ug of DNA from the canonical Plasmodium falciparum strain, 3D7, for sequencing on Cold Spring Harbor Laboratory's PacBio RS-II instrument.  This was sequenced on eight SMRT cells, which should yield approximately 100x coverage over the 23 megabase genome.  Having just received the filtered subreads (in unaligned and aligned form), we investigated some very basic properties of the data.

Read length
===========

|key    |  numReads|  minLength|  maxLength|  meanLength|  n50Value|
|:------|---------:|----------:|----------:|-----------:|---------:|
|pb0    |     67425|         50|      35806|        6197|      9225|
|pb1    |     63527|         50|      37335|        6593|      9707|
|pb2    |     70147|         50|      40116|        6174|      9320|
|pb3    |     76364|         50|      35377|        6379|      9277|
|pb4    |     49238|         50|      33801|        5700|      8719|
|pb5    |     36364|         50|      35002|        5503|      8490|
|pb6    |     33385|         50|      34139|        5461|      8389|
|pb7    |     30082|         50|      35126|        5071|      8055|
|total  |    426532|         50|      40116|        6031|      9058|


We first examined read count and read length properties, as shown in the table above.  We have labeled data from each SMRT cell with a simple identifier, pb0 .. pb7 (the full identifier is a bit long and cumbersome to use everywhere). On average, each SMRT cell yields 5.3316 &times; 10<sup>4</sup> reads, typically 5884.704 +/- 527.5589 bases long.

We next examined the read length distribution, shown below.

![plot of chunk lengthHist](figure/lengthHist.png) 


Interestingly, this distribution appears to be bimodal, with peaks at 1.1 kb and 8.0 kb.  The origin of these peaks are unclear.  Furthermore, they appear in all SMRT cells (see figure below), which could indicate a property of the library used for sequencing.  It is not clear whether this bimodality is from the original DNA sample or an artifact of the long-fragment library construction process.

![plot of chunk lengthHistPerRG](figure/lengthHistPerRG.png) 


Read alignment
==============




We next examined the alignment of the reads to the reference (as performed by CSHL using `pbalign`(?)).  To guide our expectations, we compare the PacBio alignments of the long-read data to BWA alignments of short-read data Illumina data (paired-end, 76-bp, ~200 bp fragment size) from the same parasite.  However, please note that the DNA source library for these two experiments is not the same.  The Illumina data is from sample PG0051-C, the 3D7 isolate sequenced on an Illumina GA2 for the P.f. crosses project.  The coverage metrics over the whole of chromosome 8 are listed in the table below:

|id        |  median|    mean|     sd|
|:---------|-------:|-------:|------:|
|PacBio    |      81|   80.87|  14.36|
|Illumina  |     132|  118.75|  62.52|


![plot of chunk coverageDist](figure/coverageDist.png) 


We manually inspected the alignments in IGV across the entire length of chromosome 1.

![IGV1](figure/IGV_chr1.png)

Strikingly, the PacBio data appears to have uniform coverage across the entire length of chromosome 1, while the Illumina data shows many peaks and valleys along the same chromosome.  Zooming in closer (below), we can see that the PacBio reads are replete with insertions, typically one base long, as expected.

![IGV2](figure/IGV_indelerrors.png)


Genome accessibility
--------------------

We sought to examine regions of the P. falciparum genome that are inaccessible with Illumina reads but accessible with PacBio reads.  We computed coverage profiles across each autosome by computing coverage at every nucleotide using the GATK's `DepthOfCoverage` tool, and showing the minimum coverage value found in 2,000-bp bins.  We've plotted the PacBio coverage above the appropriate ideogram and the Illumina coverage underneath it (and flipped upside-down) in the plot below.  Red and gray regions indicate areas of the genome deemed to be inaccessible (from Alistair Miles's accessibility calculations on MalariaGen datasets, personal communication).

![plot of chunk showCoverageOverIdeogram](figure/showCoverageOverIdeogram.png) 


It is evident that the PacBio coverage is roughly uniform across the entire length of the chromosome.  In contrast, the Illumina coverage spikes and dips as it moves along, reaching zero coverage in many regions (especially the biologically interesting subtelomeric repetitive regions).  Let us examine a few of these places more closely.




### Centromere on chromosome 4
![plot of chunk coverageCentromere](figure/coverageCentromere.png) 


### All 28 masked telomeric regions
![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions1.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions2.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions3.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions4.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions5.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions6.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions7.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions8.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions9.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions10.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions11.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions12.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions13.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions14.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions15.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions16.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions17.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions18.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions19.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions20.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions21.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions22.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions23.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions24.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions25.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions26.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions27.png) ![plot of chunk coverageMaskedRegions](figure/coverageMaskedRegions28.png) 


Error rate along the length of the read
---------------------------------------

Finally, we examined the error rate as a function of distance from the center of the contig (this is a slightly more meaningful measure for contigs, where it is easy to see error piling up towards the end of contigs, associated with the end of a graph traversal, but we have the code so we might as well look).  Here, we take a read, fold it in half, and start measuring the rate of insertion and deletion errors as a function of distance from the center of the read.  The result is shown below.  Insertion and deletion errors are found uniformly along the length of the reads.  Towards the tail end, we start to see a lot more variance in this rate as reads that reach these lengths become more infrequent.

![plot of chunk errorsByPosition](figure/errorsByPosition.png) 


AsmTest1, the first assembly attempt
------------------------------------

To construct an initial draft assembly, we ran the `RS_HGAP_Assembly.2` secondary analysis protocol available in SMRT Portal.  Briefly, this processing protocol performs the following steps:

1. Extract subreads (genomic sequence absent of the SMRT bell adaptors used to circularize the fragment and enable the polymerase to read it in multiple passes).
2. Filter out low quality subreads.
3. Compute a subread length threshold such that subreads greater than or equal to this length provide roughly 30x genome coverage.
4. Select "seed" reads based on the computed subread length threshold.
5. Map all of the filtered subreads to the seed reads using `BLASR`.
6. Determine a consensus sequence from the subread alignments to the seed reads and preassemble (i.e. error-correct) the reads.
7. Assemble the preassembled reads using the Celera overlap-consensus-layout assembler.
8. Refine the assembly by mapping all raw data to the new assembly using `BLASR` and trimming low-quality ends of contigs.
9. Improve the continuity of the assembly and remove errors using the quality-aware consensus algorithm, `Quiver`. 

With the exception of the estimated genome size parameter (which we set to 23,000,000 bp), we left all default settings in this protocol unchanged.  Basic metrics on the resulting assembly, hereafter referred to as "AsmTest1", are presented in the table below.  For comparison, we have also provided metrics on every *P. falciparum* assembly currently available.

|id        |  numContigs|  minLength|  maxLength|  meanLength|      n50|  totalSequence|
|:---------|-----------:|----------:|----------:|-----------:|--------:|--------------:|
|3D7       |          16|       5967|    3291936|     1458302|  1687656|       23332831|
|IT        |          17|       6616|    3219929|     1351588|  1570953|       22976997|
|AsmTest1  |          34|      11443|    3293905|      697875|  1696391|       23727741|
|IGH-CR14  |         849|       2199|     120285|       25608|    37016|       21741172|
|HB3       |        1189|        201|     377975|       20402|    96469|       24258511|
|RAJ116    |        1199|       2042|      70306|       11765|    12998|       14106529|
|DD2       |        2837|        201|     102309|        7358|    19112|       20875591|
|V34.04    |        4329|        226|      16341|        3059|     3756|       13240777|
|D10       |        4471|        259|      19127|        2992|     3707|       13375079|
|K1        |        4772|        231|      18390|        2785|     3422|       13290906|
|7G8       |        4843|        204|      19000|        2948|     3832|       14278891|
|FCC-2     |        4956|        200|      17581|        2616|     3302|       12963854|
|RO-33     |        4991|        208|      19991|        2748|     3473|       13714138|
|D6        |        5011|        266|      15451|        2638|     3231|       13216528|
|SL        |        5193|        214|      55682|        2540|     3079|       13192745|
|VS.1      |        5856|        201|      22989|        3225|     4424|       18887633|
|PFCLIN    |       18711|       1001|      33813|        2366|     2992|       44265486|


The AsmTest1 assembly compares quite favorably to the best assemblies, with 34 compared to 3D7's 16 and IT's 17.  The longest chromosome in the *P. falciparum* genome is chromosome 14 (3291936 bp).  The longest contig in the AsmTest1 assembly appears roughly this length, suggesting that we may have assemblied the majority (or the entirety) of chromosome 14 in a single contig.

We compared AsmTest1 to the 3D7 canonical genome by performing an all-by-all (contigs vs. chromosomes) alignment with `MUMmer`.  Short, spurious alignments were filtered out.  The results are presented in the table below.  The near-entirety of chromosome 14 is assembled into a single contig.  On average, each chromosome is assembled into 4.4667 +/- 3.7391 contigs.  Each chromosome seems to be more-or-less fully recovered (note that the % coverage will add up to more than 100% in some cases as contigs may overlap slightly).

|REF          |QUERY                    |       S1|       S2|       E1|       E2|    LEN_1|    LEN_2|   COV_R|
|:------------|:------------------------|--------:|--------:|--------:|--------:|--------:|--------:|-------:|
|M76611       |scf7180000000116|quiver  |        1|     3868|     5967|     9834|     5967|     5967|  100.00|
|Pf3D7_01_v3  |scf7180000000096|quiver  |        3|     5347|    39095|    44519|    39093|    39173|    6.10|
|Pf3D7_01_v3  |scf7180000000113|quiver  |    26627|        1|   313647|   287313|   287021|   287313|   44.79|
|Pf3D7_01_v3  |scf7180000000113|quiver  |   313611|   287483|   640851|   614871|   327241|   327389|   51.06|
|Pf3D7_02_v3  |scf7180000000110|quiver  |     1465|        1|   947100|   944392|   945636|   944392|   99.85|
|Pf3D7_02_v3  |scf7180000000091|quiver  |   946546|       90|   947102|      651|      557|      562|    0.06|
|Pf3D7_03_v3  |scf7180000000108|quiver  |     2077|        2|  1063160|  1061347|  1061084|  1061346|   99.36|
|Pf3D7_03_v3  |scf7180000000091|quiver  |  1063232|      393|  1063882|     1040|      651|      648|    0.06|
|Pf3D7_03_v3  |scf7180000000110|quiver  |  1063772|   943840|  1065269|   945340|     1498|     1501|    0.14|
|Pf3D7_03_v3  |scf7180000000105|quiver  |  1064820|        1|  1065596|      784|      777|      784|    0.07|
|Pf3D7_03_v3  |scf7180000000102|quiver  |  1065611|  1320645|  1065841|  1320876|      231|      232|    0.02|
|Pf3D7_03_v3  |scf7180000000090|quiver  |  1066010|  3293080|  1066830|  3293905|      821|      826|    0.08|
|Pf3D7_04_v3  |scf7180000000121|quiver  |        1|      731|   379999|   380097|   379999|   379367|   31.65|
|Pf3D7_04_v3  |scf7180000000121|quiver  |   379839|   380192|   587301|   587781|   207463|   207590|   17.28|
|Pf3D7_04_v3  |scf7180000000121|quiver  |   587292|   587917|   961405|   962088|   374114|   374172|   31.16|
|Pf3D7_04_v3  |scf7180000000123|quiver  |   946163|        1|   962879|    16711|    16717|    16711|    1.39|
|Pf3D7_04_v3  |scf7180000000122|quiver  |   962849|    23394|  1200488|   261857|   237640|   238464|   19.80|
|Pf3D7_05_v3  |scf7180000000107|quiver  |    12341|     5299|  1104063|  1097830|  1091723|  1092532|   81.26|
|Pf3D7_05_v3  |scf7180000000107|quiver  |  1103711|  1097886|  1332679|  1326910|   228969|   229025|   17.04|
|Pf3D7_05_v3  |scf7180000000107|quiver  |  1332488|  1327077|  1343557|  1338180|    11070|    11104|    0.82|
|Pf3D7_06_v3  |scf7180000000102|quiver  |        3|  1260092|    59608|  1319754|    59606|    59663|    4.20|
|Pf3D7_06_v3  |scf7180000000102|quiver  |    59608|   459639|   858954|  1259331|   799347|   799693|   56.36|
|Pf3D7_06_v3  |scf7180000000102|quiver  |   858698|        1|  1318229|   459631|   459532|   459631|   32.40|
|Pf3D7_06_v3  |scf7180000000097|quiver  |  1310157|      525|  1418240|   108657|   108084|   108133|    7.62|
|Pf3D7_07_v3  |scf7180000000106|quiver  |      111|        1|    71658|    71472|    71548|    71472|    4.95|
|Pf3D7_07_v3  |scf7180000000104|quiver  |    59901|        1|  1444792|  1400860|  1384892|  1400860|   95.83|
|Pf3D7_07_v3  |scf7180000000096|quiver  |  1444713|    44370|  1445205|    44865|      493|      496|    0.03|
|Pf3D7_08_v3  |scf7180000000105|quiver  |        1|     1624|   460322|   462086|   460322|   460463|   31.25|
|Pf3D7_08_v3  |scf7180000000105|quiver  |   443167|   463696|  1472805|  1493537|  1029639|  1029842|   69.91|
|Pf3D7_09_v3  |scf7180000000118|quiver  |        4|     1265|  1483635|  1485303|  1483632|  1484039|   96.23|
|Pf3D7_09_v3  |scf7180000000119|quiver  |  1466769|        1|  1505606|    38801|    38838|    38801|    2.52|
|Pf3D7_09_v3  |scf7180000000120|quiver  |  1493354|        1|  1541735|    48394|    48382|    48394|    3.14|
|Pf3D7_10_v3  |scf7180000000099|quiver  |        1|      526|   899716|   900953|   899716|   900428|   53.31|
|Pf3D7_10_v3  |scf7180000000099|quiver  |   899565|   900991|   936722|   938160|    37158|    37170|    2.20|
|Pf3D7_10_v3  |scf7180000000099|quiver  |   936607|   940321|  1438085|  1441932|   501479|   501612|   29.71|
|Pf3D7_10_v3  |scf7180000000099|quiver  |  1436757|  1442382|  1546357|  1552102|   109601|   109721|    6.49|
|Pf3D7_10_v3  |scf7180000000099|quiver  |  1543910|  1553786|  1584566|  1595087|    40657|    41302|    2.41|
|Pf3D7_10_v3  |scf7180000000099|quiver  |  1584690|  1595210|  1672530|  1696279|    87841|   101070|    5.20|
|Pf3D7_10_v3  |scf7180000000100|quiver  |  1659128|        1|  1686332|    27224|    27205|    27224|    1.61|
|Pf3D7_10_v3  |scf7180000000122|quiver  |  1686821|   261596|  1687129|   261895|      309|      300|    0.02|
|Pf3D7_11_v3  |scf7180000000091|quiver  |        1|      682|      811|     1479|      811|      798|    0.04|
|Pf3D7_11_v3  |scf7180000000098|quiver  |      814|        1|    87429|    86664|    86616|    86664|    4.25|
|Pf3D7_11_v3  |scf7180000000097|quiver  |    76816|    90878|  1932717|  1947480|  1855902|  1856603|   91.05|
|Pf3D7_11_v3  |scf7180000000091|quiver  |  1918006|  2799800|  2038337|  2920153|   120332|   120354|    5.90|
|Pf3D7_12_v3  |scf7180000000092|quiver  |        1|      870|   497672|   498698|   497672|   497829|   21.91|
|Pf3D7_12_v3  |scf7180000000092|quiver  |   497593|   498759|   976053|   977553|   478461|   478795|   21.06|
|Pf3D7_12_v3  |scf7180000000092|quiver  |   974230|   981025|  1030827|  1037622|    56598|    56598|    2.49|
|Pf3D7_12_v3  |scf7180000000092|quiver  |  1031055|  1037523|  1510206|  1516966|   479152|   479444|   21.09|
|Pf3D7_12_v3  |scf7180000000092|quiver  |  1510147|  1516984|  1723147|  1730017|   213001|   213034|    9.38|
|Pf3D7_12_v3  |scf7180000000092|quiver  |  1707530|  1745542|  2271389|  2309609|   563860|   564068|   24.82|
|Pf3D7_13_v3  |scf7180000000091|quiver  |        1|     1432|  1547426|  1549862|  1547426|  1548431|   52.90|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  1547233|  1550332|  2211516|  2214838|   664284|   664507|   22.71|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2211575|  2216968|  2211646|  2217038|       72|       71|    0.00|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2211612|  2216385|  2212190|  2216966|      579|      582|    0.02|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2212189|  2217003|  2212274|  2217089|       86|       87|    0.00|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2212297|  2218035|  2328325|  2334091|   116029|   116057|    3.97|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2328302|  2335084|  2495387|  2502176|   167086|   167093|    5.71|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2495388|  2502326|  2561560|  2568544|    66173|    66219|    2.26|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2562118|  2569784|  2672689|  2680392|   110572|   110609|    3.78|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2672690|  2681292|  2672861|  2681463|      172|      172|    0.01|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2672850|  2680436|  2672982|  2680569|      133|      134|    0.00|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2672939|  2681323|  2673345|  2681732|      407|      410|    0.01|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2673331|  2680590|  2673527|  2680787|      197|      198|    0.01|
|Pf3D7_13_v3  |scf7180000000091|quiver  |  2674506|  2682267|  2807398|  2815198|   132893|   132932|    4.54|
|Pf3D7_13_v3  |scf7180000000095|quiver  |  2795392|        1|  2903479|   108907|   108088|   108907|    3.70|
|Pf3D7_13_v3  |scf7180000000117|quiver  |  2908643|        1|  2921615|    13054|    12973|    13054|    0.44|
|Pf3D7_14_v3  |scf7180000000090|quiver  |        1|      391|  3291936|  3293241|  3291936|  3292851|  100.00|



We visualized these alignments as a multi-dot plot.  Nearly all of the assembly content lines up with the canonical assembly (as evidenced from the bulk of the data adhering approximately to the f(x) = x line).  Occassional spurious points do exist, likely due to some repetitive homologous regions among chromosomes.

![plot of chunk dotplot](figure/dotplot.png) 


