Lowering the barrier to entry for long-read sequencing with whole-genome amplification
==

Long-read sequencing can produce very long genomic fragments with an N50 of ~10-20 kb, with some fragments as long as ~50 kb.  However, they typically require a large amount of high molecular weight DNA (5 to 15 ug), limiting their applications.  Particular to malaria parasite sequencing, culturing a sufficient number of *Plasmodium falciparum* parasites to obtain the desired gDNA yield is enormously difficult and time-consuming.  Alternatively, one can consider whole-genome amplification on a small amount of DNA - say ~1 ng - to the requisite 15 ug level.  Standard WGA kits employing the Taq polymerase are inappropriate for such a task; Taq-produced amplicons are typically up to 3 or 4 kb, negating much of the value of long-read sequencing.  Instead, multiple displacement amplification (MDA) employs Phi-29, a polymerase with much higher replication fidelity and capability to generate amplicons greater than 20 kb in length.

We explored the use of MDA-based WGA to produce a draft-quality genome assembly of a malaria parasite genome.  Specifically, we amplified 1 ng of 3D7 to 15 ug and performed PacBio sequencing on the resulting amplicons.



![plot of chunk lengthDist](figure/lengthDist-1.png) 


|key            | numReads| minLength| maxLength| meanLength| n50Value|
|:--------------|--------:|---------:|---------:|----------:|--------:|
|3D7_CSHL_A01   |    46946|        50|     38351|   9622.495|    13685|
|3D7_CSHL_C01   |    49256|        50|     44451|  10274.992|    14412|
|3D7_CSHL_D01   |    40379|        50|     42929|  10874.385|    14128|
|3D7_PacBio_A01 |    63785|        50|     40660|   9725.101|    14767|
|3D7_PacBio_B01 |    71630|        51|     43816|   9932.797|    14851|
|3D7_PacBio_C01 |    89314|        50|     47379|  10767.940|    15354|
|3D7_PacBio_D01 |    91311|        50|     48366|  10733.769|    15334|
|3D7_PacBio_E01 |    94691|        50|     45388|  10271.377|    14274|
|3D7_PacBio_F01 |    85031|        50|     45460|  10228.038|    14368|
|3D7_PacBio_G01 |    47916|        50|     43162|  10397.394|    14460|
|3D7_PacBio_H01 |    68161|        50|     44561|  10278.915|    14452|






```
## Warning: NAs introduced by coercion
```

```
## Warning: NAs introduced by coercion
```

![plot of chunk showCoverageOverIdeogram](figure/showCoverageOverIdeogram-1.png) 


|id          | numContigs| minLength| maxLength|     n50|    ng50| totalSequence|
|:-----------|----------:|---------:|---------:|-------:|-------:|-------------:|
|3D7         |         16|      5967|   3291936| 1687656| 1687656|      23332831|
|unamplified |         34|     11443|   3293905| 1696391| 3293905|      23727741|
|IT          |         56|       580|   3219929| 1570953| 3219929|      22953932|
|amplified   |         84|       618|   1805358|  887768| 1805358|      22914917|
|IGH-CR14    |        849|      2199|    120285|   37016|  120285|      21741172|
|HB3         |       1189|       201|    377975|   96469|  377975|      24258511|
|RAJ116      |       1199|      2042|     70306|   12998|   70306|      14106529|
|DD2         |       2837|       201|    102309|   19112|  102309|      20875591|
|V34.04      |       4329|       226|     16341|    3756|   16341|      13240777|
|D10         |       4471|       259|     19127|    3707|   19127|      13375079|
|K1          |       4772|       231|     18390|    3422|   18390|      13290906|
|7G8         |       4843|       204|     19000|    3832|   19000|      14278891|
|FCC-2       |       4956|       200|     17581|    3302|   17581|      12963854|
|RO-33       |       4991|       208|     19991|    3473|   19991|      13714138|
|D6          |       5011|       266|     15451|    3231|   15451|      13216528|
|SL          |       5193|       214|     55682|    3079|   55682|      13192745|
|VS.1        |       5856|       201|     22989|    4424|   22989|      18887633|
|PFCLIN      |      18711|      1001|     33813|    2992|   33813|      44265486|

![dotplot_unamp]( figure/3D7.unamplified.filter.png )
![dotplot_amp]( figure/3D7.amplified.filter.png )







