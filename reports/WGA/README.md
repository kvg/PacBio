Lowering the barrier to entry for long-read sequencing with whole-genome amplification
==

Long-read sequencing can produce very long genomic fragments with an N50 of ~10-20 kb, with some fragments as long as ~50 kb.  However, they typically require a large amount of high molecular weight DNA (5 to 15 ug), limiting their applications.  Particular to malaria parasite sequencing, culturing a sufficient number of *Plasmodium falciparum* parasites to obtain the desired gDNA yield is enormously difficult and time-consuming.  Alternatively, one can consider whole-genome amplification on a small amount of DNA - say ~1 ng - to the requisite 15 ug level.  Standard WGA kits employing the Taq polymerase are inappropriate for such a task; Taq-produced amplicons are typically up to 3 or 4 kb, negating much of the value of long-read sequencing.  Instead, multiple displacement amplification (MDA) employs Phi-29, a polymerase with much higher replication fidelity and capability to generate amplicons greater than 20 kb in length.

We explored the use of MDA-based WGA to produce a draft-quality genome assembly of a malaria parasite genome.  Specifically, we amplified 1 ng of 3D7 to 15 ug and performed PacBio sequencing on the resulting amplicons.






```
## Warning: NAs introduced by coercion
```

```
## Warning: NAs introduced by coercion
```

![plot of chunk showCoverageOverIdeogram](figure/showCoverageOverIdeogram-1.png) 


|id          | numContigs| minLength| maxLength|    ng50| totalSequence|
|:-----------|----------:|---------:|---------:|-------:|-------------:|
|amplified   |         84|       618|   1805358|  887768|      22914917|
|unamplified |         34|     11443|   3293905| 1696391|      23727741|



