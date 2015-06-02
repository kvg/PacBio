Lowering the barrier to entry for long-read sequencing with whole-genome amplification
==

Long-read sequencing can produce very long genomic fragments with an N50 of ~10 kb, with some fragments as long as ~40 kb.  However, they typically require a large amount of high molecular weight DNA (5 to 15 ug), limiting their applications.  Particular to malaria parasite sequencing, culturing a sufficient number of Plasmodium falciparum parasites to obtain the desired gDNA yield is enormously difficult and time-consuming.  Alternatively, one can consider whole-genome amplification on a small amount of DNA - say ~1 ng - to the 15 ug level.  Standard WGA kits employing the Taq polymerase are inappropriate for such a task; Taq-produced amplicons are typically up to 3 or 4 kb, negating much of the value of long-read sequencing.  Instead, multiple displacement amplification (MDA) employs Phi-29, a polymerase with much higher replication fidelity and capable of generating amplicons greater than 20 kb in length.

We explored the use of MDA-based WGA to produce a draft-quality genome assembly of a malaria parasite genome.  Specifically, we amplified 1 ng of 3D7 to 15 ug and performed PacBio sequencing on the resulting amplicons.


```
## Warning in file(file, "rt"): cannot open file
## '../../results/wga/amplified/amplified.asmStats': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```
## Error in rbind(s.unamp, s.amp): object 's.amp' not found
```

```
## Error in is.data.frame(x): object 's' not found
```


```
## Warning in file(file, "rt"): cannot open file
## '../../results/wga/amplified/amplified.lengths.txt': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```
## Error in hist(l.amp$V1, breaks = seq(0, max(l.amp$V1) + 100, by = 100), : object 'l.amp' not found
```

```
## Error in plot(h.amp$mids, h.amp$density, col = "red", type = "l", lwd = 3, : object 'h.amp' not found
```

```
## Error in plot.xy(xy.coords(x, y), type = type, ...): plot.new has not been called yet
```

```
## Error in eval(expr, envir, enclos): object 'h.amp' not found
```

```
## Error in text(s["unamplified", "n50"], h.y, labels = paste("n50 =", s["unamplified", : object 's' not found
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): object 's' not found
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): object 's' not found
```

```
## Error in text(s["amplified", "n50"], h.y - (h.y/10), labels = paste("n50 =", : object 's' not found
```

```
## Error in strwidth(legend, units = "user", cex = cex, font = text.font): plot.new has not been called yet
```


```
## Warning in file(file, "rt"): cannot open file
## '../../results/wga/amplified/amplified.chimeras.k21.txt': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```
## Error in apply(k.amp, 1, isChimericByKmers): object 'k.amp' not found
```

```
## Warning in file(file, "rt"): cannot open file
## '../../results/wga/unamplified/unamplified.chimeras.k21.txt': No such file
## or directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```
## Error in apply(k.unamp, 1, isChimericByKmers): object 'k.unamp' not found
```


```
## Warning in file(file, "rt"): cannot open file
## '../../results/wga/amplified/amplified.sam.table': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```
## Error in k.amp$isAligned = TRUE: object 'k.amp' not found
```

```
## Error in k.amp[which(k.amp$readName %in% unique(subset(a.amp, chr == "*")$readName)), : object 'k.amp' not found
```

```
## Error in k.amp$chimericByBwa = FALSE: object 'k.amp' not found
```

```
## Error in k.amp[which(k.amp$readName %in% unique(subset(a.amp, isChimeric == : object 'k.amp' not found
```

```
## Error in k.unamp$isAligned = TRUE: object 'k.unamp' not found
```

```
## Error in k.unamp[which(k.unamp$readName %in% unique(subset(a.unamp, chr == : object 'k.unamp' not found
```

```
## Error in k.unamp$chimericByBwa = FALSE: object 'k.unamp' not found
```

```
## Error in k.unamp[which(k.unamp$readName %in% unique(subset(a.unamp, isChimeric == : object 'k.unamp' not found
```
