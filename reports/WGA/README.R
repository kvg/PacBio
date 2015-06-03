## ----setup, echo=FALSE---------------------------------------------------
library("seqinr");

results.dir = "../../results/wga";
data.dir = "../../data";

ref = "../../resources/references/3D7/PlasmoDB-9.0_Pfalciparum3D7_Genome.sorted.fasta";

## ----asmStats, echo=FALSE, results='asis'--------------------------------
asmStats = read.table(paste(results.dir, "assembly.stats", sep="/"), header=TRUE);

kable(asmStats);

## ----scaffoldLengths, echo=FALSE-----------------------------------------
n50 <- function(f) {
    lengths = sort(as.vector(unlist(lapply(f, function(x) { return(length(x)) }))));

    totalLength = sum(lengths);

    n50Length = 0;
    n50Value = 0;
    for (length in lengths) {
        n50Length = n50Length + length;
        n50Value = length;

        if (n50Length >= totalLength/2) {
            break;
        }
    }

    return(n50Value);
}

ng50 <- function(f, r) {
    lengths = sort(as.vector(unlist(lapply(f, function(x) { return(length(x)) }))));
    lengths.ref = sort(as.vector(unlist(lapply(r, function(x) { return(length(x)) }))));

    totalLength = sum(lengths.ref);

    n50Length = 0;
    n50Value = 0;
    for (length in lengths) {
        n50Length = n50Length + length;
        n50Value = length;

        if (n50Length >= totalLength/2) {
            break;
        }
    }

    return(n50Value);
}

f.ref   = read.fasta(ref);
f.unamp = read.fasta(paste(data.dir, "ASMTest1.polished_assembly.fasta", sep="/"));
f.amp   = read.fasta(paste(data.dir, "WGATest5.polished_assembly.fasta", sep="/"));

f.ref.size = sum(unlist(lapply(f.ref, function(x) { return(length(x)) })));

#plot(ecdf(as.vector(sort(unlist(lapply(f.amp, function(x) { return(length(x)) }))))))
#plot(ecdf(as.vector(sort(unlist(lapply(f.unamp, function(x) { return(length(x)) }))))))

