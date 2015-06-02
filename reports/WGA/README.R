
## ----asmStats, echo=FALSE, results='asis'--------------------------------
readAsmStats <- function(file, name) {
    d = t(read.table(file, header=FALSE, stringsAsFactors=FALSE, row.names=1));
    colnames(d) = gsub(":", "", colnames(d));
    rownames(d) = c(name);

    return(d);
}

s.unamp = readAsmStats("../../results/wga/unamplified/unamplified.asmStats", "unamplified");
s.amp   = readAsmStats("../../results/wga/amplified/amplified.asmStats", "amplified");

s = rbind(s.unamp, s.amp);

kable(s);


## ----readLengths, echo=FALSE---------------------------------------------
l.unamp = read.table("../../results/wga/unamplified/unamplified.lengths.txt", header=FALSE, stringsAsFactors=FALSE);
l.amp = read.table("../../results/wga/amplified/amplified.lengths.txt", header=FALSE, stringsAsFactors=FALSE);

h.unamp = hist(l.unamp$V1, breaks=seq(0, max(l.unamp$V1) + 100, by=100), plot=FALSE);
h.amp   = hist(l.amp$V1, breaks=seq(0, max(l.amp$V1) + 100, by=100), plot=FALSE);

plot(h.amp$mids, h.amp$density, col="red", type="l", lwd=3, bty="n", xlab="Read length (bp)", ylab="Density", cex=1.3, cex.axis=1.3, cex.lab=1.3);
points(h.unamp$mids, h.unamp$density, col="blue", type="l", lwd=3);

h.y = abs(max(h.amp$density) - max(h.unamp$density));

text(s["unamplified", "n50"], h.y, labels=paste("n50 =", s["unamplified", "n50"], "bp"), col="blue", pos=4, offset=-0.06);
abline(v=s["unamplified", "n50"], col="blue", lty=3);
abline(v=s["amplified", "n50"], col="red", lty=3);
text(s["amplified", "n50"], h.y - (h.y/10), labels=paste("n50 =", s["amplified", "n50"], "bp"), col="red", pos=4, offset=-0.06);

legend("topright", c("Amplified", "Unamplified"), col=c("red", "blue"), lwd=3, cex=1.3, bty="n");


