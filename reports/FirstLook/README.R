
## ----setup, echo=FALSE---------------------------------------------------
opts_chunk$set(dev=c('png', 'pdf'));

color.il = rgb(250, 164, 26, maxColorValue=255);
color.pb = rgb(166, 42, 66, maxColorValue=255);

rgcolors = list(
    pb0 = rgb(162, 13, 30, 100, maxColorValue=255),
    pb1 = rgb(188, 13, 53, 100, maxColorValue=255),
    pb2 = rgb(222, 38, 76, 100, maxColorValue=255),
    pb3 = rgb(240, 120, 140, 100, maxColorValue=255),
    pb4 = rgb(0, 38, 28, 100, maxColorValue=255),
    pb5 = rgb(4, 77, 41, 100, maxColorValue=255),
    pb6 = rgb(22, 128, 57, 100, maxColorValue=255),
    pb7 = rgb(69, 191, 85, 100, maxColorValue=255)
)

dir.results = "../../results/master/";
dir.resources = "../../resources/";

par(mar=c(5, 4, 4, 2)+0.1);


## ----loadGenes, echo=FALSE-----------------------------------------------
genes = read.table(paste(dir.resources, "genes.txt", sep=""), header=FALSE, stringsAsFactors=FALSE);
colnames(genes) = c("chrom", "start", "stop", "gene", "description");

for (i in 1:nrow(genes)) {
    genes$description[i] = gsub("\\+", " ", URLdecode(genes$description[i]));
}

genes$color = "gray";
genes$color[grep("VAR|erythrocyte membrane protein|PfEMP1", genes$description)] = rgb(255, 0, 71, maxColorValue=255);
genes$color[grep("RIF|rifin", genes$description)] = rgb(255, 200, 33, maxColorValue=255);
genes$color[grep("stevor", genes$description)] = rgb(0, 155, 204, maxColorValue=255);

genes$text = "";
genes$text[grep("VAR|erythrocyte membrane protein|PfEMP1", genes$description)] = "v";
genes$text[grep("RIF|rifin", genes$description)] = "r";
genes$text[grep("stevor", genes$description)] = "s";


## ----readStats, echo=FALSE, results='asis'-------------------------------
ls.pb  = read.table(paste(dir.results, "/pb/lengthStats.txt",  sep=""), header=TRUE);
ls.pb0 = read.table(paste(dir.results, "/pb0/lengthStats.txt", sep=""), header=TRUE);
ls.pb1 = read.table(paste(dir.results, "/pb1/lengthStats.txt", sep=""), header=TRUE);
ls.pb2 = read.table(paste(dir.results, "/pb2/lengthStats.txt", sep=""), header=TRUE);
ls.pb3 = read.table(paste(dir.results, "/pb3/lengthStats.txt", sep=""), header=TRUE);
ls.pb4 = read.table(paste(dir.results, "/pb4/lengthStats.txt", sep=""), header=TRUE);
ls.pb5 = read.table(paste(dir.results, "/pb5/lengthStats.txt", sep=""), header=TRUE);
ls.pb6 = read.table(paste(dir.results, "/pb6/lengthStats.txt", sep=""), header=TRUE);
ls.pb7 = read.table(paste(dir.results, "/pb7/lengthStats.txt", sep=""), header=TRUE);

ls.pb$key = "total";

ls = rbind(ls.pb0, ls.pb1, ls.pb2, ls.pb3, ls.pb4, ls.pb5, ls.pb6, ls.pb7, ls.pb);

rgcounts = read.table(paste(dir.results, "pb/rgcounts.txt", sep=""), header=FALSE, stringsAsFactors=FALSE);
rownames(rgcounts) = gsub("RG:Z:", "", rgcounts$V1);
colnames(rgcounts) = c("RG", "alignedReads");

aligned = c(
    rgcounts["80f50c5537", "alignedReads"],
    rgcounts["101bc6c4f2", "alignedReads"],
    rgcounts["13b6a42da4", "alignedReads"],
    rgcounts["db860605d6", "alignedReads"],
    rgcounts["236e875a13", "alignedReads"],
    rgcounts["eeb951b4d9", "alignedReads"],
    rgcounts["438703e3ac", "alignedReads"],
    rgcounts["06f934595b", "alignedReads"]
);

aligned = c(aligned, sum(aligned));

ls = cbind(ls, alignedReads = aligned);
ls = cbind(ls, pctAligned = 100*ls$alignedReads/ls$numReads);

kable(ls[,c("key", "numReads", "minLength", "maxLength", "meanLength", "n50Value")]);

#80f50c5537 m140912_185114_42137_c100689352550000001823145102281580_s1_p0
#101bc6c4f2 m140912_220917_42137_c100689352550000001823145102281581_s1_p0
#13b6a42da4 m140913_012852_42137_c100689352550000001823145102281582_s1_p0
#db860605d6 m140913_044743_42137_c100689352550000001823145102281583_s1_p0
#236e875a13 m140918_073419_42137_c100716502550000001823136502281504_s1_p0
#eeb951b4d9 m140918_105250_42137_c100716502550000001823136502281505_s1_p0
#438703e3ac m140918_141535_42137_c100716502550000001823136502281506_s1_p0
#06f934595b m140918_173349_42137_c100716502550000001823136502281507_s1_p0


## ----lengthHist, echo=FALSE, eval=TRUE-----------------------------------
lh.pb = read.table(paste(dir.results, "/pb/lengthHist.txt", sep=""), header=TRUE);

plot(lh.pb$length, lh.pb$pb, type="l", lwd=3, col=color.pb, bty="n", cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="Length (bp)", ylab="Counts");
labels = c(paste("reads:", ls.pb$numReads[1]),
           paste("min length:", ls.pb$minLength[1], "bp"),
           paste("max length:", ls.pb$maxLength[1], "bp"),
           paste("mean length:", as.integer(ls.pb$meanLength[1]), "bp"),
           paste("n50:", ls.pb$n50Value[1], "bp"));
legend(22000, 6000, labels, cex=1.3, bty="n");

points(1100, lh.pb$pb[which(lh.pb$length == 1100)]+100, pch=6);
text(1100, lh.pb$pb[which(lh.pb$length == 1100)]+100, labels="~ 1.1 kb", pos=4);

points(8000, lh.pb$pb[which(lh.pb$length == 8000)]+100, pch=6);
text(8000, lh.pb$pb[which(lh.pb$length == 8000)]+100, labels="~ 8.0 kb", pos=4);


## ----lengthHistPerRG, echo=FALSE, eval=TRUE, warning=FALSE---------------
plot(0, 0, type="n", bty="n", cex=1.3, cex.lab=1.3, cex.axis=1.3, xlab="Length (bp)", ylab="Counts", xlim=c(0, max(lh.pb$length)), ylim=c(0, 1500));
for (rg in c("pb0", "pb1", "pb2", "pb3", "pb4", "pb5", "pb6", "pb7")) {
    lhsub = read.table(paste(dir.results, "/", rg, "/lengthHist.txt", sep=""), header=TRUE);

    points(lhsub$length, lhsub$pb, type="l", lwd=2, lty=1, col=rgcolors[[rg]]);
}
legend("topright", c( '@m140912_185114_42137_c100689352550000001823145102281580_s1_p0', '@m140912_220917_42137_c100689352550000001823145102281581_s1_p0', '@m140913_012852_42137_c100689352550000001823145102281582_s1_p0', '@m140913_044743_42137_c100689352550000001823145102281583_s1_p0', '@m140918_073419_42137_c100716502550000001823136502281504_s1_p0', '@m140918_105250_42137_c100716502550000001823136502281505_s1_p0', '@m140918_141535_42137_c100716502550000001823136502281506_s1_p0', '@m140918_173349_42137_c100716502550000001823136502281507_s1_p0'), col=c(rgcolors$pb0, rgcolors$pb1, rgcolors$pb2, rgcolors$pb3, rgcolors$pb4, rgcolors$pb5, rgcolors$pb6, rgcolors$pb7), lwd=2, lty=1, cex=0.6, bty="n");


## ----loadCoverage, echo=FALSE, cache=TRUE--------------------------------
cov.pb = read.table(paste(dir.results, "/pb/coverage.simple.txt", sep=""), header=FALSE, stringsAsFactors=FALSE);
names(cov.pb) = c("chrom", "start", "cov");

cov.il = read.table(paste(dir.results, "/il/coverage.simple.txt", sep=""), header=FALSE, stringsAsFactors=FALSE);
names(cov.il) = c("chrom", "start", "cov");


## ----coverageTable, echo=FALSE, results='asis', cache=TRUE---------------
cov.pb.chr8 = subset(cov.pb, chrom == "Pf3D7_08_v3")$cov;
cov.il.chr8 = subset(cov.il, chrom == "Pf3D7_08_v3")$cov;

cov.pb.mean = mean(cov.pb.chr8);
cov.pb.sd = sd(cov.pb.chr8);
cov.pb.median = median(cov.pb.chr8);

cov.il.mean = mean(cov.il.chr8);
cov.il.sd = sd(cov.il.chr8);
cov.il.median = median(cov.il.chr8);

cov.table = rbind(PacBio = c(cov.pb.median, cov.pb.mean, cov.pb.sd), Illumina = c(cov.il.median, cov.il.mean, cov.il.sd));
colnames(cov.table) = c("median", "mean", "sd");

kable(cov.table);


## ----coverageDist, echo=FALSE, results='asis', cache=TRUE----------------
h.pb = hist(cov.pb.chr8, breaks=0:(max(cov.pb.chr8) + 10), plot=FALSE);
h.il = hist(cov.il.chr8, breaks=0:(max(cov.il.chr8) + 10), plot=FALSE);

plot(0, 0, type="n", xlim=c(0, 250), ylim=c(0, max(h.pb$counts)), xlab="Coverage", ylab="Frequency", cex=1.3, cex.lab=1.3, cex.axis=1.3, bty="n");
points(h.pb$mids, h.pb$counts, type="l", lwd=2, col=color.pb);
points(h.il$mids, h.il$counts, type="l", lwd=2, col=color.il);

legend("topright", c("PacBio coverage", "Illumina coverage"), lwd=2, col=c(color.pb, color.il), bty="n", cex=1.3);


## ----showCoverageOverIdeogram, echo=FALSE, fig.height=6, fig.width=12, dpi=300, cache=TRUE----
chroms = sort(na.exclude(grep("apico|M76611", unique(cov.pb$chrom), invert=TRUE, value=TRUE)));
chroms.length = list();

for (chr in chroms) {
    chroms.length[[chr]] = max(subset(cov.pb, chrom == chr)$start);
}

mask = read.table(paste(dir.resources, "/tbl_telomere.txt", sep=""), header=TRUE, stringsAsFactors=FALSE);

access = read.table(paste(dir.resources, "/regions-20130225.txt", sep=""), header=FALSE, stringsAsFactors=FALSE);
names(access) = c("chrom", "start", "stop", "type");
access = subset(access, type != "Core");

par(mar=c(5, 8, 2, 1));
plot(0, 0, type="n", xlim=c(0, max(unlist(chroms.length)) + 500000), ylim=c(0, length(chroms.length) + 1), bty="n", xlab="Length (bp)", ylab="", yaxt="n", cex=1.3, cex.axis=1.3, cex.lab=1.3);

for (chr in chroms) {
    pos = as.integer(gsub("_v3", "", gsub("Pf3D7_", "", chr)));
    chrlength = chroms.length[[chr]];

    mtext(chr, side=2, at=pos, las=1, cex=1.3);

    submask = subset(mask, chrom == chr);
    rect(submask$co_pos_min, pos - 0.1, submask$co_pos_max, pos + 0.1, col="gray", border=NA);

    acc = subset(access, chrom == chr);
    rect(acc$start, pos - 0.1, acc$stop, pos + 0.1, col="red", border=NA);

    rect(0, pos - 0.1, chrlength, pos + 0.1);

    pb.covs = subset(cov.pb, chrom == chr);
    il.covs = subset(cov.il, chrom == chr);

    pb.start = pb.covs$start;
    il.start = il.covs$start;

    pb.max = max(pb.covs$cov);
    il.max = max(il.covs$cov);

    window = 2000;
    interval = seq(1, length(pb.start), by=window);

    pb.mincov = c();
    il.mincov = c();
    for (i in interval) {
        pb.mincov = c(pb.mincov, min(pb.covs$cov[i:(i+window)]));
        il.mincov = c(il.mincov, min(il.covs$cov[i:(i+window)]));
    }

    points(pb.start[interval], pos + 0.1 +  0.3*(pb.mincov / max(pb.mincov, na.rm=TRUE)), type="l", col=color.pb, lwd=0.5);
    points(il.start[interval], pos - 0.1 + -0.3*(il.mincov / max(il.mincov, na.rm=TRUE)), type="l", col=color.il, lwd=0.5);
}


## ----coverageZoom, echo=FALSE--------------------------------------------
showRegionalCoverage <- function(region) {
    region.window = 1000;
    region.left = ifelse(region$start - region.window < 0, 0, region$start - region.window);
    region.right = ifelse(region$stop + region.window > chroms.length[[region$chrom]], chroms.length[[region$chrom]], region$stop + region.window);

    region.cov.pb = subset(cov.pb, chrom == region$chrom & start >= region.left & start <= region.right);
    region.cov.il = subset(cov.il, chrom == region$chrom & start >= region.left & start <= region.right);

    region.max = max(region.cov.pb$cov, region.cov.il$cov);

    par(mar=c(5, 4, 4, 2)+0.1);
    plot(0, 0, type="n", xlim=c(region$start - region.window, region$stop + region.window), ylim=c(-15, region.max), bty="n", xlab="Position (bp)", ylab="Coverage", cex=1.3, cex.axis=1.3, cex.lab=1.3, main=paste(region$chrom, ":", region$start, "-", region$stop, sep=""));

    acc = subset(access, chrom == region$chrom)
    rect(acc$start, -5 + -region.max/100, acc$stop, -5 + region.max/100, col=rgb(255, 0, 0, 150, maxColorValue=255), border=NA);
    rect(0, -5 + -region.max/100, chroms.length[[region$chrom]], -5 + region.max/100);

    points(region.cov.pb$start, region.cov.pb$cov, type="l", col=color.pb, lwd=1);
    points(region.cov.il$start, region.cov.il$cov, type="l", col=color.il, lwd=1);

    overlappingGenes = subset(genes, chrom == region$chrom & ( (start <= region.left & stop >= region.left) | (start <= region.right & stop >= region.right) | (start <= region.left & stop >= region.right) | (start >= region.left & stop <= region.right) ));

    for (i in 1:nrow(overlappingGenes)) {
        overlappingGene = overlappingGenes[i,];
        rect(overlappingGene$start, -15 + -region.max/100, overlappingGene$stop, -15 + region.max/100, col=overlappingGene$color);
        text(overlappingGene$start + (overlappingGene$stop - overlappingGene$start)/2, -15, labels=overlappingGene$text, cex=0.7);
    }
}


## ----coverageCentromere, echo=FALSE, cache=TRUE, fig.height=8, fig.width=12----
region.chr4_centromere = subset(access, chrom == "Pf3D7_04_v3" & type == "Centromere")[1,];
region.chr4_centromere$start = region.chr4_centromere$start - 10000;
region.chr4_centromere$stop = region.chr4_centromere$stop + 10000;
showRegionalCoverage(region.chr4_centromere);


## ----coverageMaskedRegions, echo=FALSE, cache=TRUE, fig.height=8, fig.width=12----
mask.sorted = mask[order(mask$chrom, mask$co_pos_min),];

for (i in 1:nrow(mask.sorted)) {
    region.telomeric = mask.sorted[i,];
    colnames(region.telomeric) = c( "sample", "chrom", "co_pos_mid", "start", "stop", "co_pos_range", "cross", "co_from_parent", "co_to_parent" );
    showRegionalCoverage(region.telomeric);
}


## ----errorsByPosition, echo=FALSE----------------------------------------
f = read.table(paste(dir.results, "/pb/errors.per_position.txt", sep=""), header=TRUE, stringsAsFactors=FALSE);

col.insertions = rgb(219, 56, 46, 100, maxColorValue=255);
col.deletions = rgb(95, 169, 192, 100, maxColorValue=255);

plot(f$deletionRate, pch=19, cex=0.2, col=col.deletions, bty="n", xlab="Distance from middle of contig (bp)", ylab="Error rate", cex.lab=1.3, cex.axis=1.3);
points(f$insertionRate, pch=19, cex=0.2, col=col.insertions);
legend("topleft", c("Insertions", "Deletions"), fill=c(col.insertions, col.deletions), bty="n", cex=1.3);


## ----assemblyStats, echo=FALSE, results='asis'---------------------------
asmstats = read.table(paste(dir.results, "/AsmTest1/assemblies.stats", sep=""), header=TRUE, stringsAsFactors=FALSE);
asmstats = asmstats[order(asmstats$numContigs),];

kable(asmstats, row.names=FALSE);


## ----showCoords, echo=FALSE, results='asis'------------------------------
showCoords = read.table(paste(dir.results, "/AsmTest1/AsmTest1.filter.filter.coord_summary", sep=""), header=FALSE, stringsAsFactors=FALSE);
names(showCoords) = c( "S1", "E1", "S2", "E2", "LEN_1", "LEN_2", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "REF", "QUERY" );
showCoords$QUERY = gsub("\\|quiver", "", showCoords$QUERY);

kable(showCoords[,c("REF", "QUERY", "S1", "E1", "S2", "E2", "COV_R", "COV_Q")], row.names=FALSE);


## ----dotplot, echo=FALSE, fig.width=12, fig.height=12--------------------
gplines = readLines(paste(dir.results, "/AsmTest1/AsmTest1.filter.gp", sep=""));

xaxis.labels.unparsed = strsplit(gsub(" $", "", gsub("^ ", "", gsub("[,\\\"]", "", grep("Pf3D7", gplines, value=TRUE)))), " ");
xaxis.labels = c();
xaxis.at = c();
for (i in 1:length(xaxis.labels.unparsed)) {
    label = xaxis.labels.unparsed[[i]][1];
    at = xaxis.labels.unparsed[[i]][2];

    xaxis.labels = c(xaxis.labels, label);
    xaxis.at = c(xaxis.at, at);
}

yaxis.labels.unparsed = strsplit(gsub(" $", "", gsub("^ ", "", gsub("[,\\\"]", "", grep("scf", gplines, value=TRUE)))), " ");
yaxis.labels = c();
yaxis.at = c();
for (i in 1:length(yaxis.labels.unparsed)) {
    label = yaxis.labels.unparsed[[i]][1];
    at = yaxis.labels.unparsed[[i]][2];

    yaxis.labels = c(yaxis.labels, label);
    yaxis.at = c(yaxis.at, at);
}

xrange = as.integer(gsub("^\\[1:|\\]", "", unlist(strsplit(grep("xrange", gplines, value=TRUE), " "))[3]));
yrange = as.integer(gsub("^\\[1:|\\]", "", unlist(strsplit(grep("yrange", gplines, value=TRUE), " "))[3]));

fw = read.table(paste(dir.results, "/AsmTest1/AsmTest1.filter.fplot", sep=""), header=FALSE, stringsAsFactors=FALSE);
rc = read.table(paste(dir.results, "/AsmTest1/AsmTest1.filter.rplot", sep=""), header=FALSE, stringsAsFactors=FALSE);

par(mar=c(6, 9, 2, 2));
plot(0, 0, type="n", xlim=c(0, xrange), ylim=c(0, yrange), bty="n", xlab="", ylab="", xaxt="n", yaxt="n");
axis(1, at=xaxis.at, labels=xaxis.labels, cex.axis=0.6, las=3);
axis(2, at=yaxis.at, labels=yaxis.labels, cex.axis=0.6, las=1);

abline(v=as.numeric(xaxis.at), lty=2, col="gray90");
abline(h=as.numeric(yaxis.at), lty=2, col="gray90");

for (i in seq(1, nrow(fw), by=2)) {
    points(fw$V1[i:(i+1)], fw$V2[i:(i+1)], col="red", type="p", cex=0.7);
    points(fw$V1[i:(i+1)], fw$V2[i:(i+1)], col="red", type="l", lwd=2);
}

for (i in seq(1, nrow(rc), by=2)) {
    points(rc$V1[i:(i+1)], rc$V2[i:(i+1)], col="blue", type="p", cex=0.7);
    points(rc$V1[i:(i+1)], rc$V2[i:(i+1)], col="blue", type="l", lwd=2);
}


## ----variants, echo=FALSE, results='asis'--------------------------------
variants = read.table(paste(dir.results, "/AsmTest1/AsmTest1.filter.filter.variant_summary", sep=""), header=FALSE, stringsAsFactors=FALSE);
names(variants) = c( "P1", "REF", "ALT", "P2", "BUFF", "DIST", "R", "Q", "LEN_R", "LEN_Q", "FRM", "UNKNOWN", "REFERENCE", "QUERY" );

variants.summary = cbind("SNP" = table(subset(variants, REFERENCE != "M76611" & REF != "." & ALT != ".")$REFERENCE),
                         "INS" = table(subset(variants, REFERENCE != "M76611" & REF == "." & ALT != ".")$REFERENCE),
                         "DEL" = table(subset(variants, REFERENCE != "M76611" & REF != "." & ALT == ".")$REFERENCE));
variants.summary = rbind(variants.summary, c(sum(variants.summary[,"SNP"]), sum(variants.summary[,"INS"]), sum(variants.summary[,"DEL"])));

chrom.lengths = unique(variants[,c("REFERENCE", "LEN_R")]);
rownames(chrom.lengths) = chrom.lengths$REFERENCE;

variants.pctsummary = variants.summary[1:nrow(variants.summary) - 1,];
for (chrom in chrom.lengths$REFERENCE) {
    if (chrom != "M76611") {
        variants.pctsummary[chrom,] = 100.0 * variants.summary[chrom,] / chrom.lengths[chrom, "LEN_R"];
    }
}
variants.pctsummary = rbind(variants.pctsummary, 100.0*c(sum(variants.summary[,"SNP"]), sum(variants.summary[,"INS"]), sum(variants.summary[,"DEL"]))/23332831);

variants.fullsummary = variants.summary;
variants.fullsummary[,1] = sprintf("%s (%.2f%%)", variants.summary[,1], variants.pctsummary[,1]);
variants.fullsummary[,2] = sprintf("%s (%.2f%%)", variants.summary[,2], variants.pctsummary[,2]);
variants.fullsummary[,3] = sprintf("%s (%.2f%%)", variants.summary[,3], variants.pctsummary[,3]);

kable(variants.fullsummary);


## ----varTable, echo=FALSE, results='asis'--------------------------------
varTable = read.table(paste(dir.results, "/AsmTest1/AsmTest1.varGenes.table", sep=""), header=FALSE, stringsAsFactors=FALSE);
colnames(varTable) = c( "geneName", "flags", "contig", "start", "MQ", "CIGAR", "QUAL", "UKN1", "UKN2", "NM", "MD", "AS" );
varTable$NM = as.integer(gsub("NM:i:", "", varTable$NM));

numIns = c();
numDel = c();

for (i in 1:nrow(varTable)) {
    insAndDels = table(unlist(strsplit(gsub("\\d+", "", gsub("\\d+M", "", varTable$CIGAR[i])), "")));

    numIns = c(numIns, ifelse(is.na(insAndDels["I"]), 0, insAndDels["I"]));
    numDel = c(numDel, ifelse(is.na(insAndDels["D"]), 0, insAndDels["D"]));
}

varTable = cbind(varTable, NI = numIns, ND = numDel);

kable(varTable[,c("geneName", "contig", "start", "NM", "NI", "ND", "CIGAR")]);

numPerfect = nrow(subset(varTable, NM == 0 & NI == 0 & ND == 0));


## ----varMismatchLocationSummary, echo=FALSE------------------------------
varExons = read.table(paste(dir.results, "/AsmTest1/AsmTest1.varExons.table", sep=""), header=FALSE, stringsAsFactors=FALSE);
colnames(varExons) = c( "exonName", "flags", "contig", "start", "MQ", "CIGAR", "QUAL", "UKN1", "UKN2", "NM", "MD", "AS" );
varExons$NM = as.integer(gsub("NM:i:", "", varExons$NM));
varExons$geneName = unlist(lapply(strsplit(gsub("exon_", "", varExons$exonName), "-"), function(x) { return(x[1]); }));
varExons$exonNumber = unlist(lapply(strsplit(gsub("exon_", "", varExons$exonName), "-"), function(x) { return(x[2]); }));

numIns = c();
numDel = c();

for (i in 1:nrow(varExons)) {
    insAndDels = table(unlist(strsplit(gsub("\\d+", "", gsub("\\d+M", "", varExons$CIGAR[i])), "")));

    numIns = c(numIns, ifelse(is.na(insAndDels["I"]), 0, insAndDels["I"]));
    numDel = c(numDel, ifelse(is.na(insAndDels["D"]), 0, insAndDels["D"]));
}

varExons = cbind(varExons, NI = numIns, ND = numDel);

genesToProcess = unique(varExons$geneName[!(varExons$geneName %in% subset(varExons, MQ < 10)$geneName)]);

list.nm = list();
list.ni = list();
list.nd = list();
list.ne = list();
for (gene in genesToProcess) {
    nm.total = subset(varTable, geneName == gene)$NM;
    ni.total = subset(varTable, geneName == gene)$NI;
    nd.total = subset(varTable, geneName == gene)$ND;

    nm.exon1 = subset(varExons, geneName == gene & exonNumber == 1)$NM;
    ni.exon1 = subset(varExons, geneName == gene & exonNumber == 1)$NI;
    nd.exon1 = subset(varExons, geneName == gene & exonNumber == 1)$ND;

    nm.exon2 = subset(varExons, geneName == gene & exonNumber == 2)$NM;
    ni.exon2 = subset(varExons, geneName == gene & exonNumber == 2)$NI;
    nd.exon2 = subset(varExons, geneName == gene & exonNumber == 2)$ND;

    nm.intron = nm.total - nm.exon1 - nm.exon2;
    ni.intron = ni.total - ni.exon1 - ni.exon2;
    nd.intron = nd.total - nd.exon1 - nd.exon2;

    list.nm[["total"]] = c(list.nm[["total"]], nm.total);
    list.nm[["exon1"]] = c(list.nm[["exon1"]], nm.exon1);
    list.nm[["exon2"]] = c(list.nm[["exon2"]], nm.exon2);
    list.nm[["intron"]] = c(list.nm[["intron"]], nm.intron);

    list.ni[["total"]] = c(list.ni[["total"]], ni.total);
    list.ni[["exon1"]] = c(list.ni[["exon1"]], ni.exon1);
    list.ni[["exon2"]] = c(list.ni[["exon2"]], ni.exon2);
    list.ni[["intron"]] = c(list.ni[["intron"]], ni.intron);

    list.nd[["total"]] = c(list.nd[["total"]], nd.total);
    list.nd[["exon1"]] = c(list.nd[["exon1"]], nd.exon1);
    list.nd[["exon2"]] = c(list.nd[["exon2"]], nd.exon2);
    list.nd[["intron"]] = c(list.nd[["intron"]], nd.intron);

    list.ne[["total"]] = c(list.ne[["total"]], nm.total + ni.total + nd.total);
    list.ne[["exon1"]] = c(list.ne[["exon1"]], nm.exon1 + ni.exon1 + nd.exon1);
    list.ne[["exon2"]] = c(list.ne[["exon2"]], nm.exon2 + ni.exon2 + nd.exon2);
    list.ne[["intron"]] = c(list.ne[["intron"]], nm.intron + ni.intron + nd.intron);
}

matrix.ne = cbind(total = sum(list.ne[["total"]]), intron = sum(list.ne[["intron"]]), exon1 = sum(list.ne[["exon1"]]), exon2 = sum(list.ne[["exon2"]]));
matrix.nm = cbind(total = sum(list.nm[["total"]]), intron = sum(list.nm[["intron"]]), exon1 = sum(list.nm[["exon1"]]), exon2 = sum(list.nm[["exon2"]]));
matrix.ni = cbind(total = sum(list.ni[["total"]]), intron = sum(list.ni[["intron"]]), exon1 = sum(list.ni[["exon1"]]), exon2 = sum(list.ni[["exon2"]]));
matrix.nd = cbind(total = sum(list.nd[["total"]]), intron = sum(list.nd[["intron"]]), exon1 = sum(list.nd[["exon1"]]), exon2 = sum(list.nd[["exon2"]]));

matrix.all = rbind(all = matrix.ne, mismatches = matrix.nm, insertions = matrix.ni, deletions = matrix.nd);
rownames(matrix.all) = c("all", "  mismatches", "  insertions", "  deletions");


## ----showVarMismatchLocationSummary, echo=FALSE, results='asis'----------
kable(matrix.all);


## ----sessionInfo---------------------------------------------------------
sessionInfo();


