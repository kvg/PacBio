#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use lib "$FindBin::Bin/lib";

use ParseArgs;
use DM;

my %args = &getCommandArguments(
    'ANALYSIS'     => undef,
    'DRY_RUN'    => 1,
    'NUM_JOBS'   => 1,
    'KEEP_GOING' => 0,
    'CLUSTER'    => 'localhost',
    'QUEUE'      => 'localhost',
);

my $cwd = getcwd();
my $projectName = basename($cwd);

my $binDir = "$cwd/bin";
my $dataDir = "$cwd/data";
my $listsDir = "$cwd/lists";
my $resultsDir = "$cwd/results/$args{'ANALYSIS'}";
my $reportsDir = "$cwd/reports";
my $resourcesDir = "$cwd/resources";
my $scriptsDir = "$cwd/scripts";
my $logsDir = "$resultsDir/logs";

my $dm = new DM(
    'dryRun'     => $args{'DRY_RUN'},
    'numJobs'    => $args{'NUM_JOBS'},
    'keepGoing'  => $args{'KEEP_GOING'},
    'cluster'    => $args{'CLUSTER'},
    'queue'      => $args{'QUEUE'},
    'outputFile' => "$resultsDir/dm.log",
);

my $indiana8 = "java -Xmx8g -jar $binDir/indiana.jar";
#my $gatk8 = "java -Xmx8g -jar bin/GenomeAnalysisTK.v3.2.jar";
my $gatk8 = "java -Xmx8g -jar /data1/users/kiran/repositories/GATK-Protected/protected/gatk-package-distribution/target/gatk-package-distribution-3.2.jar";
my $nucmer = "~/opt/MUMmer3.23/nucmer";
my $showaligns = "~/opt/MUMmer3.23/show-aligns";
my $showcoords = "~/opt/MUMmer3.23/show-coords";
my $showdiff = "~/opt/MUMmer3.23/show-diff";
my $showsnps = "~/opt/MUMmer3.23/show-snps";
my $showtiling = "~/opt/MUMmer3.23/show-tiling";
my $deltafilter = "~/opt/MUMmer3.23/delta-filter";
my $mummerplot = "~/opt/MUMmer3.23/mummerplot";
my $icorn2 = "~/opt/icorn2-v0.95//icorn2.sh";

my %bams = (
    'pb' => "$dataDir/Pfalc3D7/aligned_reads.bam",
    'il' => "$dataDir/PfCross/PG0051-C.ERR019061.bam",
);

my %seqs = (
    'pb' => "$dataDir/Pfalc3D7/filtered_subreads.fastq",
    'il' => "$dataDir/PfCross/ERR019061_1.fastq.gz",
);

my %fastqs = (
    'il' => {
        'end1' => "$dataDir/PfCross/ERR019061_1.fastq.gz",
        'end2' => "$dataDir/PfCross/ERR019061_2.fastq.gz"
    },
);

my %asms = (
    'AsmTest1' => "$dataDir/ASMTest1.polished_assembly.fasta",
    'ICORN2_01' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.1.fasta',
    'ICORN2_02' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.2.fasta',
    'ICORN2_03' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.3.fasta',
    'ICORN2_04' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.4.fasta',
    'ICORN2_05' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.5.fasta',
    'ICORN2_06' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.6.fasta',
    'ICORN2_07' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.7.fasta',
    'ICORN2_08' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.8.fasta',
    'ICORN2_09' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.9.fasta',
    'ICORN2_10' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.10.fasta',
    'ICORN2_11' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.11.fasta',
    'ICORN2_12' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.12.fasta',
    'ICORN2_13' => 'scratch/icorn2/ICORN2.ASMTest1.polished_assembly.13.fasta',
);

my %refs = (
    'pb' => "$resourcesDir/3D7.cshl.fasta",
    'il' => "/home/kiran/ngs/references/plasmodium/falciparum/3D7/PlasmoDB-9.0/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta",
);

my %existingAsms = (
    "3D7" => "/home/kiran/ngs/references/plasmodium/falciparum/3D7/PlasmoDB-9.0/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta",
    "7G8" => "/home/kiran/ngs/references/plasmodium/falciparum/7G8/BroadInstitute/plasmodium_falciparum__isolate_7g8__1_supercontigs.fasta",
    "D10" => "/home/kiran/ngs/references/plasmodium/falciparum/D10/BroadInstitute/plasmodium_falciparum__isolate_d10__1_supercontigs.fasta",
    "D6" => "/home/kiran/ngs/references/plasmodium/falciparum/D6/BroadInstitute/plasmodium_falciparum__isolate_d6__1_supercontigs.fasta",
    "DD2" => "/home/kiran/ngs/references/plasmodium/falciparum/Dd2/BroadInstitute/plasmodium_falciparum__isolate_dd2__1_supercontigs.fasta",
    "FCC-2" => "/home/kiran/ngs/references/plasmodium/falciparum/FCC-2.Hainan/BroadInstitute/hainan__1_supercontigs.fasta",
    "HB3" => "/home/kiran/ngs/references/plasmodium/falciparum/HB3/BroadInstitute/plasmodium_falciparum__isolate_hb3__1_supercontigs.fasta",
    "IGH-CR14" => "/home/kiran/ngs/references/plasmodium/falciparum/IGH-CR14/BroadInstitute/plasmodium_falciparum_igh-cr14_nucleus_1_supercontigs.fasta",
    "IT" => "/home/kiran/ngs/references/plasmodium/falciparum/IT/PlasmoDB-9.0/PlasmoDB-9.0_PfalciparumIT_Genome.split.fasta",
    "K1" => "/home/kiran/ngs/references/plasmodium/falciparum/K1/BroadInstitute/plasmodium_falciparum__isolate_k1__1_supercontigs.fasta",
    "PFCLIN" => "/home/kiran/ngs/references/plasmodium/falciparum/PFCLIN/SangerInstitute/PFCLIN.20080302.contigs.fasta",
    "RAJ116" => "/home/kiran/ngs/references/plasmodium/falciparum/RAJ116/BroadInstitute/plasmodium_falciparum_raj116_nucleus_1_supercontigs.fasta",
    "RO-33" => "/home/kiran/ngs/references/plasmodium/falciparum/RO-33/BroadInstitute/plasmodium_falciparum__isolate_ro-33__1_supercontigs.fasta",
    "SL" => "/home/kiran/ngs/references/plasmodium/falciparum/SL/BroadInstitute/plasmodium_falciparum__isolate_santa_lucia__1_supercontigs.fasta",
    "V34.04" => "/home/kiran/ngs/references/plasmodium/falciparum/Senegal_V34.04/BroadInstitute/plasmodium_falciparum__isolate_senegal_v34.04__1_supercontigs.fasta",
    "VS.1" => "/home/kiran/ngs/references/plasmodium/falciparum/VS.1/BroadInstitute/1__1_supercontigs.fasta",
);

my %rois = (
    'varGenes' => "$resourcesDir/var.3D7.annotated.gff",
    'varExons' => "$resourcesDir/var.3D7.annotated.gff",
    'rifinGenes' => "$resourcesDir/rifin.3D7.gff",
    'rifinExons' => "$resourcesDir/rifin.3D7.gff",
    'stevorGenes' => "$resourcesDir/stevor.3D7.gff",
    'stevorExons' => "$resourcesDir/stevor.3D7.gff",
    'allGenes' => "$resourcesDir/all.3D7.gff",
    'allExons' => "$resourcesDir/all.3D7.gff",
);

my %roiTypes = (
    'varGenes' => "gene",
    'varExons' => "exon",
    'rifinGenes' => "gene",
    'rifinExons' => "exon",
    'stevorGenes' => "gene",
    'stevorExons' => "exon",
    'allGenes' => "gene",
    'allExons' => "exon",
);

my %pbrgs = (
    'pb0' => '@m140912_185114_42137_c100689352550000001823145102281580_s1_p0',
    'pb1' => '@m140912_220917_42137_c100689352550000001823145102281581_s1_p0',
    'pb2' => '@m140913_012852_42137_c100689352550000001823145102281582_s1_p0',
    'pb3' => '@m140913_044743_42137_c100689352550000001823145102281583_s1_p0',
    'pb4' => '@m140918_073419_42137_c100716502550000001823136502281504_s1_p0',
    'pb5' => '@m140918_105250_42137_c100716502550000001823136502281505_s1_p0',
    'pb6' => '@m140918_141535_42137_c100716502550000001823136502281506_s1_p0',
    'pb7' => '@m140918_173349_42137_c100716502550000001823136502281507_s1_p0',
);

# ==============
# ANALYSIS RULES
# ==============

my %roiFastas;
foreach my $roiId (keys(%rois)) {
    my $roi = "$resultsDir/roi.$roiId.fasta";
    my $roiCmd = "$indiana8 ExtractSequenceFeature -r $refs{'il'} -g $rois{$roiId} -t $roiTypes{$roiId} -o $roi";
    $dm->addRule($roi, $rois{$roiId}, $roiCmd);

    $roiFastas{$roiId} = $roi;
}

foreach my $id (keys(%pbrgs)) {
    my $rg = $pbrgs{$id};

    my $rgsubset = "$resultsDir/$id/$id.fastq";
    my $rgsubsetCmd = "grep --no-group-separator -A3 $rg $seqs{'pb'} > $rgsubset";
    $dm->addRule($rgsubset, $seqs{'pb'}, $rgsubsetCmd, 'nopostfix' => 1);

    $bams{$id} = undef;
    $seqs{$id} = $rgsubset;
}

foreach my $id (keys(%bams)) {
    my $bam = $bams{$id};
    my $fq = $seqs{$id};

    my $outdir = "$resultsDir/$id";
    
    if (defined($bam)) {
        my $coverage = "$outdir/coverage.txt";
        my $coverageCmd = "$gatk8 -T DepthOfCoverage -R $refs{$id} -I $bam -omitSampleSummary -omitIntervals -omitLocusTable -mmq 1 -o $coverage";
        $dm->addRule($coverage, $bam, $coverageCmd);

        my $coverageSimple = "$outdir/coverage.simple.txt";
        my $coverageSimpleCmd = "grep -v Locus $coverage | sed 's/:/\\t/g' | cut -f1-3 > $coverageSimple";
        $dm->addRule($coverageSimple, $coverage, $coverageSimpleCmd);

        my $errorsPerPosition = "$outdir/errors.per_position.txt";
        my $errorsPerContig = "$outdir/errors.per_contig.txt";
        my $errorsCmd = "$indiana8 FoldedContigErrorStats -c $bam -o $errorsPerPosition -o2 $errorsPerContig";
        $dm->addRule($errorsPerPosition, $bam, $errorsCmd);

        if ($id eq 'pb') {
            my $rgCounts = "$outdir/rgcounts.txt";
            my $rgCountsCmd = "samtools view $bam | awk '{ print \$14 }' | sort | uniq -c | awk '{ print \$2, \$1 }' > $rgCounts";
            $dm->addRule($rgCounts, $bam, $rgCountsCmd);
        }
    }

    if (defined($fq)) {
        my $lengthHist = "$outdir/lengthHist.txt";
        my $lengthStats = "$outdir/lengthStats.txt";
        my $lengthHistCmd = "$indiana8 lengthdist -f $id:$fq -b 100 -o $lengthHist -so $lengthStats";
        $dm->addRule($lengthHist, $fq, $lengthHistCmd);
    }
}

if (0) {
foreach my $asmid (keys(%asms)) {
    my $contigs = $asms{$asmid};
    my $outdir = "$resultsDir/$asmid";

    my $readsSam = "$outdir/$asmid.il_reads.sam";
    my $readsSamCmd = "bwa mem -t 30 $contigs $fastqs{'il'}->{'end1'} $fastqs{'il'}->{'end2'} > $readsSam";
    $dm->addRule($readsSam, $contigs, $readsSamCmd, 'nopostfix' => 1);

    my $readsBam = "$outdir/$asmid.il_reads.bam";
    my $readsBamCmd = "java -Xmx8g -jar ~/repositories/Picard-Latest/dist/SortSam.jar I=$readsSam O=$readsBam SO=coordinate CREATE_INDEX=true";
    $dm->addRule($readsBam, $readsSam, $readsBamCmd);

    my $readsWithRGsBam = "$outdir/$asmid.il_reads.with_rgs.bam";
    my $readsWithRGsBamCmd = "java -Xmx8g -jar ~/repositories/Picard-Latest/dist/AddOrReplaceReadGroups.jar I=$readsBam O=$readsWithRGsBam RGSM=PG0051-C RGID=ERR019061 RGPL=Illumina RGCN=WTSI RGLB=unknown RGPU=unknown RGDS=unknown CREATE_INDEX=true";
    $dm->addRule($readsWithRGsBam, $readsBam, $readsWithRGsBamCmd);

    my $dedupedBam = "$outdir/$asmid.il_reads.with_rgs.deduped.bam";
    my $dedupedMetrics = "$outdir/$asmid.il_reads.with_rgs.deduped.metrics";
    my $dedupedBamCmd = "java -Xmx8g -jar ~/repositories/Picard-Latest/dist/MarkDuplicates.jar I=$readsWithRGsBam O=$dedupedBam M=$dedupedMetrics CREATE_INDEX=true";
    $dm->addRule($dedupedBam, $readsWithRGsBam, $dedupedBamCmd);

    my $rtcTargets = "$outdir/$asmid.il_reads.with_rgs.deduped.rtc.intervals";
    my $rtcTargetsCmd = "$gatk8 -T RealignerTargetCreator -R $contigs -I $dedupedBam -o $rtcTargets";
    $dm->addRule($rtcTargets, $dedupedBam, $rtcTargetsCmd);

    my $realignedBam = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.bam";
    my $realignedBamCmd = "$gatk8 -T IndelRealigner -R $contigs -I $dedupedBam -targetIntervals $rtcTargets -o $realignedBam";
    $dm->addRule($realignedBam, $rtcTargets, $realignedBamCmd);

    my $recalTable = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal_data.table";
    my $recalTableCmd = "$gatk8 -T BaseRecalibrator -R $contigs -I $realignedBam -knownSites $resourcesDir/test.gatkfakeout.vcf -nct 10 -o $recalTable";
    $dm->addRule($recalTable, $realignedBam, $recalTableCmd);

    my $recalBam = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal.bam";
    my $recalBamCmd = "$gatk8 -T PrintReads -R $contigs -I $realignedBam --BQSR $recalTable -o $recalBam -nct 30";
    $dm->addRule($recalBam, $recalTable, $recalBamCmd);

    my $annotations = join(" -A ", "AlleleBalance", "AlleleBalanceBySample", "BaseCounts", "BaseQualityRankSumTest", "ChromosomeCounts", "ClippingRankSumTest", "Coverage", "DepthPerAlleleBySample", "DepthPerSampleHC", "FisherStrand", "GCContent", "HaplotypeScore", "HardyWeinberg", "HomopolymerRun", "LikelihoodRankSumTest", "LowMQ", "MappingQualityRankSumTest", "MappingQualityZero", "MappingQualityZeroBySample", "NBaseCount", "QualByDepth", "RMSMappingQuality", "ReadPosRankSumTest", "SpanningDeletions", "StrandBiasBySample", "StrandOddsRatio", "TandemRepeatAnnotator", "VariantType" );

    my $snps = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal.snp.vcf";
    my $snpsCmd = "$gatk8 -T UnifiedGenotyper -R $contigs -I $recalBam -o $snps -nt 30 --sample_ploidy 1 -glm SNP -A $annotations";
    $dm->addRule($snps, $recalBam, $snpsCmd);

    my $snpsFiltered = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal.snp.filtered.vcf";
    my $snpsFilteredCmd = "$gatk8 -T VariantFiltration -R $contigs -V $snps --filterExpression 'QD < 2.0 || MQ < 60.0 || FS > 40.0 || MQRankSum < -7.5 || ReadPosRankSum < -4.0 || MQ0 > 0 || HRun > 3' --filterName SNPFilter -o $snpsFiltered";
    $dm->addRule($snpsFiltered, $snps, $snpsFilteredCmd);

    my $indels = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal.indel.vcf";
    my $indelsCmd = "$gatk8 -T UnifiedGenotyper -R $contigs -I $recalBam -o $indels -nt 30 --sample_ploidy 1 -glm INDEL -A $annotations";
    $dm->addRule($indels, $recalBam, $indelsCmd);

    my $indelsFiltered = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal.indel.filtered.vcf";
    my $indelsFilteredCmd = "$gatk8 -T VariantFiltration -R $contigs -V $indels --filterExpression 'QD < 2.0 || ReadPosRankSum < -8.0 || FS > 200.0' --filterName IndelFilter -o $indelsFiltered";
    $dm->addRule($indelsFiltered, $indels, $indelsFilteredCmd);

    my $variantsFiltered = "$outdir/$asmid.il_reads.with_rgs.deduped.realigned.recal.variants.filtered.vcf";
    my $variantsFilteredCmd = "$gatk8 -T CombineVariants -R $contigs -V $snpsFiltered -V $indelsFiltered -o $variantsFiltered";
    $dm->addRule($variantsFiltered, [$snpsFiltered, $indelsFiltered], $variantsFilteredCmd);

    my $altRef = "$outdir/$asmid.GATKPolished.fasta";
    my $altRefCmd = "$gatk8 -T FastaAlternateReferenceMaker -R $contigs -V $variantsFiltered -o $altRef";
    $dm->addRule($altRef, $variantsFiltered, $altRefCmd);

    #my $altRefDict = "$outdir/$asmid.GATKPolished.dict";
    #my $altRefDictCmd = "java -Xmx8g -jar ~/repositories/Picard-Latest/dist/CreateSequenceDictionary.jar R=$altRef O=$altRefDict";
    #$dm->addRule($altRefDict, $altRef, $altRefDictCmd);

    my $altRefFai = "$outdir/$asmid.GATKPolished.fasta.fai";
    my $altRefFaiCmd = "samtools faidx $altRef";
    $dm->addRule($altRefFai, $altRef, $altRefFaiCmd);

    my $altRefBwt = "$outdir/$asmid.GATKPolished.fasta.bwt";
    my $altRefBwtCmd = "bwa index $altRef";
    $dm->addRule($altRefBwt, $altRef, $altRefBwtCmd);

    $asms{"$asmid.GATKPolished"} = $altRef;

    #data/ASMTest1.polished_assembly.fasta

    my $icornPolished = "$outdir/ICORN2.$asmid.polished_assembly.fasta.5";
    my $icornPolishedCmd = "cd $outdir && $icorn2 $dataDir/PfCross/ERR019061 250 $contigs 1 5";
    #$dm->addRule($icornPolished, $contigs, $icornPolishedCmd);
}
}

foreach my $asmid (keys(%asms)) {
    my $contigs = $asms{$asmid};
    my $outdir = "$resultsDir/$asmid";

    my @tags;
    push(@tags, "$asmid:$contigs");
    foreach my $id (keys(%existingAsms)) {
        push(@tags, "$id:$existingAsms{$id}");
    }

    my $bas = "$outdir/assemblies.stats";
    my $basCmd = "$indiana8 BasicAssemblyStats -c " . join(" -c ", @tags) . " -o $bas";
    $dm->addRule($bas, $contigs, $basCmd);

    my $delta = "$outdir/$asmid.delta";
    my $deltaCmd = "$nucmer -p $outdir/$asmid $refs{'il'} $contigs";
    $dm->addRule($delta, $contigs, $deltaCmd);

    my $filter = "$outdir/$asmid.filter";
    my $filterCmd = "$deltafilter -q $delta > $filter";
    $dm->addRule($filter, $delta, $filterCmd);

    my $coords = "$outdir/$asmid.filter.coords";
    my $coordsCmd = "$showcoords -rcl $filter > $coords";
    $dm->addRule($coords, $filter, $coordsCmd);

    my $dotplot = "$outdir/$asmid.filter.ps";
    my $dotplotCmd = "$mummerplot $filter -R $refs{'il'} -Q $contigs --filter --layout --postscript --fat -p $outdir/$asmid.filter";
    $dm->addRule($dotplot, $filter, $dotplotCmd);

    my $coordSummary = "$outdir/$asmid.filter.filter.coord_summary";
    my $coordSummaryCmd = "$showcoords -b -r -l -c -T -H $outdir/$asmid.filter.filter > $coordSummary";
    $dm->addRule($coordSummary, $dotplot, $coordSummaryCmd);

    my $variantSummary = "$outdir/$asmid.filter.filter.variant_summary";
    my $variantSummaryCmd = "$showsnps -l -r -T -H $outdir/$asmid.filter.filter > $variantSummary";
    $dm->addRule($variantSummary, $dotplot, $variantSummaryCmd);

    if (0) {
    foreach my $roiId (keys(%roiFastas)) {
        my $roiSam = "$outdir/$asmid.$roiId.sam";
        my $roiSamCmd = "bwa mem $contigs $roiFastas{$roiId} | sed 's/\|quiver//g' > $roiSam";
        $dm->addRule($roiSam, $contigs, $roiSamCmd, 'nopostfix' => 1);

        my $roiBam = "$outdir/$asmid.$roiId.bam";
        my $roiBamCmd = "java -Xmx8g -jar ~/repositories/Picard-Latest/dist/SortSam.jar I=$roiSam O=$roiBam SO=coordinate CREATE_INDEX=true";
        $dm->addRule($roiBam, $roiSam, $roiBamCmd);

        my $roiTable = "$outdir/$asmid.$roiId.table";
        my $roiTableCmd = "grep -v '\@' $roiSam | cut -f1-9,12-14 | awk '{ if (NF == 12) print \$0 }' | column -t > $roiTable";
        $dm->addRule($roiTable, $roiSam, $roiTableCmd);
    }
    }
}

$dm->execute();
