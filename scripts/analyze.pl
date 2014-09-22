#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use lib "$FindBin::Bin/lib";

use ParseArgs;
use DM;

my %args = &getCommandArguments(
    'RUN_ID'     => undef,
    'DRY_RUN'    => 1,
    'NUM_JOBS'   => 1,
    'KEEP_GOING' => 1,
    'CLUSTER'    => 'localhost',
    'QUEUE'      => 'localhost',
);

my $cwd = getcwd();
my $projectName = basename($cwd);

my $binDir = "$cwd/bin";
my $dataDir = "$cwd/data";
my $listsDir = "$cwd/lists";
my $resultsDir = "$cwd/results/$args{'RUN_ID'}";
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
my $gatk8 = "java -Xmx8g -jar bin/GenomeAnalysisTK.v3.2.jar";

my %bams = (
    'pb' => "$dataDir/Pfalc3D7/aligned_reads.bam",
    'il' => "$dataDir/PfCross/PG0051-C.ERR019061.bam",
);

my %seqs = (
    'pb' => "$dataDir/Pfalc3D7/filtered_subreads.fastq",
    'il' => "$dataDir/PfCross/ERR019061_1.fastq.gz",
);

my %refs = (
    'pb' => "$resourcesDir/3D7.cshl.fasta",
    'il' => "/home/kiran/ngs/references/plasmodium/falciparum/3D7/PlasmoDB-9.0/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta",
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

$dm->execute();
