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

my %bams = (
    'pb' => "$dataDir/Pfalc3D7/aligned_reads.bam",
    'il' => "$dataDir/PfCross/PG0051-C.ERR019061.bam",
);

my %seqs = (
    'pb' => "$dataDir/Pfalc3D7/filtered_subreads.fastq",
    'il' => "$dataDir/PfCross/ERR019061_1.fastq.gz",
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
        my $coverageHist = "$outdir/coverageHist.txt";
        my $coverageHistCmd = "$indiana8 CoverageHistogram -b $bam -o $coverageHist";
        $dm->addRule($coverageHist, $bam, $coverageHistCmd);
    }

    if (defined($fq)) {
        my $lengthHist = "$outdir/lengthHist.txt";
        my $lengthStats = "$outdir/lengthStats.txt";
        my $lengthHistCmd = "$indiana8 lengthdist -f $id:$fq -b 100 -o $lengthHist -so $lengthStats";
        $dm->addRule($lengthHist, $fq, $lengthHistCmd);
    }
}

$dm->execute();
