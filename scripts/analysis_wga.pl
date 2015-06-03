#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use lib "$FindBin::Bin/lib";

use ParseArgs;
use DM;
use Data::Dumper;

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
my $indiana16 = "java -Xmx16g -jar $binDir/indiana.jar";
my $sortsam = "java -Xmx64g -jar $binDir/SortSam.jar";
my $bwa = "$binDir/bwa";
my $gatk8 = "java -Xmx8g -jar $binDir/GenomeAnalysisTK.jar";
my $nucmer = "~/opt/MUMmer3.23/nucmer";
my $showaligns = "~/opt/MUMmer3.23/show-aligns";
my $showcoords = "~/opt/MUMmer3.23/show-coords";
my $showdiff = "~/opt/MUMmer3.23/show-diff";
my $showsnps = "~/opt/MUMmer3.23/show-snps";
my $showtiling = "~/opt/MUMmer3.23/show-tiling";
my $deltafilter = "~/opt/MUMmer3.23/delta-filter";
my $mummerplot = "~/opt/MUMmer3.23/mummerplot";
my $icorn2 = "~/opt/icorn2-v0.95//icorn2.sh";

my $ref = "$resourcesDir/references/3D7/PlasmoDB-9.0_Pfalciparum3D7_Genome.sorted.fasta";

my %fqs = (
    'unamplified' => 'data/ASMTest1.filtered_subreads.fastq',
    'amplified'   => 'data/WGATest5.filtered_subreads.fastq',
);

my %asms = (
    'unamplified' => 'data/ASMTest1.polished_assembly.fasta',
    'amplified'   => 'data/WGATest5.polished_assembly.fasta',
);

# ==============
# ANALYSIS RULES
# ==============

my $asmStats = "$resultsDir/assembly.stats";
my $asmStatsCmd = "$indiana8 BasicAssemblyStats -c unamplified:$asms{'unamplified'} -c amplified:$asms{'amplified'} -o $asmStats";
$dm->addRule($asmStats, [$asms{'unamplified'}, $asms{'amplified'}], $asmStatsCmd);

my %alignedUnamp = align($dm, 'seq' => $fqs{'unamplified'}, 'sample' => '3D7', 'readgroup' => 'unamplified', 'resultsDir' => "$resultsDir/reads", 't' => 50);
my %mummerUnamp  = mummer($dm, 'seq' => $asms{'unamplified'}, 'sample' => '3D7', 'readgroup' => 'unamplified', 'resultsDir' => "$resultsDir/mummer");

my %alignedAmp = align($dm, 'seq' => $fqs{'amplified'}, 'sample' => '3D7', 'readgroup' => 'amplified', 'resultsDir' => "$resultsDir/reads", 't' => 50);
my %mummerAmp  = mummer($dm, 'seq' => $asms{'amplified'}, 'sample' => '3D7', 'readgroup' => 'amplified', 'resultsDir' => "$resultsDir/mummer");

my $coverage = "$resultsDir/coverage.txt";
my $coverageCmd = "$gatk8 -T DepthOfCoverage -R $ref -I $alignedUnamp{'bam'} -I $alignedAmp{'bam'} -pt readgroup -omitSampleSummary -omitIntervals -omitLocusTable -mmq 1 -nt 10 -o $coverage";
$dm->addRule($coverage, [$alignedUnamp{'bai'}, $alignedAmp{'bai'}], $coverageCmd);

my $coverageSimple = "$resultsDir/coverage.simple.txt";
#my $coverageSimpleCmd = "grep -v Locus $coverage | sed 's/:/\\t/g' | cut -f1-4 > $coverageSimple";
my $coverageSimpleCmd = "grep -v Locus $coverage | sed 's/:/\\t/g' | cut -f1,2,5,6 > $coverageSimple";
$dm->addRule($coverageSimple, $coverage, $coverageSimpleCmd);

$dm->execute();

sub align {
    my ($dm, %o) = @_;

    my %a = (
        't' => 4,
        %o
    );

    my $bamPrefix = "$a{'resultsDir'}/$a{'sample'}.$a{'readgroup'}";
    my $bam = "$bamPrefix.bam";
    my $bamCmd = "$bwa mem -v 1 -t $a{'t'} -R \"\@RG\\tID:$a{'readgroup'}\\tSM:$a{'sample'}\\tPL:PACBIO\" $ref $a{'seq'} | samtools view -bS - | samtools sort - $bamPrefix";
    $dm->addRule($bam, [$ref, $a{'seq'}], $bamCmd, 'nopostfix' => 1);

    my $bai = "$bam.bai";
    my $baiCmd = "samtools index $bam";
    $dm->addRule($bai, $bam, $baiCmd);

    return ('bam' => $bam, 'bai' => $bai);
}

sub mummer {
    my ($dm, %o) = @_;

    my %a = (
        %o
    );

    my $asmid = "$a{'sample'}.$a{'readgroup'}";

    my $contigs = $a{'seq'};

    my $delta = "$a{'resultsDir'}/$asmid.delta";
    my $deltaCmd = "$nucmer -p $a{'resultsDir'}/$asmid $ref $contigs";
    $dm->addRule($delta, $contigs, $deltaCmd);

    my $filter = "$a{'resultsDir'}/$asmid.filter";
    my $filterCmd = "$deltafilter -q $delta > $filter";
    $dm->addRule($filter, $delta, $filterCmd);

    my $coords = "$a{'resultsDir'}/$asmid.filter.coords";
    my $coordsCmd = "$showcoords -rcl $filter > $coords";
    $dm->addRule($coords, $filter, $coordsCmd);

    my $dotplot = "$a{'resultsDir'}/$asmid.filter.ps";
    my $dotplotCmd = "$mummerplot $filter -R $ref -Q $contigs --filter --layout --postscript --fat -p $a{'resultsDir'}/$asmid.filter";
    $dm->addRule($dotplot, $filter, $dotplotCmd);

    my $coordSummary = "$a{'resultsDir'}/$asmid.filter.filter.coord_summary";
    my $coordSummaryCmd = "$showcoords -b -r -l -c -T -H $a{'resultsDir'}/$asmid.filter.filter > $coordSummary";
    $dm->addRule($coordSummary, $dotplot, $coordSummaryCmd);

    my $variantSummary = "$a{'resultsDir'}/$asmid.filter.filter.variant_summary";
    my $variantSummaryCmd = "$showsnps -l -r -T -H $a{'resultsDir'}/$asmid.filter.filter > $variantSummary";
    $dm->addRule($variantSummary, $dotplot, $variantSummaryCmd);

    return ('delta' => $delta);
}
