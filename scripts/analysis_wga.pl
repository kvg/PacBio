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

# ==============
# ANALYSIS RULES
# ==============

my $asmStats = "$resultsDir/assembly.stats";
my $asmStatsCmd = "$indiana8 BasicAssemblyStats -c unamplified:$asms{'unamplified'} -c amplified:$asms{'amplified'} " . flattenAsmList() . " -r $ref -o $asmStats";
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

sub flattenAsmList {
    my @list;
    foreach my $key (keys(%existingAsms)) {
        push(@list, "-c $key:$existingAsms{$key}");
    }

    return join(" ", @list);
}

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
