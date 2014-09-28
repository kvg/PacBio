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
my $nucmer = "~/opt/MUMmer3.23/nucmer";
my $showaligns = "~/opt/MUMmer3.23/show-aligns";
my $showcoords = "~/opt/MUMmer3.23/show-coords";
my $showdiff = "~/opt/MUMmer3.23/show-diff";
my $showsnps = "~/opt/MUMmer3.23/show-snps";
my $showtiling = "~/opt/MUMmer3.23/show-tiling";
my $deltafilter = "~/opt/MUMmer3.23/delta-filter";
my $mummerplot = "~/opt/MUMmer3.23/mummerplot";

my %bams = (
    'pb' => "$dataDir/Pfalc3D7/aligned_reads.bam",
    'il' => "$dataDir/PfCross/PG0051-C.ERR019061.bam",
);

my %seqs = (
    'pb' => "$dataDir/Pfalc3D7/filtered_subreads.fastq",
    'il' => "$dataDir/PfCross/ERR019061_1.fastq.gz",
);

my %asms = (
    'AsmTest1' => "data/ASMTest1.polished_assembly.fasta",
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
    "IT" => "/home/kiran/ngs/references/plasmodium/falciparum/IT/PlasmoDB-9.0/PlasmoDB-9.0_PfalciparumIT_Genome.fasta",
    "K1" => "/home/kiran/ngs/references/plasmodium/falciparum/K1/BroadInstitute/plasmodium_falciparum__isolate_k1__1_supercontigs.fasta",
    "PFCLIN" => "/home/kiran/ngs/references/plasmodium/falciparum/PFCLIN/SangerInstitute/PFCLIN.20080302.contigs.fasta",
    "RAJ116" => "/home/kiran/ngs/references/plasmodium/falciparum/RAJ116/BroadInstitute/plasmodium_falciparum_raj116_nucleus_1_supercontigs.fasta",
    "RO-33" => "/home/kiran/ngs/references/plasmodium/falciparum/RO-33/BroadInstitute/plasmodium_falciparum__isolate_ro-33__1_supercontigs.fasta",
    "SL" => "/home/kiran/ngs/references/plasmodium/falciparum/SL/BroadInstitute/plasmodium_falciparum__isolate_santa_lucia__1_supercontigs.fasta",
    "V34.04" => "/home/kiran/ngs/references/plasmodium/falciparum/Senegal_V34.04/BroadInstitute/plasmodium_falciparum__isolate_senegal_v34.04__1_supercontigs.fasta",
    "VS.1" => "/home/kiran/ngs/references/plasmodium/falciparum/VS.1/BroadInstitute/1__1_supercontigs.fasta",
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
}

$dm->execute();
