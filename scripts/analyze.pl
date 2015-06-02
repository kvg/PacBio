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
    'RUN_ID'     => undef,
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
my $indiana16 = "java -Xmx16g -jar $binDir/indiana.jar";
my $sortsam = "java -Xmx64g -jar $binDir/SortSam.jar";
my $bwa = "$binDir/bwa";
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

my $ref = "/home/kiran/ngs/references/plasmodium/falciparum/3D7/PlasmoDB-9.0/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta";

my %fqs = (
    'unamplified' => {
        'unamplified' => 'data/Pfalc3D7/filtered_subreads.fastq'
    },
    'amplified'   => {
        'LID50115'        => 'data/LID50115/LID50115_filtered_subreads.fastq',
        'LID50131_032715' => 'data/LID50131/LID50131_032715_filtered_subreads.fastq',
        'LID50131_040715' => 'data/LID50131/LID50131_040715_filtered_subreads.fastq'
    }
);

print Dumper(\%fqs);

#my %bas = (
#    'unamplified' => [
#        'data/PfalcRaw/A01_1/Analysis_Results/m140912_185114_42137_c100689352550000001823145102281580_s1_p0.bas.h5',
#        'data/PfalcRaw/A01_2/Analysis_Results/m140912_220917_42137_c100689352550000001823145102281581_s1_p0.bas.h5',
#        'data/PfalcRaw/A01_3/Analysis_Results/m140913_012852_42137_c100689352550000001823145102281582_s1_p0.bas.h5',
#        'data/PfalcRaw/A01_4/Analysis_Results/m140913_044743_42137_c100689352550000001823145102281583_s1_p0.bas.h5',
#        'data/PfalcRaw/B01_1/Analysis_Results/m140918_073419_42137_c100716502550000001823136502281504_s1_p0.bas.h5',
#        'data/PfalcRaw/B01_2/Analysis_Results/m140918_105250_42137_c100716502550000001823136502281505_s1_p0.bas.h5',
#        'data/PfalcRaw/B01_3/Analysis_Results/m140918_141535_42137_c100716502550000001823136502281506_s1_p0.bas.h5',
#        'data/PfalcRaw/B01_4/Analysis_Results/m140918_173349_42137_c100716502550000001823136502281507_s1_p0.bas.h5',
#    ],
#    'amplified' => [
#        'data/LID50115/m150301_071648_42137_c100745092550000001823165807071517_s1_p0.bas.h5',
#        'data/LID50115/m150301_134425_42137_c100745502550000001823165807071550_s1_p0.bas.h5',
#        'data/LID50115/m150301_180216_42137_c100745502550000001823165807071551_s1_p0.bas.h5',
#        'data/LID50115/m150301_222135_42137_c100745502550000001823165807071552_s1_p0.bas.h5',
#        'data/LID50115/m150302_024043_42137_c100745502550000001823165807071553_s1_p0.bas.h5',
#        'data/LID50115/m150302_070004_42137_c100745502550000001823165807071554_s1_p0.bas.h5',
#        'data/LID50115/m150302_111906_42137_c100745502550000001823165807071555_s1_p0.bas.h5',
#        'data/LID50115/m150302_154108_42137_c100745502550000001823165807071556_s1_p0.bas.h5',
#    ]
#);

# ==============
# ANALYSIS RULES
# ==============

my $chrKmersK17 = "$resultsDir/kmers.k17.table";
my $chrKmersK17Cmd = "$indiana16 FindUniqueAndDistantKmers -r $ref -k 17 -o $chrKmersK17";
$dm->addRule($chrKmersK17, $ref, $chrKmersK17Cmd);

my $chrKmersK21 = "$resultsDir/kmers.k21.table";
my $chrKmersK21Cmd = "$indiana16 FindUniqueAndDistantKmers -r $ref -k 21 -o $chrKmersK21";
$dm->addRule($chrKmersK21, $ref, $chrKmersK21Cmd);

foreach my $fqkey (keys(%fqs)) {
    my %fqs = %{$fqs{$fqkey}};

    foreach my $key (keys(%fqs)) {
        my $fq = $fqs{$key};

        my $fasta = "$resultsDir/$key/$key.fasta";
        my $fastaCmd = "grep --no-group-separator -A1 '^\@m' $fq | sed 's/^\@m/>m/' > $fasta";
        $dm->addRule($fasta, $fq, $fastaCmd, 'nopostfix' => 1);

        my $asmStats = "$resultsDir/$key/$key.asmStats";
        my $asmStatsCmd = "$indiana8 AssemblyEval -f $fasta -o $asmStats";
        $dm->addRule($asmStats, $fasta, $asmStatsCmd);

        my $lengths = "$resultsDir/$key/$key.lengths.txt";
        my $lengthsCmd = "grep --no-group-separator -A1 '^\@m' $fq | grep -v '^\@m' | awk '{ print length(\$0) }' > $lengths";
        $dm->addRule($lengths, $fq, $lengthsCmd, 'nopostfix' => 1);

        my $chimeraTable = "$resultsDir/$key/$key.chimeras.k17.txt";
        my $chimeraTableCmd = "$indiana16 FindChimericLongReads -t $chrKmersK17 -l $fq -o $chimeraTable";
        $dm->addRule($chimeraTable, $chrKmersK17, $chimeraTableCmd);

        my $sam = "$resultsDir/$key/$key.sam";
        my $samCmd = "$bwa mem -x pacbio -t 10 -R '\@RG\\tID:$key\\tSM:3D7_WGA' $ref $fq > $sam";
        $dm->addRule($sam, $fq, $samCmd, 'nopostfix' => 1);

        my $bam = "$resultsDir/$key/$key.bam";
        my $bamCmd = "$sortsam I=$sam O=$bam SO=coordinate CREATE_INDEX=true";
        $dm->addRule($bam, $sam, $bamCmd);

        my $samTable = "$resultsDir/$key/$key.sam.table";
        my $samTableCmd = "grep -v '^\@' $sam | awk '{ print \$1, length(\$10), \$0 ~ \"SA:\", \$3, \$4, \$5 }' > $samTable";
        $dm->addRule($samTable, $sam, $samTableCmd, 'nopostfix' => 1);

        my $cov = "$resultsDir/$key/$key.cov";
        my $covCmd = "$gatk8 -T DepthOfCoverage -R $ref -I $bam -omitLocusTable -omitBaseOutput -pt sample -pt readgroup -o $cov";
        $dm->addRule($cov, $bam, $covCmd);
    }
}

$dm->execute();
