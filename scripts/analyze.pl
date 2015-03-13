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
    #'touch' => 1
);

my $indiana8 = "java -Xmx8g -jar $binDir/indiana.jar";
my $gatk8 = "java -Xmx8g -jar /data1/users/kiran/repositories/GATK-Protected/protected/gatk-package-distribution/target/gatk-package-distribution-3.2.jar";

my $daligner = "/home/kiran/repositories/DALIGNER/daligner";
my $HPCdaligner = "/home/kiran/repositories/DALIGNER/HPCdaligner";
my $HPCmapper = "/home/kiran/repositories/DALIGNER/HPCmapper";
my $LAsort = "/home/kiran/repositories/DALIGNER/LAsort";
my $LAmerge = "/home/kiran/repositories/DALIGNER/LAmerge";
my $LAsplit = "/home/kiran/repositories/DALIGNER/LAsplit";
my $LAcat = "/home/kiran/repositories/DALIGNER/LAcat";
my $LAshow = "/home/kiran/repositories/DALIGNER/LAshow";
my $LAcheck = "/home/kiran/repositories/DALIGNER/LAcheck";

my $fasta2DB = "/home/kiran/repositories/DAZZ_DB/fasta2DB";
my $DB2fasta = "/home/kiran/repositories/DAZZ_DB/DB2fasta";
my $quiva2DB = "/home/kiran/repositories/DAZZ_DB/quiva2DB";
my $DB2quiva = "/home/kiran/repositories/DAZZ_DB/DB2quiva";
my $DBsplit = "/home/kiran/repositories/DAZZ_DB/DBsplit";
my $DBdust = "/home/kiran/repositories/DAZZ_DB/DBdust";
my $Catrack = "/home/kiran/repositories/DAZZ_DB/Catrack";
my $DBshow = "/home/kiran/repositories/DAZZ_DB/DBshow";
my $DBstats = "/home/kiran/repositories/DAZZ_DB/DBstats";
my $DBrm = "/home/kiran/repositories/DAZZ_DB/DBrm";
my $simulator = "/home/kiran/repositories/DAZZ_DB/simulator";
my $simfromsource = "/home/kiran/repositories/DAZZ_DB/simfromsource";
my $fasta2DAM = "/home/kiran/repositories/DAZZ_DB/fasta2DAM";
my $DAM2fasta = "/home/kiran/repositories/DAZZ_DB/DAM2fasta";

my $dextract = "/home/kiran/repositories/DEXTRACTOR/dextract";
my $dexta = "/home/kiran/repositories/DEXTRACTOR/dexta";
my $undexta = "/home/kiran/repositories/DEXTRACTOR/undexta";
my $dexqv = "/home/kiran/repositories/DEXTRACTOR/dexqv";
my $undexqv = "/home/kiran/repositories/DEXTRACTOR/undexqv";

my $ref3D7 = "/home/kiran/ngs/references/plasmodium/falciparum/3D7/PlasmoDB-9.0/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta";

# ==============
# ANALYSIS RULES
# ==============

my @fastas;

#chomp(my @h5s = qx(find data/PfalcRaw/ -name \*.bax.h5));
#foreach my $h5 (@h5s) {
#    (my $basename = basename($h5)) =~ s/.bax.h5//g;
#
#    my $fasta = "$resultsDir/fastas/$basename.fasta";
#    my $fastaCmd = "$dextract -l500 -o$resultsDir/fastas/$basename $h5";
#    $dm->addRule($fasta, $h5, $fastaCmd);
#
#    push(@fastas, $fasta);
#}

my $genome = "$resultsDir/Pf3D7_01_v3.fasta";
my $genomeCmd = "~/bin/selectcontigs Pf3D7_01_v3 $ref3D7 > $genome";
$dm->addRule($genome, $ref3D7, $genomeCmd);

my $fasta = "$resultsDir/reads.fasta";
my $map = "$resultsDir/reads.map";
my $fastaCmd = "$simfromsource $genome -c5 -M$map -e0.0 > $fasta";
$dm->addRule($fasta, $genome, $fastaCmd);

push(@fastas, $fasta);

my $db = "$resultsDir/unamplified.3D7.db";
my $dbCmd = "$fasta2DB -v $db " . join(" ", @fastas);
$dm->addRule($db, \@fastas, $dbCmd);

my $dust = "$resultsDir/.unamplified.3D7.dust.data";
my $dustCmd = "$DBdust -b $db";
$dm->addRule($dust, $db, $dustCmd);

my $split = "$resultsDir/.unamplified.3D7.split";
my $splitCmd = "$DBsplit -s200 $db && touch $db $dust $split";
$dm->addRule($split, $dust, $splitCmd);

my $plan = "$resultsDir/daligner.plan";
my $planCmd = "$HPCdaligner -d -b -e0.95 -M16 $db > $plan";
$dm->addRule($plan, $split, $planCmd, 'nopostfix' => 1);

my $res = "$resultsDir/res";
my $resCmd = "./scripts/runPlan.pl PLAN=$plan DRY_RUN=0 NUM_JOBS=15 && touch $res";
$dm->addRule($res, $plan, $resCmd);

my $align = "$resultsDir/unamplified.3D7.1.las.txt";
my $alignCmd = "$LAshow $db unamplified.3D7.1.las > $align";
$dm->addRule($align, $res, $alignCmd, 'nopostfix' => 1);

$dm->execute();
