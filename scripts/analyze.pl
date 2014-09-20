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
    'KEEP_GOING' => 1
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

# ==============
# ANALYSIS RULES
# ==============

#$dm->addRule(target, dependencies, command, 'outputFile' => "$logsDir/target.log");

$dm->execute();
