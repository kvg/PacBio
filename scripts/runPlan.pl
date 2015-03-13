#!/usr/bin/perl -w

use strict;

use Cwd;
use FindBin;
use File::Basename;
use lib "$FindBin::Bin/lib";

use ParseArgs;
use DM;

my %args = &getCommandArguments(
    'PLAN'       => undef,
    'DRY_RUN'    => 1,
    'NUM_JOBS'   => 1,
    'KEEP_GOING' => 1,
    'CLUSTER'    => 'localhost',
    'QUEUE'      => 'localhost',
);

my $dm = new DM(
    'dryRun'     => $args{'DRY_RUN'},
    'numJobs'    => $args{'NUM_JOBS'},
    'keepGoing'  => $args{'KEEP_GOING'},
    'cluster'    => $args{'CLUSTER'},
    'queue'      => $args{'QUEUE'},
    'outputFile' => "$args{'PLAN'}.dm.log",
);

open(PLAN, $args{'PLAN'});
chomp(my @plan = <PLAN>);
close(PLAN);

my @jobsDaligner;
my @jobsSort;
my @jobsLevel;

my $target = '';

foreach my $plan (@plan) {
    if    ($plan =~ /^# Daligner jobs/)     { $target = "daligner"; }
    elsif ($plan =~ /^# Initial sort jobs/) { $target = "sort";     }
    elsif ($plan =~ /^# Level 1 jobs/)      { $target = "level";    }
    else {
        if ($target eq 'daligner') { push(@jobsDaligner, $plan); }
        if ($target eq 'sort')     { push(@jobsSort, $plan);     }
        if ($target eq 'level')    { push(@jobsLevel, $plan);    }
    }
}

my @douts;
for (my $i = 0; $i <= $#jobsDaligner; $i++) {
    my $out = "daligner.$i";
    my $cmd = "-$jobsDaligner[$i] > /dev/null 2>&1 && touch $out";
    $dm->addRule($out, "", $cmd);

    push(@douts, $out);
}

my @souts;
for (my $i = 0; $i <= $#jobsSort; $i++) {
    my $out = "sort.$i";
    my $cmd = "$jobsSort[$i] && touch $out";
    $dm->addRule($out, \@douts, $cmd);

    push(@souts, $out);
}

my @levels;
for (my $i = 0; $i <= $#jobsLevel; $i++) {
    my $out = "level.$i";
    my $cmd = "$jobsLevel[$i] && touch $out";
    $dm->addRule($out, \@souts, $cmd);

    push(@levels, $out);
}

$dm->execute();
