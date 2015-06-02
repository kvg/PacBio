#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use File::Basename;

sub printAnalyses {
    my %analyses = @_;

    print "\nAvailable analyses:\n";

    foreach my $analysis (sort keys(%analyses)) {
        print "  - $analysis [arguments]\n";
    }

    print "\n";
}

chomp(my @analyses = qx(find scripts/ -name \*analysis_\*.pl));
my %analyses;
foreach my $analysis (@analyses) {
    (my $name = basename($analysis)) =~ s/(analysis_|\.pl)//g;
    $analyses{$name} = $analysis;

}

if (scalar(@ARGV) == 0) {
    printAnalyses(%analyses);
} else {
    my $analysis = shift(@ARGV);

    if (!exists($analyses{$analysis})) {
        print "Analysis '$analysis' doesn't exist.\n\n";
        printAnalyses(%analyses);
    } else {
        my $cmd = "perl -w $analyses{$analysis} ANALYSIS=$analysis " . join(" ", @ARGV);
        system($cmd);
    }
}
