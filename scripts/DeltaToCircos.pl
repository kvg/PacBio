#!/usr/bin/perl -w

use strict;

my ($delta) = @ARGV;

open(DELTA, $delta);
my $chrCurrent = undef;
while (my $line = <DELTA>) {
    if ($line =~ /^>/) {
        #>Pf3D7_14_v3 Supercontig_1.1 3291936 377975
        my ($chrRef, $chrAlt, $lengthRef, $lengthAlt) = split(/\s+/, $line);

        $chrCurrent = $chrRef;
        $chrCurrent =~ s/^>//g;
    } else {
        if (defined($chrCurrent)) {
            my @pieces = split(/\s+/, $line);
            if (scalar(@pieces) == 7) {
                my ($refStart, $refEnd, $altStart, $altEnd, $junk) = @pieces;

                print join(" ", $chrCurrent, $refStart, $refEnd) . "\n";
            }
        }
    }
}
close(DELTA);
