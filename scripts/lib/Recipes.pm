package Recipes;

use strict;
use Cwd;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK); # %EXPORT_TAGS);
use POSIX qw(ceil);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = qw(addMcCortexRules addSgaRules addAlignmentRules loadMetadata);
@EXPORT_OK = qw(addMcCortexRules addSgaRules addAlignmentRules loadMetadata);

sub roundoff {
    my $num = shift;

    my $roundto = 1;
    if ($num >= 10 && $num <= 100) { $roundto = 10; }
    if ($num >= 100 && $num <= 1000) { $roundto = 100; }

    return int(ceil($num/$roundto))*$roundto;
}

sub computeMla {
    my ($end1, $end2) = @_;

    my $mlb = int( ((-e $end1 ? -s $end1 : 0) + (-e $end2 ? -s $end2 : 0)) / (1024 * 1024 * 1024) );
    my $mla = roundoff($mlb);
    if ($end1 =~ /.gz/ || $end2 =~ /.gz/) { $mla *= 3; }
    if ($mla <= 3) { $mla = 5; }

    return $mla;
}

sub addMcCortexRules {
    my ($dm, %o) = @_;

    my %a = (
        'k'         => 47,
        'fq-cutoff' => 5,
        't'         => 1,
        'l'         => 0,
        'L'         => 400,
        'end2'      => '',
        'mla'       => undef,
        'ctxbin'    => 'bin/mccortex63',
        'indiana'   => 'bin/indiana.jar',
        'build'     => 1,
        'clean'     => 1,
        'infer'     => 1,
        'thread'    => 1,
        'contigs'   => 1,
        %o,
    );

    foreach my $arg (qw/end1 resultsDir sample/) {
        die "need to define $arg" unless defined $a{$arg};
    }

    my $mla = 0;
    my @ends1;
    my @ends2;
    my $seqCode;
    if (ref($a{'end1'}) eq 'ARRAY') {
        push(@ends1, @{$a{'end1'}});
        push(@ends2, @{$a{'end2'}});

        my @seqCodes;

        for (my $i = 0; $i <= $#ends1; $i++) {
            $mla += computeMla($ends1[$i], $ends2[$i]);

            push(@seqCodes, "--seq2 $ends1[$i]:$ends2[$i]");
        }

        $seqCode = join(" ", @seqCodes);
    } else {
        push(@ends1, $a{'end1'});
        push(@ends2, $a{'end2'});

        $mla = defined($a{'mla'}) ? $a{'mla'} : computeMla($a{'end1'}, $a{'end2'});

        $seqCode = (-e $a{'end1'} && ! -e $a{'end2'}) ? "--seqi $a{'end1'}" : "--seq2 $a{'end1'}:$a{'end2'}";
    }

    # Build binary
    my $ctx = "$a{'resultsDir'}/$a{'sample'}.ctx";
    my $ctxCmd = "$a{'ctxbin'} build -f -m ${mla}G -t $a{'t'} -k $a{'k'} --sample $a{'sample'} --fq-cutoff $a{'fq-cutoff'} $seqCode $ctx";
    if ($a{'build'}) { $dm->addRule($ctx, $ends1[0], $ctxCmd, 'outputFile' => "$ctx.log", 'memRequest' => $mla); }

    # Sort binary
    my $ctxSorted = "$a{'resultsDir'}/$a{'sample'}.sorted.ctx";
    my $ctxSortedCmd = "$a{'ctxbin'} sort -f -m ${mla}G -o $ctxSorted $ctx";
    #if ($a{'build'}) { $dm->addRule($ctxSorted, $ctx, $ctxSortedCmd); }

    # Clean binary
    my $clean = "$a{'resultsDir'}/$a{'sample'}.clean.ctx";
    my $covb = "$a{'resultsDir'}/$a{'sample'}.clean.covb.csv";
    my $cova = "$a{'resultsDir'}/$a{'sample'}.clean.cova.csv";
    my $lenb = "$a{'resultsDir'}/$a{'sample'}.clean.lenb.csv";
    my $lena = "$a{'resultsDir'}/$a{'sample'}.clean.lena.csv";
    my $cleanCmd = "$a{'ctxbin'} clean -f -m ${mla}G --out $clean -B 2 --covg-before $covb --covg-after $cova --len-before $lenb --len-after $lena $ctx";
    if (-e $ends1[0] && ! -e $ends2[0]) {
        $clean = $ctx;
    } else {
        if ($a{'clean'}) { $dm->addRule($clean, $ctx, $cleanCmd, 'outputFile' => "$clean.log", 'memRequest' => $mla); }
    }

    # Infer edges
    my $infer = "$a{'resultsDir'}/$a{'sample'}.infer.ctx";
    my $inferCmd = "$a{'ctxbin'} inferedges -f -m ${mla}G --out $infer $clean";
    if ($a{'infer'}) { $dm->addRule($infer, $clean, $inferCmd, 'outputFile' => "$infer.log", 'memRequest' => $mla); }

    # Verify
    my $verify = "$a{'resultsDir'}/$a{'sample'}.verify.txt";
    my $verifyCmd = "./scripts/verifyCtxs.pl CTX=$ctx CLEAN=$clean INFER=$infer OUT=$verify";
    #if ($a{'infer'}) { $dm->addRule($verify, $infer, $verifyCmd); }

    # Sort binary
    my $inferSorted = "$a{'resultsDir'}/$a{'sample'}.infer.sorted.ctx";
    my $inferSortedCmd = "$a{'ctxbin'} sort -f -m ${mla}G -o $inferSorted $infer";
    #if ($a{'infer'}) { $dm->addRule($inferSorted, $infer, $inferSortedCmd); }

    # Thread (single-end)
    my $se = "$a{'resultsDir'}/$a{'sample'}.infer.se.ctp";
    my $seCmd = "$a{'ctxbin'} thread -f -m ${mla}G -t $a{'t'} $seqCode --out $se $infer";
    if ($a{'thread'}) { $dm->addRule($se, $infer, $seCmd, 'outputFile' => "$se.log", 'memRequest' => $mla); }

    # Thread (paired-end)
    my $pe = "$a{'resultsDir'}/$a{'sample'}.infer.pe.ctp";
    my $peCmd = "$a{'ctxbin'} thread -f -m ${mla}G -t $a{'t'} -p $se -l $a{'l'} -L $a{'L'} --one-way $seqCode --out $pe $infer";
    if (-e $a{'end1'} && ! -e $a{'end2'}) {
        $pe = $se;
    } else {
        if ($a{'thread'}) { $dm->addRule($pe, $se, $peCmd, 'outputFile' => "$pe.log", 'memRequest' => $mla); }
    }

    # Extract seeds
    my $seedKmers = "$a{'resultsDir'}/$a{'sample'}.seeds.fasta";
    my $seedKmersCmd = "$a{'ctxbin'} view --kmers $clean | awk '{ print \">kmer.\" NR \"\\n\" \$1 }' > $seedKmers";
    if ($a{'contigs'}) { $dm->addRule($seedKmers, $clean, $seedKmersCmd, 'outputFile' => "$seedKmers.log", 'memRequest' => 22); }

    # Extract contigs
    my $contigs = "$a{'resultsDir'}/$a{'sample'}.contigs.fasta";
    my $contigsCmd = "$a{'ctxbin'} contigs -f -m " . 4*${mla} . "G --seed $seedKmers --no-reseed --out $contigs -p 0:$pe $infer";
    if ($a{'contigs'}) { $dm->addRule($contigs, [$seedKmers, $pe, $infer], $contigsCmd, 'outputFile' => "$contigs.log", 'memRequest' => $mla); }

    # Remove redundant contigs
    my $nonRedundantContigs = "$a{'resultsDir'}/$a{'sample'}.non_empty.non_redundant.contigs.fasta";
    my $nonRedundantContigsCmd = "$a{'ctxbin'} rmsubstr -f -m " . 4*${mla} . "G -o $nonRedundantContigs $contigs";
    if ($a{'contigs'}) { $dm->addRule($nonRedundantContigs, $contigs, $nonRedundantContigsCmd, 'outputFile' => "$nonRedundantContigs.log", 'memRequest' => $mla); }

    # Compute stats
    my $asmStats = "$a{'resultsDir'}/$a{'sample'}.non_empty.non_redundant.contigs.stats";
    my $asmStatsCmd = "java -Xmx8g -jar $a{'indiana'} BasicAssemblyStats -c $a{'sample'}:$nonRedundantContigs -o $asmStats";
    if ($a{'contigs'}) { $dm->addRule($asmStats, $nonRedundantContigs, $asmStatsCmd); }

    return ( 'ctx' => $ctx, 'clean' => $clean, 'infer' => $infer, 'ctp' => $pe, 'contigs' => $nonRedundantContigs );
}

sub addSgaRules {
    my ($dm, %o) = @_;

    my %a = (
        'ol'         => 75,
        'mol'        => 55,
        'ck'         => 41,
        't'          => 4,
        'd'          => 4000000,
        'covFilter'  => 2,
        'R'          => 10,
        'minPairs'   => 5,
        'minLength'  => 200,
        'maxGapDiff' => 0,
        'sga'        => 'bin/sga',
        'indiana'    => 'bin/indiana.jar',
        %o,
    );

    foreach my $arg (qw/end1 end2 resultsDir sample/) {
        die "need to define $arg" unless defined $a{$arg};
    }

    my $cwd = getcwd();
    my $sga = "$cwd/$a{'sga'}";

    # ==============
    # ANALYSIS RULES
    # ==============

    # First, preprocess the data to remove ambiguous basecalls
    my $preprocess = "$a{'resultsDir'}/preprocessed.fastq";
    my $preprocessCmd = "$a{'sga'} preprocess --pe-mode 1 -o $preprocess $a{'end1'} $a{'end2'}";
    $dm->addRule($preprocess, "", $preprocessCmd, 'outputFile' => "$preprocess.log");

    # ----------------
    # Error correction
    # ----------------

    # Build the index that will be used for error correction.
    my $index = "$preprocess.sai";
    my $indexCmd = "$a{'sga'} index -a ropebwt -t $a{'t'} --no-reverse -p $preprocess $preprocess";
    $dm->addRule($index, $preprocess, $indexCmd, 'outputFile' => "$index.log");

    # Perform error correction with a 41-mer.  The k-mer cutoff parameter is learned automatically
    my $correct = "$a{'resultsDir'}/reads.k$a{'ck'}.fastq";
    my $correctCmd = "$a{'sga'} correct -k $a{'ck'} --discard --learn -t $a{'t'} -p $preprocess -o $correct $preprocess";
    $dm->addRule($correct, [$preprocess, $index], $correctCmd, 'outputFile' => "$correct.log");

    # ---------------
    # Contig assembly
    # ---------------

    # Index the corrected data.
    my $index2 = "$correct.sai";
    my $index2Cmd = "$a{'sga'} index -a ropebwt -t $a{'t'} -p $correct $correct";
    $dm->addRule($index2, $correct, $index2Cmd, 'outputFile' => "$index2.log");

    # Remove exact-match duplicates and reads with low-frequency k-mers
    my $filter = "$correct.filter.pass.fa";
    my $filterCmd = "$a{'sga'} filter -x $a{'covFilter'} -t $a{'t'} --homopolymer-check --low-complexity-check -p $correct $correct";
    $dm->addRule($filter, $index2, $filterCmd, 'outputFile' => "$filter.log");

    # Merge simple, unbranched chains of vertices
    my $merged = "$a{'resultsDir'}/merged.k$a{'ck'}.fa";
    my $mergedCmd = "$a{'sga'} fm-merge -m $a{'mol'} -t $a{'t'} -p $correct -o $merged $filter";
    $dm->addRule($merged, $filter, $mergedCmd, 'outputFile' => "$merged.log");

    # Build an index of the merged sequences
    my $index3 = "$a{'resultsDir'}/merged.k$a{'ck'}.fa.sai";
    my $index3Cmd = "$a{'sga'} index -d 1000000 -t $a{'t'} -p $merged $merged";
    $dm->addRule($index3, $merged, $index3Cmd, 'outputFile' => "$index3.log");

    # Remove any substrings that were generated from the merge process
    my $rmdup = "$a{'resultsDir'}/merged.k$a{'ck'}.rmdup.fa";
    my $rmdupCmd = "$a{'sga'} rmdup -t $a{'t'} -o $rmdup -p $merged $merged";
    $dm->addRule($rmdup, [$merged, $index3], $rmdupCmd, 'outputFile' => "$rmdup.log");

    # Compute the structure of the string graph
    my $overlap = "$a{'resultsDir'}/merged.k$a{'ck'}.rmdup.asqg.gz";
    my $overlapCmd = "cd $a{'resultsDir'} && $cwd/$a{'sga'} overlap -m $a{'mol'} -t $a{'t'} $rmdup";
    $dm->addRule($overlap, $rmdup, $overlapCmd, 'outputFile' => "$overlap.log");

    # Perform the contig assembly without bubble popping
    my $assemble = "$a{'resultsDir'}/assemble.m$a{'ol'}";
    my $assembleFasta = "$a{'resultsDir'}/assemble.m$a{'ol'}-contigs.fa";
    my $assembleGraph = "$a{'resultsDir'}/assemble.m$a{'ol'}-graph.asqg.gz";
    my $assembleCmd = "$a{'sga'} assemble -m $a{'ol'} -g $a{'maxGapDiff'} -r $a{'R'} -o $assemble $overlap";
    $dm->addRule($assembleFasta, $overlap, $assembleCmd, 'outputFile' => "$assemble.log");

    # Compute stats
    my $asmStats = "$assembleFasta.stats";
    my $asmStatsCmd = "java -Xmx8g -jar $a{'indiana'} BasicAssemblyStats -c $a{'sample'}:$assembleFasta -o $asmStats";
    $dm->addRule($asmStats, $assembleFasta, $asmStatsCmd);

    return ( 'contigs' => $assembleFasta );
}

sub addAlignmentRules {
    my ($dm, %o) = @_;

    my %a = (
        't' => 1,
        'bwa' => 'bin/bwa',
        'picard' => 'bin/picard.jar',
        %o,
    );

    foreach my $arg (qw/end1 resultsDir sample readgroup ref/) {
        die "need to define $arg" unless defined $a{$arg};
    }

    my $bamPrefix = "$a{'resultsDir'}/$a{'sample'}.$a{'readgroup'}";
    my $bam = "$bamPrefix.bam";
    my $bamCmd = "$a{'bwa'} mem -v 1 -t $a{'t'} -R \"\@RG\\tID:$a{'readgroup'}\\tSM:$a{'sample'}\\tPL:ILLUMINA\" $a{'ref'} $a{'end1'} $a{'end2'} | samtools view -bS - | samtools sort - $bamPrefix";
    $dm->addRule($bam, [$a{'ref'}, $a{'end1'}], $bamCmd, 'nopostfix' => 1);

    my $md = "$a{'resultsDir'}/$a{'sample'}.$a{'readgroup'}.md.bam";
    my $mdMetrics = "$a{'resultsDir'}/$a{'sample'}.$a{'readgroup'}.md.metrics";
    my $mdCmd = "java -Xmx8g -jar $a{'picard'} MarkDuplicates I=$bam O=$md METRICS_FILE=$mdMetrics READ_NAME_REGEX=null VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true";
    $dm->addRule($md, $bam, $mdCmd);

    return $md;
}

sub addTrustedKmerRules {
    my ($dm, %o) = @_;

    my %a = (
        'indiana' => 'bin/indiana.jar',
        'mccortex63' => 'bin/mccortex63',
        'blastn' => 'bin/blastn',
    );

    foreach my $arg (qw/joinCrossCtx/) {
        die "need to define $arg" unless defined $a{$arg};
    }

    my $indiana = "java -Xmx8g -jar $a{'indiana'}";
    my $joinCrossCtx = $a{'joinCrossCtx'};
    (my $basename = $joinCrossCtx) =~ s/.ctx//;

#    my $novelKmersCtx = "$basename.novel.ctx";
#    my $novelKmersCtxCmd = "$indiana EmitNovelKmers -g $joinCrossCtx -o $novelKmersCtx";
#    $dm->addRule($novelKmersCtx, $joinCrossCtx, $novelKmersCtxCmd);
#
#    my $novelKmersConfidentCtx = "$basename.novel.confident.ctx";
#    my $novelKmersConfidentCtxCmd = "$indiana EmitNovelKmers -g $novelKmersCtx -t $kmerCovThreshold -o $novelKmersConfidentCtx";
#    $dm->addRule($novelKmersConfidentCtx, [$novelKmersCtx, $kmerCovThreshold], $novelKmersConfidentCtxCmd);
#
#    my $novelKmersConfidentFasta = "$basename.novel.confident.fasta";
#    my $novelKmersConfidentFastaCmd = "$a{'mccortex63'} view --kmers $novelKmersConfidentCtx | awk '{ print \">\" NR \"\\n\" \$1 }' > $novelKmersConfidentFasta";
#    $dm->addRule($novelKmersConfidentFasta, $novelKmersConfidentCtx, $novelKmersConfidentFastaCmd);
#
#    my $novelKmersConfidentBlast = "$novelKmersConfidentFasta.blast";
#    my $novelKmersConfidentBlastCmd = "$a{'blastn'} -db /home/kiran/opt/rmblast-2.2.28/db/nt/nt -query $novelKmersConfidentFasta -num_alignments 1 -outfmt '6 qseqid sseqid length mismatch gaps evalue stitle' -out $novelKmersConfidentBlast";
#    $dm->addRule($novelKmersConfidentBlast, $novelKmersConfidentFasta, $novelKmersConfidentBlastCmd);
#
#    my $contaminantIds = "$novelKmersConfidentFasta.contaminantIds.txt";
#    my $contaminantIdsCmd = "grep -v 'Plasmodium' $novelKmersConfidentBlast | awk '{ print \">\" \$1 }' | sort | uniq > $contaminantIds";
#    $dm->addRule($contaminantIds, $novelKmersConfidentBlast, $contaminantIdsCmd);
#
#    my $trustedIds = "$novelKmersConfidentFasta.trustedIds.txt";
#    my $trustedIdsCmd = "grep '>' $novelKmersConfidentFasta | grep -v -x -f $contaminantIds > $trustedIds";
#    $dm->addRule($trustedIds, $contaminantIds, $trustedIdsCmd);
#
#    my $novelKmersTrusted = "$novelKmersConfidentFasta.trustedKmers.fasta";
#    my $novelKmersTrustedCmd = "$selectcontigs $trustedIds $novelKmersConfidentFasta > $novelKmersTrusted";
#    $dm->addRule($novelKmersTrusted, $trustedIds, $novelKmersTrustedCmd);
#
#    my $novelKmersTrustedList = "$novelKmersConfidentFasta.trustedKmers.list";
#    my $novelKmersTrustedListCmd = "grep -v '>' $novelKmersTrusted > $novelKmersTrustedList";
#    $dm->addRule($novelKmersTrustedList, $novelKmersTrusted, $novelKmersTrustedListCmd);
}

sub loadMetadata {
    my %o = (
        'allowCommentedOutSamples' => 0,
        'ignoreCross' => '',
        @_
    );

    # Load parent information
    my %parents;
    my %crossParents;

    open(PARENTS, "lists/parents.txt");
    chomp(my $pheader = <PARENTS>);
    my @pheader = split(/\t/, $pheader);

    while (my $line = <PARENTS>) {
        if ($line !~ /^#/) {
            my ($cross, $p1, $a1, $p2, $a2) = split(/\s+/, $line);

            next if $cross eq $o{'ignoreCross'};

            my ($c1, $c2) = split(/x/, $cross);

            $parents{$p1} = $c1;
            $parents{$p2} = $c2;

            $crossParents{$cross}->{'parent1'} = $p1;
            $crossParents{$cross}->{'parent2'} = $p2;
            $crossParents{$cross}->{'acc1'} = $a1;
            $crossParents{$cross}->{'acc2'} = $a2;
        }
    }
    close(PARENTS);

    # Load manifest file
    my %samples;

    open(MANIFEST, "lists/manifest.txt");
    chomp(my $header = <MANIFEST>);
    my @header = split(/\t/, $header);

    my $cwd = getcwd();

    my %fastas = (
        '3D7xHB3' => { '3D7' => '/home/kiran/ngs/references/plasmodium/falciparum/3D7/PlasmoDB-9.0/PlasmoDB-9.0_Pfalciparum3D7_Genome.fasta',
                       'HB3' => '/home/kiran/ngs/references/plasmodium/falciparum/HB3/BroadInstitute/plasmodium_falciparum__isolate_hb3__1_supercontigs.fasta' },
        'HB3xDD2' => { 'HB3' => '/home/kiran/ngs/references/plasmodium/falciparum/HB3/BroadInstitute/plasmodium_falciparum__isolate_hb3__1_supercontigs.fasta',
                       'DD2' => '/home/kiran/ngs/references/plasmodium/falciparum/Dd2/BroadInstitute/plasmodium_falciparum__isolate_dd2__1_supercontigs.fasta' },
        '7G8xGB4' => { '7G8' => '/home/kiran/ngs/references/plasmodium/falciparum/7G8/BroadInstitute/plasmodium_falciparum__isolate_7g8__1_supercontigs.fasta',
                       'GB4' => undef },
        '803xGB4' => { '803' => undef,
                       'GB4' => undef }
    );

    while (my $line = <MANIFEST>) {
        if ($line !~ /^#/ || $o{'allowCommentedOutSamples'}) {
            chomp($line);
            $line =~ s/^#//;

            my @line = split(/\t/, $line);
            my %entry;

            for (my $i = 0; $i < scalar(@header); $i++) {
                $entry{$header[$i]} = $line[$i];
            }

            my $cross = $entry{'cross'};

            next if $cross eq $o{'ignoreCross'};

            my ($c1, $c2) = split(/x/, $cross);
            my $id = "$entry{'sample'}.$entry{'accession'}";
            my $p1 = $crossParents{$cross}->{'parent1'};
            my $p2 = $crossParents{$cross}->{'parent2'};

            $samples{$cross}->{$id}->{'sample'} = $entry{'sample'};
            $samples{$cross}->{$id}->{'accession'} = $entry{'accession'};
            $samples{$cross}->{$id}->{'end1'} = "$cwd/$entry{'end1'}";
            $samples{$cross}->{$id}->{'end2'} = "$cwd/$entry{'end2'}";
            $samples{$cross}->{$id}->{'isParent'} = ($entry{'sample'} eq $p1 || $entry{'sample'} eq $p2) ? 1 : 0;
            $samples{$cross}->{$id}->{'parent1'} = $crossParents{$cross}->{'parent1'};
            $samples{$cross}->{$id}->{'parent2'} = $crossParents{$cross}->{'parent2'};
            $samples{$cross}->{$id}->{'acc1'} = $crossParents{$cross}->{'acc1'};
            $samples{$cross}->{$id}->{'acc2'} = $crossParents{$cross}->{'acc2'};
            $samples{$cross}->{$id}->{'fa1'} = $fastas{$cross}->{$c1};
            $samples{$cross}->{$id}->{'fa2'} = $fastas{$cross}->{$c2};
        }
    }
    close(MANIFEST);

    return %samples;
}

1;
