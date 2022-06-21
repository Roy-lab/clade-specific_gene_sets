#!/usr/bin/env perl
#usage : ./script.pl [config]
use strict;
use warnings;

my $target_anc_point = "";
my @object = ();
my @extant = ();
my @subanc = ();
my @ancestor = ();
my @negate = ();
my $gnlist_file = "";
my $matrix_file = "";

open (CONFIG, $ARGV[0]) or die  " - Check config file.\n\n";
while (<CONFIG>) {
        my $line = $_; chomp $line;
        my @line = split /=/, $line;
        if ($line =~ /^target_anc_point/) {
                $target_anc_point = $line[1];
        } elsif ($line =~ /^object/) {
                @object = split / /, $line[1];
        } elsif ($line =~ /^extant/) {
                @extant = split / /, $line[1];
        } elsif ($line =~ /^subancestor/) {
		if (exists $line[1] && $line[1] ne "") {
	                @subanc = split / /, $line[1];
		}
        } elsif ($line =~ /^ancestor/) {
                @ancestor = split / /, $line[1];
        } elsif ($line =~ /^negate/) {
                @negate = split / /, $line[1];
        } elsif ($line =~ /^gnlist/) {
                $gnlist_file = $line[1];
        } elsif ($line =~ /^matrix/) {
                $matrix_file = $line[1];
        }
}
close CONFIG;

open (GNLS, $gnlist_file) or die " - No gene list file: $gnlist_file\n\n";
my @gnls = <GNLS>; chomp @gnls;
my %gns = ();
foreach my$gn (@gnls) {$gns{$gn}++;}
close GNLS;

open (MAT, $matrix_file) or die " - No matrix file:$matrix_file\n\n";
while (<MAT>) {
	my $line = $_; chomp $line;
	my @line = split /\t/, $line;
	my $repid = shift @line;

	## Main rule
	my $target_id = $line[$target_anc_point];	
	if ($target_id == 0) {next;}
	my $more = $target_id + 2;	# * #
	my $cnt = 0;
	for (my$a=0; $a<scalar @ancestor; $a++) {
		my $comp_anc = $line[$ancestor[$a]];
        	if ($comp_anc == 0) {next;}
		if ($comp_anc >= $more) {$cnt++;}	# * #
	}
	unless ($cnt == scalar @ancestor) {next;}

	## Rule for sub-ancestor
        my $subancestor_violated=0;
        if (exists $subanc[0] && $subanc[0] ne "") {
                foreach my$sub (@subanc) {
                        if ($line[$sub] < $line[$target_anc_point]) { # * #
                                $subancestor_violated++;
                        }
                }
        }
        if ($subancestor_violated > 0) {next;}

	## Rule for negators
	$cnt = 0;
	my $cnt0 = 0;
	foreach my$i (@negate) {
		unless ($line[$i] == 0) {
			foreach my$a (@ancestor) {
				if ($line[$i] < $line[$a]) {$cnt++;}	# * #
			}
		} else {
			$cnt0++;
		}
	}
	unless ($cnt == 0) {next;}
	if ($cnt0 == scalar @negate) {next;}
	
	## Rule for object specific
	my $min = 0;
	foreach my$i (@object) {
		if ($min > $line[$i]) {$min = $line[$i];}
	}
	$cnt = 0;
	foreach my$i (@negate) {
		if ($line[$i] == 0) {next;}
		if ($line[$i] <= $min) {$cnt++;}	# * #
	}
	unless ($cnt == 0) {next;}


	print $line."\n";
}
close MAT;

exit
