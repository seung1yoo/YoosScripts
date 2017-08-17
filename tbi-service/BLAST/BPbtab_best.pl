#!/usr/bin/perl -w

use strict;

open (FH, $ARGV[0]);

my %best;
my %best_bitscore;
while(<FH>) {
	chomp;
	if ($_ =~ /^#/) { print $_."\n"; next; }
	my @t_line = split (/\t/, $_);
#	my $index = "$t_line[0]_$t_line[5]";
	my $index = "$t_line[0]";
	if ($best{$index}) {
		if ($best_bitscore{$index} < $t_line[10]) { $best{$index} = $_; $best_bitscore{$index} = $t_line[10]; }
#		elsif ($best_bitscore{$index} == $t_line[10]) { $best{$index} .= "\n".$_; }
		else { next; }
	}
	else { $best{$index} = $_; $best_bitscore{$index} = $t_line[10]; }
}

foreach my $key(keys %best) {
	print "$best{$key}\n";
}
