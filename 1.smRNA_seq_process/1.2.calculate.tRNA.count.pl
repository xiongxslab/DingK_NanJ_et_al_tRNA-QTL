#! /usr/bin/perl -w
use strict;

open IN,"samtools view -F 20 $ARGV[0] |" or die $!;
open OUT,">$ARGV[1]" or die $!;

my %count;
while(<IN>){
	chomp;
	my @sp = split /\t/;
	$sp[1] = 0;
	$count{$sp[2]} += 1;
}
close IN;

print OUT "tRNA\tcount\n";
foreach my $tRNA(sort keys %count){
	print OUT "$tRNA\t$count{$tRNA}\n";
}



