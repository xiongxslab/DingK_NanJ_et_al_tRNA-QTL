#! /usr/bin/perl -w
use strict;

my @cols = qw/pval_beta_fdr qval_beta/;
my @cutoffs = qw/0.05 0.1/;

# read in nominal p-value files
foreach my $col(@cols){
  foreach my $cutoff(@cutoffs){
    open IN,"../2.fastqtl_permute/03.permute.out.$col.$cutoff.txt" or die $!;
    my %cut;
    <IN>;
    while(<IN>){
      chomp;
      my @sp = split /\t/;
      $cut{$sp[0]} = $sp[-1];
    }
    close IN;
    # read in sumstat files
    open IN2,"gzip -dc tQTL.sumstat.txt.gz|" or die $!;
    open OUT,">tQTL.$col.$cutoff.txt" or die $!;
    print OUT "phenotype_id\tvariant_id\tdistance\tma_samples\tma_count\tmaf\tpval_nominal\tslope\tslope_se\n";
    while(<IN2>){
      chomp;
      my @sp = split /\t/;
      next if !exists $cut{$sp[0]};
      next if $sp[6] > $cut{$sp[0]};
      print OUT "$_\n";
    }
    close IN2;
    close OUT;
  }
}



