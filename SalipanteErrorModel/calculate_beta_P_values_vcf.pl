#!/usr/bin/perl

## Adam Waalkes
## Usage:
## perl calculate_beta_P_values_vcf.pl <beta matrix> <vcf2> <outfile>
## calculates and appends P value for each variant in a vcf 

use warnings;
use strict;
use Math::CDF 'pbeta';

my $beta_matrix_file = $ARGV[0];
my $vcf_file = $ARGV[1];
my $outfile = $ARGV[2];
my $one=0;
my $two=0;
my $onemax=0;
my $twomax=0;
my $oneheaderlines=0;
my $twoheaderlines=0;
my $matchcount=0;
my $mean=0.0;
my $alpha=0.0;
my $beta=0.0;
my $p_value=0.0;
my $one_minus_p_value=0.0;
my @final_split=();

open(IP, "$beta_matrix_file");
my(@line) = <IP>;
open(IP2, "$vcf_file");
my(@line2) = <IP2>;
open(OUT, '>' , $outfile);

$onemax=@line-1;
$twomax=@line2-1;

#skip vcf header
while(substr($line[$one],0,1) eq "#") {
   $one++;
   $oneheaderlines++;
}
while(substr($line2[$two],0,1) eq "#") {
   print OUT $line2[$two];
   $two++;
   $twoheaderlines++;
}

while(($one <= $onemax) && ($two <= $twomax)) {
   my @line_split=();
   my @line_split2=();
   chomp($line[$one]);
   chomp($line2[$two]);
   @line_split = split(/\t/,$line[$one]);
   @line_split2 = split(/\t/,$line2[$two]);
   if($line_split[0] ne $line_split2[0]) {
      if($line_split[0] gt $line_split2[0]) {
         $two++;
         next;
      }
      else {
         $one++;
         next;
      }
   }
   if($line_split[1] != $line_split2[1]) {
      if($line_split[1] > $line_split2[1]) {
         $two++;
         next;
      }
      else {
         $one++;
         next;
      }
   }
   else {
      @final_split = split(/\:/,$line_split2[9]);
      $mean=$final_split[3];
      if($line_split2[4] eq 'A') {
         $alpha = $line_split[2];
         $beta = $line_split[3];
      }
      elsif($line_split2[4] eq 'T') {
         $alpha = $line_split[4];
         $beta = $line_split[5];
      }
      elsif($line_split2[4] eq 'G') {
         $alpha = $line_split[6];
         $beta = $line_split[7];
      }
      else {
         $alpha = $line_split[8];
         $beta = $line_split[9];
      }
      $p_value = pbeta($mean,$alpha,$beta);
      $one_minus_p_value = 1-$p_value;
#      print OUT "mean = $mean, alpha = $alpha, beta = $beta, p_value = $p_value, 1-p_value = $one_minus_p_value\n";
      if($one_minus_p_value <= 0.005) { #write if p-value < 0.5%
         $line_split2[8]=$line_split2[8] . "\:PVAL";
         $line_split2[9]=$line_split2[9] . "\:$one_minus_p_value";
         print OUT "$line_split2[0]\t$line_split2[1]\t$line_split2[2]\t$line_split2[3]\t$line_split2[4]\t$line_split2[5]\t$line_split2[6]\t$line_split2[7]\t$line_split2[8]\t$line_split2[9]\n";
      }
   }
  $two++; #don't increment $one in case of multiple different variants at same loci
}

close IP;
close IP2;
close OUT;

exit;
