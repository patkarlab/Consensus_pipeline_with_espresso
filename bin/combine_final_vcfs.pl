#!/usr/bin/perl
## Adam Waalkes
## takes a vcf file that is sorted by chr, position, ref allele and alt allele
## then it removes duplicates, counting them and representing the count of each variant
## Usage:combine_final_vcfs.pl <sorted vcf> <file to output new simplified vcf>

use warnings;
use strict;

my $vcf_file = $ARGV[0];
my $output_file = $ARGV[1];
my @line_split=();
my $line2=();
my $current_position=0;
my $current_ref="A";
my $current_alt="C";
my $current_count=0;
my $total_count=0;
my $current_chr=0;
my $alt_frac=0.0;
my $info_field=();
my @info_field_split=();

#print "$vcf_file\n";
#open(IP, "$vcf_file"); #pre read the first line to load variable for loop comparison
#$line2 = <IP>;
#@line_split = split("\t",$line2);
#$current_chr = $line_split[0];
#$current_position = $line_split[1];
#$current_ref = $line_split[3];
#$current_alt = $line_split[4];
#print "initialize chr $current_chr pos $current_position ref $current_ref alt $current_alt\n";
#close IP;
open(IP, "$vcf_file");

open(OUT, ">", "$output_file"); #write a VCF header appropriate for our output format
print OUT "##fileformat=VCFv4.1\n";
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##FORMAT=<ID=ALT,Number=1,Type=Integer,Description=\"Number of alts found\">\n";
print OUT "##FORMAT=<ID=TOT,Number=1,Type=Integer,Description=\"total number of reads at this location\">\n";
print OUT "##FORMAT=<ID=FRAC,Number=1,Type=Float,Description=\"fraction of alt alleles compared to total read count at this location\">\n";
print OUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1\n";

while(my $line = <IP>) {
    @line_split = split("\t",$line);
    
    if ($total_count == 0){
      @line_split = split("\t",$line);
      $current_chr = $line_split[0];
      $current_position = $line_split[1];
      $current_ref = $line_split[3];
      $current_alt = $line_split[4];
      #print OUT "*****$current_chr\t$current_position"."\t\.\t$current_ref\t$current_alt"."\t\.\tPASS"."\t\.\tGT:ALT:TOT:FRAC\t1:$current_count:$total_count:$alt_frac\n";
      $total_count=0;
      $current_count=0;
      $alt_frac=0;
    }     
    if (($current_chr eq $line_split[0]) && ($current_position == $line_split[1]) && ($current_ref eq $line_split[3]) && ($current_alt eq $line_split[4])){
       @info_field_split = split(":",$line_split[9]);
#       print "current starts $current_count, total starts $total_count\n";
       $total_count += $info_field_split[2];
       $current_count += $info_field_split[1];
#       print "combining $current_chr $current_position $current_ref $current_alt current $current_count total $total_count\n";
#       print "read current $info_field_split[1]  read total $info_field_split[2]\n";
#       print "current $current_count total $total_count\n";
    }
    else {
    
       $alt_frac = $current_count/$total_count;
       print OUT "$current_chr\t$current_position"."\t\.\t$current_ref\t$current_alt"."\t\.\tPASS"."\t\.\tGT:ALT:TOT:FRAC\t1:$current_count:$total_count:$alt_frac\n";
       $current_chr = $line_split[0];
       $current_position = $line_split[1];
       $current_ref = $line_split[3];
       $current_alt = $line_split[4];
       @info_field_split = split(":",$line_split[9]);
       $current_count = $info_field_split[1];
       $total_count = $info_field_split[2];
       
    }
   $alt_frac=0; 
}

if($total_count != 0){
  $alt_frac = $current_count/$total_count;
  print OUT "$current_chr\t$current_position"."\t\.\t$current_ref\t$current_alt"."\t\.\tPASS"."\t\.\tGT:ALT:TOT:FRAC\t1:$current_count:$total_count:$alt_frac\n";
}
#print OUT "$current_chr\t$current_position"."\t\.\t$current_ref\t$current_alt"."\t\.\tPASS"."\t\.\tGT:ALT:TOT:FRAC\t1:$current_count:$total_count:$alt_frac\n";
 
close IP;
close OUT;

exit;
