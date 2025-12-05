#!/usr/bin/perl
## Adam Waalkes
## simplify_vcf2.pl was created to be used by the pipeline when we use CallMolecularConsensusReads in which case we can get multiple
## UMIDs combined into the vcf so we need to count differently than simplify_vcf
##
## takes a merged vcf file that is sorted by position and combines duplicates adding them properly it also 
## removes all variants that are not in capture area (mip arms) as defined by lowarm and higharm and them modifies
## the formating simplifying the data and representing the count of each identical entry 
## Usage: simplify_vcf2.pl <low arm coordinate> <high arm coordinate> <sorted vcf> <total number of unique molecular tags> <file to output new simplified vcf>


use warnings;
use strict;

my %reads=();
my $lowarm = $ARGV[0];
my $higharm = $ARGV[1];
my $vcf_file = $ARGV[2];
my $number_moletags = $ARGV[3];
my $output_file = $ARGV[4];
my @line_split=();
my $line2=();
my $current_name=();
my $current_position=0;
my $current_ref="A";
my $current_alt="C";
my $current_count=0;
my $alt_frac=0.0;
my @name_split=();

open(IP, "$vcf_file"); #pre read the first line to load variable for loop comparison
$line2 = <IP>;
@line_split = split("\t",$line2);
$current_name = $line_split[0];
$current_position = $line_split[1];
$current_ref = $line_split[3];
$current_alt = $line_split[4];
close IP;
open(IP, "$vcf_file");
while(my $line = <IP>) { #skip all snps in the front arm
    @line_split = split("\t",$line);
    my @allele_count_split = split("\:",$line_split[9]);
    if ($line_split[1] > $lowarm) {
       $current_name = $line_split[0];
       $current_position = $line_split[1];
       $current_ref = $line_split[3];
       $current_alt = $line_split[4];
       $current_count = $allele_count_split[1];
       last;
    }
}

open(OUT, ">", "$output_file"); #write a VCF header appropriate for our output format
print OUT "##fileformat=VCFv4.1\n";
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##FORMAT=<ID=ALT,Number=1,Type=Integer,Description=\"Number of alts found\">\n";
print OUT "##FORMAT=<ID=TOT,Number=1,Type=Integer,Description=\"total number of molecular tags at this location\">\n";
print OUT "##FORMAT=<ID=FRAC,Number=1,Type=Float,Description=\"fraction of alt alleles compared to total molecular tags at this location\">\n";
print OUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1\n";

if ($current_position >= $higharm) { #if the first variant that is bigger than the lowarm is also bigger than the higharm we are done no variants
   print "No snps detected in between $lowarm and $higharm\n";
   close IP;
   close OUT;
   exit;
}

while(my $line = <IP>) {
    @line_split = split("\t",$line);
    my @allele_count_split = split("\:",$line_split[9]);
    if ($line_split[1] >= $higharm) { #skip all snps in the high arm
       last;
    }
    if (($current_position == $line_split[1]) && ($current_ref eq $line_split[3]) && ($current_alt eq $line_split[4]) && ($line_split[1] > $lowarm)){
       $current_count+=$allele_count_split[1];
    }
    else {
       $alt_frac = $current_count/$number_moletags;
       if ($current_name =~ ":") {
          @name_split = split (":",$current_name);
       }
       else { #special case if there is only one vcf then grep returns just chr17 not 7578412.TCCGGGTG.sam.snp.vcf:chr17
          $name_split[1] = $current_name;
       }
       print OUT "$name_split[1]\t$current_position"."\t\.\t$current_ref\t$current_alt"."\t\.\tPASS"."\t\.\t"."GT:ALT:TOT:FRAC\t1:$current_count:$number_moletags:$alt_frac\n";
       $current_name = $line_split[0];
       $current_position = $line_split[1];
       $current_ref = $line_split[3];
       $current_alt = $line_split[4];
       $current_count = $allele_count_split[1];
    }
}

if ($current_count > 0) { #only print last one if there is one
   $alt_frac = $current_count/$number_moletags;
   if ($current_name =~ ":") {
      @name_split = split (":",$current_name);
   }
   else { #special case if there is only one vcf then grep returns just chr17 not 7578412.TCCGGGTG.sam.snp.vcf:chr17
      $name_split[1] = $current_name;
   }
   print OUT "$name_split[1]\t$current_position"."\t\.\t$current_ref\t$current_alt"."\t\.\tPASS"."\t\.\t"."GT:ALT:TOT:FRAC\t1:$current_count:$number_moletags:$alt_frac\n";
}
close IP;
close OUT;
exit;
