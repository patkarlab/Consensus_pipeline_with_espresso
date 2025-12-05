#!/usr/bin/env perl
## Adam Waalkes
## takes a mpileup file and calls variants, creates a vcf of the variants
## Usage: call_variants_from_mpileup.pl <mileup> <minimum base quality to call> <vcf file to output new simplified vcf>

use warnings;
use strict;
use Switch;

my %reads=();
my $mpileup_file = $ARGV[0];
my $min_qual = $ARGV[1];
my $vcf_file = $ARGV[2];
my $indel_count = 0;
my $indel_string = "";


open(IP, "$mpileup_file");
open(OUT, ">", "$vcf_file");
print OUT "##fileformat=VCFv4.1\n";
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##FORMAT=<ID=ALT,Number=1,Type=Integer,Description=\"Number of alts found\">\n";
print OUT "##FORMAT=<ID=TOT,Number=1,Type=Integer,Description=\"total number of molecular tags at this location\">\n";
print OUT "##FORMAT=<ID=FRAC,Number=1,Type=Float,Description=\"fraction of alt alleles compared to total molecular tags at this location\">\n";
print OUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1\n";

while(my $line = <IP>) {
   my @line_split = split("\t",$line);
   my $index_difference=0; #difference between the variant string and the qual score relevant for insertions and read start
   my %current_variants=();
   my %indel_variants=();
   my $alt_frac=0;
   #chr is in $line_split[0]
   #coord is in $line_split[1]
   #ref is in $line_split[2]
   #read count is in $line_split[3]
   #variant string is in $line_split[4]
   #qual for base is in $line_split[5]
   for(my $counter=0; $counter < length($line_split[4]); $counter++) {
      switch (substr($line_split[4],$counter,1)) {
         case /[.,\*]/ {next;}
         case /[\$]/ {   #$counter=$counter++;
                         $index_difference++;
                         next;}
         case /[\^]/ {   $counter=$counter+1;$index_difference=$index_difference+2;} # read start, ignore nut skip the inline qual also increment the $index_difference
         case /[gG]/ {if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                         $current_variants{"G"}++;
                      }
                      else {
                        $line_split[3]--;
                      }
         }
         case /[aA]/ {if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                         $current_variants{"A"}++;
                      }
                      else {
                         $line_split[3]--;
                      }
          }
         case /[cC]/ {if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                         $current_variants{"C"}++;
                     }
                     else {
                        $line_split[3]--;
                     }
          }
         case /[tT]/ {if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                         $current_variants{"T"}++;
                      }
                      else {
                         $line_split[3]--;
                      }
         }
         case /[-]/ {
                    if(((substr($line_split[4],$counter+2,1)) ge "0") && ((substr($line_split[4],$counter+2,1)) le "9")) {
                       $indel_count = substr($line_split[4],$counter+1,2);
                       $indel_string = "-" . $line_split[2] . substr($line_split[4],$counter+3,$indel_count);
                       if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                          $indel_variants{uc $indel_string}++;
                       }
                       else {
                          $line_split[3]--;
                       }
                       $counter=$counter+$indel_count+3;
                       $index_difference = $index_difference+$indel_count+3;
                    }
                    else {
                       $indel_count = substr($line_split[4],$counter+1,1);
                       $indel_string = "-" . $line_split[2] . substr($line_split[4],$counter+2,$indel_count);
                       if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                          $indel_variants{uc $indel_string}++;
                       }
                       else {
                          $line_split[3]--;
                       }
                       $counter=$counter+$indel_count+2;
                       $index_difference = $index_difference+$indel_count+2;

                    }

         }
         case /[+]/ { 

                    if(((substr($line_split[4],$counter+2,1)) ge "0") && ((substr($line_split[4],$counter+2,1)) le "9")) {
                       $indel_count = substr($line_split[4],$counter+1,2);
                       $indel_string = "+" . $line_split[2] . substr($line_split[4],$counter+3,$indel_count);
                       if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                          $indel_variants{uc $indel_string}++;
                       }
                       else {
                          $line_split[3]--;
                       }
                       $counter=$counter+$indel_count+3;
                       $index_difference = $index_difference+$indel_count+3;
                    }
                    else {
                       $indel_count = substr($line_split[4],$counter+1,1);
                       $indel_string = "+" . $line_split[2] . substr($line_split[4],$counter+2,$indel_count);
                       if((substr($line_split[5],$counter-$index_difference,1) gt $min_qual)) {
                          $indel_variants{uc $indel_string}++;
                       }
                       else {
                          $line_split[3]--;
                       }
                       $counter=$counter+$indel_count+2;
                       $index_difference = $index_difference+$indel_count+2;

                    }
         }
      }
   }
   foreach my $hdr(sort {$current_variants{$b} <=> $current_variants{$a}} keys %current_variants){
      $alt_frac=$current_variants{$hdr}/$line_split[3];
      print OUT "$line_split[0]\t$line_split[1]\t.\t$line_split[2]\t$hdr\t.\tPASS\t.\tGT:ALT:TOT:FRAC\t1:$current_variants{$hdr}:$line_split[3]:$alt_frac\n";
    }

   foreach my $hdr(sort {$indel_variants{$b} <=> $indel_variants{$a}} keys %indel_variants){
      $alt_frac=$indel_variants{$hdr}/$line_split[3];
      if((substr($hdr,0,1)) eq "-") {
         my $indel_out = substr($hdr,1);
         print OUT "$line_split[0]\t$line_split[1]\t.\t$indel_out\t$line_split[2]\t.\tPASS\t.\tGT:ALT:TOT:FRAC\t1:$indel_variants{$hdr}:$line_split[3]:$alt_frac\n";
      }
      else {
         my $indel_out = substr($hdr,1);
         print OUT "$line_split[0]\t$line_split[1]\t.\t$line_split[2]\t$indel_out\t.\tPASS\t.\tGT:ALT:TOT:FRAC\t1:$indel_variants{$hdr}:$line_split[3]:$alt_frac\n";
      }
   }
}


close IP;
close OUT;
exit;
