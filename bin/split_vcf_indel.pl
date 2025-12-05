#!/usr/bin/perl

## perl split_vcf_indel.pl <infile>
## input is a vcf file this routine splits out indels from point mutations
## creating $infile.indel and $infile.point

use warnings;
use strict;

my $infile =$ARGV[0];
my @line_split=();
my $line="";

open(IP, "$infile");
open(OUT, ">", "$infile.point");
open(OUT2, ">", "$infile.indel");

while($line = <IP>){
   chomp $line;
   @line_split = split(/\t/,$line);
   if (substr($line,0,1) eq "#") {
       print OUT "$line\n";
       print OUT2 "$line\n";
   }
   elsif (((length $line_split[3])==1) && ((length $line_split[4])==1)){
      print OUT "$line\n";      
   }
   else {
      print OUT2 "$line\n";
   }
}
close IP;
close OUT;
close OUT2;

exit;

