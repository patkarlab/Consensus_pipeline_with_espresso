#!/usr/bin/perl

## perl combine_columns.pl <infile1> <infile2> <infile3>
## reads in two files of identical length and writes to STDOUT line1\tline2\tline3 for all lines
## used to combine coordinates and counts or percentages in the error analysis process for smmips data

use warnings;
use strict;

my $infile = $ARGV[0];
my $infile2 = $ARGV[1];
my $infile3 = $ARGV[2];
my $line="";
my $line2="";
my $line3="";


open(IP, "$infile");
open(IP2, "$infile2");
open(IP3, "$infile3");
while($line = <IP>){
   $line2 = <IP2>;
   $line3 = <IP3>;
   chomp $line;
   chomp $line2;
   print "$line\t$line2\t$line3";
}

close IP;
close IP2;
close IP3;
exit;

