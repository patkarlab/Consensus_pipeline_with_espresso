#!/usr/bin/perl

## perl add_RX_field.pl <sam file>
## copies the UMID into a new RX field 

use warnings;
use strict;

my $file =$ARGV[0];
my @first_linearray=();
my $current_tag="";
my $single_count=0;
my $non_single_count=0;
my $previous_line="";
my $previous_read="";
my $current_tag_count=1;
my $single_percent=0;
my $last_tag="";
my $first_line="";
my $UMID_count=0;
my @first_field2="";

open(IP, "$file");
open(OUT, ">", "$file.RX");
$first_line= <IP>;
@first_linearray = split("\t", $first_line);
@first_field2 = split("\:", $first_linearray[0]);
$last_tag = $first_field2[1];
close IP;
open(IP, "$file");

while(my $line = <IP>){
   my @linearray=();
   my @first_field=();
   my $first_half="";
   my $second_half="";

   chomp $line;
   @linearray = split("\t", $line);
   @first_field = split("\:", $linearray[0]);
   $first_half = substr $first_field[7],0,4;
   $second_half= substr $first_field[7],4,4;
   if($first_field[7] eq $last_tag) {
      print OUT "$line\tRG:Z:runzzz.1\tMI:Z:$UMID_count\tQX:Z:AAAA-AAAA\tRX:Z:$first_half-$second_half\n";
   }
   else {
      $last_tag=$first_field[7];
      $UMID_count++;
      print OUT "$line\tRG:Z:runzzz.1\tMI:Z:$UMID_count\tQX:Z:AAAA-AAAA\tRX:Z:$first_half-$second_half\n";

   }
}

close IP;
close OUT;

exit;

