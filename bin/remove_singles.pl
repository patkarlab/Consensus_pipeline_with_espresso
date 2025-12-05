#!/usr/bin/perl

## perl remove_singles.pl <sam file>
## splits a sam file with molecular tags into a file with singletons <sam file>.singles and non singletons <sam file>.nosingles 

use warnings;
use strict;
use Switch;

my $file = $ARGV[0];
my $outfile = $ARGV[1];
my @first_linearray=();
my $current_tag="";
#my $single_count=0;
#my $non_single_count=0;
my $previous_line="";
my $previous_read="";
my $current_tag_count=1;
my $single_percent=0;
my %counts=(
     "1",0,
     "2",0,
     "3",0,
     "4",0,
     "5",0,
     "6",0,
     "7",0,
     "8",0,
     "9",0,
     "10",0,
     "11",0,
     "12",0,
     "13",0,
     "14",0,
     "15",0,
     "16",0,
     "17",0,
     "18",0,
     "19",0,
     "20",0,
     "21",0,
     "22",0,
     "23",0,
     "24",0,
     "25",0,
     "26",0,
     "27",0,
     "28",0,
     "29",0,
     "30",0,
     "31",0,
     "32",0,
     "33",0,
     "34",0,
     "35",0,
     "36",0,
     "37",0,
     "38",0,
     "39",0,
     "40",0,
     "41",0,
     "42",0,
     "43",0,
     "44",0,
     "45",0,
     "46",0,
     "47",0,
     "48",0,
     "49",0,
     "50",0,
);

open(IP, "$file");
open(OUT, ">", "$file.nosingles") or die "Can't write to file '$file.nosingles' [$!]\n";;
open(S_OUT, ">", "$file.singles") or die "Can't write to file '$file.singles' [$!]\n";;
open(OUT_COUNT, ">>", "$outfile") or die "Can't write to file '$outfile' [$!]\n";

$previous_line = <IP>;
@first_linearray = split("\t", $previous_line);
$current_tag = (split(":",$first_linearray[0]))[-1];

while(my $line = <IP>){
   my @linearray=();

   @linearray = split("\t", $line);
   if((split(":",$linearray[0]))[-1] eq $current_tag) {
      print OUT $previous_line;
#      $counts{"1"}++;
      $previous_line=$line;
      $current_tag_count++;
   }
   elsif($current_tag_count>1) {
      print OUT $previous_line;
      if($current_tag_count > 50) {
         print "Current tag count of $current_tag_count is larger than 50.  Trimming to 50\n";
         $current_tag_count=50;
      }
      $counts{$current_tag_count}++;
      $previous_line=$line;
      $current_tag_count=1;
      $current_tag=(split(":",$linearray[0]))[-1];
   }
   else {
      $counts{"1"}++;
      print S_OUT $previous_line;
      $previous_line=$line;
      $current_tag_count=1;
      $current_tag=(split(":",$linearray[0]))[-1];
   }

}
if($current_tag_count>1) {
   print OUT $previous_line;
   $counts{$current_tag_count}++;
}
else {
   print S_OUT $previous_line;
   $counts{"1"}++;
}
#$single_percent=$single_count/($single_count+$non_single_count);
#print "for $file the singleton count was $single_count, non singleton count was $non_single_count, percent single was $single_percent\n";
my $mip_name = (split("/",$file))[1];
my $tot_count = 0;
my $counter = 0;
my $non_singleton_count = 0;
print OUT_COUNT "COUNT\t$mip_name";
for my $hdr (sort {$a <=> $b} keys %counts) {
    print OUT_COUNT "\t$counts{$hdr}";
	$tot_count = $tot_count + $counts{$hdr};
	$counter ++;
	if ($counter > 1){
		$non_singleton_count = $non_singleton_count + $counts{$hdr};
	}
}
my $column_count = scalar(keys %counts);
if ( $column_count == 50 ) {
	print OUT_COUNT "\t0";
}

print OUT_COUNT "\t";
print OUT_COUNT "\t$tot_count";
print OUT_COUNT "\t$non_singleton_count";
print OUT_COUNT "\n";

close IP;
close OUT;
close S_OUT;
close OUT_COUNT;

exit;

