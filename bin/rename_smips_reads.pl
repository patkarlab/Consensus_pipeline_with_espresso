#!/usr/bin/env perl
## takes a fastq file of smmips generated dna and parses the entire file into subfiles based on the molecular tag 
## Usage: rename_smips_reads.pl <*.fastq>  

use warnings;
use strict;
use Switch;

my %reads=();
my $fastq_file = $ARGV[0];	
my $fasta_file = $fastq_file;
$fasta_file =~  s/\.gz$//;	# Removing the gz file extension
my $output_file = $ARGV[1];
my @name_array=();
my $chomped_name=();
my $first_half_tag=();
my $second_half_tag=();
my $nucleotides=();
my $sep="";
my $qual="";
my @args=();

open(IP, "zcat $fastq_file |");
open(AIP, ">", "$fasta_file.A");
open(CIP, ">", "$fasta_file.C");
open(GIP, ">", "$fasta_file.G");
open(TIP, ">", "$fasta_file.T");

while(my $line = <IP>) {
    $nucleotides=<IP>;
    $sep=<IP>;
    $qual=<IP>;
    switch (substr($nucleotides,0,1)) {
       case "A" {
          print AIP $line;
          print AIP $nucleotides;
          print AIP $sep;
          print AIP $qual;
       }
       case "C" {
          print CIP $line;
          print CIP $nucleotides;
          print CIP $sep;
          print CIP $qual;
       }
       case "G" {
          print GIP $line;
          print GIP $nucleotides;
          print GIP $sep;
          print GIP $qual;
       }
       case "T" {
          print TIP $line;
          print TIP $nucleotides;
          print TIP $sep;
          print TIP $qual;
       }
       else {print "Non nucleotide value in tag\n"}
    }
}

close IP;

close AIP;
close CIP;
close GIP;
close TIP;

open(IP,">","$output_file");
open(AIP, "$fasta_file.A");

while(my $line = <AIP>) {
    $reads{$line}{nucleotides}=<AIP>;
    $reads{$line}{sep}=<AIP>;
    $reads{$line}{qual_score}=<AIP>;
    $first_half_tag=substr($reads{$line}{nucleotides},0,4);
    $second_half_tag=substr($reads{$line}{nucleotides},length($reads{$line}{nucleotides})-5,4);
    $reads{$line}{mole_tag}=$first_half_tag . $second_half_tag . "\n";

}

foreach my $name(sort {$reads{$a}{mole_tag} cmp $reads{$b}{mole_tag}} keys %reads) {
   @name_array = split(' ',$name);
   $chomped_name=$name_array[0];
   chomp $chomped_name;
   if (length($reads{$name}{nucleotides}) > 53) { # for cellfree mips panel arm is 42, min capture is 8, still have bardcode which is 8 so smallest we want is 58
      print IP "$chomped_name:";                  # but we want to allow for indels to cut the number to 53
      print IP $reads{$name}{mole_tag};
      print IP substr($reads{$name}{nucleotides},4,length($reads{$name}{nucleotides})-9);
      print IP "\n$reads{$name}{sep}";
      print IP substr($reads{$name}{qual_score},4,length($reads{$name}{qual_score})-9);
      print IP "\n";
   }
} 

close AIP;
%reads=();
open(CIP, "$fasta_file.C");
while(my $line = <CIP>) {
    $reads{$line}{nucleotides}=<CIP>;
    $reads{$line}{sep}=<CIP>;
    $reads{$line}{qual_score}=<CIP>;
    $first_half_tag=substr($reads{$line}{nucleotides},0,4);
    $second_half_tag=substr($reads{$line}{nucleotides},length($reads{$line}{nucleotides})-5,4);
    $reads{$line}{mole_tag}=$first_half_tag . $second_half_tag . "\n";

}

foreach my $name(sort {$reads{$a}{mole_tag} cmp $reads{$b}{mole_tag}} keys %reads) {
   @name_array = split(' ',$name);
   $chomped_name=$name_array[0];
   chomp $chomped_name;
   if (length($reads{$name}{nucleotides}) > 53) { # for cellfree mips panel arm is 42, min capture is 8, still have bardcode which is 8 so smallest we want is 58
      print IP "$chomped_name:";                  # but we want to allow for indels to cut the number to 53
      print IP $reads{$name}{mole_tag};
      print IP substr($reads{$name}{nucleotides},4,length($reads{$name}{nucleotides})-9);
      print IP "\n$reads{$name}{sep}";
      print IP substr($reads{$name}{qual_score},4,length($reads{$name}{qual_score})-9);
      print IP "\n";
   }
}

close CIP;
%reads=();
open(GIP, "$fasta_file.G");
while(my $line = <GIP>) {
    $reads{$line}{nucleotides}=<GIP>;
    $reads{$line}{sep}=<GIP>;
    $reads{$line}{qual_score}=<GIP>;
    $first_half_tag=substr($reads{$line}{nucleotides},0,4);
    $second_half_tag=substr($reads{$line}{nucleotides},length($reads{$line}{nucleotides})-5,4);
    $reads{$line}{mole_tag}=$first_half_tag . $second_half_tag . "\n";

}

foreach my $name(sort {$reads{$a}{mole_tag} cmp $reads{$b}{mole_tag}} keys %reads) {
   @name_array = split(' ',$name);
   $chomped_name=$name_array[0];
   chomp $chomped_name;
   if (length($reads{$name}{nucleotides}) > 53) { # for cellfree mips panel arm is 42, min capture is 8, still have bardcode which is 8 so smallest we want is 58
      print IP "$chomped_name:";                  # but we want to allow for indels to cut the number to 53
      print IP $reads{$name}{mole_tag};
      print IP substr($reads{$name}{nucleotides},4,length($reads{$name}{nucleotides})-9);
      print IP "\n$reads{$name}{sep}";
      print IP substr($reads{$name}{qual_score},4,length($reads{$name}{qual_score})-9);
      print IP "\n";
   }
}
close GIP;
%reads=();
open(TIP, "$fasta_file.T");
while(my $line = <TIP>) {
    $reads{$line}{nucleotides}=<TIP>;
    $reads{$line}{sep}=<TIP>;
    $reads{$line}{qual_score}=<TIP>;
    $first_half_tag=substr($reads{$line}{nucleotides},0,4);
    $second_half_tag=substr($reads{$line}{nucleotides},length($reads{$line}{nucleotides})-5,4);
    $reads{$line}{mole_tag}=$first_half_tag . $second_half_tag . "\n";

}

foreach my $name(sort {$reads{$a}{mole_tag} cmp $reads{$b}{mole_tag}} keys %reads) {
   @name_array = split(' ',$name);
   $chomped_name=$name_array[0];
   chomp $chomped_name;
   if (length($reads{$name}{nucleotides}) > 53) { # for cellfree mips panel arm is 42, min capture is 8, still have bardcode which is 8 so smallest we want is 58
      print IP "$chomped_name:";                  # but we want to allow for indels to cut the number to 53
      print IP $reads{$name}{mole_tag};
      print IP substr($reads{$name}{nucleotides},4,length($reads{$name}{nucleotides})-9);
      print IP "\n$reads{$name}{sep}";
      print IP substr($reads{$name}{qual_score},4,length($reads{$name}{qual_score})-9);
      print IP "\n";
   }
}
close TIP;

@args = ("rm", "-f", "$fasta_file.A");
system(@args) == 0 or die "system @args failed: $?";
@args = ("rm", "-f", "$fasta_file.C");
system(@args) == 0 or die "system @args failed: $?";
@args = ("rm", "-f", "$fasta_file.G");
system(@args) == 0 or die "system @args failed: $?";
@args = ("rm", "-f", "$fasta_file.T");
system(@args) == 0 or die "system @args failed: $?";

close IP;
exit;
