#!/usr/bin/perl
## Adam Waalkes
## takes a vcf file and counts the different kinds of mutations outputing to standard out total
## counts and percentage of total also creates a vcf file for each type of mutation and indels 
## Usage:$perl mutation_types_and_split_vcf_error_rates.pl <vcf> 
## Output:
## <vcf>.AC <vcf>.GT split vcfs by mutation type
## <vcf>.count.AC <vcf>.count.GT alt allele count by mutation type
## writes to stdout count and % of each type of mutaion, total mutations and mean and stdev of alt allele frequencies by type

use warnings;
use strict;

sub average{
   my $data = @_;
   if ( $data == 0 ) {
      die("Empty arrayn");
   }
   my $total = 0;
   foreach (@_) {
      $total += $_;
   }
   my $average = $total / $data;
   return $average;
}

sub stdev{
   my $data = @_;
   if($data == 1){
           return 0;
   }
   my $average = &average( @_ );
   my $sqtotal = 0;
   foreach(@_) {
      $sqtotal += ($average-$_) ** 2;
   }
   my $std = ($sqtotal / ($data-1)) ** 0.5;
   return $std;
}

my $vcf_file = $ARGV[0];
my @line_split=();
my @sample_split=();
my $total_count=0;
my $AC=0;
my $AG=0;
my $AT=0;
my $A_total=0;
my $CA=0;
my $CG=0;
my $CT=0;
my $C_total=0;
my $GA=0;
my $GC=0;
my $GT=0;
my $G_total=0;
my $TA=0;
my $TC=0;
my $TG=0;
my $T_total=0;
my $INDEL=0;
my $AC_per=0;
my $AG_per=0;
my $AT_per=0;
my $A_per=0;
my $CA_per=0;
my $CG_per=0;
my $CT_per=0;
my $C_per=0;
my $GA_per=0;
my $GC_per=0;
my $GT_per=0;
my $G_per=0;
my $TA_per=0;
my $TC_per=0;
my $TG_per=0;
my $T_per=0;
my $INDEL_per=0;
my @A_err=();
my @AC_err=();
my @AG_err=();
my @AT_err=();
my @C_err=();
my @CA_err=();
my @CG_err=();
my @CT_err=();
my @G_err=();
my @GA_err=();
my @GC_err=();
my @GT_err=();
my @T_err=();
my @TA_err=();
my @TC_err=();
my @TG_err=();
my @INDEL_err=();
my $het_floor=0.18;
my $AC_mean=0;
my $AG_mean=0;
my $AT_mean=0;
my $A_mean=0;
my $CA_mean=0;
my $CG_mean=0;
my $CT_mean=0;
my $C_mean=0;
my $GA_mean=0;
my $GC_mean=0;
my $GT_mean=0;
my $G_mean=0;
my $TA_mean=0;
my $TC_mean=0;
my $TG_mean=0;
my $T_mean=0;
my $INDEL_mean=0;
my $AC_stdev=0;
my $AG_stdev=0;
my $AT_stdev=0;
my $A_stdev=0;
my $CA_stdev=0;
my $CG_stdev=0;
my $CT_stdev=0;
my $C_stdev=0;
my $GA_stdev=0;
my $GC_stdev=0;
my $GT_stdev=0;
my $G_stdev=0;
my $TA_stdev=0;
my $TC_stdev=0;
my $TG_stdev=0;
my $T_stdev=0;
my $INDEL_stdev=0;
my $A_count_in_MIPS=4467; #from previous count of number of As in the area covered by panel
my $C_count_in_MIPS=3510;
my $G_count_in_MIPS=3710;
my $T_count_in_MIPS=4626;

open(IP, "$vcf_file");
open(ACIP,">","$vcf_file.AC");
open(AGIP,">","$vcf_file.AG");
open(ATIP,">","$vcf_file.AT");
open(CAIP,">","$vcf_file.CA");
open(CGIP,">","$vcf_file.CG");
open(CTIP,">","$vcf_file.CT");
open(GCIP,">","$vcf_file.GC");
open(GAIP,">","$vcf_file.GA");
open(GTIP,">","$vcf_file.GT");
open(TCIP,">","$vcf_file.TC");
open(TGIP,">","$vcf_file.TG");
open(TAIP,">","$vcf_file.TA");
open(INDELIP,">","$vcf_file.INDEL");
open(ACIPC,">","$vcf_file.count.AC");
open(AGIPC,">","$vcf_file.count.AG");
open(ATIPC,">","$vcf_file.count.AT");
open(CAIPC,">","$vcf_file.count.CA");
open(CGIPC,">","$vcf_file.count.CG");
open(CTIPC,">","$vcf_file.count.CT");
open(GCIPC,">","$vcf_file.count.GC");
open(GAIPC,">","$vcf_file.count.GA");
open(GTIPC,">","$vcf_file.count.GT");
open(TCIPC,">","$vcf_file.count.TC");
open(TGIPC,">","$vcf_file.count.TG");
open(TAIPC,">","$vcf_file.count.TA");
open(INDELIPC,">","$vcf_file.count.INDEL");

while(my $line = <IP>) {
   if (substr($line,0,1) ne "#") {
      @line_split = split("\t",$line);
      if (($line_split[3] eq "A") && ($line_split[4] eq "C")){
         $AC++;
         print ACIP $line;
         @sample_split = split(":",$line_split[9]);
         print ACIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @AC_err,$sample_split[2];
         } 
      }
      elsif (($line_split[3] eq "A") && ($line_split[4] eq "G")){
         $AG++;
         print AGIP $line;
         @sample_split = split(":",$line_split[9]);
         print AGIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @AG_err,$sample_split[2];
         }
      }
      elsif (($line_split[3] eq "A") && ($line_split[4] eq "T")){
         $AT++;
         print ATIP $line;
         @sample_split = split(":",$line_split[9]);
         print ATIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @AT_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "C") && ($line_split[4] eq "A")){
         $CA++;
         print CAIP $line;
         @sample_split = split(":",$line_split[9]);
         print CAIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @CA_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "C") && ($line_split[4] eq "G")){
         $CG++;
         print CGIP $line;
         @sample_split = split(":",$line_split[9]);
         print CGIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @CG_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "C") && ($line_split[4] eq "T")){
         print CTIP $line;
         $CT++;
         @sample_split = split(":",$line_split[9]);
         print CTIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @CT_err,$sample_split[2];
         }
      }
      elsif (($line_split[3] eq "G") && ($line_split[4] eq "A")){
         $GA++;
         print GAIP $line;
         @sample_split = split(":",$line_split[9]);
         print GAIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @GA_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "G") && ($line_split[4] eq "C")){
         $GC++;
         print GCIP $line;
         @sample_split = split(":",$line_split[9]);
         print GCIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @GC_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "G") && ($line_split[4] eq "T")){
         $GT++;
         print GTIP $line;
         @sample_split = split(":",$line_split[9]);
         print GTIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @GT_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "T") && ($line_split[4] eq "A")){
         $TA++;
         print TAIP $line;
         @sample_split = split(":",$line_split[9]);
         print TAIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @TA_err,$sample_split[2];
         }
      }
      elsif (($line_split[3] eq "T") && ($line_split[4] eq "C")){
         $TC++;
         print TCIP $line;
         @sample_split = split(":",$line_split[9]);
         print TCIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @TC_err,$sample_split[2];
         }
       }
      elsif (($line_split[3] eq "T") && ($line_split[4] eq "G")){
         $TG++;
         print TGIP $line;
         @sample_split = split(":",$line_split[9]);
         print TGIPC "$sample_split[1]\n";
         if ($sample_split[2] <= $het_floor) {
            push @TG_err,$sample_split[2];
         }
       }
      else {
         $INDEL++;
         print INDELIP $line;
       }
   }
}

$A_total=$AC+$AG+$AT;
$C_total=$CA+$CG+$CT;
$G_total=$GC+$GA+$GT;
$T_total=$TC+$TG+$TA;
$total_count=$AC+$AG+$AT+$CA+$CG+$CT+$GA+$GC+$GT+$TA+$TC+$TG+$INDEL;
$AC_per=100*$AC/$total_count;
$AG_per=100*$AG/$total_count;
$AT_per=100*$AT/$total_count;
$A_per=100*$A_total/$total_count;
$CA_per=100*$CA/$total_count;
$CG_per=100*$CG/$total_count;
$CT_per=100*$CT/$total_count;
$C_per=100*$C_total/$total_count;
$GA_per=100*$GA/$total_count;
$GC_per=100*$GC/$total_count;
$GT_per=100*$GT/$total_count;
$G_per=100*$G_total/$total_count;
$TA_per=100*$TA/$total_count;
$TC_per=100*$TC/$total_count;
$TG_per=100*$TG/$total_count;
$T_per=100*$T_total/$total_count;
$INDEL_per=100*$INDEL/$total_count;

@A_err = @AC_err;
push(@A_err, @AG_err);
push(@A_err, @AT_err);
$A_total=$AC+$AG+$AT;
for (my $i=$A_total;$i<$A_count_in_MIPS;$i++) {
   push @A_err,"0";
}
$A_mean = &average(@A_err);
$A_stdev = &stdev(@A_err);

@C_err = @CA_err;
push(@C_err, @CG_err);
push(@C_err, @CT_err);
$C_total=$CA+$CG+$CT;
for (my $i=$C_total;$i<$C_count_in_MIPS;$i++) {
   push @C_err,"0";
}
$C_mean = &average(@C_err);
$C_stdev = &stdev(@C_err);

@G_err = @GA_err;
push(@G_err, @GC_err);
push(@G_err, @GT_err);
$G_total=$GC+$GT+$GA;
for (my $i=$G_total;$i<$G_count_in_MIPS;$i++) {
   push @G_err,"0";
}
$G_mean = &average(@G_err);
$G_stdev = &stdev(@G_err);

@T_err = @TA_err;
push(@T_err, @TC_err);
push(@T_err, @TG_err);
$T_total=$TC+$TG+$TA;
for (my $i=$T_total;$i<$T_count_in_MIPS;$i++) {
   push @T_err,"0";
}
$T_mean = &average(@T_err);
$T_stdev = &stdev(@T_err);

for (my $i=$AC;$i<$A_count_in_MIPS;$i++) {
   push @AC_err,"0";
}
$AC_mean = &average(@AC_err);
$AC_stdev = &stdev(@AC_err);
for (my $i=$AG;$i<$A_count_in_MIPS;$i++) {
   push @AG_err,"0";
}
$AG_mean = &average(@AG_err);
$AG_stdev = &stdev(@AG_err);
for (my $i=$AT;$i<$A_count_in_MIPS;$i++) {
   push @AT_err,"0";
}
$AT_mean = &average(@AT_err);
$AT_stdev = &stdev(@AT_err);
for (my $i=$CA;$i<$C_count_in_MIPS;$i++) {
   push @CA_err,"0";
}
$CA_mean = &average(@CA_err);
$CA_stdev = &stdev(@CA_err);
for (my $i=$CG;$i<$C_count_in_MIPS;$i++) {
   push @CG_err,"0";
}
$CG_mean = &average(@CG_err);
$CG_stdev = &stdev(@CG_err);
for (my $i=$CT;$i<$C_count_in_MIPS;$i++) {
   push @CT_err,"0";
}
$CT_mean = &average(@CT_err);
$CT_stdev = &stdev(@CT_err);
for (my $i=$GA;$i<$G_count_in_MIPS;$i++) {
   push @GA_err,"0";
}
$GA_mean = &average(@GA_err);
$GA_stdev = &stdev(@GA_err);
for (my $i=$GC;$i<$G_count_in_MIPS;$i++) {
   push @GC_err,"0";
}
$GC_mean = &average(@GC_err);
$GC_stdev = &stdev(@GC_err);
for (my $i=$GT;$i<$G_count_in_MIPS;$i++) {
   push @GT_err,"0";
}
$GT_mean = &average(@GT_err);
$GT_stdev = &stdev(@GT_err);
for (my $i=$TA;$i<$T_count_in_MIPS;$i++) {
   push @TA_err,"0";
}
$TA_mean = &average(@TA_err);
$TA_stdev = &stdev(@TA_err);
for (my $i=$TC;$i<$T_count_in_MIPS;$i++) {
   push @TC_err,"0";
}
$TC_mean = &average(@TC_err);
$TC_stdev = &stdev(@TC_err);
for (my $i=$TG;$i<$T_count_in_MIPS;$i++) {
   push @TG_err,"0";
}
$TG_mean = &average(@TG_err);
$TG_stdev = &stdev(@TG_err);

printf("Type\tcount\t%%\tmean per base error\tstdev\n");
printf("A->*\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$A_total,$A_per,$A_mean,$A_stdev);
printf("C->*\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$C_total,$C_per,$C_mean,$C_stdev);
printf("G->*\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$G_total,$G_per,$G_mean,$G_stdev);
printf("T->*\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$T_total,$T_per,$T_mean,$T_stdev);
printf("Indel\t%d\t%.1f\tNA\t\t\tNA\n",$INDEL,$INDEL_per);
print "Total variants:\t$total_count\n\n";
printf("A->C\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$AC,$AC_per,$AC_mean,$AC_stdev);
printf("A->G\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$AG,$AG_per,$AG_mean,$AG_stdev);
printf("A->T\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$AT,$AT_per,$AT_mean,$AT_stdev);
printf("C->A\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$CA,$CA_per,$CA_mean,$CA_stdev);
printf("C->G\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$CG,$CG_per,$CG_mean,$CG_stdev);
printf("C->T\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$CT,$CT_per,$CT_mean,$CT_stdev);
printf("G->A\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$GA,$GA_per,$GA_mean,$GA_stdev);
printf("G->C\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$GC,$GC_per,$GC_mean,$GC_stdev);
printf("G->T\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$GT,$GT_per,$GT_mean,$GT_stdev);
printf("T->A\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$TA,$TA_per,$TA_mean,$TA_stdev);
printf("T->C\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$TC,$TC_per,$TC_mean,$TC_stdev);
printf("T->G\t%d\t%.1f\t%1.8f\t\t%1.8f\n",$TG,$TG_per,$TG_mean,$TG_stdev);
printf("Indel\t%d\t%.1f\tNA\t\t\tNA\n",$INDEL,$INDEL_per);
print "Total variants:\t$total_count\n";

close IP;
close ACIP;
close AGIP;
close ATIP;
close CAIP;
close CGIP;
close CTIP;
close GAIP;
close GCIP;
close GTIP;
close TAIP;
close TGIP;
close TCIP;
close INDELIP;
close ACIPC;
close AGIPC;
close ATIPC;
close CAIPC;
close CGIPC;
close CTIPC;
close GAIPC;
close GCIPC;
close GTIPC;
close TAIPC;
close TGIPC;
close TCIPC;
close INDELIPC;

exit;

