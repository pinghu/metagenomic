#! /usr/bin/perl -w
###illumina ID is located at the column 14. need to unique the data 
use strict;
#my $usage = "usage: $0 <csv> \n";
#die $usage unless @ARGV == 1;
my ($filename, $col) = @ARGV;


my %ID2Name;

open (B, "</mnt/G6_2D/db/Brenda/brenda.name")||die "could not open brenda.name\n";
while (my $x=<B>){
    chomp $x;
    my @tmp=split /\t/, $x;
    if(scalar @tmp>=2){
     $ID2Name{$tmp[0]}=$tmp[1]."\t".$tmp[2];
    }else{
	print  $x, "\n";
    }
}
close B;
 
open (A, "<$filename")||die "could not open $filename\n";
open (B, ">$filename.ECName")||die "could not open $filename.ECName\n";
while (my $xx= <A>){
    chomp $xx;
    for($xx){s/\r//gi;}
    my @tmp=split /\t/, $xx;
    if (defined $ID2Name{$tmp[$col-1]}){
	print B  $xx, "\t", $ID2Name{$tmp[$col-1]}, "\n";
    }else{
	print B $xx, "\tNA\tNA\n";
    }
	    
}
close A;
close B;
