#! /usr/bin/perl -w
###illumina ID is located at the column 14. need to unique the data 
use strict;
#my $usage = "usage: $0 <csv> \n";
#die $usage unless @ARGV == 1;
my ($filename, $col) = @ARGV;

my %ID;
my %ID2Name;
open (A, "<$filename")||die "could not open $filename\n";

while (my $xx= <A>){
    chomp $xx;
    for($xx){s/\r//gi;}
    my @tmp=split /\t/, $xx;
    $ID{$tmp[$col-1]}=$xx;
    $ID2Name{$tmp[$col-1]}="NA";
    
}
close A;

open (B, "</mnt/G6_2D/db/Brenda/brenda.name")||die "could not open brenda.name\n";
while (my $x=<B>){
    chomp $x;
    my @tmp=split /\t/, $x;
    if(defined $ID{$tmp[0]}){
	#print $tmp[0], "\t",$tmp[1], "\t",  $tmp[2], "\n";
	$ID2Name{$tmp[0]}=$tmp[1]."\t".$tmp[2];
    }
}
close B;
 


open (A, ">$filename.ECName")||die "could not open $filename.ECName\n";
foreach my $i (keys %ID2Name){
    print A  $ID{$i},"\t", $ID2Name{$i}, "\n";
}
close A;
