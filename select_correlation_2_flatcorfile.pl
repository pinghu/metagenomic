#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  <file> \n"; 
die $usage unless @ARGV == 3;
my ($spearmanfile, $pearsonfile, $outname) = @ARGV;
my %result1P;
my %result1Cor;
#my %result;
open (A, "<$spearmanfile")||die "could not open $spearmanfile\n";
my $tt=<A>;
while (my $a=<A>){
    
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    if($tmp[4]<=0.05){
	$result1P{$tmp[1]}{$tmp[2]}=$tmp[4];
	$result1Cor{$tmp[1]}{$tmp[2]}=$tmp[3]; 

    }
    
}
close A;

open (A, "<$pearsonfile")||die "could not open $pearsonfile\n";
my $tt=<A>;
open (B, ">$outname.corp05")||die "could not open $outname.corp05\n";
print B "variable1\tvariable1\tP.$spearmanfile\tP.$pearsonfile\tCor.$spearmanfile\tCor.$pearsonfile\n";
open (C, ">$outname.cor5p05")||die "could not open $outname.cor5p05\n";
print C "variable1\tvariable1\tP.$spearmanfile\tP.$pearsonfile\tCor.$spearmanfile\tCor.$pearsonfile\n";

while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
    if($tmp[4]<=0.05){
	if(defined $result1P{$tmp[1]}{$tmp[2]}){
	    print B $tmp[1], "\t", $tmp[2], "\t", $result1P{$tmp[1]}{$tmp[2]}, "\t",$tmp[4], "\t",
		$result1Cor{$tmp[1]}{$tmp[2]}, "\t", $tmp[3], "\n";
	    if((abs($tmp[3]) >=0.5) &&(abs(	$result1Cor{$tmp[1]}{$tmp[2]})>=0.5)&& ($tmp[3] *$result1Cor{$tmp[1]}{$tmp[2]} >0) ){
		print C $tmp[1], "\t", $tmp[2], "\t", $result1P{$tmp[1]}{$tmp[2]}, "\t",$tmp[4], "\t",
		$result1Cor{$tmp[1]}{$tmp[2]}, "\t", $tmp[3], "\n";
	    }
	}
    }
}
close A;
close B;
close C;
