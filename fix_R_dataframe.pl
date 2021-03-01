#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV ==1;
my ($file) = @ARGV;

open (A, "<$file")||die "could not open $file\n";
my $x=<A>;
chomp $x;
open (B, ">$file.xls")||die "could not open $file.xls\n";
for($x){s/\"//gi;s/\r//gi;}
print B "ID\t$x\n";
while (my $a=<A>){
    chomp $a;
    for($a){s/\"//gi;s/\r//gi;}
    print B $a, "\n";
}
close A;
close B;
