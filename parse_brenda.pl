#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### Download brenda database from: https://www.brenda-enzymes.org/download_brenda_without_registration.php
##########################################################################
use strict;
#my $usage="$0  [p][col][file] first line is title\n"; 
#die $usage unless @ARGV ==3;
#my ($p, $col,$file) = @ARGV;
my $ID="";
my $name="";
my $def="";
open (A, "<brenda_download.txt")||die "could not open ko\n";
while (my $a=<A>){
    chomp $a;
   
    if ($a =~ /ID\s+(\S+)/){
	if($ID ne ""){
	    print $ID, "\t", $name, "\t", $def, "\n";
	    $ID=""; $name=""; $def="";
        }
	$ID=$1;
    }elsif($a =~ /RN\s+(\S+.*)$/){
	$name=$1;
    }elsif($a =~ /SN\s+(\S+.*)$/){
	$def=$1;
    }

}
close A;
