#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### 
##########################################################################
use strict;
my $usage="$0  [file] first line is title\n"; 
die $usage unless @ARGV==1;
my ($file, $col, $metafile) = @ARGV;
my $metafile="MDEWhiteLipidMeta2.txt";
#my $col=7; #SG area um2
#lipid column 22-38

open (A, "<$metafile")||die "could not open $metafile\n";
my %site;
my %ageGrp;
my %gender;
my %age;
my %race;
my %skin;
my %id3;
my %lipid;
my $tt=<A>;
chomp $tt;
my @title=split /\t/, $tt;
my %lipidT;
my %cutC;
#my %shannon;
my %appear_age;
for(my $j=22; $j<38; $j++){
    $lipidT{$j}=$title[$j];
    $cutC{$j}="cut -f1";
}

while (my $a=<A>){
    chomp $a;
    for($a){s/\r//gi;}
    my @tmp=split /\t/, $a;
   
    $site{$tmp[0]}=$tmp[4];
    $id3{$tmp[0]}=$tmp[5];
    $age{$tmp[0]}=$tmp[6];
    $appear_age{$tmp[0]}=$tmp[8];
    $ageGrp{$tmp[0]}=$tmp[9];
    $skin{$tmp[0]}=$tmp[10];
    $race{$tmp[0]}=$tmp[11];
    for(my $j=22; $j<38; $j++){
	$lipid{$tmp[0]}{$j}=$tmp[$j];
    } 
}
close A;

open (B, "<$file")||die "could not open $file\n";
open (C, ">$file.lipid")||die "could not open $file.lipid\n";
my $x=<B>;
chomp $x;
chomp $x;for($x){s/\r//gi;}
my @aa=split /\t/, $x;
my $aaa="age";
my $cutS="cut -f1";
for(my $i=1; $i<@aa; $i++){
    my @tmp=split /\./, $aa[$i];
    my $oldid=$tmp[5]."-".$tmp[6];
    print STDERR $oldid,"\n";
    if(defined $age{$oldid}){
        $aaa=$aaa."\t".$age{$oldid};
        if($age{$oldid} ne "NA"){
            $cutS=$cutS.",".($i+1);
        }
    }else{
	print STDERR "could not match age  $oldid\n";
        $aaa=$aaa."\tNA";       
    }

    for(my $j=22; $j<38; $j++){
	if(defined $lipid{$oldid}{$j}){
	    $lipidT{$j}=$lipidT{$j}."\t".$lipid{$oldid}{$j};
	    if( $lipid{$oldid}{$j} ne "NA"){
		$cutC{$j}=$cutC{$j}.",".($i+1);
	    }
	}else{
	    print STDERR "could not match $j  $oldid\n";
	    $lipidT{$j}=$lipidT{$j}."\tNA";	
	}
    }
}

print C "$x\n";
for(my $j=22; $j<38; $j++){
    print C $lipidT{$j}, "\n";
}  

print C $aaa, "\n";
while (my $b=<B>){
    print C $b;
}
close B;
close C;
for(my $j=22; $j<38; $j++){
    print($cutC{$j}, " $file.lipid >$file.lipid.clean\n");
}

