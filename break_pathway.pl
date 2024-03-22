#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### need to remove all " in the file to match correctly
##########################################################################
use strict;
#my $usage="$0  [file] first line is title\n"; 
#die $usage unless @ARGV ==1;
#my ($list) = @ARGV;
#my $list="compiled_microprofile.xls";
my $list ="pstat_microprofile";
#my $list ="test4";
open (A, "<$list")||die "could not open $list\n";
my $title=<A>;
my %content;



while (my $a = <A>) {
    chomp $a;
    $a =~ s/\r//g;
    my @tmp = split /\t/, $a;
    
    if (scalar @tmp < 6) {
        print STDERR $a, "\n";
        next;
    }

    my $file = $tmp[0] . "." . $tmp[1] . "." . $tmp[2];
    #my $ko_name = $tmp[4] . "|" . $tmp[5];
    my $ko_name = $tmp[4];
    my $count = $tmp[12];
    my $direction = $tmp[13];
    my $p = $tmp[8];
    my $pa = $tmp[9];
    $content{$ko_name}{$a}=$tmp[4] . "|" . $tmp[5];
        
}

close A;

foreach my $i (sort keys %content){
    open (B, ">$i.koenrich")||die "could not open $i.koenrich\n";
    print B $title;
    foreach my $j (sort keys %{$content{$i}}){
	if(defined $content{$i}{$j}){
	    print B $j, "\n";
	}
    }
    close B;
}
