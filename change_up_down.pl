#! /usr/bin/perl -w
##########################################################################
#### Author: Ping HU
#### need to remove all " in the file to match correctly
##########################################################################
use strict;
#my $usage="$0  [file] first line is title\n"; 
#die $usage unless @ARGV ==1;
#my ($list) = @ARGV;
my $list="all_Microprofile.xls";
#my $list ="test4";
open (A, "<$list")||die "could not open $list\n";
my $title=<A>;
my %content;
my %number;
my %Dir;
my %Kname;
my %Fname;
my %totalCount;

# Assuming open filehandle A beforehand, and declaration of %Kname, %Fname, %Dir, %number, %content
while (my $a = <A>) {
    chomp $a;
    # Removing carriage return characters; no need to do it twice
    $a =~ s/\r//g;
    
    # Splitting the line into fields based on tab delimiter
    my @tmp = split /\t/, $a;
    
    # Checking if there are less than 6 fields, if so, print to STDERR and skip to the next iteration
    if (scalar @tmp < 6) {
        print STDERR $a, "\n";
        next;
    }
    
    # Constructing various strings from parts of the line
    my $file = $tmp[0] . "." . $tmp[1] . "." . $tmp[2];
    my $ko_name = $tmp[4] . "|" . $tmp[5];
    my $count = $tmp[12];
    my $direction = $tmp[13];
    my $p = $tmp[8];
    my $pa = $tmp[9];
    
    # Populating hashes
    $Kname{$ko_name} = 1;
    $Fname{$file} = 1;
    print STDERR $ko_name, $file, $direction, $count, "\n";
    if($direction eq "All"){
	$totalCount{$file}{$ko_name}=$count;
    }
    if (!defined $Dir{$file}{$ko_name}) {
        # Initial population if not already defined
        $Dir{$file}{$ko_name} = $direction;
        $number{$file}{$ko_name} = $count;
        $content{$file}{$ko_name} = $a;
        print STDERR "1\n";
    }elsif(($Dir{$file}{$ko_name} eq "All")&&($direction ne "All")&&($count*2 != $number{$file}{$ko_name})){
	# Logic when direction is "All"; assuming this is a special case
	$Dir{$file}{$ko_name} = $direction;
	$number{$file}{$ko_name} = $count;
	$content{$file}{$ko_name} = $a;
	print STDERR "2\n";
    }elsif (($Dir{$file}{$ko_name} ne "All") &&($direction ne "All")&&($count == $number{$file}{$ko_name})){
	# Update if current count is greater than stored count
	$Dir{$file}{$ko_name} = "All";
	$number{$file}{$ko_name} = $count;
	$content{$file}{$ko_name} = $a;
	print STDERR "3\n";
    }elsif (($Dir{$file}{$ko_name} ne "All") &&($direction ne "All")&&($count > $number{$file}{$ko_name})){
	# Update if current count is greater than stored count
	$Dir{$file}{$ko_name} = $direction;
	$number{$file}{$ko_name} = $count;
	$content{$file}{$ko_name} = $a;
	print STDERR "3\n";
    }
        
}

close A;
print $title;
foreach my $i (sort keys %Fname){
    foreach my $j (sort keys %Kname){
	if(defined $content{$i}{$j}){
	    print $content{$i}{$j}, "\n";
	}
    }
}
