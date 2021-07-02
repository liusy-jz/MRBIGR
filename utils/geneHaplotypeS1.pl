#!/usr/bin/perl
use strict;
use warnings;


my %hash;
open IN, "$ARGV[0]" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /\t/, $head;
shift(@head);
shift(@head);

while(<IN>){
	chomp;
	my @tmp=split /\t/;
	if(exists $hash{$tmp[0]}{'pos'}){
		$hash{$tmp[0]}{'pos'} .= ','.$tmp[1];
	}
	else{
		$hash{$tmp[0]}{'pos'}=$tmp[1];
	}
	for(my $i=0; $i<=$#head;$i++){
		my $s=$head[$i];
		if(exists $hash{$tmp[0]}{$s}){
			$hash{$tmp[0]}{$s} .= ','.$tmp[$i+2];
		}
		else{
			$hash{$tmp[0]}{$s} = $tmp[$i+2];
		}
	}
}
close IN;

print "$head\n";
foreach my $k(sort keys %hash){
	print "$k\t$hash{$k}{'pos'}";
	foreach my $k2(@head){
		print "\t$hash{$k}{$k2}";
	}
	print "\n";
	close IN;
}
close IN;
