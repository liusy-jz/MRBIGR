#!/usr/bin/perl
use strict;
use warnings;

my %hash;
open IN, "$ARGV[0]" or die $!;
while(<IN>){
	chomp;
	my @tmp=split /,/;
	my @gene=split /;/, $tmp[8];
	foreach my $g(@gene){
		$hash{$g}='';
	}
}
close IN;

my %ctrl;
open IN, "$ARGV[1]" or die $!;
while(<IN>){
        chomp;
	my @tmp=split /\t/;
	my @temp=split /,/, $tmp[6];
	$tmp[0] =~ s/Chr//;
	$tmp[0] =~ s/chr//;
	my $ctrl="chr".$tmp[0].".s_".$tmp[1];
	
	foreach my $e(@temp){
		if(exists $hash{$e}){
			unless(exists $ctrl{$ctrl}){
				print "$_\n";
				$ctrl{$ctrl}= '';
			}
		}
	}
}
close IN;

