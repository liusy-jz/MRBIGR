#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my %gene;
while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my $chr=$tmp[0];
	my $pos=$tmp[1];
	my $genes=$tmp[6];
	$gene{$chr}{$pos}=$genes;
}
close IN;

open IN, "$ARGV[1]" or die $!;
while(<IN>){
	chomp;
	next if(/^##/);
	my @tmp=split /\t/;
	next if($tmp[4] =~ /,/);
	if(/^#/){
		print "#CHROM\tPOS";
		for(my $i=9; $i<=$#tmp; $i++){
			print "\t$tmp[$i]";
		}
		print "\n";
		next;
	}
	if(exists $gene{$tmp[0]}{$tmp[1]}){
		my $chr= $tmp[0];
		$chr=~ s/Chr//;
		$chr=~ s/chr//;
		my $genes = $gene{$tmp[0]}{$tmp[1]};
		my @gene = split /,/, $genes;
		foreach my $g(@gene){
			print "$g\tchr"."$chr".".s_"."$tmp[1]";
			for(my $i=9; $i<=$#tmp; $i++){
				my $geno=(split(':',$tmp[$i]))[0];
				my $num='NA';
				if($geno=~/\d/){
					$num = (split('/',$geno))[0]+(split('/',$geno))[1];
				}
				print "\t$num";
			}
			print "\n";
		}
	}
}
close IN;
