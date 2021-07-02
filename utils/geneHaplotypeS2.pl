#!/usr/bin/perl
use strict;
use warnings;

open IN, "$ARGV[0]" or die $!;
my %geno;
while(<IN>){
        chomp;
        my @tmp=split /\t/;
	$tmp[0]=~ s/Chr//;
	$tmp[0]=~ s/chr//;
	my $id="chr".$tmp[0].".s_".$tmp[1];
	$geno{$id}=$tmp[3]."/".$tmp[4];
}
close IN;

open IN, "$ARGV[1]" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /\t/, $head;
shift(@head);
shift(@head);

while(<IN>){
	chomp;
	my @tmp=split /\t/;
	my %hash=();
	my %count=();
	for(my $i=2; $i<=$#tmp; $i++){
		my $k=$tmp[$i];
		next if($k =~ /NA/);
		my $s=$head[$i-2];
		if(exists $hash{$k}){
			$hash{$k} .= ','.$s;
		}
		else{
			$hash{$k}=$s;
		}	
	}
	
	foreach my $k(keys %hash){
		my @test=split /,/,$hash{$k};
		$count{$k}=@test;
		if(@test < 5){
			delete($hash{$k});
			delete($count{$k});
		}
	}
        my $n=keys %hash;
	next if($n >= 9 ||  $n < 2);
	my @ss=split /,/, $tmp[1];
	
	my $ss='';
	foreach my $s(@ss){
		$ss .= ','.$geno{$s};
	}
	$ss=~s/^,//;
	print "$tmp[0]\t$tmp[1]\t$ss\t$n";

	my $j=1;
	
	foreach my $k(sort keys %hash){
		
		print "\tHap$j|$count{$k}|$k|$hash{$k}";
		$j++;
	}
	print "\n";
}
close IN;

