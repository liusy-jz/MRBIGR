#!/usr/bin/perl
use strict;
use warnings;

die "perl $0 <pheno.csv> <hap.tsv> <qtl.anno.csv> <outdir>" unless @ARGV==4;

my $outdir=$ARGV[3];
`mkdir -p $outdir`;
open RST, ">$outdir/results.tsv" or die $!;

open IN, "$ARGV[0]" or die $!;
my $head=<IN>;
chomp($head);
my @head=split /,/, $head;
my %hash;
while(<IN>){
	chomp;
	my @tmp=split /,/;
	for(my $i=1; $i<=$#tmp; $i++){
		my $t=$head[$i];
		$hash{$t}{$tmp[0]}=$tmp[$i];
	}
}
close IN;

my %hap;
open IN, "$ARGV[1]" or die $!;
while(<IN>){
        chomp;
        my @tmp=split /\t/;
	my @gene=split /,/, $tmp[0];
	foreach my $g(@gene){
		$hap{$g}=$_;
	}
}
close IN;	

open IN, "$ARGV[2]" or die $!;
<IN>;
while(<IN>){
        chomp;
	my $seq=$_;
        my @tmp=split /,/;
	my $trait=$tmp[7];
	my $genes=$tmp[8];
	my @gene=split /;/, $genes;
	foreach my $gene(@gene){
		next unless exists($hap{$gene});
		my @haplo=split /\t/, $hap{$gene};
		my $haplo_info="$haplo[1]\t$haplo[2]\t$haplo[3]";
		my $haplo_detail='';
	
	
		my $otfile="$outdir/$trait.$gene.ttest.csv";
		my $outpfile="$outdir/$trait.$gene.pval.txt";	
		open OUT, ">$otfile" or die $!;
		for(my $i=4; $i<=$#haplo; $i++){
			my $haplow=$haplo[$i];
			my @haplow=split /\|/, $haplow;
			my $hap_name=$haplow[0];
			my $hap_num=$haplow[1];
			my $hap_type=$haplow[2];
			my $sample=$haplow[3];
			my @samples=split /,/, $sample;
			$haplo_detail .= '|'.$hap_name.'|'.$hap_num.'|'.$hap_type;

			foreach my $s(@samples){
				if(exists $hash{$trait}{$s}){
					print OUT "$s,$hash{$trait}{$s},$hap_name\n";
				}
			}
		}
	
		close OUT;
		$haplo_detail =~s/^\|//;

my $otfileprefix=$otfile;
$otfileprefix=~s/\.csv$//;
my $ttest=<<"END";
pwd<-getwd()
setwd(pwd)
data<-read.csv(\"$otfile\",head=F)
df<-data.frame(Value=data[,2],Group=factor(data[,3]))
len<-length(levels(df\$Group))
#Welch’s t-test
if(len==2){
	a1<-t.test(df\$Value~df\$Group, df)
	pval<-a1\$p.value
	write.table(pval,\"$outpfile\",col.name=F,quote=F)
}
#Tukey’s test
if(len>2){
	a1<-aov(df\$Value~df\$Group)
	pval<-TukeyHSD(x=a1,\'df\$Group\', conf.level=0.95)\$\`df\$Group\`[,4]
	write.table(pval,\"$outpfile\",col.name=F,quote=F)
}

library(ggplot2)
require(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, \"Set1\"))
p<-ggplot(df, aes(Group, Value,colour=Group,fill=NA))  + geom_jitter(size=0.5, colour="grey")+geom_boxplot(outlier.colour=NA) + labs(title=\"$trait $gene\") + theme(panel.background = element_rect(fill=\'white\', colour=\'black\'),panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + scale_fill_manual(values = getPalette(10)) 
p<-p + theme_bw()+ theme(legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size = 4),
          legend.key.size = unit(3,'mm'),
          legend.position = 'none',
          panel.grid = element_blank(),
          axis.line = element_line(colour = 'black',size=0.4),
          axis.text = element_text(size = 6,color = 'black'),
          axis.ticks.length=unit(.1, 'cm'))
pdf(\"$otfileprefix.pdf\",width=3.5,height=4.5)
p
dev.off()

END
	
		open R,"|R --vanilla --slave" or die $!;
		print R $ttest;
		close R;
	
		my $minP=1;
		my $pval='';
		open PVAL, "$outpfile" or die $!;
		while(<PVAL>){
			chomp;
			my @tmp=split /\s+/;
			if($tmp[0] eq '1'){
				$pval=$tmp[1];
				$minP=$tmp[1];
			}
			else {
				$pval .= ','.$tmp[0].":".$tmp[1];
				$minP=$tmp[1] if($tmp[1] < $minP);
			}
		}
		close PVAL;
		$pval =~s/^,//;
		print RST "$trait\t$gene\t$haplo_info\t$haplo_detail\t$pval\t$minP\n";
	}
}
close IN;	
close RST;
		

