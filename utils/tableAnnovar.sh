#!/bin/bash

vcf=$1
dbdir=$2
dbprefix=$3
outprefix=$4

table_annovar.pl $vcf $dbdir -buildver $dbprefix -out $outprefix -remove -protocol refGene -operation g -nastring . -vcfinput
cut -f 1-10 $outprefix.$dbprefix\_multianno.txt >$outprefix.$dbprefix\_multianno.bed
rm $outprefix.avinput $outprefix.$dbprefix\_multianno.txt

cat $outprefix.$dbprefix\_multianno.bed|grep -v 'unknown'|grep -v 'UTR'|grep -v 'ncRNA'|grep -v 'stream'|grep -v 'intergenic'|grep -v 'intronic'|grep -v 'nonframeshift'|grep -vP '\tsynonymous SNV' >$outprefix.$dbprefix\_largeEffect.bed

