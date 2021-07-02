#!/bin/bash

gtf=$1
fa=$2
outprefix=$3

gtfToGenePred -genePredExt $gtf $outprefix\_refGene.txt
retrieve_seq_from_fasta.pl --format refGene --seqfile $fa $outprefix\_refGene.txt --out $outprefix\_refGeneMrna.fa



