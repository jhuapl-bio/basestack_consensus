#!/bin/bash

rawdir=/uru/Data/Nanopore/projects/ncov/dRNA/raw
datadir=~/data/ncov

if [ $1 == untar ] ; then
    for i in HP4 HP20 ;
    do
	mkdir -p $datadir/$i
	mkdir -p $datadir/$i/raw
	tar -xzf $rawdir/200601_ncov_$i.tar.gz -C $datadir/$i/raw/
    done
fi

	
if [ $1 == call ] ; then
    for i in HP4 HP20 ;
    do
	mkdir -p $datadir/$i/called

	guppy_basecaller -r \
	    -i $datadir/$i/raw \
	    -s $datadir/$i/called \
	    --flowcell FLO-MIN106 \
	    --kit SQK-RNA002 \
	    --device "cuda:0"
    done
fi
	
if [ $1 == fq ] ; then
    for i in HP4 HP20 ;
    do
	mkdir -p $datadir/$i/fastqs

	cat $datadir/$i/called/*fastq > $datadir/$i/fastqs/$i.fq
    done
fi

if [ $1 == align ] ; then
    ref=~/data/ncov/ref/nCoV-2019.reference.fasta    
    for i in HP4 HP20 ;
    do
	mkdir -p $datadir/$i/align

	minimap2 \
	    -ax splice -k14 -uf -t 36 \
	    $ref \
	    $datadir/$i/fastqs/$i.fq |\
	    samtools view -@ 36 -b |\
	    samtools sort -@ 36 -o $datadir/$i/align/$i.sorted.bam -T $datadir/$i/align/reads.tmp -

	samtools index $datadir/$i/align/$i.sorted.bam
    done
fi

	    
if [ $1 == UT ] ; then
    for i in HP4 HP20 ;
    do
	seqkit seq --rna2dna $datadir/$i/fastqs/$i.fq > $datadir/$i/fastqs/${i}_dna.fq
    done
fi
