#!/bin/bash

##datadir=/uru/Data/Nanopore/projects/ncov/illumina
##rawdir=/uru/Data/NGS/Raw/200421_ncov

datadir=/uru/Data/Nanopore/projects/ncov/nextera_v2
rawdir=/uru/Data/NGS/Raw/200807_ncov_nextflex

if [ $1 == trimmomatic ] ; then
    ##stick with trimmo for getting rid of illumina adapter for now

    mkdir -p $datadir/trimmed
    cat ~/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa ~/software/Trimmomatic-0.39/adapters/TruSeq2-PE.fa > ~/Code/ncov/analysis/TruSeq_all-PE.fa
    
    for i in $rawdir/*R1*fastq.gz ;
    do
	prefix=`basename $i _L001_R1_001.fastq.gz`
	echo $prefix
	java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
	     $rawdir/${prefix}_L001_R1_001.fastq.gz $rawdir/${prefix}_L001_R2_001.fastq.gz \
	     $datadir/trimmed/${prefix}_fwd_paired.fq.gz $datadir/trimmed/${prefix}_fwd_unpaired.fq.gz \
	     $datadir/trimmed/${prefix}_rev_paired.fq.gz $datadir/trimmed/${prefix}_rev_unpaired.fq.gz \
	     ILLUMINACLIP:TruSeq_all-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
    done
fi


    
if [ $1 == align ] ; then
    mkdir -p $datadir/align_trimmed

    ##aligning to artic network reference for trimming step
    ref=~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta

    bwa index $ref
    
    for i in $datadir/trimmed/*_fwd_paired.fq.gz ;
    do
	prefix=`basename $i _fwd_paired.fq.gz`
	
	bwa mem -t 36 \
	    $ref \
	    $datadir/trimmed/${prefix}_fwd_paired.fq.gz \
	    $datadir/trimmed/${prefix}_rev_paired.fq.gz | \
	    samtools view -@ 36 -bS - | \
	    samtools sort -@ 36 -T $prefix -o $datadir/align_trimmed/$prefix.sorted.bam

	samtools index $datadir/align_trimmed/$prefix.sorted.bam
    done
fi
	    
if [ $1 == trim_primers ] ; then
    bed=~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed
    
    for i in $datadir/align_trimmed/*.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`

	align_trim \
	    --start \
	    --report $datadir/align_trimmed/$prefix.report.txt \
	    $bed \
	    < $i | \
	    samtools sort -@ 36 -T $prefix.tmp -o $datadir/align_trimmed/$prefix.trimmed.sorted.bam

	samtools index $datadir/align_trimmed/$prefix.trimmed.sorted.bam
    done
fi

if [ $1 == trim_ends ] ; then
    bed=~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed
    
    for i in $datadir/align_trimmed/*.sorted.bam ;
    do
	prefix=`basename $i .sorted.bam`

	align_trim \
	    --report $datadir/align_trimmed/$prefix.report_end.txt \
	    $bed \
	    < $i | \
	    samtools sort -@ 36 -T $prefix.tmp -o $datadir/align_trimmed/$prefix.trimmed_end.sorted.bam

	samtools index $datadir/align_trimmed/$prefix.trimmed_end.sorted.bam
    done
fi

	    
if [ $1 == trim_pairs ] ; then
    bed=~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed

    for i in $datadir/align_trimmed/*report.txt ;
    do
	prefix=`basename $i .report.txt`
	align_trim \
	    --remove-incorrect-pairs \
	    --report $datadir/align_trimmed/$prefix.report_pairs.txt \
	    $bed \
	    < $datadir/align_trimmed/$prefix.sorted.bam | \
	    samtools sort -@ 36 -T $prefix.tmp -o $datadir/align_trimmed/$prefix.trimmed_pairs.sorted.bam
	
	samtools index $datadir/align_trimmed/$prefix.trimmed_pairs.sorted.bam
    done
fi

	     

    
