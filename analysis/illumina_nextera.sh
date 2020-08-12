#!/bin/bash

##datadir=/uru/Data/Nanopore/projects/ncov/nextera
##rawdir=/uru/Data/NGS/Raw/200722_ncov

datadir=/uru/Data/Nanopore/projects/ncov/nextera_v2
rawdir=/uru/Data/NGS/Raw/200807_ncov_nextflex 

if [ $1 == quick_alignment ] ; then
    mkdir -p $datadir/align

    ##reference included in artic
    ref=~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta

    bwa index $ref

    for i in $rawdir/*_L001_R1_001.fastq.gz  ;
    do
        prefix=`basename $i _L001_R1_001.fastq.gz `

        bwa mem -t 36 \
            $ref \
            $rawdir/${prefix}_L001_R1_001.fastq.gz \
            $rawdir/${prefix}_L001_R2_001.fastq.gz | \
            samtools view -@ 36 -bS - | \
            samtools sort -@ 36 -o $datadir/align/$prefix.sorted.bam

        samtools index $datadir/align/$prefix.sorted.bam
    done
fi

if [ $1 == genomecov ] ; then
   mkdir -p $datadir/cov

   for i in $rawdir/*_L001_R1_001.fastq.gz ;
   do
       prefix=`basename $i _L001_R1_001.fastq.gz `
       bedtools genomecov -d -ibam $datadir/align/$prefix.sorted.bam > $datadir/cov/$prefix.cov &
   done
fi
   
if [ $1 == trimmomatic ] ; then
    ##stick with trimmo for getting rid of illumina adapter for now
    mkdir -p $datadir/trimmed
    cp ~/software/Trimmomatic-0.39/adapters/all_adapters.fa ~/Code/ncov/analysis/

    for i in $rawdir/*R1*fastq.gz ;
    do
        prefix=`basename $i _L001_R1_001.fastq.gz`
        echo $prefix
        java -jar ~/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 36 -phred33 \
             $rawdir/${prefix}_L001_R1_001.fastq.gz $rawdir/${prefix}_L001_R2_001.fastq.gz \
             $datadir/trimmed/${prefix}_fwd_paired.fq.gz $datadir/trimmed/${prefix}_fwd_unpaired.fq.gz \
             $datadir/trimmed/${prefix}_rev_paired.fq.gz $datadir/trimmed/${prefix}_rev_unpaired.fq.gz \
             ILLUMINACLIP:all_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36
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
            samtools sort -@ 36 -o $datadir/align_trimmed/$prefix.sorted.bam
        samtools index $datadir/align_trimmed/$prefix.sorted.bam
    done
fi

if [ $1 == trim_pairs ] ; then
    bed=~/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed

    for i in $datadir/align_trimmed/*.sorted.bam ;
    do
        prefix=`basename $i .sorted.bam`
        align_trim \
            --remove-incorrect-pairs \
            --report $datadir/align_trimmed/$prefix.report_pairs.txt \
            $bed \
            < $datadir/align_trimmed/$prefix.sorted.bam | \
            samtools sort -@ 36 -T $prefix.tmp -o $datadir/align_trimmed/$prefix.sorted.trimmed_pairs.bam

        samtools index $datadir/align_trimmed/$prefix.sorted.trimmed_pairs.bam
    done
fi


if [ $1 == trimmed_genomecov ] ; then
   mkdir -p $datadir/cov

   for i in $datadir/align_trimmed/*.sorted.trimmed_pairs.bam ;	    
   do
       prefix=`basename $i .sorted.trimmed_pairs.bam `
       bedtools genomecov -d -ibam $datadir/align_trimmed/$prefix.sorted.trimmed_pairs.bam > $datadir/cov/trimmed/$prefix.trimmed.cov
   done
fi
