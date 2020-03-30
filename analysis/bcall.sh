#!/bin/bash

datadir=/uru/Data/Nanopore/projects/ncov

##needs basecalling docker:
##docker run --runtime=nvidia --name yfan -i -t -v ~/data/ncov:/media yfan2012/whack:guppy_bcall /bin/bash
if [ $1 == init_accurate ] ; then
    ##defaults to hac (high accuracy)

    ##check if the fast5 data is there
    rawdir=/media/initial-clinical/$2/fast5_pass
    if [ -d $rawdir ] ; then
	##call
	calldir=/media/initial-clinical/$2/called_accurate 
	mkdir -p $calldir
	guppy_basecaller -i $rawdir -s $calldir --flowcell FLO-MIN106 --kit SQK-LSK109 -x "cuda:0"
    fi
    
fi

if [ $1 == demux ] ; then
    ##cat stuff into one file
    calldir=/media/initial-clinical/$2/called_accurate

    guppy_barcoder -i /media/initial-clinical/$2/called_accurate \
		   -s /media/initial-clinical/$2/guppy_demux \
		   -t 4 \
		   --barcode_kits "EXP-PBC001" \
		   --trim_barcodes \
		   --compress_fastq \
    		   -x "cuda:0" 
fi

if [ $1 == demux_single ] ; then
    ##cat stuff into one file
    calldir=/media/initial-clinical/$2/called_accurate

    guppy_barcoder -i /media/initial-clinical/$2/called_accurate \
		   -s /media/initial-clinical/$2/guppy_demux_front \
		   -t 4 \
		   --barcode_kits "EXP-PBC001" \
		   --trim_barcodes \
		   --min_score 60 \
		   --min_score_rear_override 0 \
		   --compress_fastq \
    		   -x "cuda:0"
    guppy_barcoder -i /media/initial-clinical/$2/called_accurate \
		   -s /media/initial-clinical/$2/guppy_demux_rear \
		   -t 4 \
		   --barcode_kits "EXP-PBC001" \
		   --trim_barcodes \
		   --min_score 0 \
		   --min_score_rear_override 60 \
		   --compress_fastq \
    		   -x "cuda:0" 
fi

    
if [ $1 == cat_barcodes ] ; then
    mkdir -p $datadir/initial-clinical/$2/guppy_fastqs
    for i in 01 02 03 04 05 06 07 08 09 10 11 12 ;
    do
	cat $datadir/initial-clinical/$2/guppy_demux_rear/barcode$i/*fastq.gz > $datadir/initial-clinical/$2/guppy_fastqs/barcode${i}_rear.fq.gz
	cat $datadir/initial-clinical/$2/guppy_demux_front/barcode$i/*fastq.gz > $datadir/initial-clinical/$2/guppy_fastqs/barcode${i}_front.fq.gz
    done
fi
	
