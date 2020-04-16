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
	
if [ $1 == demux_v3 ] ; then
    calldir=~/data/ncov/initial-clinical/v3_run1/20200405_0422_GA30000_FAN30842_00fb1614

    mkdir -p $calldir/barcodes
    
    guppy_barcoder \
	-i $calldir/fastq_pass \
	-s $calldir/guppy_demux \
	-t 4 \
	-q 0 \
	--trim_barcodes \
	--barcode_kits "EXP-NBD104" \
	-x "cuda:0"
fi

if [ $1 == demux_count ] ; then
    calldir=~/data/ncov/initial-clinical/v3_run1/20200405_0422_GA30000_FAN30842_00fb1614
    touch $calldir/guppy_demux/counts.csv
    for i in 01 02 03 04 05 06 07 08 09 10 11 12 ;
    do
	numlines=`wc -l $calldir/guppy_demux/barcode$i/*fastq | cut -d ' ' -f 1 `
	numreads=`expr $numlines / 4`
	echo barcode$i,$numreads >> $calldir/guppy_demux/counts.csv
    done
    uclines=`wc -l $calldir/guppy_demux/unclassified/*fastq | cut -d ' ' -f 1 `
    ucreads=`expr $uclines / 4`
    echo unclassified,$ucreads >> $calldir/guppy_demux/counts.csv
fi

if [ $1 == dmux_test ] ; then
    calldir=~/data/ncov/initial-clinical

    for i in 20200405_0422_GA30000_FAN30842_00fb1614  20200410_2018_X4_FAN32204_327837a0 ;
    do
	guppy_barcoder \
	    -i $calldir/$i/fastq_pass \
	    -s $calldir/$i/guppy_demux \
	    -t 4 \
	    -q 0 \
	    --trim_barcodes \
	    --barcode_kits "EXP-NBD104" \
	    -x "cuda:0"

	guppy_barcoder \
	    -i $calldir/$i/fastq_pass \
	    -s $calldir/$i/guppy_demux_both \
	    -t 4 \
	    -q 0 \
	    --trim_barcodes \
	    --require_barcodes_both_ends \
	    --barcode_kits "EXP-NBD104" \
	    -x "cuda:0"

    done
fi

if [ $1 == dmux_test_count ] ; then
    calldir=~/data/ncov/initial-clinical

    for run in 20200405_0422_GA30000_FAN30842_00fb1614  20200410_2018_X4_FAN32204_327837a0 ;
    do
	for mode in guppy_demux guppy_demux_both
	do
	    touch $calldir/$run/$mode/counts.csv
	    for i in 01 02 03 04 05 06 07 08 09 10 11 12 ;
	    do
		numlines=`wc -l $calldir/$run/$mode/barcode$i/*fastq | cut -d ' ' -f 1 `
		numreads=`expr $numlines / 4`
		echo barcode$i,$numreads >> $calldir/$run/$mode/counts.csv
	    done
	    uclines=`wc -l $calldir/$run/$mode/unclassified/*fastq | cut -d ' ' -f 1 `
	    ucreads=`expr $uclines / 4`
	    echo unclassified,$ucreads >> $calldir/$run/$mode/counts.csv
	done
    done
fi
