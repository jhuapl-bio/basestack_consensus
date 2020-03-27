#!/bin/bash

datadir=/uru/Data/Nanopore/projects/ncov

##needs basecalling docker:
##docker run --runtime=nvidia --name yfan -i -t -v ~/data/ncov:/media yfan2012/whack:guppy_bcall /bin/bash
if [ $1 == init_accurate ] ; then
    ##defaults to hac (high accuracy)
    for i in 4 5 6 8 ;
    do
	##check if the fast5 data is there
	rawdir=/media/initial-clinical/JHU-0000$i/fast5_pass
	if [ -d $rawdir ] ; then

	    ##call
	    calldir=/media/initial-clinical/JHU-0000$i/called_accurate 
	    mkdir -p $calldir
	    guppy_basecaller -i $rawdir -s $calldir --flowcell FLO-MIN106 --kit SQK-LSK109 -x "cuda:0"
	fi
    done
fi
