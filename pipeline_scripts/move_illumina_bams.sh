#!/bin/bash
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus"

REF="/home/idies/workspace/covid19/ncov_reference/sequence.fasta"
OUTDIR="/home/idies/workspace/covid19/curation/*$RUN/bams"
OUTDIR=`echo $OUTDIR`
OUTDIR=$OUTDIR/illumina
if [ ! -d $OUTDIR ]
then
    mkdir $OUTDIR
fi
for bamfile in `ls $DIR/*.nanopolish.primertrimmed.rg.sorted.bam`; do

    sample=${bamfile##*/}
    samplename=${sample%%_*}

    if [ ! "$samplename" = "NTC" ]; then

        echo $samplename

        prefix=`echo $bamfile`
        prefix=${prefix##*/}
        prefix=${prefix%%.*}
	
	samplemap=$BINDIR/samplenamemap.txt
	
	illuminasample=`cat $samplemap | grep $samplename | cut -f 1`
	if [ -z $illuminasample ]
        then
            continue
        fi

	echo 'Illumina sample: '$illuminasample
	illuminabam=`ls /home/idies/workspace/covid19/illumina/*/trimmedvariantcalling/$illuminasample'_'*'.norm200.bam' 2>/dev/null`
	
        if [ -z $illuminabam ]
        then
            continue
        fi
        echo 'Illumina bam: '$illuminabam
	samplewithbarcode=${sample%%.*}
        OUTPATH=$OUTDIR/$samplewithbarcode.illumina.bam
	if [ ! -r $OUTPATH ]
        then
	    ln -s $illuminabam $OUTPATH
	    ln -s $illuminabam'.bai' $OUTPATH'.bai'
        fi
    fi
done

