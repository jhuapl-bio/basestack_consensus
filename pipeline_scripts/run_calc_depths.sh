#!/bin/bash

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus"

for consfile in $DIR/*nanopolish.consensus.fasta; do

	sample=${consfile##*/}
	samplename=${sample%%.*}

	# loop through all samples including NTC
	echo $samplename

	bamfile=$DIR/$samplename.nanopolish.primertrimmed.rg.sorted.bam
	outfile="$DIR/$samplename.nanopolish.primertrimmed.rg.sorted.del.depth"
	
	# run script
	if [ ! -f $outfile ]; then
		python /home/idies/workspace/covid19/code/ncov/pipeline_scripts/calc_sample_depths.py $bamfile $outfile
	fi

done