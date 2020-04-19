#!/bin/bash

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200"

# make and save output directory
if [ ! -d /home/idies/workspace/Storage/swohl/persistent/covid/$RUN ]; then
	mkdir /home/idies/workspace/Storage/swohl/persistent/covid/$RUN
fi
outdir="/home/idies/workspace/Storage/swohl/persistent/covid/$RUN/"

# save path to NTC bamfile
ntc_bamfile="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200/NTC*primertrimmed.rg.sorted.bam"

for consfile in $DIR/*.consensus.fasta; do

	sample=${consfile##*/}
	samplename=${sample%%_*}

	# loop through all NTC samples
	if [ ! "$samplename" = "NTC" ]; then

		echo $samplename

		vcffile="$DIR/$samplename*.pass.vcf.gz"
		bamfile="$DIR/$samplename*.primertrimmed.rg.sorted.bam"
		consensus="$DIR/$samplename*.consensus.fasta"
		prefix=`echo $consensus`
		prefix=${prefix##*/}
		prefix=${prefix%%.*}

		# run script
		python /home/idies/workspace/Storage/swohl/persistent/covid/vcf_postfilter.py --vcffile $vcffile --bamfile $bamfile --consensus $consensus --ntc-bamfile $ntc_bamfile -o $outdir --prefix $prefix
	fi
done

