#!/bin/bash

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200/"
NTC="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200/NTC/"

# make and save output directory
if [ ! -d /home/idies/workspace/Storage/swohl/persistent/covid/$RUN ]; then
	mkdir /home/idies/workspace/Storage/swohl/persistent/covid/$RUN
fi
outdir="/home/idies/workspace/Storage/swohl/persistent/covid/$RUN/"

# save path to NTC bamfile
ntc_bamfile="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200/NTC/*primertrimmed.rg.sorted.bam"

for sampledir in $DIR*/; do

	sampledir=${sampledir%/} # remove trailing slash
	samplename=${sampledir##*/}

	# loop through all NTC samples
	if [ ! "$samplename" = "NTC" ]; then

		echo $samplename

		vcffile="$sampledir/*.pass.vcf.gz"
		bamfile="$sampledir/*.primertrimmed.rg.sorted.bam"
		consensus="$sampledir/*.consensus.fasta"
		prefix=`echo $consensus`
		prefix=${prefix##*/}
		prefix=${prefix%%.*}

		# run script
		python /home/idies/workspace/Storage/swohl/persistent/covid/vcf_postfilter.py --vcffile $vcffile --bamfile $bamfile --consensus $consensus --ntc-bamfile $ntc_bamfile -o $outdir --prefix $prefix
	fi
done

