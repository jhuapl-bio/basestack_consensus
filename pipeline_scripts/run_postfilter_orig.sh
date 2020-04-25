#!/bin/bash

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200-orig/"
NTC="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200-orig/NTC/"

# make and save output directory
outdir="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-post-filter-orig/"

# save path to NTC bamfile
ntc_bamfile="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/3-hac-medaka-norm200-orig/NTC/*primertrimmed.rg.sorted.bam"

# save path to nextstrain vcf
vcf_next="/home/idies/workspace/covid19/nextstrain/archive/2020-04-14/alignments.vcf"

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
		python /home/idies/workspace/Storage/swohl/persistent/covid/vcf_postfilter.py --vcffile $vcffile --bamfile $bamfile --consensus $consensus --ntc-bamfile $ntc_bamfile --vcf-nextstrain $vcf_next --ns-snp-threshold 2 -o $outdir --prefix $prefix
	fi
done

