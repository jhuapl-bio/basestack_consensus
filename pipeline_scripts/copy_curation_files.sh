#!/bin/bash

# specify a sequencing run name
RUN=$1

# specify the main destination directory
OUT="/home/idies/workspace/covid19/curation"

# get the report name
# and extract the plate name prefix
report="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/run_stats/plate*report.pdf"
prefix=`echo $report`
prefix=${prefix##*/}
prefix=$(echo $prefix | cut -f1,2 -d'-')

# get the full subdirectory name including the plate information
fullname=$(echo $prefix-$RUN)
DIR="/home/idies/workspace/covid19/curation/$fullname"

# make this curation run directory if it does not exist
if [ ! -d $DIR ]; then
	mkdir $DIR
	mkdir $DIR/bams $DIR/bams/ont
fi

# path to postfilter output
postfilt="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter/postfilt_all.txt"

# path to ont bam files
bamdir="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus"

# copy in all the files
cp $report $DIR
cp $postfilt $DIR/$prefix-$RUN-variants.txt

for bamfile in `ls $bamdir/*nanopolish.primertrimmed.rg.sorted.ba*`; do
	bamname=${bamfile##*/}
	if [ ! -f $DIR/bams/ont/$bamname ]; then
		ln -s $bamfile $DIR/bams/ont/$bamname
	fi
done