#!/bin/bash

if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

#Run directory
RUN=$1

if [ $RUN == "test_data" ]; then
   DIR="/home/idies/workspace/covid19/sequencing_runs/test_data/artic-pipeline/5-post-filter/"
   echo "Running SnpEff on test data "
else
   DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter"
fi   
OUTDIR=$DIR

#SCRIPTS PATH
SCRIPT_DIR="/home/idies/workspace/covid19/code/ncov/pipeline_scripts" 
DATA="/home/idies/workspace/covid19/ncov_reference/lineages/lineages/data/"
TMPDIR=$OUTDIR
source /home/idies/workspace/covid19/bashrc
conda activate pangolin
echo "Making Pangolin lineages for consensus sequences in ${DIR}"
cat ${DIR}/MDHP*.complete.fasta > $OUTDIR/postfilt_consensus_all.fasta
pangolin $OUTDIR/postfilt_consensus_all.fasta -f -d ${DATA} -o ${OUTDIR} --tempdir $TMPDIR
echo "DONE"
