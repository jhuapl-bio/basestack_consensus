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

SCRIPT_DIR="/home/idies/workspace/covid19/code/ncov/pipeline_scripts" 
REF_GB="/home/idies/workspace/covid19/ncov_reference/reference_seq.gb"
NEXTSTRAIN_CLADES="/home/idies/workspace/covid19/ncov_reference/clades.tsv"

source /home/idies/workspace/covid19/bashrc
conda activate nextstrain

echo "Assigning nextstrain clades for consensus sequences in ${DIR}"

CONS_FASTA=$OUTDIR/postfilt_consensus_all.fasta

if [ ! -f "$CONS_FASTA" ]; then
    echo " File : $CONS_FASTA does not exist... Making "
    cat ${DIR}/MDHP*.complete.fasta > $CONS_FASTA
fi

#usage: assign_clades.py [-h] --sequences SEQUENCES --clade CLADE --gbk GBK [--output OUTPUT] [--keep-temporary-files] [--chunk-size CHUNK_SIZE]
#                        [--nthreads NTHREADS]

${SCRIPT_DIR}/assign_clades.py --sequences ${CONS_FASTA} --output ${OUTDIR}/nextstrain_clades.tsv --gbk ${REF_GB} --clade ${NEXTSTRAIN_CLADES}
echo "DONE"
