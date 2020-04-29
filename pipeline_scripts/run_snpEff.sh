#!/bin/bash

if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

#Run directory
RUN=$1

if [ $RUN == "test_data" ]; then
   DIR="/home/idies/workspace/covid19/sequencing_runs/test_data/artic-pipeline/5-post-filter_nanopolish/"
   echo "Running SnpEff on test data "
else
   DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter_nanopolish"
   
OUTDIR=$DIR

#SCRIPTS PATH
SCRIPT_DIR="/home/idies/workspace/covid19/code/ncov/pipeline_scripts" 

CONFIG="/home/idies/workspace/covid19/ncov_reference/snpEff.config"
DBNAME="ncov"

#mkdir -p $OUTDIR

#for vcf in $TEST_DIR/*.consensus.combined.vcf; do
for vcf in $DIR/*.consensus.combined.vcf; do
    ${SCRIPT_DIR}/annotate_variants.sh ${vcf} ${CONFIG} ${DBNAME} ${OUTDIR}
    echo "SnpEff completed on run ${DIR}"
done
