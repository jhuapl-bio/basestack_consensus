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
   
OUTDIR=$DIR

#SCRIPTS PATH
SCRIPT_DIR="/home/idies/workspace/covid19/code/ncov/pipeline_scripts" 

CONFIG="/home/idies/workspace/covid19/ncov_reference/snpEff.config"
DBNAME="ncov"

#mkdir -p $OUTDIR

#for vcf in $TEST_DIR/*.consensus.combined.vcf; do
for vcf in $DIR/*.allsnps.combined.vcf; do
    ${SCRIPT_DIR}/annotate_variants.sh ${vcf} ${CONFIG} ${DBNAME} ${OUTDIR}
    echo "SnpEff completed on run ${DIR}"
done

echo "Making final reports on run ${DIR}"
cat ${OUTDIR}/MDHP*_ann_report.txt  | awk '$4 != "N" { print $0}'  | awk '!seen[$0]++' > ${OUTDIR}/final_snpEff_report.txt
cat ${OUTDIR}/MDHP*_ann_report.txt  | awk '!seen[$0]++' | awk 'NR == 1  || $4 == "N" { print $0}'  > ${OUTDIR}/snpEff_report_with_Ns.txt
echo "DONE"
