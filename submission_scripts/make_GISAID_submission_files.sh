#!/bin/bash

if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

GOOGLE_API_KEY="/home/idies/workspace/covid19/config/keys/covid19_google_api_key.json"
LINK="https://docs.google.com/spreadsheets/d/11pXGdGqhe1R408HnJCzfSEN0oxSOGuXiWIX_rHNCmYQ/edit?usp=sharing"
SEQ_PATH="/home/idies/workspace/covid19/sequencing_runs/"
FASTA_PATH="artic-pipeline/5-post-filter"
SUB_MANIFEST="/home/idies/workspace/covid19/jhu_sequences/master/submission_manifest.tsv"
SUBMISSIONS_PATH="/home/idies/workspace/covid19/jhu_sequences/submissions"
SUBMITTED_SEQ=${SUBMISSIONS_PATH}"/all_submitted_jhu_sequences.fasta"
MASTER_META="/home/idies/workspace/covid19/jhu_sequences/master/jhu_sequences_metadata_master.tsv"

SCRIPTS_PATH="/home/idies/workspace/covid19/code/ncov/submission_scripts"
source /home/idies/workspace/covid19/bashrc
conda activate nextstrain

#usage: prepare_for_GISAID_submission.py [-h] -in EXCEL_LINK --api_key API_KEY --submission_manifest SUB_MANIFEST --seq-path SEQ_PATH --fasta-path FASTA_PATH
#                                        --submission-path SUBMISSIONS_PATH --submitted-fasta SUBMITTED_SEQ --master-meta MASTER_META

${SCRIPTS_PATH}/prepare_for_GISAID_submission.py -in ${LINK} --api_key ${GOOGLE_API_KEY} --submission_manifest ${SUB_MANIFEST} --seq-path ${SEQ_PATH} --fasta-path ${FASTA_PATH} --submission-path ${SUBMISSIONS_PATH} --submitted-fasta ${SUBMITTED_SEQ} --master-meta ${MASTER_META}
echo "DONE"
