#!/usr/bin/env bash


###
# Script to build new databases for SnpEff
###

FASTA=$1
GFF3=$2
DBNAME=$3
SNPEFF_DIR=$4

# CHECK 1: INPUT FILES

usage() {
    cat <<EOM
    Usage: $(basename $0) <fasta> <gff> <dbname> <path to snpeff config dir>

EOM
    exit 0
}

PATHENV="/home/idies/workspace/Storage/sramakr4/persistent/miniconda3/etc/profile.d/conda.sh"
JARPATH="/home/idies/workspace/Storage/sramakr4/persistent/miniconda3/envs/nextstrain/share/snpeff-4.3.1t-3/"

source ${PATHENV}
conda activate nextstrain

[ -z $1 ] && { usage; }

if [[ ! -r "$FASTA" ]]
then
        echo "$0: input $FASTA not found"
        exit 1
fi

if [[ ! -r "$GFF3" ]]
then
        echo "$0: input $GFF3 not found"
        exit 1
fi

# CHECK 2: LOAD ENVIRONMENT
. "/home/idies/workspace/Storage/sramakr4/persistent/miniconda3/etc/profile.d/conda.sh"
conda activate nextstrain

if ! [ -x "$(command -v snpEff)" ]; then
  echo 'Error: snpEff is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v java)" ]; then
  echo 'Error: java is not installed.' >&2
  exit 1
fi

# MAKE DIRS
#Go into the snpEff directory and create a directory for your files
cd $SNPEFF_DIR
mkdir -p data/$DBNAME

#Copy the files into snpEff's directory structure
cp $GFF3 data/$DBNAME/genes.gff
cp $FASTA data/$DBNAME/sequences.fa

#Edit snpEff.config and insert your specific database information:
echo "$DBNAME.genome : $DBNAME" > snpEff.config

#Build the database
java -jar ${JARPATH}/snpEff.jar build -gff3 -v $DBNAME
