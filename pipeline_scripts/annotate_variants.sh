#!/usr/bin/env bash

###

# Script to snpeff annotations for all vcf

###

VCF=$1
CONFIG=$2
DBNAME=$3
OUT_DIR=$4

# CHECK 1: INPUT FILES

usage() {
    cat <<EOM
    Usage: $(basename $0) <vcf> <config> <dbname> <path to snpeff output dir>

EOM
    exit 0
}

[ -z $1 ] && { usage; }

PATHENV="/home/idies/workspace/Storage/sramakr4/persistent/miniconda3/etc/profile.d/conda.sh"
JARPATH="/home/idies/workspace/Storage/sramakr4/persistent/miniconda3/envs/nextstrain/share/snpeff-4.3.1t-3/"

VCF_BASE=$( basename ${VCF} ".vcf")
CONFIG_DIR=$( dirname ${CONFIG} )
CONFIG_DATA=${CONFIG_DIR}/data/
AA_DATA=${CONFIG_DIR}/amino_acid_codes.txt

source ${PATHENV}
conda activate nextstrain


if [[ ! -r "$CONFIG" ]]
then
        echo "$0: input $CONFIG not found"
        exit 1
fi

if [[ ! -r "$VCF" ]]
then
        echo "$0: input $VCF not found"
        exit 1
fi

# CHECK 2: LOAD ENVIRONMENT

if ! [ -x "$(command -v snpEff)" ]; then
  echo 'Error: snpEff is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v java)" ]; then
  echo 'Error: java is not installed.' >&2
  exit 1
fi

# MAKE DIRS
java -Xmx4g -jar $JARPATH/snpEff.jar eff -c ${CONFIG} -dataDir ${CONFIG_DATA} ncov ${VCF} > ${OUT_DIR}/${VCF_BASE}_ann.vcf

if [[ ! -r "${OUT_DIR}/${VCF_BASE}_ann.vcf" ]]
then
        echo "$0: ERROR while running SnpEff"
        exit 1
fi

if grep -q "ERROR_" ${OUT_DIR}/${VCF_BASE}_ann.vcf; 
then
   error=$( grep -n "ERROR_" ${OUT_DIR}/${VCF_BASE}_ann.vcf )
   echo "$0: ERROR while running SnpEff ${error} "
   exit 1
fi
  
# Make report
# Report of all 3 letter amino acide codes
awk '/^#/ { if ( $1 == "#CHROM" ) { new_fields="GENE\tANN\tAA_MUT" ; OFS="\t"; print $0"\t"new_fields ; next } ; print ; next } {  split($8,a,"|") ; split(a[11],m,".") ; o_cols="."; ann=a[4]"\t"a[2]"\t"m[2] ; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"o_cols"\t"ann }' ${OUT_DIR}/${VCF_BASE}_ann.vcf > ${OUT_DIR}/${VCF_BASE}_ann_3letter_code.vcf


# Report of all clean vcf
awk '/^#/ { if ( $1 == "#CHROM" ) { OFS="\t" ; print } else { print }}' ${OUT_DIR}/${VCF_BASE}_ann_3letter_code.vcf > ${OUT_DIR}/${VCF_BASE}_ann_clean.vcf 


### Modify this for all mutations
awk 'NR==FNR { a[$2] = $3 ; next } { FS="\t" ; if ( !/^#/ ) {  OFS="\t" ; last=$NF; gsub(/[0-9]*/,"",last) ; for(j=1;j<length(last);j+=3) { k=substr(last,j,3) ; gsub(k,a[k],$NF) } ; print }}' ${AA_DATA} ${OUT_DIR}/${VCF_BASE}_ann_3letter_code.vcf >> ${OUT_DIR}/${VCF_BASE}_ann_clean.vcf

# Tab demilited report
echo -e "#CHROM\tPOS\tREF\tALT\tGENE\tANN\tAA_MUT"  > ${OUT_DIR}/${VCF_BASE}_ann_report.txt 
awk '!/^#/{ print $1"\t"$2"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11}' ${OUT_DIR}/${VCF_BASE}_ann_clean.vcf >> ${OUT_DIR}/${VCF_BASE}_ann_report.txt 


exit 0
## DONE
