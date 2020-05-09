#!/bin/bash
source /home/idies/workspace/covid19/bashrc
conda activate artic-ncov2019-medaka

#---------------------------------------------------------------------------------------------------

# set default values here

# define colors for error messages
red='\033[0;31m'
RED='\033[1;31m'
green='\033[0;32m'
GREEN='\033[1;32m'
yellow='\033[0;33m'
YELLOW='\033[1;33m'
blue='\033[0;34m'
BLUE='\033[1;34m'
purple='\033[0;35m'
PURPLE='\033[1;35m'
cyan='\033[0;36m'
CYAN='\033[1;36m'
NC='\033[0m'

# usage function
usage() {
	echo -e "usage: ${YELLOW}$0${NC} [options]"
	echo -e ""
	echo -e "OPTIONS:"
	echo -e "   -h      show this message"
	echo -e "   -i      path to 2-length-filter fastq file (i.e. full/path/sequencing_run/artic-pipeline/2-length-filter/<sample>.fq)"
    echo -e "   -t      number of threads"
	echo -e ""
}

#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:t:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) fastq=$OPTARG ;;
        t) threads=$OPTARG ;;
		?) usage; exit ;;
	esac
done

if ((OPTIND == 1))
then
    usage
    exit
fi

#===================================================================================================
# DEFINE FUNCTIONS
#===================================================================================================

echo_log() {
	input="$*"
	# if input is non-empty string, prepend initial space
	if [[ -n "$input" ]]; then
		input=" $input"
	fi
	# print to STDOUT
	#echo -e "[$(date +"%F %T")]$input"
	# print to log file (after removing color strings)
	echo -e "[$(date +"%F %T")]$input\r" | sed -r 's/\x1b\[[0-9;]*m?//g' >> "$logfile"
}

#===================================================================================================
# QUALITY CHECKING
#===================================================================================================

sequencing_run=$(dirname $(dirname $(dirname "$fastq")))

if [ ! -s ${fastq} ];then
    >&2 echo_log "Error Sequencing run ${sequencing_run} does not exist"
    exit 1
fi

if [ ! -s ${sequencing_run}/run_config.txt ];then
    >&2 echo_log "Error Require a run_config.txt file in the sequencing run directory"
    exit 1
fi

if [ ! -s ${sequencing_run}/manifest.txt ];then
    >&2 echo_log "Error Require a manifest.txt file in the sequencing run directory"
    exit 1
fi

if [ ! -s $fastq ];then
    >&2 echo_log "Error: Module 2 output '${fastq}' not found"
    exit 1
fi

#===================================================================================================
# Default values
#===================================================================================================

# location of programs used by pipeline
software_path=/home/idies/workspace/covid19/code
JAVA_PATH="${software_path}/jdk-14.0.1/bin"
samtools_path="${software_path}/samtools-1.10/bin"
NormalizeCoveragePath="${software_path}/CoverageNormalization"

# input files and directories
base=$(basename "${fastq%.fastq}")
manifest=${sequencing_run}/manifest.txt
run_configuration="${sequencing_run}/run_config.txt"
gather_dir=${sequencing_run}/artic-pipeline/2-length-filter

# reference sequence
reference="$scheme_dir/$protocol/nCoV-2019.reference.fasta"

# location for primer schemes
scheme_dir="$software_path/artic-ncov2019/primer_schemes"

# primer protocol
protocol=$(awk '/primers/{ print $2 }' "${sequencing_run}/run_config.txt")

# reference fasta
reference="$scheme_dir/$protocol/nCoV-2019.reference.fasta"

# Output directories
normalize_dir="${sequencing_run}/artic-pipeline/3-normalization/$base"


# Optional program parameters - check with Tom , even_strand or median_strand
norm_parameters="coverage_threshold=150 --qual_sort --even_strand" 

# create output directories, need separate directories for each sample 
mkdir -p $normalize_dir

# log file
logfile="${sequencing_run}/artic-pipeline/3-normalization/$(date +"%F-%H%M%s"-module3.log"


#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "SAMPLE ${base}:  sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "SAMPLE ${base}: recording software version numbers"
echo_log "SAMPLE ${base}: Coverage Normalization from: https://github.com/mkirsche/CoverageNormalization"
echo_log "SAMPLE ${base}: Coverage Normalization parameters: $norm_parameters"
echo_log "SAMPLE ${base}: run configuration file: ${sequencing_run}/run_config.txt"
echo_log "SAMPLE ${base}: run manifest file: ${manifest}"
echo_log "SAMPLE ${base}: inputs: Length directory: ${gather_dir}"
echo_log "SAMPLE ${base}: output normalization directory: ${normalize_dir}"
echo_log "SAMPLE ${base}: ------ processing pipeline output ------"

#---------------------------------------------------------------------------------------------------
# module 3
#---------------------------------------------------------------------------------------------------

echo_log "SAMPLE ${base}: Starting normalize module 3"

cd $normalize_dir

# create alignment file
align_out="$normalize_dir/$base.sam"
echo_log "SAMPLE ${base}: output file = $align_out" 

minimap2 -a \
	-x map-ont \
	-t 32 \
	"$reference" \
	"$fastq" > "$align_out" 2>> "$logfile"

samtools sort $align_out > ${align_out%.sam}.bam 2>> "$logfile"
samtools depth -a -d 0 "${align_out%.sam}.bam" > "${align_out%.sam}.depth" 2>> "$logfile"

# normalization, txt file output went to working directory
out_sam="${align_out%.sam}.covfiltered.sam"

$JAVA_PATH/java \
	-cp $NormalizeCoveragePath/src \
	NormalizeCoverage \
	input="$align_out" \
	$norm_parameters 2>> "$logfile"

# fastq conversion
samtools fastq "${out_sam}" > "${out_sam%.sam}.fq" 2>> "$logfile"

samtools sort "${out_sam}" > "${out_sam%.sam}.bam" 2>> "$logfile"
samtools depth -a -d 0 "${out_sam%.sam}.bam" > "${out_sam%.sam}.depth" 2>> "$logfile"

#---------------------------------------------------------------------------------------------------

#chgrp -R 5102 $demux_dir

#===================================================================================================
# QUALITY CHECKING AND MODULE 4 JOB SUBMISSION
#===================================================================================================

if [ ! -f "${out_sam%.sam}.fq" ];then
    >&2 echo_log "SAMPLE ${base}: Error: Module 3 output "${out_sam%.sam}.fq" not found"
    exit 1
else
	echo_log "SAMPLE ${base}: Module 3 complete for sample '${base}'"
	submit_sciserver_ont_job.py -m 3 -i "${out_sam%.sam}.fq -t 5 2>> "$logfile" 
fi



