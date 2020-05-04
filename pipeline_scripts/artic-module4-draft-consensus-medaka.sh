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
        echo -e "   -i      /full/path/to/normalizd_sample.fq"
        echo -e "   -t      number of threads (default: 6)"
        echo -e ""
}

#---------------------------------------------------------------------------------------------------
#default threads
threads=6
#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:t:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) normalized_fastq=$OPTARG ;;
                t) threads=$OPTARG ;;
                ?) usage; exit ;;
       esac
done

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
# Sequencing run directory
#===================================================================================================

# sequencing run directory
sequencing_run=$(dirname $(dirname $(dirname $(dirname "$normalized_fastq"))))

#===================================================================================================
# QUALITY CHECKING
#===================================================================================================

if [ ! -f ${normalized_fastq} ];then
    >&2 echo "Fastq file ${normalized_fastq} does not exist"
    exit 1
fi

if [ ! -d ${sequencing_run} ];then
    >&2 echo "Error Sequencing run ${sequencing_run} does not exist"
    exit 1
fi

if [ ! -s ${sequencing_run}/run_config.txt ];then
    # test why this is getting triggered
    >&2 echo "Error Require a run_config.txt file in the sequencing run directory"
    >&2 echo "${sequencing_run}/run_config.txt does not exist"
    exit 1
fi

if [ ! -s ${sequencing_run}/manifest.txt ];then
    >&2 echo "Error Require a manifest.txt file in the sequencing run directory"
    >&2 echo "${sequencing_run}/manifest.txt does not exist"
    exit 1
fi

if [ ! -d ${sequencing_run}/artic-pipeline/3-normalization ];then
    >&2 echo "Error Require module 3 normalization output"
    >&2 echo "${sequencing_run}/artic-pipeline/3-normalization does not exist"
    exit 1
fi

#===================================================================================================
# Default values
#===================================================================================================

# location of programs used by pipeline
software_path=/home/idies/workspace/covid19/code

# log file
logfile=${sequencing_run}/pipeline.log

# input files, these files should be in the sequencing run directory
manifest=${sequencing_run}/manifest.txt
run_configuration="${sequencing_run}/run_config.txt"

# location for primer schemes
scheme_dir=${software_path}/artic-ncov2019/primer_schemes

# primer protocol
protocol=$(awk '/primers/{ print $2 }' "${run_configuration}")

# Output directories
consensus_dir=${sequencing_run}/artic-pipeline/4-draft-consensus

# Optional program parameters
out_prefix="$consensus_dir/$(basename ${normalized_fastq%.covfiltered.fq}.medaka)"


#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "------ Medaka Paramters:"
echo_log "sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "recording software version numbers..."
echo_log "Software version: $(medaka --version)"
echo_log "run configuration file: ${sequencing_run}/run_config.txt"
echo_log "run manifest file: ${manifest}"
echo_log "fasta file: ${normalized_fastq}"
echo_log "output medaka directory: ${consensus_dir}"
echo_log "------ processing pipeline output ------"

#---------------------------------------------------------------------------------------------------
# module 4
#---------------------------------------------------------------------------------------------------

echo_log "Starting Module 4 Medaka on $normalized_fastq"
mkdir -p $consensus_dir

artic minion \
        --medaka \
        --normalise 1000000 \
        --threads $threads \
        --scheme-directory "$scheme_dir" \
        --read-file $normalized_fastq \
        "$protocol" "$out_prefix"



#---------------------------------------------------------------------------------------------------

echo_log "Module 4 Medaka: processing complete"
#chgrp -R 5102 $demux_dir

