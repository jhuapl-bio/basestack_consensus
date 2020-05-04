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
        echo -e ""
}

#---------------------------------------------------------------------------------------------------

# parse input arguments
# parse input arguments
while getopts "hi:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) sequencing_run=$OPTARG ;;
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
# QUALITY CHECKING
#===================================================================================================

if [ ! -d ${sequencing_run} ];then
    >&2 echo "Error Sequencing run ${sequencing_run} does not exist"
    exit 1
fi

if [ ! -s ${sequencing_run}/run_config.txt ];then
    >&2 echo "Error Require a run_config.txt file in the sequencing run directory"
    exit 1
fi

if [ ! -s ${sequencing_run}/manifest.txt ];then
    >&2 echo "Error Require a manifest.txt file in the sequencing run directory"
    exit 1
fi

if [ ! -d ${sequencing_run}/artic-pipeline/1-barcode-demux ];then
    >&2 echo "Error Require module 1-barcode-demux output"
    exit 1
fi

#===================================================================================================
# Default values
#===================================================================================================

# location of programs used by pipeline

# log file
logfile=${sequencing_run}/pipeline.log

# input files, these files should be in the sequencing run directory
manifest=${sequencing_run}/manifest.txt
demux_dir=${sequencing_run}/artic-pipeline/1-barcode-demux

# Output directories
gather_dir=${sequencing_run}/artic-pipeline/2-length-filter

# Optional program parameters


#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "  sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "recording software version numbers"
echo_log "Artic guppyplex from: $(artic --version)"
echo_log "run configuration file: ${sequencing_run}/run_config.txt"
echo_log "run manifest file: ${manifest}"
echo_log "inputs: fastq_directory: ${demux_dir}"
echo_log "output gather directory: ${gather_dir}"
echo_log "------ processing pipeline output ------"

#---------------------------------------------------------------------------------------------------
# module 2
#---------------------------------------------------------------------------------------------------

echo_log "Starting artic guppyplex module 2 on $sequencing_run"
mkdir -p $gather_dir

while read barcode name; do
echo_log "${name}_${barcode}"
artic guppyplex \
	--skip-quality-check \
	--min-length 400 \
	--max-length 700 \
	--directory "$demux_dir"/"${barcode}" \
    --prefix "$gather_dir"/"${name}"
      
done < "$manifest"
    
#---------------------------------------------------------------------------------------------------

echo_log "run complete"
#chgrp -R 5102 $demux_dir

#===================================================================================================
# QUALITY CHECKING AND MODULE 3 JOB SUBMISSION
#===================================================================================================

if [ ! -d $gather_dir ];then
    >&2 echo_log "Error $gather_dir not created"
    exit 1
fi

#if find "$gather_dir" -maxdepth 0 -empty | read;
#    echo_log "$gather_dir empty."
#else
#    echo_log "Begin submitting module 3"
#    python <submit_module2>.py
#fi

