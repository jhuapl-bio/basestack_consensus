#!/bin/bash
source /home/idies/workspace/covid19/bashrc

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
	echo -e "   -i      path/to/sequencing_run_folder"
	echo -e ""
}

#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) sequencing_run=$OPTARG ;;
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

# check for the existence of the sequencing run directory
if [ ! -d ${sequencing_run} ];then
    >&2 echo "Error Sequencing run ${sequencing_run} does not exist"
    exit 1
fi

# check for existence of run_config.txt and for barcoding 
if [ ! -s ${sequencing_run}/run_config.txt ];then
    >&2 echo "Error Require a run_config.txt file in the sequencing run directory"
    exit 1
else
    if ! grep -q "barcoding" ${sequencing_run}/run_config.txt;then 
        echo "require barcoding file"
        >&2 echo "Error barcode file not found"
        exit 1
    fi
fi

# check for existence of manifest.txt and that it has two columns
if [ ! -s ${sequencing_run}/manifest.txt ];then 
    >&2 echo "Error Require a manifest.txt file in the sequencing run directory"
    exit 1
else
    columns=$( awk -F' ' '{print NF}' ${sequencing_run}/manifest3.txt )
    if [ $columns -ne 2 ];then 
        >&2 echo "Error manifest.txt file does not have two columns"
        exit 1
    fi
fi

# check for existence of fastq_pass directory
if [ ! -d ${sequencing_run}/fastq_pass ];then
    >&2 echo "Error Require fastq_pass directory in the sequencing run directory"
    exit 1
fi

#===================================================================================================
# Default values
#===================================================================================================

# location of programs used by pipeline - double check if source bashrc don't have to hardcode paths
software_path=/home/idies/workspace/covid19/code
guppy_barcoder_path=${software_path}/ont-guppy-cpu/bin

# log file
logfile=${sequencing_run}/artic-pipeline/pipeline.log

# input files, these files should be in the sequencing run directory
run_configuration="${sequencing_run}/run_config.txt"
barcode_file=$(awk '/barcoding/{ print $2 }' "${run_configuration}")
fastq_dir=${sequencing_run}/fastq_pass
manifest=${sequencing_run}/manifest.txt

# Output directories
demux_dir=${sequencing_run}/artic-pipeline/1-barcode-demux

#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "  sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "recording software version numbers"
echo_log $(${guppy_barcoder_path}/guppy_barcoder --version)
echo_log "run configuration file: ${sequencing_run}/run_config.txt"
echo_log "run manifest file: ${sequencing_run}/manifest.txt"
echo_log "inputs: fastq_directory: ${fastq_dir}, arrangements files: ${barcode_file}"
echo_log "output demultiplex directory: ${demux_dir}"
echo_log "------ processing pipeline output ------"

#---------------------------------------------------------------------------------------------------
# module 1 
#---------------------------------------------------------------------------------------------------

echo_log "Starting guppy demux module 1 $sequencing_run"

$guppy_barcoder_path/guppy_barcoder \
	--require_barcodes_both_ends \
	-i "$fastq_dir" \
	-s "$demux_dir" \
	--arrangements_files $barcode_file

while read barcode name; do
mv "$demux_dir"/"${barcode/NB/barcode}" "$demux_dir"/"${barcode}"
done < "$manifest"
    
#---------------------------------------------------------------------------------------------------

#chgrp -R 5102 $demux_dir

#===================================================================================================
# QUALITY CHECKING AND MODULE 2 JOB SUBMISSION
#===================================================================================================

if [ ! -d $demux_dir ];then
    >&2 echo_log "Error $demux_dir not created"
    exit 1
fi

while read name barcode; do
    if [ -d "$demux_dir"/"$name" ]; then
        complete=TRUE
    else
        complete=FALSE
        echo_log Error "$demux_dir"/"$name" does not exist
        exit 1
    fi
done < $manifest

if complete=TRUE; then
   touch $demux_dir/1-barcode-demux.complete
   echo_log "run complete"
fi
    
if [ -s $demux_dir/1-barcode-demux.complete ]; then
    while read name barcode; do
        echo_log executing <submit_module2>.py artic-module2-length-filter.sh "$demux_dir"/"$name"       
        #python <submit_module2.py>
    done < $manifest
fi




