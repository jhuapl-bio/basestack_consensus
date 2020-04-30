#!/bin/bash
. "/home/idies/workspace/Storage/ernluaw1/persistent/Miniconda3/etc/profile.d/conda.sh"
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
	echo -e "   -i      path/to/sequencing run folder"
    echo -e "   -t      number of threads"
	#echo -e "   -b      repository path (working on calculating this directly)"
	#echo -e "   -1      barcode demux folder (default: ${CYAN}<run-folder>/artic-pipeline/3-normalization${NC})"
	echo -e ""
}

#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:t:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) sequencing_run=$OPTARG ;;
        t) threads=$OPTARG ;;
		#b) bin_path=$OPTARG ;;
		#1) demux_path=$OPTARG ;;
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
	echo -e "[$(date +"%F %T")]$input\r" | sed -r 's/\x1b\[[0-9;]*m?//g'
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

if [ ! -d ${sequencing_run}/artic-pipeline/2-length-filter ];then
    >&2 echo "Error Require module 2-length-filter output"
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

# input files, these files should be in the sequencing run directory
manifest=${sequencing_run}/manifest.txt
run_configuration="${sequencing_run}/run_config.txt"
gather_dir=${sequencing_run}/artic-pipeline/2-length-filter

# location for primer schemes
scheme_dir=${software_path}/artic-ncov2019/primer_schemes

# primer protocol
protocol=$(awk '/primers/{ print $2 }' "${run_configuration}")

# Output directories
normalize_dir=${sequencing_run}/artic-pipeline/3-normalization

# Optional program parameters
norm_parameters="coverage_threshold=150 --qual_sort" 

#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "  sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "recording software version numbers"
echo_log "Coverage Normalization from: https://github.com/mkirsche/CoverageNormalization"
echo_log "Coverage Normalization parameters: $norm_parameters"
echo_log "run configuration file: ${sequencing_run}/run_config.txt"
echo_log "run manifest file: ${manifest}"
echo_log "inputs: Length directory: ${gather_dir}"
echo_log "output normalization directory: ${normalize_dir}"
echo_log "------ processing pipeline output ------"

#---------------------------------------------------------------------------------------------------
# module 3
#---------------------------------------------------------------------------------------------------

echo_log "Starting normalize module 3 $sequencing_run"

while read barcode name; do

# create output directories, need separate directories for each sample 
mkdir -p $normalize_dir/${name}_${barcode}
cd $normalize_dir/${name}_${barcode}

# create alignment file
align_out=${name}_${barcode}.sam

minimap2 -a \
-x map-ont \
-t $threads \
$scheme_dir/$protocol/nCoV-2019.reference.fasta \
$gather_dir/${name}_${barcode}.fastq > $normalize_dir/${name}_${barcode}/$align_out

# normalization, txt file output went to working directory
out_sam=$normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.sam

$JAVA_PATH/java -cp $NormalizeCoveragePath/src NormalizeCoverage input=$normalize_dir/${name}_${barcode}/$align_out $norm_parameters

# fastq conversion
$samtools_path/samtools fastq $out_sam > $normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.fq

done < "${manifest}"
    
#---------------------------------------------------------------------------------------------------

echo_log "run complete"
#chgrp -R 5102 $demux_dir

