#!/bin/bash

#---------------------------------------------------------------------------------------------------

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
	echo -e "   -i      sequencing run folder"
	echo -e "   -b      repository path (working on calculating this directly)"
	echo -e "   -p      protocol (default: )"
	echo -e "   -r      reference FASTA (default: )"
#	echo -e "   -o      output folder (default: ${CYAN}<run-folder>/artic-pipeline/run_stats${NC})"
	echo -e "   -1      barcode demux folder (default: ${CYAN}<run-folder>/artic-pipeline/1-barcode-deumux${NC})"
	echo -e "   -2      length filter folder (default: ${CYAN}<run-folder>/artic-pipeline/2-length-filter{NC})"
	echo -e "   -3      normalization folder (default: ${CYAN}<run-folder>/artic-pipeline/3-normalize${NC})"
	echo -e "   -4      draft consensus folder (default: ${CYAN}<run-folder>/artic-pipeline/4-draft-consensus${NC})"
	echo -e "   -5      nextstrain folder (default: ${CYAN}not used yet${NC})"
	echo -e "   -6      post-filter folder (default: ${CYAN}not used yet${NC})"
	echo -e ""
}
#---------------------------------------------------------------------------------------------------
# set default values here - Tom's defaults
logfile="/dev/null"
tempdir="/tmp"

# report current hash of miseq-analysis git repo
#script_path=$(dirname $(readlink -f $0))
#GIT_DIR="$script_path/.git"
#export GIT_DIR
#hash=$(git rev-parse --short HEAD)

#echo $GIT_DIR

# the sequencing run ID - update*
runID=$1

logfile=/home/idies/workspace/Temporary/ernluaw1/scratch/log_${runID}-submit.txt

# input sequecing directory
sequencing_run="/home/idies/workspace/covid19/sequencing_runs/${runID}"

# location of programs used by pipeline
software_path=/home/idies/workspace/covid19/code

# location for primer schemes
scheme_dir=${software_path}/artic-ncov2019/primer_schemes

# primer protocol
protocol=$(awk '/primers/{ print $2 }' "${sequencing_run}/run_config.txt")

# Have these files be in the sequencing run directory - update*
barcode_file=$(awk '/barcoding/{ print $2 }' "${sequencing_run}/run_config.txt")
manifest=${sequencing_run}/manifest.txt
fastq_dir=${sequencing_run}/fastq_pass

#---------------------------------------------------------------------------------------------------
# parse input arguments
while getopts "hi:b:p:r:o:1:2:3:4:5:6:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) run_path=$OPTARG ;;
		b) bin_path=$OPTARG ;;
		p) run_path=$OPTARG ;;
		r) reference=$OPTARG ;;
		o) stats_path=$OPTARG ;;
		1) demux_path=$OPTARG ;;
		2) lengthfilter_path=$OPTARG ;;
		3) normalize_path=$OPTARG ;;
		4) draftconsensus_path=$OPTARG ;;
		5) nextstrain_path=$OPTARG ;;
		6) postfilter_path=$OPTARG ;;
		?) usage; exit ;;
	esac
done

#===================================================================================================
# QUALITY CHECKING
#===================================================================================================

# Output directories
pipeline_label=hac-medaka-norm200
demux_dir=${sequencing_run}/artic/1-barcode-demux
gather_dir=${sequencing_run}/artic/2-length-filter
normalize_dir=${sequencing_run}/artic/3-normalization
consensus_dir=${sequencing_run}/artic/4-draft-consensus

#===================================================================================================
# DEFINE FUNCTIONS
#===================================================================================================

#---------------------------------------------------------------------------------------------------
# log function that will output to STDOUT and a log file
echo_log() {

	input="$*"

	# if input is non-empty string, prepend initial space
	if [[ -n "$input" ]]; then
		input=" $input"
	fi

	# print to STDOUT
	echo -e "[$(date +"%F %T")]$prefix$input"

	# print to log file (after removing color strings)
	echo -e "[$(date +"%F %T")]$prefix$input" | gawk '{ printf("%s\n", gensub(/\x1b\[[0-9;]*m?/, "", "g", $0)); }' >> "$logfile"
}

#---------------------------------------------------------------------------------------------------

#===================================================================================================
# MAIN BODY
#===================================================================================================

echo -e "$(date +"%F %T") sequencing run start $runID"
#echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

# create directory to hold temporary files
#runtime=$(date +"%Y%m%d%H%M%S%N")
#workdir="$tempdir/build_krakendb-$runtime"
#mkdir -m 775 -p "$workdir"

#echo_log "recording software version numbers"
#echo_log "current git hash: $hash"
#echo_log "  guppy barcoder: "
#echo_log "input arguments"
#echo_log "  sequencing run folder: ${CYAN}$run_path${NC}"
#echo_log "  working directory: ${CYAN}$workdir${NC}"
#echo_log "  threads: ${CYAN}1${NC}"
#echo_log "output arguments"
#echo_log "  log file: ${CYAN}$logfile${NC}"
#echo_log "  summary file: ${CYAN}$summary${NC}"
#echo_log "  depth file: ${CYAN}$summary${NC}"
#echo_log "  mutation file: ${CYAN}$summary${NC}"
#echo_log "------ processing pipeline output ------"

#exit

#---------------------------------------------------------------------------------------------------
# module 1
#---------------------------------------------------------------------------------------------------

# need to fix hardcoded path to software
guppy_barcoder_path=${software_path}/ont-guppy-cpu/bin

echo_log "Starting guppy demux"

$guppy_barcoder_path/guppy_barcoder \
	--require_barcodes_both_ends \
	-i "$fastq_dir" \
	-s "$demux_dir" \
	--arrangements_files $barcode_file
	

#---------------------------------------------------------------------------------------------------
# module 2
#---------------------------------------------------------------------------------------------------

echo_log "Starting artic guppyplex"
mkdir -p $gather_dir

while read barcode name; do
echo_log "${name}_${barcode}"
artic guppyplex \
	--skip-quality-check \
	--min-length 400 \
	--max-length 700 \
	--directory "$demux_dir"/"$barcode" \
    --prefix "$gather_dir"/"$name"
done < "$manifest"

#---------------------------------------------------------------------------------------------------
# module 3
#---------------------------------------------------------------------------------------------------

echo_log "Starting normalize"

echo_log "Starting normalize"

# software - need to move java software and samtools over to code directory
JAVA_PATH="${software_path}/jdk-14.0.1/bin"
samtools_path="${software_path}/samtools-1.10/bin"
NormalizeCoveragePath="${software_path}/CoverageNormalization"

while read barcode name; do

# create output directories, need separate directories for each sample 
mkdir -p $normalize_dir/${name}_${barcode}
cd $normalize_dir/${name}_${barcode}

# create alignment file
align_out=${name}_${barcode}.sam

minimap2 -a \
-x map-ont \
-t 32 \
$scheme_dir/$protocol/nCoV-2019.reference.fasta \
$gather_dir/${name}_${barcode}.fastq > $normalize_dir/${name}_${barcode}/$align_out

# normalization,.txt file output went to working directory
out_sam=$normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.sam

$JAVA_PATH/java -cp $NormalizeCoveragePath/src NormalizeCoverage input=$normalize_dir/${name}_${barcode}/$align_out coverage_threshold=150 --qual_sort

# fastq conversion
$samtools_path/samtools fastq $out_sam > $normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.fq

done < "${manifest}"

#---------------------------------------------------------------------------------------------------
# module 4
#---------------------------------------------------------------------------------------------------

echo_log "Starting medaka"
mkdir -p $consensus_dir

while read barcode name; do
echo_log "${name}_${barcode}"
artic minion \
	--medaka \
	--normalise 1000000 \
	--threads 32 \
	--scheme-directory "$scheme_dir" \
	--read-file $normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.fq \
	"$protocol" "$consensus_dir"/${name}_${barcode}-${pipeline_label}
done < "$manifest"

#---------------------------------------------------------------------------------------------------
# Summary file creation
#---------------------------------------------------------------------------------------------------

# count the reads in fastq files
read_count () {
  fastq=$1
  total=$(wc -l $fastq | awk -F' ' '{print $1}')
  expr $total / 4
}

# create summary file
summary_csv="${sequencing_run}/artic-pipeline/summary.csv"


# header for summary file
echo "#SAMPLE,1-barcode-dumux reads,2-single-fastq reads,3-normalize reads,4-consensus reads" > "$summary_csv"


# populate summary file
while read barcode name; do

sample=$read_$barcode

# read count 1-barcode-demux
one_total=0        
for i in $demux_dir/$barcode/*.fastq; do
    fastq_count=$(read_count $i)
    let one_total=one_total+$fastq_count
done

# read count 2-single-fastq
two_total=$(read_count $gather_dir/${name}_${barcode}.fastq)

# read count 3-normalize
three_total=$(read_count $normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.fq)

# read count 4-consensus
bam_prefix="$consensus_dir"/${name}_${barcode}-${pipeline_label}

# counting only mapped (primary aligned) reads
four_total=$($samtools_path/samtools view -c -F 260 ${bam_prefix}.primertrimmed.rg.sorted.bam)

echo "${sample},${one_total},${two_total},${three_total},${four_total}" >> "$summary_csv"

done < "$manifest"


#---------------------------------------------------------------------------------------------------
# End
#---------------------------------------------------------------------------------------------------

echo_log "run complete"
chgrp -R 5102 $demux_dir $gather_dir $normalize_dir $consensus_dir 
