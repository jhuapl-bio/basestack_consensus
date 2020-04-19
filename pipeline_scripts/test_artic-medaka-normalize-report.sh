. "/home/idies/workspace/Storage/ernluaw1/persistent/Miniconda3/etc/profile.d/conda.sh"
conda activate artic-ncov2019-medaka

echo_log() {
	input="$*"
	# if input is non-empty string, prepend initial space
	if [[ -n "$input" ]]; then
		input=" $input"
	fi
	# print to STDOUT
	echo -e "[$(date +"%F %T")]$input"
	# print to log file (after removing color strings)
	echo -e "[$(date +"%F %T")]$input\r" | sed -r 's/\x1b\[[0-9;]*m?//g' >> "$logfile"
}



# Command Line Arguments

# the sequencing run ID - update*
runID="20200410_2018_X5_FAN29376_72d1837b"
#runID=$1

logfile=/home/idies/workspace/Temporary/ernluaw1/scratch/log_${runID}-test-report.txt

# location for the output of the medaka pipeline - update*
#working_dir="/home/idies/workspace/Temporary/ernluaw1/scratch/${runID}_analysis-submit"
working_dir="/home/idies/workspace/Temporary/ernluaw1/scratch/${runID}_test-normalize"

# input sequecing directory
sequencing_run="/home/idies/workspace/covid19/sequencing_runs/${runID}"

# location for primer schemes
scheme_dir=/home/idies/workspace/Storage/ernluaw1/persistent/artic-ncov2019/primer_schemes

# primer protocol
protocol=nCoV-2019/V3

mkdir -p $working_dir

# Have these files be in the sequencing run directory - update*
barcode_file=$(sed -n 's/^barcoding//p' "${sequencing_run}/run_config.txt")
#manifest=${sequencing_run}/manifest.txt
manifest="/home/idies/workspace/Temporary/ernluaw1/scratch/test_${runID}/manifest.txt"
fastq_dir=${sequencing_run}/fastq_pass

# Output directories
pipeline_label=hac-medaka-norm200
demux_dir=${working_dir}/1-barcode-demux
gather_dir=${working_dir}/2-single-fastq
normalize_dir=${working_dir}/3-normalize
consensus_dir=${working_dir}/4-consensus

echo -e "$(date +"%F %T") sequencing run start $runID"

# module 1 ################################################################################

# need to fix hardcoded path to software
guppy_barcoder_path=/home/idies/workspace/Storage/ernluaw1/persistent/bin/ont-guppy-cpu/bin

echo_log "Starting guppy demux"

$guppy_barcoder_path/guppy_barcoder \
	--require_barcodes_both_ends \
	-i "$fastq_dir" \
	-s "$demux_dir" \
	--arrangements_files $barcode_file
	


# module 2 #################################################################################

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


# module 3 ################################################################################

echo_log "Starting normalize"
while read barcode name; do

# create output directories, need separate directories for each sample 
mkdir -p $normalize_dir/${name}_${barcode}

# create alignment file
align_out=${name}_${barcode}.sam
minimap2 -a \
-x map-ont \
-t 32 \
$scheme_dir/$protocol/nCoV-2019.reference.fasta \
$gather_dir/${name}_${barcode}.fastq > $working_dir/$align_out

# normalization,.txt file output went to working directory
out_sam=$working_dir/${name}_${barcode}.covfiltered.sam
$JAVA_PATH/java -cp $NormalizeCoveragePath/src NormalizeCoverage input=$working_dir/$align_out > $normalize_dir/${name}_${barcode}

# fastq conversion
samtools fastq $out_sam > $working_dir/${name}_${barcode}.covfiltered.fq

done < "$manifest"


# module 4 ##################################################################################

echo_log "Starting medaka"
mkdir -p $consensus_dir

while read barcode name; do
echo_log "${name}_${barcode}"
artic minion \
	--medaka \
	--normalise 1000000 \
	--threads 32 \
	--scheme-directory "$scheme_dir" \
	--read-file $working_dir/${name}_${barcode}.covfiltered.fq \
	"$protocol" "$consensus_dir"/${name}_${barcode}-${pipeline_label}
done < "$manifest"


# End ######################################################################################

echo_log "run complete"
