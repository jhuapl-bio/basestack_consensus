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

# the sequencing run ID
#runID="20200410_2018_X4_FAN32204_327837a0"
runID=$1

logfile=/home/idies/workspace/Temporary/ernluaw1/scratch/log_${runID}-submit.txt

# location for the output of the medaka pipeline
working_dir="/home/idies/workspace/Temporary/ernluaw1/scratch/${runID}_analysis-submit"

# input sequecing directory
sequencing_run="/home/idies/workspace/covid19/sequencing_runs/${runID}"

# location for primer schemes
scheme_dir=/home/idies/workspace/Storage/ernluaw1/persistent/artic-ncov2019/primer_schemes

# primer protocol
protocol=nCoV-2019/V3

mkdir -p $working_dir

# Have these files be in the sequencing run directory
barcode_file=$(sed -n 's/^barcoding//p' "${sequencing_run}/run_config.txt")
manifest=${sequencing_run}/manifest.txt
fastq_dir=${sequencing_run}/fastq_pass

# Output directories
pipeline_label=hac-medaka-norm200
demux_dir=${working_dir}/1-barcode-demux
gather_dir=${working_dir}/2-single-fastq
consensus_dir=${working_dir}/3-consensus

echo -e "$(date +"%F %T") sequencing run start $runID"

# need to fix
guppy_barcoder_path=/home/idies/workspace/Storage/ernluaw1/persistent/bin/ont-guppy-cpu/bin

echo_log "Starting guppy demux"
$guppy_barcoder_path/guppy_barcoder \
	--require_barcodes_both_ends \
	-i "$fastq_dir" \
	-s "$demux_dir" \
	--arrangements_files $barcode_file
	

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

echo_log "Starting medaka"
mkdir -p $consensus_dir
while read barcode name; do
echo_log "${name}_${barcode}"
artic minion \
	--medaka \
	--normalise 200 \
	--threads 32 \
	--scheme-directory "$scheme_dir" \
	--read-file "$gather_dir"/${name}_${barcode}.fastq \
	"$protocol" "$consensus_dir"/${name}_${barcode}-${pipeline_label}
done < "$manifest"

echo_log "run complete"
