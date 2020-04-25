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

# Output directories
pipeline_label=hac-medaka-norm200
demux_dir=${sequencing_run}/artic/1-barcode-demux
gather_dir=${sequencing_run}/artic/2-length-filter
normalize_dir=${sequencing_run}/artic/3-normalization
consensus_dir=${sequencing_run}/artic/4-draft-consensus

echo -e "$(date +"%F %T") sequencing run start $runID"

# module 1 ################################################################################

# need to fix hardcoded path to software
guppy_barcoder_path=${software_path}/ont-guppy-cpu/bin

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
	--read-file $normalize_dir/${name}_${barcode}/${name}_${barcode}.covfiltered.fq \
	"$protocol" "$consensus_dir"/${name}_${barcode}-${pipeline_label}
done < "$manifest"


# Summary file creation ####################################################################

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


# End ######################################################################################

echo_log "run complete"
chgrp -R 5102 $demux_dir $gather_dir $normalize_dir $consensus_dir 
