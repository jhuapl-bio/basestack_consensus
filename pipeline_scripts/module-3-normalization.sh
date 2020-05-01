#!/bin/bash

. "/home/idies/workspace/covid19/miniconda3/etc/profile.d/conda.sh"
conda activate artic-ncov2019

# read input arguments
fastq=$1
out_base=$2

# set up output directory
base=$(basename "${fastq%.fastq}")
sequencing_run=$(dirname $(dirname $(dirname "$fastq")))
out_dir="$out_base/$(basename $sequencing_run)/3-normalization/$base"

mkdir -p "$out_dir"

# set up environment invo
software_path=/home/idies/workspace/covid19/code
JAVA_PATH="$software_path/jdk-14.0.1/bin"
samtools_path="$software_path/samtools-1.10/bin"
NormalizeCoveragePath="$software_path/CoverageNormalization"

# get reference sequence
scheme_dir="$software_path/artic-ncov2019/primer_schemes"
protocol=$(awk '/primers/{ print $2 }' "${sequencing_run}/run_config.txt")
reference="$scheme_dir/$protocol/nCoV-2019.reference.fasta"

# align all data to reference
align_out="$out_dir/$base.sam"
minimap2 -a \
	-x map-ont \
	-t 32 \
	"$reference" \
	"$fastq" > "$align_out"

# normalization,.txt file output went to working directory
out_sam="${align_out%.sam}.covfiltered.sam"

$JAVA_PATH/java \
	-cp $NormalizeCoveragePath/src \
	NormalizeCoverage \
	input="$align_out" \
	coverage_threshold=150 \
	--qual_sort \
	--median_strand

# fastq conversion
$samtools_path/samtools fastq "$out_sam" > "${out_sam%.sam}.fq"

# BAM creation
$samtools_path/samtools sort "$out_sam" > "${out_sam%.sam}.bam"
$samtools_path/samtools index "${out_sam%.sam}.bam"
