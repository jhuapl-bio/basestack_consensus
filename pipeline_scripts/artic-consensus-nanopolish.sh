# activate conda env
. "/home/idies/workspace/covid19/miniconda3/etc/profile.d/conda.sh"
conda activate artic-ncov2019

# grab input file from 3-normalization
fastq=$1

# grab environment and run info
run_path=$(dirname $(dirname $(dirname "$fastq")))
summary=$(find "$run_path" -max-depth 2 -name "*sequencing_summary*.txt")
scheme_dir=/home/idies/workspace/covid19/code/artic-ncov2019/primer_schemes
protocol=$(grep primers "$run_path/run_config.txt" | cut -f2)

# set up output folder
consensus_dir=$run_path/artic-pipeline/4-draft-consensus
mkdir -p "$consensus_dir"

out_prefix="$consensus_dir/$(basename ${fastq%.covfiltered.fq})"

# run ARTIC pipeline
artic minion \
    --normalise 1000000 \
    --threads 32 \
    --scheme-directory "$scheme_dir" \
    --read-file "$fastq" \
    --fast5-directory "$run_path/fast5_pass" \
    --sequencing-summary "$summary" \
    "$protocol" "$out_prefix"
