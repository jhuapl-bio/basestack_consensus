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
        echo -e "   -i      /full/path/to/fastq_of_interest.fq"
        echo -e "   -t      number of threads (default: 6)"
        echo -e "   -r      reference_genome.fasta"
        echo -e "OUTPUT:"
        echo -e "human_subtracted.fastq, human_subtracted.sam, human_subtracted.fast5"
}

#---------------------------------------------------------------------------------------------------
#default threads
threads=6
#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:t:r:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) fastq=$OPTARG ;;
                t) threads=$OPTARG ;;
                r) reference=$OPTARG ;;
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
sequencing_run=$(dirname $(dirname $(dirname "$fastq")))

#===================================================================================================
# Default values
#===================================================================================================

# input files, these files should be in the sequencing run directory, leaving samfile quoted for future updates 
base=$(basename $fastq)
sample_name=${base%.fastq}
#sample_name=$(basename "${samfile%.covfiltered.sam}")

# Output files and directory
align_out="${sample_name}.human_aligned.sam"
read_ids="${sample_name}-read_ids.txt"
out_sam="${sample_name}.human_removed.sam"

outdir="${sequencing_run}/artic-pipeline/fast5-subset-human-filtered/{sample_name}"
mkdir -p outdir

# log file
logfile="${outdir}"/logs/fast5-subset-human-filtered"${sample_name}"-$(date +"%F-%H%M%S").log

#===================================================================================================
# QUALITY CHECKING
#===================================================================================================

if [ ! -d "${sequencing_run}" ];then
    >&2 echo "Error: Sequencing run ${sequencing_run} does not exist"
    exit 1
fi

if [ -s "$outdir/$name-human-filtered-subset.fast5" ];then
    >&2 echo "Fast5 subset already exists for this sample: $outdir/${sample_name}-human-filtered-subset.fast5"
    >&2 echo "    Archive previous fast5 subset processing before rerunning."
    exit 1
else
    mkdir -p "$(dirname $outdir)/logs"
fi

#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "SAMPLE ${sample_name}: ------ Fast5 Subset Human Filter Paramters:"
echo_log "SAMPLE ${sample_name}: samtools version $(samtools --version)"
echo_log "SAMPLE ${sample_name}: minimap version $(minimap2 --version)"
echo_log "SAMPLE ${sample_name}: sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "SAMPLE ${sample_name}: sample fastq: ${fastq}"
echo_log "SAMPLE ${sample_name}: output directory: ${outdir}"
echo_log "SAMPLE ${sample_name}: ------ processing fast5 subset output ------"

#---------------------------------------------------------------------------------------------------
# module 3 Human Filter
#---------------------------------------------------------------------------------------------------

# if input is samfile, convert to fastq
#samtools fastq "${out_sam}" > "${out_sam%.sam}.fq"

echo_log "starting minimap2"

minimap2 -a \
	-x map-ont \
	-t 32 \
	"${reference}" \
	"${fastq}" > "${outdir}"/"${align_out}"

echo_log "extracting unmapped reads"
# extract unmapped reads
samtools view -f 4 "${outdir}"/"${align_out}" > "${outdir}"/"${out_sam}"

echo_log "converting to fastq"
# fastq conversion - if we want to output a fastq
samtools fastq "${outdir}"/"${out_sam}" > "${outdir}"/"${out_sam%.sam}.fq"

echo_log "retrieving read ids"
# retrieve read ids that pass
awk '{if ( $1 ~ "^@" ){}else{print $1}}' "${outdir}"/"${out_sam}" > "${outdir}"/"${read_ids}"

echo_log "performing fast5_subset"
fast5_subset \
--input "${sequencing_run}/fast5_pass" \
--save_path "${outdir}" \
--read_id_list "${outdir}"/"${read_ids}" \
--batch_size 100 \
-t $threads \
--recursive

#---------------------------------------------------------------------------------------------------

echo_log "SAMPLE ${sample_name}: Module 3 post-normalization human filter fast5 subsetting complete"

