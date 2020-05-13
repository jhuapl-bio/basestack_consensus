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
        echo -e "   -i      /full/path/to/sequencing_directory"
        echo -e "   -c      name of control sample provided in manifest (Default: 'NTC')"
        echo -e ""
}

#---------------------------------------------------------------------------------------------------
# default control name
control_name="NTC"
#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:c:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) sequencing_run=$OPTARG ;;
		c) control_name=$OPTARG ;;
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
# Default values
#===================================================================================================

# input files
sequencing_run_name=$(basename $sequencing_run)
consensus_dir=${sequencing_run}/artic-pipeline/4-draft-consensus
manifest=${sequencing_run}/manifest.txt

# output directory
postfilter_dir="${sequencing_run}/artic-pipeline/5-post-filter"

# Postfiltering reference files
vcf_next="/home/idies/workspace/covid19/nextstrain/latest/alpha/alignments.vcf"
case_defs="/home/idies/workspace/covid19/code/ncov/pipeline_scripts/variant_case_definitions.csv"


# Posterfiltering NTC files for baselining
control_barcode=$(awk -F $'\t' -v control_name=${control_name} '$2 == control_name { print $1 }' ${manifest})
ntc_depthfile="${consensus_dir}/${control_name}_${control_barcode}.nanopolish.primertrimmed.rg.sorted.depth"
ntc_bamfile="${consensus_dir}/${control_name}_${control_barcode}.nanopolish.primertrimmed.rg.sorted.bam"

# CombineVariants reference files
reference="/home/idies/workspace/covid19/ncov_reference/sequence.fasta"
reference_annotation="/home/idies/workspace/covid19/ncov_reference/genes.gff3"

# Pangolin reference directory
pangolin_data="/home/idies/workspace/covid19/ncov_reference/lineages/lineages/data

# snpEff reference file
snpEff_config="/home/idies/workspace/covid19/ncov_reference/snpEff.config"

# location of programs used by pipeline
software_path=/home/idies/workspace/covid19/code
JAVA_PATH="${software_path}/jdk-14.0.1/bin"

# log file
logfile=${postfilter_dir}/logs/module5-postfilter-${sequencing_run_name}-$(date +"%F-%H%M%S").log

#===================================================================================================
# QUALITY CHECKING
#===================================================================================================

if [ ! -d ${sequencing_run} ];then
    >&2 echo "Error: Sequencing run ${sequencing_run} does not exist"
    exit 1
fi

if [ ! -d ${consensus_dir} ];then
    >&2 echo "Error: Require Module 4 draft consensus output"
    >&2 echo "    ${consensus_dir} does not exist"
    exit 1
fi

if [ ! -s ${manifest} ];then
    >&2 echo "Error: Require Run Manifest for processing"
    >&2 echo "    ${manifest} does not exist"
    exit 1
fi

module4_complete_flag="TRUE"
while IFS=$'\t' read barcode name; do
	if [ ! -f ${consensus_dir}/module4-${name}_${barcode}.all_callers.complete ];then
	    >&2 echo "Error: Module 4 has not completed for sample ${name}. Module 4 must be completed for all samples before running Module 5"
	    >&2 echo "    Details: Module 4 .complete file not found: ${consensus_dir}/module4-${name}_${barcode}.all_callers.complete"
            module4_complete_flag="FALSE"
            exit 1
	fi
done < ${manifest}

ref_files_found_flag="TRUE"
if [[ -s ${vcf_next} ]] 
	>&2 echo "Error: Module 5 reference file not found: ${vcf_next}"
	ref_files_found_flag="FALSE"
	exit 1
fi

if [[ -s ${case_defs} ]] 
	>&2 echo "Error: Module 5 reference file not found: ${case_defs}"
	ref_files_found_flag="FALSE"
	exit 1
fi

if [[ -f ${ntc_depthfile} ]] 
	>&2 echo "Error: Control sample depth file not found: ${ntc_depthfile}"
	ref_files_found_flag="FALSE"
	exit 1
fi

if [[ -f ${ntc_bamfile} ]] 
	>&2 echo "Error: Control Sample BAM file not found: ${ntc_bamfile}"
	ref_files_found_flag="FALSE"
	exit 1
fi

# create logs folder after successful completion of Module 4 for all samples
if [[ ${module4_complete_flag} == "TRUE" ]] && [[ ${ref_files_found_flag} == "TRUE" ]]; then
	mkdir -p ${postfilter_dir}/logs

#===================================================================================================
# MAIN BODY
#===================================================================================================

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

echo_log "RUN ${sequencing_run_name}: ------ Samtools / Merge Paramters:"
echo_log "RUN ${sequencing_run_name}: sequencing run folder: ${CYAN}$sequencing_run${NC}"
echo_log "RUN ${sequencing_run_name}: recording software version numbers..."
echo_log "RUN ${sequencing_run_name}: Software version: $(samtools --version)"
echo_log "RUN ${sequencing_run_name}: Reference fasta: ${reference}"
echo_log "RUN ${sequencing_run_name}: input fasta file: ${normalized_fastq}"
echo_log "RUN ${sequencing_run_name}: input vcfs: ${input_nanopolish_vcf}, ${input_medaka_vcf_zip}"
echo_log "RUN ${sequencing_run_name}: iput bam file: ${input_nanopolish_bamfile}"
echo_log "RUN ${sequencing_run_name}: output consensus directory: ${consensus_dir}"
echo_log "RUN ${sequencing_run_name}: ------ processing Samtools/Merge ------"

#---------------------------------------------------------------------------------------------------
# Module 5 Postfiler
#---------------------------------------------------------------------------------------------------

# Run module
echo_log "RUN ${sequencing_run_name}: Starting Module 5 Postfilter Submodule 1 on ${sequencing_run}"

bash -x ${software_path}/run_postfilter.sh ${sequencing_run_name} 2>> ${logfile}

#---------------------------------------------------------------------------------------------------
# Module 5 Postfiler Summarization
#---------------------------------------------------------------------------------------------------

post_filtering_complete_flag="TRUE"
while read barcode name; do
	if  [[ name != ${control_name} ]]; then
		if [ ! -s ${postfilter_dir}/${name}_${barcode}*variant_data.txt ]; then
			echo_log "RUN ${sequencing_run_name}: Error: Postfiltering must be completed for all samples prior to summarization."
			echo_log "RUN ${sequencing_run_name}:     {postfilter_dir}/${name}_${barcode}*variant_data.txt does not exist"
			post_filtering_complete_flag="FALSE"
		fi
	fi
done < ${manifest}

# Run module
if [[ ${post_filtering_complete_flag} == "TRUE" ]]; then
	echo_log "RUN ${sequencing_run_name}: Module 5 Postfiltering completed for ${sequencing_run}"
	echo_log "RUN ${sequencing_run_name}: Starting Module 5 Postfilter Summarization on ${sequencing_run}"

	python ${software_path}/postfilter_summary.py $postfilter_dir 2>> ${logfile}

	echo_log "RUN ${sequencing_run_name}: Module 5 Postfilter Summarization completed for ${sequencing_run}"
else
	echo_log "RUN ${sequencing_run_name}: Error: Module 5 Summarization not performed."
	exit 1
fi
#---------------------------------------------------------------------------------------------------
# module 5 Combine Variants
#---------------------------------------------------------------------------------------------------

if [[ ! -s ${reference_annotation} ]]; then
	echo_log "RUN ${sequencing_run_name}: Error: nCoV reference GFF not found: ${reference_annotation}"
	exit 1
fi

if [[ ! -s ${reference} ]]; then
	echo_log "RUN ${sequencing_run_name}: Error: nCoV reference FASTA not found: ${reference}"
	exit 1
fi

echo_log "RUN ${sequencing_run_name}: Combing variants for each sample..."
bash -x ${software_path}/run_combine.sh ${sequencing_run_name} 2>> ${logfile}

#---------------------------------------------------------------------------------------------------
# Module 5 Pangolin and snpEff
#---------------------------------------------------------------------------------------------------

if [[ ! -s ${snpEff_config} ]]; then
	echo_log "RUN ${sequencing_run_name}: Error: snpEff Config file not found: ${snpEff_config}"
	exit 1
fi

if [[ ! -d ${pangolin_data} ]]; then
	echo_log "RUN ${sequencing_run_name}: Error: Pangolin data directory not found: ${pangolin_data}"
	exit 1
fi

combine_variants_complete_flag="TRUE"
while read barcode name; do
	if  [[ name != ${control_name} ]]; then
		if [ ! -s ${postfilter_dir}/${name}_${barcode}.consensus.combined.vcf ]; then
			echo_log "RUN ${sequencing_run_name}:Error: Variants must be combined for all samples prior to running Pangolin and snpEff."
			echo_log "RUN ${sequencing_run_name}:     {postfilter_dir}/${name}_${barcode}.consensus.combined.vcf does not exist"
			combine_variants_complete_flag="FALSE"
		fi
	fi
done < ${manifest}

# Run pangolin and snpEff
if [[ ${combine_variants_complete_flag} == "TRUE" ]]; then
	echo_log "RUN ${sequencing_run_name}: Starting Module 5 Pangolin and snpEff on ${sequencing_run}"

	bash -x ${software_path}/run_pangolin.sh 2>> ${logfile}
	bash -x ${software_path}/run_snpEff.sh 2>> ${logfile}

else
	echo_log "RUN ${sequencing_run_name}: Error: Module 5 Pangolin and snpEff not performed."
	exit 1
fi

#---------------------------------------------------------------------------------------------------
# Post-processing quality checking

module5_complete_flag="TRUE"
if [[ ! -s ${postfiler_dir}/postfilt_consensus_all.fasta ]]; then
	echo_log "RUN ${sequencing_run_name}: Error: Panglolin output file - ${postfiler_dir}/postfilt_consensus_all.fasta not detected."
	module5_complete_flag="FALSE"
fi

if [[ ! -s ${postfiler_dir}/final_snpEff_report.txt ]]; then
	echo_log "RUN ${sequencing_run_name}: Error: snpEff output file - ${postfiler_dir}/final_snpEff_report.txt not detected."
	module5_complete_flag="FALSE"
fi

if [[ ${module5_complete} == "TRUE" ]]; then 
	echo_log "RUN ${sequencing_run_name}: Module 5 Postfilter completed for ${sequencing_run}"
	echo_log "RUN ${sequencing_run_name}: Creating ${postfilter_dir}/module5-${sequencing_run_name}.complete"
	touch ${postfilter_dir}/module5-${sequencing_run_name}.complete
else
	echo_log "RUN ${sequencing_run_name}: Error: Module 5 did not complete."
	exit 1
fi
