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
	echo -e "   -o      output folder (default: ${CYAN}<run-folder>/artic-pipeline/run_stats${NC})"
	echo -e "   -1      barcode demux folder (default: ${CYAN}<run-folder>/artic-pipeline/1-guppy-barcoder${NC})"
	echo -e "   -2      length filter folder (default: ${CYAN}<run-folder>/artic-pipeline/2-guppyplex${NC})"
	echo -e "   -3      normalization folder (default: ${CYAN}not used yet${NC})"
	echo -e "   -4      draft consensus folder (default: ${CYAN}<run-folder>/artic-pipeline/3-hac-medaka-norm200${NC})"
	echo -e "   -5      nextstrain folder (default: ${CYAN}not used yet${NC})"
	echo -e "   -6      post-filter folder (default: ${CYAN}<run-folder>/artic-pipeline/4-post-filter${NC})"
	echo -e ""
}
#---------------------------------------------------------------------------------------------------
# set default values here
logfile="/dev/null"
tempdir="/tmp"

# report current hash of miseq-analysis git repo
#script_path=$(dirname $(readlink -f $0))
#GIT_DIR="$script_path/.git"
#export GIT_DIR
#hash=$(git rev-parse --short HEAD)

#echo $GIT_DIR

#base_path="/sciserver/vc_crypt/covid19/vc1/sequencing_runs"
#base_path="/home/idies/workspace/covid19/sequencing_runs"
#bin_path="$base_path/../code/ncov/pipeline_scripts"

primerscheme_path="$base_path/../code/artic-ncov2019/primer_schemes"
protocol="nCoV-2019/V3"
reference="$primerscheme_path/$protocol/nCoV-2019.reference.fasta"
ref_header=$(head -n1 "$reference" | cut -c2-)
ref_length=$("$bin_path/fix_fasta.sh" "$reference" | tail -n1 | awk '{print length($1)}')

stats_base="artic-pipeline/run_stats"
demux_base="artic-pipeline/1-guppy-barcoder"
lengthfilter_base="artic-pipeline/2-guppyplex"
draftconsensus_base="artic-pipeline/3-hac-medaka-norm200"
postfilter_base="artic-pipeline/4-post-filter"

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

stats_base="artic-pipeline/run_stats"
demux_base="artic-pipeline/1-guppy-barcoder"
lengthfilter_base="artic-pipeline/2-guppyplex"
draftconsensus_base="artic-pipeline/3-hac-medaka-norm200"
postfilter_base="artic-pipeline/4-post-filter"

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

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

# create directory to hold temporary files
runtime=$(date +"%Y%m%d%H%M%S%N")
workdir="$tempdir/build_krakendb-$runtime"
mkdir -m 775 -p "$workdir"

echo_log "recording software version numbers"
echo_log "current git hash: $hash"
echo_log "  guppy barcoder: "
echo_log "input arguments"
echo_log "  sequencing run folder: ${CYAN}$run_path${NC}"
echo_log "  working directory: ${CYAN}$workdir${NC}"
echo_log "  threads: ${CYAN}1${NC}"
echo_log "output arguments"
echo_log "  log file: ${CYAN}$logfile${NC}"
echo_log "  summary file: ${CYAN}$summary${NC}"
echo_log "  depth file: ${CYAN}$summary${NC}"
echo_log "  mutation file: ${CYAN}$summary${NC}"
echo_log "------ processing pipeline output ------"

exit

outfile="$run_path/$stats_path/summary.txt"

if ! [[ -s "$outfile" ]]; then
	make_new_outfile="true"
else
	make_new_outfile="false"
fi
if [[ "$make_new_outfile" == "true" ]]; then
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		"Sample" \
		"Barcode" \
		"Demultiplexed Reads" \
		"Length filter" \
		"SARS-CoV-2 Aligned" \
		"Consensus coverage" \
		"Consensus" > "$outfile"
fi

if ! [[ -s "$run_path/$stats_path/demux_count.txt" ]]; then
	tail -n+2 "$run_path/$demux_path/barcoding_summary.txt" | cut -f2 | sort | uniq -c > "$run_path/$stats_path/demux_count.txt"
fi

while read barcode label; do

	echo_log "  $barcode"

	demux_reads=$(grep "$barcode" "$run_path/$stats_path/demux_count.txt" | sed 's/^ \+//' | cut -d" " -f1)
	gather=$(find "$run_path/$lengthfilter_path" -name "*$barcode.fastq")
	length_filter=$(($(wc -l < "$gather") / 4))
	alignment=$(find "$run_path/$draftconsensus_path" -name "*$barcode*.sorted.bam" ! -name "*trimmed*")
	aligned_reads=$(samtools idxstats "$alignment" | grep "$ref_header" | cut -f3)
	draft_consensus=$(find "$run_path/$draftconsensus_path" -name "*$barcode*.consensus.fasta")
	consensus_length=$(($ref_length - $(tail -n+2 "$draft_consensus" | grep -o N | wc -l)))
	post_filter=$(find "$run_path/$postfilter_path" -name "*$barcode*.variant_data.txt")
	if [[ "$consensus_length" -lt 25000 ]]; then
		flag="No"
	elif ! [[ -s "$post_filter" ]]; then
		flag="No"
	else
		flag=$(tail -n+2 "$post_filter" | awk -F $'\t' '
			BEGIN{out="Yes"} {
				if($8 == "MAF>0.10" && out=="Yes"){ out="Yes*"; }
				if(length($8) > 0 && $8 != "MAF>0.10"){out="Maybe";}
			} END {
				print out;
			}')
	fi

	if [[ "$make_new_outfile" == "true" ]]; then
		printf "%s\t%s\t%'d\t%'d\t%'d\t%'d (%s %%)\t%s\n" \
			"$label" \
			"$barcode" \
			"$demux_reads" \
			"$length_filter" \
			"$aligned_reads" \
			"$consensus_length" \
			$(echo "$consensus_length" | awk -v L="$ref_length" '{printf("%0.1f", (100*$1/L))}') \
			"$flag" >> "$outfile"
	fi

	depth_outfile="$run_path/$stats_path/depth-${label}-${barcode}.txt"
	if ! [[ -s "$depth_outfile" ]]; then
		samtools depth -d 0 -a "$alignment" > "$depth_outfile"
	fi

	mutations_outfile="$run_path/$stats_path/mutations-${label}-${barcode}.txt"
	if ! [[ -s "$mutations_outfile" ]]; then
		final_consensus=$(find "$run_path/$postfilter_path" -name "*$barcode*.complete.fasta")
		if [[ -s "$final_consensus" ]]; then
			"$bin_path/mutations.sh" \
				$("$bin_path/fix_fasta.sh" "$reference" | tail -n1) \
				$("$bin_path/fix_fasta.sh" "$final_consensus" | tail -n1) \
				| tail -n+55 | head -n-67 > "$mutations_outfile"
		fi
	fi

done < "$run_path/manifest.txt"

if [[ "$make_new_outfile" == "true" ]]; then
	printf "\tuncalled\t%'d\tNA\tNA\tNA\n" $(grep unclassified "$run_path/$stats_path/demux_count.txt" | sed 's/^ \+//' | cut -d" " -f1) >> "$outfile"
fi

echo_log "Calculating depth"
find "$run_path/$stats_path" -name "depth-*.txt" -print0 | while read -d $'\0' f; do
	base="${f%.txt}"
	awk -v BASE="${base#depth-}" '{printf("%s\t%s\n", BASE, $0);}' "$f"
done > "$run_path/$stats_path/depth-all.txt"

echo_log "Identifying mutations"
mutations_pos="$run_path/$stats_path/mutations-pos.txt"
mutations_all="$run_path/$stats_path/mutations-all.txt"
mutations_table="$run_path/$stats_path/mutations-table.txt"
find "$run_path/$stats_path" -name "mutations-*.txt" | while read fn; do
	base=$(basename "${fn%.txt}")
	base="${base#mutations-}"
	awk -v BASE="$base" 'BEGIN{ printf("%s", BASE); } {
		printf("\t%s", $1);
	} END { printf("\n"); }' "$fn" >> "$mutations_all"
	cut -c2- "$fn" | rev | cut -c2- | rev
done | sort | uniq > "$mutations_pos"

awk '{
	if(NR==FNR) {
		a[$1];
	} else {
		for(i=2; i<=NF; i++) {
			pos=substr($i,2,length($i)-2);
			m[$1][pos] = $i;
		}
	}
} END {
	printf("virus");
	PROCINFO["sorted_in"] = "@ind_num_asc";
	for(i in a) { printf("\t%s", i); }
	for(sample in m) {
		printf("\n%s", sample);
		for(i in a){
			printf("\t");
			if(i in m[sample]) {
				printf("%s", m[sample][i]);
			}
		}
	}
	printf("\n");
}' "$mutations_pos" "$mutations_all" > "$mutations_table"

echo_log "${GREEN}Done${NC}"
