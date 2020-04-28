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

#===================================================================================================
# SET DEFAULT VALUES
#===================================================================================================

# set default values here
logfile="/dev/null"
tempdir="/tmp"

# report current hash of miseq-analysis git repo
bin_path="$(dirname $0)"
#script_path=$(dirname $(readlink -f $0))
GIT_DIR="$script_path/../.git"
export GIT_DIR
hash=$(git rev-parse --short HEAD)

primerscheme_path="$bin_path/../../artic-ncov2019/primer_schemes"
protocol="nCoV-2019/V3"
reference="$primerscheme_path/$protocol/nCoV-2019.reference.fasta"

stats_base="artic-pipeline/run_stats"
demux_base="artic-pipeline/1-guppy-barcoder"
lengthfilter_base="artic-pipeline/2-guppyplex"
normalize_base="artic-pipeline/3-normalization_update"
draftconsensus_base="artic-pipeline/4-draft-consensus_update"
postfilter_base="artic-pipeline/5-post-filter_update"

#===================================================================================================
# PARSE INPUT ARGUMENTS
#===================================================================================================

# parse input arguments
while getopts "hi:m:b:p:r:o:1:2:3:4:5:6:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) run_path=$OPTARG ;;
		m) manifest=$OPTARG ;;
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

if ! [[ -d "$run_path" ]]; then
	echo -e "${RED}Error: run path ${CYAN}$run_path${RED} does not exist.${NC}"
	usage
	exit
fi

if ! [[ -s "$manifest" ]]; then
	manifest="$run_path/manifest.txt"
	if ! [[ -s "$manifest" ]]; then
		echo -e "${RED}Error: manifest file ${CYAN}$manifest${RED} does not exist.${NC}"
		usage
		exit
	fi
fi

if ! [[ -s "$reference" ]]; then
	echo -e "${RED}Error: reference sequence ${CYAN}$reference${RED} does not exist.${NC}"
	usage
	exit
else
	ref_header=$(head -n1 "$reference" | cut -c2-)
	ref_length=$("$bin_path/fix_fasta.sh" "$reference" | tail -n1 | awk '{print length($1)}')
fi

if ! [[ -d "$bin_path" ]]; then
	echo -e "${RED}Error: reference sequence ${CYAN}$reference${RED} does not exist.${NC}"
	usage
	exit
fi
vcfigv_repo_path="$bin_path/../../vcfigv"
if ! [[ -d "$vcfigv_repo_path" ]]; then
	echo -e "${RED}Error: vcfigv repository ${CYAN}$vcfigv_repo_path${RED} does not exist.${NC}"
	usage
	exit
fi

if [[ -z "$stats_path" ]]; then
	stats_path="$run_path/$stats_base"
fi
if ! [[ -d "$stats_path" ]]; then
	echo -e "${RED}Error: output path ${CYAN}$stats_path${RED} does not exist.${NC}"
	usage
	exit
fi
if [[ -z "$demux_path" ]]; then
	demux_path="$run_path/$demux_base"
fi
if ! [[ -d "$demux_path" ]]; then
	echo -e "${RED}Error: demux path ${CYAN}$demux_path${RED} does not exist.${NC}"
	usage
	exit
fi
if ! [[ -s "$demux_path/barcoding_summary.txt" ]]; then
	echo -e "${RED}Error: demux summary ${CYAN}$demux_path/barcoding_summary.txt${RED} does not exist.${NC}"
	usage
	exit
fi
if [[ -z "$lengthfilter_path" ]]; then
	lengthfilter_path="$run_path/$lengthfilter_base"
fi
if ! [[ -d "$lengthfilter_path" ]]; then
	echo -e "${RED}Error: length filter path ${CYAN}$lengthfilter_path${RED} does not exist.${NC}"
	usage
	exit
fi
if [[ -z "$normalize_path" ]]; then
	normalize_path="$run_path/$normalize_base"
fi
if ! [[ -d "$normalize_path" ]]; then
	echo -e "${RED}Error: normalization path ${CYAN}$normalize_path${RED} does not exist.${NC}"
	usage
	exit
fi
if [[ -z "$draftconsensus_path" ]]; then
	draftconsensus_path="$run_path/$draftconsensus_base"
fi
if ! [[ -d "$draftconsensus_path" ]]; then
	echo -e "${RED}Error: draft consensus path ${CYAN}$draftconsensus_path${RED} does not exist.${NC}"
	usage
	exit
fi
if [[ -z "$postfilter_path" ]]; then
	postfilter_path="$run_path/$postfilter_base"
fi
if ! [[ -d "$postfilter_path" ]]; then
	echo -e "${RED}Error: post-filter path ${CYAN}$postfilter_path${RED} does not exist.${NC}"
	usage
	exit
fi
postfilt_summary="$postfilter_path/postfilt_summary.txt"
if ! [[ -s "$postfilt_summary" ]]; then
	echo -e "${RED}Error: post-filter summary ${CYAN}$postfilt_summary${RED} does not exist.${NC}"
	usage
	exit
fi
postfilt_all="$postfilter_path/postfilt_all.txt"
if ! [[ -s "$postfilt_summary" ]]; then
	echo -e "${RED}Error: post-filter full report ${CYAN}$postfilt_all${RED} does not exist.${NC}"
	usage
	exit
fi

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

outfile="$stats_path/summary.txt"
demuxfile="$stats_path/demux_count.txt"
depthfile="$stats_path/depth-all.txt"
mutations_pos="$stats_path/mutations-pos.txt"
mutations_all="$stats_path/mutations-all.txt"
mutations_table="$stats_path/mutations-table.txt"

echo_log "====== Call to ${YELLOW}"$(basename $0)"${NC} from ${GREEN}"$(hostname)"${NC} ======"

# create directory to hold temporary files
runtime=$(date +"%Y%m%d%H%M%S%N")
workdir="$tempdir/report_summary_table-$runtime"
mkdir -m 775 -p "$workdir"

echo_log "recording software version numbers"
echo_log "current git hash: $hash"
echo_log "  guppy barcoder: "
echo_log "input arguments"
echo_log "  sequencing run folder: ${CYAN}$run_path${NC}"
echo_log "    1)    barcode demux: ├── ${CYAN}${demux_path#$run_path}${NC}"
echo_log "    2)    length filter: ├── ${CYAN}${lengthfilter_path#$run_path}${NC}"
echo_log "    3)    normalization: ├── ${CYAN}${normalize_path#$run_path}${NC}"
echo_log "    4)  draft consensus: ├── ${CYAN}${draftconsensus_path#$run_path}${NC}"
echo_log "    5)       nextstrain: ├── ${CYAN}${nextstrain_path#$run_path}${NC}"
echo_log "    6)      post-filter: └── ${CYAN}${postfilter_path#$run_path}${NC}"
echo_log "  manifest: ${CYAN}$manifest${NC}"
echo_log "  reference sequence: ${CYAN}$reference${NC}"
echo_log "  working directory: ${CYAN}$workdir${NC}"
echo_log "  threads: ${CYAN}1${NC}"
echo_log "output arguments"
echo_log "  log file: ${CYAN}$logfile${NC}"
echo_log "  summary file: ${CYAN}${summary#$run_path}${NC}"
echo_log "  demux file: ${CYAN}${demuxfile#$run_path}${NC}"
echo_log "  depth file: ${CYAN}${depthfile#$run_path}${NC}"
echo_log "  mutation file: ${CYAN}${mutations_table#$run_path}${NC}"
echo_log "------ processing pipeline output ------"

if ! [[ -s "$outfile" ]]; then
	make_new_outfile="true"
else
	echo_log "Summary file already present - will not overwrite it"
	make_new_outfile="false"
fi
if [[ "$make_new_outfile" == "true" ]]; then
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		"Sample" \
		"Barcode" \
		"Raw Reads" \
		"Length filter" \
		"SARS-CoV-2 Aligned" \
		"Coverage" \
		"Consensus" > "$outfile"
fi

if ! [[ -s "$stats_path/demux_count.txt" ]]; then
	echo_log "Pulling barcode demux stats"
	tail -n+2 "$demux_path/barcoding_summary.txt" | cut -f2 | sort | uniq -c > "$demuxfile"
else
	echo_log "Demux count already present - will not overwrite it"
fi

echo_log "Reading manifest"
while read barcode label; do

	echo_log "  $barcode"

	demux_reads=$(grep "$barcode" "$stats_path/demux_count.txt" | sed 's/^ \+//' | cut -d" " -f1)
	gather=$(find "$lengthfilter_path" -name "*$barcode.fastq")
	echo_log "    FASTQ file: $gather"
	if [[ -s "$gather" ]]; then
		length_filter=$(($(wc -l < "$gather") / 4))
	else
		length_filter=0
	fi
	echo_log "      number of reads: $length_filter"
#	alignment=$(find "$draftconsensus_path" -name "*${barcode}*.sorted.bam" ! -name "*trimmed*")
	alignment=$(find "$normalize_path" -name "*$barcode.sam")
	if [[ -s "$gather" ]]; then
		alignment_depth=$(wc -l "$alignment")
	else
		alignment_depth=0
	fi
	echo_log "    full alignment: $alignment"
	echo_log "      aligned reads: $alignment_depth"
	normalized_alignment=$(find "$draftconsensus_path" -name "*$barcode*.sorted.bam" ! -name "*trimmed*")
	normalized_reads=$(samtools idxstats "$normalized_alignment" | grep "$ref_header" | cut -f3)
	echo_log "    normalized alignment: $alignment"
	echo_log "      normalized reads: $alignment_depth"
	draft_consensus=$(find "$draftconsensus_path" -name "*$barcode*.consensus.fasta")
	echo_log "    draft consensus: $draft_consensus"
	consensus_length=$(($ref_length - $(tail -n+2 "$draft_consensus" | grep -o N | wc -l)))
	echo_log "      unambiguous consensus: $consensus_length"
	post_filter=$(find "$postfilter_path" -name "*$barcode*.variant_data.txt")
	echo_log "    variant file: $post_filter"

	if [[ "$make_new_outfile" == "true" ]]; then
		echo_log "  adding line to summary file"
		flag=$(grep "^$label" "$postfilt_summary" | cut -d$'\t' -f6)
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

	depth_outfile="$stats_path/depth-${label}_${barcode}.txt"
	if ! [[ -s "$depth_outfile" ]]; then
		echo_log "  creating depth file"
		samtools sort "$alignment" | samtools depth -a -d 10000000 - > "$depth_outfile"
	fi

	normalized_depth_outfile="$stats_path/depth-${label}_${barcode}.txt"
	if ! [[ -s "$normalized_depth_outfile" ]]; then
		echo_log "  creating normalized depth file"
		samtools depth -d 0 -a "$normalized_alignment" > "$normalized_depth_outfile"
	fi

	mutations_outfile="$stats_path/mutations-${label}_${barcode}.txt"
	if ! [[ -s "$mutations_outfile" ]]; then
		echo_log "  creating mutations file"
		final_consensus=$(find "$postfilter_path" -name "*$barcode*.complete.fasta")
		echo_log "    final consensus: $final_consensus"
		if [[ -s "$final_consensus" ]]; then
			"$bin_path/mutations.sh" \
				$("$bin_path/fix_fasta.sh" "$reference" | tail -n1) \
				$("$bin_path/fix_fasta.sh" "$final_consensus" | tail -n1) \
				| tail -n+55 | head -n-67 > "$mutations_outfile"
		fi
	fi

	trimmed_alignment=$(find "$draftconsensus_path" -name "*${barcode}*.primertrimmed.rg.sorted.bam")
	echo_log "    primer trimmed alignment: $trimmed_alignment"
	vcf=$(find "$draftconsensus_path/VariantValidator" -name "*${barcode}*.merged.vcf" ! -name "*medaka*")
	outPrefix=$(basename "${vcf%.merged.vcf}")

	if [[ -s "$vcf" && -s "$trimmed_alignment" ]]; then

		igv_out_path="$stats_path/igv"
		mkdir -p "$igv_out_path"

		java -cp "$vcfigv_repo_path/src" \
			Vcf2Bat \
			--squish \
			--nocombine \
			--svg \
			aln="$trimmed_alignment" \
			var="$vcf" \
			genome="$reference" \
			outprefix="$outPrefix"

		mv "$outPrefix.bat" "$outPrefix"
		mv "$outPrefix" "$igv_out_path"

		"$vcfigv_repo_path/xvfb-run" \
			--auto-servernum "$vcfigv_repo_path/IGV_2.3.98/igv.sh" \
			-b "$pipelinepath/$outPrefix/$outPrefix.bat"

	fi

done < "$manifest"

if [[ "$make_new_outfile" == "true" ]]; then
	printf "\tuncalled\t%'d\tNA\tNA\tNA\tNA\n" $(grep unclassified "$stats_path/demux_count.txt" | sed 's/^ \+//' | cut -d" " -f1) >> "$outfile"
fi

echo_log "Calculating depth"
find "$stats_path" -name "depth-*.txt" ! -name "depth-all.txt" -print0 | while read -d $'\0' f; do
	base="${f%.txt}"
	awk -v BASE=$(basename "${base#depth-}") '{printf("%s\t%s\n", BASE, $0);}' "$f"
done > "$depthfile"

echo_log "Identifying mutations"
rm "$mutations_all"
find "$stats_path" -name "mutations-*.txt" ! -name "mutations-pos.txt" ! -name "mutations-all.txt" ! -name "mutations-table.txt" | while read fn; do
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

cp "$postfilt_summary" "$stats_path"
cp "$postfilt_all" "$stats_path"

echo_log "${GREEN}Done${NC}"
