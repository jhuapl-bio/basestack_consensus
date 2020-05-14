#!/bin/bash

# usage function
usage() {
        echo -e "usage: ${YELLOW}$0${NC} [options]"
        echo -e ""
        echo -e "OPTIONS:"
        echo -e "   -h      show this message"
        echo -e "   -i      /full/path/to/sequencing_run/artic-pipeline/4-draft-consensus"
        echo -e "   -d      /full/path/to/<control>nanopolish.primertrimmed.rg.sorted.depth"
        echo -e "   -b      /full/path/to/<control>nanopolish.primertrimmed.rg.sorted.bam"
        echo -e "   -v      /full/path/to/nextstrain_alignments.vcf"
        echo -e "   -c      /full/path/to/variant_case_definitions.csv"
        echo -e ""
}

#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:d:b:v:c:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) consensus_dir=$OPTARG ;;
                d) depthfile=$OPTARG ;;
                b) bamfile=$OPTARG ;;
                v) vcf_next=$OPTARG ;;
                c) case_defs=$OPTARG ;;
                ?) usage; exit ;;
       esac
done

postfilter_dir="$(dirname ${consensus_dir})/5-post-filter"

# make and save output directory
if [ ! -d $postfilter_dir ]; then
        mkdir $postfilter_dir
fi

# save path to NTC depthfile and mpileup
ntc_depthfile="${depthfile}"
ntc_bamfile="${bamfile}"

# save path to nextstrain vcf
vcf_next="${vcf_next}"

# save path to case definitions
case_defs="${case_defs}"

for consfile in "${consensus_dir}/*nanopolish.consensus.fasta"; do

	sample=${consfile##*/}
	samplename=${sample%%_*}

	# loop through all NTC samples
	if [ ! "$samplename" = "NTC" ]; then

		echo $samplename

		vcffile=$DIR/$samplename*all_callers.combined.vcf
		mpileup="$DIR/$samplename*mpileup"
		depth="$DIR/$samplename*nanopolish.primertrimmed.rg.sorted.depth"
		consensus="$DIR/$samplename*nanopolish.consensus.fasta"
		prefix=`echo $consensus`
		prefix=${prefix##*/}
		prefix=${prefix%%.*}

		# run script
		python /home/idies/workspace/covid19/code/ncov/pipeline_scripts/vcf_postfilter.py \
		--vcffile $vcffile \
		--mpileup $mpileup \
		--depthfile $depth \
		--consensus $consensus \
		--ntc-bamfile $ntc_bamfile \
		--ntc-depthfile $ntc_depthfile \
		--vcf-nextstrain $vcf_next \
		--case-defs $case_defs \
		--ns-snp-threshold 2 \
		--maf-flag 15 \
		--outdir $postfilter_dir \
		--prefix $prefix
	fi
done

