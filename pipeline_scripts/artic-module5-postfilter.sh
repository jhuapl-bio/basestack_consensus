#!/bin/bash
source /opt/basestack_consensus/bashrc
conda activate jhu-ncov

# usage function
usage() {
        echo -e "usage: ${YELLOW}$0${NC} [options]"
        echo -e ""
        echo -e "OPTIONS:"
        echo -e "   -h      show this message"
        echo -e "   -i      /full/path/to/sequencing_run/artic-pipeline/4-draft-consensus"
        echo -e "   -d      /full/path/to/<control>nanopolish.primertrimmed.rg.sorted.del.depth"
        echo -e "   -b      /full/path/to/<control>nanopolish.primertrimmed.rg.sorted.bam"
        echo -e "   -v      /full/path/to/approx_global_diversity.tsv"
        echo -e "   -c      /full/path/to/variant_case_definitions.txt"
        echo -e "   -r      /full/path/to/<reference>.fasta"
        echo -e "   -a      /full/path/to/amplicons"
        echo -e "   -p      /full/path/to/homopolymer_positions.txt"
        echo -e "   -k      /full/path/to/key_positions.txt"
        echo -e "   -m      /full/path/to/sequencing_run/manifest.txt"
        echo -e "   -n      name of control sample in manifest (Default = 'NTC')"
        echo -e ""
}

#---------------------------------------------------------------------------------------------------
control_name="NTC"
#---------------------------------------------------------------------------------------------------


# parse input arguments
while getopts "hi:d:b:v:c:r:a:p:k:m:n:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) consensus_dir=$OPTARG ;;
                d) depthfile=$OPTARG ;;
                b) bamfile=$OPTARG ;;
                v) global_vars=$OPTARG ;;
                c) case_defs=$OPTARG ;;
                r) reference=$OPTARG ;;
                a) amplicons=$OPTARG ;;
                p) homopolymers=$OPTARG ;;
                k) key_positions=$OPTARG ;;
        m) manifest=$OPTARG ;;
        n) control_name=$OPTARG;;
                ?) usage; exit ;;
       esac
done

postfilter_dir="$(dirname ${consensus_dir})/5-post-filter"

# make and save output directory
if [ ! -d "$postfilter_dir" ]; then
        mkdir "$postfilter_dir"
fi

# save path to reference genome and consensus
reference="${reference}"
consensus_dir="${consensus_dir}"

# save path to NTC depthfile and mpileup
ntc_depthfile="${depthfile}"
ntc_bamfile="${bamfile}"

# save path to global vars
global_vars="${global_vars}"

# save path to case definitions
case_defs="${case_defs}"

# save path to amplicon sites file
amplicons="${amplicons}"

# save path to homopolymers file
homopolymers="${homopolymers}"

# save path to key positions file
keypos="${key_positions}"


while read barcode name; do

    # loop through all NTC samples
    if [[  "$name" != "$control_name" ]]; then

        # align sample to reference genome
        echo "SAMPLE $name: aligning to reference genome"
        cat "${reference}" "${consensus_dir}/${name}_${barcode}.nanopolish.consensus.fasta" > "${postfilter_dir}/${name}_${barcode}.ref.fasta"
        mafft --preservecase "${postfilter_dir}/${name}_${barcode}.ref.fasta" > "${postfilter_dir}/${name}_${barcode}.align.ref.fasta"


        echo "SAMPLE $name: running vcf_postfilter.py"
        vcffile="${consensus_dir}/${name}_${barcode}.all_callers.combined.vcf"
        depth="${consensus_dir}/${name}_${barcode}.nanopolish.primertrimmed.rg.sorted.del.depth"
        alignment="${postfilter_dir}/${name}_${barcode}.align.ref.fasta"

        # run script
        vcf_postfilter.py \
        --vcffile "$vcffile" \
        --depthfile "$depth" \
        --aln-to-ref "$alignment" \
        --ntc-bamfile "$ntc_bamfile" \
        --ntc-depthfile "$ntc_depthfile" \
        --global-vars "$global_vars" \
        --key-vars "$keypos" \
        --homopolymers "$homopolymers" \
        --case-defs "$case_defs" \
        --amplicons "$amplicons" \
        --outdir "$postfilter_dir" \
        --samplename "${name}_${barcode}"
    fi

done < "${manifest}"

