#!/bin/bash
source /home/idies/workspace/covid19/bashrc
conda activate jhu-ncov

# usage function
usage() {
        echo -e "usage: ${YELLOW}$0${NC} [options]"
        echo -e ""
        echo -e "OPTIONS:"
        echo -e "   -h      show this message"
        echo -e "   -i      /full/path/to/sequencing_run/artic-pipeline/5-postfilter"
        echo -e "   -c      /full/path/to/<snpEff_configuration>"
        echo -e "   -m      /full/path/to/manifest.txt"
        echo -e ""
}

#---------------------------------------------------------------------------------------------------

# parse input arguments
while getopts "hi:c:m:" OPTION
do
       case $OPTION in
                h) usage; exit 1 ;;
                i) postfilter_dir=$OPTARG ;;
                c) snpEff_config=$OPTARG ;;
                m) manifest=$OPTARG ;;
                ?) usage; exit ;;
       esac
done

DBNAME="ncov"

annotate=$(which annotate_variants.sh)
if [[ -z "$annotate" ]]; then
	echo "Error: annotate_variants.sh is not in your path! Please place in your path before rerunning module."
	exit 1
fi

while read barcode name; do
    vcf="${postfilter_dir}"/"${name}"_"${barcode}".allsnps.combined.vcf
    if [[ -s "$vcf"  ]]; then
	    bash -x "${annotate}" "${vcf}" "${snpEff_config}" "${DBNAME}" "${postfilter_dir}"
	    echo "SnpEff completed on run ${postfilter_dir}"
	    echo "Making final reports on run ${postfilter_dir}"
	    cat "${postfilter_dir}"/"${name}"_"${barcode}"_ann_report.txt  | awk '$4 != "N" { print $0}'  | awk '!seen[$0]++' >> "${postfilter_dir}/final_snpEff_report.txt"
	    cat "${postfilter_dir}"/"${name}"_"${barcode}"_ann_report.txt  | awk '!seen[$0]++' | awk 'NR == 1  || $4 == "N" { print $0}'  >> "${postfilter_dir}/snpEff_report_with_Ns.txt"
    else
	    echo "File not found (snpEff not run): $vcf"
    fi
done < "$manifest"


