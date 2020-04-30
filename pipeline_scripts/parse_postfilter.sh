#!/bin/bash
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

REF="/home/idies/workspace/covid19/ncov_reference/sequence.fasta"
GENES="/home/idies/workspace/covid19/ncov_reference/genes.gff3"

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter"

# make and save output directory
outdir=$DIR

for table in $DIR/*.variant_data.txt; do

    sample=${table##*/}
    samplename=${sample%%.*}

    # loop through all non-NTC samples
    if [ ! "$samplename" = "NTC" ]; then

        echo 'Table file: '$table
        echo 'Out prefix: '$outdir/$samplename

	# run script
        $BINDIR/VariantValidator/parsetable.sh $table $GENES $REF $outdir/$samplename
    fi
done

