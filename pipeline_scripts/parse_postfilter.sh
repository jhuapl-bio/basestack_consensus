#!/bin/bash
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

# specify a sequencing run directory
RUN=$1

# make and save output directory
outdir="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter_nanopolish"
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

table=`ls /home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter_nanopolish/postfilt_all.txt`

# run script
prefix='postfilt'
$BINDIR/VariantValidator/parsetable.sh $table $outdir/$prefix

