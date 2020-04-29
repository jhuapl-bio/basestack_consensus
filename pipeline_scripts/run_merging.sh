#!/bin/bash
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus"

REF="$BINDIR/VariantValidator/nCoV-2019.reference.fasta"

for consfile in $DIR/*.consensus.fasta; do

    sample=${consfile##*/}
    samplename=${sample%%_*}

    # loop through all non-NTC samples
    if [ ! "$samplename" = "NTC" ]; then

        echo $samplename

        nanopolishvcfunzipped=`ls /home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus/$samplename*.nanopolish.merged.vcf`
        #nanopolishvcfunzipped=${nanopolishvcffile::-3}
        #if [ ! -r $nanopolishvcfunzipped ]
        #then
        #    echo 'Unzipping '$nanopolishvcffile
        #    gunzip -c $nanopolishvcffile > $nanopolishvcfunzipped          
        #fi

        medakavcffile=`ls /home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus/$samplename*.medaka.merged.vcf.gz`
        medakavcfunzipped=${medakavcffile::-3}
        if [ ! -r $medakavcfunzipped ]
        then
            echo 'Unzipping '$medakavcffile
            gunzip -c $medakavcffile > $medakavcfunzipped          
        fi

        bamfile=`ls $DIR/$samplename*.medaka.primertrimmed.rg.sorted.bam`
        consensus=`ls $DIR/$samplename*.medaka.consensus.fasta`
        prefix=`echo $consensus`
        prefix=${prefix##*/}
        prefix=${prefix%%.*}

        echo 'Running merging'
        echo 'Reference: '$REF
        echo 'Bam file: '$bamfile
        echo 'Nanopolish VCF: '$nanopolishvcfunzipped
        echo 'Medaka VCF: '$medakavcfunzipped
        echo 'Out prefix: '$DIR/$prefix

	# run script
        $BINDIR/VariantValidator/run.sh $REF $bamfile $nanopolishvcfunzipped,$medakavcfunzipped $DIR/$prefix
    fi
done

