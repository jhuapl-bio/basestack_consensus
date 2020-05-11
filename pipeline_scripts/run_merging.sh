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

REF="/home/idies/workspace/covid19/ncov_reference/sequence.fasta"

for bamfile in `ls $DIR/*.nanopolish.primertrimmed.rg.sorted.bam`; do

    sample=${bamfile##*/}
    samplename=${sample%%_*}

    # loop through all non-NTC samples
    if [ ! "$samplename" = "NTC" ]; then
	continue
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

        prefix=`echo $bamfile`
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
    else
	prefix=`echo $bamfile`
        prefix=${prefix##*/}
        prefix=${prefix%%.*}
        samtools mpileup --reference $REF $bamfile -o $DIR/$prefix.mpileup
    fi
done

