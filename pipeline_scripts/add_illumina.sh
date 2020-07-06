#!/bin/bash
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

javac $BINDIR/VariantValidator/src/*.java

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/4-draft-consensus"

REF="/home/idies/workspace/covid19/ncov_reference/sequence.fasta"

for bamfile in `ls $DIR/*.nanopolish.primertrimmed.rg.sorted.bam`; do

    sample=${bamfile##*/}
    samplename=${sample%%_*}

    if [ ! "$samplename" = "NTC" ]; then

        echo $samplename

        nanopolishvcf=`ls $DIR/$samplename*.nanopolish.merged.vcf`

        medakavcf=`ls $DIR/$samplename*.medaka.merged.vcf`

        prefix=`echo $bamfile`
        prefix=${prefix##*/}
        prefix=${prefix%%.*}
	samtoolsvcf=`ls $DIR/$samplename*.samtools.vcf`
	
	samplemap=$BINDIR/samplenamemap.txt
	
	illuminasample=`cat $samplemap | grep $samplename | cut -f 1`
	if [ -z $illuminasample ]
        then
            continue
        fi

	echo 'Illumina sample: '$illuminasample
	freebayesvcf=`ls /home/idies/workspace/covid19/illumina/*/trimmedv2variantcalling/$illuminasample'_'*'.freebayes.vcf' 2>/dev/null`
        ivarvcf=`ls /home/idies/workspace/covid19/illumina/*/trimmedv2variantcalling/$illuminasample'_'*'.ivar.vcf' 2>/dev/null`
	illuminasamtoolsvcf=`ls /home/idies/workspace/covid19/illumina/*/trimmedv2variantcalling/$illuminasample'_'*'.samtools.vcf' 2>/dev/null`
	illuminabam=`ls /home/idies/workspace/covid19/illumina/*/trimmedv2variantcalling/$illuminasample'_'*'.norm200.bam' 2>/dev/null`
	illuminampileup=`ls /home/idies/workspace/covid19/illumina/*/trimmedv2variantcalling/$illuminasample'_'*'.mpileup' 2>/dev/null`
	
	ontmpileup=`ls $DIR/$samplename*.mpileup 2>/dev/null`
        if [ -z $freebayesvcf ] || [ -z $ivarvcf ] || [ -z $illuminasamtoolsvcf ]
        then
            continue
        fi
	echo 'Freebayes vcf: '$freebayesvcf
	echo 'ivar vcf: '$ivarvcf
	echo 'Illumina samtools vcf: '$illuminasamtoolsvcf
        echo 'Illumina bam: '$illuminabam
	echo 'Illumina mpileup: '$illuminampileup
	echo 'ONT mpileup: '$ontmpileup

	vcffilelist=$DIR/$prefix.filelist.txt
	echo 'File list: '$vcffilelist
	echo $nanopolishvcf > $vcffilelist
	echo $medakavcf >> $vcffilelist
	echo $samtoolsvcf >> $vcffilelist

	echo $freebayesvcf >> $vcffilelist
	echo $ivarvcf >> $vcffilelist
	echo $illuminasamtoolsvcf >> $vcffilelist

	mergedvcf=$DIR/$prefix.all_callers.combined.noallelefreqs.vcf
	finalvcf=$DIR/$prefix.all_callers.combined.vcf
	# Run merging
	java -cp $BINDIR/VariantValidator/src MergeVariants illumina_bam=$illuminabam file_list=$vcffilelist out_file=$mergedvcf

	java -cp $BINDIR/VariantValidator/src AddAlleleFrequencies vcf_file=$mergedvcf illumina_mpileup=$illuminampileup  ont_mpileup=$ontmpileup out_file=$finalvcf

    fi
done

