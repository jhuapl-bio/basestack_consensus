#!/bin/bash
if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

# specify a sequencing run directory
RUN=$1

# get known directory names
DIR="/home/idies/workspace/covid19/sequencing_runs/$RUN/artic-pipeline/5-post-filter"
REF="/home/idies/workspace/covid19/ncov_reference/sequence.fasta"
GENES="/home/idies/workspace/covid19/ncov_reference/genes.gff3"
for i in `ls $DIR/*variant_data.txt`
do
    sample=${i##*/}
    samplename=${sample%%.*}
    echo $i
    echo $samplename

    consensusvcf=$DIR/$samplename.consensus.vcf
    allvcf=$DIR/$samplename.allsnps.vcf
    java -cp $BINDIR/VariantValidator/src TableToVcf table_file=$i consensus_file=$consensusvcf all_file=$allvcf

    consensuscombinedvcf=$DIR/$samplename.consensus.combined.vcf
    allcombinedvcf=$DIR/$samplename.allsnps.combined.vcf
    java -cp $BINDIR/VariantValidator/src CombineVariants vcf_file=$consensusvcf out_file=$consensuscombinedvcf genome_file=$REF gene_file=$GENES
    java -cp $BINDIR/VariantValidator/src CombineVariants vcf_file=$allvcf out_file=$allcombinedvcf genome_file=$REF gene_file=$GENES
done


