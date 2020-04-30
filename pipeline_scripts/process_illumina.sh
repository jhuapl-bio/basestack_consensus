if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

RUN=$1
DIR=/home/idies/workspace/covid19/illumina/$RUN/trimmed
REF=/home/idies/workspace/covid19/ncov_reference/sequence.fasta
OUTDIR=$DIR/../trimmedvariantcalling

if [ ! -d $OUTDIR ]
then
    mkdir $OUTDIR
fi

#bwa index $REF

for i in `ls $DIR/*.bam`
do
    echo $i
    prefix=`echo $i`
    prefix=${prefix##*/}
    prefix=${prefix%%.*}
    #prefix=`echo "${prefix/_R1_/_}"`
    echo 'Prefix: '$prefix

    #samunsorted=$OUTDIR/$prefix.sam
    #echo 'Unsorted sam file: '$samunsorted
    #if [ ! -r $samunsorted ]
    #then
    #	echo 'Aligning reads'    
    #    bwa mem $REF $i $pair > $samunsorted
    #fi

    #bamunsorted=$OUTDIR/$prefix.unsorted.bam
    #echo 'Unsorted bam file: '$bamunsorted
    #if [ ! -r $bamunsorted ]
    #then
    #    echo 'Converting alignments to BAM'
    #    samtools view -b $samunsorted > $bamunsorted	    
    #fi
    bamfile=$OUTDIR/$prefix.reheader.bam
    echo 'Bam file: '$bamfile
    if [ ! -r $bamfile ]
    then
        samtools view -H $i | sed 's/@RG/@RG\tSM:sample/g' | samtools reheader - $i > $bamfile
    fi
    bamfile=$OUTDIR/$prefix.bam
    
    if [ ! -r $bamfile.bai ]
    then
        echo 'Indexing alignments'
        samtools index $bamfile
    fi

    mpileup=$OUTDIR/$prefix.mpileup
    if [ ! -r $mpileup ]
    then
	echo 'Computing mpileup'
        samtools mpileup --reference $REF $bamfile -o $mpileup
    fi

    # Run samtools-based variant calling
    javac $BINDIR/VariantValidator/src/*.java
    samtoolsvcf=$OUTDIR/$prefix.samtools.vcf
    echo 'Samtools vcf: '$samtoolsvcf
    if [ ! -r $samtoolsvcf ]
    then
        echo 'Calling samtools variants'
        java -cp $BINDIR/VariantValidator/src CallVariants flag_prefix=ILLUMINA_ pileup_file=$mpileup out_file=$samtoolsvcf
    fi

    freebayesunfiltered=$OUTDIR/$prefix.freebayes.unfiltered.vcf
    echo 'Freebayes unfiltered vcf: '$freebayesunfiltered
    if [ ! -r $freebayesunfiltered ]
    then
        echo 'Calling freebayes variants'
        freebayes -f $REF $bamfile > $freebayesunfiltered
    fi

    freebayesvcf=$OUTDIR/$prefix.freebayes.vcf
    echo 'Freebayes vcf: '$freebayesvcf
    if [ ! -r $freebayesvcf ]
    then
        echo 'Filtering Freebayes vcf'
        vcftools --vcf $freebayesunfiltered --minQ 10 --min-meanDP 20 --recode --recode-INFO-all --out $OUTDIR/$prefix.freebayes
        mv $OUTDIR/$prefix.freebayes.recode.vcf $freebayesvcf
    fi

    shortreadsample=$prefix
    shortreadsample=${shortreadsample%_*}
    
    echo 'Short read sample: '$shortreadsample
    longreadsample=`cat $BINDIR/samplenamemap.txt | awk -v srs="$shortreadsample" '{ if ($1 == srs) { print $2; } }'`
    echo 'Long read sample: '$longreadsample

    if [ "$longreadsample" = 'NTC' ]
    then
        continue
    fi
    longreaddir='/home/idies/workspace/covid19/sequencing_runs'
    medakavcfzipped=`ls $longreaddir/*/artic-pipeline/4-draft-consensus/$longreadsample*.medaka.merged.vcf.gz`
    echo 'Medaka vcf (zipped): '$medakavcfzipped

    nanopolishvcf=`ls $longreaddir/*/artic-pipeline/4-draft-consensus/$longreadsample*.nanopolish.merged.vcf`
    echo 'Nanopolish vcf: '$nanopolishvcf

    medakavcf=${medakavcfzipped::-3}
    if [ ! -r $medakavcf ]
    then
        echo 'Unzipping '$medakavcfzipped
        gunzip -c $medakavcfzipped > $medakavcf
    fi

    longsamtoolsvcf=`ls $longreaddir/*/artic-pipeline/4-draft-consensus/$longreadsample*.samtools.vcf`
    echo 'Long-read samtoolsvcf: '$longsamtoolsvcf

    samplewithbarcode=$medakavcf
    samplewithbarcode=${samplewithbarcode##*/}
    samplewithbarcode=${samplewithbarcode%%.*}
    echo 'Sample with barcode: '$samplewithbarcode

    # Create symlinks for easier lookups
    ln -s $freebayesvcf $OUTDIR/$samplewithbarcode.freebayes.vcf
    ln -s $freebayesunfiltered $OUTDIR/$samplewithbarcode.freebayes.unfiltered.vcf
    ln -s $samtoolsvcf $OUTDIR/$samplewithbarcode.illumina.samtools.vcf

    vcffilelist=$OUTDIR/$samplewithbarcode.filelist.txt
    echo 'Vcf file list: '$vcffilelist
    if [ ! -r $vcffilelist ]
    then
        echo 'Making VCF file list'
	echo $nanopolishvcf > $vcffilelist
	echo $medakavcf >> $vcffilelist
        echo $longsamtoolsvcf >> $vcffilelist
	echo $freebayesvcf >> $vcffilelist
	echo $samtoolsvcf >> $vcffilelist
    fi

    allcallersvcf=$OUTDIR/$samplewithbarcode.all_callers.combined.vcf
    echo 'Merged vcf: '$allcallersvcf
    if [ ! -r $allcallersvcf ]
    then
        echo 'Merging calls'
	java -cp $BINDIR/VariantValidator/src MergeVariants file_list=$vcffilelist out_file=$allcallersvcf
    fi

done
