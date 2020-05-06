if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

RUN=$1
DIR=/home/idies/workspace/covid19/illumina/$RUN/trimmed
REF=/home/idies/workspace/covid19/ncov_reference/sequence.fasta
OUTDIR=/home/idies/workspace/covid19/illumina/$RUN/trimmedvariantcalling

javac $BINDIR/CoverageNormalization/src/*.java

if [ ! -d $OUTDIR ]
then
    mkdir $OUTDIR
fi

#bwa index $REF

for i in `ls $DIR/*.bam`
do
    # Get prefix without directory/file extensions
    echo $i
    prefix=`echo $i`
    prefix=${prefix##*/}
    prefix=${prefix%%.*}
    echo 'Prefix: '$prefix

    # Reheader the bam with SM added into the read group lines so freebayes won't complain
    bamreheader=$OUTDIR/$prefix.reheader.bam
    echo 'Bam file with fixed header: '$bamreheader
    if [ ! -r $bamreheader ]
    then
        samtools view -H $i | sed 's/@RG/@RG\tSM:sample/g' | samtools reheader - $i > $bamreheader
    fi

    # Convert reheadered bam to sam to pass to normalization
    samfull=$OUTDIR/$prefix.unnormalized.sam
    echo 'Unnormalized sam file: '$samfull
    if [ ! -r $samfull ]
    then
        samtools view -h $bamreheader > $samfull
    fi

    # Normalize sam file
    samfile=$OUTDIR/$prefix.sam
    echo 'Normalized SAM file: '$samfile
    if [ ! -r $samfile ]
    then
        java -cp $BINDIR/CoverageNormalization/src NormalizePairedReads input=$samfull output=$samfile coverage_threshold=200 --qual_sort 
    fi

    bamfileunsorted=$OUTDIR/$prefix.norm200.unsorted.bam
    echo 'Normalized unsorted BAM file: '$bamfileunsorted
    if [ ! -r $bamfileunsorted ]
    then
        samtools view -h -b $samfile > $bamfileunsorted
    fi

    bamfile=$OUTDIR/$prefix.norm200.bam
    echo 'Normalized BAM file: '$bamfile
    if [ ! -r $bamfile ]
    then
        samtools sort $bamfileunsorted > $bamfile
    fi
    
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

    ivartsv=$OUTDIR/$prefix.ivar.tsv
    echo 'ivar tsv: '$ivartsv
    if [ ! -r $ivartsv ]
    then
        echo 'Calling ivar variants'
        samtools mpileup -A -d 0 --reference $REF -B -Q 0 $bamfile | ivar variants -p $OUTDIR/$prefix.ivar -t 0.15
    fi

    ivarvcf=$OUTDIR/$prefix.ivar.vcf
    echo 'ivar vcf: '$ivarvcf
    if [ ! -r $ivarvcf ]
    then
        echo 'Converting ivar tsv to vcf'
	java -cp $BINDIR/VariantValidator/src IvarToVcf table_file=$ivartsv out_file=$ivarvcf
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

    longreadmpileup=`ls $longreaddir/*/artic-pipeline/4-draft-consensus/$longreadsample*.mpileup`
    echo 'Long read mpileup: '$longreadmpileup

    samplewithbarcode=$medakavcf
    samplewithbarcode=${samplewithbarcode##*/}
    samplewithbarcode=${samplewithbarcode%%.*}
    echo 'Sample with barcode: '$samplewithbarcode

    # Create symlinks for easier lookups
    fblink=$OUTDIR/$samplewithbarcode.freebayes.vcf
    fbulink=$OUTDIR/$samplewithbarcode.freebayes.unfiltered.vcf
    slink=$OUTDIR/$samplewithbarcode.samtools_illumina.vcf
    ivlink=$OUTDIR/$samplewithbarcode.ivar.vcf

    if [ -r $fblink ]
    then
        rm $fblink
    fi

    if [ -r $fbulink ]
    then
        rm $fbulink
    fi

    if [ -r $slink ]
    then
        rm $slink
    fi
    if [ -r $ivlink ]
    then
        rm $ivlink
    fi
    
    ln -s $freebayesvcf $fblink
    ln -s $freebayesunfiltered $fbulink
    ln -s $samtoolsvcf $slink
    ln -s $ivarvcf $ivlink

    vcffilelist=$OUTDIR/$samplewithbarcode.filelist.txt
    echo 'Vcf file list: '$vcffilelist
    if [ ! -r $vcffilelist ]
    then
        echo 'Making VCF file list'
	echo $nanopolishvcf > $vcffilelist
	echo $medakavcf >> $vcffilelist
        echo $longsamtoolsvcf >> $vcffilelist
	echo $freebayesvcf >> $vcffilelist
	echo $ivarvcf >> $vcffilelist
	echo $samtoolsvcf >> $vcffilelist
    fi

    allcallersvcf=$OUTDIR/$samplewithbarcode.all_callers.combined.vcf
    echo 'Merged vcf: '$allcallersvcf
    if [ ! -r $allcallersvcf ]
    then
        echo 'Merging calls'
	java -cp $BINDIR/VariantValidator/src MergeVariants file_list=$vcffilelist out_file=$allcallersvcf illumina_bam=$bamfile
    fi

    mergedallelefreqvcf=$OUTDIR/$samplewithbarcode.combined_allele_freqs.vcf
    echo 'VCF with combined allele frequencies: '$mergedallelefreqvcf
    java -cp $BINDIR/VariantValidator/src AddAlleleFrequencies vcf_file=$allcallersvcf illumina_mpileup=$mpileup ont_mpileup=$longreadmpileup out_file=$mergedallelefreqvcf

    if [ -r $freebayesvcf ] && [ -r $ivarvcf ] && [ -r $samtoolsvcf ] && [ -r $nanopolishvcf ] && [ -r $medakavcf ] && [ -r $longsamtoolsvcf ] && [ -r $allcallersvcf ]
    then
        newpath=${medakavcf::-18}.all_callers.combined.vcf
	echo 'Copying combined vcf to: '$newpath
        cp $mergedallelefreqvcf $newpath
    fi

done
