#!/bin/bash

basepath="/sciserver/vc_crypt/covid19/vc1"
repopath="$basepath/code/vcfigv"

reference="$basepath/code/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"

runpath="$basepath/sequencing_runs/20200405_0424_GA40000_FAN28191_4c470e7c"
pipelinepath="$runpath/artic-pipeline/3-hac-medaka-norm200"

find "$pipelinepath" -name "*.pass.vcf.gz" -print0 | while read -d $'\0' vcf; do

	gzip -d "$vcf"
	vcf="${vcf%.gz}"
	bam="${vcf%.pass.vcf}.primertrimmed.rg.sorted.bam"
	outPrefix=$(basename "${vcf%.pass.vcf}")

	if [[ -s "$vcf" && -s "$bam" ]]; then

		mkdir "$outPrefix"

		java -cp "$repopath/src" \
			Vcf2Bat \
			--squish \
			--nocombine \
			aln="$bam" \
			var="$vcf" \
			genome="$reference" \
			outprefix="$outPrefix"

		mv "$outPrefix.bat" "$outPrefix"
		mv "$outPrefix" "$pipelinepath"

		"$repopath/xvfb-run" \
			--auto-servernum "$repopath/IGV_2.3.98/igv.sh" \
			-b "$pipelinepath/$outPrefix/$outPrefix.bat"

	fi

	gzip "$vcf"

done

