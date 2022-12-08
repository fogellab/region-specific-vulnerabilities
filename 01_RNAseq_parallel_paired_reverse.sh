#!/bin/sh

#for i in SRR6819407 SRR6819410; do echo $i; done | xargs -I{} -n 1 -P 10 prefetch {}; /home/dtommasini/scripts/RNAseq_parallel_paired_reverse.sh 

for accession in SRR*
do
	#Fastq
	fasterq-dump --split-files $accession;
	rm -r $accession
	
	#Fastqc
        fastqc ${accession}_1.fastq
        fastqc ${accession}_2.fastq

	#STAR
	mkdir ALIGNMENT
	/opt/STAR/bin/Linux_x86_64/STAR --runThreadN 12 --genomeDir /saveddata/shared/Reference_Genomes/STAR_indices_overhang75_mm10 --readFilesIn ${accession}_1.fastq ${accession}_2.fastq  --sjdbGTFfile /saveddata/shared/Reference_Genomes/STAR_indices_overhang75_mm10/genes.gtf --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 10 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFilterMismatchNmax 3 --twopassMode Basic --outFileNamePrefix ALIGNMENT/${accession} --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15 --quantMode TranscriptomeSAM	
	rm *fastq
	rm ALIGNMENT/*Aligned.toTranscriptome.out.bam
done

#exit here for just bam files
#exit 1

#HTSEQ
mkdir HTSEQ/
mv ALIGNMENT/*Aligned.sortedByCoord.out.bam ./

#parallelize this slow step
ls *Aligned.sortedByCoord.out.bam | xargs -I{} -n 1 -P 12 samtools index {} {}.bai
ls *Aligned.sortedByCoord.out.bam | xargs -I{} -n 1 -P 12 sh -c "~/.local/bin/htseq-count -f bam -m intersection-strict -r pos -i gene_name -s reverse {} /saveddata/shared/Reference_Genomes/STAR_indices_overhang75_mm10/genes.gtf > {}.txt"
	
for file in `ls *.txt`
do
	newname=`basename $file | sed -e "s/.txt/.LastLinesRem.txt/"`
	head -n -5 $file > "$newname"
	echo $file
	echo $newname
done
	
#move htseq results to HTSEQ folder
#move back mapped .bam files back to ALIGNMENT folder
mkdir HTSEQ/LastLinesRem/
mv *.LastLinesRem.txt HTSEQ/LastLinesRem/
mv *.txt HTSEQ/
rm *Aligned.sortedByCoord.out.bam	
mkdir FASTQC
mv *.zip FASTQC
mv *.html FASTQC/

#remove files to save space
rm *.bai

