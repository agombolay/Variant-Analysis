#!/usr/bin/env bash

#Author: Alli Gombolay
#This script aligns cases and controls to reference genome

#Usage statement
function usage () {
        echo "Usage: Align.sh [options]
		-s Sample names (YS486 CM3 CM6 CM9 CM10 CM11 CM12 CM41)
	      	-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
	      	-i Basename of Bowtie2 index (e.g., sacCer3, pombe, ecoli, mm9, or hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:p:i:d:h" opt; do
    case $opt in
	s ) samples=$OPTARG ;;
        p ) path=$OPTARG ;;
	i ) index=$OPTARG ;;
	d ) directory=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Determine coordinates
for sample in ${samples[@]}; do

	#Input files
	trimmomatic=$path/trimmomatic-0.36.jar
	adapters=$path/adapters/NexteraPE-PE.fa
	bed=$directory/Variant-Calling/References/sacCer3.bed
	read1=$directory/Variant-Calling/FASTQ-Files/$sample-R1.fastq
	read2=$directory/Variant-Calling/FASTQ-Files/$sample-R2.fastq
	
	#Output directory
	output=$directory/Variant-Calling/Alignment
	
	#Create Directory
	mkdir -p $output

#############################################################################################################################
	#STEP 1: Trim FASTQ files based on quality and Illumina adapter content
	java -jar $trimmomatic PE -phred33 $read1 $read2 $output/Paired1.fq $output/Unpaired1.fq \
	$output/Paired2.fq $output/Unpaired2.fq ILLUMINACLIP:$adapters:2:30:10 MINLEN:75
		
	#STEP 2: Align pairs of reads to reference genome and save Bowtie2 log file
	bowtie2 -x $index -1 $output/Paired1.fq -2 $output/Paired2.fq --no-mixed --no-discordant \
	2> $output/$sample-Bowtie2.log -S $output/temporary.sam

	#STEP 3: Extract mapped reads, convert SAM file to BAM, and sort BAM file
	samtools view -bS -f3 -F260 $output/temporary.sam | samtools sort - -o $output/$sample.bam
	
	#Index BAM file
	samtools index $output/$sample.bam

	#STEP 4: Calculate genome coverage
	bedtools genomecov -d -ibam $output/$sample.bam -g $bed > $output/$sample.bed
	
	#Remove temporary files
	rm -f $output/Paired*.fq $output/Unpaired*.fq $output/temporary.sam

done
