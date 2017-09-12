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
	#STEP 1: Trim based on quality and adapters
	trim_galore --paired --length 50 $read1 $read2 -o $output
	
	#STEP 2: Align pairs of reads to reference genome and save log file
	bowtie2 -x $index -1 $output/$sample*_val_1.fq -2 $output/$sample*_val_2.fq \
	--no-mixed --no-discordant 2> $output/$sample-Bowtie2.log -S $output/temporary.sam

	#STEP 3: Extract mapped reads, convert SAM file to BAM, and sort BAM file
	samtools view -bS -f3 -F260 $output/temporary.sam | samtools sort - -o $output/$sample.bam
	
	#Index BAM file
	samtools index $output/$sample.bam

	#Remove temporary files
	rm -f $output/$sample*.fq $output/temporary.sam

done
