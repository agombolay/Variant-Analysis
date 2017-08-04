#!/usr/bin/env bash

#Author: Alli Gombolay
#This script aligns controls to reference genome

#Usage statement
function usage () {
        echo "Usage: Align.sh [options]
		-s Sample name(YS486 CM3 CM6 CM9 CM10 CM11 CM12 CM41)
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

	#Create Directory
	mkdir -p $directory/Variant-Calling/Alignment

	#Input files
	trimmomatic=$path/trimmomatic-0.36.jar
	adapters=$path/adapters/NexteraPE-PE.fa
	read1=$directory/Variant-Calling/Sequencing/$sample-R1.fastq
	read2=$directory/Variant-Calling/Sequencing/$sample-R2.fastq
	
	#Output files
	mapped=$directory/Variant-Calling/Alignment/$sample.bam
	statistics=$directory/Variant-Calling/Alignment/$sample-Bowtie2.log

	#STEP 1: Trim FASTQ files based on quality and Illumina adapter content
	java -jar $trimmomatic PE -phred33 $read1 $read2 R1Paired.fq R1Unpaired.fq
	R2Paired.fq R2Unpaired.fq ILLUMINACLIP:$adapters:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75

	#STEP 2: Align pairs of reads to reference genome and save Bowtie2 log file
	bowtie2 -x $index -1 Paired1-Output.fq -2 Paired2-Output.fq 2> $statistics -S temp.sam

	#STEP 3: Extract mapped reads, convert SAM file to BAM, and sort/index BAM file
	samtools view -bSf3 -F256 temp.sam | samtools sort - -o $mapped; samtools index $mapped

	#Remove temporary files
	rm -f R1Paired.fq R1Unpaired.fq R2Paired.fq R2Unpaired.fq temp.sam

done
