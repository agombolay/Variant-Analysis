#!/usr/bin/env bash

#Author: Alli Gombolay
#This script aligns controls to reference genome

#Usage statement
function usage () {
        echo "Usage: Alignment.sh [-s] 'Sample(s)' [-abcd] 'FASTQ.GZ' [-p] 'Path' [-i] 'Index' [-d] 'Directory' [-h]
		-s Sample name(s) (e.g., FS1, FS2, FS3)
		-w Input Read 1, lane 1 FASTQ.GZ filename
	      	-x Input Read 1, lane 2 FASTQ.GZ filename
	      	-y Input Read 2, lane 1 FASTQ.GZ filename
	      	-z Input Read 2, lane 2 FASTQ.GZ filename
	      	-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
	      	-i Basename of Bowtie2 index (e.g., sacCer2, pombe, ecoli, mm9, or hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:w:x:y:z:p:i:d:h" opt; do
    case $opt in
	s ) sample=$OPTARG ;;
	w ) forward1=$OPTARG ;;
	x ) forward2=$OPTARG ;;
        y ) reverse1=$OPTARG ;;
	z ) reverse2=$OPTARG ;;
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

#############################################################################################################################
#Input files
input1=$directory/Variant-Calling/Sequencing/$forward1
input2=$directory/Variant-Calling/Sequencing/$forward2
input3=$directory/Variant-Calling/Sequencing/$reverse1
input4=$directory/Variant-Calling/Sequencing/$reverse2

#Directory
mkdir -p $directory/Variant-Calling/Alignment

#Output files
statistics=$directory/Variant-Calling/Alignment/Bowtie2.log
mapped=$directory/Variant-Calling/Alignment/$sample-MappedReads.bam

#############################################################################################################################
#STEP 1: Concatenate FASTQ files from lanes 1 and 2
cat $input1 $input2 > R1.fq.gz; cat $input3 $input4 > R2.fq.gz

#############################################################################################################################
#STEP 2: Trim FASTQ files based on quality and Illumina adapter content
java -jar $path/trimmomatic-0.36.jar PE -phred33 R1.fq.gz R2.fq.gz R1Paired.fq.gz R1Unpaired.fq.gz \
R2Paired.fq.gz R2Unpaired.fq.gz ILLUMINACLIP:$path/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75

#Unzip files
gunzip R1Paired.fq.gz R2Paired.fq.gz

#############################################################################################################################
#STEP 3: Align pairs of reads to reference genome and save Bowtie2 log file
bowtie2 -x $index -1 R1Paired.fq -2 R2Paired.fq 2> $statistics -S temp.sam

#############################################################################################################################
#STEP 4: Extract mapped reads, convert SAM file to BAM, and sort/index BAM file
samtools view -bSf3 -F256 temp.sam | samtools sort - -o $mapped; samtools index $mapped

#Remove temporary files
rm -f R1.fq.gz R2.fq.gz R1Paired.fq R1Unpaired.fq.gz R2Paired.fq R2Unpaired.fq.gz temp.sam
