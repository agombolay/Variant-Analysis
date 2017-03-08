#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: March 7, 2017
#This script uses Trimmomatic program to remove Illumina adapter sequences and trim reads based on quality
#The goal of this pre-processing step is to increase mapping percentage and also decrease incorrect mappings

#Usage statement
function usage () {
        echo "Usage: call-variants.sh [-s] 'sample' [-abcd] 'FASTQ.GZ' [-p] '/path/to/Trimmomatic/' [-h]
	      -a Path to input R1 FASTQ.GZ files, lane 1 ('/path/to/R1.fastq.gz')
	      -b Path to input R1 FASTQ.GZ files, lane 2 ('/path/to/R1.fastq.gz')
	      -c Path to input R2 FASTQ.GZ files, lane 1 ('/path/to/R2.fastq.gz')
	      -d Path to input R2 FASTQ.GZ files, lane 2 ('/path/to/R2.fastq.gz')
	      -s Sample name of sequenced library; used to name output files
	      -p '/projects/home/agombolay3/data/bin/Trimmomatic-0.36'"
}

#Command-line options
while getopts "a:b:c:d:s:p:h" opt;
do
    case $opt in
	a ) inputForward1=$OPTARG ;;
	b ) inputForward2=$OPTARG ;;
        c ) inputReverse1=$OPTARG ;;
	d ) inputReverse2=$OPTARG ;;
	s ) sample=$OPTARG ;;
        p ) path=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Concatenate FASTQ files from lanes 1 and 2
#cat $inputForward1 $inputForward2 > $sample-R1.fq.gz
#cat $inputReverse1 $inputReverse2 > $sample-R2.fq.gz

#STEP 1
#Trim FASTQ files based on quality and Illumina adapter content
#java -jar $path/trimmomatic-0.36.jar PE -phred33 $sample-R1.fq.gz $sample-R2.fq.gz \
#$sample-R1Paired.fq.gz $sample-R1Unpaired.fq.gz $sample-R2Paired.fq.gz $sample-R2Unpaired.fq.gz \
#ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

#Unzip files
#gunzip $sample-R1Paired.fq.gz $sample-R1Unpaired.fq.gz $sample-R2Paired.fq.gz $sample-R2Unpaired.fq.gz

#STEP 2
#Align trimmed read pairs to reference genome of interest
#bowtie2 -x sacCer2 -1 $sample-R1Paired.fq -2 $sample-R2Paired.fq \
#-U $sample-R1Unpaired.fq,$sample-R2Unpaired.fq -S $sample.sam

#SAM to BAM
samtools view -b -S $sample.sam > $sample.bam

#Sort BAM file
samtools sort -o $sample-sorted.bam $sample.bam

#Create index file
samtools index $sample-sorted.bam
