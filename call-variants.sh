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
while getopts "a:b:c:d:s:p:l:u:n:h" opt;
do
    case $opt in
	a ) inputForward1=$OPTARG ;;
	b ) inputForward2=$OPTARG ;;
        c ) inputReverse1=$OPTARG ;;
	d ) inputReverse2=$OPTARG ;;
	s ) sample=$OPTARG ;;
        p ) path=$OPTARG ;;
	l ) RGLB=$OPTARG ;;
	u ) RGPU=$OPTARG ;;
	n ) RGSM=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Remove old version of output
rm $sample-R1.fq.gz $sample-R2.fq.gz $sample-R1Paired.fq \
$sample-R2Paired.fq $sample.sam $sample.bam $sample.bam.bai

#Concatenate FASTQ files from lanes 1 and 2
cat $inputForward1 $inputForward2 > $sample-R1.fq.gz
cat $inputReverse1 $inputReverse2 > $sample-R2.fq.gz

#STEP 1
#Trim FASTQ files based on quality and Illumina adapter content
java -jar $path/trimmomatic-0.36.jar PE -phred33 $sample-R1.fq.gz $sample-R2.fq.gz \
$sample-R1Paired.fq.gz $sample-R1Unpaired.fq.gz $sample-R2Paired.fq.gz $sample-R2Unpaired.fq.gz \
ILLUMINACLIP:$path/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:15 MINLEN:75

#Unzip files
gunzip $sample-R1Paired.fq.gz $sample-R2Paired.fq.gz

#STEP 2
#Align trimmed forward and reverse reads to reference genome of interest
bowtie2 -x sacCer2 -1 $sample-R1Paired.fq -2 $sample-R2Paired.fq -S $sample.sam

#Convert SAM file to BAM file format and sort BAM file
samtools view -b -S $sample.sam | samtools sort -o $sample.bam -

#Create index file
samtools index $sample.bam

#Add read groups to alignment file
java -jar picard.jar AddOrReplaceReadGroups I=$sample.bam \
O=$sample-v1.bam RGLB=$RGLB RGPL=illumina RGPU=$RGPU RGSM=$RGSM
      
#Mark duplicates (account for PCR duplicates)
#java -jar picard.jar MarkDuplicates I=$sample-v1.bam \
#O=$sample-final.bam M=$sample.duplication-metrics.txt
      
#GATK
