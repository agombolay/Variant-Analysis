#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: March 7, 2017
#This script uses Trimmomatic program to remove Illumina adapter sequences and trim reads based on quality
#The goal of this pre-processing step is to increase mapping percentage and also decrease incorrect mappings

#Usage statement
function usage () {
        echo "Usage: call-variants.sh [-s] 'sample' [-f] 'FASTQ.GZ' [-r] '.FASTQ.GZ' [-t] '/path/to/Trimmomatic/' [-h]
              -s Sample name of sequenced library; used to name output files
	      -f1 Path to input forward FASTQ.GZ files, lane 1 ('/path/to/forward.fastq.gz')
	      -f2 Path to input forward FASTQ.GZ files, lane 2 ('/path/to/forward.fastq.gz')
	      -r1 Path to input reverse FASTQ.GZ files, lane 1 ('/path/to/reverse.fastq.gz')
	      -r2 Path to input reverse FASTQ.GZ files, lane 2 ('/path/to/reverse.fastq.gz')
	      -t Path to Trimmomatic ('/projects/home/agombolay3/data/bin/Trimmomatic-0.36')"
}

#Command-line options
while getopts "f1:f2:r1:r2:t:h" opt;
do
    case $opt in
        f1 ) inputForward1=$OPTARG ;;
	f2 ) inputForward2=$OPTARG ;;
        r1 ) inputReverse1=$OPTARG ;;
	r2 ) inputReverse2=$OPTARG ;;
        t ) path=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#STEP 1
#Trim FASTQ files based on quality and Illumina adapter content (adapter read-through)
java -jar $path/trimmomatic-0.36.jar PE -phred33 $inputForward $inputReverse $sample-ForwardPaired.fastq.gz \
$sample-ForwardUnpaired.fastq.gz $sample-ReversePaired.fastq.gz $sample-ReverseUnpaired.fastq.gz \
ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75

#Unzip files
gunzip $sample-ForwardPaired.fastq.gz $sample-ForwardUnpaired.fastq.gz \
$sample-ReversePaired.fastq.gz $sample-ReverseUnpaired.fastq.gz

#STEP 2
#Align trimmed read pairs to reference genome of interest
bowtie2 -x $index -1 $sample-ForwardPaired.fastq -2 $sample-ReversePaired.fastq \
-U $sample-ForwardUnpaired.fastq $sample-ReverseUnpaired.fastq -S $sample.sam

#SAM to BAM
samtools view -b -S $sample.sam > $sample.bam

#Sort BAM file
samtools sort -o $sample-sorted.bam $sample.bam

#Create index file
samtools index $sample-sorted.bam
