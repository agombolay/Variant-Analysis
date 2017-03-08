#!/usr/bin/env bash

#Author: Alli Gombolay
#Date: March 7, 2017
#This script uses Trimmomatic program to remove Illumina adapter sequences and trim reads based on quality
#The goal of this pre-processing step is to increase mapping percentage and also decrease incorrect mappings

#Usage statement
function usage () {
        echo "Usage: call-variants.sh [-s] 'sample' [-f] 'FASTQ.GZ' [-r] '.FASTQ.GZ' [-t] '/path/to/Trimmomatic/' [-h]
              -s Sample name of sequenced library; used to name output files
	      -f Path to input forward FASTQ.GZ files ('/path/to/forward.fastq.gz')
	      -r Path to input reverse FASTQ.GZ files ('/path/to/reverse.fastq.gz')
	      -t '/projects/home/agombolay3/data/bin/Trimmomatic-0.36'"
}

#Command-line options
while getopts "f:r:t:h" opt;
do
    case $opt in
        f ) inputForward=$OPTARG ;;
        r ) inputReverse=$OPTARG ;;
        t ) path=$OPTARG ;;
        #Print usage statement
        h ) usage ;;
    esac
done

#Exit program if [-h]
if [ "$1" == "-h" ]; then
        exit
fi

#Trim FASTQ files based on quality and Illumina adapter content (adapter read-through)
java -jar $path/trimmomatic-0.36.jar PE -phred33 $inputForward $inputReverse $sample-ForwardPaired.fastq.gz \
$sample-ForwardUnpaired.fastq.gz $sample-ReversePaired.fastq.gz $sample-ReverseUnpaired.fastq.gz \
ILLUMINACLIP:$path/adapters/TruSeq3-SE.fa:2:30:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
