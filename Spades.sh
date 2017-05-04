#!/bin/bash

#Run SPADES

#File paths
path=/projects/home/agombolay3/data/repository/Variant-Calling-Project/Sequencing
output=/projects/home/agombolay3/data/repository/Variant-Calling-Project/Assemblies

#Merge files of Library 1
zcat $path/YS486-1_S49_L001_R1_001.fastq.gz $path/YS486-1_S49_L002_R1_001.fastq.gz > $path/YS486-1_R1.fastq.gz
zcat $path/YS486-1_S49_L001_R2_001.fastq.gz $path/YS486-1_S49_L002_R2_001.fastq.gz > $path/YS486-1_R2.fastq.gz

#Merge files of Library 2
zcat $path/YS486-2_S58_L001_R1_001.fastq.gz $path/YS486-2_S58_L002_R1_001.fastq.gz > $path/YS486-2_R1.fastq.gz
zcat $path/YS486-2_S58_L001_R2_001.fastq.gz $path/YS486-2_S58_L002_R2_001.fastq.gz > $path/YS486-2_R2.fastq.gz

#Unzip files
gunzip $path/YS486-1_R1.fastq.gz $path/YS486-1_R2.fastq.gz $path/YS486-2_R1.fastq.gz $path/YS486-2_R2.fastq.gz

#Run Spades
spades.py --pe<1>-1 $path/YS486-1_R1.fastq.gz --pe<1>-2 $path/YS486-1_R2.fastq.gz --pe<2>-1 $path/YS486-2_R1.fastq.gz --pe<2>-2 $path/YS486-2_R2.fastq.gz --careful -o $output
