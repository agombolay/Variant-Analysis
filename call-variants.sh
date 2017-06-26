#!/usr/bin/env bash

#Author: Alli Gombolay
#This script identifies variants among cases and controls

#Usage statement
function usage () {
        echo "Usage: Variants.sh [options]
		-s Sample name(s)
	      	-f Input Read 1 FASTQ.GZ filename
	      	-r Input Read 2 FASTQ.GZ filename
	      	-p Path (e.g., /projects/home/agombolay3/data/bin/Trimmomatic-0.36)
	      	-i Basename of Bowtie2 index (e.g., sacCer2, pombe, ecoli, mm9, or hg38)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:f:r:p:i:d:h" opt; do
    case $opt in
	  s ) sample=$OPTARG ;;
    f ) forward=$OPTARG ;;
	  r ) reverse=$OPTARG ;;
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

#Path to Picard java files
picard=/projects/home/agombolay3/data/bin/picard.jar
gatk=/projects/home/agombolay3/data/bin/GenomeAnalysisTK.jar

#Reference
reference=/projects/home/agombolay3/data/repository/Variant-Calling/References/sgdModified.fa

#Create FASTA dictionary file
java -jar $picard CreateSequenceDictionary R=$reference O=sgdModified.dict

samples=("YS486" "CM3" "CM6" "CM9" "CM10" "CM11" "CM12" "CM41")

#Determine coordinates
for sample in ${samples[@]}; do

	#Add read groups to alignment file
  	java -jar $picard AddOrReplaceReadGroups I=$sample.bam O=temp1.bam RGLB=$sample RGPL=Illumina RGPU=HiSeq RGSM=$sample 

  	#Mark duplicates (account for PCR duplicates)
  	java -jar $picard MarkDuplicates I=temp1.bam O=temp2.bam M=$sample.metrics.txt; samtools index temp2.bam

	#Call variants
	java -jar $gatk -I temp2.bam -ERC GVCF -o $sample.g.vcf -T HaplotypeCaller -R $reference

done
  
#Joint genotyping
java -jar $gatk -T GenotypeGVCFs --variant YS486.g.vcf --variant CM3.g.vcf --variant CM6.g.vcf --variant CM9.g.vcf \
--variant CM10.g.vcf --variant CM11.g.vcf --variant CM12.g.vcf --variant CM41.g.vcf -R $reference -o Variants.vcf    

#Annotate VCF file based on sacCer2 reference
#java -Xmx4g -jar /projects/home/agombolay3/data/bin/snpEff/snpEff.jar sacCer2 Variants1.vcf  > Variants1-Annotated.vcf 

#Filter variants in VCF file by quality score
#cat Variants1-Annotated.vcf | java -jar $bin/snpEff/SnpSift.jar filter "( QUAL >= 30 )" > Variants-Filtered.vcf
