#!/usr/bin/env bash

#Author: Alli Gombolay
#This script identifies variants among cases and controls

#Usage statement
function usage () {
        echo "Usage: Variants.sh [options]
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "d:h" opt; do
    case $opt in
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

#Input file
mapped=$directory/Variant-Calling/Alignment/$sample.bam

#Reference
reference=$directory/Variant-Calling/References/sgdModified.fa

#Create FASTA dictionary file
java -jar $picard CreateSequenceDictionary R=$reference O=sgdModified.dict

samples=("YS486" "CM3" "CM6" "CM9" "CM10" "CM11" "CM12" "CM41")

#Determine coordinates
for sample in ${samples[@]}; do

	#Add read groups to alignment file
  	java -jar $picard AddOrReplaceReadGroups I=$sample.bam O=temp1.bam \
	RGLB=$sample RGPL=Illumina RGPU=HiSeq RGSM=$sample 

  	#Mark duplicates (account for PCR duplicates)
  	java -jar $picard MarkDuplicates I=temp1.bam O=temp2.bam M=$sample.metrics.txt

	#Call variants with GATK's HaplotypeCaller tool
	java -jar $gatk -I temp2.bam -ERC GVCF -o $sample.g.vcf -T HaplotypeCaller -R $reference

done
  
#Joint genotyping with GATK's GenotypeGVCFs tool
java -jar $gatk -T GenotypeGVCFs --variant YS486.g.vcf --variant CM3.g.vcf --variant CM6.g.vcf --variant CM9.g.vcf \
--variant CM10.g.vcf --variant CM11.g.vcf --variant CM12.g.vcf --variant CM41.g.vcf -R $reference -o Variants.vcf

#Filter variants in VCF file by quality score with SnpSift
cat Variants.vcf | java -jar $bin/snpEff/SnpSift.jar filter "( QUAL >= 30 )" > Variants-Filtered.vcf

#Remove temporary files
temp1.bam temp2.bam Variants.vcf
