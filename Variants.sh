#!/usr/bin/env bash

#Author: Alli Gombolay
#This script identifies variants among cases and controls

#Usage statement
function usage () {
        echo "Usage: Variants.sh [options]
		-s Sample name(YS486, CM1, CM2, etc.)
		-d Local user directory (e.g., /projects/home/agombolay3/data/repository)"
}

#Command-line options
while getopts "s:d:h" opt; do
    case $opt in
	#Allow multiple input arguments
	s ) samples=($OPTARG) ;;
	#Allow only one input argument
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
snpSift=/projects/home/agombolay3/data/bin/snpEff/SnpSift.jar

#Reference
reference=$directory/Variant-Calling/References/sgdModified.fa
dictionary=$directory/Variant-Calling/References/sgdModified.dict

#Create FASTA dictionary file
#java -jar $picard CreateSequenceDictionary R=$reference O=$dictionary

#Determine coordinates
#for sample in ${samples[@]}; do

	#Input file
#	mapped=$directory/Variant-Calling/Alignment/$sample.bam

	#Add read groups to alignment file
#  	java -jar $picard AddOrReplaceReadGroups I=$mapped O=$sample-AddRG.bam \
#	RGLB=$sample-library RGPL=Illumina RGPU=HiSeq RGSM=$sample-sample 

#	samtools sort $sample-AddRG.bam -o $sample-AddRGSorted.bam; samtools index $sample-AddRGSorted.bam
	
  	#Mark duplicates (account for PCR duplicates)
 # 	java -jar $picard MarkDuplicates I=$sample-AddRGSorted.bam O=$sample-MarkDups.bam M=$sample.metrics.txt

#	samtools sort $sample-MarkDups.bam -o $sample-MarkDupsSorted.bam; samtools index $sample-MarkDupsSorted.bam

	#Call variants with GATK's HaplotypeCaller tool
#	java -jar $gatk -I $sample-MarkDupsSorted.bam -ERC GVCF -o $sample-raw.g.vcf -T HaplotypeCaller -R $reference

#done
  
#Joint genotyping with GATK's GenotypeGVCFs tool
#java -jar $gatk -T GenotypeGVCFs --variant YS486.g.vcf --variant CM3.g.vcf --variant CM6.g.vcf --variant CM9.g.vcf \
#--variant CM10.g.vcf --variant CM11.g.vcf --variant CM12.g.vcf --variant CM41.g.vcf -R $reference -o Variants.vcf

java -jar $gatk -T GenotypeGVCFs --variant YS486-raw.g.vcf --variant CM3-raw.g.vcf -R $reference -o Variants.vcf

#Filter variants in VCF file by quality score with SnpSift
#cat Variants.vcf | java -jar $snpSift filter "((QUAL >= 30) && (DP >= 25))" > Variants-Filtered.vcf

#java -jar $snpSift extractFields Variants-Filtered.vcf "CHROM" "POS" "REF" "ALT" "GEN[*].GT" > out.vcf

#awk -F'\t' '$6!=$7 {print $0}' out.vcf

#Create tab-delimited file of variants
#java -jar $gatk -R $reference -T VariantsToTable -V Variants-Filtered.vcf -F CHROM -F POS -F QUAL -o Variants.table
     
#Remove temporary files
#rm -f YS486-raw.g.vcf CM3-raw.g.vcf Variants.vcf
