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
reference=$directory/Variant-Calling/References/sacCer2.fa
dictionary=$directory/Variant-Calling/References/sacCer2.dict

#Create FASTA index
#samtools faidx sacCer2.fa

#Create FASTA dictionary file
#java -jar $picard CreateSequenceDictionary R=$reference O=$dictionary

#Determine coordinates
for sample in ${samples[@]}; do

	#Input file
	mapped=$directory/Variant-Calling/Alignment/$sample.bam

	#Add read groups
  	java -jar $picard AddOrReplaceReadGroups I=$mapped O=$sample-AddRG.bam \
	RGLB=$sample-library RGPL=Illumina RGPU=HiSeq RGSM=$sample-sample 

	samtools sort $sample-AddRG.bam -o $sample-AddRGSort.bam; samtools index $sample-AddRGSort.bam
	
  	#Mark duplicates
  	java -jar $picard MarkDuplicates I=$sample-AddRGSort.bam O=$sample-MarkDups.bam M=$sample.metrics.txt

	samtools sort $sample-MarkDups.bam -o $sample-MarkDupsSort.bam; samtools index $sample-MarkDupsSort.bam
	
	#Base quality score recalibration
	java -jar $gatk -T BaseRecalibrator -R $reference -I $sample-MarkDupsSort.bam -knownSites sacCer3.vcf -o table
   
   	java -jar GenomeAnalysisTK.jar -T PrintReads -R reference.fasta -I input.bam -BQSR recalibration_report.grp -o output.bam
   
	#Call variants with GATK's HaplotypeCaller tool
	java -jar $gatk -T HaplotypeCaller -R $reference -I $sample-MarkDupsSort.bam -ERC GVCF -o $sample.g.vcf -ploidy 1

done
  
#Joint genotyping with GATK's GenotypeGVCFs tool
java -jar $gatk -T GenotypeGVCFs --variant YS486.g.vcf --variant CM3.g.vcf --variant CM6.g.vcf --variant CM9.g.vcf \
--variant CM10.g.vcf --variant CM11.g.vcf --variant CM12.g.vcf --variant CM41.g.vcf -R $reference -o Variants.vcf

#Remove temporary files
