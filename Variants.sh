#!/usr/bin/env bash

#Author: Alli Gombolay
#This script identifies variants among cases and controls

#Usage statement
function usage () {
        echo "Usage: Variants.sh [options]
		-s Sample names (YS486, CM1, CM2, etc.)
		-d /projects/home/agombolay3/data/repository"
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

#Path to java files
picard=/projects/home/agombolay3/data/bin/picard.jar
gatk=/projects/home/agombolay3/data/bin/GenomeAnalysisTK.jar
snpSift=/projects/home/agombolay3/data/bin/snpEff/SnpSift.jar

#Reference input files
reference=$directory/Variant-Calling/References/sacCer3.fa
dictionary=$directory/Variant-Calling/References/sacCer3.dict

#Create FASTA index
samtools faidx $reference

#Create FASTA dictionary file
java -jar $picard CreateSequenceDictionary R=$reference O=$dictionary

#Determine coordinates
for sample in ${samples[@]}; do

	#Input file
	VCF=$directory/Variant-Calling/References/sacCer3.vcf
	mapped=$directory/Variant-Calling/Alignment/$sample.bam

	#Add read groups
  	java -jar $picard AddOrReplaceReadGroups I=$mapped O=$sample-RG.bam \
	RGLB=$sample-library RGPL=Illumina RGPU=HiSeq RGSM=$sample-sample 

	#Sort and index BAM
	samtools sort $sample-RG.bam -o $sample-RGSort.bam; samtools index $sample-RGSort.bam
	
  	#Mark duplicate reads
  	java -jar $picard MarkDuplicates I=$sample-RGSort.bam O=$sample-DeDup.bam M=$sample.txt

	#Sort and index BAM
	samtools sort $sample-DeDup.bam -o $sample-DeDupSort.bam; samtools index $sample-DeDupSort.bam
	
	#Base quality score recalibration
	java -jar $gatk -T BaseRecalibrator -R $reference -I $sample-DeDupSort.bam -knownSites $VCF -o recal.grp
   
   	#Create a recalibrated BAM with print reads
   	#java -jar $gatk -T PrintReads -R $reference -I $sample-DeDupSort.bam -BQSR recal.grp -o $sample-Recalibrated.bam
   
	#Call variants with HaplotypeCaller (ploidy=1)
	#java -jar $gatk -T HaplotypeCaller -R $reference -I $sample-Recalibrated.bam -ERC GVCF -o $sample.g.vcf -ploidy 1

done
  
#Joint genotyping with GATK's GenotypeGVCFs tool
#java -jar $gatk -T GenotypeGVCFs --variant YS486.g.vcf --variant CM3.g.vcf --variant CM6.g.vcf --variant CM9.g.vcf \
#--variant CM10.g.vcf --variant CM11.g.vcf --variant CM12.g.vcf --variant CM41.g.vcf -R $reference -o Variants.vcf

#Remove temporary files
#rm -f *-RG.bam *-RGSort.bam *-DeDup.bam *-DeDupSort.bam *-Recalibrated.bam
