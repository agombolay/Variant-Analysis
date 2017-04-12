#!/usr/bin/env bash

#Author: Alli Gombolay
#This script identifies variants among cases and controls

#Path to bin folder
bin=/projects/home/agombolay3/data/bin

#Path to reference file
reference=/projects/home/agombolay3/data/repository/Variant-Calling-Project/Variant-Calling/sacCer2.fa

#Create FASTA index file
samtools faidx sacCer2.fa

#Create FASTA dictionary file
java -jar $bin/picard.jar CreateSequenceDictionary R=$reference O=sacCer2.dict

samples=("YS486-1" "CM3" "CM6" "CM9" "CM10" "CM11" "CM12" "CM41")

#Determine coordinates
for sample in ${samples[@]}; do

#Add read groups to alignment file
java -jar $bin/picard.jar AddOrReplaceReadGroups I=$sample.bam O=$sample-AddReadGroups.bam \
RGLB=$sample RGPL=Illumina RGPU=HiSeq RGSM=$sample #LB = library, PL = platform, SM = sample

#Mark duplicates (account for PCR duplicates)
#java -jar $bin/picard.jar MarkDuplicates I=$sample-AddReadGroups.bam O=$sample-MarkDups.bam M=$sample.metrics.txt

#Create index file
#samtools index $sample-MarkDups.bam

#Call variants
#java -jar $bin/GenomeAnalysisTK.jar -I $sample-MarkDups.bam -ERC GVCF -o $sample.g.vcf -T HaplotypeCaller -R $reference

#Joint genotyping
#java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs --variant YS486-1.g.vcf \
#--variant CM3.g.vcf --variant CM6.g.vcf --variant CM9.g.vcf --variant CM10.g.vcf \
#--variant CM11.g.vcf --variant CM12.g.vcf --variant CM41.g.vcf -R $reference -o Variants1.vcf    

#Annotate VCF file based on sacCer2 reference
#java -Xmx4g -jar /projects/home/agombolay3/data/bin/snpEff/snpEff.jar sacCer2 Variants1.vcf  > Variants1-Annotated.vcf 

#Filter variants in VCF file by quality score
#cat Variants1-Annotated.vcf | java -jar $bin/snpEff/SnpSift.jar filter "( QUAL >= 30 )" > Variants1-Filtered.vcf
