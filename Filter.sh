#!/usr/bin/env bash

#Author: Alli Gombolay
#This script annotates and filters variants in VCF file

SnpEff=/projects/home/agombolay3/data/bin/snpEff/snpEff.jar
SnpSift=/projects/home/agombolay3/data/bin/snpEff/SnpSift.jar

java -Xmx4g -jar $SnpEff Saccharomyces_cerevisiae Variants.vcf > Variants-Annotated.vcf

#Filter by quality and depth
cat Variants-Annotated.vcf | java -jar $SnpSift filter "((QUAL >= 30) && (DP >= 25))" > Variants1.tab

#Extract fields from VCF file
java -jar $SnpSift extractFields Variants1.tab "CHROM" "POS" "REF" "ALT" "GEN[*].GT" "ANN[*].IMPACT:" > Variants2.tab

#Remove variants where genotype of control = cases
awk -F'\t' '$12!=$5||$12!=$6||$12!=$7||$12!=$8||$12!=$9||$12!=$10||$12!=$11 {print $0}' Variants2.tab > Variants3.tab

#Remove variants that are synonymous for all samples
awk -F'\t' '$13!='synonymous_variant'||$14!='synonymous_variant'||$15!='synonymous_variant'||$16!='synonymous_variant'\
||$17!='synonymous_variant' {print $0}' Variants3.tab > Variants.tab

#Remove temporary files
rm -f Variants1.tab Variants2.tab Variants3.tab
