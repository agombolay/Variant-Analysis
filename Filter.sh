#!/usr/bin/env bash

#Author: Alli Gombolay
#This script annotates and filters variants in VCF file

SnpEff=/projects/home/agombolay3/data/bin/snpEff/SnpEff.jar
SnpSift=/projects/home/agombolay3/data/bin/snpEff/SnpSift.jar

java -Xmx4g -jar $SnpEff Saccharomyces_cerevisiae Variants.vcf > Variants-Annotated.vcf

#Filter by quality and depth, extract fields, and remove variants where GT of control and cases are same
cat Variants-Annotated.vcf | java -jar $SnpSift filter "((QUAL >= 30) && (DP >= 25))" | \

java -jar $SnpSift extractFields - "CHROM" "POS" "REF" "ALT" "GEN[*].GT" | \

awk -F'\t' '$12!=$5 || $12!=$6 || $12!=$7 || $12!=$8 || $12!=$9 || $12!=$10 || $12!=$11 {print $0}' - > Variants-Filtered.tab
