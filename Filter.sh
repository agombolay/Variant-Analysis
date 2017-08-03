#!/usr/bin/env bash

#Author: Alli Gombolay
#This script annotates and filters variants in VCF file

SnpEff=/projects/home/agombolay3/data/bin/snpEff/snpEff.jar
SnpSift=/projects/home/agombolay3/data/bin/snpEff/SnpSift.jar

java -Xmx4g -jar $SnpEff Saccharomyces_cerevisiae Variants.vcf > Variants-Annotated.vcf

#Filter by quality and depth and remove variants where GT of control and cases are same
cat Variants-Annotated.vcf | java -jar $SnpSift filter "((QUAL >= 30) && (DP >= 25))" > Variants.tab

java -jar $SnpSift extractFields Variants.tab "CHROM" "POS" "REF" "ALT" "GEN[*].GT" "ANN[*].EFFECT" > Variants2.tab

#awk -F'\t' '$12!=$5 || $12!=$6 || $12!=$7 || $12!=$8 || $12!=$9 || $12!=$10 || $12!=$11 {print $0}' > Variants.tab
