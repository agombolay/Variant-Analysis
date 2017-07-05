#!/usr/bin/env bash

#Author: Alli Gombolay
#This script filter variants in VCF file

#Filter by quality and depth, extract fields, and remove variants where GT of control and cases are same
cat Variants.vcf | java -jar /projects/home/agombolay3/data/bin/snpEff/SnpSift.jar filter "((QUAL >= 30) && (DP >= 25))" | \
java -jar /projects/home/agombolay3/data/bin/snpEff/SnpSift.jar extractFields - "CHROM" "POS" "REF" "ALT" "GEN[*].GT" | \
awk -F'\t' '$12!=$5 || $12!=$6 || $12!=$7 || $12!=$8 || $12!=$9 || $12!=$10 || $12!=$11 {print $0}' - > Variants.tab

