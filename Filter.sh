#!/usr/bin/env bash

#Author: Alli Gombolay
#This script filter variants in VCF file

snpSift=/projects/home/agombolay3/data/bin/snpEff/SnpSift.jar

#Filter variants in VCF file by quality score
cat Variants.vcf | java -jar $snpSift filter "((QUAL >= 30) && (DP >= 25))" > Variants-Filtered1.vcf

#Extract chrom, position, ref, alt, and GT from VCF file
java -jar $snpSift extractFields Variants-Filtered1.vcf "CHROM" "POS" "REF" "ALT" "GEN[*].GT" > Variants-Filtered2.vcf

#Remove variants where GT of controls and cases are same
awk -F'\t' '$5!=$6 && $5!=$7 && $5!=$8 && $5!=$9 && $5!=$10 && $5!=$11 && $5!=$12 {print $0}' Variants-Filtered2.vcf > Variants-Filtered3.vcf

#Remove temporary files
rm -f Variants-Filtered1.vcf Variants-Filtered2.vcf Variants-Filtered3.vcf
