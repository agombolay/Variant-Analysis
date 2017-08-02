#!/usr/bin/env bash

#Author: Alli Gombolay
#This script annotates variants in VCF file

java -Xmx4g -jar snpEff.jar Saccharomyces_cerevisiae Variants.vcf > Variants-Annotated.vcf
