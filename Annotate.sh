#!/usr/bin/env bash

#Author: Alli Gombolay
#This script annotates variants in VCF file

java -Xmx4g -jar snpEff.jar Saccharomyces_cerevisiae examples/test.chr22.vcf > test.chr22.ann.vcf
