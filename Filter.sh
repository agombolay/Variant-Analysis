#Filter variants in VCF file by quality score
cat Variants.vcf | java -jar $snpSift filter "((QUAL >= 30) && (DP >= 25))" > Variants-Filtered.vcf

java -jar $snpSift extractFields Variants-Filtered.vcf "CHROM" "POS" "REF" "ALT" "GEN[*].GT" > out.vcf

awk -F'\t' '$6!=$7 {print $0}' out.vcf

awk -F'\t' '($6!= "./.") {print $0}' out2.vcf
