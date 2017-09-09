# Variant Analysis

1. Trim reads: Trimmomatic
2. Alignment to reference: Bowtie2
3. Examine unmapped reads with FastQC
4. Determine alignment coverage of genome
5. Add read groups and mark duplicates: Picard Tools
6. Variant calling and joint genotyping (Halotype GVCF): GATK Tools
7. Filter VCF file based on quality, depth, and mutation type: SnpSift


# Reference Genome
* [sacCer3 (2011 Yeast Genome)](http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/)
