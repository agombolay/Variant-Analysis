# Genome Assembly
1. Spades, Abyss, and Velvet

2. CISA

3. [Quast](http://quast.sourceforge.net/download.html)

4. GeneMark

Scaffolder: [Bambus](http://amos.sourceforge.net/wiki/index.php/Bambus)
Align large genomes: [LASTZ](http://www.geneious.com/plugins/lastz-plugin)
NGS analysis tools: [Geneious](http://www.geneious.com/)
Genome annotation pipeline: [Maker](http://www.yandell-lab.org/software/maker.html)
Visualizing de novo assembly graphs: [Bandage](https://rrwick.github.io/Bandage/)

# Variant-Analysis

1. Trim reads: Trimmomatic

2. Alignment to reference: Bowtie2

3. Examine unmapped reads with FastQC to determine why they did not align

4. Determine alignment coverage and identify any regions with no coverage

5. Add read groups and mark duplicates: Picard Tools

6. Variant calling and joint genotyping: GATK Tools

7. Filter VCF file: SnpSift
