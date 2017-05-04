# Genome Assembly
1. Spades, Abyss, and Velvet
2. CISA
3. [Quast](http://quast.sourceforge.net/download.html)
4. GeneMark

Helpful tools:
* Scaffolder: [Bambus](http://amos.sourceforge.net/wiki/index.php/Bambus)  
* Genome annotation pipeline: [Maker](http://www.yandell-lab.org/software/maker.html)  
* Visualize de novo assembly graphs: [Bandage](https://rrwick.github.io/Bandage/)  
* Biology and NGS data analysis tools: [Geneious](http://www.geneious.com/)  
* Align large genomes: [LASTZ (Geneious plugin)](http://www.geneious.com/plugins/lastz-plugin)  

# Variant-Analysis

1. Trim reads: Trimmomatic

2. Alignment to reference: Bowtie2

3. Examine unmapped reads with FastQC to determine why they did not align

4. Determine alignment coverage and identify any regions with no coverage

5. Add read groups and mark duplicates: Picard Tools

6. Variant calling and joint genotyping: GATK Tools

7. Filter VCF file: SnpSift

## Install Spades:
1. Download the software
```
wget http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz
```
2. Unpack the software
```
tar -xzf SPAdes-3.10.1-Linux.tar.gz
```
3. Move to directory
```
cd SPAdes-3.10.1-Linux/bin/
```

## Install Abyss:
1. Create my own abyss folder in which to save executables
```
mkdir abyss
```

2. Configure program and save the executables to abyss folder
```
./configure --prefix=/projects/home/agombolay3/data/bin/abyss
```

3. Build the software
```
make
```

4. Install the software
```
make install
```

## Install Velvet
1. Download the software
```
wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
```
2. Unpack the software
```
tar -zxvf velvet/velvet_1.2.10.tgz
```
3. Build the software
```
make
```

## Install CISA
1. Download the software
```
wget http://sb.nhri.org.tw/CISA/upload/en/2014/3/CISA_20140304-05194132.tar
```
2. Unpack the software
```
tar -xvf CISA_20140304-05194132.tar
```
3. Change permissions
```
chmod 755 -R CISA1.3
```

## Run SPADES
```
path=/projects/home/agombolay3/data/repository/Variant-Calling-Project/Sequencing
output=/projects/home/agombolay3/data/repository/Variant-Calling-Project/Assemblies
```

Merge files of Library 1
```
zcat $path/YS486-1_S49_L001_R1_001.fastq.gz $path/YS486-1_S49_L002_R1_001.fastq.gz > $path/YS486-1_R1.fastq.gz
zcat $path/YS486-1_S49_L001_R2_001.fastq.gz $path/YS486-1_S49_L002_R2_001.fastq.gz > $path/YS486-1_R2.fastq.gz
```
Merge files of Library 2
```
zcat $path/YS486-2_S58_L001_R1_001.fastq.gz $path/YS486-2_S58_L002_R1_001.fastq.gz > $path/YS486-2_R1.fastq.gz
zcat $path/YS486-2_S58_L001_R2_001.fastq.gz $path/YS486-2_S58_L002_R2_001.fastq.gz > $path/YS486-2_R2.fastq.gz
```
Run Spades
```
spades.py --pe<1>-1 $path/YS486-1_R1.fastq.gz --pe<1>-2 $path/YS486-1_R2.fastq.gz --pe<2>-1 $path/YS486-2_R1.fastq.gz --pe<2>-2 $path/YS486-2_R2.fastq.gz --careful -o $output
```
