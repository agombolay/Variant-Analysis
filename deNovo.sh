#Assemble reference with de novo assembler
path=/projects/home/agombolay3/data/repository/Variant-Calling-Project/Sequencing-Results/

/projects/home/agombolay3/data/bin/SPAdes-3.10.1-Linux/bin/spades.py \
--pe1-1 $path/YS486-1_S49_L001_R1_001.fastq.gz --pe1-2 $path/YS486-1_S49_L001_R2_001.fastq.gz \
--pe1-1 $path/YS486-1_S49_L002_R1_001.fastq.gz --pe1-2 $path/YS486-1_S49_L002_R2_001.fastq.gz \
-o YS486-1-Spades-Output
