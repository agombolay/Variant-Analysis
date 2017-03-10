#Alignment Results

##__YS486-1__:  
99.24% overall alignment rate  
Number of mapped reads: 11,946,179  
Number of unmapped reads: 91,457

###Time
Trimming and alignment: 92m27.523s

###Examine unmapped reads
```
samtools view -f 4 YS486-1.bam > YS486-1-unmapped.bam
samtools bam2fq YS486-1-unmapped.bam > YS486-1-unmapped.fastq
```

##__YS486-2__:  
% overall alignment rate  
Number of mapped reads:
Number of unmapped reads:

```
samtools view -Sb -f 4 YS486-2.sam > YS486-2-unmapped.bam
samtools bam2fq YS486-2-unmapped.bam > YS486-2-unmapped.fastq
```
