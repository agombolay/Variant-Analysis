#Alignment Results

##__YS486-1__:  
99.24% overall alignment rate  
Number of mapped reads: 11,946,179  
Number of unmapped reads: 91,457

```
samtools view -Sb -f 4 YS486-1.sam > YS486-1-unmapped.bam
samtools bam2fq YS486-1-unmapped.bam > YS486-1-unmapped.fastq
```

