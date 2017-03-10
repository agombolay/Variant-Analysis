#Alignment Results

###Count # of mapped reads
```
samtools view -c -F 4 YS486-1.bam
```

###Count # of unmapped reads
```
samtools view -c -f 4 YS486-2.bam
```

##__YS486-1__:  
99.24% overall alignment rate  
Number of mapped reads: 11,946,179  
Number of unmapped reads: 91,457

###Examine unmapped reads
```
samtools view -b -f 4 YS486-1.bam > YS486-1-unmapped.bam
samtools bam2fq YS486-1-unmapped.bam > YS486-1-unmapped.fastq
```

##__YS486-2__:  
99.27% overall alignment rate  
Number of mapped reads: 11,477,747  
Number of unmapped reads: 84,002  

```
samtools view -b -f 4 YS486-2.bam > YS486-2-unmapped.bam
samtools bam2fq YS486-2-unmapped.bam > YS486-2-unmapped.fastq
```
