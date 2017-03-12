#Alignment Results

###Count # of mapped reads
```
samtools view -c -F 4 sample.bam
```

###Count # of unmapped reads
```
samtools view -c -f 4 sample.bam
```

###Examine unmapped reads
```
samtools view -b -f 4 sample.bam > sample-unmapped.bam
samtools bam2fq sample-unmapped.bam > sample-unmapped.fastq
```

##__YS486-1__:  
99.23% overall alignment rate  
Number of mapped reads: 11,919,403  
Number of unmapped reads: 92,433

##__YS486-2__:  
99.26% overall alignment rate  
Number of mapped reads: 11,449,813  
Number of unmapped reads: 84,981
