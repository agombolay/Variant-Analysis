#Alignment Results

###Time:
About 1 hour per sample

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

##__Controls__:

###YS486-1:  
99.23% overall alignment rate  
Number of mapped reads: 11,919,403  
Number of unmapped reads: 92,433

###YS486-2:  
99.26% overall alignment rate  
Number of mapped reads: 11,449,813  
Number of unmapped reads: 84,981

##__Cases__:
###CM 3:

###CM 6:

###CM 9:

###CM 10:

###CM 11:

###CM 12:

###CM 41:
