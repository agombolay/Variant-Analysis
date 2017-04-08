# Alignment Results

### Time:
About 1 hour per sample

### Count # of mapped reads
```
samtools view -c -F 4 sample.bam
```

### Count # of unmapped reads
```
samtools view -c -f 4 sample.bam
```

### Examine unmapped reads
```
samtools view -b -f 4 sample.bam > sample-unmapped.bam
samtools bam2fq sample-unmapped.bam > sample-unmapped.fastq
```

## __Controls__:

### YS486-1:  
99.23% overall alignment rate  
Number of mapped reads: 11,919,403  
Number of unmapped reads: 92,433

About 0.6% of positions have 0 coverage
```
bedtools genomecov -d -ibam YS486-1.bam -g sacCer2.bed > YS486-1-Coverage.bed
grep -w 0$ YS486-1-Coverage.bed | wc -l
wc -l  YS486-2-Coverage.bed
```

### YS486-2:  
99.26% overall alignment rate  
Number of mapped reads: 11,449,813  
Number of unmapped reads: 84,981

About 0.6% of positions have 0 coverage
```
bedtools genomecov -d -ibam YS486-2.bam -g sacCer2.bed > YS486-2-Coverage.bed
grep -w 0$ YS486-2-Coverage.bed | wc -l
wc -l  YS486-2-Coverage.bed
```

## __Cases__:
### CM 3:
99.27% overall alignment rate

### CM 6:
99.09% overall alignment rate

### CM 9:
99.13% overall alignment rate

### CM 10:
99.17% overall alignment rate

### CM 11:
99.26% overall alignment rate

### CM 12:

### CM 41:
