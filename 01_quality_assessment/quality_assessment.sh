#!/bin/env bash
## 
## Quality assessment was performed using FastQC for raw read sequences. 
## Summary report was generated using MultiQC.

## Ensure that FastQC and MultiQC in your PATH:
## echo "$PATH"
for name in *.fastq.gz; 
do 
fastqc -t 10 "$name";
done

## Run MultiQC:
multiqc ./
