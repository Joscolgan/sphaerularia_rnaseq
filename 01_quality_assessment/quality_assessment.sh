#!/bin/env bash
########################################################################
##
## Author: Joe Colgan                     Name: quality_assessment.sh
##
##
## Purpose:
## Quality assessment was performed using FastQC for raw read sequences.
## Summary report was generated using MultiQC.
## This script takes compressed fastq files as input.
## The script expects the files to be in the same directory as where
## the script is run.
##
########################################################################
## Define the number of threads:
threads=4

## Ensure that FastQC and MultiQC in your PATH:
## echo "$PATH"
for name in *.fastq.gz;
do
fastqc -t "$threads" "$name";
done

## Run MultiQC:
multiqc ./
