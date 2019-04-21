#!/bin/env bash
########################################################################
##
## Author: Joe Colgan                     Name: run_star.sh
##
## Purpose:
## This script takes pairs of compressed fastq files and performs 
## alignment against a database of indexed transcripts using the 
## RNA-Seq aligner, STAR. The output sam file is converted in bam file
## format and sorted by read coordinates.
## The final output is a bam file
##
########################################################################

for name in ./input/*R1*fastq.gz;
do
new_name="$(echo "$name" | cut -d '/' -f 3 | cut -d '_' -f 1,2,3,4,5 )"; 
## Prin to console the name of the sample being processed:
echo "$new_name";
./src/STAR/bin/Linux_x86_64_static/STAR  \
--genomeDir ./star_indexes/ \
--runThreadN 10 \
--readFilesCommand zcat \
--readFilesIn "$name" ./input/"$new_name"R2.*fastq \
--outFileNamePrefix results/"$new_name";
## Convert sam to bam and sort by read name:
../src/samtools/samtools view -bS "$name" | ../src/samtools/samtools sort -o "$new_name".bam;
done
