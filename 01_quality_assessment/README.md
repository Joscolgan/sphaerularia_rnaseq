### TJ Colgan 2019

The present folder contains scripts involved in the quality assessment of Illumina sequences.  

The methods of quality assessment included:  
1. Quality assessment of raw Illumina sequences using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
2. Alignent of raw sequences against the [bumblebee reference genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000214255.1) using the RNA-Seq aligner, [STAR](https://github.com/alexdobin/STAR).  
3. Resultant alignment BAM files were indexed using [samtools](https://github.com/samtools/samtools) and assessed using [Qualimap](http://qualimap.bioinfo.cipf.es/).  

For both methods, visualisation of summary statistics for quality was performed using [MultiQC](https://multiqc.info/).  
