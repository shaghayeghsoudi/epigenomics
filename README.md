# Epigenomics
## This repository contains scripts and pipelines for processign and analyzing DNA methylation.
### Sample processing
RRBS libraries were prepared using the Premium Reduced Representation Bisulfite Sequencing (RRBS) Kit (Diagenode Cat# C02030033), according to the manufacturer’s protocol. 100 ng of genomic DNA were used to start library preparation for each sample.
Following library preparation, samples were pooled together by 6.
RRBS library pools quality control was performed by measuring DNA concentration of the pools using the Qubit. the profile of the poolswas checked using theHigh SensitivityDNA
chip for 2100 Bioanalyzer.

### Bioinformatics workflow
Raw RRBS FASTQ files were trimmed using Trimegalore RRBS mode (https://github.com/FelixKrueger/TrimGalore). 
Trimmed RRBS FASTQ files were mapped to NCBI Human Reference Genome Build GRCh37 (hg19) using Bismark pbat mode.
(https://github.com/FelixKrueger/Bismark)

DNA methylation ratio and differential methylated cytosine (DMCs/DMRs) were analyzed by using DSS (https://github.com/haowulab/DSS/tree/master) and edgeR (https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

