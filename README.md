# RNAseq-basics
Scripts for preliminary analysis of RNA-seq data.
Parameters_file.txt is required for pipernas.sh running parameters, and pheno_data.csv is required for downstream_analysis.R functions.
PID2022-136633OA

#### Author
Rodrigo Bedera García

## 1. Summary
## 2. Installation
## 3. Dependencies
## 4. Input and usage
## 5. Parameters

## 1. Summary

The pipeline consists of 3 scripts that run subsequently, namely pipernas.sh, sample_processing.sh and transcriptome_merging.sh. Pipernas.sh prepares the workspace and builds the index. When done, runs sample_processing.sh for each replicate. Sample_processing.sh performs quality control, sequence alignment and gene expression quantification. Transcriptome_merging.sh merges the transcriptome of each sample and compares it to the annotated transcriptome. The output is prepared for downstream analysis with ballgown. Downstream_analysis.R file shows an example of a downstream analysis, for which the pheno_data.csv file is required

## 2. Installation

To install this software, follow these steps:

1. Download scripts folder and parameters_file.txt (located at example)
2. Place the scripts folder and parameters_file.txt wherever is desired
3. In the parameteres_file.txt, write the path to the "scripts" folder you downloaded, in path_to_scripts: . Final dash ("/") must not be written at the end.
4. Installation done!

## 3. Dependencies

Tools needed to run the bash scripts: fastqc (https://howtoinstall.co/en/fastqc), samtools (https://howtoinstall.co/en/samtools), hisat2 (https://ubuntu.pkgs.org/20.04/ubuntu-universe-arm64/hisat2_2.1.0-4_arm64.deb.html)and stringtie (https://ccb.jhu.edu/software/stringtie/index.shtml).

Packages needed to run the R script: ballgown, FactoMineR, factoextra, limma, DESeq2, biomaRt, httr, jsonlite, xml2, pheatmap, RColorBrewer, dplyr, ggplot2, ggpubrand plyr. These packages can be installed from Bioconductor or CRAN.

## 4. Input and usage

To run the bash scripts, set the parameters_file.txt accordingly to your environment, and run pipernas.sh:

pipernas.sh /full/path/to/parameters_file.txt

The R script needs to be modified manually and is presented only as a demonstration of downstream analysis.

## 5. Parameters

working-directory: The main directory where the working folder will be created

folder-name: The directory name of the workspace that will be created

genome: Path to genome file (.fa)

annotation: Path to annotation file (.gtf)

number-samples: Number of samples in the analysis

sampleX: Number of sample, and path to the fq.gz file. Add or remove lines to match your sample number

scripts-folder: Path to the folder where the bash scripts are located.

