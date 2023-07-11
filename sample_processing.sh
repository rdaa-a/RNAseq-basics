#!/bin/bash


## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$


## Downloading sample file

SAMPLEDIR=$1
i=$2
NUMBERSAMPLES=$3
PIPER=$4

cd $SAMPLEDIR 

## Sample quality control and read mapping to reference genome

fastqc sample${i}.fastaq.gz
hisat2 --dta -x ../../genome/index -U sample${i}.fastaq.gz -S sample${i}.sam


## Generting sorted bam file
samtools sort -o sample${i}.bam sample${i}.sam
rm sample${i}.sam

samtools index sample${i}.bam


## Transcript assembly
stringtie -G ../../annotation/annotation.gtf -o sample${i}.gtf -l sample${i} sample${i}.bam

## Preparing merge list file for transcriptome merging
echo ${SAMPLEDIR}/sample${i}.gtf >> ../../results/merge_list.txt

## Gene Expression Quantification
stringtie -e -B -G ../../annotation/annotation.gtf -o sample${i}.gtf sample${i}.bam
rm sample${i}.bam sample${i}.bam.bai

## Blackboard

echo "sample${i} is finished" >> ../../results/blackboard.txt

## Read Blackboard

NSAMPLE=$( wc -l ../../results/blackboard.txt | awk '{ print $1 }' )

if [ $NUMBERSAMPLES -eq $NSAMPLE ]
then
	cd ../../logs

	sbatch ${PIPER}/transcriptome_merging.sh ${SAMPLEDIR}

fi




