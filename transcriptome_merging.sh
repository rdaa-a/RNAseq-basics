#!/bin/bash

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0$

## Accessing results folder

SAMPLEFOLDER=$1

cd $SAMPLEFOLDER
cd ../../results/
echo ""
echo "merged transcriptome en proceso"

## Merging sample transcriptomes
stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf merge_list.txt

## Comparing our assembly with the reference
gffcompare -r ../annotation/annotation.gtf -G -o comparison stringtie_merged.gtf

echo ""
echo "merged"
