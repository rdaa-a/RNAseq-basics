#! /bin/bash

## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$

## Help message

if [ $# -ne 1 ]
then
	echo ""
	echo "usage: pipernas <param-file> "
	echo ""
	echo "param-file: file where the parametres are specified "
	echo ""
	exit

fi

## Param reading

PFILE=$1

echo "loading parameter file"

WD=$( grep "working-directory:" $PFILE | awk '{ print $2 }' )

echo "working directory:"$WD

FD=$( grep "folder-name" $PFILE | awk '{ print $2 }' )
GNM=$( grep "genome:" $PFILE | awk '{ print $2 }' )
ANN=$( grep "annotation:" $PFILE | awk '{ print $2 }' )
PIPER=$( grep "scripts-folder:" $PFILE | awk '{ print $2 }' )

echo ""
echo "Parametres have been loaded"

cd $WD
mkdir $FD
cd $FD
mkdir genome annotation samples results

cd genome
cp  $GNM genome.fa.gz
gunzip genome.fa.gz

cd ../annotation/
cp $ANN annotation.gtf.gz
gunzip annotation.gtf.gz

## Reference index building


extract_splice_sites.py  annotation.gtf > annot_splice.ss
extract_exons.py annotation.gtf > annot_exons.exon

cd ../genome

hisat2-build --ss ../annotation/annot_splice.ss --exon ../annotation/annot_exons.exon genome.fa index

## Generating samples folder

cd $WD

NUMSAM=$( grep "number-samples:" $PFILE | awk '{ print $2 }' )

echo "number of samples: ${NUMSAM}"

cd ${FD}/samples

IND=1

while [ $IND -le $NUMSAM ]
do
	mkdir sample${IND}
	cd sample${IND}
	
	echo ""
	echo "copying sample${IND}"

	ROUTE=$( grep "sample${IND}:" $PFILE | awk '{ print $2 }' )
	cp $ROUTE sample${IND}.fastaq.gz
	
	IND=$(($IND + 1))
	cd ../
done

cd ..
mkdir logs
cd logs

j=1

while [ $j -le $NUMSAM ]
do
	sbatch $PIPER/sample_processing.sh ${WD}${FD}/samples/sample${j} $j $NUMSAM $PIPER
	j=$(( $j + 1 ))
done
