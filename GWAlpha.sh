#!/bin/bash

DIR=`dirname $1`
FILE=`basename $1`

Pheno_Name=`echo $FILE | cut -f1 -d.`

rm -f ${DIR}/*${Pheno_Name}_tmp*;

nSNP=`cat $1 | wc -l`
step=$(($nSNP/(`nproc`-1)))

k=0
sed -n $(($k*$step+1)),$(($k*$step+$step))p $1 > ${DIR}/${Pheno_Name}_tmp${k}.sync

while ((`cat ${DIR}/${Pheno_Name}_tmp${k} | wc -l` > 0)); do
	k=$(($k+1))
	sed -n $(($k*$step+1)),$(($k*$step+$step))p $1 > ${DIR}/${Pheno_Name}_tmp${k}.sync
done

echo "GWAlpha will be split into $k parallel jobs of $step SNPs maximum" 

rm $DIR${Pheno_Name}_tmp${k}

parallel --gnu -j `nproc` 'GWAlpha.py {1} {2} {3} {4} {5}' ::: $(ls ${DIR}/${Pheno_Name}_tmp*) ::: ${2:-ML} ::: ${3:-""} ::: ${4:-""} ::: ${5:-""}

printf Chromosome,Position,Allele,Alpha,MAF\\n > ${DIR}/GWAlpha_${Pheno_Name}_out.csv

for FILE in ${DIR}/GWAlpha_${Pheno_Name}_tmp*; do
	sed '1d' $FILE >> ${DIR}/GWAlpha_${Pheno_Name}_out.csv
done

sort -t\, -k 1,1n -k 2,2n ${DIR}/GWAlpha_${Pheno_Name}_out.csv -o ${DIR}/GWAlpha_${Pheno_Name}_out.csv

rm ${DIR}/*${Pheno_Name}_tmp*

GWAlphaPlot.r ${DIR}/GWAlpha_${Pheno_Name}_out.csv pval

