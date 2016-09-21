#!usr/bin/bash

DIR=`dirname $1`
FILE=`basename $1`

Pheno_Name=`echo $FILE | cut -f1 -d.`

rm -f ${DIR}/*${Pheno_Name}_tmp*;

nSNP=`cat $1 | wc -l`
step=$(($nSNP/(`nproc`-1)))

k=0
sed -n $(($k*$step+1)),$(($k*$step+$step))p $1 > ${DIR}/${Pheno_Name}_tmp${k}.sync

while ((`cat ${DIR}/${Pheno_Name}_tmp${k}.sync | wc -l` > 0)); do
	k=$(($k+1))
	sed -n $(($k*$step+1)),$(($k*$step+$step))p $1 > ${DIR}/${Pheno_Name}_tmp${k}.sync
done

echo "GWAlpha will be split into $k parallel jobs of $step SNPs maximum" 

rm $DIR/${Pheno_Name}_tmp${k}.sync

parallel --gnu -j `nproc` 'GWAlpha.py {1} {2} {3} {4} {5}' ::: $(ls ${DIR}/${Pheno_Name}_tmp*) ::: ${2:-ML} ::: ${3:-""} ::: ${4:-""} ::: ${5:-""}

printf Chromosome,Position,Allele,Alpha,MAF\\n > ${DIR}/GWAlpha_${Pheno_Name}_out.csv

for FILE in ${DIR}/GWAlpha_${Pheno_Name}_tmp*; do
	sed '1d' $FILE >> ${DIR}/GWAlpha_${Pheno_Name}_out.csv
done

grep Chromosome, ${DIR}/GWAlpha_${Pheno_Name}_out.csv > HD
grep -v Chromosome, ${DIR}/GWAlpha_${Pheno_Name}_out.csv | sort -t\, -k 1,1 -k 2,2n -o ${DIR}/GWAlpha_${Pheno_Name}_out.csv
cat HD ${DIR}/GWAlpha_${Pheno_Name}_out.csv > tmp && mv tmp ${DIR}/GWAlpha_${Pheno_Name}_out.csv

rm HD ${DIR}/*${Pheno_Name}_tmp*

GWAlphaPlot.r ${DIR}/GWAlpha_${Pheno_Name}_out.csv

