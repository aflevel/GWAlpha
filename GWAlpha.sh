#!usr/bin/bash

Pheno_Name=`echo $1 | cut -f1 -d.`

rm -f *${Pheno_Name}_tmp*;

nSNP=`cat $1 | wc -l`
step=$(($nSNP/(`nproc`-1)))

k=0
sed -n $(($k*$step+1)),$(($k*$step+$step))p $1 > ${Pheno_Name}_tmp${k}

while ((`cat ${Pheno_Name}_tmp${k} | wc -l` > 0)); do
        k=$(($k+1))
        sed -n $(($k*$step+1)),$(($k*$step+$step))p $1 > ${Pheno_Name}_tmp${k}
done

echo "GWAlpha will be split into $k parallel jobs of $step SNPs maximum" 

rm ${Pheno_Name}_tmp${k}

parallel --gnu -j `nproc` 'python GWAlpha.py {} ML' ::: $(ls ${Pheno_Name}_tmp*)

printf Chromosome,Position,Allele,Alpha,MAF\\n > GWAlpha_${Pheno_Name}_out.csv

for FILE in GWAlpha_${Pheno_Name}_tmp*; do
	sed '1d' $FILE >> GWAlpha_${Pheno_Name}_out.csv
done

sort -t\, -k 1,1n -k 2,2n GWAlpha_${Pheno_Name}_out.csv -o GWAlpha_${Pheno_Name}_out.csv

rm *${Pheno_Name}_tmp*

Rscript GWAlphaPlot.r ${Pheno_Name}

