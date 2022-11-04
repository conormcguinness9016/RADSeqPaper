#usr/bin/env bash
rm uniqueSNPsnoheader.txt
rm AllSNPsnoheader.txt
rm SNPsbygroupnoheader.txt
cd SNPs/uniqueSNPs
for i in `ls`
do
  grep -v ^# $i | awk -v myvar=${i%%.vcf*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../uniqueSNPsnoheader.txt
done
cd ../../
cd SNPs/tmp
for i in `ls *gz`
do
  bcftools view -i "INFO/DP>20" $i | grep -v ^# | awk -v myvar=${i%%SNPS.vcf.gz*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../AllSNPsnoheader.txt
done

cd ../SNPsbygroup

for i in `ls`
do
  grep -v ^# $i | awk -v myvar=${i%%.vcf*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../SNPsbygroupnoheader.txt
done

cd ../../


