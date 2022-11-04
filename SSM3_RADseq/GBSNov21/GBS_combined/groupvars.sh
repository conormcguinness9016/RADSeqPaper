#usr/bin/env bash
cd SNPs/tmp
rm -r ../SNPsbygroup
mkdir ../SNPsbygroup
for i in `ls A*gz`
do
  bcftools isec -w1 -n ~1000000 $i `ls E*gz` | bcftools view -i "INFO/DP>20" > ../SNPsbygroup/${i%%SNPS.vcf.gz*}.vcf
done

for i in `ls E*gz`
do
  bcftools isec -w1 -n ~100000 $i `ls A*gz` | bcftools view -i "INFO/DP>20" > ../SNPsbygroup/${i%%SNPS.vcf.gz*}.vcf
done
cd ../../