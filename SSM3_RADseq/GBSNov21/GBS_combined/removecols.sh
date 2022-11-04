#usr/bin/env bash

cd SNPs
echo moved to $PWD
rm -r tmp
mkdir tmp
for i in `ls *vcf`
do
  x=$i
  y=${x%%.vcf*}
  z=${y##*R_}
  echo z is $z
  cut -f 1-10 $i > tmp/$z.vcf
  TMPVCF=tmp/$z.vcf
  echo tmpvcf is $TMPVCF
  bgzip $TMPVCF
  if [ $? = 0 ]
  then
    echo "Successfully bgzipped $z (time: `date +"%T"`)"
  else
    echo "Failed to bgzip for $z (time: `date +"%T"`)"
  fi
  tabix $TMPVCF.gz
  if [ $? = 0 ]
  then
    echo "Successfully tabixed $z (time: `date +"%T"`)"
  else
    echo "Failed to tabixed for $z (time: `date +"%T"`)"
  fi
done
 cd tmp
bcftools isec -n +6 -w1 `ls *vcf.gz` > ../../commonSNPsbcftools.vcf
bcftools isec -n +6 `ls *vcf.gz` > ../../commonSNPsbcftoolsSITESonly.bed
bgzip ../../commonSNPsbcftools.vcf
tabix ../../commonSNPsbcftools.vcf.gz


cd ..
cd ..
mkdir SNPs/uniqueSNPs
for i in `ls SNPs/tmp/*.gz`
do
  x=$i
  y=${x%%SNPS.vcf.gz*}
  z=${y##*/}
  echo $z
  bcftools isec -w3 -n ~001 commonSNPsbcftools.vcf.gz 0D5data/postBQSR_0D5SNPS.vcf.gz $i | bcftools view -i "INFO/DP>20" > SNPs/uniqueSNPs/$z.vcf
done

#>
