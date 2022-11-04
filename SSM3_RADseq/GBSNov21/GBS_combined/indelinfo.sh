#usr/bin/env bash

cd indels
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
  rm ../../commonindels*
  bcftools isec -n +6 -w1 `ls *vcf.gz` > ../../commonindelsbcftools.vcf
  bcftools isec -n +6 `ls *vcf.gz` > ../../commonindelsbcftoolsSITESonly.bed
  bgzip ../../commonindelsbcftools.vcf
  tabix ../../commonindelsbcftools.vcf.gz
  
  
  cd ..
  cd ..
  rm -r indels/uniqueindels
  mkdir indels/uniqueindels
  for i in `ls indels/tmp/*.gz`
  do
  x=$i
  y=${x%%indels.vcf.gz*}
  z=${y##*/}
  echo $z
    bcftools isec -w2 -n ~01 commonindelsbcftools.vcf.gz $i | bcftools view -i "INFO/DP>20" > indels/uniqueindels/$z.vcf
  done

rm uniqueindelsnoheader.txt
rm Allindelsnoheader.txt
rm indelsbygroupnoheader.txt
cd indels/uniqueindels
for i in `ls`
do
  grep -v ^# $i | awk -v myvar=${i%%.vcf*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../uniqueindelsnoheader.txt
done
cd ../../
cd indels/tmp
for i in `ls *gz`
do
  bcftools view -i "INFO/DP>20" $i | grep -v ^# | awk -v myvar=${i%%indels.vcf.gz*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../Allindelsnoheader.txt
done

# cd ../indelsbygroup
# 
# for i in `ls`
# do
#   grep -v ^# $i | awk -v myvar=${i%%.vcf*} '{print $1, $2, $4, $5, myvar, $10 }' >> ../../indelsbygroupnoheader.txt
# done

cd ../../

    