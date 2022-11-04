#!/usr/bin/env bash
PATH2BAMS=${1?Provide BAM file folder}
PATH2REF=${2?Error: no reference file}
mkdir $PATH2BAMS/BAMQC

for i in `ls $PATH2BAMS`
do
  x=$i
  echo x is $x
  y=${x%%.bam}
  echo y is $y
  z=${y##*/}
  echo z is $z
  gatk DepthOfCoverage \
  -R $PATH2REF \
  -o $PATH2BAMS/BAMQC/$z \
  -I $PATH2BAMS/$z \
  -ct 30 \
  -ct 50 \
  -ct 100
done
