#!/bin/bash

PATH2VCFs=${1?No VCFs found}
PATH2PILEUPS=${2?No pileup files found}
PATH2REF=${3?No reference provided}
OUTPUTDIR=${4:-"."}
FILTERMIN=20
FILTERBY=20
FILTERMAX=20
i=20
mkdir $OUTPUTDIR/filtercounts/
FILTERCOUNTS=$OUTPUTDIR/filtercounts/
mkdir $PATH2VCFs/GATKmutectfiltered/
for i in `seq $FILTERMIN $FILTERBY $FILTERMAX`
do
 if [ -d $FILTERCOUNTS/$i ]
  then
    echo WARNING: $i filter folder already exists, assuming finished. Check all samples have been filtered!
  else
    echo filtering threshold: $i
    mkdir $FILTERCOUNTS/$i
    for sample in `find $PATH2VCFs -name "*.vcf"`
    do
      x=$sample
      y=${x%.*}
      z=${y##*/}
      echo found file $z, filtering by $i reads
      bcftools view -O v -f .,PASS $x | bcftools query -i "INFO/DP>$i" -f '%CHROM %POS %REF %ALT %QUAL %DP\n' >$FILTERCOUNTS/$i/$z.filtered.vcf
    done
    for sample in `ls $PATH2PILEUPS/*mpileup`
    do
      x=$sample
      y=${x%.*}
      z=${y##*/}
      echo pileups for $z
      bcftools query -i "INFO/DP>$i" -f '%CHROM %POS %REF %ALT %QUAL %DP\n' $PATH2PILEUPS/$z.mpileup > $FILTERCOUNTS/$i/$z.filtered.mpileup
    done
  fi
  if [ -e $FILTERCOUNTS/$i/variantscalledfiltered$i.txt ]
  then
    echo already counted for $i move on son
  else
    wc -l $FILTERCOUNTS/$i/*.filtered.mpileup > $FILTERCOUNTS/$i/basesmappedfiltered$i.txt
    awk -v var=$i 'BEGIN{ OFS=","; print "Bases Mapped, File, Filter threshold"}; NR > 0{print $1, $2, var;}' $FILTERCOUNTS/$i/basesmappedfiltered$i.txt > $FILTERCOUNTS/$i/basesmappedfiltered$i.csv
    wc -l $FILTERCOUNTS/$i/*.filtered.vcf > $FILTERCOUNTS/$i/variantscalledfiltered$i.txt
    awk -v var=$i 'BEGIN{ OFS=","; print "Variants called, File, Filter threshold"}; NR > 0{print $1, $2, var;}' $FILTERCOUNTS/$i/variantscalledfiltered$i.txt > $FILTERCOUNTS/$i/variantscalledfiltered$i.csv
    find $FILTERCOUNTS/$i/ -name "*.csv" -exec sed -i "s&$FILTERCOUNTS/$i/&&g" '{}' +
    find $FILTERCOUNTS/$i/ -name "*.csv" -exec sed -i "s&.filtered.vcf&&g" '{}' +
    find $FILTERCOUNTS/$i/ -name "*.csv" -exec sed -i "s&.filtered.mpileup&&g" '{}' +
    echo counted $i
  fi
done
if [ -d $OUTPUTDIR/tables/variants ]
then
  echo
else
  mkdir $OUTPUTDIR/tables
  mkdir $OUTPUTDIR/tables/coverage
  mkdir $OUTPUTDIR/tables/variants
fi
find $FILTERCOUNTS -name "bases*.csv" -exec mv {} $OUTPUTDIR/tables/coverage \;
find $FILTERCOUNTS -name "variants*.csv" -exec mv {} $OUTPUTDIR/tables/variants \;
