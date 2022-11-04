#!/bin/bash
#use GATK to call variants in GBS data
PATH2BAMS=BAMfiles/polishedBAMs
PATH2REF=/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa
FILTERMIN=${3:-20}
FILTERBY=${4:-20}
FILTERMAX=${5:-20}
mkdir GATKtumouronly
GATK=GATKtumouronly
mkdir $GATK/VCFs
mkdir $GATK/mpileups
mkdir $GATK/VCFs/cleanVCFs
PATH2VCFs=$GATK/VCFs
PATH2MPILEUPS=$GATK/mpileups
CLEANVCFs=$GATK/VCFs/cleanVCFs
##note gatk must be in path
for sample in `ls $PATH2BAMS/*.bam`
do
  x=$sample
  y=${x%%.bam*}
  z=${y##*/}
  echo "starting variant call $z (time: `date +"%T"`) ">> variantcalllog.txt
  if [ ! -f $GATK/mpileups/$z.mpileup ]
  then
    echo "starting mpileup for $z (time: `date +"%T"`)" >> variantcalllog.txt
    bcftools mpileup -f $PATH2REF $PATH2BAMS/$z.bam | bcftools call -mv -Ob -o $GATK/mpileups/$z.mpileup
    if [ $? = 0 ]
    then
      echo "Successfully completed mpileup for $z (time: `date +"%T"`)"
    else
      echo "Failed to complete mpileup for $z (time: `date +"%T"`)"
    fi >> variantcalllog.txt
  fi
  if [ ! -f $GATK/VCFs/$z.vcf ]
  then
      echo "starting Mutect2 for $z (time: `date +"%T"`)" >> variantcalllog.txt
      gatk Mutect2 \
      -R $PATH2REF \
      -I $PATH2BAMS/$z.bam \
      -O $GATK/VCFs/raw_$z.vcf
      if [ $? = 0 ]
      then
        echo "Successfully completed Mutect2 for $z (time: `date +"%T"`)"
      else
        "echo Failed to complete Mutect2 for $z (time: `date +"%T"`)"
      fi >> variantcalllog.txt
  fi
  if [ ! -f $PATH2VCFs/GATKfiltered_$z.vcf ]
  then
      echo "Filtering Mutect for $z (time: `date +"%T"`)" >> variantcalllog.txt
      gatk FilterMutectCalls \
      -V $GATK/VCFs/raw_$z.vcf \
      -R $PATH2REF \
      -O $CLEANVCFs/GATKfiltered_$z.vcf
      if [ $? = 0 ]
      then
        echo "Successfully filtered variants for $z (time: `date +"%T"`)"
      else
        echo "Failed to complete filtering variants for $z (time: `date +"%T"`)"
      fi >> variantcalllog.txt
  fi
  echo $z finito >> variantcalllog.txt
done
if [ ! -d $GATK/VCFs/SNPs ]
then
    sepvars.sh $CLEANVCFs
fi
if [ ! -f $PATH2VCFs/*total.vcf ]
then
for sample in `ls $PATH2VCFs/*vcf`
do
  # x=$sample
  # y=${x%%.vcf*}
  # z=${y##*/}
  # mv $PATH2VCFs/$z.vcf $PATH2VCFs/$z"total".vcf
done
fi
echo $FILTERMIN, $FILTERBY, $FILTERMAX
./multiplefilters.sh $PATH2VCFs $PATH2MPILEUPS $PATH2REF $GATK $FILTERMIN $FILTERBY $FILTERMAX
