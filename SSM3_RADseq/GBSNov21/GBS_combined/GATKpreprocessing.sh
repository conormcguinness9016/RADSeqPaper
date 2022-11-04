#/!/bin/bash

echo is this GBS data? Type Y or N
#read GBS_DATA
GBS_DATA=Y
if [ $GBS_DATA = Y ]
then
  echo Duplicate reads will not be marked!
fi
if [ $GBS_DATA = N ]
then
  echo Duplicate reads will be marked.
fi
PATH2REF=/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa
PATH2BAM=BAMfiles
echo bampath is $PATH2BAM
#PATH2FQ=${2?Error: no fastq file provided}
KNOWNSITES=/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/00-All.vcf
mkdir /Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/tmp
TEMP_DIR=/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/tmp
#mkdir $PATH2FQ/BAMfiles
#
mkdir $PATH2BAM/polishedBAMs
mkdir $PATH2BAM/RGBAMS
for sample in `ls $PATH2BAM/*.bam`
do
  x=$sample
  y=${x%%.bam*}
  z=${y##*/}
  echo x is $x y is $y z is $z
  BAMfile=$PATH2BAM/$z.bam
  if [ -f $BAMfile ]
  then
    echo sorted $z.bam and index already exists
  fi
  RG_BAM=$PATH2BAM/RGBAMS/RG_added_$z.bam
  if [ -f $RG_BAM ]
  then
    echo read groups already added for $z
  else
    echo "Adding readgroups for $z" >> errorlog.txt
    gatk AddOrReplaceReadGroups \
    -I $BAMfile \
    -O $PATH2BAM/RGBAMS/RG_added_$z.bam \
    -RGID 4 \
    -RGLB lib1 \
    -RGPL illumina \
    -RGPU unit1 \
    -RGSM $z
    if [ $? = 0 ]
    then
      echo Successfully added readgroups for $z
    else
      echo Failed to add RG for $z
    fi >> errorlog.txt
  fi
  if [ $GBS_DATA = N ]
  then
    MARKED_DUPS=$PATH2BAM/markeddups_$z.bam
    if [ ! -f $MARKED_DUPS ]
    then
    echo marking duplicates in $z as this is not GBS data
    gatk MarkDuplicates -I $RG_BAM -O $MARKED_DUPS -M $PATH2BAM/marked_dup_metrics_$z.txt --TMP_DIR $TEMP_DIR
    cp $MARKED_DUPS $PATH2BAM/preBQSR_$z.bam
    fi
  else
    echo this is GBS data, so duplicates will not be marked.
  fi
  if [ $GBS_DATA = Y ]
  then
    cp $RG_BAM $PATH2BAM/preBQSR_$z.bam
  fi
  RECAL_DATA=$PATH2BAM/recal_data_$z.table
  echo $RECAL_DATA
  if [ ! -f $RECAL_DATA ]
  then
    echo "Running BaseRecalibrator for $z" >> errorlog.txt
    gatk BaseRecalibrator -R $PATH2REF -I $PATH2BAM/preBQSR_$z.bam \
    --known-sites $KNOWNSITES \
    -O $RECAL_DATA
    if [ $? = 0 ]
    then
      echo Successfully ran baserecalibrator for $z
    else
      echo Failed to ran baserecalibrator for $z
    fi >> errorlog.txt
  fi
  POSTBQSR=$PATH2BAM/polishedBAMs/postBQSR_$z.bam
  if [ ! -f $POSTBQSR ]
  then
  echo "Running ApplyBQSR for $z" >> errorlog.txt
  gatk ApplyBQSR --bqsr-recal-file $RECAL_DATA -I $PATH2BAM/preBQSR_$z.bam -O $POSTBQSR
  if [ $? = 0 ]
  then
    echo Successfully ran ApplyBQSR for $z
  else
    echo Failed to ran ApplyBQSR for $z
  fi >> errorlog.txt
  fi
 done
ANREADYBAMS=$PATH2BAM/polishedBAMs
ls $ANREADYBAMS
./checkbams.sh $ANREADYBAMS $PATH2REF
