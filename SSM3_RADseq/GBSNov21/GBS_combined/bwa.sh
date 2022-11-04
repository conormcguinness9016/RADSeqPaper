#!/usr/bin/env bash
#Script to align reads to index
PATH2REF=/Volumes/archive/cancergeneticslab/ConorM_scratch_Nov21/ref_genomes/mouse2/mouse.fa
PATH2FQ=./trimmedfiles
mkdir BAMfiles
for sample in `ls $PATH2FQ/*.fq.gz`
do
  x=$sample
  y=${x%.*}
  z=${y##*/}
  p=${z%.*}
  echo x is $x
  echo y is $y
  echo z is $z
  echo p is $p
  echo processing $z with BWA
  bwa mem -t 4 $PATH2REF $x | samtools sort -o BAMfiles/$p.bam
  BAM=BAMfiles/$p.bam
  samtools index $BAM
  echo bam file generated for $z
  echo the BAM is called $BAM
done
#echo “SAM files created, moved to folder SAMfiles within $PATH2FQ.”
