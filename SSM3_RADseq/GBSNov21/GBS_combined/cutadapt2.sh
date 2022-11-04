#usr/bin/env bash
mkdir cutadapt2
mkdir cutadapt2/lane1
mkdir cutadapt2/lane2
mkdir adaptertrimmed
cutadapt --no-indels -e 0 -g file:bc.fasta -o "cutadapt2/lane1/{name}.fastq.gz" Illumina/SQ1724_HHYT5DRXY_s_1.fastq
cutadapt --no-indels -e 0 -g file:bc.fasta -o "cutadapt2/lane2/{name}.fastq.gz" Illumina/SQ1724_HHYT5DRXY_s_2.fastq
for i in `cat samples.txt`
do
  cat cutadapt2/lane1/$i*fastq.gz > cutadapt2/lane1_$i.fastq.gz
  cat cutadapt2/lane2/$i*fastq.gz > cutadapt2/lane2_$i.fastq.gz
  cat cutadapt2/lane*_$i.fastq.gz > adaptertrimmed/$i.fastq.gz
  done
