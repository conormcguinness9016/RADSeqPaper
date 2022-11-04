#usr/bin/env bash
###cut the universal adapter off each read

cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -o Adaptercut_1.fastq Illumina/SQ1724_HKMH2DRXY_s_1.fastq
cutadapt -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGCA -o Adaptercut_1_2.fastq Adaptercut_1.fastq
cutadapt -a CGAGATCGGAAGAGCGGACTTTAAGC -o Adaptercut_1_3.fastq Adaptercut_1_2.fastq
cutadapt -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -o Adaptercut_1_4.fastq Adaptercut_1_3.fastq

cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG -o Adaptercut_2.fastq Illumina/SQ1724_HKMH2DRXY_s_2.fastq
cutadapt -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGCA -o Adaptercut_2_2.fastq Adaptercut_2.fastq
cutadapt -a CGAGATCGGAAGAGCGGACTTTAAGC -o Adaptercut_2_3.fastq Adaptercut_2_2.fastq
cutadapt -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -o Adaptercut_2_4.fastq Adaptercut_2_3.fastq

mkdir cutadapt2
cutadapt --no-indels -e 0 -g file:bc.fasta -o "cutadapt2/{name}.fastq.gz" Adaptercut.fastq

echo cutadapt.sh executed successfully

