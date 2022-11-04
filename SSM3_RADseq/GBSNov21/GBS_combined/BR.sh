#/!/bin/bash
gatk BaseRecalibrator -R ../../ref_genomes/mouse2/mouse.fa -I test/preBQSR_A1A12.bam \
--known-sites /Volumes/scratch/cancergeneticslab/ConorM/ref_genomes/mouse2/00-All.vcf \
-O recal_data_A1A12.table
