#usr/bin/env bash

./cutadapt.sh
R
source("barcodekeyfilegeneration.R")
# ./demultiplex_files.sh
# ./combine_cleanup_files.sh
./cutadapt2.sh
./trim.sh
mv cutadapt2/*trimmed* trimmedfiles/
cd trimmedfiles
for sample in `ls *trimmed.fq`; do mv $sample ${sample%.fastq*}.fq; pigz ${sample%.fastq*}.fq; done
cd ..
./bwa.sh
./GATKpreprocessing.sh
./GATKvariantcalltumouronly.sh
cd GATKtumouronly/
mv mpileups/ ..
cd VCFs/cleanVCFs/
mkdir ../../../totalVCFs
mv *totals*vcf ../../../totalVCFs/

mv SNPs ../../../
mv indels/ ../../../
cd ../../../SNPs
mkdir controls
mv *NEG* controls/
cd ..
./removecols.sh
mkdir mpileups/varsonly
mv mpileups/* mpileups/varsonly
cd BAMfiles/polishedBAMs
samtools merge ../../bigbam.bam `ls *bam`
cd ../../
./groupvars.sh
./SNPinfo.sh
./depth.sh
./depthsum.R
./SigAnalysis.R
