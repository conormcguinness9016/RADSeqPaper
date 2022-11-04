
#/!/bin/bash

PATH2REF=../../ref_genomes/mouse2/mouse.fa
#KNOWNSITES=${3?Error: need the known sites vcf}
PATH2BAM=test
ANREADYBAMS=$PATH2BAM/polishedBAMs
echo 1st argument is $1
echo 2nd argumnt is $2
echo path2BAM is $PATH2BAM

./GATKpreprocessing.sh
./GATKvariantcalltumouronly.sh $ANREADYBAMS $PATH2REF
