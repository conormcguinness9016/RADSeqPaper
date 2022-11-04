#usr/bin/env bash
mkdir split_files

##nonamers
PATH2BARCODES=/Volumes/scratch/cancergeneticslab/ConorM/Oct2020_SSM3Data/MyAnalysis/barcodes/nonabarcodes.txt
PATH2FASTQ=/Volumes/scratch/cancergeneticslab/ConorM/Oct2020_SSM3Data/MyAnalysis/Adaptercut.fastq

cat $PATH2FASTQ | fastx_barcode_splitter.pl --bcfile $PATH2BARCODES --prefix split_files/nona_ --suffix .fq --bol

#octamers
PATH2BARCODES=/Volumes/scratch/cancergeneticslab/ConorM/Oct2020_SSM3Data/MyAnalysis/barcodes/octabarcodes.txt
PATH2FASTQ=/Volumes/scratch/cancergeneticslab/ConorM/Oct2020_SSM3Data/MyAnalysis/split_files/nona_unmatched.fq
##make new directory split_files
## demultiplexes fastq file (fastq file must be present
cat $PATH2FASTQ | fastx_barcode_splitter.pl --bcfile $PATH2BARCODES --prefix split_files/octa_ --suffix .fq --bol


###heptamers
PATH2BARCODES=/Volumes/scratch/cancergeneticslab/ConorM/Oct2020_SSM3Data/MyAnalysis/barcodes/heptabarcodes.txt
PATH2FASTQ=/Volumes/scratch/cancergeneticslab/ConorM/Oct2020_SSM3Data/MyAnalysis/split_files/octa_unmatched.fq
##make new directory split_files
## demultiplexes fastq file (fastq file must be present
cat $PATH2FASTQ | fastx_barcode_splitter.pl --bcfile $PATH2BARCODES --prefix split_files/hepta_ --suffix .fq --bol

echo demultiplex_files.sh executed successfully


