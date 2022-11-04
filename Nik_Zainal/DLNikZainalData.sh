#!/usr/bin/env bash
echo downloading Nik-Zainal data
wget ftp://ftp.sanger.ac.uk/pub/cancer/Nik-ZainalEtAl-560BreastGenomes/Caveman_560_20Nov14_clean.txt
echo take header off and create allmutsnohder.txt file
grep -v "##" Caveman_560_20Nov14_clean.txt > geninfo560genomes.txt
