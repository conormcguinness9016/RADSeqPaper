#usr/bin/env bash

for i in `cat diffmers.txt`
do
ls split_files/*_$i.fq
done