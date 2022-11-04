#usr/bin/env bash
mkdir depth
rm -r depthlong
mkdir depthlong
rm depthlong.txt
for i in `ls BAMfiles/polishedBAMs/*bam`
do
  BAM=$i
  x=$i
  y=${x%%.bam*}
  z=${y##*R_}
  p=${z#*.bam}
  echo $y
  echo $z
  echo $p
  # samtools depth $BAM  > depth/$z.txt
  # BP20=$(awk '{ if ($3 > 19) print $0}'< depth/$z.txt | wc -l)
  # echo $z $BP20 >> depthtab.txt
  awk -v var=$z '{print var, $0}' OFS='\t' depth/$z.txt >> depthlong.txt
done
rm chroms.txt
seq 1 19 >> chroms.txt
echo X >> chroms.txt
echo Y >> chroms.txt
for i in `cat chroms.txt`
do
  awk -v var=$i '{ if ($2 ==var) print $0 }' depthlong.txt >> depthlong/depth_$i.txt
done
  


