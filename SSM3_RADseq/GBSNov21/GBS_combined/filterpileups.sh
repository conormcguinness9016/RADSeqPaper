#usr/bin/env bash
bcftools view -i "%QUAL>=30 & INFO/DP>20" postBQSR_A1A2.mpileup | wc -l
