#!/bin/bash
# 2023-03-08 11:21:00
conda activate bedtools # also have bedops

VCF='/central/groups/carnegie_poc/lczech/1001g/greneNet_final.vcf.gz'
SINK='/central/groups/carnegie_poc/meixilin/grenenet/analyses/data-raw/seedmix_vcf/'

rsync -ahv $VCF ${SINK}

cd $SINK
zcat greneNet_final.vcf.gz | vcf2bed | cut -f 1-3 > greneNet_final.vcf.bed

# yep! still the file is not outputting the right files
