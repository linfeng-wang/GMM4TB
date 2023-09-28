#!/bin/bash
set -e
set -u
set -o pipefail

RAW='/mnt/storage7/jody/tb_ena/per_sample'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/QuantTB/in_vitro_samples/fastq'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/Artificial_mixture_creation/'


cd $PROCESSED

# cat $PROCESSED/mix_clinical_name.txt | parallel -j 30 "samtools sort -n $RAW/{}.bqsr.cram | samtools fastq --reference $REFGENOME -1 {}_1.fastq.gz -2 {}_2.fastq.gz -"
cat '/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable_eval/invitro_mix_name.csv' | parallel -j 30 "samtools sort -n $RAW/{}.bqsr.cram | samtools fastq --reference $REFGENOME -1 {}_1.fastq.gz -2 {}_2.fastq.gz -"

# < - > hyphen means to take input from standard in

# samtools sort -n /mnt/storage7/jody/tb_ena/per_sample/ERR2864244.bqsr.cram | samtools fastq --reference /mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa -1 ERR2864244_1.fastq.gz -2 ERR2864244_2.fastq.gz -

# samtools sort -n /mnt/storage7/jody/tb_ena/per_sample/ERR751749.bqsr.cram | samtools fastq --reference /mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa -1 ERR751749_1.fastq.gz -2 ERR751749_2.fastq.gz -