#!/bin/bash
set -e
set -u
set -o pipefail

declare sample1_name= $2
declare sample2_name= $3

declare proportion= $1 #change output file name

declare i num_read_sample1=3610000 #change number of reads to get from sample1
declare i num_read_sample2=190000 #change number of reads to get from sample2

#file paths
RAW='/mnt/storage7/lwang/trial_tb_philippines/data/raw_fastq'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/data/seqtk/'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis'

#Sample paths
SAMPLE1_READ1='/mnt/storage7/lwang/trial_tb_philippines/data/raw_fastq/ERR6634978_1.fastq.gz'
SAMPLE1_READ2='/mnt/storage7/lwang/trial_tb_philippines/data/raw_fastq/ERR6634978_2.fastq.gz'
SAMPLE2_READ1='/mnt/storage7/lwang/trial_tb_philippines/data/raw_fastq/ERR6635032_1.fastq.gz'
SAMPLE2_READ2='/mnt/storage7/lwang/trial_tb_philippines/data/raw_fastq/ERR6635032_2.fastq.gz'

#activating fast2matrix enviroment
echo "===activating variant_detection enviroment==="
eval "$(conda shell.bash hook)"
conda activate tb-profiler

mkdir -p $PROCESSED
cd $PROCESSED

# echo '====number of reads===='
# echo 'SAMPLE1_READ1' >> log.txt
# echo $(zcat $SAMPLE1_READ1|wc -l)/4|bc >> log.txt
# echo 'SAMPLE1_READ2' >> log.txt
# echo $(zcat $SAMPLE1_READ2|wc -l)/4|bc >> log.txt
# echo 'SAMPLE2_READ1' >> log.txt
# echo $(zcat $SAMPLE2_READ1|wc -l)/4|bc >> log.txt
# echo 'SAMPLE2_READ2' >> log.txt
# echo $(zcat $SAMPLE2_READ2|wc -l)/4|bc >> log.txt

echo '====subsampling===='
seqtk sample -s100 $SAMPLE1_READ1 $num_read_sample1 > $sample1_name-$sample2_name-${proportion}_1.fastq
seqtk sample -s100 $SAMPLE2_READ1 $num_read_sample2 >> $sample1_name-$sample2_name-${proportion}_1.fastq
seqtk sample -s100 $SAMPLE1_READ2 $num_read_sample1 > $sample1_name-$sample2_name-${proportion}_2.fastq
seqtk sample -s100 $SAMPLE2_READ2 $num_read_sample2 >> $sample1_name-$sample2_name-${proportion}_2.fastq

pigz $sample1_name-$sample2_name-${proportion}_1.fastq
pigz $sample1_name-$sample2_name-${proportion}_2.fastq


echo 'Input: num_read_sample1=' $num_read_sample1 >> log.txt
echo 'Input: num_read_sample2=' $num_read_sample2 >> log.txt

echo 'Ouput:' $sample1_name-$sample2_name-${proportion}_1.fastq >> log.txt
echo $(zcat $sample1_name-$sample2_name-${proportion}_1.fastq.gz|wc -l)/4|bc >> log.txt

echo 'Ouput:' $sample1_name-$sample2_name-${proportion}_2.fastq >> log.txt
echo $(zcat $sample1_name-$sample2_name-${proportion}_2.fastq.gz|wc -l)/4|bc >> log.txt
echo ' ' >> log.txt

echo '***Programme finished***' ${proportion}
