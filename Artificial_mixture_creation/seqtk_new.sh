#!/bin/bash
set -e
set -u
set -o pipefail

# test run
# ./seqtk_new.sh 50 50 ERR2503407 ERR2864238
echo '====Load variables====='

#file paths
RAW='/mnt/storage7/lwang/trial_tb_philippines/data/insilico_mix_prep'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/data/insilico_mix_mix'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/Artificial_mixture_creation/'

sample1_proportion=$1 #change output file name
sample2_proportion=$2
proportion=$1$2

sample1_name=$3
sample2_name=$4

#Sample paths
SAMPLE1_READ1=$RAW/${sample1_name}_1.fastq.gz
SAMPLE1_READ2=$RAW/${sample1_name}_2.fastq.gz
SAMPLE2_READ1=$RAW/${sample2_name}_1.fastq.gz
SAMPLE2_READ2=$RAW/${sample2_name}_2.fastq.gz


# echo '====get depth====='

# depth1=`gunzip -c $SAMPLE1_READ1 | wc -l`
# depth2=`gunzip -c $SAMPLE2_READ1 | wc -l`

# let num_read_sample1=$depth1/100*$1 #change number of reads to get from sample1
# let num_read_sample2=$depth2/100*$2 #change number of reads to get from sample2

let num_read_sample1=3000000/100*$1 #change number of reads to get from sample1
let num_read_sample2=3000000/100*$2 #change number of reads to get from sample2

#activating fast2matrix enviroment
# echo "===activating variant_detection enviroment==="
# eval "$(conda shell.bash hook)"
# conda activate tb-profiler

# echo 'environment activated'

# mkdir -p $PROCESSED
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

echo '***Programme finished***' $sample1_name-$sample2_name-${proportion}