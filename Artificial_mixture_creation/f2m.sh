#fast2matrix analysis
#local path: cd ~/trial_tb_philippines/pipelines/Genomic_data_analysis/trial_pipe


RAW='/mnt/storage7/lwang/trial_tb_philippines/data/seqtk'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis'


##activating fast2matrix enviroment
eval "$(conda shell.bash hook)"
conda activate fastq2matrix

mkdir -p $PROCESSED
cd $PROCESSED
#Maping out variant
echo "***Maping out variant***"
fastq2vcf.py all --read1 $RAW/$1_1.fastq.gz --read2 $RAW/$1_2.fastq.gz --ref $REFGENOME --prefix $1

#End of process for this file
echo "======End of process for $1======"
# cp $PROCESSED/sep_f4m/$1/* $PROCESSED/
