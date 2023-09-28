#freebayes bam->vcf conversion parallel file running script
#local path: cd ~/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis


RAW='/mnt/storage7/lwang/trial_tb_philippines/data/seqtk'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis'

#reference file
REFGENOME='refgenome/MTB-h37rv_asm19595v2-eg18.fa'

#activating variant_detection enviroment
echo '===activating variant_detection enviroment==='
eval "$(conda shell.bash hook)"
conda activate variant_detection 

mkdir -p $PROCESSED/freebayesVCF
cd $PROCESSED/freebayesVCF


cat $RAW/sample_name.txt | parallel -j 5 "$PIPELINE/freebayes.sh {}"
