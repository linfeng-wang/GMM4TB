#generate strain variant json file using tb-profiler
#local path: cd ~/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis

RAW='/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/json_result'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis'

#reference file
REFGENOME='refgenome/MTB-h37rv_asm19595v2-eg18.fa'

#activating variant_detection enviroment
echo '===activating tb-profiler enviroment==='
eval "$(conda shell.bash hook)"
conda activate tb-profiler

cd $RAW

#generating json file
echo '***generating json file from *tb-profiler vcf_profile* command***'
# find . -name '*.recode.vcf.gz' -exec tb-profiler vcf_profile {} \;
find . -name '*.vcf.gz' -exec tb-profiler vcf_profile {} \;
