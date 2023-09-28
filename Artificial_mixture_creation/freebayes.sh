#freebayes bam to vcf conversion
#local path: cd ~/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis

RAW='/mnt/storage7/lwang/trial_tb_philippines/data/seqtk'
PROCESSED='/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk'
REFGENOME='/mnt/storage7/lwang/trial_tb_philippines/refgenome/MTB-h37rv_asm19595v2-eg18.fa'
PIPELINE='/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis'

#activating variant_detection enviroment
echo 'activating variant_detection enviroment'
eval "$(conda shell.bash hook)"
conda activate variant_detection 

mkdir -p $PROCESSED/freebayesVCF
cd $PROCESSED/freebayesVCF

#Converting bam to vcf file
echo 'converting .bam to .vcf for sample' $1 
freebayes -f $REFGENOME $PROCESSED/$1.mkdup.bam > $PROCESSED/freebayesVCF/$1.vcf

#freebayes -f $REFGENOME $PROCESSED/pool_f4m/ERR6634978.mkdup.bam > $PROCESSED/freebayesVCF/ERR6634978.vcf

#Quality filtering
# echo 'Quality filtering' $1 
# vcftools --vcf $PROCESSED/freebayesVCF/q20/$1.vcf --minQ 20 --recode --recode-INFO-all --out $1_q20
bgzip $1.vcf

# vcftools --vcf $PROCESSED/freebayesVCF/q20/$1.vcf --minQ 20 --recode --recode-INFO-all --out $1_q20
# bgzip $1_q20.recode.vcf
