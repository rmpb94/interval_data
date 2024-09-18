## BCFtools liftover tool test ####
## Dowload precompiled binaries
## Working directory
base_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/"
cd ${base_dir}

mkdir ${base_dir}imputation ${base_dir}interval

##  Copy files that are needed across from /home/rmb218/rds/hpc-work/
mkdir /home/rmb218/hcp-work/imputation/
cp /home/rmb218/rds/hpc-work/interval/imputation/liftover/testing/* 



module load bcftools-1.9-gcc-5.4.0-b2hdt5n

## Unzip download
cd bcftools_liftover
unzip /home/rmb218/hpc-work/imputation/bcftools_liftover/score_1.20-20240505.zip

## Define path to plugin
BCFTOOLS_PLUGINS="/home/rmb218/hpc-work/imputation/bcftools_liftover"

bcftools plugin $BCFTOOLS_PLUGINS/liftover.so