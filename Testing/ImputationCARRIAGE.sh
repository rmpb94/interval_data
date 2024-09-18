############################################
##  INTERVAL CARRIAGE TOPMed imputation ####
############################################
module load plink/2.00-alpha

##  Data directories    ####
qc_dir="/home/rmb218/rds/hpc-work/interval/imputation/qc_files/"
##  sample_ids.id, variant_ids.snplist
raw_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
###################################
##  Pre-imputation files       ####
###################################
##  Split samples into 2 batches
shuf ${qc_dir}sample_ids.id -o ${qc_dir}sample_ids_shuffled.id
total_lines=$(wc -l < ${qc_dir}sample_ids_shuffled.id)
midpoint=$(( (total_lines + 1) / 2))
split -l $midpoint ${qc_dir}sample_ids_shuffled.id ${qc_dir}samples_

interval_pre_impute_files.sh
#!/bin/bash
#SBATCH --job-name=interval_pre_impute
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputation_jobs/out/interval_pre_impute.out
#SBATCH --time=05:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load plink/2.00-alpha
raw_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
qc_dir="/home/rmb218/rds/hpc-work/interval/imputation/qc_files/"
batch_out="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_imputation/pre_impute/"
##  Batch A
plink2 \
    --bfile ${raw_data_dir}interval_qced_24.8.18 \
    --keep ${qc_dir}samples_aa \
    --extract ${qc_dir}variant_ids.snplist \
    --chr ${SLURM_ARRAY_TASK_ID} \
    --recode vcf id-paste=iid bgz \
    --snps-only just-acgt \
    --out ${batch_out}A/interval_batch_a_chr${SLURM_ARRAY_TASK_ID}

##  Batch B
plink2 \
    --bfile ${raw_data_dir}interval_qced_24.8.18 \
    --keep ${qc_dir}samples_ab \
    --extract ${qc_dir}variant_ids.snplist \
    --chr ${SLURM_ARRAY_TASK_ID} \
    --recode vcf id-paste=iid bgz \
    --snps-only just-acgt \
    --out ${batch_out}B/interval_batch_b_chr${SLURM_ARRAY_TASK_ID}
