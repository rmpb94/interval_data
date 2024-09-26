#######################################################################################
## PRE IMPUTATION QC                                                                 ##
#######################################################################################
# Interactive node: sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
# Original Data: /home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/
# This should have been through QC as described in /home/rmb218/rds/rds-jmmh2-post_qc_data/interval/reference_files_genotyped_imputed/QC_process_log.txt

# Confirming original QC has been sucessful

module load plink/1.9

interval_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/qc/pre_impute/"

plink --bfile ${interval_data_dir}interval_qced_24.8.18 --hardy --out ${out_dir}interval_preimpute
plink --bfile ${interval_data_dir}interval_qced_24.8.18 --freq --out ${out_dir}interval_preimpute
plink --bfile ${interval_data_dir}interval_qced_24.8.18 --missing --out ${out_dir}interval_preimpute
plink --bfile ${interval_data_dir}interval_qced_24.8.18 --het --out ${out_dir}interval_preimpute
plink  --bfile ${interval_data_dir}interval_qced_24.8.18 --indep 50 5 2 --out ${out_dir}interval_preimpute
plink --bfile ${interval_data_dir}interval_qced_24.8.18  --extract ${out_dir}interval_preimpute.prune.in --make-bed --out ${out_dir}interval_preimpute_pruned

#!/bin/bash
#SBATCH --job-name=pre_impute_pca
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/qc/out/pre_impute_pca.out
#SBATCH --time=01:00:00
#SBATCH --mem=30000
#SBATCH --nodes=8
#SBATCH --tasks=8
#SBATCH --partition=icelake-himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load plink/1.9
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/qc/pre_impute/"
plink --bfile ${out_dir}interval_preimpute_pruned --Z-genome --memory 30000 --out ${out_dir}interval_preimpute_relationships

