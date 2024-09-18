##############################################
##  Lift over                               ##  
##############################################
##  Interval data is GRCh37, needs to be 38, before imputation
## sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
## module load plink/2.00-alpha
## Jobs are in /home/rmb218/hpc-work/jobs/imputation
sbatch interval_liftover.sh

#!/bin/bash
#SBATCH --job-name=interval_liftover
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/out/interval_liftover.out
#SBATCH --time=01:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
############################################
##  Make pre imputation files             ##
############################################
##  From previous QC performed in 2018, list 655,045 variants and 42,396 individuals
module load plink/2.00-alpha
## These file paths are ultimately the ones that should be used, but currently there's no space in the prject folder
## For now, put these in my space /home/rmb218/hpc-work/
work_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/"
info_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/"
interval_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/pre_impute/"
mkdir ${interval_dir}
clean_interval_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/"
## raw_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
## data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/"
#work_dir="/home/rmb218/hpc-work/"
#info_dir="/home/rmb218/hpc-work/interval/info/"
#interval_dir="/home/rmb218/hpc-work/interval/"
raw_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/"

plink2 --bfile ${raw_data_dir}interval_qced_24.8.18 --write-samples --out ${info_dir}sample_ids_2018
plink2 --bfile ${raw_data_dir}merged_imputation --write-snplist --out ${info_dir}variant_ids_2018

##  Pre-imputation files
plink2 \
    --bfile ${raw_data_dir}interval_qced_24.8.18 \
    --keep ${info_dir}sample_ids_2018.id \
    --extract ${info_dir}variant_ids_2018.snplist \
    --make-bed \
    --out ${interval_dir}pre_impute_interval

############################################
##  Matching SNPs                         ##
############################################
## Locally, made 2 lists of all SNPs in INTERVAL data
## my_snps_rsids.txt, my_snps_affx.txt
## Extract these from the big UKBB file
## Update reference and alternative alleles 
## Update SNP names
## Make all names
## Transfer files from local R project into /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/

## Update SNP names from old names to chr:poshg38 names
plink2 \
    --bfile ${interval_dir}pre_impute_interval \
    --update-name ${info_dir}all_snps_data.txt 5 1 \
    --make-bed \
    --out ${interval_dir}pre_impute_temp

## Write out list of updated variants    
plink2 --bfile ${interval_dir}pre_impute_temp --write-snplist --out ${info_dir}updated_snps

## Make pgen files of the same data
plink2 \
    --bfile ${interval_dir}pre_impute_temp \
    --make-pgen \
    --out ${interval_dir}pre_impute_temp

##  Keep only SNPs which I have update info for
plink2 \
    --pfile ${interval_dir}pre_impute_temp \
    --extract ${info_dir}rb_snplist2024.txt \
    --make-bed \
    --out ${interval_dir}pre_impute_interval_hg19_trimmed

##  Update positions
plink2 \
    --bfile ${interval_dir}pre_impute_interval_hg19_trimmed \
    --update-map ${info_dir}all_snps_data.txt 4 5 \
    --sort-vars 'n' \
    --make-pgen \
    --out ${clean_interval_dir}pre_impute_interval_hg38

## Remove all intermediate files to free up space
rm -R /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/intermediate/

plink2 \
    --pfile ${clean_interval_dir}pre_impute_interval_hg38 \
    --make-bed \
    --out ${clean_interval_dir}pre_impute_interval_hg38
