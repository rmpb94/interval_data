##############################################
##  Lift over                               ##  
##############################################
##  Interval data is GRCh37, needs to be 38, before imputation
##  Nick Gleadhall was using this file ~/rds/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed
## UKBB_v2_contents_20180502.hg19_to_hg38.txt
##  Olga had noticed some errors with file, so I'll need to check it

############################################
##  Make some checks on the command line  ##
############################################
##   Print the first line of the file to understand what columns are in it
data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/reference_files/genetic/reference_files_genotyped_imputed/"
head -1 ${data_dir}UKBB_v2_contents_20180502.hg19_to_hg38.txt 
##  probeset_id     affy_snp_id     rsid    chr     pos_hg19_vcf    ref_hg19_vcf    alt_hg19_vcf    pos_hg38_vcf    ref_hg38_vcf alt_hg38_vcf    source  category


## Put all of the following into a job script
## Jobs are in /home/rmb218/hpc-work/jobs/imputation
interval_liftover.sh

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
work_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/"
info_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/"
interval_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/"
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
plink2 \
    --bfile ${interval_dir}pre_impute_interval \
    --update-name ${info_dir}all_snps_data.txt 5 1 \
    --make-bed \
    --out ${interval_dir}pre_impute_temp
    
plink2 --bfile ${interval_dir}pre_impute_temp --write-snplist --out ${info_dir}updated_snps

##  Keep only SNPs which I have update info for
plink2 \
    --bfile ${interval_dir}pre_impute_temp \
    --extract ${info_dir}rb_snplist2024.txt \
    --make-bed \
    --out ${interval_dir}pre_impute_interval_hg19_trimmed

##  Update positions
plink2 \
    --bfile ${interval_dir}pre_impute_interval_hg19_trimmed \
    --update-map ${info_dir}all_snps_data.txt 4 5 \
    --sort-vars 'n' \
    --make-pgen \
    --out ${interval_dir}pre_impute_interval_hg38

