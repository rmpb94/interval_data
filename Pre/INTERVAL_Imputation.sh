##  Modules
sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
module load plink/2.00-alpha
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy

#!/bin/bash
#SBATCH --job-name=pre_impute_processing
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/out/pre_impute_processing.out
#SBATCH --time=01:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load plink/2.00-alpha

## Re-align to reference genome before doing anything
#cd /home/rmb218/hpc-work/ref_files
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# use shared drive paths to store files
#interval_dir="/home/rmb218/hpc-work/interval/"
interval_shared_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/"
plink2 \
  --bfile ${interval_shared_dir}pre_impute_interval_hg38 \
  --fa /home/rmb218/rds/hpc-work/ref_files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
  --ref-from-fa 'force' \
  --make-bed \
  --out ${interval_shared_dir}interval_refhg38

## Split the cohort in half (batch a and batch b)
info_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/"
plink2 \
    --bfile ${interval_shared_dir}interval_refhg38 \
    --write-samples \
    --out ${info_dir}full_sample_list

##  Generate first random number list
shuf -i 1-42396 -n 21198 > ${info_dir}interval_random_split.txt
##  Use numbers to extract individuals
awk 'NR==FNR {lines[$1]; next} FNR in lines' ${info_dir}interval_random_split.txt ${info_dir}full_sample_list.id > ${info_dir}batch_a_individuals.txt

mkdir ${interval_shared_dir}batch_a
mkdir ${interval_shared_dir}batch_b
## Batch A
plink2 \
    --bfile ${interval_shared_dir}interval_refhg38 \
    --keep ${info_dir}batch_a_individuals.txt \
    --make-bed \
    --out ${interval_shared_dir}batch_a/pre_impute_interval_batch_a
## Batch B
plink2 \
    --bfile ${interval_shared_dir}interval_refhg38 \
    --remove ${info_dir}batch_a_individuals.txt \
    --make-bed \
    --out ${interval_shared_dir}batch_b/pre_impute_interval_batch_b

batch_a_split_chromosomes.sh
#!/bin/bash
#SBATCH --job-name=split_batch_a_chr
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/split_batch_a_chr.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load plink/2.00-alpha
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy

interval_batch_a="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/"
info_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/"
interval_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/"
mkdir ${interval_batch_a}chromosomes/

## Split by chromosomes
plink2 \
  --bfile ${interval_batch_a}pre_impute_interval_batch_a \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --make-bed \
  --out ${interval_batch_a}chromosomes/batch_a_${SLURM_ARRAY_TASK_ID}
plink2 \
  --bfile ${interval_batch_a}chromosomes/batch_a_${SLURM_ARRAY_TASK_ID} \
  --recode vcf id-paste=iid bgz \
  --snps-only just-acgt \
  --out ${interval_batch_a}chromosomes/batch_a_${SLURM_ARRAY_TASK_ID}
bcftools index ${interval_batch_a}chromosomes/batch_a_${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "${SLURM_ARRAY_TASK_ID} chr${SLURM_ARRAY_TASK_ID}" > ${info_dir}chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt

bcftools annotate --rename-chrs ${info_dir}chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt ${interval_batch_a}chromosomes/batch_a_${SLURM_ARRAY_TASK_ID}.vcf.gz -o ${interval_batch_a}chromosomes/batch_a_chr${SLURM_ARRAY_TASK_ID}.vcf

## Remove old vcf files
rm ${interval_batch_a}chromosomes/batch_a_${SLURM_ARRAY_TASK_ID}.vcf.gz

## Zip and index new vcf file
bgzip ${interval_batch_a}chromosomes/batch_a_chr${SLURM_ARRAY_TASK_ID}.vcf
bcftools index ${interval_batch_a}chromosomes/batch_a_chr${SLURM_ARRAY_TASK_ID}.vcf.gz
## Remove update file
rm ${info_dir}chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt


############################
## Batch B
############################

batch_b_split_chromosomes.sh
#!/bin/bash
#SBATCH --job-name=split_batch_b_chr
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/split_batch_b_chr.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load plink/2.00-alpha
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy

interval_batch_b="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/"
info_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/"
interval_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/"
mkdir ${interval_batch_b}chromosomes/

## Split by chromosomes
plink2 \
  --bfile ${interval_batch_b}pre_impute_interval_batch_b \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --make-bed \
  --out ${interval_batch_b}chromosomes/batch_b_${SLURM_ARRAY_TASK_ID}
plink2 \
  --bfile ${interval_batch_b}chromosomes/batch_b_${SLURM_ARRAY_TASK_ID} \
  --recode vcf id-paste=iid bgz \
  --snps-only just-acgt \
  --out ${interval_batch_b}chromosomes/batch_b_${SLURM_ARRAY_TASK_ID}
bcftools index ${interval_batch_b}chromosomes/batch_b_${SLURM_ARRAY_TASK_ID}.vcf.gz

echo "${SLURM_ARRAY_TASK_ID} chr${SLURM_ARRAY_TASK_ID}" > ${info_dir}chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt

bcftools annotate --rename-chrs ${info_dir}chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt ${interval_batch_b}chromosomes/batch_b_${SLURM_ARRAY_TASK_ID}.vcf.gz -o ${interval_batch_b}chromosomes/batch_b_chr${SLURM_ARRAY_TASK_ID}.vcf

## Remove old vcf files
rm ${interval_batch_b}chromosomes/batch_b_${SLURM_ARRAY_TASK_ID}.vcf.gz

## Zip and index new vcf file
bgzip ${interval_batch_b}chromosomes/batch_b_chr${SLURM_ARRAY_TASK_ID}.vcf
bcftools index ${interval_batch_b}chromosomes/batch_b_chr${SLURM_ARRAY_TASK_ID}.vcf.gz
## Remove update file
rm ${info_dir}chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt