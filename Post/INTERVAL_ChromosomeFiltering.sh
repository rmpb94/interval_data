## Filter chromosome 22 in batches separately
batch_a_filtering._sh
#!/bin/bash
#SBATCH --job-name=batch_a_filtering
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/filtering/out/batch_a_filtering.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/filtering/out/batch_filtering.out
#SBATCH --time=24:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL2-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
batch_a_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/imputed/"
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/imputed/"
bcftools view -i 'INFO/MAF > 0.005 && INFO/R2 > 0.7' ${batch_a_dir}chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${out_dir}batch_a_chr${SLURM_ARRAY_TASK_ID}.dose.vcf
bgzip ${out_dir}batch_a_chr${SLURM_ARRAY_TASK_ID}.dose.vcf
bcftools index ${out_dir}batch_a_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz 

batch_b_filtering.sh
#!/bin/bash
#SBATCH --job-name=batch_b_filtering
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/filtering/out/batch_b_filtering.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/filtering/out/batch_b_filtering.out
#SBATCH --time=24:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL2-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
batch_b_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/"
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/imputed/"
bcftools view -i 'INFO/MAF > 0.005 && INFO/R2 > 0.7' ${batch_b_dir}chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${out_dir}batch_b_chr${SLURM_ARRAY_TASK_ID}.dose.vcf
bgzip -f ${out_dir}batch_b_chr${SLURM_ARRAY_TASK_ID}.dose.vcf
bcftools index -f ${out_dir}batch_b_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz 


## Create QC files ####
batch_qc_files.sh
#!/bin/bash
#SBATCH --job-name=batch_qc_files
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/filtering/out/batch_qc_files.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/filtering/out/batch_qc_files.out
#SBATCH --time=10:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
data_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/imputed/"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/ER2\t%INFO/IMPUTED\t%INFO/TYPED\n' ${data_dir}batch_a_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${data_dir}ba_chr${SLURM_ARRAY_TASK_ID}_qc_stats.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/ER2\t%INFO/IMPUTED\t%INFO/TYPED\n' ${data_dir}batch_b_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${data_dir}bb_chr${SLURM_ARRAY_TASK_ID}_qc_stats.txt

## Using the above files, plot MAF vs R2 graphs
