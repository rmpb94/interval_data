## Merge batch_a and batch_b files into one chromosome
## Output QC files

#!/bin/bash
#SBATCH --job-name=chr_merge
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/merges/out/chr_merge_all_out.out
#SBATCH --time=08:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
input_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/imputed/"
output_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/post_impute/"

bcftools merge -O z ${input_dir}batch_a_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz ${input_dir}batch_b_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz -o ${output_dir}interval_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz
bcftools index ${output_dir}interval_chr${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools query -f '%CHROM\t%POS\t%INFO/AF\t%INFO/MAF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/ER2\t%INFO/IMPUTED\t%INFO/TYPED\n' ${output_dir}interval_chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${output_dir}qc/chr${SLURM_ARRAY_TASK_ID}_qc_stats.txt
