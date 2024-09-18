## Indexing ####
## Index all batch A
batch_a_index.sh
#!/bin/bash
#SBATCH --job-name=batch_a_index
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/out/batch_a_index_%a.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/out/batch_a_index_error_%a.out
#SBATCH --time=12:00:00
#SBATCH --partition=icelake-himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
batch_a_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/imputed/"
bcftools index ${batch_a_dir}chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz 

batch_b_index.sh
#!/bin/bash
#SBATCH --job-name=batch_a_index
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/out/batch_b_index_%a.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/out/batch_b_index_error_%a.out
#SBATCH --time=12:00:00
#SBATCH --partition=icelake-himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
batch_b_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/"
bcftools index ${batch_b_dir}chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz 

## Quality control filtering ####
## What information is in the vcf file?
zcat ${batch_a_dir}chr22.dose.vcf.gz | more

## Make initial QC files from raw imputed batch files
batch_qc_files.sh
#!/bin/bash
#SBATCH --job-name=batch_qc_files
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/qc/out/batches_qc.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/qc/out/batches_qc.out
#SBATCH --time=05:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy

batch_a_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/imputed/"
a_out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/"

batch_b_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/"
b_out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/"

## Make stripped back files for QC from both batch a and batch b
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/ER2\t%INFO/IMPUTED\t%INFO/TYPED\n' ${batch_a_dir}chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${a_out_dir}ba_chr${SLURM_ARRAY_TASK_ID}_qc_stats.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/ER2\t%INFO/IMPUTED\t%INFO/TYPED\n' ${batch_b_dir}chr${SLURM_ARRAY_TASK_ID}.dose.vcf.gz > ${b_out_dir}bb_chr${SLURM_ARRAY_TASK_ID}_qc_stats.txt

## Look at QC output to decide thresholds for analysis
post_qc_plots.sh
#!/bin/bash
#SBATCH --job-name=post_imputation_qc_plot
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/qc/out/post_imputation_qc_plot.out
#SBATCH --error=/home/rmb218/hpc-work/jobs/imputation/qc/out/post_imputation_qc_plot.out
#SBATCH --time=00:30:00
#SBATCH --partition=icelake-himem
#SBATCH --mem=30000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load R
cd /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/
Rscript /home/rmb218/hpc-work/r_scripts/interval_qc/INTERVALChromosomalQC.R