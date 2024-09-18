## Merge batches together, output stats file
chr_22_merge.sh
#!/bin/bash
#SBATCH --job-name=chr22_merge
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/qc/out/chr22_merge.out
#SBATCH --time=10:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
input_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/imputed/qc_test/"
output_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/imputed/"
bcftools merge -O z ${input_dir}batch_a_chr22.dose.vcf.gz ${input_dir}batch_b_chr22.dose.vcf.gz --output ${output_dir}interval_chr22.dose.vcf.gz
bcftools index ${output_dir}interval_chr22.vcf.gz
bcftools query -f '%CHROM\t%POS\t%INFO/AF\t%INFO/MAF\t%INFO/AVG_CS\t%INFO/R2\t%INFO/ER2\t%INFO/IMPUTED\t%INFO/TYPED\n' ${output_dir}interval_chr22.dose.vcf.gz > ${input_dir}chr22_qc_stats.txt

## Assess merge ####
## To look at stats, easier to convert this vcf into a plink file
## Then calculate missingness etc from there
module load plink/2.00-alpha
## Paths
chrom_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/imputed/"
out_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/imputed/qc_out/"

## Convert file into plink binary files
plink2 --vcf ${chrom_dir}interval_chr22.dose.vcf.gz dosage=DS --make-bed --out ${out_dir}interval_chr22

## Then lpoad plink 1.9 to run QC commands
sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
out_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/imputed/qc_out/"
module load plink/1.9

plink --bfile ${out_dir}interval_chr22 --hardy --out ${out_dir}interval_chr22
plink --bfile ${out_dir}interval_chr22 --freq --out ${out_dir}interval_chr22
plink --bfile ${out_dir}interval_chr22 --missing --out ${out_dir}interval_chr22
plink --bfile ${out_dir}interval_chr22 --het --out ${out_dir}interval_chr22
plink --bfile ${out_dir}interval_chr22 --indep 50 5 2 --out ${out_dir}interval_chr22
plink --bfile ${out_dir}interval_chr22  --extract ${out_dir}interval_chr22.prune.in --make-bed --out ${out_dir}interval_chr22_pruned

plink --bfile ${out_dir}interval_chr22 --geno 0.03 --maf 0.005 --hwe 0.00000001 --make-bed --out ${out_dir}interval_chr22_qced
    ##  Total genotyping rate is 0.978789.
    ##  15664 variants removed due to missing genotype data (--geno).
    ##  --hwe: 2 variants removed due to Hardy-Weinberg exact test.
    ##  2855 variants removed due to minor allele threshold(s)
    ##  (--maf/--max-maf/--mac/--max-mac).
    ##  108193 variants and 42396 people pass filters and QC.