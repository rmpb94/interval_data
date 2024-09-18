## Create `.snpstats` files of current INTERVAL data
## Using QCtool
sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
module load qctool/v2.0.5
module load plink/2.00-alpha
## Make PLINK files for batch_a, batch_b, and chr22 merged

## Directories ####
## This is to the .vcf files post-filtering for MAF >= 0.005 and R2 >= 0.7
batch_data_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/imputed/qc_test/"
## This is merged batch_a and batch_b into one VCF file
chrom_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/imputed/"
out_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval/data/imputed/qc_out/"

plink2 --vcf ${batch_data_dir}batch_a_chr22.dose.vcf.gz dosage=DS --make-bed --out ${out_dir}batch_a_chr22
plink2 --vcf ${batch_data_dir}batch_b_chr22.dose.vcf.gz dosage=DS --make-bed --out ${out_dir}batch_b_chr22
plink2 --vcf ${chrom_dir}interval_chr22.dose.vcf.gz dosage=DS --make-bed --out ${out_dir}interval_chr22

qctool -g ${out_dir}batch_a_chr22.bed -snp-stats -osnp ${out_dir}batch_a_chr22.snpstats
qctool -g ${out_dir}batch_b_chr22.bed -snp-stats -osnp ${out_dir}batch_b_chr22.snpstats
qctool -g ${out_dir}interval_chr22.bed -snp-stats -osnp ${out_dir}interval_chr22.snpstats

## Using previous (hg19) imputation stats to compare new imputed data
## Access `.snpstats` files from previous imputation
file_dir="/home/rmb218/rds/rds-post_qc_data-pNR2rM6BWWA/interval/reference_files/genetic/reference_files_genotyped_imputed/"