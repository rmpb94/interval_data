#######################################################################################
## POST IMPUTATION QC                                                                ##
#######################################################################################
## Following: https://www.biostars.org/p/335605/
## INTERVAL SNP-genotype data has been lifted to hg38, imputed to TOPMed r3, and merged per chromosome
## Now this needs to go through several QC stages before I can carry it forward to CARRIAGE
sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR

## Full INTERVAL population PCA ####

dir_path="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_analysis/pca/"

ref_download.sh
#!/bin/bash
#SBATCH --job-name=ref_download
#SBATCH --output=/home/rmb218/hpc-work/jobs/interval/pca/out/ref_download.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
cd /rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_analysis/pca/ref_files/
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr" ;
suffix=".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" ;
for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done

## Step 2
## Download 1k genomes ped file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;

## This is the link for downloading the hg38 reference genome https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22)) 

## Convert 1k genome reference files into PLINK format
ref_plink_format.sh
#!/bin/bash
#SBATCH --job-name=ref_plink_format
#SBATCH --output=/home/rmb218/hpc-work/jobs/interval/pca/out/ref_plink_format.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
cd /rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_analysis/pca/ref_files/
module load plink/2.00-alpha
module load bcftools 
## Make sure variants have IDs
# Annotate the VCF file and output to a temporary file
bcftools annotate ALL.chr${SLURM_ARRAY_TASK_ID}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz -x ID -I +'%CHROM:%POS:%REF:%ALT' -O z -o ALL.chr${SLURM_ARRAY_TASK_ID}.annotated.vcf.gz

# Normalize the annotated VCF and output to a temporary file
bcftools norm --rm-dup both ALL.chr${SLURM_ARRAY_TASK_ID}.annotated.vcf.gz \
    -O z -o ALL.chr${SLURM_ARRAY_TASK_ID}.normalized.vcf.gz

# Index the normalized VCF file
bcftools index ALL.chr${SLURM_ARRAY_TASK_ID}.normalized.vcf.gz

plink2 \
    --vcf ALL.chr${SLURM_ARRAY_TASK_ID}.normalized.vcf.gz \
    --make-pgen \
    --out 1kg_chr${SLURM_ARRAY_TASK_ID}_ref_genome_hg38

## Prune variants
ref_prune_variants.sh
#!/bin/bash
#SBATCH --job-name=ref_prune_variants
#SBATCH --output=/home/rmb218/hpc-work/jobs/interval/pca/out/ref_prune_variants.out
#SBATCH --time=01:00:00
#SBATCH --partition=icelake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#SBATCH --array=1-22
cd /rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_analysis/pca/ref_files/
module load plink/2.00-alpha

plink2 \
    --pfile 1kg_chr${SLURM_ARRAY_TASK_ID}_ref_genome_hg38 \
    --maf 0.10 --indep-pairwise 50 5 0.5 \
    --out Pruned/1kg_chr${SLURM_ARRAY_TASK_ID}_ref_genome_hg38 

plink2 \
    --pfile 1kg_chr${SLURM_ARRAY_TASK_ID}_ref_genome_hg38 \
    --extract Pruned/1kg_chr${SLURM_ARRAY_TASK_ID}_ref_genome_hg38.prune.in \
    --make-bed \
    --out Pruned/1kg_chr${SLURM_ARRAY_TASK_ID}_ref_genome_hg38

## List all files for PLINK
cd /rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_analysis/pca/ref_files/
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;
sed -i 's/.bim//g' ForMerge.list

ref_merge_files.sh
#!/bin/bash
#SBATCH --job-name=ref_merge_files
#SBATCH --output=/home/rmb218/hpc-work/jobs/interval/pca/out/ref_merge_files.out
#SBATCH --time=01:00:00
#SBATCH --partition=icelake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
cd /rds/project/rds-zuZwCZMsS0w/carriage_gwas/interval_analysis/pca/ref_files/
module load plink/1.9
plink --merge-list ForMerge.list --out Merge
## Write out list of SNPs
plink2 --bfile Merge --write-snplist --out reference

## Prep INTERVAL data ####

#!/bin/bash
#SBATCH --job-name=pca_interval_plink
#SBATCH --output=/home/rmb218/hpc-work/jobs/interval/pca/out/pca_interval_plink.out
#SBATCH --time=01:00:00
#SBATCH --partition=icelake
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load plink/2.00-alpha
int_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval_data/pre_impute/"
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval_analysis/pca/int_files/"

plink2 --bfile ${int_dir}interval_refhg38 --maf 0.1 --make-bed --out ${out_dir}interval_refhg38
plink2 --bfile ${out_dir}interval_refhg38 --indep-pairwise 50 5 0.5 --out ${out_dir}interval
plink2 --bfile ${out_dir}interval_refhg38 --extract ${out_dir}interval.prune.in --make-bed --out ${out_dir}interval_pruned_hg38

plink2 --bfile ${out_dir}interval_pruned_hg38 --write-snplist --out ${out_dir}interval

## Combine Data ####
## Extract SNPs from both reference and INTERVAL cohorts
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval_analysis/pca/"

plink2 --bfile ${out_dir}int_files/interval_pruned_hg38 --extract ${out_dir}pca_snps.txt --make-bed --out ${out_dir}combined_files/interval
plink2 --bfile ${out_dir}ref_files/Merge --extract ${out_dir}ref_pca_snps.txt --make-bed --out ${out_dir}combined_files/1kg

plink2 --bfile ${out_dir}combined_files/1kg --update-name ${out_dir}pca_snps_rename.txt 1 2 --make-bed --out ${out_dir}combined_files/1kg_renamed
plink2 --bfile ${out_dir}combined_files/1kg_renamed --write-snplist --out ${out_dir}combined_files/test

## Files for merge
/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval_analysis/pca/combined_files/interval
/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval_analysis/pca/combined_files/1kg_renamed
module unload plink
module load plink/1.9
plink --bfile ${out_dir}combined_files/interval --bmerge ${out_dir}combined_files/1kg_renamed --out ${out_dir}combined_files/int_pca_combined

pca.sh
#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --output=/home/rmb218/hpc-work/jobs/interval/pca/out/pca.out
#SBATCH --time=05:00:00
#SBATCH --partition=icelake-himem
#SBATCH --mem=30000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load plink/1.9
out_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval_analysis/pca/"
plink --bfile ${out_dir}combined_files/int_pca_combined --pca --out ${out_dir}int_pca