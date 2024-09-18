# Make directories
mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/imputed/
mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/topmed_qc/
mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/
mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/topmed_qc/

# Download job
b.sh
#!/bin/bash
#SBATCH --job-name=b1
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/downloads/out/b1.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
#a_imputed_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/imputed/"
#a_qc_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/topmed_qc/"
b_imputed_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/"
b_qc_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/topmed_qc/"

cd ${b_qc_dir}logs/
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/33593fbc519587ac7f04d57dee272f0d23055dcc7fea95b2471191773914b96d/chr_1.log
mkdir ${b_qc_dir}chr1
cd ${b_qc_dir}chr1/
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/c6b29efde3362a9d86c6ba9b473f670fcf476739d53d5bc8c4237251ff3d7cd5/qcreport.html
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/2f5358f82bf0f8bf8deaf19ac88e9d52b19c506865df040f01bcc2d3acfef0bb/snps-excluded.txt
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/84eb5b0390ed1a7c1058f469c800320e63dd7e72580814648b9eb4e9b5b01d94/typed-only.txt
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/6ad6ffe32efc4fa4b937b65cce14f24d3ca77b3200dc0cbc6daf0eb73f2f4732/results.md5
cd ${b_imputed_dir}
wget https://imputation.biodatacatalyst.nhlbi.nih.gov/share/results/23446580e184934875b7027813273bf2a9b869e25042cec31524934ca03bc7bb/chr_1.zip

## Unzipping files
#!/bin/bash
#SBATCH --job-name=unzip_chr#b
#SBATCH --output=/home/rmb218/hpc-work/jobs/imputation/downloads/out/unzip_chr21b.out
#SBATCH --time=01:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
ENCRYPTED_ZIP_FILE="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/chr_21.zip"
OUTPUT_DIR="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/imputed/"
# Define the password for the encrypted zip file
PASSWORD="ToIgxX4$Y6tVWz"
# Unzip the encrypted .zip file
unzip -P "$PASSWORD" "$ENCRYPTED_ZIP_FILE" -d "$OUTPUT_DIR"
# Print completion message
echo "Unzipping completed"
