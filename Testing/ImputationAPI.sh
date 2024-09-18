###################################
##  INTERVAL TOPMed imputation ####
###################################
## Beginning work to impute INTERVAL data with TOPMed r3
##  Run an interactive node
sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
module load plink/2.00-alpha

###################################
##  Directories                ####
###################################
##  Data locations
geno_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
working_dir="/home/rmb218/rds/hpc-work/interval/imputation/"

###################################
##  Variant and sample lists   ####
###################################
##  This info comes from AB and the README.txt file in the related directory
##  List of sampleIDs that comes from interval_qced_24.8.18 (42,396 samples, 743,989 variants)
##  List of variants from more stringent QC merged_imputation (43,059 samples, 655,045 variants)
plink2 --bfile ${geno_data_dir}interval_qced_24.8.18 --write-samples --out ${working_dir}qc_files/sample_ids
plink2 --bfile ${geno_data_dir}merged_imputation --write-snplist --out ${working_dir}qc_files/variant_ids
##  Check numbers
wc -l ${working_dir}/qc_files/sample_ids.id ## 42397 (file has header)
wc -l ${working_dir}/qc_files/variant_ids.snplist ## 665045

###################################
##  Pre-imputation files       ####
###################################
##  Make new files using sample_ids.id and variant_ids.snplist
plink2 --bfile ${geno_data_dir}interval_qced_24.8.18 --keep ${working_dir}qc_files/sample_ids.id --extract ${working_dir}qc_files/variant_ids.snplist --make-bed --out ${working_dir}intermediate_files/pre_impute

###################################
##  Pre-imputation files       ####
###################################
##  Make 2 sets of test files
##  The below pulls out 1000 lines at random, and uses a random seed to minimise chances of overlap
shuf --random-source=/dev/random -n 1000 ${working_dir}qc_files/sample_ids.id > ${working_dir}testing_files/random_sample1.txt
shuf --random-source=/dev/random -n 1000 ${working_dir}qc_files/sample_ids.id > ${working_dir}testing_files/random_sample2.txt
##  Extract the randomised individuals from plink files, make new plink files
plink2 --bfile ${working_dir}intermediate_files/pre_impute --keep ${working_dir}testing_files/random_sample1.txt --make-bed --out ${working_dir}testing_files/random_sample1
plink2 --bfile ${working_dir}intermediate_files/pre_impute --keep ${working_dir}testing_files/random_sample2.txt --make-bed --out ${working_dir}testing_files/random_sample2

###################################
##  Prepaere imputation files  ####
###################################
plink2 --bfile ${working_dir}testing_files/random_sample1 --keep-allele-order --recode vcf id-paste=iid bgz --snps-only just-acgt --out ${working_dir}testing_files/random_sample1
plink2 --bfile ${working_dir}testing_files/random_sample2 --keep-allele-order --recode vcf id-paste=iid bgz --snps-only just-acgt --out ${working_dir}testing_files/random_sample2

## TOPMed needs individual chromosomes
## Run this as a job
random_sample_chromosome_vcfs.sh
#!/bin/bash
#SBATCH --job-name=random_sample_vcfs
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputation_jobs/out/random_sample_vcfs_%j.out
#SBATCH --time=01:00:00
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
geno_data_dir="/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/"
working_dir="/home/rmb218/rds/hpc-work/interval/imputation/"

plink2 \
  --bfile ${working_dir}testing_files/random_sample1 \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --recode vcf id-paste=iid bgz \
  --snps-only just-acgt \
  --out ${working_dir}testing_files/randomsample1_chromosomevcf/rs1_chr${SLURM_ARRAY_TASK_ID}

#bcftools index ${working_dir}testing_files/randomsample1_chromosomevcf/rs1_chr${SLURM_ARRAY_TASK_ID}.vcf.gz
plink2 \
  --bfile ${working_dir}testing_files/random_sample2 \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --recode vcf id-paste=iid bgz \
  --snps-only just-acgt \
  --out ${working_dir}testing_files/randomsample2_chromosomevcf/rs2_chr${SLURM_ARRAY_TASK_ID}

#bcftools index ${working_dir}testing_files/randomsample2_chromosomevcf/rs2_chr${SLURM_ARRAY_TASK_ID}.vcf.gz

###################################
##  Access TOPMed API          ####
###################################
##  Follow guidence from: https://topmedimpute.readthedocs.io/en/latest/api/
python_path="/home/rmb218/rds/hpc-work/python_env/"
module load python/3.8
#virtualenv ${python_path}rb_python
source ${python_path}rb_python/bin/activate
#deactivate
#pip install requests
###################################
##  Random sample 1            ####
###################################
##  Python TOPMed API access
python
import requests

##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoiSk1IektKZ0k3MDBTS3FvMU5iNzNTZVYwVTNFTjc5IiwiZXhwaXJlIjoxNzE3MDc3MDgwMjAxLCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.cWUytpgbX2OTGca4cmoOhzcAK_DgFzt9PgFs2ebCAlc'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'random_sample_1',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg19',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job. Multiple files per chromosome
file_paths = ['/home/rmb218/rds/hpc-work/interval/imputation/testing_files/randomsample1_chromosomevcf/rs1_chr{}.vcf.gz'.format(i) for i in range(1,23)]
files = [('files', open(file_path, 'rb')) for file_path in file_paths]

endpoint = "/jobs/submit/imputationserver"
resp = requests.post(base + endpoint, files=files, data=data, headers=headers)

output = resp.json()

if resp.status_code != 200:
  print(output['message'])
  raise Exception('POST {} {}'.format(endpoint, resp.status_code))
else:
    # print message
    print(output['message'])
    print(output['id'])

###################################
##  RS1 download and format    ####
###################################

rs1_download.sh
#!/bin/bash
#SBATCH --job-name=rs1_download
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputation_jobs/out/rs1_download.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU

cd /home/rmb218/rds/hpc-work/interval/imputation/testing_files/rs1_imputed/
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1252372/7c67397123b80383ca08921c62849f0ea78dcbcff5ed828b3a0c024aa43b0dd6 | bash

##  Unzip all files
unzip -P '<nl3yTEXeeiQ9O' '/home/rmb218/rds/hpc-work/interval/imputation/testing_files/rs1_imputed/*.zip'

##  Combine files into one vcf, create plink file
rs1_combine_post_impute.sh
#!/bin/bash
#SBATCH --job-name=rs1_combine_post_impute
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputation_jobs/out/rs1_combine_post_impute.out
#SBATCH --time=10:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
data_dir="/home/rmb218/rds/hpc-work/interval/imputation/testing_files/rs1_imputed/"

bcftools concat ${data_dir}chr{1..22}.dose.vcf.gz \
    -O z \
    -o ${data_dir}merged_rs1_imputed_dose.vcf.gz 
bcftools index ${data_dir}/merged_rs1_imputed_dose.vcf.gz 

###################################
##  Random sample 2            ####
###################################
module load python/3.8
python_path="/home/rmb218/rds/hpc-work/python_env/"
source ${python_path}rb_python/bin/activate
##  Python TOPMed API access
import requests

##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoiSk1IektKZ0k3MDBTS3FvMU5iNzNTZVYwVTNFTjc5IiwiZXhwaXJlIjoxNzE3MDc3MDgwMjAxLCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.cWUytpgbX2OTGca4cmoOhzcAK_DgFzt9PgFs2ebCAlc'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'random_sample_2',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg19',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job. Multiple files per chromosome
file_paths = ['/home/rmb218/rds/hpc-work/interval/imputation/testing_files/randomsample2_chromosomevcf/rs2_chr{}.vcf.gz'.format(i) for i in range(1,23)]
files = [('files', open(file_path, 'rb')) for file_path in file_paths]

endpoint = "/jobs/submit/imputationserver"
resp = requests.post(base + endpoint, files=files, data=data, headers=headers)

output = resp.json()

if resp.status_code != 200:
  print(output['message'])
  raise Exception('POST {} {}'.format(endpoint, resp.status_code))
else:
    # print message
    print(output['message'])
    print(output['id'])

###################################
##  RS2 download and format    ####
###################################
rs2_download.sh
#!/bin/bash
#SBATCH --job-name=rs2_download
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/invterval_jobs/imputation_jobs/out/rs2_download.out
#SBATCH --time=02:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
cd /home/rmb218/rds/hpc-work/interval/imputation/testing_files/rs2_imputed/
curl -sL https://imputation.biodatacatalyst.nhlbi.nih.gov/get/1257168/020d7c17bdc7daa3f2ba72a8675d7e40f5f6f0f92d6d68de3463207c1e57ce38 | bash

##  Unzip all files
##  Make sure to unzip separately from download
unzip -P 'jb4CG1rPdHp9BK' '/home/rmb218/rds/hpc-work/interval/imputation/testing_files/rs2_imputed/*.zip'

##  Combine files into one vcf, create plink file
rs2_combined_post_impute.sh
#!/bin/bash
#SBATCH --job-name=rs1_vcfmerge
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputatiobn_jobs/out/rs2_combine_post_impute.out
#SBATCH --time=10:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
data_dir="/home/rmb218/rds/hpc-work/interval/imputation/testing_files/rs2_imputed/"

bcftools concat ${data_dir}chr{1..22}.dose.vcf.gz \
    -O z \
    -o ${data_dir}merged_rs2_imputed_dose.vcf.gz 
bcftools index ${data_dir}/merged_rs2_imputed_dose.vcf.gz 

###################################
##  Combine RS1 and RS2        ####
###################################
rs1_combine_rs2.sh
#!/bin/bash
#SBATCH --job-name=sort_random_samples
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputation_jobs/out/sort_random_samples.out
#SBATCH --time=10:00:00
#SBATCH --partition=icelake
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmb218@cam.ac.uk
#SBATCH --account=BUTTERWORTH-SL3-CPU
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
data_dir="/home/rds/project/rds-zuZwCZMsS0w/carriage_gwas/first_pass/results/"
out_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/imputation/testing/"
bcftools sort ${data_dir}rs1_imputed/merged_rs1_imputed_dose.vcf.gz -O z -o ${out_dir}rs1_sorted.vcf.gz
bcftools index ${out_dir}rs1_sorted.vcf.gz
bcftools sort ${data_dir}rs2_imputed/merged_rs2_imputed_dose.vcf.gz -O z -o ${out_dir}rs2_sorted.vcf.gz
bcftools index ${out_dir}rs2_sorted.vcf.gz

## Merge
bcftools merge ${out_dir}rs1_sorted.vcf.gz ${out_dir}rs2_sorted.vcf.gz -O z -o ${out_dir} interval_merged_rs1rs2