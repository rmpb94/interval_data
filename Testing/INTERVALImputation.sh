###################################
##  INTERVAL TOPMed imputation ####
###################################
## Beginning work to impute INTERVAL data with TOPMed r3

###################################
##  Directories                ####
###################################
##  Data locations
shared_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/"
pre_data_dir="/home/rmb218/rds/hpc-work/interval/imputation/liftover/"
##  pre_impute_interval_hg38.pgen pre_impute_interval_hg38.pvar pre_impute_interval_hg38.psam

###################################
##  Pre-imputation             ####
###################################
##  TOPMed won't accept all samples at once, need to split the data into 2
##  Do this randomly
working_dir="/home/rmb218/rds/hpc-work/interval/imputation/"
##  Write out a list of sample IDs
plink2 \
    --pfile ${pre_data_dir}pre_impute_interval_hg38 \
    --write-samples \
    --out ${working_dir}info_files/full_sample_list

##  Generate first random number list
shuf -i 1-42396 -n 21198 > ${working_dir}info_files/interval_part1.txt
##  Use numbers to extract individuals
awk 'NR==FNR {lines[$1]; next} FNR in lines' ${working_dir}info_files/interval_part1.txt ${working_dir}info_files/full_sample_list.id > ${working_dir}info_files/batch_a_individuals.txt

##  Make 2 data sets, one for each list of participants
module load plink/2.00-alpha
plink2 \
    --pfile ${pre_data_dir}pre_impute_interval_hg38 \
    --keep ${working_dir}info_files/batch_a_individuals.txt \
    --make-bed \
    --out ${working_dir}intermediate_files/pre_impute_interval_batch_a
plink2 \
    --pfile ${pre_data_dir}pre_impute_interval_hg38 \
    --remove ${working_dir}info_files/batch_a_individuals.txt \
    --make-bed \
    --out ${working_dir}intermediate_files/pre_impute_interval_batch_b

###################################
##  Prepare imputation files   ####
###################################
plink2 \
    --bfile ${working_dir}intermediate_files/pre_impute_interval_batch_a \
    --recode vcf id-paste=iid bgz \
    --snps-only just-acgt \
    --out ${working_dir}intermediate_files/batch_a/interval_a_hg38
plink2 \
    --bfile ${working_dir}intermediate_files/pre_impute_interval_batch_b \
    --recode vcf id-paste=iid bgz \
    --snps-only just-acgt \
    --out ${working_dir}intermediate_files/batch_b/interval_b_hg38

###################################
##  Chromosome files           ####
###################################
interval_a_chromosome_vcf
#!/bin/bash
#SBATCH --job-name=interval_a_chromosome_vcf
#SBATCH --output=/home/rmb218/rds/hpc-work/jobs/interval_jobs/imputation_jobs/out/interval_a_chromosome_vcf-%j.out
#SBATCH --time=00:30:00
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
working_dir="/home/rmb218/rds/hpc-work/interval/imputation/intermediate_files/"
parent_dir="/home/rmb218/rds/hpc-work/interval/imputation/"
plink2 \
  --bfile ${working_dir}pre_impute_interval_batch_a \
  --chr ${SLURM_ARRAY_TASK_ID} \
  --recode vcf id-paste=iid bgz \
  --snps-only just-acgt \
  --out ${working_dir}batch_a/batch_a_${SLURM_ARRAY_TASK_ID}
bcftools index ${working_dir}batch_a/batch_a_${SLURM_ARRAY_TASK_ID}.vcf.gz
##  Recode chromosomes
echo "${SLURM_ARRAY_TASK_ID} chr${SLURM_ARRAY_TASK_ID}" > ${parent_dir}info_files/a_chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt
bcftools annotate --rename-chrs ${parent_dir}info_files/a_chr_${SLURM_ARRAY_TASK_ID}_name_conv.txt ${working_dir}batch_a/batch_a_${SLURM_ARRAY_TASK_ID}.vcf.gz -o ${working_dir}batch_a/batch_a_chr${SLURM_ARRAY_TASK_ID}.vcf
bgzip ${working_dir}batch_a/batch_a_chr${SLURM_ARRAY_TASK_ID}.vcf
bcftools index ${working_dir}batch_a/batch_a_chr${SLURM_ARRAY_TASK_ID}.vcf.gz

###################################
##  Imputation                 ####
###################################
##  Submit to TOPMedr3
python_path="/home/rmb218/rds/hpc-work/python_env/"
module load python/3.8
source ${python_path}topmed/bin/activate
python
import requests
##  Batch A
##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoiWlVMc0JDRTZrVmNrY0xDYXV4UWJjdkhDRXdyNGo4IiwiZXhwaXJlIjoxNzIwMDAwODA3MzA2LCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.ZFhGOkNbkN-KlhhfoRc_cO6t2quj1fbkaack2RksF8g'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'interval_batch_a_1_4',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg38',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job. Multiple files per chromosome
file_paths = ['/home/rmb218/rds/hpc-work/interval/imputation/intermediate_files/batch_a/batch_a_chr{}.vcf.gz'.format(i) for i in range(1,5)]

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

##  Batch B
##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoiS05IMmRFUmhGWFpsMEE0d3RQYk1jT2NJalNhSGE3IiwiZXhwaXJlIjoxNzE5NzM3MjcwNDU5LCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.1q_Dh45h3XR7nR0bfwWRWuSwnpyULkHLyLfqa1lD-Gs'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'interval_batch_b',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg38',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job. Multiple files per chromosome
file_paths = ['/home/rmb218/rds/hpc-work/interval/imputation/intermediate_files/batch_b/batch_b_chr{}.vcf.gz'.format(i) for i in range(1,23)]

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
