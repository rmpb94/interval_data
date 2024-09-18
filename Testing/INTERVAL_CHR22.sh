## QC of chr22, but with maf filtering before imputation
## How does maf filter affect imputation quality? How does it affect the MAF distribution of the imputed data?
mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/chr_22/
mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/chr_22/input

working_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/chr_22/"

# Move chromosome 22 files into directory
cp /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/clean/batch_a_chr22.vcf.gz ${working_dir}input/
cp /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/clean/batch_b_chr22.vcf.gz ${working_dir}input/

sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR

module load plink/2.00-alpha
module load bcftools
working_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/chr_22/"

## Make plink bed files
plink2 --vcf ${working_dir}input/batch_a_chr22.vcf.gz --make-bed --out ${working_dir}input/batch_a_chr22
plink2 --vcf ${working_dir}input/batch_b_chr22.vcf.gz --make-bed --out ${working_dir}input/batch_b_chr22

## Filter for MAF > 0.01
plink2 --bfile ${working_dir}input/batch_a_chr22 --maf 0.01 --make-bed --out ${working_dir}input/batch_a_chr22_maf
plink2 --bfile ${working_dir}input/batch_b_chr22 --maf 0.01 --make-bed --out ${working_dir}input/batch_b_chr22_maf

plink2 --bfile ${working_dir}input/batch_a_chr22_maf --recode vcf id-paste=iid bgz --out ${working_dir}new_batch_a_chr22
plink2 --bfile ${working_dir}input/batch_b_chr22_maf --recode vcf id-paste=iid bgz --out ${working_dir}new_batch_b_chr22

bcftools index ${working_dir}new_batch_a_chr22.vcf.gz
bcftools index ${working_dir}new_batch_b_chr22.vcf.gz

echo "22 chr22" > ${working_dir}input/chr_22_name_conv.txt
bcftools annotate --rename-chrs ${working_dir}input/chr_22_name_conv.txt ${working_dir}new_batch_b_chr22.vcf.gz -o ${working_dir}new_batch_b_chr22.vcf.gz
bcftools annotate --rename-chrs ${working_dir}input/chr_22_name_conv.txt ${working_dir}new_batch_a_chr22.vcf.gz -o ${working_dir}new_batch_a_chr22.vcf
bgzip ${working_dir}new_batch_a_chr22.vcf
bgzip ${working_dir}new_batch_b_chr22.vcf

## Submit to imputation server
source /home/rmb218/rds/hpc-work/python_env/topmed/bin/activate
import requests
##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoidG8xS2F4d0k3RDZYeU45QmJ4cWNoWGpuTkRSQUczIiwiZXhwaXJlIjoxNzI4OTAwODU2MTM5LCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.daOixkQMBJ3Yf7Wh1utE6awmpQeafMjgKeE6dN3feDU'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'batch_b_chr22',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg38',
  'phasing': 'eagle',
  'r2Filter': 0.3,
  'meta': 'yes',
}
##  Submit job
file_paths = ['/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/chr_22/new_batch_b_chr{}.vcf.gz'.format(i) for i in range(22,23)]
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

