#################################################
## Test TOPMed liftover                        ##
#################################################
## I want to know how manual lift-over performs against the liftover inbuilt into TOPMed
cd /home/rmb218/rds/hpc-work/interval/imputation/liftover
module load plink/2.00-alpha
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy

rb_dir="/home/rmb218/rds/hpc-work/interval/imputation/liftover/"
info_dir="/home/rmb218/rds/hpc-work/interval/imputation/"

## Use pre-liftover file, but just trim to chromosome 22
plink2 \
  --pfile ${rb_dir}rb_pre_impute_interval \
  --chr 22 \
  --make-bed \
  --out ${rb_dir}testing/batch_a_22_hg19

## Split cohort into batch a
plink2 \
    --bfile ${rb_dir}testing/batch_a_22_hg19 \
    --keep ${info_dir}info_files/batch_a_individuals.txt \
    --make-bed \
    --out ${rb_dir}testing/batch_a_22_hg19

## Test with udating ref and alt alleles
cd /home/rmb218/rds/hpc-work/ref_files/
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
## Try reference fasta for just chr22
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr22.fa.gz

cd ${rb_dir}

plink2 \
  --bfile ${rb_dir}testing/batch_a_22_hg19 \
  --fa /home/rmb218/rds/hpc-work/ref_files/chr22.fa.gz \
  --ref-from-fa 'force' \
  --make-bed \
  --out ${rb_dir}testing/batch_a_22_hg19

plink2 \
    --bfile ${rb_dir}testing/batch_a_22_hg19 \
    --recode vcf id-paste=iid bgz \
    --snps-only just-acgt \
    --out ${rb_dir}testing/batch_a_22_hg19

bcftools index ${rb_dir}testing/batch_a_22_hg19.vcf.gz

## Submit job to TOPMed
source /home/rmb218/rds/hpc-work/python_env/topmed/bin/activate
import requests
##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoiWlVMc0JDRTZrVmNrY0xDYXV4UWJjdkhDRXdyNGo4IiwiZXhwaXJlIjoxNzIwMDAwODA3MzA2LCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.ZFhGOkNbkN-KlhhfoRc_cO6t2quj1fbkaack2RksF8g'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'chr22_test_hg19_ref_update',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg19',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job. Multiple files per chromosome
file_paths = ['/home/rmb218/rds/hpc-work/interval/imputation/liftover/testing/batch_a_{}_hg19.vcf.gz'.format(i) for i in range(22,23)]
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
