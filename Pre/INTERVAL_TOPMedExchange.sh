mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/clean
mv /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/chromosomes/batch_a_chr*.vcf.gz /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/clean/.
mv /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/chromosomes/batch_a_chr*.vcf.gz.csi /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/clean/.
rm -R /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_a/chromosomes/


mkdir /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/clean
mv /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/chromosomes/batch_b_chr*.vcf.gz /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/clean/.
mv /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/chromosomes/batch_b_chr*.vcf.gz.csi /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/clean/.
rm -R /home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/chromosomes/


## Files for imputation need to be uploaded in groups
## This needs to be done for batch_a and batch_b
#  Python virtual environment for TOPMed
#virtualenv /home/rmb218/rds/hpc-work/python_env/topmed
source /home/rmb218/rds/hpc-work/python_env/topmed/bin/activate

import requests
##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoicHJMZGJTUTlXY0dkRERQTzlYOGdtRG8xdFE2dW1BIiwiZXhwaXJlIjoxNzIyNjkyODA4MjU0LCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.MDR9AbMNvnGPvswqIQusYwqrp-3QymtVRbsgxVIO9BM'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'batch_b_chr1',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg38',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job
file_paths = ['/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/data/batch_b/clean/batch_b_chr{}.vcf.gz'.format(i) for i in range(1,2)]
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


