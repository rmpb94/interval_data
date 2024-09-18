##########################################
##  INTERVAL chr22                      ## 
##########################################
##  Testing imputation with chromosome 22

##  Modules
module load plink/2.00-alpha
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy

cd /home/rmb218/hpc-work/ref_files
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

##  Data    ####
##  Lifted INTERVAL data
interval_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/"
info_dir="/home/rmb218/rds/rds-jmmh2-projects/carriage_gwas/interval/info/"

##  Some initial tests
plink2 \
    --pfile ${interval_dir}pre_impute_interval_hg38 \
    --write-snplist \
    --write-samples \
    --out ${info_dir}interval_hg38
wc -l ${info_dir}interval_hg38.snplist ## 646,726 variants in the hg38 liftover file
## 42,396 samples
## pvar file shows positional info

##  Split cohort in half
plink2 \
    --bfile ${interval_dir}pre_impute_interval_hg38 \
    --write-samples \
    --out ${info_dir}/full_sample_list

##  Generate first random number list
shuf -i 1-42396 -n 21198 > ${info_dir}interval_random_split.txt
##  Use numbers to extract individuals
awk 'NR==FNR {lines[$1]; next} FNR in lines' ${info_dir}/interval_random_split.txt ${info_dir}full_sample_list.id > ${info_dir}batch_a_individuals.txt


plink2 \
    --pfile ${rb_dir}pre_impute_interval_hg38 \
    --keep ${info_dir}info_files/batch_a_individuals.txt \
    --make-bed \
    --out ${rb_dir}testing/pre_impute_interval_batch_a
##  Prepare imputation file

##  Testing re-alignment of ref/alt alleles with hg38 .fa file
sintr -A BUTTERWORTH-SL3-CPU -p icelake -N2 -n2 -t 1:0:0 --qos=INTR
module load plink/2.00-alpha
module load bcftools/1.14/gcc/intel-oneapi-mkl/xrcih6dy
rb_dir="/home/rmb218/rds/hpc-work/interval/imputation/liftover/"

plink2 \
  --bfile ${rb_dir}testing/pre_impute_interval_batch_a \
  --chr 22 \
  --make-bed \
  --out ${rb_dir}testing/batch_a_22_ref

plink2 \
  --bfile ${rb_dir}testing/batch_a_22_ref \
  --fa /home/rmb218/rds/hpc-work/ref_files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
  --ref-from-fa 'force' \
  --make-bed \
  --out ${rb_dir}testing/batch_22_a_refhg38
 
plink2 \
  --bfile ${rb_dir}testing/batch_22_a_refhg38 \
  --recode vcf id-paste=iid bgz \
  --snps-only just-acgt \
  --out ${rb_dir}testing/batch_a_22

bcftools index ${rb_dir}testing/batch_a_22.vcf.gz

##  Chromosomes must be coded chr#
echo "22 chr22" > ${rb_dir}testing/chr_22_name_conv.txt
bcftools annotate --rename-chrs ${rb_dir}testing/chr_22_name_conv.txt ${rb_dir}testing/batch_a_22.vcf.gz -o ${rb_dir}testing/batch_a_chr22.vcf
bgzip ${rb_dir}testing/batch_a_chr22.vcf
bcftools index ${rb_dir}testing/batch_a_chr22.vcf.gz

##  Python virtual environment for TOPMed
#virtualenv /home/rmb218/rds/hpc-work/python_env/topmed
source /home/rmb218/rds/hpc-work/python_env/topmed/bin/activate

import requests
##   Imputation server URL
base = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
##  NB. Tokens are only valid for 30 days
token = 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoicm1iMjE4QGNhbS5hYy51ayIsImFwaV9oYXNoIjoiWlVMc0JDRTZrVmNrY0xDYXV4UWJjdkhDRXdyNGo4IiwiZXhwaXJlIjoxNzIwMDAwODA3MzA2LCJuYW1lIjoiUm9pc2luIEJvZ2dhbiIsImFwaSI6dHJ1ZSwidXNlcm5hbWUiOiJyb2lzaW5ib2dnYW4ifQ.ZFhGOkNbkN-KlhhfoRc_cO6t2quj1fbkaack2RksF8g'
##  Add API token to headers
headers = {'X-Auth-Token' : token }
data = {
  'job-name': 'chr22_test',
  'refpanel': 'apps@topmed-r3',
  'population': 'all',
  'build': 'hg38',
  'phasing': 'eagle',
  'r2Filter': 0,
  'meta': 'yes',
}
##  Submit job. Multiple files per chromosome
file_paths = ['/home/rmb218/rds/hpc-work/interval/imputation/liftover/testing/batch_a_chr{}.vcf.gz'.format(i) for i in range(22,23)]
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


