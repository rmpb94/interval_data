main_dir="/rds/project/rds-zuZwCZMsS0w/carriage_gwas/"

mkdir ${main_dir}interval_data
mkdir ${main_dir}interval_data/bed
mv ${main_dir}interval/data/pre_impute ${main_dir}interval_data
mv ${main_dir}interval/data/imputed ${main_dir}interval_data
