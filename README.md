# The following is the workflow for INTERVAL

## Structure

All work for imputation of INTERVAL individuals (~42K) is organised in the following directory structure:
```
└── 📁interval_data

    └── 📁Post
        └── INTERVAL_ChromosomeFiltering.sh
        └── INTERVAL_ImputationProcessing.sh
        └── INTERVAL_Merge.sh
        └── INTERVALDownloadTemplate.sh
        └── snpStats.sh
    └── 📁Pre
        └── INTERVAL_Imputation.sh
        └── INTERVAL_Liftover.sh
        └── INTERVAL_TOPMedExchange.sh
    └── 📁QC
        └── .DS_Store
        └── INTERVAL_PostImputationQC.sh
        └── INTERVAL_QC_PCA.sh
        └── INTERVAL_QC_PRE.sh
        └── INTERVALChromosomalQC.R
    └── 📁Testing
        └── BCFtoolsLiftover.sh
        └── Chr22Merge.sh
        └── ImputationAPI.sh
        └── ImputationCARRIAGE.sh
        └── INTERVAL_CHR22.sh
        └── INTERVAL_Imputation_CHR22.sh
        └── INTERVALImputation.sh
        └── test_INTERVALLiftOver.sh
        └── TestingLiftover.sh
        └── TOPMedSubmission.txt
    └── .DS_Store
    └── .gitattributes
    └── INTERVAL_Organisation.sh
    └── README.md
```

## Order of operations

1. INTERVAL_Liftover.sh

This uses data from `/home/rmb218/rds/rds-jmmh2-post_qc_data/interval/genotype/affy_ukbiobank_array/genotyped/` which contains the merged_imputation files.

Script which details the manual liftover process, matching SNP IDs and positions, updates from GRCh37 to GRCh38.
From previous QC performed in 2018, list 655,045 variants and 42,396 individuals.
Trims to recommended number of SNPs, and individuals.
Create pre_impute_interval_hg38.bed/bim/fam files.

The reference file here comes form Nick Gleadall.
Comparisons between manual liftover and using inbuilt TOPMed lisftover did not show any significant differences.

Information from Adam Butterworth and Tao Jiang outines that MAF filtering was not performed on this whole cohort, as it was decided that rarer typed SNPs would help with the imputation of rare variants. For additional information see the readme file in the relevant directory.

2. INTERVAL_Imputation.sh

Takes in lifted over files, aligns to reference .fasta to create interva_refhg38.
Splits INTERVAL into two batches (a and b), random process.
Splits batches into 22 chromosomes, makes vcf files, and indexes them.

3. INTERVAL_TOPMedExchange.sh

This submits each chromosome file to the TOPMed server for imputation.

4. INTERVALDownloadTemplate.sh

Download files from the TOPMed imputation server.

5. INTERVAL_ImputationProcessing.sh

Unzips downloaded files with password provided by TOPMed.
Indexes .dose.vcf.gz files.

6. INTERVAL_PostImputationQC.sh

Makes initial QC files from each batch.
Extracts information from the info files for each chromosome in each batch.
Creates plots of MAF and R2 outliers, also creates a file of identified SNP outliers.
Uses R script to make plots for all chromosomes.

6.1. snpStats.sh

This looks at comparing the previous imputation of the INTERVAL data on hg19 with my current data imputed on hg38

7. INTERVAL_ChromosomeFiltering.sh

Filters batch a and batch b files to maf >= 0.005 and R2 >= 0.7

8. INTERVAL_Merge.sh

Merge all batch_a and batch_b files per chromosome. 
Output QC files per chromosome.

9. INTERVAL_Organisation.sh

Organise all files related to the INTERVAL cohort created in this process. 
This is really just me organising the output, deleting temporary files, and zipping, ready for other people to access.

10. INTERVAL_PCA.sh

Details the process for an ancestry PCA analysis with 1KG as the reference population (hg38). This is done for the >42K samples.