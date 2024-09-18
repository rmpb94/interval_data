## INTERVAL QC ####

## Libraries ####
library(tidyverse)
library(patchwork)
library(MetBrewer)

pal = met.brewer("Thomas")

## Data ####
batch_a_file_list = list.files("batch_a/", pattern = ".txt")
batch_b_file_list = list.files("batch_b/", pattern = ".txt")

maf_outlier_snps = c()
r2_outlier_snps = c()
total_snps_all = c()

for(f in 1:length(batch_a_file_list)){
  batch_a = read_table(paste0("batch_a/", batch_a_file_list[f]), 
                    col_names = c("chr", "bp", "ref", "alt", "a_af", "a_maf", "a_avg_cs", "a_r2", "a_er2", "a_imputed", "a_typed"))
  batch_b = read_table(paste0("batch_b/", batch_b_file_list[f]), 
                       col_names = c("chr", "bp", "ref", "alt", "b_af", "b_maf", "b_avg_cs", "b_r2", "b_er2", "b_imputed", "b_typed"))
  
  ## Pre-filtering
  batch_a = subset(batch_a, batch_a$a_maf >= 0.005)
  batch_a = subset(batch_a, batch_a$a_r2 >= 0.7)
  batch_b = subset(batch_b, batch_b$b_maf >= 0.005)
  batch_b = subset(batch_b, batch_b$b_r2 >= 0.7)
 
  print(dim(batch_a))
  print(dim(batch_b))
 
  chrom_dat = merge(batch_a, batch_b, by = c("chr", "bp", "ref", "alt"))
  chrom_dat$maf_diff = abs(chrom_dat$a_maf - chrom_dat$b_maf)
  chrom_dat$r2_diff = abs(chrom_dat$a_r2 - chrom_dat$b_r2)
  

  maf_dat_1 = subset(chrom_dat, chrom_dat$maf_diff < 0.01)
  maf_dat_2 = subset(chrom_dat, chrom_dat$maf_diff > 0.01)

  r2_dat_1 = subset(chrom_dat, chrom_dat$r2_diff < 0.05)
  r2_dat_2 = subset(chrom_dat, chrom_dat$r2_diff > 0.05)

  maf_cor_plot = 
    ggplot() + theme_light() +
    geom_point(data = maf_dat_1, aes(x = a_maf, y = b_maf), colour = pal[8], alpha = 0.5) +
    geom_point(data = maf_dat_2, aes(x = a_maf, y = b_maf), colour = pal[2], alpha = 0.5) +
    geom_abline(colour = pal[7], alpha = 0.7) +
    labs(x = "Batch A MAF", y = "Batch B MAF", title = paste0("Chromosome ", f))
  
  r2_cor_plot = 
    ggplot() + theme_light() +
    geom_point(data = r2_dat_1, aes(x = a_r2, y = b_r2), alpha = 0.5, colour = pal[8]) +
    geom_point(data = r2_dat_2, aes(x = a_r2, y = b_r2), colour = pal[2], alpha = 0.5) +
    geom_abline(colour = pal[7], alpha = 0.7) +
    labs(x = "Batch A R2", y = "Batch B R2", title = paste0("Chromosome ", f))

chrom_dat$colour = ifelse(chrom_dat$maf_diff >= 0.01 & chrom_dat$r2_diff >= 0.05, "R2 and MAF",
                    ifelse(chrom_dat$maf_diff >= 0.01 & chrom_dat$r2_diff < 0.05, "MAF",
                           ifelse(chrom_dat$maf_diff < 0.01 & chrom_dat$r2_diff >= 0.05, "R2",
                                  ifelse(chrom_dat$maf_dif < 0.01 & chrom_dat$r2_diff < 0.05, "Pass", NA))))

total_snps = data.frame(batch_a_file_list[f], dim(chrom_dat)[1])

dat1 = subset(chrom_dat, chrom_dat$colour == "Pass")
dat2 = subset(chrom_dat, chrom_dat$colour != "Pass")

  maf_chromosome_plot = 
    ggplot() + theme_light() +
      geom_point(data = dat1, aes(x = bp, y = maf_diff), colour = "grey80") +
      geom_point(data = dat2, aes(x = bp, y = maf_diff, colour = colour), size = 2) +
      scale_colour_manual(values = c(pal[1], pal[3], pal[8])) +
      labs(colour = " ", x = "Base pair position (bp)", y = "MAF Difference")

all_plots = (maf_cor_plot + r2_cor_plot) / maf_chromosome_plot

  ## Save plots
  jpeg(paste0("qc/post_impute/graphs/filtered_chromosome_", f, ".jpeg"), units = "in", height = 6, width = 10, res = 200)
  print(all_plots)
  dev.off()

  print("Plot saved")

  temp_maf_outliers = subset(chrom_dat, chrom_dat$maf_diff >= 0.01)
  temp_r2_outliers = subset(chrom_dat, chrom_dat$r2_diff >= 0.05)
  
  maf_outlier_snps = rbind(maf_outlier_snps, temp_maf_outliers)
  r2_outlier_snps = rbind(r2_outlier_snps, temp_r2_outliers)

  total_snps_all = rbind(total_snps_all, total_snps)
}
write.table(maf_outlier_snps, "qc/post_impute/high_r2_maf_outlier_snps.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(r2_outlier_snps, "qc/post_impute/high_r2_r2_outlier_snps.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(total_snps_all, "qc/post_impute/total_snps.txt", col.names = T, row.names = F, quote = F, sep = "\t")