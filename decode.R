
rm(list=ls())
library(TwoSampleMR)
library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)

#### 
decode <- fread("../../../proteomics/decode_result/decode_mendelian_5000_0.1.csv")
stomach <- fread("../finngen_R11_C3_STOMACH_EXALLC.gz")
write.csv(stomach,"stomach.csv")

#### mr
exposure <- read_exposure_data(filename = "../../../proteomics/decode_result/decode_mendelian_5000_0.1.csv",
                               clump = F,sep = ",",phenotype_col = "gene"
                               ,snp_col = "rsids"
                               ,beta_col = "Beta"
                               ,se_col = "SE"
                               ,eaf_col = "effectAlleleFreq"
                               ,effect_allele_col = "effectAllele"
                               ,other_allele_col = "otherAllele"
                               ,pval_col = "Pval"
                               ,samplesize_col = "N"
)
outcome <- read_outcome_data(filename = "stomach.csv"
                             ,snps = exposure$SNP
                             ,sep = ","
                             ,snp_col = "rsids"
                             ,beta_col = "beta"
                             ,se_col = "sebeta"
                             ,eaf_col = "af_alt"
                             ,effect_allele_col = "alt"
                             ,other_allele_col = "ref"
                             ,pval_col = "pval"
)
exposure[3449,2] <- "T"
exposure[6125,2] <- "T"
dat <- harmonise_data(exposure_dat = exposure,outcome_dat = outcome)

####  
cl <- makeCluster(8)
registerDoParallel(cl)
gene <- unique(dat$exposure)
res <- foreach(gene2 = gene,.combine = rbind,.packages = c("TwoSampleMR","tidyverse"),.verbose = T) %dopar% {
  element <- dat %>% filter(dat$exposure %in% gene2)
  mr(element)
}
stopCluster(cl)
res2 <- res %>% filter(method=="Wald ratio"|method=="Inverse variance weighted")
res2 <- res2 %>% filter(pval < 0.05)
write.csv(test,"res_final.csv")
