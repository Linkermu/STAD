
rm(list=ls())
library(TwoSampleMR)
library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)

ukb <- fread("../../../proteomics/UKB-PPP_result/ukbppp_mendelian_10000_0.001.csv")

exposure <- read_exposure_data(filename = "../../../proteomics/UKB-PPP_result/ukbppp_mendelian_10000_0.001.csv",
                               clump = F,sep = ",",phenotype_col = "gene"
                               ,snp_col = "rsids"
                               ,beta_col = "Beta"
                               ,se_col = "se"
                               ,eaf_col = "eaf"
                               ,effect_allele_col = "alt"
                               ,other_allele_col = "ref"
                               ,pval_col = "Pval"
                               ,samplesize_col = "N"
)
outcome <- read_outcome_data(filename = "../decode/stomach.csv"
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
for (i in seq_along(exposure$effect_allele.exposure)) {
  if(is.na(exposure[i,5])){
    snp <- exposure[i,1]
    if(snp %in% outcome$SNP){
      jianji <- outcome %>% filter(SNP %in% snp) %>% pull(effect_allele.outcome)
      exposure[i,5] <- jianji
    }
  }
}

exposure <- exposure %>% filter(!is.na(effect_allele.exposure))

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
res2$p.adjust <- p.adjust(res2$pval,method = "fdr")
res3 <- res2[res2$pval <0.05,]



