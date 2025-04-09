
# BiocManager::install("infercnv")
library(infercnv)
library(tidyverse)
library(ggplot2)
options(stringsAsFactors = F)
#  
epi_sce <-  subset(epi, downsample=1000)  
#  
Idents(epi_sce) = paste0("C", Idents(epi_sce))
table(Idents(epi_sce))
sum(table(Idents(epi_sce)))
epi_name <- colnames(epi_sce)

## 
Endothelial <-  sce.all.int[, sce.all.int$celltype == c( 'Endothelial' )]
endo_name  <-  sample(row.names(Endothelial@meta.data),2000) 
Mast <-  sce.all.int[, sce.all.int$celltype == c( 'Mast' )]
Mast_name <-  sample(row.names(Mast@meta.data),2000)  

all_cell <- c(epi_name,endo_name,Mast_name)

### 
sce_cnv <- subset(sce.all.int,cell = all_cell)

# 
sce_cnv$cellname <- colnames(sce_cnv)
epi_sce$cellname <- colnames(epi_sce)

sce_cnv@meta.data <- sce_cnv@meta.data %>% 
  mutate(celltype2 = case_when(
    cellname %in% epi_sce$cellname ~ epi_sce$celltype[match(cellname, epi_sce$cellname)],
    celltype == "Mast" ~ "Mast",
    celltype == "Endothelial" ~ "Endothelial",
    TRUE ~ NA_character_
  ))

#### 
library(AnnoProbe) # devtools::install_github("jmzeng1314/AnnoProbe")  
geneInfor <- annoGene(rownames(sce_cnv), "SYMBOL",'human')  #(GRCh38)
head(geneInfor)
#
geneInfor <- geneInfor[with(geneInfor,order(chr, start)),c(1,4:6)]
#
geneInfor <- geneInfor[!duplicated(geneInfor[,1]),] 
#
write.table(geneInfor ,file = "geneFile.txt",sep = '\t',quote = F,col.names = F,row.names = F)

### 
dat <- GetAssayData(sce_cnv,layer ='counts',assay='RNA')
dat=dat[match(geneInfor[,1], rownames(dat)),]
dim(dat)
identical(rownames(dat),geneInfor$SYMBOL)
##
write.table(dat,file = 'expFile.txt',row.names = T ,sep = '\t',quote = FALSE)

##
groupinfo <-  data.frame(v1=rownames(sce_cnv@meta.data),v2= sce_cnv@meta.data$celltype2 )
groupinfo$v1 <- gsub("-",".",groupinfo$v1)  
##
groupinfo$v2[sample(which(groupinfo$v2 == "Mast"), 1000)] <- "Ref_mast"
groupinfo$v2[sample(which(groupinfo$v2 == "Endothelial"), 1000)] <- "Ref_endo"
#
write.table(groupinfo,file = "groupFiles.txt",           
            sep = '\t',quote = F,col.names = F,row.names = F)


# 
infercnv_obj <-  CreateInfercnvObject(raw_counts_matrix='expFile.txt',
                                      annotations_file='groupFiles.txt',
                                      delim="\t",
                                      gene_order_file= 'geneFile.txt',
                                      ref_group_names=c("Mast","Endothelial")) 
#
dir.create("SNV")
infercnv_sc <-  infercnv::run(infercnv_obj,
                              cutoff=0.1, 
                              out_dir="SNV/", 
                              cluster_by_groups=TRUE, 
                              denoise=T,
                              HMM=T,
                              write_expr_matrix=T,
                              output_format = "pdf")

##cnv_score
data <- read.table("SNV/infercnv.observations.txt", header=T)
expr <- data %>% as.matrix()
expr.scale <- scale(t(expr))
tmp1 <- sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2 <- apply(expr.scale, 2, max) - apply(expr.scale,2,min)
expr_1 <- t(2*sweep(tmp1, 2, tmp2, "/")-1)
cnv_score <- as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score) <- "cnv_score"
cnv_score <- rownames_to_column(cnv_score, var='cell')
cnv_score$cell <- gsub("\\.","-",cnv_score$cell)  
##
colnames(groupinfo)[1] <- "cell"
groupinfo$cell <- gsub("\\.","-",groupinfo$cell)
#
meta <- merge(cnv_score,groupinfo,by="cell", all=F)
meta <- rename(meta, 'celltype' = 'V2')
##
ggplot2::ggplot(meta,aes(x=celltype,y=cnv_score))+
  geom_violin(aes(fill=celltype),cex=1.2)+  
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none')
##
ggplot(meta, aes(x=celltype  , y=cnv_score, fill=celltype  )) +
  geom_boxplot()

