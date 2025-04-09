
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

library(stringr)
library(Seurat)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)
library(miloR)
library(SingleCellExperiment)
library(ggbeeswarm)
library(scales)
library(forcats)
library(data.table)
library(Matrix)

register(MulticoreParam(workers = 8, progressbar = TRUE))

sce_immu = subset(sce_immu2,celltype2 %in% c('B','DC','NK','Mast','Mono_c1','Mono_c2' ,'CD8_T','CD8_Treg' ,'CD4/Treg','Mac_c3','Mac_c1_M1',
                                             'Mac_c2_M2','B_naive','Plasma_c1','Plasma_c2_IgG','Plasma_c3','Plasma_c4'))

immu_tumor = sce_immu2
immu_tumor$group <- factor(immu_tumor$group, levels=c("low","high"))
immu_tumor$orig.ident=as.character(immu_tumor$orig.ident)

immu_tumor <- as.SingleCellExperiment(immu_tumor)
# 
sc_milo <- miloR::Milo(immu_tumor) 
#
sc_milo <- miloR::buildGraph(sc_milo, k = 10, d = 20) 
#
sc_milo <- makeNhoods(sc_milo, 
                      prop = 0.05, 
                      k = 10, 
                      d=20, 
                      refined = TRUE,
                      refinement_scheme="graph")
#
sc_milo <- countCells(sc_milo, meta.data = as.data.frame(colData(sc_milo)), 
                      sample="orig.ident") 

#
traj_design <- data.frame(colData(sc_milo))[,c("orig.ident", "group")]
traj_design$orig.ident <- as.factor(traj_design$orig.ident)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident

#
da_results <- testNhoods(sc_milo, 
                         design = ~ group, 
                         design.df = traj_design,
                         fdr.weighting="graph-overlap")

## 
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) 
ggplot(da_results, aes(logFC, -log10(SpatialFDR)))+ 
  geom_point() +
  geom_hline(yintercept = 1)

## 
sc_milo <- buildNhoodGraph(sc_milo)
# 
plotNhoodGraphDA(sc_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) 

da_results <- annotateNhoods(sc_milo, da_results, coldata_col = "celltype2")
ggplot(da_results,aes(celltype2_fraction)) + geom_histogram(bins=50)+theme_classic()

plotDAbeeswarm(da_results, group.by = "celltype2") +
  scale_color_gradient2(low="#3D98D3FF",
                        mid="#CCCCCC",
                        high="#910000",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black'))

plotDAbeeswarm(da_results,group.by = "celltype2", alpha=1)+
  geom_boxplot(aes(group =celltype2),outlier.shape = NA,color =  "#666666")+
  scale_color_gradient2(midpoint=0, low='#3F51B5FF', mid="white",high='#96281BFF', space ="Lab" )+
  theme_bw(base_size=11)+
  geom_hline(yintercept=0, linetype="dashed")
