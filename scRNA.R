library(Seurat)
library(patchwork)
library(tidyverse)
library(multtest)
library(metap)
library(tibble)
library(DoubletFinder)
library(data.table)
library(harmony)
library(ggsci)
library(future)
library(ggplot2)
library(Matrix)
library(clustree)
library(cowplot)
library(data.table)

## scRNA2
matrix_data <- fread("GSE206785_scgex.txt.gz",data.table = F)
meta <- fread("GSE206785_metadata.txt.gz")
rownames(meta) <- rownames(matrix_data)
meta$cell <- rownames(meta)
matrix_data <- t(matrix_data)
sce.all=CreateSeuratObject(counts = matrix_data,meta.data = meta)
sce.all$orig.ident <- sce.all$cell
save.image("scrna2.RData")

## scRNA1
samples <- list.files()
sceList = lapply(samples,function(pro){ 
  tmp = Read10X(file.path(".",pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
}) 

##
sce.all_1 <- merge(sceList[[1]],sceList[-1],add.cell.ids = samples)
sceList2 <- list(sce.all,sce.all_1)
sce.all_2 <- merge(sceList2[[1]],sceList2[-1])
sce.all <- JoinLayers(sce.all_2)
sce.all@meta.data <- sce.all@meta.data[,1:3]
sce.all@meta.data$sample <- ifelse(
  str_detect(sce.all@meta.data[["orig.ident"]], "T") | 
    str_detect(sce.all@meta.data[["orig.ident"]], "tumor"),
  "tumor",
  "normal"
)

#######  
sp='human'
setwd("stomach/sc rna")
source('qc.R')
sce.all.filt = basic_qc(sce.all)

###### harmony ######
library(harmony)
if(T){
  dir.create("2-harmony")
  getwd()
  setwd("2-harmony")
  source('harmony.R')
  sce.all.int = run_harmony(sce.all.filt)
}

sce.all.int  
DimPlot(sce.all.int, reduction = "umap",label = T,raster=FALSE) 
DimPlot(sce.all.int, reduction = "umap",raster=FALSE,group.by = "orig.ident",cols = col) + NoLegend()

sce.all.int <- SetIdent(sce.all.int,value = "RNA_snn_res.0.05")
DotPlot(sce.all.int, features = unique(genes_to_check),
        assay='RNA' )  + coord_flip()

genes_to_check = c('EPCAM','KRT19','KRT8',  #epi
                   'COL1A1', 'COL1A2', 'DCN',  #str
                   'CD3D', 'CD3E', 'CD3G', #T 
                   'CD79A', 'IGKC','MZB1',  #b
                   'CD14', 'MS4A7', 'FCGR3A', #mye
                   'VWF','PECAM1','ENG',  #endo
                   'TPSAB1','CPA3','MS4A2', #mast
                   "MKI67","PCNA","TOP2A", # prof
                   "RGS5","ACTA2","FRZB" # Mural Cells
)


celltype=data.frame(ClusterID=0:8,
                    celltype= 0:8) 
celltype[celltype$ClusterID %in% c( 1),2]='Epithelial'
celltype[celltype$ClusterID %in% c( 4),2]='Stromal'
celltype[celltype$ClusterID %in% c(0),2]='T/NK'
celltype[celltype$ClusterID %in% c( 2 ),2]='B'
celltype[celltype$ClusterID %in% c( 3 ),2]='Myeloid'
celltype[celltype$ClusterID %in% c(5 ),2]='Endothelial'
celltype[celltype$ClusterID %in% c( 6 ),2]='Mural'
celltype[celltype$ClusterID %in% c( 7 ),2]='Mast'
celltype[celltype$ClusterID %in% c( 8 ),2]='Cycling'

sce.all.int@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce.all.int@meta.data[which(sce.all.int@meta.data$RNA_snn_res.0.05 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all.int@meta.data$celltype)

DimPlot(sce.all.int, reduction = "umap",label = T,group.by = "celltype",
        cols = c('#E64A35','#4DBBD4' ,'#01A187'  ,'#AF9E85','#3C5588'  ,'#F29F80'  ,
                 '#8491B6','#91D0C1','#BB6239') ) 

DimPlot(sce.all.int, reduction = "umap",group.by = "sample")  


pbmc.markers <- FindAllMarkers(sce.all.int, only.pos = TRUE,
                               min.pct = 0.25,
                               test.use = "wilcox",
                               logfc.threshold = 0.25) 

# top20
top20 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.csv(top20,"top20.csv")

#####epi
epi <- subset(sce.all.int,idents = c(1))
cellname <- colnames(epi)
epi <- subset(sce.all.filt,cells = cellname)
epi_tumor <- epi[,epi$sample %in% "tumor"]
epi_normal <- epi[,epi$sample %in% "normal"]

## 
sce_modify <- RunTSNE(sce_modify, dims = 1:20)
DimPlot(sce_modify, reduction = "tsne",label = T) 
DimPlot(sce_modify, reduction = "tsne",label = T,group.by = "sample") 
saveRDS(sce_modify,"epi.rds")

# tumor
FeaturePlot(sce_modify,features = "CLDN4",reduction = "tsne")
FeaturePlot(sce_modify,features = "CLDN7",reduction = "tsne")
FeaturePlot(sce_modify,features = "TFF3",reduction = "tsne")
# normal
FeaturePlot(sce_modify,features = "MUC5AC",reduction = "tsne")
FeaturePlot(sce_modify,features = "GKN1",reduction = "tsne")
FeaturePlot(sce_modify,features = "PGC",reduction = "tsne")
FeaturePlot(sce_modify,features = "LIPF",reduction = "tsne")

markers <- FindAllMarkers(epi,only.pos = TRUE,
                          min.pct = 0.25,
                          test.use = "wilcox", 
                          logfc.threshold = 0.25) 

# top20
top20_gene <- markers %>% filter(cluster %in% c('C2','C6','C10','C14','C15','C16','C19')) %>% 
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)
gene_name <- unique(top20_gene$gene)
gene_name <- c(gene_name,"TPST2","NUCB1","SLURP1","PCBD1","PDE5A")
gene_name <- unique(gene_name)

## lasso  cox
library(tidyr)
library(reshape2)
library (gplots) 

tb=table(scRNA$sample, scRNA$celltype)
balloonplot(tb)
bar_data <- as.data.frame(tb)
bar_per <- bar_data %>% 
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)

#write.csv(bar_per,file = "celltype_by_group_percent.csv")
col =c("#3176B7","#F78000","#3FA116","#CE2820","#9265C1",
       "#885649","#DD76C5","#BBBE00","#41BED1")

#####
library(ggthemes)
p1 = ggplot(bar_per, aes(x = percent, y = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=col)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))
p1


##
epi_markers <- FindAllMarkers(epi, min.pct = 0.25, 
                              logfc.threshold = 0.25)
epi_markers = epi_markers %>% filter(cluster %in% c('C2','C6','C10','C14','C15','C16','C19'))

colnames(epi_markers)[6] = "celltype"
k = epi_markers$p_val_adj<0.05;table(k)
epi_markers = epi_markers[k,]

##
epi_markers$label <- ifelse(epi_markers$avg_log2FC<0,"sigDown","sigUp")
topgene <- epi_markers %>%
  group_by(celltype) %>%
  top_n(n = 7, wt = avg_log2FC) %>%
  bind_rows(group_by(epi_markers, celltype) %>%
              top_n(n = 7, wt = -avg_log2FC))
##
dfbar = epi_markers %>%
  group_by(celltype) %>%
  summarise(low = round(min(avg_log2FC)-0.5),
            up = round(max(avg_log2FC)+0.5))

##
p1 <- ggplot()+
  geom_col(aes(x = celltype ,y = low),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(aes(x = celltype ,y = up),dfbar,
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(aes(x = celltype, y = avg_log2FC, color = label),epi_markers,
              width =0.4,size = 1)+
  scale_color_manual(values = c("#80B1D3",'#E26565'))+
  scale_y_continuous(
    breaks = seq(-10, 10, by = 5),  # 从 -6 到 8，间隔为 2
  ) +
  theme_classic()
p1
##
library(RColorBrewer)
ctys = unique(topgene$celltype)
mycol <- colorRampPalette(rev(brewer.pal(n = 7, name ="Set1")))(length(ctys))
p2 <- p1 + 
  geom_tile(aes(x = ctys,y = 0),
            height = 0.6,fill = mycol, show.legend = F)+
  geom_text(aes(x= ctys, y = 0, label = ctys),
            size = 3,fontface = "bold")
p2
library(ggrepel)
##
p3 <- p2 + 
  geom_text_repel(aes(x = celltype,y = avg_log2FC,label = gene),
                  topgene,size = 3 )+
  labs(x = "CellType",y = "Average log2FoldChange",
       title = "Differential expression genes")+
  theme(
    plot.title = element_text(size = 14,color = "black",face = "bold"),
    axis.title = element_text(size = 12,color = "black",face = "bold"),
    axis.line.y = element_line(color = "black",linewidth = 0.3),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position  = c(0.98,0.96),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )+
  guides(color = guide_legend(override.aes = list(size = 4)))  
p3




