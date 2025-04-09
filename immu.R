mycolors <-c("#658E67","#E1C62F",
             "#BEAED4","#7570B3",'#446983',
             "#8f657d","#E5C494",'#739B57',"#FBB4AE"
             
             ,"#C0645A","#64495D","#8DD3C7","#B3DE69",
             "#1B9E77","#BBBE00",'#CE3D33',"#FDB462",
             '#5DB1DC','#AF9E85',"#F1E2CC")
sce_immu <- subset(sce.all.int,idents=c(0,2,3,7))
sce_immu <- NormalizeData(sce_immu, normalization.method = "LogNormalize", scale.factor = 10000)
sce_immu <- FindVariableFeatures(sce_immu, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(sce_immu)
sce_immu <- ScaleData(sce_immu, features = scale.genes)
sce_immu <- RunPCA(sce_immu, features = VariableFeatures(sce_immu))
ElbowPlot(sce_immu)
sce_immu <- FindNeighbors(sce_immu, dims = 1:10,reduction = "pca")  
sce_immu <- FindClusters(sce_immu, resolution = 0.5)  
sce_immu <- RunUMAP(sce_immu, dims = 1:10,reduction = "pca",n.neighbors = 20,
                    min.dist = 0.5,spread = 3)
DimPlot(sce_immu, reduction = "umap",label = T,cols =mycolors)  
DimPlot(sce_immu,group.by = "orig.ident",cols = mycolors,label = F,)

marker <- FindAllMarkers(sce_immu, only.pos = TRUE,
                         min.pct = 0.25,
                         test.use = "wilcox", 
                         logfc.threshold = 0.25) 
top20 <- marker %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

top20 %>% filter(cluster == 15) %>% pull(gene)

genes_to_check <- c("CD8B","TRAC","GZMK", #cd8
                    "KLRK1","FYB1","GIMAP5" ,#NK
                    "CD3D","CD3E","CD4", #cd4
                    "FOXP3","CTLA4", #Treg
                    "MS4A1","CD19","CD79B", #B
                    "PAX5","CXCR5", # b_naive
                    "JCHAIN", "MZB1", "XBP1" , #Plasma
                    'IGHG4',# lg G
                    'MUC5AC', #lg_muc5ac
                    'IGKV1-5', # lg κ
                    'IGHM', # lg M 
                    "IGHA1", # lg A
                    "FCN1", "S100A9", "TREM1", # mono
                    "TPSAB1","CPA3","CD9", #mast
                    'C1QA','APOE',"IL2RA", # M1 
                    "CD163","MRC1", #M2
                    "ITGAX","CD83","CD86" #dc
)

DotPlot(sce_immu, features = genes_to_check)+RotatedAxis()

###定义细胞类型
celltype=data.frame(ClusterID=0:15,
                    celltype= 0:15) 
celltype[celltype$ClusterID %in% c( 0),2]='CD8+_T'    
celltype[celltype$ClusterID %in% c( 1),2]='NK'    
celltype[celltype$ClusterID %in% c(2),2]='Plasma_c1'
celltype[celltype$ClusterID %in% c(3 ),2]='Monocyte'   
celltype[celltype$ClusterID %in% c( 4 ),2]='B_naive'
celltype[celltype$ClusterID %in% c(5 ),2]='Plasma_IgG'
celltype[celltype$ClusterID %in% c( 6 ),2]='Plasma_c2'
celltype[celltype$ClusterID %in% c( 7 ),2]='B'
celltype[celltype$ClusterID %in% c( 8 ),2]='Plasma_IgM/A'
celltype[celltype$ClusterID %in% c( 9 ),2]='CD4+_T'     
celltype[celltype$ClusterID %in% c( 10 ),2]='CD4+_Treg'   
celltype[celltype$ClusterID %in% c( 11 ),2]='Mast'      
celltype[celltype$ClusterID %in% c( 12 ),2]='DC'       
celltype[celltype$ClusterID %in% c( 13 ),2]="Mac_M1"        
celltype[celltype$ClusterID %in% c( 14 ),2]='Mac_M2'               
celltype[celltype$ClusterID %in% c( 15 ),2]='Plasma_MUC5AC'      

sce_immu@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce_immu@meta.data[which(sce_immu@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce_immu@meta.data$celltype)
sce_immu <- SetIdent(sce_immu,value = "celltype")

DimPlot(sce_immu, reduction = "umap",label = T,cols = mycolors) 

sce_immu@meta.data[["celltype1"]] = factor(sce_immu@meta.data[["celltype"]],levels = c(
  "B" ,"DC","NK" ,"Mast"  ,"Monocyte" ,"CD8+_T" ,"CD4+_Treg" ,"CD4+_T","Mac_M1" , 
  "Mac_M2","B_naive","Plasma_c1" , "Plasma_c2",
  "Plasma_IgG" ,"Plasma_IgM/A", "Plasma_MUC5AC"
))

DimPlot(sce_immu, reduction = "umap",label = T,cols = mycolors,group.by = "celltype1") 








