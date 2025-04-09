meta = sce_immu@meta.data
meta = left_join(meta,group,by = c("orig.ident"= "id"))
sce_immu$group = meta$group

library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)

##
immu_tumor  =  subset(sce_immu, sample == "tumor")
data <- sce_immu2@meta.data
# Ro/e 
Roe <- calTissueDist(data,
                     byPatient = F,
                     colname.cluster = "celltype2", # 不同细胞亚群
                     colname.patient = "orig.ident", # 不同样本,这里用orig信息来替代
                     colname.tissue = "group", # 根据自己的选择不同分组
                     method = "chisq", # "chisq", "fisher", and "freq" 
                     min.rowSum = 0) 
Roe
##
col_fun <-colorRamp2(c(min(Roe,na.rm = TRUE),1,max(Roe,na.rm = TRUE)),
                     c('#3D98D3FF', "white", '#800202'))
Heatmap(as.matrix(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e Index",
          at = seq(0.5, 1.5, by = 0.5), 
          labels = seq(0.5, 1.5, by = 0.5) 
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 8, col = "black"))
        }
)

####  cluster

library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

sce_high_immu = subset(sce_immu2, group == 'high')
test <- subset(sce_high_immu,downsample = 1000)
# 
mat.pca=test@reductions$pca@cell.embeddings
mat.pca =mat.pca[rownames(test@meta.data),]
mat.pca =as.data.frame(mat.pca)
mat.pca$celltype<-test@meta.data$celltype2

result <- mat.pca %>% group_by(celltype) %>% summarise(across(PC_1:PC_40, ~mean(.x, na.rm = TRUE)))
result<-result%>% column_to_rownames(var="celltype")
result<-t(result)

library(ggplot2)
library(ggdendro)
hc = hclust(dist(t(as.matrix((result)))),method = 'ward.D')
hc = as.dendrogram(hc)
hc_dat<-dendro_data(hc,type="triangle")

ggplot() + 
  geom_segment(data= hc_dat$segments, aes(x = x, y =-y, xend = xend, yend =-yend),size=1, color='#8c6bb1') +
  geom_text(data= hc_dat$labels, aes(x = x, y =-y+5, label = label), hjust =0, angle =0) +
  geom_point(data= hc_dat$labels, aes(x = x, y =-y, color = label),size=6)+
  geom_text(data= hc_dat$labels, aes(x = x, y =-y, label = x))+
  theme_dendro()+ 
  coord_flip()+ 
  ylim(-60,45)+ 
  NoLegend()+
  scale_color_manual(values = c("#d2981a","#a53e1f","#457277","#8f657d","#42819F","#86AA7D","#CBB396",
                                "#A65628" ,"#F781BF", "#999999", "#1B9E77", "#D95F02",
                                "#7570B3" , "#66A61E" ,"#E6AB02" ,"#A6761D","#AA3B50"))
