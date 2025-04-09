data1$color <- ifelse(data1$b > 0,"right","left")
data1$color <- ifelse(data1$pval > 0.05,"bottom",data1$color)
data1$color <- as.factor(data1$color)

ggplot(data1,aes(b, -log10(pval)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "#999999")+
  geom_point(aes(size=-log10(pval), color= color,alpha = 1))+
  scale_color_manual(values = c("left" = "#4F4FFF", "right" = "#E31A1C", "bottom" = "grey"))+
  scale_size_continuous(range = c(1.2,2.2))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9),
        legend.justification = c(0,1),
        legend.text = element_text(size = 11),  
        legend.key.size = unit(1.0, "lines")
  )+
  guides( color = "none", alpha="none", 
          size = "none"  
  )+
  geom_label_repel(aes(label = label, color = color), size = 3)+
  xlab("In(OR)")+
  ylab("-Log10(p-value)")
