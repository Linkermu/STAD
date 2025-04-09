
##
result <- nmf(exp,
              rank = 2:6,
              method = "lee",
              nrun = 50,
              seed = 123)

plot(result)

#  
nmf_rank <- nmf(exp,
                rank = 3,
                method = "lee",
                nrun = 50,
                seed = 123)

group  <- predict(nmf_rank) 
table(group)
jco <- c("#2874C5","#EABF00","#C6524A", "#C38E2B")
pdf("geo.pdf", width = 8, height = 8)
consensusmap(nmf_rank,labRow = NA,labCol = NA,
             annCol = data.frame("cluster"=group[colnames(exp)]),
             annColors = list(cluster=jco))
dev.off()
pdf("geo1.pdf", width = 5, height = 5)
basismap(nmf_rank,
         cexCol = 1.5,
         cexRow = 1.5,
         annColors=list(c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
dev.off()
