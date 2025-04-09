
library(survival)
library(dplyr)
library(glmnet)

data <- fread("../exp_data/clinical_log_scale.csv") 
gene <- all_gene
# 
data <- as.data.frame(data)
x <- as.matrix( data[,colnames(data)%in%gene] )
y <- as.matrix( data[c("time","state")] )
colnames(y) <- c("time","status")

fit <- glmnet(x,y,family="cox",alpha = 1,nlambda = 100)
plot(fit,label = T, xvar = "norm")  
plot(fit,label = T, xvar = "lambda")
print(fit)

##
cvfit <- cv.glmnet(x,y,family="cox",type.measure = "deviance")
plot(cvfit)

##  
cvfit$lambda.min
cvfit$lambda.1se
coef <- coef(cvfit,s="lambda.min")
coef2 <- coef(cvfit,s="lambda.1se")
beta <-  coef[which(coef != 0)]
name <-  row.names(coef)[which(coef != 0)]
geneCoef <- cbind(Gene=name,Coef=beta)
geneCoef
