
library(Mime1)

Training <- fread("/TCGA/exp_data/clinical_log.csv")
Training <- Training[,c(2,11,12,13:19950)]
colnames( Training )[1] <- "ID"
colnames( Training )[2] <- "OS.time"
colnames( Training )[3] <- "OS"
Training$OS.time <- as.numeric(Training$OS.time*30)

Validation <- read.csv("/GEO/GSE66229/GSE66229.csv")
Validation$X <- NULL
Validation <- Validation %>%  relocate(time, .before = state)
colnames(Validation)[2] <- "OS.time"
colnames(Validation)[3] <- "OS"
Validation$OS.time <- as.numeric(Validation$OS.time *30)


gene %in% colnames(Training)
gene %in% colnames(Validation)

list_train_vali_Data <- list(Training=Training,Validation=Validation)

genelist <- gene

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Training,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize =5,
                       seed = 12345)

###  
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1], order = names(list_train_vali_Data),width = 0.35)

#
cindex_dis_select(res,model="StepCox[forward] + Ridge", order= names(list_train_vali_Data))

#
survplot <- vector("list",2) 
for (i in c(1:2)) {  
  print(survplot[[i]] <- rs_sur(res, 
                                model_name = "StepCox[forward] + Ridge", 
                                dataset = names(list_train_vali_Data)[i], 
                                color=c("blue","red"), 
                                median.line = "hv",
                                cutoff = 0.5,
                                conf.int = T,
                                xlab="Days",
                                pval.coord=c(2200,0.9))  
  )
}
aplot::plot_list(gglist=survplot,ncol=2)

# 
head(res$riskscore$`StepCox[forward] + Ridge`[[1]])


# 
unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,
                                   optimal.model = "StepCox[forward] + Ridge",
                                   type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))

all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                             train_data = list_train_vali_Data[["Training"]], 
                             inputmatrix.list = list_train_vali_Data,
                             mode = 'all',
                             AUC_time = 1,  
                             auc_cal_method="KM")

auc_dis_all(all.auc.1y,       
            dataset = names(list_train_vali_Data),     
            validate_set=names(list_train_vali_Data)[-1],     
            order= names(list_train_vali_Data),    
            width = 0.35,      
            year=1)

roc_vis(all.auc.1y,   
        model_name = "StepCox[forward] + Ridge",     
        dataset = names(list_train_vali_Data),      
        order= names(list_train_vali_Data),  
        anno_position=c(0.65,0.55),     
        year=1)

#
res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Training,
                                              candidate_genes = genelist,
                                              mode = "all",nodesize =5,seed = 5201314 )
# 
core_feature_select(res.feature.all)

#
core_feature_rank(res.feature.all, top=10)

#
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                             train_data = list_train_vali_Data[["Training"]],  
                             inputmatrix.list = list_train_vali_Data,
                             mode = 'all',
                             AUC_time = 3,                      
                             auc_cal_method="KM")

all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,
                             train_data = list_train_vali_Data[["Training"]],     
                             inputmatrix.list = list_train_vali_Data,
                             mode = 'all',
                             AUC_time = 5,                   
                             auc_cal_method="KM")

roc_vis(all.auc.3y,   
        model_name = "StepCox[forward] + Ridge",     
        dataset = names(list_train_vali_Data),      
        order= names(list_train_vali_Data),  
        anno_position=c(0.65,0.55),     
        year=3)

roc_vis(all.auc.5y,   
        model_name = "StepCox[forward] + Ridge",     
        dataset = names(list_train_vali_Data),      
        order= names(list_train_vali_Data),  
        anno_position=c(0.65,0.55),     
        year=5)

