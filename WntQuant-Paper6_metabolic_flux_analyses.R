# ==============================================================================
# WntQuant-Paper6 -- Metabolic flux analyses
# ------------------------------------------------------------------------------
# Purpose:    Infer metabolic flux (METAFlux) and link Wnt activity to metabolic reprogramming across cell lines and tumor cohorts.
# Inputs:     Expression matrices (TCGA, TARGET, GDSC, CCLE, CTRP) and Wnt scores.
# Outputs:    Metabolic flux scores and pathway-level activity estimates.
# NOTE:       Part of the code involves the core interests of the laboratory and is not open source at this time.
# ==============================================================================


rm(list=ls())
gc()
library(METAFlux)
setwd("input")
load("input/pancancer_exp.Rdata")
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
pancancer_exp=pancancer_exp[,rownames(pancancer_group)]
scores=calculate_reaction_score(pancancer_exp)
metabolism_scores=scores
data("human_blood")
flux=compute_flux(mras=scores,medium=human_blood)
metabolism_flux=flux
load("input/metabolism_flux.Rdata")
flux=metabolism_flux
cbrt=function(x){
  sign(x)*abs(x)^(1/3)
}
flux1=cbrt(flux)
data("nutrient_lookup_files")
pathway=unique(unlist(human_gem$SUBSYSTEM))
pathway_score=list()
for(i in pathway){
  path=i
  activity_score=c()
  for (d in 1:ncol(flux1)){
    activity_score[d]=mean(abs(flux1[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    names(activity_score)[d]=colnames(flux1)[d]
  }
  pathway_score[[i]]=activity_score
}
all_pathway_score=as.data.frame(do.call(rbind,pathway_score))
TCGA_metabolism=all_pathway_score
save(TCGA_metabolism,file="TCGA_metabolism.Rdata")
# ------------------------------------------------------------------------------
# TARGET
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(METAFlux)
setwd("input")
load("input/TARGET.Rdata")
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
TARGET=TARGET[,rownames(TARGET_G)]
scores=calculate_reaction_score(TARGET)
data("human_blood")
flux=compute_flux(mras=scores,medium=human_blood)
cbrt=function(x){
  sign(x)*abs(x)^(1/3)
}
flux1=cbrt(flux)
data("nutrient_lookup_files")
pathway=unique(unlist(human_gem$SUBSYSTEM))
pathway_score=list()
for(i in pathway){
  path=i
  activity_score=c()
  for (d in 1:ncol(flux1)){
    activity_score[d]=mean(abs(flux1[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    names(activity_score)[d]=colnames(flux1)[d]
  }
  pathway_score[[i]]=activity_score
}
all_pathway_score=as.data.frame(do.call(rbind,pathway_score))
TARGET_metabolism=all_pathway_score
save(TARGET_metabolism,file="TARGET_metabolism.Rdata")
# ------------------------------------------------------------------------------
# GDSC1
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(METAFlux)
setwd("input")
GDSC1_exp=readRDS("input/DataFiles/DataFiles/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
scores=calculate_reaction_score(GDSC1_exp)
metabolism_scores=scores
data("cell_medium")
flux=compute_flux(mras=scores,medium=cell_medium)
cbrt=function(x){
  sign(x)*abs(x)^(1/3)
}
flux1=cbrt(flux)
data("nutrient_lookup_files")
pathway=unique(unlist(human_gem$SUBSYSTEM))
pathway_score=list()
for(i in pathway){
  path=i
  activity_score=c()
  for (d in 1:ncol(flux1)){
    activity_score[d]=mean(abs(flux1[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    names(activity_score)[d]=colnames(flux1)[d]
  }
  pathway_score[[i]]=activity_score
}
all_pathway_score=as.data.frame(do.call(rbind,pathway_score))
GDSC1_metabolism=all_pathway_score
save(GDSC1_metabolism,file="GDSC1_metabolism.Rdata")
# ------------------------------------------------------------------------------
# GDSC2
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(METAFlux)
setwd("input")
GDSC2_exp=readRDS("input/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
scores=calculate_reaction_score(GDSC2_exp)
metabolism_scores=scores
data("cell_medium")
flux=compute_flux(mras=scores,medium=cell_medium)
cbrt=function(x){
  sign(x)*abs(x)^(1/3)
}
flux1=cbrt(flux)
data("nutrient_lookup_files")
pathway=unique(unlist(human_gem$SUBSYSTEM))
pathway_score=list()
for(i in pathway){
  path=i
  activity_score=c()
  for (d in 1:ncol(flux1)){
    activity_score[d]=mean(abs(flux1[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    names(activity_score)[d]=colnames(flux1)[d]
  }
  pathway_score[[i]]=activity_score
}
all_pathway_score=as.data.frame(do.call(rbind,pathway_score))
GDSC2_metabolism=all_pathway_score
save(GDSC2_metabolism,file="GDSC2_metabolism.Rdata")
# ------------------------------------------------------------------------------
# CTRP
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(METAFlux)
setwd("input")
CTRP_exp=readRDS("input/DataFiles/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP_exp=log2(CTRP_exp+1)
scores=calculate_reaction_score(CTRP_exp)
metabolism_scores=scores
data("cell_medium")
flux=compute_flux(mras=scores,medium=cell_medium)
cbrt=function(x){
  sign(x)*abs(x)^(1/3)
}
flux1=cbrt(flux)
data("nutrient_lookup_files")
pathway=unique(unlist(human_gem$SUBSYSTEM))
pathway_score=list()
for(i in pathway){
  path=i
  activity_score=c()
  for (d in 1:ncol(flux1)){
    activity_score[d]=mean(abs(flux1[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    names(activity_score)[d]=colnames(flux1)[d]
  }
  pathway_score[[i]]=activity_score
}
all_pathway_score=as.data.frame(do.call(rbind,pathway_score))
CTRP_metabolism=all_pathway_score
save(CTRP_metabolism,file="CTRP_metabolism.Rdata")
# ------------------------------------------------------------------------------
# CCLE
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(METAFlux)
setwd("input")
load("input/CCLE_exp.Rdata")
scores=calculate_reaction_score(CCLE_exp)
metabolism_scores=scores
data("cell_medium")
flux=compute_flux(mras=scores,medium=cell_medium)
cbrt=function(x){
  sign(x)*abs(x)^(1/3)
}
flux1=cbrt(flux)
data("nutrient_lookup_files")
pathway=unique(unlist(human_gem$SUBSYSTEM))
pathway_score=list()
for(i in pathway){
  path=i
  activity_score=c()
  for (d in 1:ncol(flux1)){
    activity_score[d]=mean(abs(flux1[which(unlist(human_gem$SUBSYSTEM)==i),d]))
    names(activity_score)[d]=colnames(flux1)[d]
  }
  pathway_score[[i]]=activity_score
}
all_pathway_score=as.data.frame(do.call(rbind,pathway_score))
CCLE_metabolism=all_pathway_score
save(CCLE_metabolism,file="CCLE_metabolism.Rdata")
# ------------------------------------------------------------------------------
# xgboost
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(shapviz)
library(xgboost)
library(ggplot2)
library(caret)
library(Matrix)
load("input/TCGA_metabolism.Rdata")
load("input/TARGET_metabolism.Rdata")
load("input/pancancer_WNT_Score.Rdata")
load("input/TARGET_WNT_Score.Rdata")
TCGA=merge(pancancer_WNT_Score$final_activity_score[,c(1),drop=F], 
           as.data.frame(t(TCGA_metabolism)), by="row.names")
rownames(TCGA)=TCGA[,1]
TCGA=TCGA[,-1]
set.seed(123)
xgb_grid = expand.grid(
  nrounds = c(50, 100, 150),
  max_depth = c(3, 6, 9),
  eta = c(0.01, 0.1, 0.3),
  gamma = c(0, 0.1, 0.2),
  colsample_bytree = c(0.6, 0.8, 1),
  min_child_weight = c(1, 3, 5),
  subsample = c(0.5, 0.75, 1)
)
xgb_trcontrol=trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  returnData = FALSE,
  returnResamp = "all",
  allowParallel = F
)
xgb_train = train(
  x = data.matrix(TCGA[, -1]),
  y = TCGA[, 1],
  trControl = xgb_trcontrol,
  tuneGrid = xgb_grid[sample(1:nrow(xgb_grid), 10), ],
  method = "xgbTree",
  nthread = 1
)
print(xgb_train$bestTune)
best_params = xgb_train$bestTune
fit = xgboost(
  data = data.matrix(TCGA[, -1]),
  label = TCGA[, 1],
  nrounds = best_params$nrounds,
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  gamma = best_params$gamma,
  colsample_bytree = best_params$colsample_bytree,
  min_child_weight = best_params$min_child_weight,
  subsample = best_params$subsample,
  nthread = 1,
  verbose = 0
)
results=xgb_train$results
TCGA_pred = TCGA
TCGA_pred$pred = predict(fit, as.matrix(TCGA_pred[, -1]))
library(ggplot2)
library(ggpubr)
p=ggscatter(TCGA_pred, x = "pred", y = "activity_score",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "#d71345", fill = "#BBBDBE"),
            color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -2,
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TCGA-XGBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
shap_xgboost = shapviz(fit, X_pred = data.matrix(TCGA[, -1]), X = TCGA[, -1])
sv_importance(shap_xgboost, kind = "beeswarm", max_display = 20) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"))

impor = as.data.frame(sort(colMeans(abs(shap_xgboost$S)), decreasing = TRUE))
colnames(impor) = "XGBoost"
xgboost_impor = impor
save(xgboost_impor, file = "xgboost_impor.Rdata")
# ------------------------------------------------------------------------------
# TARGET
# ------------------------------------------------------------------------------
load("input/TARGET_metabolism.Rdata")
load("input/TARGET_WNT_Score.Rdata")
TARGET=merge(TARGET_WNT_Score$final_activity_score[,c(1),drop=F], 
             as.data.frame(t(TARGET_metabolism)), by="row.names")
rownames(TARGET)=TARGET[,1]
TARGET=TARGET[,-1]
TARGET_pred=TARGET
TARGET_pred$pred = predict(fit, as.matrix(TARGET_pred[, -1]))
p=ggscatter(TARGET_pred, x = "pred", y = "activity_score",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "#d71345", fill = "#BBBDBE"),
            color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TARGET-XGBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC1
# ------------------------------------------------------------------------------
load("input/GDSC1_metabolism.Rdata")
load("input/GDSC1_WNT_score.Rdata")
GDSC1=merge(GDSC1_WNT_score[,c(1),drop=F], 
            as.data.frame(t(GDSC1_metabolism)), by="row.names")
rownames(GDSC1)=GDSC1[,1]
GDSC1=GDSC1[,-1]
GDSC1_pred=GDSC1
GDSC1_pred$pred=predict(fit, as.matrix(GDSC1_pred[, -1]))
p=ggscatter(GDSC1_pred, x = "pred", y = "activity_score",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "#d71345", fill = "#BBBDBE"),
            color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC1-XGBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC2
# ------------------------------------------------------------------------------
load("input/GDSC2_metabolism.Rdata")
load("input/GDSC2_WNT_score.Rdata")
GDSC2=merge(GDSC2_WNT_score[,c(1),drop=F], 
            as.data.frame(t(GDSC2_metabolism)), by="row.names")
rownames(GDSC2)=GDSC2[,1]
GDSC2=GDSC2[,-1]
GDSC2_pred=GDSC2
GDSC2_pred$pred=predict(fit, as.matrix(GDSC2_pred[, -1]))
p=ggscatter(GDSC2_pred, x = "pred", y = "activity_score",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "#d71345", fill = "#BBBDBE"),
            color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC2-XGBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CCLE
# ------------------------------------------------------------------------------
load("input/CCLE_metabolism.Rdata")
load("input/CCLE_WNT_score.Rdata")
CCLE=merge(CCLE_WNT_score[,c(1),drop=F], 
           as.data.frame(t(CCLE_metabolism)), by="row.names")
rownames(CCLE)=CCLE[,1]
CCLE=CCLE[,-1]
CCLE_pred=CCLE
CCLE_pred$pred=predict(fit, as.matrix(CCLE_pred[, -1]))
p=ggscatter(CCLE_pred, x = "pred", y = "activity_score",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "#d71345", fill = "#BBBDBE"),
            color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CCLE-XGBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CTRP
# ------------------------------------------------------------------------------
load("input/CTRP_metabolism.Rdata")
load("input/CTRP_WNT_score.Rdata")
CTRP=merge(CTRP_WNT_score[,c(1),drop=F], 
           as.data.frame(t(CTRP_metabolism)), by="row.names")
rownames(CTRP)=CTRP[,1]
CTRP=CTRP[,-1]
CTRP_pred=CTRP
CTRP_pred$pred=predict(fit, as.matrix(CTRP_pred[, -1]))
p=ggscatter(CTRP_pred, x = "pred", y = "activity_score",
            add = "reg.line",
            conf.int = TRUE, 
            add.params = list(color = "#d71345", fill = "#BBBDBE"),
            color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CTRP-XGBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
xgboost_fit=fit
save(xgboost_fit,file="xgboost_fit.Rdata")
# ------------------------------------------------------------------------------
# lightgbm
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(shapviz)
library(lightgbm)
library(ggplot2)
library(caret)
library(Matrix)
load("TCGA_metabolism.Rdata")
name=rownames(TCGA_metabolism)
load("pancancer_WNT_Score.Rdata")
TCGA = merge(pancancer_WNT_Score$final_activity_score[,c(1),drop=F], 
              as.data.frame(t(TCGA_metabolism)), by="row.names")
rownames(TCGA) = TCGA[,1]
TCGA = TCGA[,-1]
TCGA = as.matrix(TCGA)
colnames(TCGA) = gsub("[^a-zA-Z0-9_.]", "_", colnames(TCGA))
load("TARGET_metabolism.Rdata")
load("TARGET_WNT_Score.Rdata")
TARGET = merge(TARGET_WNT_Score$final_activity_score[,c(1),drop=F], 
                as.data.frame(t(TARGET_metabolism)), by="row.names")
rownames(TARGET) = TARGET[,1]
TARGET = TARGET[,-1]
TARGET = as.matrix(TARGET)
colnames(TARGET) = gsub("[^a-zA-Z0-9_.]", "_", colnames(TARGET))
set.seed(123)
best_score = Inf
best_params = list()
num_folds = 5
folds = createFolds(TCGA[,1], k = num_folds)
param_grid = expand.grid(
  num_leaves = c(31, 63),
  learning_rate = c(0.05, 0.1),
  feature_fraction = c(0.7, 0.9),
  min_data_in_leaf = c(20, 50)
)
for(i in 1:nrow(param_grid)) {
  cat("Testing parameter set", i, "of", nrow(param_grid), "\n")
  cv_scores = numeric(num_folds)
  for(fold in 1:num_folds) {
    train_idx = unlist(folds[-fold])
    valid_idx = folds[[fold]]
    dtrain = lgb.Dataset(TCGA[train_idx, -1], label = TCGA[train_idx, 1])
    dvalid = lgb.Dataset(TCGA[valid_idx, -1], label = TCGA[valid_idx, 1])
    model = lgb.train(
      params = list(
        objective = "regression",
        metric = "rmse",
        num_leaves = param_grid$num_leaves[i],
        learning_rate = param_grid$learning_rate[i],
        feature_fraction = param_grid$feature_fraction[i],
        min_data_in_leaf = param_grid$min_data_in_leaf[i]
      ),
      data = dtrain,
      valids = list(valid = dvalid),
      nrounds = 100,
      early_stopping_rounds = 10,
      verbose = -1
    )
    cv_scores[fold] = model$best_score
  }
  mean_score = mean(cv_scores)
  cat("Mean RMSE:", mean_score, "\n")
  
  if(mean_score < best_score) {
    best_score = mean_score
    best_params = as.list(param_grid[i,])
    cat("New best parameters found!\n")
  }
}
dtrain=lgb.Dataset(TCGA[, -1], label = TCGA[, 1])
dtest=lgb.Dataset.create.valid(dtrain, data = TARGET[, -1], label = TARGET[, 1])
fit=lgb.train(
  params = c(
    list(
      objective = "regression",
      metric = "rmse",
      boosting_type = "gbdt"
    ),
    best_params
  ),
  data = dtrain,
  valids = list(test = dtest),
  nrounds = 100,
  early_stopping_rounds = 10,
  verbose = 1
)
TCGA_pred = as.data.frame(TCGA)
TCGA_pred$pred = predict(fit, as.matrix(TCGA_pred[, -1]))
p = ggscatter(TCGA_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0,
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TCGA-LightGBM (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
shap_lightgbm = shapviz(fit, X_pred = TCGA[, -1])
colnames(shap_lightgbm$X) = name
colnames(shap_lightgbm$S) = name
sv_importance(shap_lightgbm, kind = "beeswarm", max_display = 20) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"))

impor = as.data.frame(sort(colMeans(abs(shap_lightgbm$S)), decreasing = TRUE))
colnames(impor) = "LightGBM"
lightgbm_impor = impor
save(lightgbm_impor, file = "lightgbm_impor.Rdata")
TARGET_pred = as.data.frame(TARGET)
TARGET_pred$pred = predict(fit, as.matrix(TARGET_pred[, -1]))
p = ggscatter(TARGET_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TARGET-LightGBM (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC1
# ------------------------------------------------------------------------------
load("GDSC1_metabolism.Rdata")
load("GDSC1_WNT_score.Rdata")
GDSC1 = merge(GDSC1_WNT_score[,c(1),drop=F], 
               as.data.frame(t(GDSC1_metabolism)), by="row.names")
rownames(GDSC1) = GDSC1[,1]
GDSC1 = GDSC1[,-1]
GDSC1 = as.matrix(GDSC1)
colnames(GDSC1) = gsub("[^a-zA-Z0-9_.]", "_", colnames(GDSC1))
GDSC1_pred = as.data.frame(GDSC1)
GDSC1_pred$pred = predict(fit, as.matrix(GDSC1_pred[, -1]))
p = ggscatter(GDSC1_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC1-LightGBM (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC2
# ------------------------------------------------------------------------------
load("GDSC2_metabolism.Rdata")
load("GDSC2_WNT_score.Rdata")
GDSC2 = merge(GDSC2_WNT_score[,c(1),drop=F], 
               as.data.frame(t(GDSC2_metabolism)), by="row.names")
rownames(GDSC2) = GDSC2[,1]
GDSC2 = GDSC2[,-1]
GDSC2 = as.matrix(GDSC2)
colnames(GDSC2) = gsub("[^a-zA-Z0-9_.]", "_", colnames(GDSC2))
GDSC2_pred = as.data.frame(GDSC2)
GDSC2_pred$pred = predict(fit, as.matrix(GDSC2_pred[, -1]))
p = ggscatter(GDSC2_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC2-LightGBM (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CCLE
# ------------------------------------------------------------------------------
load("CCLE_metabolism.Rdata")
load("CCLE_WNT_score.Rdata")
CCLE = merge(CCLE_WNT_score[,c(1),drop=F], 
              as.data.frame(t(CCLE_metabolism)), by="row.names")
rownames(CCLE) = CCLE[,1]
CCLE = CCLE[,-1]
CCLE = as.matrix(CCLE)
colnames(CCLE) = gsub("[^a-zA-Z0-9_.]", "_", colnames(CCLE))
CCLE_pred = as.data.frame(CCLE)
CCLE_pred$pred = predict(fit, as.matrix(CCLE_pred[, -1]))
p = ggscatter(CCLE_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CCLE-LightGBM (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CTRP
# ------------------------------------------------------------------------------
load("CTRP_metabolism.Rdata")
load("CTRP_WNT_score.Rdata")
CTRP = merge(CTRP_WNT_score[,c(1),drop=F], 
              as.data.frame(t(CTRP_metabolism)), by="row.names")
rownames(CTRP) = CTRP[,1]
CTRP = CTRP[,-1]
CTRP = as.matrix(CTRP)
colnames(CTRP) = gsub("[^a-zA-Z0-9_.]", "_", colnames(CTRP))
CTRP_pred = as.data.frame(CTRP)
CTRP_pred$pred = predict(fit, as.matrix(CTRP_pred[, -1]))
p = ggscatter(CTRP_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CTRP-LightGBM (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
lightgbm_fit=fit
save(lightgbm_fit,file="lightgbm_fit.Rdata")
# ------------------------------------------------------------------------------
# rf
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(shapviz)
library(treeshap)
library(ranger)
library(ggplot2)
library(caret)
library(doParallel)
set.seed(123)
load("TCGA_metabolism.Rdata")
load("pancancer_WNT_Score.Rdata")
TCGA = merge(pancancer_WNT_Score$final_activity_score[,c(1),drop=F], 
              as.data.frame(t(TCGA_metabolism)), by="row.names")
rownames(TCGA) = TCGA[,1]
TCGA = TCGA[,-1]
cl = makePSOCKcluster(4)
registerDoParallel(cl)
rf_grid = expand.grid(
  mtry = seq(10, 50, by = 10),
  min.node.size = c(1, 5, 10),
  splitrule = "variance"
)
ctrl = trainControl(
  method = "cv",
  number = 5,
  verboseIter = TRUE,
  allowParallel = TRUE
)
rf_model = train(
  x = TCGA[, -1],
  y = TCGA[, 1],
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = ctrl,
  num.trees = 100,
  importance = "permutation",
  verbose = TRUE
)
stopCluster(cl)
print(rf_model$bestTune)
fit = ranger(
  y = TCGA$activity_score,
  x = TCGA[, -1],
  mtry = rf_model$bestTune$mtry,
  min.node.size = rf_model$bestTune$min.node.size,
  splitrule = rf_model$bestTune$splitrule,
  num.trees = 100,
  importance = "permutation"
)
TCGA_pred = TCGA
TCGA_pred$pred = predict(fit, data = TCGA_pred[, -1])$predictions
p = ggscatter(TCGA_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0,
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TCGA-RF (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
unified_model = ranger.unify(fit, TCGA[, -1])
shaps = treeshap(unified_model, TCGA[, -1])
shap_rf = shapviz(shaps, X = TCGA[, -1])
sv_importance(shap_rf, kind = "beeswarm", max_display = 20) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"))

impor = as.data.frame(sort(colMeans(abs(shap_rf$S)), decreasing = TRUE))
colnames(impor) = "RF"
rf_impor = impor
save(rf_impor, file = "rf_impor.Rdata")
# ------------------------------------------------------------------------------
# TARGET
# ------------------------------------------------------------------------------
load("TARGET_metabolism.Rdata")
load("TARGET_WNT_Score.Rdata")
TARGET = merge(TARGET_WNT_Score$final_activity_score[,c(1),drop=F], 
                as.data.frame(t(TARGET_metabolism)), by="row.names")
rownames(TARGET) = TARGET[,1]
TARGET = TARGET[,-1]
TARGET_pred = TARGET
TARGET_pred$pred = predict(fit, data = TARGET_pred[, -1])$predictions
p = ggscatter(TARGET_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TARGET-RF (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC1
# ------------------------------------------------------------------------------
load("GDSC1_metabolism.Rdata")
load("GDSC1_WNT_score.Rdata")
GDSC1 = merge(GDSC1_WNT_score[,c(1),drop=F], 
               as.data.frame(t(GDSC1_metabolism)), by="row.names")
rownames(GDSC1) = GDSC1[,1]
GDSC1 = GDSC1[,-1]
GDSC1_pred = GDSC1
GDSC1_pred$pred = predict(fit, data = GDSC1_pred[, -1])$predictions
p = ggscatter(GDSC1_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC1-RF (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC2
# ------------------------------------------------------------------------------
load("GDSC2_metabolism.Rdata")
load("GDSC2_WNT_score.Rdata")
GDSC2 = merge(GDSC2_WNT_score[,c(1),drop=F], 
               as.data.frame(t(GDSC2_metabolism)), by="row.names")
rownames(GDSC2) = GDSC2[,1]
GDSC2 = GDSC2[,-1]
GDSC2_pred = GDSC2
GDSC2_pred$pred = predict(fit, data = GDSC2_pred[, -1])$predictions
p = ggscatter(GDSC2_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC2-RF (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CCLE
# ------------------------------------------------------------------------------
load("CCLE_metabolism.Rdata")
load("CCLE_WNT_score.Rdata")
CCLE = merge(CCLE_WNT_score[,c(1),drop=F], 
              as.data.frame(t(CCLE_metabolism)), by="row.names")
rownames(CCLE) = CCLE[,1]
CCLE = CCLE[,-1]
CCLE_pred = CCLE
CCLE_pred$pred = predict(fit, data = CCLE_pred[, -1])$predictions
p = ggscatter(CCLE_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CCLE-RF (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CTRP
# ------------------------------------------------------------------------------
load("CTRP_metabolism.Rdata")
load("CTRP_WNT_score.Rdata")
CTRP = merge(CTRP_WNT_score[,c(1),drop=F], 
              as.data.frame(t(CTRP_metabolism)), by="row.names")
rownames(CTRP) = CTRP[,1]
CTRP = CTRP[,-1]
CTRP_pred = CTRP
CTRP_pred$pred = predict(fit, data = CTRP_pred[, -1])$predictions
p = ggscatter(CTRP_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0, 
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CTRP-RF (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
rf_fit=fit
save(rf_fit,file="rf_fit.Rdata")
# ------------------------------------------------------------------------------
# catboost
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(shapviz)
library(catboost)
library(ggplot2)
library(caret)
library(doParallel)
shapviz.catboost.Model = function(object, X_pred, X = X_pred, collapse = NULL, ...) {  
  if (!requireNamespace("catboost", quietly = TRUE)) {    
    stop("Package 'catboost' not installed")  
  }  
  stopifnot(    
    "X must be a matrix or data.frame. It can't be an object of class catboost.Pool" = 
      is.matrix(X) || is.data.frame(X),    
    "X_pred must be a matrix, a data.frame, or a catboost.Pool" = 
      is.matrix(X_pred) || is.data.frame(X_pred) || inherits(X_pred, "catboost.Pool"),    
    "X_pred must have column names" = !is.null(colnames(X_pred))  
  )    
  if (!inherits(X_pred, "catboost.Pool")) {    
    X_pred = catboost.load_pool(X_pred)  
  }
  S = catboost.get_feature_importance(object, X_pred, type = "ShapValues", ...)
  pp = ncol(X_pred) + 1L  
  baseline = S[1L, pp]  
  S = S[, -pp, drop = FALSE]  
  colnames(S) = colnames(X_pred)  
  shapviz(S, X = X, baseline = baseline, collapse = collapse)
}
load("TCGA_metabolism.Rdata")
load("pancancer_WNT_Score.Rdata")
TCGA = merge(pancancer_WNT_Score$final_activity_score[,c(1),drop=F], 
              as.data.frame(t(TCGA_metabolism)), by="row.names")
rownames(TCGA) = TCGA[,1]
TCGA = TCGA[,-1]
param_grid = expand.grid(
  depth = c(4, 6, 8),
  learning_rate = c(0.01, 0.05, 0.1),
  l2_leaf_reg = c(1, 3, 5)
)
set.seed(123)
folds = createFolds(TCGA[,1], k = 5)
best_score = Inf
best_params = NULL
for(i in 1:nrow(param_grid)) {
  cat("Testing parameter set", i, "of", nrow(param_grid), "\n")
  cv_scores = numeric(length(folds))
  for(fold in seq_along(folds)) {
    train_idx = unlist(folds[-fold])
    valid_idx = folds[[fold]]
    train_pool = catboost.load_pool(TCGA[train_idx, -1], label = TCGA[train_idx, 1])
    valid_pool = catboost.load_pool(TCGA[valid_idx, -1], label = TCGA[valid_idx, 1])
    model = catboost.train(
      train_pool,
      params = list(
        loss_function = "RMSE",
        iterations = 500,
        depth = param_grid$depth[i],
        learning_rate = param_grid$learning_rate[i],
        l2_leaf_reg = param_grid$l2_leaf_reg[i],
        logging_level = "Silent"
      )
    )
    pred = catboost.predict(model, valid_pool)
    cv_scores[fold] = sqrt(mean((pred - TCGA[valid_idx, 1])^2))
  }
  mean_score = mean(cv_scores)
  cat("Mean RMSE:", mean_score, "\n")
  if(mean_score < best_score) {
    best_score = mean_score
    best_params = param_grid[i,]
    cat("New best parameters found!\n")
  }
}
train_pool = catboost.load_pool(TCGA[, -1], label = TCGA[, 1])
fit = catboost.train(
  train_pool,
  params = list(
    loss_function = "RMSE",
    iterations = 500,
    depth = best_params$depth,
    learning_rate = best_params$learning_rate,
    l2_leaf_reg = best_params$l2_leaf_reg,
    logging_level = "Silent"
  )
)
TCGA_pred = TCGA
TCGA_pool = catboost.load_pool(TCGA_pred[, -1])
TCGA_pred$pred = catboost.predict(fit, TCGA_pool)
p = ggscatter(TCGA_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -2,
           label.y = -2,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TCGA-CatBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
shp = shapviz(fit, X_pred = TCGA[,-1])
sv_importance(shp, kind = "beeswarm", max_display = 20) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"))

impor = as.data.frame(sort(colMeans(abs(shp$S)), decreasing = TRUE))
colnames(impor) = "CatBoost"
catboost_impor = impor
save(catboost_impor, file = "catboost_impor.Rdata")
# ------------------------------------------------------------------------------
# TARGET
# ------------------------------------------------------------------------------
load("TARGET_metabolism.Rdata")
load("TARGET_WNT_Score.Rdata")
TARGET = merge(TARGET_WNT_Score$final_activity_score[,c(1),drop=F], 
                as.data.frame(t(TARGET_metabolism)), by="row.names")
rownames(TARGET) = TARGET[,1]
TARGET = TARGET[,-1]
TARGET_pred = TARGET
TARGET_pool = catboost.load_pool(TARGET_pred[, -1])
TARGET_pred$pred = catboost.predict(fit, TARGET_pool)
p = ggscatter(TARGET_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2.5,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("TARGET-CatBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC1
# ------------------------------------------------------------------------------
load("GDSC1_metabolism.Rdata")
load("GDSC1_WNT_score.Rdata")
GDSC1 = merge(GDSC1_WNT_score[,c(1),drop=F], 
               as.data.frame(t(GDSC1_metabolism)), by="row.names")
rownames(GDSC1) = GDSC1[,1]
GDSC1 = GDSC1[,-1]
GDSC1_pred = GDSC1
GDSC1_pool = catboost.load_pool(GDSC1_pred[, -1])
GDSC1_pred$pred = catboost.predict(fit, GDSC1_pool)
p = ggscatter(GDSC1_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2.5,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC1-CatBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# GDSC2
# ------------------------------------------------------------------------------
load("GDSC2_metabolism.Rdata")
load("GDSC2_WNT_score.Rdata")
GDSC2 = merge(GDSC2_WNT_score[,c(1),drop=F], 
               as.data.frame(t(GDSC2_metabolism)), by="row.names")
rownames(GDSC2) = GDSC2[,1]
GDSC2 = GDSC2[,-1]
GDSC2_pred = GDSC2
GDSC2_pool = catboost.load_pool(GDSC2_pred[, -1])
GDSC2_pred$pred = catboost.predict(fit, GDSC2_pool)
p = ggscatter(GDSC2_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2.5,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("GDSC2-CatBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CCLE
# ------------------------------------------------------------------------------
load("CCLE_metabolism.Rdata")
load("CCLE_WNT_score.Rdata")
CCLE = merge(CCLE_WNT_score[,c(1),drop=F], 
              as.data.frame(t(CCLE_metabolism)), by="row.names")
rownames(CCLE) = CCLE[,1]
CCLE = CCLE[,-1]
CCLE_pred = CCLE
CCLE_pool = catboost.load_pool(CCLE_pred[, -1])
CCLE_pred$pred = catboost.predict(fit, CCLE_pool)
p = ggscatter(CCLE_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2.5,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CCLE-CatBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
# ------------------------------------------------------------------------------
# CTRP
# ------------------------------------------------------------------------------
load("CTRP_metabolism.Rdata")
load("CTRP_WNT_score.Rdata")
CTRP = merge(CTRP_WNT_score[,c(1),drop=F], 
              as.data.frame(t(CTRP_metabolism)), by="row.names")
rownames(CTRP) = CTRP[,1]
CTRP = CTRP[,-1]
CTRP_pred = CTRP
CTRP_pool = catboost.load_pool(CTRP_pred[, -1])
CTRP_pred$pred = catboost.predict(fit, CTRP_pool)
p = ggscatter(CTRP_pred, x = "pred", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE, 
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = -1, 
           label.y = -2.5,
           size = 6) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("CTRP-CatBoost (Tuned)") +
  xlab("Model prediction score") +
  ylab("Wnt/β-catenin pathway activity score")
print(p)
catboost_fit=fit
save(catboost_fit,file="catboost_fit.Rdata")
# ------------------------------------------------------------------------------
# visualisation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/xgboost_impor.Rdata")
load("input/lightgbm_impor.Rdata")
load("input/rf_impor.Rdata")
load("input/catboost_impor.Rdata")
ML_importance=merge(xgboost_impor,lightgbm_impor,by="row.names")
ML_importance=merge(ML_importance,rf_impor,by.x="Row.names",by.y="row.names")
ML_importance=merge(ML_importance,catboost_impor,by.x="Row.names",by.y="row.names")
rownames(ML_importance)=ML_importance$Row.names
ML_importance=as.data.frame(scale(ML_importance[,-1]))
ML_importance$importance=rowMeans(ML_importance)
ML_importance=ML_importance[order(-ML_importance$importance),]
library(openxlsx)
setwd(".")
write.xlsx(ML_importance, "Table S8.xlsx", rowNames = TRUE, colNames = TRUE)
heatmap_data=as.matrix(ML_importance[, -which(colnames(ML_importance) == "importance")])
ML_importance$feature=rownames(heatmap_data)
ML_importance=head(ML_importance,20)
heatmap_data=heatmap_data[rownames(ML_importance),]
library(dplyr)  
library(ggplot2)  
library(patchwork)  
library(tidyr)  
library(tibble)  
sorted_features = ML_importance %>%  
  arrange(desc(importance)) %>%  
  pull(feature)
heatmap_df = as.data.frame(heatmap_data) %>%  
  tibble::rownames_to_column("feature") %>%  
  mutate(feature = factor(feature, levels = sorted_features)) %>%  
  tidyr::pivot_longer(-feature)
heatmap_plot = ggplot(heatmap_df, aes(x = name, y = feature, fill = value)) +  
  geom_tile(color = "gray70", size = 0.5) + 
  scale_fill_gradient2(low = "#1E88E5", mid = "white", high = "#D81B60") +  
  theme_minimal(base_size = 12) +
  labs(x = "", y = "", fill = "Value") +  
  theme(  
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "none",
    text = element_text(color = "black")
  )
barplot_plot = ML_importance %>%  
  mutate(feature = factor(feature, levels = sorted_features)) %>%  
  ggplot(aes(x = feature, y = importance)) +  
  geom_col(fill = "#9b95c9", width = 0.8) +  
  coord_flip() +  
  theme_minimal(base_size = 12) +
  labs(x = "", y = "Importance") +  
  theme(  
    axis.text.y = element_blank(),  
    axis.text.x = element_text(color = "black"),
    legend.position = "none",
    text = element_text(color = "black")
  )  
combined_plot = heatmap_plot + barplot_plot +   
  plot_layout(widths = c(2, 1))  
print(combined_plot)