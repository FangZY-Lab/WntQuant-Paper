# ==============================================================================
# WntQuant revision -- Pan-cancer Wnt score comparison
# ------------------------------------------------------------------------------
# Purpose:    Compare Wnt activity across the TCGA and TARGET pan-cancer
#             cohorts and correlate it with immune, stemness and Hallmark
#             programmes.
# Inputs:     pancancer_WNT_Score.Rdata, pancancer_group.Rdata and the
#             Hallmark matrices under input/.
# Outputs:    Boxplots, correlation heatmaps and the CD274 relationship
#             plots.
# ==============================================================================

################################################################################Comparison of pan-cancer WNT scores across TCGA####
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
WNT_Score_group=merge(pancancer_WNT_Score$final_activity_score,pancancer_group,by="row.names",all=T)
p = ggplot(WNT_Score_group, aes(Project, pos_score, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +
  labs(x = " ", y = "Wnt/β-catenin pathway activation score") + #activity, activation, inhibition
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(linewidth = 1),
    panel.grid.major.y = element_line(linewidth = 1),
    panel.grid.minor.x = element_line(linewidth = 0.8),
    panel.grid.minor.y = element_line(linewidth = 0.8),
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),
    axis.line.y = element_line(colour = 'black', linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 16),  # Add this line to set y-axis text size
    axis.title.x = element_text(colour = "black", size = 16),  # Add this line to set x-axis title size
    axis.title.y = element_text(angle = 90, colour = "black", size = 16)  # Update to set y-axis title size
  ) +
  scale_fill_manual(values = c("#90d7ec", "#f58f98")) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
################################################################################Matching TCGA and TARGET pan-cancer WNT scores####
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
library(limma)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_score.Rdata")
pancancer_WNT_score=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
WNT_Score_group=merge(pancancer_WNT_score,pancancer_group,by="row.names",all=T)
WNT_Score_group_normal=WNT_Score_group[WNT_Score_group$Group=="Normal",]
WNT_Score_group_normal=WNT_Score_group_normal[,c(1,2,5)]
WNT_Score_group_normal$Row.names=ifelse(substr(WNT_Score_group_normal$Row.names,1,4)=="TCGA",
                                        substr(WNT_Score_group_normal$Row.names, 1, 12),
                                        substr(WNT_Score_group_normal$Row.names,1,(nchar(WNT_Score_group_normal$Row.names)-4)))
WNT_Score_group_normal=as.matrix(WNT_Score_group_normal)
rownames(WNT_Score_group_normal)=WNT_Score_group_normal[,1]
WNT_Score_group_normal=WNT_Score_group_normal[,-1]
WNT_Score_group_normal=as.data.frame(avereps(WNT_Score_group_normal))
WNT_Score_group_tumor=WNT_Score_group[WNT_Score_group$Group=="Tumor",]
WNT_Score_group_tumor=WNT_Score_group_tumor[,c(1,2,5)]
WNT_Score_group_tumor$Row.names=ifelse(substr(WNT_Score_group_tumor$Row.names,1,4)=="TCGA",
                                       substr(WNT_Score_group_tumor$Row.names, 1, 12),
                                       substr(WNT_Score_group_tumor$Row.names,1,(nchar(WNT_Score_group_tumor$Row.names)-4)))
WNT_Score_group_tumor=as.matrix(WNT_Score_group_tumor)
rownames(WNT_Score_group_tumor)=WNT_Score_group_tumor[,1]
WNT_Score_group_tumor=WNT_Score_group_tumor[,-1]
WNT_Score_group_tumor=as.data.frame(avereps(WNT_Score_group_tumor))
WNT_Score_group_tumor$ID=rownames(WNT_Score_group_tumor)
WNT_Score_group_normal$ID=rownames(WNT_Score_group_normal)
rownames(WNT_Score_group_tumor)=NULL
rownames(WNT_Score_group_normal)=NULL
WNT_Score_group_tumor$type="Tumor"
WNT_Score_group_normal$type="Normal"
pairs_data=rbind(WNT_Score_group_tumor,WNT_Score_group_normal)
library(dplyr)
pairs_data <- pairs_data %>%
  group_by(ID) %>%
  filter(n() > 1)
pairs_data$activity_score=as.numeric(pairs_data$activity_score)
pairs_data=pairs_data[!pairs_data$Project%in%c("SARC","SKCM","THYM"),]
mytheme <- theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
                 axis.line = element_line(color = "black", linewidth = 0.4),
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_line(size = 0.1, color = "#e5e5e5"),
                 axis.text.y = element_text(color = "black", size = 16),
                 axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
                 axis.title.y = element_text(color = "black", size = 16),
                 axis.title.x = element_text(color = "black", size = 16),
                 legend.position = "none",
                 panel.spacing = unit(0, "lines"),
                 plot.title = element_text(color = "black", size = 16))
fill.color = c(Normal = "#90d7ec", Tumor = "#f58f98")
p = ggplot(pairs_data, aes(x = type, y = activity_score, fill = type)) +
  geom_line(aes(group = ID), position = position_dodge(0.2), color = "grey90") +
  geom_point(aes(group = ID, size = activity_score), pch = 21,
             size = 1.8, position = position_jitter(w = 0.1)) +
  stat_summary(fun.data = 'mean_se', geom = "errorbar", color = "red",
               width = 0.25, size = 0.8, position = position_dodge(0.9)) +
  facet_wrap(. ~ Project, scales = "free_y", ncol = 10) +
  scale_fill_manual(values = fill.color) +
  labs(x = NULL, y = "Wnt/β-catenin pathway activity score") +
  stat_compare_means(aes(group = type),
                     paired = TRUE,
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 3) +
  theme_bw() + mytheme +
  ggtitle("TCGA")
p
################################################################################Changes in Wnt-related genes####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
wntgene=read.delim("WNT_GENES.txt",header=TRUE,sep='\t',check.names=F)

load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
wntgene_exp=pancancer_exp[intersect(wntgene$WNT_GENES,rownames(pancancer_exp)),]
wntgene_exp=as.data.frame(t(wntgene_exp))


load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
pancancer_WNT_Score=pancancer_WNT_Score[["final_activity_score"]]
pancancer_wnt=pancancer_WNT_Score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]

fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tme_all_fs=wntgene_exp[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tme_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tme_all_fs[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[i],tme=cycle[j],cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_tme_all=rbind(cancer_tme_all,cancer_tme)
  }
  cancer_cor=rbind(cancer_cor,cancer_tme_all)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "",
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*",
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
library(ggplot2)
library(dplyr)
cancer_cor$tme=factor(cancer_cor$tme, levels = colnames(tme_all_fs))
ggplot(cancer_cor, aes(cancer, tme)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 6) + # Adjusted size for text
  theme_minimal() + # Use minimal theme
  theme(
    axis.title.x = element_blank(), # Remove x-axis title
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Adjust x-axis text size
    axis.text.y = element_text(size = 16) # Adjust y-axis text size
  ) +
  # Adjust legend
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +
  theme(
    legend.text = element_text(size = 16), # Set legend text size
    legend.title = element_text(size = 16)  # Set legend title size (if applicable)
  )
################################################################################Validation of the gene sets in pan-cancer####
rm(list=ls())
gc()
library(limma)
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
pancancer_wnt=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
pancancer_wnt=pancancer_wnt[rownames(pancancer_group),]
rm(pancancer_WNT_Score)
vector=unique(pancancer_group$Project)
limma.one.sided = function(fit, lower = FALSE){
  se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total=fit$df.prior + fit$df.residual
  pt(fit$t, df = df.total, lower.tail = lower)
}
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(H-L, levels=design)
  fit=lmFit(gene, design)
  fit=contrasts.fit(fit, compare)
  fit=eBayes(fit)
  if(alternative_limma=="greater"){
    result=limma.one.sided(fit, lower=FALSE)
    result=as.data.frame(result)
    colnames(result)="P.Value"
    result$gene_id=rownames(result)
    topTable=topTable(fit, coef=1, number=Inf)
    topTable$gene_id=rownames(topTable)
    topTable=topTable[,c("gene_id","logFC")]
    result=merge(result,topTable,by="gene_id")
  }else if(alternative_limma=="less"){
    result=limma.one.sided(fit, lower=TRUE)
    result=as.data.frame(result)
    colnames(result)="P.Value"
    result$gene_id=rownames(result)
    topTable=topTable(fit, coef=1, number=Inf)
    topTable$gene_id=rownames(topTable)
    topTable=topTable[,c("gene_id","logFC")]
    result=merge(result,topTable,by="gene_id")
  }else if(alternative_limma=="two.sided"){
    result=topTable(fit, coef=1, number=Inf)
    result=na.omit(result)
    result$gene_id=rownames(result)
  }else{
    print("Please enter a correct method.(greater/less/two.sided)")
  }
  return(result)
}
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
for (i in 1:length(vector)){
  pancancer_group_fs=pancancer_group[pancancer_group$Project==vector[i],]
  exp=pancancer_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  pancancer_group_fs=pancancer_group[pancancer_group$Project==vector[i],]
  group=pancancer_wnt[rownames(pancancer_group_fs),]
  library(dplyr)
  group = group %>%
    mutate(group = ntile(activity_score,2))
  group$group <- factor(group$group, levels =1:2, labels = c("L","H"))
  group$Tag=rownames(group)
  group=group[,c("Tag","group")]
  result=Get_limma(gene=exp,
                   group=group,
                   alternative_limma="greater")
  assign(paste0(vector[i],"_result"),result)
}
geneSets=getGmt("wpags_wpigs.gmt")
library(fgsea)
for (i in 1:length(vector)){
  alldiff=get(paste0(vector[i],"_result"))
  alldiff=alldiff[order(alldiff$logFC,decreasing = T),]
  id=alldiff$logFC
  names(id)=alldiff$gene_id
  canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {
    data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)
  }))
  canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
  canonical_pathways$term=as.character(canonical_pathways$term)
  canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
  fgseaRes=fgsea(pathways=canonical_pathways,
                 stats=id,
                 minSize=1,
                 maxSize=10000,
                 nperm=1000)
  fgsea_result=data.frame(fgseaRes)
  fgsea_result=fgsea_result[,c("pathway","NES")]
  colnames(fgsea_result)[2]=paste0(vector[i])
  if(i==1){
    fgsea_all=fgsea_result
  }else{
    fgsea_all=merge(fgsea_all,fgsea_result,by="pathway",all=T)
  }
}
rownames(fgsea_all)=fgsea_all[,1]
fgsea_all=fgsea_all[,-1]
fgsea_all=as.data.frame(t(fgsea_all))
library(ComplexHeatmap)
library(circlize)
ha = rowAnnotation(Group = factor(c(rep("WPAGSs", 36), rep("WPIGSs", 27))),
                   col = list(Group = c(WPAGSs = "#f15b6c", WPIGSs = "#009ad6")),
                   show_annotation_name = FALSE)
fgsea_all_transposed = t(as.matrix(fgsea_all))
Heatmap(fgsea_all_transposed,
        name = "NES",
        column_title = "Wnt/β-catenin pathway activity score",  # Set the title correctly
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),
                         c("#90d7ec", "white", "#f58f98")),
        left_annotation = ha,
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
        heatmap_legend_param = list(
          title = "NES",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 12)
        ),
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
################################################################################Validating the gene sets with marker genes####
rm(list=ls())
gc()
library(limma)
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
gene_group=pancancer_exp[c("AXIN2","ZNRF3","RNF43","LGR5","TCF7"),]
gene_group=as.data.frame(t(gene_group))
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
vector=unique(pancancer_group$Project)
limma.one.sided = function(fit, lower = FALSE){
  se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total=fit$df.prior + fit$df.residual
  pt(fit$t, df = df.total, lower.tail = lower)
}
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(H-L, levels=design)
  fit=lmFit(gene, design)
  fit=contrasts.fit(fit, compare)
  fit=eBayes(fit)
  if(alternative_limma=="greater"){
    result=limma.one.sided(fit, lower=FALSE)
    result=as.data.frame(result)
    colnames(result)="P.Value"
    result$gene_id=rownames(result)
    topTable=topTable(fit, coef=1, number=Inf)
    topTable$gene_id=rownames(topTable)
    topTable=topTable[,c("gene_id","logFC")]
    result=merge(result,topTable,by="gene_id")
  }else if(alternative_limma=="less"){
    result=limma.one.sided(fit, lower=TRUE)
    result=as.data.frame(result)
    colnames(result)="P.Value"
    result$gene_id=rownames(result)
    topTable=topTable(fit, coef=1, number=Inf)
    topTable$gene_id=rownames(topTable)
    topTable=topTable[,c("gene_id","logFC")]
    result=merge(result,topTable,by="gene_id")
  }else if(alternative_limma=="two.sided"){
    result=topTable(fit, coef=1, number=Inf)
    result=na.omit(result)
    result$gene_id=rownames(result)
  }else{
    print("Please enter a correct method.(greater/less/two.sided)")
  }
  return(result)
}
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
for (i in 1:length(vector)){
  pancancer_group_fs=pancancer_group[pancancer_group$Project==vector[i],]
  exp=pancancer_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  pancancer_group_fs=pancancer_group[pancancer_group$Project==vector[i],]
  group=gene_group[rownames(pancancer_group_fs),]
  group=group[,"AXIN2",drop=F]# Replacement point
  library(dplyr)
  group = group %>%
    mutate(group = ntile(AXIN2,2))# Replacement point
  group$group <- factor(group$group, levels =1:2, labels = c("L","H"))
  group$Tag=rownames(group)
  group=group[,c("Tag","group")]
  identical(group$Tag,colnames(exp))
  result=Get_limma(gene=exp,
                   group=group,
                   alternative_limma="greater")
  assign(paste0(vector[i],"_result"),result)
}
geneSets=getGmt("wpags_wpigs.gmt")
library(fgsea)
for (i in 1:length(vector)){
  alldiff=get(paste0(vector[i],"_result"))
  alldiff=alldiff[order(alldiff$logFC,decreasing = T),]
  id=alldiff$logFC
  names(id)=alldiff$gene_id
  canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {
    data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)
  }))
  canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
  canonical_pathways$term=as.character(canonical_pathways$term)
  canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
  fgseaRes=fgsea(pathways=canonical_pathways,
                 stats=id,
                 minSize=1,
                 maxSize=10000,
                 nperm=1000)
  fgsea_result=data.frame(fgseaRes)
  fgsea_result=fgsea_result[,c("pathway","NES")]
  colnames(fgsea_result)[2]=paste0(vector[i])
  if(i==1){
    fgsea_all=fgsea_result
  }else{
    fgsea_all=merge(fgsea_all,fgsea_result,by="pathway",all=T)
  }
}
rownames(fgsea_all)=fgsea_all[,1]
fgsea_all=fgsea_all[,-1]
fgsea_all=as.data.frame(t(fgsea_all))
library(ComplexHeatmap)
library(circlize)
ha = rowAnnotation(Group = factor(c(rep("WPAGSs", 36), rep("WPIGSs", 27))),
                   col = list(Group = c(WPAGSs = "#f15b6c", WPIGSs = "#009ad6")),
                   show_annotation_name = FALSE)
fgsea_all_transposed = t(as.matrix(fgsea_all))
Heatmap(fgsea_all_transposed,
        name = "NES",
        column_title = "AXIN2",  # Set the title correctly
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),
                         c("#90d7ec", "white", "#f58f98")),
        left_annotation = ha,
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
        heatmap_legend_param = list(
          title = "NES",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 12)
        ),
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
################################################################################Hallmark correlation analysis####
rm(list=ls())
gc()
# load("C:/Users/.../input/tme_all_hallmarker.Rdata")
# load("C:/Users/.../input/TARGET_all_hallmarker.Rdata")
load("C:/Users/user/Desktop/input/hallmarker_activity.Rdata")
tme_all_hallmarker=hallmarker_activity
tme_all=as.data.frame(tme_all_hallmarker)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
pancancer_wnt=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  tme_all_fs=tme_all[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tme_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tme_all_fs[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[i],tme=cycle[j],cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_tme_all=rbind(cancer_tme_all,cancer_tme)
  }
  cancer_cor=rbind(cancer_cor,cancer_tme_all)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "",
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*",
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
library(ggplot2)
library(dplyr)
cancer_cor$tme=factor(cancer_cor$tme, levels = colnames(tme_all_fs))
cancer_cor=cancer_cor[cancer_cor$tme!="WNT_β_catenin_signaling",]
ggplot(cancer_cor, aes(cancer, tme)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 6) + # Adjusted size for text
  theme_minimal() + # Use minimal theme
  theme(
    axis.title.x = element_blank(), # Remove x-axis title
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16), # Adjust x-axis text size
    axis.text.y = element_text(size = 16) # Adjust y-axis text size
  ) +
  # Adjust legend
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +
  theme(
    legend.text = element_text(size = 16), # Set legend text size
    legend.title = element_text(size = 16)  # Set legend title size (if applicable)
  )
################################################################################Relationship between positive and negative cases across TCGA pan-cancer####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
group=pancancer_group[pancancer_group$Group=="Tumor",]
WNT_Score_ssGSEA=pancancer_WNT_Score$final_activity_score[rownames(group),]
WNT_Score_ssGSEA=WNT_Score_ssGSEA[,-1]
wnt_group=merge(WNT_Score_ssGSEA,group,by="row.names")
wnt_group=wnt_group[,-c(1,5)]
wnt_group_avg = aggregate(cbind(pos_score, neg_score) ~ Project, data = wnt_group, FUN = mean)
colnames(wnt_group_avg)[2:3] = c("avg_pos_score", "avg_neg_score")
library(ggplot2)
library(ggpubr)
library(ggrepel)
ggplot(wnt_group_avg, aes(x = avg_pos_score, y = avg_neg_score)) +
  geom_point(size = 2.8, shape = 21, alpha = 1, color = "#90d7d8", fill = "#90d7d8") +
  geom_text_repel(data = wnt_group_avg, aes(label = Project), size = 3.8,
                  box.padding = unit(0.8, 'lines'), show.legend = FALSE, max.overlaps = 500) +
  geom_abline(slope = -1, intercept = 0, color = "#bed742", linewidth = 2, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "#f58f98", linewidth = 2, alpha = 0.5) +
  geom_vline(xintercept = 0, color = "#d3d7d4", linewidth = 0.5, alpha = 0.3) +
  geom_hline(yintercept = 0, color = "#d3d7d4", linewidth = 0.5, alpha = 0.3) +
  labs(x = 'Wnt/β-catenin pathway activation score (Mean)',
       y = 'Wnt/β-catenin pathway inhibition score (Mean)') +
  theme_bw() +
  theme(
    text = element_text(color = "black", size = 16), # Set text color and size
    axis.title = element_text(size = 16),           # Axis title size
    axis.text = element_text(size = 16),             # Axis text size
    plot.title = element_text(size = 16)             # Title size if you add a plot title
  )
library(reshape2)
wnt_group=melt(wnt_group)
wnt_group$variable=ifelse(wnt_group$variable=="pos_score","Activation","Inhibition")
colnames(wnt_group)=c("Project","Group","value")
p = ggplot(wnt_group, aes(Project, value, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +
  labs(x = " ", y = "Wnt/β-catenin pathway activation/inhibition score") +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(linewidth = 1),
    panel.grid.major.y = element_line(linewidth = 1),
    panel.grid.minor.x = element_line(linewidth = 0.8),
    panel.grid.minor.y = element_line(linewidth = 0.8),
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),
    axis.line.y = element_line(colour = 'black', linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16),  # Set font size for x-axis text
    axis.text.y = element_text(colour = "black", size = 16),         # Set font size for y-axis text
    axis.title.y = element_text(angle = 90, colour = "black", size = 16),  # Set font size and color for y-axis title
    axis.title.x = element_text(colour = "black", size = 16)         # Set font size and color for x-axis title
  ) +
  scale_fill_manual(values = c("#f58f98", "#90d7ec")) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test")
p
################################################################################WNT-immune infiltration correlation analysis####
rm(list=ls())
gc()
library(IOBR)
load("C:/Users/user/Desktop/input/tme_all_literatures.Rdata")
tme_all=as.data.frame(tme_all_literatures)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
sig=signature_collection_citation[c(19:54,127:141,282:314),]
tme_all=tme_all[,intersect(sig$Signatures,colnames(tme_all))]
sig=sig[sig$Signatures%in%colnames(tme_all),]
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
pancancer_wnt=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tme_all_fs=tme_all[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tme_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity,tme_all_fs[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[i],tme=cycle[j],cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_tme_all=rbind(cancer_tme_all,cancer_tme)
  }
  cancer_cor=rbind(cancer_cor,cancer_tme_all)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "",
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*",
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
library(ggplot2)
library(dplyr)
cancer_cor$tme=factor(cancer_cor$tme, levels = colnames(tme_all_fs))
ggplot(cancer_cor, aes(cancer, tme)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 6) + # Adjusted size for text
  theme_minimal() + # Use minimal theme
  theme(
    axis.title.x = element_blank(), # Remove x-axis title
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16), # Adjust x-axis text size
    axis.text.y = element_text(size = 16) # Adjust y-axis text size
  ) +
  # Adjust legend
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +
  theme(
    legend.text = element_text(size = 16), # Set legend text size
    legend.title = element_text(size = 16)  # Set legend title size (if applicable)
  )
################################################################################WNT-immunomodulator/chemokine correlation analysis####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
Chemokinegene=read.delim("Chemokine.txt",header=TRUE,sep='\t',check.names=F)#Chemokine
colnames(Chemokinegene)="gene"
Immunomodulatorgene=read.delim("Immunomodulator.txt",header=TRUE,sep='\t',check.names=F)#Immunomodulator
colnames(Immunomodulatorgene)="gene"
Chemokine_Immunomodulator=rbind(Immunomodulatorgene,Chemokinegene)
Chemokine_Immunomodulator=as.data.frame(unique(Chemokine_Immunomodulator$gene))
Chemokine_Immunomodulator=Chemokine_Immunomodulator[order(Chemokine_Immunomodulator$`unique(Chemokine_Immunomodulator$gene)`),,drop=F]
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
wntgene_exp=pancancer_exp[intersect(Chemokine_Immunomodulator$`unique(Chemokine_Immunomodulator$gene)`,rownames(pancancer_exp)),]
wntgene_exp=as.data.frame(t(wntgene_exp))
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
pancancer_wnt=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tme_all_fs=wntgene_exp[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tme_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tme_all_fs[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[i],tme=cycle[j],cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_tme_all=rbind(cancer_tme_all,cancer_tme)
  }
  cancer_cor=rbind(cancer_cor,cancer_tme_all)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "",
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*",
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
library(ggplot2)
library(dplyr)
cancer_cor$tme=factor(cancer_cor$tme, levels = colnames(tme_all_fs))
ggplot(cancer_cor, aes(cancer, tme)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 6) + # Adjusted size for text
  theme_minimal() + # Use minimal theme
  theme(
    axis.title.x = element_blank(), # Remove x-axis title
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16), # Adjust x-axis text size
    axis.text.y = element_text(size = 16) # Adjust y-axis text size
  ) +
  # Adjust legend
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +
  theme(
    legend.text = element_text(size = 16), # Set legend text size
    legend.title = element_text(size = 16)  # Set legend title size (if applicable)
  )
################################################################################Stemness correlation analysis####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/tme_all_stemness.Rdata")
tme_all=as.data.frame(tme_all_stemness)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
pancancer_wnt=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  tme_all_fs=tme_all[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tme_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tme_all_fs[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[i],tme=cycle[j],cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_tme_all=rbind(cancer_tme_all,cancer_tme)
  }
  cancer_cor=rbind(cancer_cor,cancer_tme_all)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "",
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*",
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
library(ggplot2)
library(dplyr)
cancer_cor$tme=factor(cancer_cor$tme, levels = colnames(tme_all_fs))
ggplot(cancer_cor, aes(cancer, tme)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 6) + # Adjusted size for text
  theme_minimal() + # Use minimal theme
  theme(
    axis.title.x = element_blank(), # Remove x-axis title
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.title.y = element_blank(), # Remove y-axis title
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16), # Adjust x-axis text size
    axis.text.y = element_text(size = 16) # Adjust y-axis text size
  ) +
  # Adjust legend
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +
  theme(
    legend.text = element_text(size = 16), # Set legend text size
    legend.title = element_text(size = 16)  # Set legend title size (if applicable)
  )
################################################################################Relationship between CD274 and Wnt####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
rt1=read.delim("EXP0011_PCG_log2_afterQC.txt")[-1,]
rt2=read.delim("EXP0012_PCG_log2_afterQC.txt")[-1,]
rt3=read.delim("EXP0051_PCG_log2_afterQC.txt")[-1,]
################################################################################BCR/TCR####
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/hippo")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
pancancer_cluster$Tag = substr(pancancer_cluster$Tag, 1, 12)
library(readxl)
geneset = readxl::read_excel("TCR.xlsx")
geneset = geneset[, c(1, 21:26)]
colnames(geneset)[1] = "Tag"
geneset_hippo = merge(geneset, pancancer_cluster, by = "Tag")
# Shannon
geneset_hippo1 = geneset_hippo[, c(3, 8)]
geneset_hippo1 = na.omit(geneset_hippo1)

# TCR Shannon
geneset_hippo2 = geneset_hippo[, c(5, 8)]
geneset_hippo2 = na.omit(geneset_hippo2)

# BCR Evenness
geneset_hippo3 = geneset_hippo[, c(2, 8)]
geneset_hippo3 = na.omit(geneset_hippo3)

# TCR Evenness
geneset_hippo4 = geneset_hippo[, c(7, 8)]
geneset_hippo4 = na.omit(geneset_hippo4)

# BCR Richness
geneset_hippo5 = geneset_hippo[, c(4, 8)]
geneset_hippo5 = na.omit(geneset_hippo5)

# TCR Richness
geneset_hippo6 = geneset_hippo[, c(6, 8)]
geneset_hippo6 = na.omit(geneset_hippo6)
library(ggplot2)
library(ggpubr)
library(patchwork)
# 8.1 BCR Shannon
p1 <- ggplot(geneset_hippo1, aes(x = cluster, y = `BCR Shannon`, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "BCR Shannon", x = "", y = "Shannon index") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "plain")) +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  scale_fill_manual(values = c("#d71345", "#f47920", "#145b7d", "#6950a1")) +
  coord_flip(ylim = c(0, 6))  # Widen the range so the boxes stand out

# 8.2 TCR Shannon
p2 <- ggplot(geneset_hippo2, aes(x = cluster, y = `TCR Shannon`, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "TCR Shannon", x = "", y = "Shannon index") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "plain")) +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  scale_fill_manual(values = c("#d71345", "#f47920", "#145b7d", "#6950a1")) +
  coord_flip(ylim = c(0, 6))

# 8.3 BCR Evenness
p3 <- ggplot(geneset_hippo3, aes(x = cluster, y = `BCR Evenness`, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "BCR Evenness", x = "", y = "Evenness index") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "plain")) +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  scale_fill_manual(values = c("#d71345", "#f47920", "#145b7d", "#6950a1")) +
  coord_flip(ylim = c(0.5, 1))

# 8.4 TCR Evenness
p4 <- ggplot(geneset_hippo4, aes(x = cluster, y = `TCR Evenness`, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "TCR Evenness", x = "", y = "Evenness index") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "plain")) +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  scale_fill_manual(values = c("#d71345", "#f47920", "#145b7d", "#6950a1")) +
  coord_flip(ylim = c(0.9, 1))

# 8.5 BCR Richness
p5 <- ggplot(geneset_hippo5, aes(x = cluster, y = `BCR Richness`, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "BCR Richness", x = "", y = "Richness index") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "plain")) +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  scale_fill_manual(values = c("#d71345", "#f47920", "#145b7d", "#6950a1")) +
  coord_flip(ylim = c(0, 100))

# 8.6 TCR Richness
p6 <- ggplot(geneset_hippo6, aes(x = cluster, y = `TCR Richness`, fill = cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  labs(title = "TCR Richness", x = "", y = "Richness index") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(color = "black", size = 14),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(color = "black", size = 14),
        plot.title = element_text(color = "black", size = 14, face = "plain")) +
  stat_compare_means(method = "kruskal.test", label = "p.signif") +
  scale_fill_manual(values = c("#d71345", "#f47920", "#145b7d", "#6950a1")) +
  coord_flip(ylim = c(0, 100))
combined_plot <- (p1 | p3 | p5) / (p2 | p4 | p6)
print(combined_plot)
