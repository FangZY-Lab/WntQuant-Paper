# ==============================================================================
# WntQuant-Paper3 -- Initial pan-cancer explorations
# ------------------------------------------------------------------------------
# Purpose:    Compute Wnt activity scores across TCGA and TARGET pan-cancer atlases and perform initial cross-cancer comparisons.
# Inputs:     Pan-cancer expression matrices (TCGA / TARGET) and annotations.
# Outputs:    Pan-cancer Wnt activity score objects.
# NOTE:       Part of the code involves the core interests of the laboratory and is not open source at this time.
# ==============================================================================


rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/pancancer_exp.Rdata")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=pancancer_exp,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
pancancer_WNT_Score=WNT_Score
save(pancancer_WNT_Score,file="pancancer_WNT_Score.Rdata")
# ------------------------------------------------------------------------------
# TARGET pan-cancer WNT score calculation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/TARGET.Rdata")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=TARGET,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
TARGET_WNT_Score=WNT_Score
save(TARGET_WNT_Score,file="TARGET_WNT_Score.Rdata")
# ------------------------------------------------------------------------------
# Comparison of TCGA pan-cancer WNT scores
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
WNT_Score_group=merge(pancancer_WNT_Score$final_activity_score,pancancer_group,by="row.names",all=T)
p = ggplot(WNT_Score_group, aes(Project, activity_score, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") + #activity, activation, inhibition  
  theme_bw() +  
  theme(  
    legend.position = "right",  
    panel.grid.major.x = element_line(linewidth = 1),  
    panel.grid.major.y = element_line(linewidth = 1),  
    panel.grid.minor.x = element_line(linewidth = 0.8),  
    panel.grid.minor.y = element_line(linewidth = 0.8),  
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
    axis.line.y = element_line(colour = 'black', linewidth = 1),  
    axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black", size = 16),  
    axis.text.y = element_text(colour = "black", size = 16),
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, colour = "black", size = 16)
  ) +  
  scale_fill_manual(values = c("#90d7ec", "#f58f98")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") 
p
# ------------------------------------------------------------------------------
# Comparison of TARGET pan-cancer WNT scores
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
setwd("input")
load("input/TARGET_WNT_Score.Rdata")
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
WNT_Score_group=merge(TARGET_WNT_Score$final_activity_score,TARGET_G,by="row.names",all=T)
p = ggplot(WNT_Score_group, aes(Project, activity_score, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +  #activity, activation, inhibition  
  theme_bw() +  
  theme(  
    legend.position = "right",  
    panel.grid.major.x = element_line(linewidth = 1),  
    panel.grid.major.y = element_line(linewidth = 1),  
    panel.grid.minor.x = element_line(linewidth = 0.8),  
    panel.grid.minor.y = element_line(linewidth = 0.8),  
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
    axis.line.y = element_line(colour = 'black', linewidth = 1),  
    axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 16),
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, colour = "black", size = 16)
  ) +  
  scale_fill_manual(values = c("#90d7ec", "#f58f98")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") 
p
# ------------------------------------------------------------------------------
# Matching comparison of TCGA-TARGET pan-cancer WNT scores
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
library(limma)
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
load("input/TARGET_WNT_Score.Rdata")
pancancer_WNT_score=rbind(pancancer_WNT_Score$final_activity_score,TARGET_WNT_Score$final_activity_score)
load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=rbind(pancancer_group,TARGET_G)

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
pairs_data=pairs_data[!pairs_data$Project%in%c("TCGA_SARC","TCGA_SKCM","TCGA_THYM"),]
mytheme <- theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),  
                 axis.line = element_line(color = "black", linewidth = 0.4),  
                 panel.grid.minor = element_blank(),  
                 panel.grid.major = element_line(size = 0.1, color = "#e5e5e5"),  
                 axis.text.y = element_text(color = "black", size = 16),  
                 axis.text.x = element_text(angle = 60, hjust = 1, color = "black", size = 16),  
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
  ggtitle("TCGA+TARGET") 
p
# ------------------------------------------------------------------------------
# Changes in Wnt-related genes
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
wntgene=read.delim("WNT_GENES.txt",header=TRUE,sep='\t',check.names=F)

load("input/pancancer_exp.Rdata")
wntgene_exp=pancancer_exp[intersect(wntgene$WNT_GENES,rownames(pancancer_exp)),]
wntgene_exp=as.data.frame(t(wntgene_exp))
load("input/TARGET.Rdata")
wntgene_TARGET=TARGET[intersect(wntgene$WNT_GENES,rownames(TARGET)),]
wntgene_TARGET=as.data.frame(t(wntgene_TARGET))
wntgene_exp=rbind(wntgene_exp,wntgene_TARGET)

load("input/pancancer_WNT_Score.Rdata")
pancancer_WNT_Score=pancancer_WNT_Score[["final_activity_score"]]
load("input/TARGET_WNT_Score.Rdata")
TARGET_WNT_Score=TARGET_WNT_Score[["final_activity_score"]]
pancancer_wnt=rbind(pancancer_WNT_Score,TARGET_WNT_Score)

load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
pancancer_group=rbind(pancancer_group,TARGET_G)

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
  geom_text(aes(label = pstar), color = "black", size = 6) +
  theme_minimal() +
  theme(  
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +  
  theme(  
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )  
# ------------------------------------------------------------------------------
# Characterised gene validation gene set
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
setwd("input")
load("input/pancancer_exp.Rdata")
load("input/TARGET.Rdata")
TARGET=TARGET[rownames(pancancer_exp),]
identical(rownames(pancancer_exp),rownames(TARGET))
pancancer_exp=cbind(pancancer_exp,TARGET)
gene_group=pancancer_exp[c("AXIN2","ZNRF3","RNF43","LGR5","TCF7"),]
gene_group=as.data.frame(t(gene_group))
load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
pancancer_group=rbind(pancancer_group,TARGET_G)
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
load("input/pancancer_exp.Rdata")
load("input/TARGET.Rdata")
identical(rownames(pancancer_exp),rownames(TARGET))
pancancer_exp=cbind(pancancer_exp,TARGET)
for (i in 1:length(vector)){
  pancancer_group_fs=pancancer_group[pancancer_group$Project==vector[i],]
  exp=pancancer_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  pancancer_group_fs=pancancer_group[pancancer_group$Project==vector[i],]
  group=gene_group[rownames(pancancer_group_fs),]
  group=group[,"AXIN2",drop=F]#Location of replacement genes
  library(dplyr)
  group = group %>%  
    mutate(group = ntile(AXIN2,2))#Location of replacement genes
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha, 
        show_row_names = F,  
        show_column_names = T, 
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,  
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1) 
)
grid.text("AXIN2", x = 0.1, y = unit(0.97, "npc"), gp = gpar(fontsize = 16, fontface = "bold"))#Location of replacement genes
# ------------------------------------------------------------------------------
# Validation of gene sets in pan-cancer
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
load("input/TARGET_WNT_Score.Rdata")
pancancer_wnt=rbind(pancancer_WNT_Score$final_activity_score,TARGET_WNT_Score$final_activity_score)
load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
pancancer_group=rbind(pancancer_group,TARGET_G)
pancancer_wnt=pancancer_wnt[rownames(pancancer_group),]
rm(TARGET_WNT_score)
rm(TARGET_G)
rm(pancancer_WNT_score)
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
load("input/pancancer_exp.Rdata")
load("input/TARGET.Rdata")
identical(rownames(pancancer_exp),rownames(TARGET))
pancancer_exp=cbind(pancancer_exp,TARGET)
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T,
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
grid.text("", x = 0.2, y = unit(0.99, "npc"), gp = gpar(fontsize = 12, fontface = "bold")) 
# ------------------------------------------------------------------------------
# TCGA hallmark ssGSEA
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
library(GSEABase)
hallmarker=getGmt("geneset_hallmark.gmt")
hallmarker_gmt=lapply(hallmarker, function(x){ x@geneIds })
names(hallmarker_gmt)=names(hallmarker)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_hallmarker=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_hallmarker=calculate_sig_score(eset = pancancer_exp_fs,signature = hallmarker_gmt,method="ssgsea",mini_gene_count=1)
  tme_all_hallmarker=rbind(tme_all_hallmarker,im_hallmarker)
}
save(tme_all_hallmarker,file="tme_all_hallmarker.Rdata")
# ------------------------------------------------------------------------------
# TARGET hallmark ssGSEA
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
library(GSEABase)
hallmarker=getGmt("geneset_hallmark.gmt")
hallmarker_gmt=lapply(hallmarker, function(x){ x@geneIds })
names(hallmarker_gmt)=names(hallmarker)
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_all_hallmarker=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_hallmarker=calculate_sig_score(eset = pancancer_exp_fs,signature = hallmarker_gmt,method="ssgsea",mini_gene_count=1)
  TARGET_all_hallmarker=rbind(TARGET_all_hallmarker,im_hallmarker)
}
save(TARGET_all_hallmarker,file="TARGET_all_hallmarker.Rdata")
# ------------------------------------------------------------------------------
# Hallmark correlation analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/tme_all_hallmarker.Rdata")
load("input/TARGET_all_hallmarker.Rdata")
tme_all=as.data.frame(tme_all_hallmarker)
tme_all=rbind(tme_all,as.data.frame(TARGET_all_hallmarker))
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("input/pancancer_WNT_Score.Rdata")
load("input/TARGET_WNT_Score.Rdata")
pancancer_wnt=rbind(pancancer_WNT_Score$final_activity_score,TARGET_WNT_Score$final_activity_score)
load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
pancancer_group=rbind(pancancer_group,TARGET_G)
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
  geom_text(aes(label = pstar), color = "black", size = 6) +
  theme_minimal() +
  theme(  
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.y = element_blank(), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +  
  theme(  
    legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16)
  )  
# ------------------------------------------------------------------------------
# TCGA pan-cancer pos vs. neg
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/pancancer_WNT_Score.Rdata")
load("input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
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
    text = element_text(color = "black", size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16)
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
    axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, colour = "black", size = 16),
    axis.title.x = element_text(colour = "black", size = 16)
  ) +  
  scale_fill_manual(values = c("#f58f98", "#90d7ec")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") 
p
# ------------------------------------------------------------------------------
# TARGET pan-cancer pos vs neg
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/TARGET_WNT_Score.Rdata")
load("input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
group=TARGET_G[TARGET_G$Group=="Tumor",]
WNT_Score_ssGSEA=TARGET_WNT_Score$final_activity_score[rownames(group),]
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
    text = element_text(color = "black", size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 16)
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
    axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, colour = "black", size = 16),
    axis.title.x = element_text(colour = "black", size = 16)
  ) +  
  scale_fill_manual(values = c("#f58f98", "#90d7ec")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") 
p
# ------------------------------------------------------------------------------
# Calculation of WNT scores for GTEx
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
load("input/GTEx_exp.Rdata")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=GTEx_exp,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
GTEx_WNT_score=WNT_Score$final_activity_score
save(GTEx_WNT_score,file="GTEx_WNT_score.Rdata")
rm(list=ls())
gc()
load("input/GTEx_pancancer_mrna_pheno.rdata")
rownames(gtex_mrna_pheno)=gtex_mrna_pheno[,1]
gtex_mrna_pheno=gtex_mrna_pheno[,-1]
load("input/GTEx_WNT_score.Rdata")
tissue=gtex_mrna_pheno[,1,drop=F]
WNT_Score_tissue=merge(GTEx_WNT_score,tissue,by="row.names")
library(RColorBrewer)
library(ggplot2)
median_scores=aggregate(activity_score ~ primary_site, data = WNT_Score_tissue, FUN = median)
sorted_sites=median_scores$primary_site[order(median_scores$activity_score, decreasing = TRUE)]
colors=colorRampPalette(c("#f58f98", "#90d7ec"))(length(sorted_sites))
WNT_Score_tissue$primary_site=factor(WNT_Score_tissue$primary_site, levels = sorted_sites)
p = ggplot(WNT_Score_tissue, aes(x = primary_site, y = activity_score, fill = primary_site)) +  
  geom_boxplot(color = "black") +  
  theme_bw() +  
  theme(  
    axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16),  
    axis.text.y = element_text(color = "black", size = 16),  
    axis.title = element_text(color = "black", size = 16),  
    plot.title = element_text(color = "black", size = 16),  
    legend.position = 'none'  
  ) +  
  scale_fill_manual(values = colors) +  
  labs(  
    y = "Wnt/β-catenin pathway activity score", #activation/inhibition 
    x = ""  
  ) +   
  ggtitle("GTEx")  
p
# ------------------------------------------------------------------------------
# Validation of gene sets in GTEx
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
library(dplyr)
setwd("input")
load("input/GTEx_pancancer_mrna_pheno.rdata")
rownames(gtex_mrna_pheno)=gtex_mrna_pheno[,1]
gtex_mrna_pheno=gtex_mrna_pheno[,-1]
load("input/GTEx_WNT_score.Rdata")
tissue=gtex_mrna_pheno[,1,drop=F]
WNT_Score_tissue=merge(GTEx_WNT_score,tissue,by="row.names")
vector=unique(WNT_Score_tissue$primary_site)
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
load("input/GTEx_exp.Rdata")
rownames(WNT_Score_tissue)=WNT_Score_tissue[,1]
WNT_Score_tissue=WNT_Score_tissue[,-1]
for (i in 1:length(vector)){
  pancancer_group_fs=WNT_Score_tissue[WNT_Score_tissue$primary_site==vector[i],]
  exp=GTEx_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  group=pancancer_group_fs
  library(dplyr)
  group = group %>%  
    mutate(group = ntile(activity_score,2))
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T,
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
grid.text(" ", x = 0.3, y = unit(1, "npc"), gp = gpar(fontsize = 12, fontface = "bold")) 
# ------------------------------------------------------------------------------
# Calculation of WNT scores for the CCLE
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
load("input/CCLE_exp.Rdata")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=CCLE_exp,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
CCLE_WNT_score=WNT_Score$final_activity_score
save(CCLE_WNT_score,file="CCLE_WNT_score.Rdata")
rm(list=ls())
gc()
load("input/CCLE_WNT_score.Rdata")
CCLE_WNT_score$tissue=rownames(CCLE_WNT_score)
CCLE_WNT_score$tissue=gsub("^[^_]*_", "", CCLE_WNT_score$tissue)
CCLE_WNT_score$tissue=tolower(CCLE_WNT_score$tissue)
CCLE_WNT_score$tissue=gsub("_", " ", CCLE_WNT_score$tissue)
CCLE_WNT_score$tissue <- sapply(strsplit(CCLE_WNT_score$tissue, " "), function(x) {
  sapply(x, function(y) {
    paste(toupper(substring(y, 1, 1)), substring(y, 2), sep = "")
  }) %>%
    paste(collapse = " ")
})
median_scores=aggregate(neg_score ~ tissue, data = CCLE_WNT_score, FUN = median)
sorted_sites=median_scores$tissue[order(median_scores$neg_score, decreasing = TRUE)]
colors=colorRampPalette(c("#f58f98", "#90d7ec"))(length(sorted_sites))
CCLE_WNT_score$tissue=factor(CCLE_WNT_score$tissue, levels = sorted_sites)
p = ggplot(CCLE_WNT_score, aes(x = tissue, y = neg_score, fill = tissue)) +  
  geom_boxplot(color = "black") +  
  theme_bw() +  
  theme(  
    axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16),  
    axis.text.y = element_text(color = "black", size = 16),  
    axis.title = element_text(color = "black", size = 16),  
    plot.title = element_text(color = "black", size = 16),  
    legend.position = 'none'  
  ) +  
  scale_fill_manual(values = colors) +  
  labs(  
    y = "Wnt/β-catenin pathway inhibition score", #activation/inhibition 
    x = ""  
  ) +   
  ggtitle("CCLE")  
p
# ------------------------------------------------------------------------------
# Validation of gene sets in CCLE
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
library(dplyr)
setwd("input")
load("input/CCLE_WNT_score.Rdata")
CCLE_WNT_score$tissue=rownames(CCLE_WNT_score)
CCLE_WNT_score$tissue=gsub("^[^_]*_", "", CCLE_WNT_score$tissue)
CCLE_WNT_score$tissue=tolower(CCLE_WNT_score$tissue)
CCLE_WNT_score$tissue=gsub("_", " ", CCLE_WNT_score$tissue)
CCLE_WNT_score$tissue <- sapply(strsplit(CCLE_WNT_score$tissue, " "), function(x) {
  sapply(x, function(y) {
    paste(toupper(substring(y, 1, 1)), substring(y, 2), sep = "")
  }) %>%
    paste(collapse = " ")
})
ccle_g=as.data.frame(table(CCLE_WNT_score$tissue))
vector=unique(CCLE_WNT_score$tissue)
vector=vector[!vector %in% c("Small Intestine", "Salivary Gland", "Cervix")] 
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
load("input/CCLE_exp.Rdata")
for (i in 1:length(vector)){
  pancancer_group_fs=CCLE_WNT_score[CCLE_WNT_score$tissue==vector[i],]
  exp=CCLE_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  group=pancancer_group_fs
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T,
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
grid.text(" ", x = 0.3, y = unit(1, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
# ------------------------------------------------------------------------------
# Calculation of WNT scores for GDSC version 1
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
GDSC1_exp=readRDS("input/DataFiles/DataFiles/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=as.data.frame(GDSC1_exp),
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
GDSC1_WNT_score=WNT_Score$final_activity_score
save(GDSC1_WNT_score,file="GDSC1_WNT_score.Rdata")
rm(list=ls())
gc()
library(ggplot2)
library(readxl)
setwd("input/DataFiles/DataFiles/GLDS/GDSCv1")
load("input/GDSC1_WNT_score.Rdata")
group=read_excel("Cell_Lines_Details.xlsx")
group=na.omit(group[,c(2,9)])
group$`COSMIC identifier`=paste0("COSMIC_",group$`COSMIC identifier`)
GDSC1_WNT_score=merge(GDSC1_WNT_score,group,by.x="row.names",by.y="COSMIC identifier")
colnames(GDSC1_WNT_score)[5]="tissue"
process_text=function(x) {
  x = gsub(" ", "_", x)
  x = sapply(x, function(s) {
    if (nchar(s) > 0) {
      paste0(toupper(substring(s, 1, 1)), 
             substring(s, 2))
    } else {
      s
    }
  }, USE.NAMES = FALSE)
  return(x)
}
GDSC1_WNT_score$tissue=process_text(GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("Lung", GDSC1_WNT_score$tissue), "Lung", GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("le[ui]k[a]?emia", GDSC1_WNT_score$tissue, ignore.case = TRUE), "Leukaemia", GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("lymphoma", GDSC1_WNT_score$tissue, ignore.case = TRUE), "Lymphoma", GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("neoplasm", GDSC1_WNT_score$tissue, ignore.case = TRUE), "Neoplasm", GDSC1_WNT_score$tissue)
median_scores=aggregate(neg_score ~ tissue, data = GDSC1_WNT_score, FUN = median)
sorted_sites=median_scores$tissue[order(median_scores$neg_score, decreasing = TRUE)]
colors=colorRampPalette(c("#f58f98", "#90d7ec"))(length(sorted_sites))
GDSC1_WNT_score$tissue=factor(GDSC1_WNT_score$tissue, levels = sorted_sites)
p = ggplot(GDSC1_WNT_score, aes(x = tissue, y = neg_score, fill = tissue)) +  
  geom_boxplot(color = "black") +  
  theme_bw() +  
  theme(  
    axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16),  
    axis.text.y = element_text(color = "black", size = 16),  
    axis.title = element_text(color = "black", size = 16),  
    plot.title = element_text(color = "black", size = 16),  
    legend.position = 'none'  
  ) +  
  scale_fill_manual(values = colors) +  
  labs(  
    y = "Wnt/β-catenin pathway inhibition score", #activation/inhibition 
    x = ""  
  ) +   
  ggtitle("GDSC Version1")  
p
# ------------------------------------------------------------------------------
# Validation of gene sets in GDSC version1
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
library(dplyr)
setwd("input/DataFiles/DataFiles/GLDS/GDSCv1")
load("input/GDSC1_WNT_score.Rdata")
group=read_excel("Cell_Lines_Details.xlsx")
group=na.omit(group[,c(2,9)])
group$`COSMIC identifier`=paste0("COSMIC_",group$`COSMIC identifier`)
GDSC1_WNT_score=merge(GDSC1_WNT_score,group,by.x="row.names",by.y="COSMIC identifier")
colnames(GDSC1_WNT_score)[5]="tissue"
process_text=function(x) {
  x = gsub(" ", "_", x)
  x = sapply(x, function(s) {
    if (nchar(s) > 0) {
      paste0(toupper(substring(s, 1, 1)), 
             substring(s, 2))
    } else {
      s
    }
  }, USE.NAMES = FALSE)
  return(x)
}
GDSC1_WNT_score$tissue=process_text(GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("Lung", GDSC1_WNT_score$tissue), "Lung", GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("le[ui]k[a]?emia", GDSC1_WNT_score$tissue, ignore.case = TRUE), "Leukaemia", GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("lymphoma", GDSC1_WNT_score$tissue, ignore.case = TRUE), "Lymphoma", GDSC1_WNT_score$tissue)
GDSC1_WNT_score$tissue=ifelse(grepl("neoplasm", GDSC1_WNT_score$tissue, ignore.case = TRUE), "Neoplasm", GDSC1_WNT_score$tissue)
number=as.data.frame(table(GDSC1_WNT_score$tissue))
number=number[number$Freq>=6,]
vector=unique(number$Var1)
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
rownames(GDSC1_WNT_score)=GDSC1_WNT_score[,1]
GDSC1_WNT_score=GDSC1_WNT_score[,-1]
GDSC1_exp=readRDS("input/DataFiles/DataFiles/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
for (i in 1:length(vector)){
  pancancer_group_fs=GDSC1_WNT_score[GDSC1_WNT_score$tissue==vector[i],]
  exp=GDSC1_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  group=pancancer_group_fs
  library(dplyr)
  group = group %>%  
    mutate(group = ntile(activity_score,2))
  group$group = factor(group$group, levels =1:2, labels = c("L","H"))
  group$Tag=rownames(group)
  group=group[,c("Tag","group")]
  result=Get_limma(gene=exp,
                   group=group,
                   alternative_limma="greater")
  assign(paste0(vector[i],"_result"),result)
}
setwd("input")
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T,
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
grid.text(" ", x = 0.3, y = unit(1, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
# ------------------------------------------------------------------------------
# Calculation of WNT scores for GDSC version 2
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
GDSC2_exp=readRDS("input/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=as.data.frame(GDSC2_exp),
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
GDSC2_WNT_score=WNT_Score$final_activity_score
save(GDSC2_WNT_score,file="GDSC2_WNT_score.Rdata")
rm(list=ls())
gc()
library(readxl)
setwd("input/DataFiles/DataFiles/GLDS/GDSCv2")
load("input/GDSC2_WNT_score.Rdata")
group=read_excel("Cell_Lines_Details.xlsx")
group=na.omit(group[,c(2,9)])
group$`COSMIC identifier`=paste0("COSMIC_",group$`COSMIC identifier`)
GDSC2_WNT_score=merge(GDSC2_WNT_score,group,by.x="row.names",by.y="COSMIC identifier")
colnames(GDSC2_WNT_score)[5]="tissue"
process_text=function(x) {
  x <- gsub(" ", "_", x)
  x <- sapply(x, function(s) {
    if (nchar(s) > 0) {
      paste0(toupper(substring(s, 1, 1)), 
             substring(s, 2))
    } else {
      s
    }
  }, USE.NAMES = FALSE)
  return(x)
}
GDSC2_WNT_score$tissue=process_text(GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("Lung", GDSC2_WNT_score$tissue), "Lung", GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("le[ui]k[a]?emia", GDSC2_WNT_score$tissue, ignore.case = TRUE), "Leukaemia", GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("lymphoma", GDSC2_WNT_score$tissue, ignore.case = TRUE), "Lymphoma", GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("neoplasm", GDSC2_WNT_score$tissue, ignore.case = TRUE), "Neoplasm", GDSC2_WNT_score$tissue)
median_scores=aggregate(neg_score ~ tissue, data = GDSC2_WNT_score, FUN = median)
sorted_sites=median_scores$tissue[order(median_scores$neg_score, decreasing = TRUE)]
colors=colorRampPalette(c("#f58f98", "#90d7ec"))(length(sorted_sites))
GDSC2_WNT_score$tissue=factor(GDSC2_WNT_score$tissue, levels = sorted_sites)
p = ggplot(GDSC2_WNT_score, aes(x = tissue, y = neg_score, fill = tissue)) +  
  geom_boxplot(color = "black") +  
  theme_bw() +  
  theme(  
    axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16),  
    axis.text.y = element_text(color = "black", size = 16),  
    axis.title = element_text(color = "black", size = 16),  
    plot.title = element_text(color = "black", size = 16),  
    legend.position = 'none'  
  ) +  
  scale_fill_manual(values = colors) +  
  labs(  
    y = "Wnt/β-catenin pathway inhibition score", #activation/inhibition 
    x = ""  
  ) +   
  ggtitle("GDSC Version2")  
p
# ------------------------------------------------------------------------------
# Validation of gene sets in GDSC version2
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
library(dplyr)
setwd("input/DataFiles/DataFiles/GLDS/GDSCv2")
load("input/GDSC2_WNT_score.Rdata")
group=read_excel("Cell_Lines_Details.xlsx")
group=na.omit(group[,c(2,9)])
group$`COSMIC identifier`=paste0("COSMIC_",group$`COSMIC identifier`)
GDSC2_WNT_score=merge(GDSC2_WNT_score,group,by.x="row.names",by.y="COSMIC identifier")
colnames(GDSC2_WNT_score)[5]="tissue"
process_text=function(x) {
  x = gsub(" ", "_", x)
  x = sapply(x, function(s) {
    if (nchar(s) > 0) {
      paste0(toupper(substring(s, 1, 1)), 
             substring(s, 2))
    } else {
      s
    }
  }, USE.NAMES = FALSE)
  return(x)
}
GDSC2_WNT_score$tissue=process_text(GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("Lung", GDSC2_WNT_score$tissue), "Lung", GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("le[ui]k[a]?emia", GDSC2_WNT_score$tissue, ignore.case = TRUE), "Leukaemia", GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("lymphoma", GDSC2_WNT_score$tissue, ignore.case = TRUE), "Lymphoma", GDSC2_WNT_score$tissue)
GDSC2_WNT_score$tissue=ifelse(grepl("neoplasm", GDSC2_WNT_score$tissue, ignore.case = TRUE), "Neoplasm", GDSC2_WNT_score$tissue)
number=as.data.frame(table(GDSC2_WNT_score$tissue))
number=number[number$Freq>=6,]
vector=unique(number$Var1)
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
rownames(GDSC2_WNT_score)=GDSC2_WNT_score[,1]
GDSC2_WNT_score=GDSC2_WNT_score[,-1]
GDSC2_exp=readRDS("input/DataFiles/DataFiles/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
for (i in 1:length(vector)){
  pancancer_group_fs=GDSC2_WNT_score[GDSC2_WNT_score$tissue==vector[i],]
  exp=GDSC2_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  group=pancancer_group_fs
  library(dplyr)
  group = group %>%  
    mutate(group = ntile(activity_score,2))
  group$group = factor(group$group, levels =1:2, labels = c("L","H"))
  group$Tag=rownames(group)
  group=group[,c("Tag","group")]
  result=Get_limma(gene=exp,
                   group=group,
                   alternative_limma="greater")
  assign(paste0(vector[i],"_result"),result)
}
setwd("input")
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T,
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
grid.text(" ", x = 0.3, y = unit(1, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
# ------------------------------------------------------------------------------
# Calculation of WNT scores for CTRP
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input/")
library(GSEABase)
CTRP_exp=readRDS("input/DataFiles/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=as.data.frame(log2(CTRP_exp+1)),
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
CTRP_WNT_score=WNT_Score$final_activity_score
save(CTRP_WNT_score,file="CTRP_WNT_score.Rdata")
rm(list=ls())
gc()
library(readxl)
library(ggplot2)
setwd("input/DataFiles/DataFiles/GLDS/CTRPv2")
load("input/CTRP_WNT_score.Rdata")
group=read_excel("Harmonized_CCL_Data_v1.0.xlsx")
group=na.omit(group[,c(1,4)])
CTRP_WNT_score=merge(CTRP_WNT_score,group,by.x="row.names",by.y="cvcl cell line name")
colnames(CTRP_WNT_score)[5]="tissue"
CTRP_WNT_score=CTRP_WNT_score[CTRP_WNT_score$tissue!="NA",]
median_scores=aggregate(neg_score ~ tissue, data = CTRP_WNT_score, FUN = median)
sorted_sites=median_scores$tissue[order(median_scores$neg_score, decreasing = TRUE)]
colors=colorRampPalette(c("#f58f98", "#90d7ec"))(length(sorted_sites))
CTRP_WNT_score$tissue=factor(CTRP_WNT_score$tissue, levels = sorted_sites)
p = ggplot(CTRP_WNT_score, aes(x = tissue, y = neg_score, fill = tissue)) +  
  geom_boxplot(color = "black") +  
  theme_bw() +  
  theme(  
    axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16),  
    axis.text.y = element_text(color = "black", size = 16),  
    axis.title = element_text(color = "black", size = 16),  
    plot.title = element_text(color = "black", size = 16),  
    legend.position = 'none'  
  ) +  
  scale_fill_manual(values = colors) +  
  labs(  
    y = "Wnt/β-catenin pathway inhibition score", #activation/inhibition 
    x = ""  
  ) +   
  ggtitle("CTRP")  
p
# ------------------------------------------------------------------------------
# Validation of gene sets in CTRP
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
library(GSEABase)
library(dplyr)
setwd("input/DataFiles/DataFiles/GLDS/CTRPv2")
load("input/CTRP_WNT_score.Rdata")
group=read_excel("Harmonized_CCL_Data_v1.0.xlsx")
group=na.omit(group[,c(1,4)])
CTRP_WNT_score=merge(CTRP_WNT_score,group,by.x="row.names",by.y="cvcl cell line name")
colnames(CTRP_WNT_score)[5]="tissue"
CTRP_WNT_score=CTRP_WNT_score[CTRP_WNT_score$tissue!="NA",]
tissue_df=CTRP_WNT_score[, c("Row.names", "tissue")]  
agg_df=aggregate(. ~ Row.names, data = CTRP_WNT_score[, !(names(CTRP_WNT_score) %in% "tissue")], FUN = mean)  
tissue_unique=unique(tissue_df)
CTRP_WNT_score=merge(agg_df, tissue_unique, by = "Row.names", all.x = TRUE)
number=as.data.frame(table(CTRP_WNT_score$tissue))
number=number[number$Freq>=6,]
vector=unique(number$Var1)
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
rownames(CTRP_WNT_score)=CTRP_WNT_score[,1]
CTRP_WNT_score=CTRP_WNT_score[,-1]
CTRP_exp=readRDS("input/DataFiles/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP_exp=as.data.frame(log2(CTRP_exp+1))
for (i in 1:length(vector)){
  pancancer_group_fs=CTRP_WNT_score[CTRP_WNT_score$tissue==vector[i],]
  exp=CTRP_exp[,rownames(pancancer_group_fs)]
  exp=exp[rowSums(exp)>0,]
  group=pancancer_group_fs
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
setwd("input")
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
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T,
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60,
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1)
)
grid.text(" ", x = 0.3, y = unit(1, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))