# ==============================================================================
# WntQuant-Paper4 -- Deep pan-cancer explorations
# ------------------------------------------------------------------------------
# Purpose:    Pan-cancer subtyping by k-means clustering of Wnt activity scores and deep multi-omic characterization of the resulting subtypes.
# Inputs:     Pan-cancer Wnt activity scores and group annotations.
# Outputs:    Subtype assignments and downstream comparative analyses.
# NOTE:       Part of the code involves the core interests of the laboratory and is not open source at this time.
# ==============================================================================


rm(list=ls())
gc()
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
genesets=as.matrix(pancancer_WNT_Score$single_score_matrix)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
genesets=genesets[rownames(pancancer_group),]
library(mlbench)
library(factoextra)
library(ggplot2)
library (cluster)
library(fpc)
set.seed(123)
fviz_nbclust(genesets, kmeans, method = "wss") +   
  geom_vline(xintercept = 4, linetype = 2) +
  ggtitle("Optimal number of clusters") + 
  labs(x = "Number of clusters", y = "Total within sum of square") +
  theme_minimal() +
  theme(  
    text = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black")
  )
km_result = kmeans(genesets, centers=4,nstart = 1000,algorithm="Hartigan-Wong",iter.max=100)
panduan=as.data.frame(km_result$cluster)
panduan$Tag=rownames(panduan)
colnames(panduan)=c("cluster","Tag")
rank1=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==1,]$Tag,]$activity_score)
rank1
rank2=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==2,]$Tag,]$activity_score)
rank2
rank3=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==3,]$Tag,]$activity_score)
rank3
rank4=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==4,]$Tag,]$activity_score)
rank4
rank=data.frame(Tag=c(1,2,3,4),value=c(rank1, rank2, rank3, rank4))
rank=rank[order(-rank$value),]
km_result$cluster=ifelse(km_result$cluster==rank$Tag[1],1,ifelse(km_result$cluster==rank$Tag[2],2,ifelse(km_result$cluster==rank$Tag[3],3,4)))
dd=cbind(genesets, cluster = km_result$cluster)
fviz_cluster(km_result, data = genesets,  
             palette = c("#d71345","#f47920","#145b7d","#6950a1"),  
             ellipse.type = "euclid",  
             star.plot = TRUE,   
             repel = TRUE,  
             geom = "point",
             pointsize = 0.5) +   
  ggtitle("K-means") +   
  theme_minimal() +  
  theme(  
    text = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black")
  )  
cluster=data.frame(Tag = rownames(as.data.frame(dd)),cluster=as.data.frame(dd)$cluster,stringsAsFactors = FALSE)
cluster$cluster = paste("cluster",cluster$cluster, sep = "")
cluster$cluster=ifelse(cluster$cluster=="cluster1","panWSS1",
                       ifelse(cluster$cluster=="cluster2","panWSS2",
                              ifelse(cluster$cluster=="cluster3","panWSS3","panWSS4")))
WNT=merge(pancancer_WNT_Score$final_activity_score,cluster,by.x="row.names",by.y="Tag")
WNT$project="TCGA"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = activity_score, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c(paste0("panWSS1\n(n=",table(WNT$cluster)[1],")"),   
                              paste0("panWSS2\n(n=",table(WNT$cluster)[2],")"),  
                              paste0("panWSS3\n(n=",table(WNT$cluster)[3],")"),  
                              paste0("panWSS4\n(n=",table(WNT$cluster)[4],")"))) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16,  color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
  ggtitle("TCGA-Pancancer")  
print(p)
pancancer_cluster=cluster
save(pancancer_cluster,file="pancancer_cluster.Rdata")
# ------------------------------------------------------------------------------
# TARGET pan-cancer typing matrix
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/TARGET_WNT_Score.Rdata")
genesets=as.matrix(TARGET_WNT_Score$single_score_matrix)
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
genesets=genesets[rownames(TARGET_G),]
library(mlbench)
library(factoextra)
library(ggplot2)
library (cluster)
library(fpc)
set.seed(123)
fviz_nbclust(genesets, kmeans, method = "wss") +   
  geom_vline(xintercept = 4, linetype = 2) +
  ggtitle("Optimal number of clusters") +
  labs(x = "Number of clusters", y = "Total within sum of square") + 
  theme_minimal() +
  theme(  
    text = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black")
  )
km_result = kmeans(genesets, centers=4,nstart = 1000,algorithm="Hartigan-Wong",iter.max=100)
panduan=as.data.frame(km_result$cluster)
panduan$Tag=rownames(panduan)
colnames(panduan)=c("cluster","Tag")
rank1=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==1,]$Tag,]$activity_score)
rank1
rank2=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==2,]$Tag,]$activity_score)
rank2
rank3=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==3,]$Tag,]$activity_score)
rank3
rank4=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==4,]$Tag,]$activity_score)
rank4
rank=data.frame(Tag=c(1,2,3,4),value=c(rank1, rank2, rank3, rank4))
rank=rank[order(-rank$value),]
km_result$cluster=ifelse(km_result$cluster==rank$Tag[1],1,ifelse(km_result$cluster==rank$Tag[2],2,ifelse(km_result$cluster==rank$Tag[3],3,4)))
dd=cbind(genesets, cluster = km_result$cluster)
fviz_cluster(km_result, data = genesets,  
             palette = c("#d71345","#f47920","#145b7d","#6950a1"),  
             ellipse.type = "euclid",  
             star.plot = TRUE,   
             repel = TRUE,  
             geom = "point",
             pointsize = 0.5) +   
  ggtitle("K-means") +   
  theme_minimal() +  
  theme(  
    text = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black")
  )  
cluster=data.frame(Tag = rownames(as.data.frame(dd)),cluster=as.data.frame(dd)$cluster,stringsAsFactors = FALSE)
cluster$cluster = paste("cluster",cluster$cluster, sep = "")
cluster$cluster=ifelse(cluster$cluster=="cluster1","panWSS1",
                       ifelse(cluster$cluster=="cluster2","panWSS2",
                              ifelse(cluster$cluster=="cluster3","panWSS3","panWSS4")))
WNT=merge(TARGET_WNT_Score$final_activity_score,cluster,by.x="row.names",by.y="Tag")
WNT$project="TARGET"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = activity_score, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c(paste0("panWSS1\n(n=",table(WNT$cluster)[1],")"),   
                              paste0("panWSS2\n(n=",table(WNT$cluster)[2],")"),  
                              paste0("panWSS3\n(n=",table(WNT$cluster)[3],")"),  
                              paste0("panWSS4\n(n=",table(WNT$cluster)[4],")"))) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16,  color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
  ggtitle("TARGET-Pancancer")  
print(p)
TARGET_cluster=cluster
save(TARGET_cluster,file="TARGET_cluster.Rdata")
# ------------------------------------------------------------------------------
# TARGET+TCGA cumulative bar chart
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/TARGET_cluster.Rdata")
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
TARGET_G$Project=paste0("TARGET-",TARGET_G$Project)
load("input/pancancer_cluster.Rdata")
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
pancancer_group$Project=paste0("TCGA-",pancancer_group$Project)
pancancer_group=rbind(pancancer_group,TARGET_G)
pancancer_cluster=rbind(pancancer_cluster,TARGET_cluster)
WNT=merge(pancancer_cluster,pancancer_group,by.x="Tag",by.y="row.names")
df_counts = as.data.frame(table(WNT$Project, WNT$cluster))
colnames(df_counts) <- c("project", "cluster", "count")
df_percent = df_counts %>%
  group_by(project) %>%
  mutate(percent = count / sum(count) * 100)
ggplot(df_percent, aes(x = project, y = percent, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("#f7acbc", "#fedcbd", "#90d7ec", "#afb4db")) +
  labs(x = "", y = "Percentage (%)", fill = "Cluster") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
# ------------------------------------------------------------------------------
# Expression of the AXIN2 gene in a new typology of TCGA
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(reshape2)
setwd("input")
load("input/pancancer_exp.Rdata")
load("input/pancancer_cluster.Rdata")
group=pancancer_cluster
colnames(group)=c("Tag","group")
co_sample=intersect(colnames(pancancer_exp),group$Tag)
pancancer_exp=pancancer_exp[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(pancancer_exp)
exp=as.data.frame(t(exp[intersect(c("AXIN2","ZNRF3","RNF43","LGR5","TCF7"),rownames(exp)),]))
biomarker_group=merge(exp,group,by.x="row.names",by.y="Tag",all=T)
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("AXIN2","group")]
colnames(biomarker_group)=c("AXIN2","cluster")
WNT=biomarker_group
WNT$project="TCGA"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = AXIN2, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("AXIN2 expression") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c(paste0("panWSS1\n(n=",table(WNT$cluster)[1],")"),   
                              paste0("panWSS2\n(n=",table(WNT$cluster)[2],")"),  
                              paste0("panWSS3\n(n=",table(WNT$cluster)[3],")"),  
                              paste0("panWSS4\n(n=",table(WNT$cluster)[4],")"))) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16,  color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
  ggtitle("TCGA-Pancancer")  
print(p)
# ------------------------------------------------------------------------------
# Expression of the AXIN2 gene in a new typology of TARGET
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(reshape2)
setwd("input")
load("input/TARGET.Rdata")
load("input/TARGET_cluster.Rdata")
group=TARGET_cluster
colnames(group)=c("Tag","group")
co_sample=intersect(colnames(TARGET),group$Tag)
TARGET=TARGET[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(TARGET)
exp=as.data.frame(t(exp[intersect(c("AXIN2","ZNRF3","RNF43","LGR5","TCF7"),rownames(exp)),]))
biomarker_group=merge(exp,group,by.x="row.names",by.y="Tag",all=T)
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("AXIN2","group")]
colnames(biomarker_group)=c("AXIN2","cluster")
WNT=biomarker_group
WNT$project="TARGET"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = AXIN2, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("AXIN2 expression") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c(paste0("panWSS1\n(n=",table(WNT$cluster)[1],")"),   
                              paste0("panWSS2\n(n=",table(WNT$cluster)[2],")"),  
                              paste0("panWSS3\n(n=",table(WNT$cluster)[3],")"),  
                              paste0("panWSS4\n(n=",table(WNT$cluster)[4],")"))) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16,  color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
  ggtitle("TARGET-Pancancer")  
print(p)
# ------------------------------------------------------------------------------
# TCGA-OS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
survival_data=survival_data[,c("sample","OS","OS.time")]
survival_data_cms=merge(survival_data,pancancer_cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
survival_data_cms$OS.time=survival_data_cms$OS.time/30
survival_data_cms=na.omit(survival_data_cms)
library(survival)
library(survminer)
fitd = survdiff(Surv(OS.time, OS) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(OS.time, OS)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(OS.time, OS)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p=ggsurvplot(  
  fit = fit,  
  conf.int = FALSE,
  risk.table = F,
  risk.table.col = "strata",  
  palette = mycol,
  data = survival_data_cms,  
  xlim = c(0, 376),
  size = 0.1,  
  break.time.by = 48,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25, 
  title = paste0("OS(TCGA-Pancancer,n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) +
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# TARGET-OS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/TARGET_cluster.Rdata")
survival_data1=read.delim("TARGET-AML.survival.tsv")
survival_data2=read.delim("TARGET-ALL-P1.survival.tsv")
survival_data3=read.delim("TARGET-ALL-P2.survival.tsv")
survival_data4=read.delim("TARGET-ALL-P3.survival.tsv")
survival_data5=read.delim("TARGET-CCSK.survival.tsv")
survival_data6=read.delim("TARGET-NBL.survival.tsv")
survival_data7=read.delim("TARGET-RT.survival.tsv")
survival_data8=read.delim("TARGET-WT.survival.tsv")
survival_data9=read.delim("TARGET-OS.survival.tsv")
survival_data=rbind(survival_data1,survival_data2)
survival_data=rbind(survival_data,survival_data3)
survival_data=rbind(survival_data,survival_data4)
survival_data=rbind(survival_data,survival_data5)
survival_data=rbind(survival_data,survival_data6)
survival_data=rbind(survival_data,survival_data7)
survival_data=rbind(survival_data,survival_data8)
survival_data=rbind(survival_data,survival_data9)
survival_data=survival_data[,c(1,2,3)]
survival_data=na.omit(survival_data)
survival_data=survival_data[,c("sample","OS","OS.time")]
survival_data_cms=merge(survival_data,TARGET_cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
survival_data_cms$OS.time=survival_data_cms$OS.time/30
survival_data_cms=na.omit(survival_data_cms)
library(survival)
library(survminer)
fitd = survdiff(Surv(OS.time, OS) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(OS.time, OS)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(OS.time, OS)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p=ggsurvplot(  
  fit = fit,  
  conf.int = FALSE,
  risk.table = F,
  risk.table.col = "strata",  
  palette = mycol,
  data = survival_data_cms,  
  xlim = c(0, 153),
  size = 0.9,  
  break.time.by = 48,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = paste0("OS(TARGET-Pancancer,n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) +
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# DSS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
survival_data=survival_data[,c("sample","DSS","DSS.time")]
survival_data_cms=merge(survival_data,pancancer_cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
survival_data_cms$DSS.time=survival_data_cms$DSS.time/30
survival_data_cms=na.omit(survival_data_cms)
library(survival)
library(survminer)
fitd = survdiff(Surv(DSS.time, DSS) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(DSS.time, DSS)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(DSS.time, DSS)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p=ggsurvplot(  
  fit = fit,  
  conf.int = FALSE,
  risk.table = F,
  risk.table.col = "strata",  
  palette = mycol,
  data = survival_data_cms,  
  xlim = c(0, 376),
  size = 0.1,  
  break.time.by = 48,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = paste0("DSS(TCGA-Pancancer,n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) +
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# DFI
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
survival_data=survival_data[,c("sample","DFI","DFI.time")]
survival_data_cms=merge(survival_data,pancancer_cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
survival_data_cms$DFI.time=survival_data_cms$DFI.time/30
survival_data_cms=na.omit(survival_data_cms)
library(survival)
library(survminer)
fitd = survdiff(Surv(DFI.time, DFI) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(DFI.time, DFI)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(DFI.time, DFI)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p=ggsurvplot(  
  fit = fit,  
  conf.int = FALSE,
  risk.table = F,
  risk.table.col = "strata",  
  palette = mycol,
  data = survival_data_cms,  
  xlim = c(0, 376),
  size = 0.1,  
  break.time.by = 48,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = paste0("DFI(TCGA-Pancancer,n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) + 
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# PFI
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
survival_data=survival_data[,c("sample","PFI","PFI.time")]
survival_data_cms=merge(survival_data,pancancer_cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
survival_data_cms$PFI.time=survival_data_cms$PFI.time/30
survival_data_cms=na.omit(survival_data_cms)
library(survival)
library(survminer)
fitd = survdiff(Surv(PFI.time, PFI) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(PFI.time, PFI)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(PFI.time, PFI)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p=ggsurvplot(  
  fit = fit,  
  conf.int = FALSE,
  risk.table = F,
  risk.table.col = "strata",  
  palette = mycol,
  data = survival_data_cms,  
  xlim = c(0, 376),
  size = 0.1,  
  break.time.by = 48,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = paste0("PFI(TCGA-Pancancer,n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) +
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# TCGA-age
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
age=survival_data[,c(1,4)]
colnames(age)=c("sample","age")
age=na.omit(age)
age$age=ifelse(age$age<=40,"age≤40",
               ifelse(age$age<=60,"40<age≤60","age>60"))
age=merge(pancancer_cluster,age,by.x="Tag",by.y="sample")
age$age=factor(age$age,levels=c("age≤40","40<age≤60","age>60"))
age$status="Age"
WNT_cluster1=age[age$cluster=="panWSS1",]
WNT_cluster2=age[age$cluster=="panWSS2",]
WNT_cluster3=age[age$cluster=="panWSS3",]
WNT_cluster4=age[age$cluster=="panWSS4",]
age=age[,-c(1,4)]
age=as.data.frame(table(age))
contingency_table=acast(age, cluster ~ age, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤40" = "#90d7ec",   
                               "40<age≤60" = "#009ad6",   
                               "age>60" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=10364, Age)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = "none", 
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤40" = "#90d7ec",   
                               "40<age≤60" = "#009ad6",   
                               "age>60" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤40" = "#90d7ec",   
                               "40<age≤60" = "#009ad6",   
                               "age>60" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤40" = "#90d7ec",   
                               "40<age≤60" = "#009ad6",   
                               "age>60" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-gender
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
gender=survival_data[,c(1,5)]
gender$gender=ifelse(gender$gender=="MALE","Male","Female")
gender=merge(pancancer_cluster,gender,by.x="Tag",by.y="sample")
gender=gender[,c("cluster","gender")]
gender=na.omit(gender)
gender$status="gender"
WNT_cluster1=gender[gender$cluster=="panWSS1",]
WNT_cluster2=gender[gender$cluster=="panWSS2",]
WNT_cluster3=gender[gender$cluster=="panWSS3",]
WNT_cluster4=gender[gender$cluster=="panWSS4",]
gender=gender[,-3]
gender=as.data.frame(table(gender))
contingency_table=acast(gender, cluster ~ gender, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=10415, Gender)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none", 
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16) 
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot)  
# ------------------------------------------------------------------------------
# TCGA-race
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
race=survival_data[,c(1,6)]
table(race$race)
race=race[!race$race%in%c("","[Not Evaluated]","[Unknown]"),]
race$race=ifelse(race$race=="AMERICAN INDIAN OR ALASKA NATIVE","Indian or Alaskan",
                 ifelse(race$race=="ASIAN","Asian",
                        ifelse(race$race=="BLACK OR AFRICAN AMERICAN","Black",
                               ifelse(race$race=="NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER","Hawaiian",
                                      "White"))))
race=race[!race$race%in%c("Hawaiian","Indian or Alaskan"),]
race=merge(pancancer_cluster,race,by.x="Tag",by.y="sample")
race=race[,c("cluster","race")]
race=na.omit(race)
race$status="race"
WNT_cluster1=race[race$cluster=="panWSS1",]
WNT_cluster2=race[race$cluster=="panWSS2",]
WNT_cluster3=race[race$cluster=="panWSS3",]
WNT_cluster4=race[race$cluster=="panWSS4",]
race=race[,-3]
race=as.data.frame(table(race))
contingency_table=acast(race, cluster ~ race, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=9113, Race)") +
  theme(  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16) 
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-stage
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
stage=survival_data[,c(1,8)]
colnames(stage)=c("sample","stage")
table(stage$stage)
stage=stage[!stage$stage%in%c("","[Discrepancy]","Stage IS"),]
stage$stage=ifelse(stage$stage%in%c("I","Stage I","Stage IA","Stage IA1","Stage IA2","Stage IB","Stage IB1","Stage IB2","Stage IC"),"Stage1",
                   ifelse(stage$stage%in%c("IIa","IIb","Stage II","Stage IIA","Stage IIA1","Stage IIA2","Stage IIB","Stage IIC"),"Stage2",
                          ifelse(stage$stage%in%c("III","Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IIIC1","Stage IIIC2"),"Stage3","Stage4")))
stage=na.omit(stage)
stage=merge(pancancer_cluster,stage,by.x="Tag",by.y="sample")
stage$status="stage"
WNT_cluster1=stage[stage$cluster=="panWSS1",]
WNT_cluster2=stage[stage$cluster=="panWSS2",]
WNT_cluster3=stage[stage$cluster=="panWSS3",]
WNT_cluster4=stage[stage$cluster=="panWSS4",]
stage=stage[,-c(1,4)]
stage=as.data.frame(table(stage))
contingency_table=acast(stage, cluster ~ stage, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = stage)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Stage0" = "#cde6c7", "Stage1" = "#abc88b",   
                               "Stage2" = "#84bf96", "Stage3" = "#45b97c",   
                               "Stage4" = "#005831")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=2423, Clinical stage)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = stage)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Stage0" = "#cde6c7", "Stage1" = "#abc88b",   
                               "Stage2" = "#84bf96", "Stage3" = "#45b97c",   
                               "Stage4" = "#005831")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = stage)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Stage0" = "#cde6c7", "Stage1" = "#abc88b",   
                               "Stage2" = "#84bf96", "Stage3" = "#45b97c",   
                               "Stage4" = "#005831")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = stage)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Stage0" = "#cde6c7", "Stage1" = "#abc88b",   
                               "Stage2" = "#84bf96", "Stage3" = "#45b97c",   
                               "Stage4" = "#005831")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-grade
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
grade=survival_data[,c(1,10)]
colnames(grade)=c("sample","grade")
table(grade$grade)
grade=grade[!grade$grade%in%c("","[Discrepancy]","[Unknown]","GB","GX"),]
table(grade$grade)
grade$grade=ifelse(grade$grade%in%c("G1","G2","Low Grade"),"Low","High")
grade=na.omit(grade)
grade=merge(pancancer_cluster,grade,by.x="Tag",by.y="sample")
grade$status="grade"
WNT_cluster1=grade[grade$cluster=="panWSS1",]
WNT_cluster2=grade[grade$cluster=="panWSS2",]
WNT_cluster3=grade[grade$cluster=="panWSS3",]
WNT_cluster4=grade[grade$cluster=="panWSS4",]
grade=grade[,-c(1,4)]
grade=as.data.frame(table(grade))
contingency_table=acast(grade, cluster ~ grade, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Low" = "#afb4db", "High" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=4322, Grade)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none", 
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Low" = "#afb4db", "High" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Low" = "#afb4db", "High" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Low" = "#afb4db", "High" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# TCGA-status
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
status=survival_data[,c(1,15)]
colnames(status)=c("sample","tumor_status")
table(status$tumor_status)
status=status[!status$tumor_status%in%c("","[Discrepancy]"),]
status$tumor_status=ifelse(status$tumor_status=="TUMOR FREE","Tumor free","With tumor")
status=na.omit(status)
status=merge(pancancer_cluster,status,by.x="Tag",by.y="sample")
status$status="status"
WNT_cluster1=status[status$cluster=="panWSS1",]
WNT_cluster2=status[status$cluster=="panWSS2",]
WNT_cluster3=status[status$cluster=="panWSS3",]
WNT_cluster4=status[status$cluster=="panWSS4",]
status=status[,-c(1,4)]
status=as.data.frame(table(status))
contingency_table=acast(status, cluster ~ tumor_status, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = tumor_status)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Tumor free" = "#f7acbc", "With tumor" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=9601, Status)") + 
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = tumor_status)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Tumor free" = "#f7acbc", "With tumor" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = tumor_status)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Tumor free" = "#f7acbc", "With tumor" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = tumor_status)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Tumor free" = "#f7acbc", "With tumor" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-outcome
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("clinical.tsv")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
outcome=survival_data[,c(1,23)]
colnames(outcome)=c("sample","outcome")
table(outcome$outcome)
outcome=outcome[!outcome$outcome%in%c("","[Discrepancy]","[Unknown]","[Not Evaluated]"),]
outcome$outcome=ifelse(outcome$outcome%in%c("Complete Remission/Response","No Measureable Tumor or Tumor Markers","Normalization of Tumor Markers, but Residual Tumor Mass","Partial Remission/Response"),"CR/PR","SD/PD")
outcome=na.omit(outcome)
outcome=merge(pancancer_cluster,outcome,by.x="Tag",by.y="sample")
outcome$status="status"
WNT_cluster1=outcome[outcome$cluster=="panWSS1",]
WNT_cluster2=outcome[outcome$cluster=="panWSS2",]
WNT_cluster3=outcome[outcome$cluster=="panWSS3",]
WNT_cluster4=outcome[outcome$cluster=="panWSS4",]
outcome=outcome[,-c(1,4)]
outcome=as.data.frame(table(outcome))
contingency_table=acast(outcome, cluster ~ outcome, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = outcome)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CR/PR" = "#009ad6", "SD/PD" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=5250, Outcome)") + 
  theme(  
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic") 
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = outcome)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CR/PR" = "#009ad6", "SD/PD" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = "none", 
    text = element_text(color = "black", size = 16) 
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = outcome)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CR/PR" = "#009ad6", "SD/PD" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = outcome)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CR/PR" = "#009ad6", "SD/PD" = "#d71345")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-margin
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
margin=survival_data[,c(1,24)]
colnames(margin)=c("sample","margin")
table(margin$margin)
margin=margin[!margin$margin%in%c("","[Unknown]"),]
margin=na.omit(margin)
margin=merge(pancancer_cluster,margin,by.x="Tag",by.y="sample")
margin$status="status"
WNT_cluster1=margin[margin$cluster=="panWSS1",]
WNT_cluster2=margin[margin$cluster=="panWSS2",]
WNT_cluster3=margin[margin$cluster=="panWSS3",]
WNT_cluster4=margin[margin$cluster=="panWSS4",]
margin=margin[,-c(1,4)]
margin=as.data.frame(table(margin))
contingency_table=acast(margin, cluster ~ margin, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = margin)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Negative" = "#94d6da", "Positive" = "#f15b6c","Close"="#d3d7d4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=1245, Margin status)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = margin)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Negative" = "#94d6da", "Positive" = "#f15b6c","Close"="#d3d7d4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = margin)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Negative" = "#94d6da", "Positive" = "#f15b6c","Close"="#d3d7d4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = margin)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Negative" = "#94d6da", "Positive" = "#f15b6c","Close"="#d3d7d4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-residual
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
residual=survival_data[,c(1,25)]
colnames(residual)=c("sample","residual")
table(residual$residual)
residual=residual[!residual$residual%in%c("","[Unknown]"),]
residual=na.omit(residual)
residual=merge(pancancer_cluster,residual,by.x="Tag",by.y="sample")
residual$status="status"
WNT_cluster1=residual[residual$cluster=="panWSS1",]
WNT_cluster2=residual[residual$cluster=="panWSS2",]
WNT_cluster3=residual[residual$cluster=="panWSS3",]
WNT_cluster4=residual[residual$cluster=="panWSS4",]
residual=residual[,-c(1,4)]
residual=as.data.frame(table(residual))
contingency_table=acast(residual, cluster ~ residual, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = residual)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("R0" = "#f173ac", "R1" = "#c77eb5","R2"="#8552a1","RX"="#f2eada")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=1190, Residual tumor)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = residual)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("R0" = "#f173ac", "R1" = "#c77eb5","R2"="#8552a1","RX"="#f2eada")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = residual)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("R0" = "#f173ac", "R1" = "#c77eb5","R2"="#8552a1","RX"="#f2eada")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = residual)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("R0" = "#f173ac", "R1" = "#c77eb5","R2"="#8552a1","RX"="#f2eada")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA-mutation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
mut_data=read.delim("mc3.v0.2.8.PUBLIC.nonsilentGene.xena",check.names = F)
mut_data=mut_data[mut_data$sample%in%c("APC","AXIN1","AXIN2","RNF43","DKK1","DKK3","TCF7L2","GSK3B","CTNNB1"),]
rownames(mut_data)=mut_data[,1]
mut_data=mut_data[,-1]
mut_data=t(mut_data)
mut_data=ifelse(mut_data==1,"MT","WT")
mut_data=as.data.frame(mut_data)
save(mut_data,file="mut_data.Rdata")
load("input/mut_data.Rdata")
gene=merge(pancancer_cluster,mut_data,by.x="Tag",by.y="row.names")
gene=gene[,c("cluster","APC")]
colnames(gene)=c("cluster","gene")
gene=na.omit(gene)
gene$status="gene"
WNT_cluster1=gene[gene$cluster=="panWSS1",]
WNT_cluster2=gene[gene$cluster=="panWSS2",]
WNT_cluster3=gene[gene$cluster=="panWSS3",]
WNT_cluster4=gene[gene$cluster=="panWSS4",]
gene=gene[,-3]
gene=as.data.frame(table(gene))
contingency_table=acast(gene, cluster ~ gene, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = gene)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TCGA(n=8797, APC)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = gene)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = gene)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = gene)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# Protein labelling standardisation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
protein_data=read.delim("TCGA-RPPA-pancan-clean.xena",check.names = F)
rownames(protein_data)=protein_data[,1]
protein_data=protein_data[,-1]
protein_data=as.data.frame(t(protein_data))
protein_data=merge(protein_data,pancancer_cluster,by.x="row.names",by.y="Tag")
rownames(protein_data)=protein_data[,1]
protein_data=protein_data[,-1]
protein_data=protein_data[, colSums(is.na(protein_data)) < nrow(protein_data)]
protein_data=melt(protein_data)
colnames(protein_data)=c("Cluster","protein","value")
protein_data=na.omit(protein_data)
proteins=unique(protein_data$protein)
print(unique(protein_data$protein))
protein_names=c(
  "14-3-3ε", "4E-BP1", "4E-BP1_pS65", "4E-BP1_pT37/T46", "53BP1", "ACC_pS79", "ACC1", "AKT", "AKT_pS473", "AKT_pT308",
  "AMPKα", "AMPKα_pT172", "AR", "ASNS", "ATM", "BAK", "BAX", "BCL-2", "BCL-xL", "Beclin-1",
  "β-catenin", "BID", "BIM", "c-Jun_pS73", "c-KIT", "c-MET_pY1235", "c-MYC", "CRAF", "CRAF_pS338", "Caspase-7 (cleaved D198)",
  "Caveolin-1", "CD31 (PECAM-1)", "CD49b", "CDK1", "CHK1", "CHK1_pS345", "CHK2", "CHK2_pT68", "cIAP", "Claudin-7",
  "Collagen VI", "Cyclin B1", "Cyclin D1", "Cyclin E1", "DJ-1", "Dvl3", "E-cadherin", "eEF2", "eEF2K", "EGFR",
  "EGFR_pY1068", "EGFR_pY1173", "eIF4E", "ERα", "ERα_pS118", "ERK2", "Fibronectin", "FOXO3A", "GAB2", "GATA3",
  "GSK-3α/β", "GSK-3α/β_pS21/S9", "HER2", "HER2_pY1248", "HER3", "HER3_pY1289", "HSP70", "IGFBP2", "INPP4B", "IRS-1",
  "JNK_pT183/Y185", "JNK2", "Ku80", "LCK", "LKB1", "MAPK_pT202/Y204", "MEK1", "MEK1_pS217/S221", "MIG6", "MRE11",
  "mTOR", "mTOR_pS2448", "N-cadherin", "NF-κB p65_pS536", "NF2 (Merlin)", "Notch1", "P-cadherin", "p27", "p27_pT157", "p38 MAPK",
  "p38_pT180/Y182", "p53", "p70 S6K", "p70 S6K_pT389", "p90 RSK_pT359/S363", "PAI-1", "Paxillin", "PCNA", "PDK1_pS241", "PEA15",
  "PI3K p110α", "PKCα", "PKCα_pS657", "PKCδ_pS664", "PR", "PRAS40_pT246", "PTEN", "RAD50", "RAD51", "Rb_pS807/S811",
  "Ribosomal protein S6", "S6_pS235/S236", "S6_pS240/S244", "Shc_pY317", "SMAD1", "SMAD3", "SMAD4", "Src", "Src_pY416", "Src_pY527",
  "STAT3_pY705", "STAT5α", "Stathmin", "SYK", "Tuberin (TSC2)", "VEGFR2", "XRCC1", "YAP", "YAP_pS127", "YB-1",
  "YB-1_pS102", "4E-BP1_pT70", "ARAF_pS299", "Annexin VII", "ARID1A", "BRAF", "BAD_pS112", "BAP1", "BRCA2", "CD20",
  "Cyclin E2", "ETS-1", "eIF4G", "FASN", "FOXO3A_pS318/S321", "FOXM1", "G6PD", "GAPDH", "GSK-3_pS9", "Heregulin",
  "Myosin-11", "Myosin IIA_pS1943", "NRAS", "NDRG1_pT346", "p21", "p27_pT198", "p90 RSK", "PDCD4", "PDK1", "PEA15_pS116",
  "PI3K p85", "PKCβII_pS660", "PRDX1", "Rab11", "Rab25", "Raptor", "RBM15", "Rictor", "Rictor_pT1135", "SCD1",
  "SF2 (SRSF1)", "TAZ", "TIGAR", "Transglutaminase", "Transferrin receptor (TfR)", "TSC1", "Tuberin_pT1462", "EPPK1", "XBP1", "Acetyl-α-Tubulin (K40)",
  "CD62L ligand", "14-3-3β", "14-3-3ζ", "ACVRL1", "DIRAS3", "Annexin A1", "PREX1", "ENY2", "GCN5L2 (KAT2B)", "ADAR1",
  "JAB1 (COPS5)", "c-MET", "Caspase-8", "ERCC1", "MSH2", "MSH6", "PARP (cleaved)", "Rb", "SETD2", "SMAC (DIABLO)",
  "Snail", "AXL", "Myosin IIA (MYH9)", "SLC1A5 (ASCT2)", "GATA6", "BRD4", "CDK1_pY15", "ARAF", "BRAF_pS445", "BCL2A1",
  "c-ABL", "Caspase-3", "CD26 (DPP4)", "CHK1_pS296", "COG3", "DUSP4", "ERCC5 (XPG)", "IGF1R_pY1135/Y1136", "IRF1", "JAK2",
  "p16 INK4A", "SHP2_pY542", "PD-L1 (CD274)", "PARP1", "CA9 (Carbonic anhydrase IX)", "Complex II subunit 30 (SDHB)", "Glycogenin-1", "Glycogen synthase (GYS)", "GYS_pS641", "HIF-1α",
  "LDHA", "LDHB", "Mitochondrial marker", "ATP synthase β (OXPHOS)", "PKM2", "PYGB", "PYGB (ab2)", "PYGL", "PYGM", "CTLA-4",
  "PD-1 (PDCD1)", "Caspase-9", "E2F1", "EZH2", "KEAP1", "Lipocalin-2", "MACC1", "NRF2", "PARP3", "Thymidylate synthase",
  "TTF-1", "Chromogranin A (N-term)", "CK5", "Napsin A", "p63", "RET_pY905", "Synaptophysin", "α-catenin"
)
protein_tag=data.frame(Tag=unique(protein_data$protein),Names=protein_names)
save(protein_tag,file="protein_tag.Rdata")
# ------------------------------------------------------------------------------
# Expression of β-catenin protein
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
protein_data=read.delim("TCGA-RPPA-pancan-clean.xena",check.names = F)
rownames(protein_data)=protein_data[,1]
protein_data=protein_data[,-1]
protein_data=as.data.frame(t(protein_data))
protein_data=merge(protein_data,pancancer_cluster,by.x="row.names",by.y="Tag")
WNT=protein_data[,c("BETACATENIN","cluster")]
WNT$project="TCGA"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = BETACATENIN, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("β-catenin protein level") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c(paste0("panWSS1\n(n=",table(WNT$cluster)[1],")"),   
                              paste0("panWSS2\n(n=",table(WNT$cluster)[2],")"),  
                              paste0("panWSS3\n(n=",table(WNT$cluster)[3],")"),  
                              paste0("panWSS4\n(n=",table(WNT$cluster)[4],")"))) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16,  color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
  ggtitle("TCGA-Pancancer")  
print(p)
# ------------------------------------------------------------------------------
# Protein differential expression
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
library(limma)
library(ggrepel)
limma.one.sided = function(fit, lower = FALSE){
  se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total=fit$df.prior + fit$df.residual
  pt(fit$t, df = df.total, lower.tail = lower)
}
load("input/pancancer_cluster.Rdata")
protein_data=read.delim("TCGA-RPPA-pancan-clean.xena",check.names = F)
rownames(protein_data)=protein_data[,1]
protein_data=protein_data[,-1]
# ------------------------------------------------------------------------------
# panWSS1 vs. Other
# ------------------------------------------------------------------------------
pancancer_cluster_G=pancancer_cluster
pancancer_cluster_G$cluster=ifelse(pancancer_cluster_G$cluster=="panWSS1","panWSS1","Other")
sample=intersect(pancancer_cluster_G$Tag,colnames(protein_data))
pancancer_cluster_G=pancancer_cluster_G[pancancer_cluster_G$Tag%in%sample,]
protein_data_exp=protein_data[,pancancer_cluster_G$Tag]
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(panWSS1-Other, levels=design)
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
colnames(pancancer_cluster_G)=c("Tag","group")
result=Get_limma(gene=protein_data_exp,
                 group=pancancer_cluster_G,
                 alternative_limma="two.sided")
load("protein_tag.Rdata")
result = merge(result, protein_tag, by.x = "gene_id", by.y = "Tag")
result$Sig = "No diff"
result$Sig[result$logFC >= 0.1 & result$P.Value < 1e-10] = "Up"
result$Sig[result$logFC <= -0.1 & result$P.Value < 1e-10] = "Down"
result_up = result[result$logFC > 0, ]
result_dn = result[result$logFC < 0, ]
result_up = result_up[order(result_up$P.Value), ]
result_dn = result_dn[order(result_dn$P.Value), ]
top_genes = rbind(head(result_up, 8), head(result_dn, 8))
volcano_p = ggplot(result, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = Sig, size = -log10(P.Value)), alpha = 0.5) +
  scale_color_manual(values = c("Down" = "#90d7ec", 
                                "No diff" = "grey30", 
                                "Up" = "#f7acbc")) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  labs(x = expression(log[2]("FC")), 
       y = expression(-log[10]("p"))) +
  xlim(-1.7, 1.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = Names),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.5,
    force = 1,
    nudge_x = 0.2,
    nudge_y = 0.2
  )+ggtitle("panWSS1")
print(volcano_p)
# ------------------------------------------------------------------------------
# panWSS2 vs. Other
# ------------------------------------------------------------------------------
pancancer_cluster_G=pancancer_cluster
pancancer_cluster_G$cluster=ifelse(pancancer_cluster_G$cluster=="panWSS2","panWSS2","Other")
sample=intersect(pancancer_cluster_G$Tag,colnames(protein_data))
pancancer_cluster_G=pancancer_cluster_G[pancancer_cluster_G$Tag%in%sample,]
protein_data_exp=protein_data[,pancancer_cluster_G$Tag]
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(panWSS2-Other, levels=design)
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
colnames(pancancer_cluster_G)=c("Tag","group")
result=Get_limma(gene=protein_data_exp,
                 group=pancancer_cluster_G,
                 alternative_limma="two.sided")
result$P.Value=ifelse(result$P.Value==0,2.094838e-321*0.1,result$P.Value)
load("protein_tag.Rdata")
result = merge(result, protein_tag, by.x = "gene_id", by.y = "Tag")
result$Sig = "No diff"
result$Sig[result$logFC >= 0.1 & result$P.Value < 1e-10] = "Up"
result$Sig[result$logFC <= -0.1 & result$P.Value < 1e-10] = "Down"
result_up = result[result$logFC > 0, ]
result_dn = result[result$logFC < 0, ]
result_up = result_up[order(result_up$P.Value), ]
result_dn = result_dn[order(result_dn$P.Value), ]
top_genes = rbind(head(result_up, 8), head(result_dn, 8))
volcano_p = ggplot(result, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = Sig, size = -log10(P.Value)), alpha = 0.5) +
  scale_color_manual(values = c("Down" = "#90d7ec", 
                                "No diff" = "grey30", 
                                "Up" = "#f7acbc")) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"), # Keep the axis lines
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  labs(x = expression(log[2]("FC")), 
       y = expression(-log[10]("p"))) +
  xlim(-2, 2) +
  geom_text_repel(
    data = top_genes,
    aes(label = Names),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.5,
    force = 1,
    nudge_x = 0.2,
    nudge_y = 0.2
  )+ggtitle("panWSS2")
print(volcano_p)
# ------------------------------------------------------------------------------
# panWSS3 vs. Other
# ------------------------------------------------------------------------------
pancancer_cluster_G=pancancer_cluster
pancancer_cluster_G$cluster=ifelse(pancancer_cluster_G$cluster=="panWSS3","panWSS3","Other")
sample=intersect(pancancer_cluster_G$Tag,colnames(protein_data))
pancancer_cluster_G=pancancer_cluster_G[pancancer_cluster_G$Tag%in%sample,]
protein_data_exp=protein_data[,pancancer_cluster_G$Tag]
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(panWSS3-Other, levels=design)
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
colnames(pancancer_cluster_G)=c("Tag","group")
result=Get_limma(gene=protein_data_exp,
                 group=pancancer_cluster_G,
                 alternative_limma="two.sided")
load("protein_tag.Rdata")
result = merge(result, protein_tag, by.x = "gene_id", by.y = "Tag")
result$Sig = "No diff"
result$Sig[result$logFC >= 0.1 & result$P.Value < 1e-10] = "Up"
result$Sig[result$logFC <= -0.1 & result$P.Value < 1e-10] = "Down"
result_up = result[result$logFC > 0, ]
result_dn = result[result$logFC < 0, ]
result_up = result_up[order(result_up$P.Value), ]
result_dn = result_dn[order(result_dn$P.Value), ]
top_genes = rbind(head(result_up, 8), head(result_dn, 8))
volcano_p = ggplot(result, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = Sig, size = -log10(P.Value)), alpha = 0.5) +
  scale_color_manual(values = c("Down" = "#90d7ec", 
                                "No diff" = "grey30", 
                                "Up" = "#f7acbc")) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  labs(x = expression(log[2]("FC")), 
       y = expression(-log[10]("p"))) +
  xlim(-1.7, 1.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = Names),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.5,
    force = 1,
    nudge_x = 0.2,
    nudge_y = 0.2 
  )+ggtitle("panWSS3")
print(volcano_p)
# ------------------------------------------------------------------------------
# panWSS4 vs. Other
# ------------------------------------------------------------------------------
pancancer_cluster_G=pancancer_cluster
pancancer_cluster_G$cluster=ifelse(pancancer_cluster_G$cluster=="panWSS4","panWSS4","Other")
sample=intersect(pancancer_cluster_G$Tag,colnames(protein_data))
pancancer_cluster_G=pancancer_cluster_G[pancancer_cluster_G$Tag%in%sample,]
protein_data_exp=protein_data[,pancancer_cluster_G$Tag]
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(panWSS4-Other, levels=design)
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
colnames(pancancer_cluster_G)=c("Tag","group")
result=Get_limma(gene=protein_data_exp,
                 group=pancancer_cluster_G,
                 alternative_limma="two.sided")
load("protein_tag.Rdata")
result = merge(result, protein_tag, by.x = "gene_id", by.y = "Tag")
result$Sig = "No diff"
result$Sig[result$logFC >= 0.1 & result$P.Value < 1e-10] = "Up"
result$Sig[result$logFC <= -0.1 & result$P.Value < 1e-10] = "Down"
result_up = result[result$logFC > 0, ]
result_dn = result[result$logFC < 0, ]
result_up = result_up[order(result_up$P.Value), ]
result_dn = result_dn[order(result_dn$P.Value), ]
top_genes = rbind(head(result_up, 8), head(result_dn, 8))
volcano_p = ggplot(result, aes(logFC, -log10(P.Value))) +
  geom_point(aes(color = Sig, size = -log10(P.Value)), alpha = 0.5) +
  scale_color_manual(values = c("Down" = "#90d7ec", 
                                "No diff" = "grey30", 
                                "Up" = "#f7acbc")) +
  theme_bw() +
  theme(
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 15, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15)
  ) +
  labs(x = expression(log[2]("FC")), 
       y = expression(-log[10]("p"))) +
  xlim(-1.7, 1.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = Names),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.2,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.size = 0.5,
    force = 1,
    nudge_x = 0.2,
    nudge_y = 0.2
  )+ggtitle("panWSS4")
print(volcano_p)
# ------------------------------------------------------------------------------
# Identification of subtypes at the mRNA level of biological identity
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
pancancer_exp=pancancer_exp[,rownames(pancancer_group)]
hallmarker_activity=calculate_sig_score(eset = pancancer_exp,signature = hallmarker_gmt,method="ssgsea",mini_gene_count=1)
save(hallmarker_activity,file="hallmarker_activity.Rdata")
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
library(limma)
library(GSEABase)
load("input/pancancer_cluster.Rdata")
load("input/pancancer_exp.Rdata")
pancancer_exp=pancancer_exp[,pancancer_cluster$Tag]
load("input/hallmarker_activity.Rdata")
hallmarker_activity=as.data.frame(hallmarker_activity)
rownames(hallmarker_activity)=hallmarker_activity[,1]
hallmarker_activity=hallmarker_activity[,-1]
hallmarker_activity=as.data.frame(scale(hallmarker_activity))
tme_all=merge(hallmarker_activity,pancancer_cluster,by.x="row.names",by.y="Tag")
tme_all=melt(tme_all)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
tme_all$cluster <- factor(tme_all$cluster, levels = c("panWSS4", "panWSS3", "panWSS2", "panWSS1"))
p = ggplot(tme_all, aes(x = variable, y = value, fill = cluster)) +
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    linewidth = 0.5,
    alpha = 0.8,
    position = position_dodge(width = 0.8)
  ) +
  labs(x = "", y = "ssGSEA score") +
  theme_bw() +
  theme(
    legend.position = "left",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
  ) +
  scale_fill_manual(
    values = c("#d71345","#f47920","#145b7d","#6950a1"),
    breaks = c("panWSS1", "panWSS2", "panWSS3", "panWSS4"),
    labels = c("panWSS1", "panWSS2", "panWSS3", "panWSS4")
  ) +
  stat_compare_means(aes(group = cluster, label = ..p.signif..)) +
  ggtitle("TCGA") +
  coord_flip()
p
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
TARGET=TARGET[,rownames(TARGET_G)]
TARGET_hallmarker_activity=calculate_sig_score(eset = TARGET,signature = hallmarker_gmt,method="ssgsea",mini_gene_count=1)
save(TARGET_hallmarker_activity,file="TARGET_hallmarker_activity.Rdata")
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
library(limma)
library(GSEABase)
load("input/TARGET_cluster.Rdata")
load("input/TARGET.Rdata")
TARGET=TARGET[,TARGET_cluster$Tag]
load("input/TARGET_hallmarker_activity.Rdata")
TARGET_hallmarker_activity=as.data.frame(TARGET_hallmarker_activity)
rownames(TARGET_hallmarker_activity)=TARGET_hallmarker_activity[,1]
TARGET_hallmarker_activity=TARGET_hallmarker_activity[,-1]
TARGET_hallmarker_activity=as.data.frame(scale(TARGET_hallmarker_activity))
tme_all=merge(TARGET_hallmarker_activity,TARGET_cluster,by.x="row.names",by.y="Tag")
tme_all=melt(tme_all)
tme_all=tme_all[tme_all$variable!="WNT_β_catenin_signaling",]
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
tme_all$cluster=factor(tme_all$cluster,levels = c("panWSS1","panWSS2","panWSS3","panWSS4"))
p = ggplot(tme_all, aes(variable, value, fill = cluster, color = cluster)) +
  geom_boxplot(outlier.shape = NA,color = "black",linewidth = 0.5,alpha = 0.8) +
  labs(x = " ", y = "ssGSEA score") +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(linewidth = 1),
    panel.grid.major.y = element_line(linewidth = 1),
    panel.grid.minor.x = element_line(linewidth = 0.8),
    panel.grid.minor.y = element_line(linewidth = 0.8),
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),
    axis.line.y = element_line(colour = 'black', linewidth = 1),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.y = element_text(angle = 90, size = 12, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 12, colour = "black", hjust = 0)
  ) +
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +
  stat_compare_means(aes(group = cluster, label = ..p.signif..)) +
  ggtitle("TARGET")
p
# ------------------------------------------------------------------------------
# HRD
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/pancancer_cluster.Rdata")
HRD_data=read.delim("TCGA.HRD_withSampleID.txt",check.names = F)
rownames(HRD_data)=HRD_data[,1]
HRD_data=HRD_data[,-1]
HRD_data=as.data.frame(t(HRD_data))
colnames(HRD_data)=c("TAI","LST","HRD-LOH","HRD")
WNT=merge(HRD_data,pancancer_cluster,by.x="row.names",by.y="Tag")
WNT$project="TCGA"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = HRD, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Homologous recombination deficiency") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c(paste0("panWSS1\n(n=",table(WNT$cluster)[1],")"),   
                              paste0("panWSS2\n(n=",table(WNT$cluster)[2],")"),  
                              paste0("panWSS3\n(n=",table(WNT$cluster)[3],")"),  
                              paste0("panWSS4\n(n=",table(WNT$cluster)[4],")"))) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16,  color = "black"), 
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black") 
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
  ggtitle("TCGA-Pancancer")  
print(p)
# ------------------------------------------------------------------------------
# TARGET-age
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/TARGET_cluster.Rdata")
rt1_G=read.table(file = "TARGET-AML.clinical.tsv",header = T,sep = "\t",check.names = F)
rt1_G=rt1_G[,c("sample","age_at_index.demographic")]
rt2_G=read.table(file = "TARGET-ALL-P1.clinical.tsv",header = T,sep = "\t",check.names = F)
rt2_G=rt2_G[,c("sample","age_at_index.demographic")]
rt3_G=read.table(file = "TARGET-ALL-P2.clinical.tsv",header = T,sep = "\t",check.names = F)
rt3_G=rt3_G[,c("sample","age_at_index.demographic")]
rt4_G=read.table(file = "TARGET-ALL-P3.clinical.tsv",header = T,sep = "\t",check.names = F)
rt4_G=rt4_G[,c("sample","age_at_index.demographic")]
rt5_G=read.table(file = "TARGET-CCSK.clinical.tsv",header = T,sep = "\t",check.names = F)
rt5_G=rt5_G[,c("sample","age_at_index.demographic")]
rt6_G=read.table(file = "TARGET-NBL.clinical.tsv",header = T,sep = "\t",check.names = F)
rt6_G=rt6_G[,c("sample","age_at_index.demographic")]
rt7_G=read.table(file = "TARGET-OS.clinical.tsv",header = T,sep = "\t",check.names = F)
rt7_G=rt7_G[,c("sample","age_at_index.demographic")]
rt8_G=read.table(file = "TARGET-RT.clinical.tsv",header = T,sep = "\t",check.names = F)
rt8_G=rt8_G[,c("sample","age_at_index.demographic")]
rt9_G=read.table(file = "TARGET-WT.clinical.tsv",header = T,sep = "\t",check.names = F)
rt9_G=rt9_G[,c("sample","age_at_index.demographic")]
rt_G=rbind(rt1_G,rt2_G)
rt_G=rbind(rt_G,rt3_G)
rt_G=rbind(rt_G,rt4_G)
rt_G=rbind(rt_G,rt5_G)
rt_G=rbind(rt_G,rt6_G)
rt_G=rbind(rt_G,rt7_G)
rt_G=rbind(rt_G,rt8_G)
rt_G=rbind(rt_G,rt9_G)
survival_data=na.omit(rt_G)
colnames(survival_data)=c("sample","age")
age=survival_data
age=na.omit(age)
age$age=ifelse(age$age<=6,"age≤6",
               ifelse(age$age<=12,"6<age≤12","age>12"))
age=merge(TARGET_cluster,age,by.x="Tag",by.y="sample")
age$age=factor(age$age,levels=c("age≤6","6<age≤12","age>12"))
age$status="Age"
WNT_cluster1=age[age$cluster=="panWSS1",]
WNT_cluster2=age[age$cluster=="panWSS2",]
WNT_cluster3=age[age$cluster=="panWSS3",]
WNT_cluster4=age[age$cluster=="panWSS4",]
age=age[,-c(1,4)]
age=as.data.frame(table(age))
contingency_table=acast(age, cluster ~ age, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤6" = "#90d7ec",   
                               "6<age≤12" = "#009ad6",   
                               "age>12" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TARGET(n=3467, Age)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = "none", 
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤6" = "#90d7ec",   
                               "6<age≤12" = "#009ad6",   
                               "age>12" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤6" = "#90d7ec",   
                               "6<age≤12" = "#009ad6",   
                               "age>12" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = age)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("age≤6" = "#90d7ec",   
                               "6<age≤12" = "#009ad6",   
                               "age>12" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# TARGET-gender
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/TARGET_cluster.Rdata")
rt1_G=read.table(file = "TARGET-AML.clinical.tsv",header = T,sep = "\t",check.names = F)
rt1_G=rt1_G[,c("sample","gender.demographic")]
rt2_G=read.table(file = "TARGET-ALL-P1.clinical.tsv",header = T,sep = "\t",check.names = F)
rt2_G=rt2_G[,c("sample","gender.demographic")]
rt3_G=read.table(file = "TARGET-ALL-P2.clinical.tsv",header = T,sep = "\t",check.names = F)
rt3_G=rt3_G[,c("sample","gender.demographic")]
rt4_G=read.table(file = "TARGET-ALL-P3.clinical.tsv",header = T,sep = "\t",check.names = F)
rt4_G=rt4_G[,c("sample","gender.demographic")]
rt5_G=read.table(file = "TARGET-CCSK.clinical.tsv",header = T,sep = "\t",check.names = F)
rt5_G=rt5_G[,c("sample","gender.demographic")]
rt6_G=read.table(file = "TARGET-NBL.clinical.tsv",header = T,sep = "\t",check.names = F)
rt6_G=rt6_G[,c("sample","gender.demographic")]
rt7_G=read.table(file = "TARGET-OS.clinical.tsv",header = T,sep = "\t",check.names = F)
rt7_G=rt7_G[,c("sample","gender.demographic")]
rt8_G=read.table(file = "TARGET-RT.clinical.tsv",header = T,sep = "\t",check.names = F)
rt8_G=rt8_G[,c("sample","gender.demographic")]
rt9_G=read.table(file = "TARGET-WT.clinical.tsv",header = T,sep = "\t",check.names = F)
rt9_G=rt9_G[,c("sample","gender.demographic")]
rt_G=rbind(rt1_G,rt2_G)
rt_G=rbind(rt_G,rt3_G)
rt_G=rbind(rt_G,rt4_G)
rt_G=rbind(rt_G,rt5_G)
rt_G=rbind(rt_G,rt6_G)
rt_G=rbind(rt_G,rt7_G)
rt_G=rbind(rt_G,rt8_G)
rt_G=rbind(rt_G,rt9_G)
survival_data=na.omit(rt_G)
survival_data=survival_data[!survival_data$gender.demographic%in%c("","not reported","unknown"),]
survival_data$gender.demographic=ifelse(survival_data$gender.demographic=="male","Male","Female")
colnames(survival_data)=c("sample","gender")
gender=survival_data
gender=merge(TARGET_cluster,gender,by.x="Tag",by.y="sample")
gender=gender[,c("cluster","gender")]
gender=na.omit(gender)
gender$status="gender"
WNT_cluster1=gender[gender$cluster=="panWSS1",]
WNT_cluster2=gender[gender$cluster=="panWSS2",]
WNT_cluster3=gender[gender$cluster=="panWSS3",]
WNT_cluster4=gender[gender$cluster=="panWSS4",]
gender=gender[,-3]
gender=as.data.frame(table(gender))
contingency_table=acast(gender, cluster ~ gender, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TARGET(n=3504, Gender)") + 
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = gender)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Male" = "#6EB9C3", "Female" = "#C98B88")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TARGET-race
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
library(reshape2)
load("input/TARGET_cluster.Rdata")
rt1_G=read.table(file = "TARGET-AML.clinical.tsv",header = T,sep = "\t",check.names = F)
rt1_G=rt1_G[,c("sample","race.demographic")]
rt2_G=read.table(file = "TARGET-ALL-P1.clinical.tsv",header = T,sep = "\t",check.names = F)
rt2_G=rt2_G[,c("sample","race.demographic")]
rt3_G=read.table(file = "TARGET-ALL-P2.clinical.tsv",header = T,sep = "\t",check.names = F)
rt3_G=rt3_G[,c("sample","race.demographic")]
rt4_G=read.table(file = "TARGET-ALL-P3.clinical.tsv",header = T,sep = "\t",check.names = F)
rt4_G=rt4_G[,c("sample","race.demographic")]
rt5_G=read.table(file = "TARGET-CCSK.clinical.tsv",header = T,sep = "\t",check.names = F)
rt5_G=rt5_G[,c("sample","race.demographic")]
rt6_G=read.table(file = "TARGET-NBL.clinical.tsv",header = T,sep = "\t",check.names = F)
rt6_G=rt6_G[,c("sample","race.demographic")]
rt7_G=read.table(file = "TARGET-OS.clinical.tsv",header = T,sep = "\t",check.names = F)
rt7_G=rt7_G[,c("sample","race.demographic")]
rt8_G=read.table(file = "TARGET-RT.clinical.tsv",header = T,sep = "\t",check.names = F)
rt8_G=rt8_G[,c("sample","race.demographic")]
rt9_G=read.table(file = "TARGET-WT.clinical.tsv",header = T,sep = "\t",check.names = F)
rt9_G=rt9_G[,c("sample","race.demographic")]
rt_G=rbind(rt1_G,rt2_G)
rt_G=rbind(rt_G,rt3_G)
rt_G=rbind(rt_G,rt4_G)
rt_G=rbind(rt_G,rt5_G)
rt_G=rbind(rt_G,rt6_G)
rt_G=rbind(rt_G,rt7_G)
rt_G=rbind(rt_G,rt8_G)
rt_G=rbind(rt_G,rt9_G)
survival_data=na.omit(rt_G)
survival_data=survival_data[survival_data$race.demographic%in%c("asian","black or african american","white"),]
race=survival_data
colnames(race)=c("sample","race")
race$race=ifelse(race$race=="white","White",ifelse(race$race=="asian","Asian","Black"))
race=merge(TARGET_cluster,race,by.x="Tag",by.y="sample")
race=race[,c("cluster","race")]
race=na.omit(race)
race$status="race"
WNT_cluster1=race[race$cluster=="panWSS1",]
WNT_cluster2=race[race$cluster=="panWSS2",]
WNT_cluster3=race[race$cluster=="panWSS3",]
WNT_cluster4=race[race$cluster=="panWSS4",]
race=race[,-3]
race=as.data.frame(table(race))
contingency_table=acast(race, cluster ~ race, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("TARGET(n=2968, Race)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = race)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Asian" = "#fedcbd","Black"="#d3d7d4","White"="#cde6c7")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(panWSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# TCGA pan-cancer for single cancers
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
genesets=as.matrix(pancancer_WNT_Score$single_score_matrix)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
genesets=genesets[rownames(pancancer_group),]
fs=unique(pancancer_group$Project)
for(i in 1:length(fs)){
  setwd("input")
  pancancer_group_single=pancancer_group[pancancer_group$Project==fs[i],]
  genesets_single=genesets[rownames(pancancer_group_single),]
  library(mlbench)
  library(factoextra)
  library(ggplot2)
  library (cluster)
  library(fpc)
  set.seed(123)
  km_result = kmeans(genesets_single, centers=4,nstart = 1000,algorithm="Hartigan-Wong",iter.max=100)
  panduan=as.data.frame(km_result$cluster)
  panduan$Tag=rownames(panduan)
  colnames(panduan)=c("cluster","Tag")
  rank1=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==1,]$Tag,]$activity_score)
  rank1
  rank2=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==2,]$Tag,]$activity_score)
  rank2
  rank3=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==3,]$Tag,]$activity_score)
  rank3
  rank4=median(pancancer_WNT_Score$final_activity_score[panduan[panduan$cluster==4,]$Tag,]$activity_score)
  rank4
  rank=data.frame(Tag=c(1,2,3,4),value=c(rank1, rank2, rank3, rank4))
  rank=rank[order(-rank$value),]
  km_result$cluster=ifelse(km_result$cluster==rank$Tag[1],1,ifelse(km_result$cluster==rank$Tag[2],2,ifelse(km_result$cluster==rank$Tag[3],3,4)))
  dd=cbind(genesets_single, cluster = km_result$cluster)
  p0=fviz_cluster(km_result, data = genesets_single,  
                  palette = c("#d71345","#f47920","#145b7d","#6950a1"),  
                  ellipse.type = "euclid",  
                  star.plot = TRUE,   
                  repel = TRUE,  
                  geom = "point",
                  pointsize = 0.5) +   
    ggtitle(paste0("TCGA_", fs[i])) +   
    theme_minimal() +  
    theme(  
      text = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black")
    )
  setwd("input")
  cluster=data.frame(Tag = rownames(as.data.frame(dd)),cluster=as.data.frame(dd)$cluster,stringsAsFactors = FALSE)
  cluster$cluster = paste("cluster",cluster$cluster, sep = "")
  cluster$cluster=ifelse(cluster$cluster=="cluster1",paste0(tolower(fs[i]),"WSS1"),
                         ifelse(cluster$cluster=="cluster2",paste0(tolower(fs[i]),"WSS2"),
                                ifelse(cluster$cluster=="cluster3",paste0(tolower(fs[i]),"WSS3"),paste0(tolower(fs[i]),"WSS4"))))
  WNT=merge(pancancer_WNT_Score$final_activity_score,cluster,by.x="row.names",by.y="Tag")
  WNT$project="TCGA"
  table(WNT$cluster)
  WNT$cluster=factor(WNT$cluster,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                          paste0(tolower(fs[i]),"WSS2"),
                                          paste0(tolower(fs[i]),"WSS3"),
                                          paste0(tolower(fs[i]),"WSS4")))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c(paste0(tolower(fs[i]),"WSS1"),
                       paste0(tolower(fs[i]),"WSS2"),
                       paste0(tolower(fs[i]),"WSS3"),
                       paste0(tolower(fs[i]),"WSS4")), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p1=ggplot(WNT, aes(x = cluster, y = activity_score, fill = cluster)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c(paste0(tolower(fs[i]),"WSS1\n(n=",table(WNT$cluster)[1],")"),   
                                paste0(tolower(fs[i]),"WSS2\n(n=",table(WNT$cluster)[2],")"),  
                                paste0(tolower(fs[i]),"WSS3\n(n=",table(WNT$cluster)[3],")"),  
                                paste0(tolower(fs[i]),"WSS4\n(n=",table(WNT$cluster)[4],")"))) +  
    theme_bw() +  
    theme(  
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16,  color = "black"),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, color = "black")
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
    ggtitle(paste0("TCGA-",fs[i]))  
  print(p1)
  setwd("input")
  assign(paste0(fs[i],"_cluster"),cluster)
  save(list=paste0(fs[i],"_cluster"),file=paste0(paste0(fs[i],"_cluster.Rdata")))
}
# ------------------------------------------------------------------------------
# TARGET pan-cancer for single cancers
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/TARGET_WNT_Score.Rdata")
genesets=as.matrix(TARGET_WNT_Score$single_score_matrix)
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
genesets=genesets[rownames(TARGET_G),]
fs=unique(TARGET_G$Project)
for(i in 1:length(fs)){
  setwd("input")
  TARGET_group_single=TARGET_G[TARGET_G$Project==fs[i],]
  genesets_single=genesets[rownames(TARGET_group_single),]
  library(mlbench)
  library(factoextra)
  library(ggplot2)
  library (cluster)
  library(fpc)
  set.seed(123)
  km_result = kmeans(genesets_single, centers=4,nstart = 1000,algorithm="Hartigan-Wong",iter.max=100)
  panduan=as.data.frame(km_result$cluster)
  panduan$Tag=rownames(panduan)
  colnames(panduan)=c("cluster","Tag")
  rank1=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==1,]$Tag,]$activity_score)
  rank1
  rank2=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==2,]$Tag,]$activity_score)
  rank2
  rank3=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==3,]$Tag,]$activity_score)
  rank3
  rank4=median(TARGET_WNT_Score$final_activity_score[panduan[panduan$cluster==4,]$Tag,]$activity_score)
  rank4
  rank=data.frame(Tag=c(1,2,3,4),value=c(rank1, rank2, rank3, rank4))
  rank=rank[order(-rank$value),]
  km_result$cluster=ifelse(km_result$cluster==rank$Tag[1],1,ifelse(km_result$cluster==rank$Tag[2],2,ifelse(km_result$cluster==rank$Tag[3],3,4)))
  dd=cbind(genesets_single, cluster = km_result$cluster)
  p0=fviz_cluster(km_result, data = genesets_single,  
                  palette = c("#d71345","#f47920","#145b7d","#6950a1"),  
                  ellipse.type = "euclid",  
                  star.plot = TRUE,   
                  repel = TRUE,  
                  geom = "point",
                  pointsize = 0.5) +   
    ggtitle(paste0("TARGET_", fs[i])) +   
    theme_minimal() +  
    theme(  
      text = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black")
    )
  setwd("input")
  cluster=data.frame(Tag = rownames(as.data.frame(dd)),cluster=as.data.frame(dd)$cluster,stringsAsFactors = FALSE)
  cluster$cluster = paste("cluster",cluster$cluster, sep = "")
  cluster$cluster=ifelse(cluster$cluster=="cluster1",paste0(tolower(fs[i]),"WSS1"),
                         ifelse(cluster$cluster=="cluster2",paste0(tolower(fs[i]),"WSS2"),
                                ifelse(cluster$cluster=="cluster3",paste0(tolower(fs[i]),"WSS3"),paste0(tolower(fs[i]),"WSS4"))))
  WNT=merge(TARGET_WNT_Score$final_activity_score,cluster,by.x="row.names",by.y="Tag")
  WNT$project="TCGA"
  table(WNT$cluster)
  WNT$cluster=factor(WNT$cluster,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                          paste0(tolower(fs[i]),"WSS2"),
                                          paste0(tolower(fs[i]),"WSS3"),
                                          paste0(tolower(fs[i]),"WSS4")))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c(paste0(tolower(fs[i]),"WSS1"),
                       paste0(tolower(fs[i]),"WSS2"),
                       paste0(tolower(fs[i]),"WSS3"),
                       paste0(tolower(fs[i]),"WSS4")), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p1=ggplot(WNT, aes(x = cluster, y = activity_score, fill = cluster)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c(paste0(tolower(fs[i]),"WSS1\n(n=",table(WNT$cluster)[1],")"),   
                                paste0(tolower(fs[i]),"WSS2\n(n=",table(WNT$cluster)[2],")"),  
                                paste0(tolower(fs[i]),"WSS3\n(n=",table(WNT$cluster)[3],")"),  
                                paste0(tolower(fs[i]),"WSS4\n(n=",table(WNT$cluster)[4],")"))) +  
    theme_bw() +  
    theme(  
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16,  color = "black"),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, color = "black")
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +  
    ggtitle(paste0("TARGET-",fs[i]))  
  print(p1)
  setwd("input")
  assign(paste0(fs[i],"_cluster"),cluster)
  save(list=paste0(fs[i],"_cluster"),file=paste0(paste0(fs[i],"_cluster.Rdata")))
}
# ------------------------------------------------------------------------------
# TCGA single cancer OS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)
p_TCGA=data.frame()
for(i in 1:length(fs)){
  setwd("input")
  load(paste0("input/",fs[i],"_cluster.Rdata"))
  cluster=get(paste0(fs[i],"_cluster"))
  survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
  survival_data=survival_data[,c("sample","OS","OS.time")]
  survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
  survival_data_cms$level=survival_data_cms$cluster
  survival_data_cms$level=factor(survival_data_cms$level,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                                                  paste0(tolower(fs[i]),"WSS2"),
                                                                  paste0(tolower(fs[i]),"WSS3"),
                                                                  paste0(tolower(fs[i]),"WSS4")))
  survival_data_cms$OS.time=survival_data_cms$OS.time/30
  survival_data_cms=na.omit(survival_data_cms)
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(OS.time, OS) ~ level,data= survival_data_cms,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p_TCGA_single=data.frame(Cancer=fs[i],p_value=fitd[["pvalue"]])
  p_TCGA=rbind(p_TCGA,p_TCGA_single)
  fit = survfit(Surv(OS.time, OS)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
  ps = pairwise_survdiff(Surv(OS.time, OS)~ level,data= survival_data_cms,p.adjust.method = "none")
  mycol=c("#d71345","#f47920","#145b7d","#6950a1")
  names(fit$strata)=gsub("group=", "", names(fit$strata))
  library(RColorBrewer)
  library(tibble)
  library(ggpp)
  p=ggsurvplot(  
    fit = fit,  
    conf.int = FALSE,
    risk.table = F,
    risk.table.col = "strata",  
    palette = mycol,
    data = survival_data_cms,  
    xlim = c(0, max(survival_data_cms$OS.time)),
    size = 0.9,  
    break.time.by = 48,
    legend.title = "",
    xlab = "Time (months)",
    ylab = "Probability",
    risk.table.y.text = FALSE,
    tables.height = 0.25,
    title = paste0("OS(TCGA-",fs[i],",n=",length(survival_data_cms$sample),")"),
    ylim = c(0, 1)
  )  
  p.lab = paste0("log-rank test p",  
                  ifelse(p.val < 0.001, " < 0.001",
                         paste0(" = ", round(p.val, 3))))  
  p$plot = p$plot +   
    annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
             fontface = "italic", label = p.lab, size =6) +
    theme(text = element_text(size = 16))
  setwd("./supplementary_figures")
  ggsave(filename = paste0("TCGA-", fs[i], "_os.pdf"), plot = p$plot, width=7, height=4.5)
}
# ------------------------------------------------------------------------------
# TARGET single cancer OS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
fs=unique(TARGET_G$Project)
p_TARGET=data.frame()
for(i in 1:length(fs)){
  setwd("input")
  load(paste0("input/",fs[i],"_cluster.Rdata"))
  cluster=get(paste0(fs[i],"_cluster"))
  survival_data1=read.delim("TARGET-AML.survival.tsv")
  survival_data2=read.delim("TARGET-ALL-P1.survival.tsv")
  survival_data3=read.delim("TARGET-ALL-P2.survival.tsv")
  survival_data4=read.delim("TARGET-ALL-P3.survival.tsv")
  survival_data5=read.delim("TARGET-CCSK.survival.tsv")
  survival_data6=read.delim("TARGET-NBL.survival.tsv")
  survival_data7=read.delim("TARGET-RT.survival.tsv")
  survival_data8=read.delim("TARGET-WT.survival.tsv")
  survival_data9=read.delim("TARGET-OS.survival.tsv")
  survival_data=rbind(survival_data1,survival_data2)
  survival_data=rbind(survival_data,survival_data3)
  survival_data=rbind(survival_data,survival_data4)
  survival_data=rbind(survival_data,survival_data5)
  survival_data=rbind(survival_data,survival_data6)
  survival_data=rbind(survival_data,survival_data7)
  survival_data=rbind(survival_data,survival_data8)
  survival_data=rbind(survival_data,survival_data9)
  survival_data=survival_data[,c(1,2,3)]
  survival_data=na.omit(survival_data)
  survival_data=survival_data[,c("sample","OS","OS.time")]
  survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
  survival_data_cms$level=survival_data_cms$cluster
  survival_data_cms$level=factor(survival_data_cms$level,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                                                  paste0(tolower(fs[i]),"WSS2"),
                                                                  paste0(tolower(fs[i]),"WSS3"),
                                                                  paste0(tolower(fs[i]),"WSS4")))
  survival_data_cms$OS.time=survival_data_cms$OS.time/30
  survival_data_cms=na.omit(survival_data_cms)
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(OS.time, OS) ~ level,data= survival_data_cms,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p_TARGET_single=data.frame(Cancer=fs[i],p_value=fitd[["pvalue"]])
  p_TARGET=rbind(p_TARGET,p_TARGET_single)
  fit = survfit(Surv(OS.time, OS)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
  ps = pairwise_survdiff(Surv(OS.time, OS)~ level,data= survival_data_cms,p.adjust.method = "none")
  mycol=c("#d71345","#f47920","#145b7d","#6950a1")
  names(fit$strata)=gsub("group=", "", names(fit$strata))
  library(RColorBrewer)
  library(tibble)
  library(ggpp)
  p=ggsurvplot(  
    fit = fit,  
    conf.int = FALSE, 
    risk.table = F, 
    risk.table.col = "strata",  
    palette = mycol, 
    data = survival_data_cms,  
    xlim = c(0, max(survival_data_cms$OS.time)),
    size = 0.9,  
    break.time.by = 48,  
    legend.title = "",
    xlab = "Time (months)",
    ylab = "Probability",
    risk.table.y.text = FALSE, 
    tables.height = 0.25,
    title = paste0("OS(TARGET-",fs[i],",n=",length(survival_data_cms$sample),")"),
    ylim = c(0, 1)
  )  
  p.lab = paste0("log-rank test p",  
                  ifelse(p.val < 0.001, " < 0.001",
                         paste0(" = ", round(p.val, 3))))  
  p$plot = p$plot +   
    annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
             fontface = "italic", label = p.lab, size =6) +
    theme(text = element_text(size = 16))
  setwd("./supplementary_figures")
  ggsave(filename = paste0("TARGET-", fs[i], "_os.pdf"), plot = p$plot, width=7, height=4.5)
}
# ------------------------------------------------------------------------------
# TCGA single cancer DSS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)[-c(6)]
p_TCGA=data.frame()
for(i in 1:length(fs)){
  setwd("input")
  load(paste0("input/",fs[i],"_cluster.Rdata"))
  cluster=get(paste0(fs[i],"_cluster"))
  survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
  survival_data=survival_data[,c("sample","DSS","DSS.time")]
  survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
  survival_data_cms$level=survival_data_cms$cluster
  survival_data_cms$level=factor(survival_data_cms$level,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                                                  paste0(tolower(fs[i]),"WSS2"),
                                                                  paste0(tolower(fs[i]),"WSS3"),
                                                                  paste0(tolower(fs[i]),"WSS4")))
  survival_data_cms$DSS.time=survival_data_cms$DSS.time/30
  survival_data_cms=na.omit(survival_data_cms)
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(DSS.time, DSS) ~ level,data= survival_data_cms,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p_TCGA_single=data.frame(Cancer=fs[i],p_value=fitd[["pvalue"]])
  p_TCGA=rbind(p_TCGA,p_TCGA_single)
  fit = survfit(Surv(DSS.time, DSS)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
  ps = pairwise_survdiff(Surv(DSS.time, DSS)~ level,data= survival_data_cms,p.adjust.method = "none")
  mycol=c("#d71345","#f47920","#145b7d","#6950a1")
  names(fit$strata)=gsub("group=", "", names(fit$strata))
  library(RColorBrewer)
  library(tibble)
  library(ggpp)
  p=ggsurvplot(  
    fit = fit,  
    conf.int = FALSE,
    risk.table = F,
    risk.table.col = "strata",  
    palette = mycol,
    data = survival_data_cms,  
    xlim = c(0, max(survival_data_cms$DSS.time)),
    size = 0.9,  
    break.time.by = 48,
    legend.title = "",
    xlab = "Time (months)",
    ylab = "Probability",
    risk.table.y.text = FALSE,
    tables.height = 0.25,
    title = paste0("DSS(TCGA-",fs[i],",n=",length(survival_data_cms$sample),")"),
    ylim = c(0, 1)
  )  
  p.lab = paste0("log-rank test p",  
                  ifelse(p.val < 0.001, " < 0.001",
                         paste0(" = ", round(p.val, 3))))  
  p$plot = p$plot +   
    annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
             fontface = "italic", label = p.lab, size =6) +
    theme(text = element_text(size = 16))
  setwd("./supplementary_figures")
  ggsave(filename = paste0("TCGA-", fs[i], "_dss.pdf"), plot = p$plot, width=7, height=4.5)
}
# ------------------------------------------------------------------------------
# TCGA single cancer DFI
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)[-c(6,7,17,20,28)]
p_TCGA=data.frame()
for(i in 1:length(fs)){
  setwd("input")
  load(paste0("input/",fs[i],"_cluster.Rdata"))
  cluster=get(paste0(fs[i],"_cluster"))
  survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
  survival_data=survival_data[,c("sample","DFI","DFI.time")]
  survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
  survival_data_cms$level=survival_data_cms$cluster
  survival_data_cms$level=factor(survival_data_cms$level,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                                                  paste0(tolower(fs[i]),"WSS2"),
                                                                  paste0(tolower(fs[i]),"WSS3"),
                                                                  paste0(tolower(fs[i]),"WSS4")))
  survival_data_cms$DFI.time=survival_data_cms$DFI.time/30
  survival_data_cms=na.omit(survival_data_cms)
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(DFI.time, DFI) ~ level,data= survival_data_cms,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p_TCGA_single=data.frame(Cancer=fs[i],p_value=fitd[["pvalue"]])
  p_TCGA=rbind(p_TCGA,p_TCGA_single)
  fit = survfit(Surv(DFI.time, DFI)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
  ps = pairwise_survdiff(Surv(DFI.time, DFI)~ level,data= survival_data_cms,p.adjust.method = "none")
  mycol=c("#d71345","#f47920","#145b7d","#6950a1")
  names(fit$strata)=gsub("group=", "", names(fit$strata))
  library(RColorBrewer)
  library(tibble)
  library(ggpp)
  p=ggsurvplot(  
    fit = fit,  
    conf.int = FALSE,
    risk.table = F,
    risk.table.col = "strata",  
    palette = mycol,
    data = survival_data_cms,  
    xlim = c(0, max(survival_data_cms$DFI.time)),
    size = 0.9,  
    break.time.by = 48,
    legend.title = "",
    xlab = "Time (months)",
    ylab = "Probability",
    risk.table.y.text = FALSE,
    tables.height = 0.25,
    title = paste0("DFI(TCGA-",fs[i],",n=",length(survival_data_cms$sample),")"),
    ylim = c(0, 1)
  )  
  p.lab = paste0("log-rank test p",  
                  ifelse(p.val < 0.001, " < 0.001",
                         paste0(" = ", round(p.val, 3))))  
  p$plot = p$plot +   
    annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
             fontface = "italic", label = p.lab, size =6) +
    theme(text = element_text(size = 16))
  setwd("./supplementary_figures")
  ggsave(filename = paste0("TCGA-", fs[i], "_dfi.pdf"), plot = p$plot, width=7, height=4.5)
}
# ------------------------------------------------------------------------------
# TCGA single cancer PFI
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(dplyr)
library(ggplot2)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)[-c(6)]
p_TCGA=data.frame()
for(i in 1:length(fs)){
  setwd("input")
  load(paste0("input/",fs[i],"_cluster.Rdata"))
  cluster=get(paste0(fs[i],"_cluster"))
  survival_data=read.delim("Survival_SupplementalTable_S1_20171025_xena_sp")
  survival_data=survival_data[,c("sample","PFI","PFI.time")]
  survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
  survival_data_cms$level=survival_data_cms$cluster
  survival_data_cms$level=factor(survival_data_cms$level,levels=c(paste0(tolower(fs[i]),"WSS1"),
                                                                  paste0(tolower(fs[i]),"WSS2"),
                                                                  paste0(tolower(fs[i]),"WSS3"),
                                                                  paste0(tolower(fs[i]),"WSS4")))
  survival_data_cms$PFI.time=survival_data_cms$PFI.time/30
  survival_data_cms=na.omit(survival_data_cms)
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(PFI.time, PFI) ~ level,data= survival_data_cms,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p_TCGA_single=data.frame(Cancer=fs[i],p_value=fitd[["pvalue"]])
  p_TCGA=rbind(p_TCGA,p_TCGA_single)
  fit = survfit(Surv(PFI.time, PFI)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
  ps = pairwise_survdiff(Surv(PFI.time, PFI)~ level,data= survival_data_cms,p.adjust.method = "none")
  mycol=c("#d71345","#f47920","#145b7d","#6950a1")
  names(fit$strata)=gsub("group=", "", names(fit$strata))
  library(RColorBrewer)
  library(tibble)
  library(ggpp)
  p=ggsurvplot(  
    fit = fit,  
    conf.int = FALSE,
    risk.table = F,
    risk.table.col = "strata",  
    palette = mycol,
    data = survival_data_cms,  
    xlim = c(0, max(survival_data_cms$PFI.time)),
    size = 0.9,  
    break.time.by = 48,
    legend.title = "",
    xlab = "Time (months)",
    ylab = "Probability",
    risk.table.y.text = FALSE,
    tables.height = 0.25,
    title = paste0("PFI(TCGA-",fs[i],",n=",length(survival_data_cms$sample),")"),
    ylim = c(0, 1)
  )  
  p.lab = paste0("log-rank test p",  
                  ifelse(p.val < 0.001, " < 0.001",
                         paste0(" = ", round(p.val, 3))))  
  p$plot = p$plot +   
    annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,   
             fontface = "italic", label = p.lab, size =6) +
    theme(text = element_text(size = 16))
  setwd("./supplementary_figures")
  ggsave(filename = paste0("TCGA-", fs[i], "_pfi.pdf"), plot = p$plot, width=7, height=4.5)
}