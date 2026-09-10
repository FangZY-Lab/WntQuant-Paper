# ==============================================================================
# WntQuant-Paper5 -- Immunological and stemness analyses
# ------------------------------------------------------------------------------
# Purpose:    Link Wnt activity to immune subtypes, MSI status, stemness signatures, and the tumor microenvironment.
# Inputs:     Pan-cancer Wnt activity scores and immune/stemness annotations.
# Outputs:    Immunological and stemness association results.
# NOTE:       Part of the code involves the core interests of the laboratory and is not open source at this time.
# ==============================================================================


rm(list=ls())
gc()
setwd("input")
Immune_subtype = read.table("Subtype_Immune_Model_Based.txt",header = T,sep = "\t",fill = F)
load("input/pancancer_WNT_Score.Rdata")
immune_exp=merge(pancancer_WNT_Score$final_activity_score,Immune_subtype,by.x="row.names",by.y = "sample")
immune_exp=immune_exp[,c("activity_score","Subtype_Immune_Model_Based")]
table(immune_exp$Subtype_Immune_Model_Based)
immune_exp$Subtype_Immune_Model_Based=factor(immune_exp$Subtype_Immune_Model_Based,
                                             levels=c("Wound Healing (Immune C1)","IFN-gamma Dominant (Immune C2)",
                                                      "Inflammatory (Immune C3)","Lymphocyte Depleted (Immune C4)",
                                                      "Immunologically Quiet (Immune C5)","TGF-beta Dominant (Immune C6)"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("Wound Healing (Immune C1)","IFN-gamma Dominant (Immune C2)",
                     "Inflammatory (Immune C3)","Lymphocyte Depleted (Immune C4)",
                     "Immunologically Quiet (Immune C5)","TGF-beta Dominant (Immune C6)"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggboxplot(immune_exp, x="Subtype_Immune_Model_Based", y="activity_score", 
            fill = "Subtype_Immune_Model_Based", width = 0.5,
            xlab = "",
            palette = c('#f58f98','#b2d235',"#1d953f","#90d7ec","#009ad6","#ef5b9c")) +
  ylab(label = "Wnt/β-catenin pathway activity score") +
  guides(fill="none") +
  scale_x_discrete(labels = c("Wound healing(Immune C1,n=2414)", 
                              "IFN-γ dominant(Immune C2,n=2588)",
                              "Inflammatory(Immune C3,n=2395)",
                              "Lymphocyte depleted(Immune C4,n=1156)",
                              "Immunologically quiet(Immune C5,n=385)",
                              "TGF-β dominant(Immune C6,n=180)")) +
  theme(axis.text = element_text(size=16, color="black"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, color="black"),
        axis.title = element_text(size=16, color="black"),
        plot.title = element_text(size=16, color="black")) +
  stat_compare_means(comparisons=my.comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 4) +
  ggtitle("TCGA") +
  coord_flip()
p
# ------------------------------------------------------------------------------
# MSI
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
WNT_Score=pancancer_WNT_Score$final_activity_score
WNT_Score=WNT_Score[as.numeric(substr(rownames(WNT_Score), 14, 15)) <= 10, ]
WNT_Score=as.matrix(WNT_Score)
rownames(WNT_Score)=substr(rownames(WNT_Score), 1, 12)
library(limma)
WNT_Score=as.data.frame(avereps(WNT_Score))
library(BiocOncoTK)
data(MSIsensor.10k)
data(patient_to_tumor_code)
names(MSIsensor.10k)[1] = "patient_barcode"
save(MSIsensor.10k,patient_to_tumor_code,file = "MSIsensor10k.Rdata")
load("MSIsensor10k.Rdata")
MSI=merge(MSIsensor.10k,patient_to_tumor_code,by = "patient_barcode")
MSI=MSI[,c(1,2,4)]
WNT_Score_MSI=merge(WNT_Score,MSI,by.x="row.names",by.y="patient_barcode")
WNT_Score_MSI$tumor_code=ifelse(WNT_Score_MSI$tumor_code%in%c("COAD","READ"),"CRC",WNT_Score_MSI$tumor_code)
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","CRC","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
cancer_all=data.frame()
for(i in 1:length(fs)){
  WNT_Score_fs=WNT_Score_MSI[WNT_Score_MSI$tumor_code==fs[i],]
  cor_test=cor.test(WNT_Score_fs$activity_score,WNT_Score_fs$MSIsensor.score,method="pearson")
  cancer=data.frame(cancer=fs[i],cor=cor_test$estimate,p.value=cor_test$p.value)
  cancer_all=rbind(cancer_all,cancer)
}
library(ggsci)
library(ggplot2)
Color = pal_npg(alpha = 0.7)(6)
cancer_all=cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05, 
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$cancer = paste0(cancer_all$Signature,cancer_all$cancer)
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) + 
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("MSI")
# ------------------------------------------------------------------------------
# TMB
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(TCGAmutations)
tcga_available()
rm(list=ls())
setwd("input/TMB")
fs = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
        "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV",
        "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
        "UCEC","UCS","UVM")
lapply(fs,function(x){
  library(TCGAmutations)
  print(x)
  maf = TCGAmutations::tcga_load(study = x)
  save(maf, file = paste0("input/TMB/TCGA_",x,"_MC3_maf.Rdata"))
})
rm(list=ls())
setwd("input/TMB")
fs = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
        "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV",
        "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
        "UCEC","UCS","UVM")
TCGAtmb = lapply(fs,function(x){
  library(maftools) 
  load(paste0("input/TMB/TCGA_",x,"_MC3_maf.Rdata"))
  TMB = data.frame(tmb(maf = maf))
  TMB$Tumor_Sample_Barcode <- as.character(TMB$Tumor_Sample_Barcode)
  TMB$Group = x
  print(x)
  return(TMB)
})
TCGAtmb = do.call(rbind,TCGAtmb)
head(TCGAtmb)
TCGAtmb$Tumor_Sample_Barcode=substr(TCGAtmb$Tumor_Sample_Barcode, 1, 15)
TMB=TCGAtmb[,c(1,3,5)]
rownames(TMB)=TMB[,1]
TMB=TMB[,-1]
save(TMB,file="TMB.Rdata")
load("input/pancancer_WNT_Score.Rdata")
WNT_Score=pancancer_WNT_Score$final_activity_score
WNT_Score_TMB=merge(WNT_Score,TMB,by = "row.names")
WNT_Score_TMB$Group=ifelse(WNT_Score_TMB$Group%in%c("COAD","READ"),"CRC",WNT_Score_TMB$Group)
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","CRC","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
cancer_all=data.frame()
for(i in 1:length(fs)){
  WNT_Score_fs=WNT_Score_TMB[WNT_Score_TMB$Group==fs[i],]
  cor_test=cor.test(WNT_Score_fs$activity_score,WNT_Score_fs$total_perMB,method="pearson")
  cancer=data.frame(cancer=fs[i],cor=cor_test$estimate,p.value=cor_test$p.value)
  cancer_all=rbind(cancer_all,cancer)
}
library(ggsci)
library(ggplot2)
Color = pal_npg(alpha = 0.7)(6)
cancer_all=cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05, 
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$cancer = paste0(cancer_all$Signature,cancer_all$cancer)
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
library(GGally)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) + 
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("TMB")
# ------------------------------------------------------------------------------
# FGA
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
WNT_Score=pancancer_WNT_Score$final_activity_score
library(cBioPortalData)
cbio=cBioPortal()
studies=getStudies(cbio)
study=getStudies(api = cbio)
study=study$studyId[grep("tcga", study$studyId, ignore.case = TRUE)]
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COADREAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
pancancer_df=data.frame()
for (i in 1:length(fs)) {
  clinical=clinicalData(api = cbio, studyId = paste0(tolower(fs[i]),"_tcga"))
  df=clinical[,c("sampleId","FRACTION_GENOME_ALTERED")]
  df$Group=fs[i]
  pancancer_df=rbind(pancancer_df,df)
}
pancancer_df=na.omit(pancancer_df)
save(pancancer_df,file="FGA.Rdata")
colnames(pancancer_df)=c("sample","friction_genome_altered","Group")
pancancer_df$Group=ifelse(pancancer_df$Group=="COADREAD","CRC",pancancer_df$Group)
WNT_Score_TMB=merge(WNT_Score,pancancer_df,by.x="row.names",by.y = "sample")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","CRC","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
WNT_Score_TMB$friction_genome_altered=as.numeric(WNT_Score_TMB$friction_genome_altered)
cancer_all=data.frame()
for(i in 1:length(fs)){
  WNT_Score_fs=WNT_Score_TMB[WNT_Score_TMB$Group==fs[i],]
  cor_test=cor.test(WNT_Score_fs$activity_score,WNT_Score_fs$friction_genome_altered,method="pearson")
  cancer=data.frame(cancer=fs[i],cor=cor_test$estimate,p.value=cor_test$p.value)
  cancer_all=rbind(cancer_all,cancer)
}
library(ggsci)
library(ggplot2)
Color = pal_npg(alpha = 0.7)(6)
cancer_all=cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05, 
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$cancer = paste0(cancer_all$Signature,cancer_all$cancer)
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
library(GGally)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) + 
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("FGA")
# ------------------------------------------------------------------------------
# mutation count
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/pancancer_WNT_Score.Rdata")
WNT_Score=pancancer_WNT_Score$final_activity_score
library(cBioPortalData)
cbio=cBioPortal()
studies=getStudies(cbio)
study=getStudies(api = cbio)
study=study$studyId[grep("tcga", study$studyId, ignore.case = TRUE)]
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COADREAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
pancancer_df=data.frame()
for (i in 1:length(fs)) {
  clinical=clinicalData(api = cbio, studyId = paste0(tolower(fs[i]),"_tcga"))
  df=clinical[,c("sampleId","MUTATION_COUNT")]
  df$Group=fs[i]
  pancancer_df=rbind(pancancer_df,df)
}
pancancer_df=na.omit(pancancer_df)
save(pancancer_df,file="mutation_count.Rdata")
colnames(pancancer_df)=c("sample","MUTATION_COUNT","Group")
pancancer_df$Group=ifelse(pancancer_df$Group=="COADREAD","CRC",pancancer_df$Group)
WNT_Score_MT=merge(WNT_Score,pancancer_df,by.x="row.names",by.y = "sample")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","CRC","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
WNT_Score_MT$MUTATION_COUNT=log2(as.numeric(WNT_Score_MT$MUTATION_COUNT))
cancer_all=data.frame()
for(i in 1:length(fs)){
  WNT_Score_fs=WNT_Score_MT[WNT_Score_MT$Group==fs[i],]
  cor_test=cor.test(WNT_Score_fs$activity_score,WNT_Score_fs$MUTATION_COUNT,method="pearson")
  cancer=data.frame(cancer=fs[i],cor=cor_test$estimate,p.value=cor_test$p.value)
  cancer_all=rbind(cancer_all,cancer)
}
library(ggsci)
library(ggplot2)
Color = pal_npg(alpha = 0.7)(6)
cancer_all=cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05, 
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$cancer = paste0(cancer_all$Signature,cancer_all$cancer)
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
library(GGally)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) + 
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("Mutation count", subtitle = expression("(log"[2] * "-transformation)"))
# ------------------------------------------------------------------------------
# Immune gene collation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
signature_tme_literatures=signature_tme
load("input/tme_all_literatures.Rdata")
tme_all=as.data.frame(tme_all_literatures)
load("input/TARGET_G_all_literatures.Rdata")
tme_all=rbind(as.data.frame(TARGET_G_all_literatures),tme_all)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
sig=signature_collection_citation[c(19:54,127:141,282:314),]
tme_all=tme_all[,intersect(sig$Signatures,colnames(tme_all))]
sig=sig[sig$Signatures%in%colnames(tme_all),]
signature_tme_literatures=signature_tme_literatures[sig$Signatures]
for(i in 1:length(signature_tme_literatures)) {  
  signature_tme_literatures[[i]] <- c(paste0("PMID:",sig$PMID[i]), signature_tme_literatures[[i]])  
}  
load("input/ssGSEA28.Rdata")
for(i in 1:length(cellMarker)) {  
  cellMarker[[i]] <- c(paste0("PMID:28052254"), cellMarker[[i]])  
} 
write.gmt=function (gs.list, file) {
  if (length(grep(".gmt$", file)) == 0) 
    stop("write.gmt(): file should be .gmt")
  gs.names=names(gs.list)
  gs.desc=sapply(gs.list, function(x) x[1])
  gs.lines=unlist(lapply(gs.list, function(x) {
    paste(x[-1], collapse = "\t")
  }))
  gs.lines=paste(gs.names, gs.desc, gs.lines, sep = "\t")
  writeLines(gs.lines, con = file)
}
write.gmt(signature_tme_literatures,file="./signature_tme_literatures.gmt")
write.gmt(cellMarker,file="./cellMarker.gmt")
# ------------------------------------------------------------------------------
# TCGA immunoinfiltration analysis (literatures)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
load("input/pancancer_group.Rdata")
signature_tme_literatures=signature_tme
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_literatures=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_literatures=calculate_sig_score(eset = pancancer_exp_fs,signature = signature_tme_literatures,method="ssgsea",mini_gene_count=1)
  tme_all_literatures=rbind(tme_all_literatures,im_literatures)
}
save(tme_all_literatures,file="tme_all_literatures.Rdata")
# ------------------------------------------------------------------------------
# TARGET immunoinfiltration analysis (literatures)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
load("input/TARGET_G.Rdata")
signature_tme_literatures=signature_tme
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_G_all_literatures=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_literatures=calculate_sig_score(eset = pancancer_exp_fs,signature = signature_tme_literatures,method="ssgsea",mini_gene_count=1)
  TARGET_G_all_literatures=rbind(TARGET_G_all_literatures,im_literatures)
}
save(TARGET_G_all_literatures,file="TARGET_G_all_literatures.Rdata")
# ------------------------------------------------------------------------------
# WNT-immune infiltration correlation analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(IOBR)
load("input/tme_all_literatures.Rdata")
tme_all=as.data.frame(tme_all_literatures)
load("input/TARGET_G_all_literatures.Rdata")
tme_all=rbind(as.data.frame(TARGET_G_all_literatures),tme_all)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
sig=signature_collection_citation[c(19:54,127:141,282:314),]
tme_all=tme_all[,intersect(sig$Signatures,colnames(tme_all))]
sig=sig[sig$Signatures%in%colnames(tme_all),]
load("input/pancancer_WNT_score.Rdata")
load("input/TARGET_WNT_score.Rdata")
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
  group_fs=pancancer_group[pancancer_group==fs[i],]
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
# TCGA immunoinfiltration analysis (cellmarker)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
load("input/ssGSEA28.Rdata")
names(cellMarker)=paste0(gsub(" ", "_", names(cellMarker)), "_cellMarker")
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_cellMarker=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_cellMarker=calculate_sig_score(eset = pancancer_exp_fs,signature = cellMarker,method="ssgsea",mini_gene_count=1)
  tme_all_cellMarker=rbind(tme_all_cellMarker,im_cellMarker)
}
save(tme_all_cellMarker,file="tme_all_cellMarker.Rdata")
# ------------------------------------------------------------------------------
# TARGET immunoinfiltration analysis (cellmarker)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
load("input/ssGSEA28.Rdata")
names(cellMarker)=paste0(gsub(" ", "_", names(cellMarker)), "_cellMarker")
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_all_cellMarker=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_cellMarker=calculate_sig_score(eset = pancancer_exp_fs,signature = cellMarker,method="ssgsea",mini_gene_count=1)
  TARGET_all_cellMarker=rbind(TARGET_all_cellMarker,im_cellMarker)
}
save(TARGET_all_cellMarker,file="TARGET_all_cellMarker.Rdata")
# ------------------------------------------------------------------------------
# WNT-immunocellmarker correlation analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/tme_all_cellMarker.Rdata")
tme_all=as.data.frame(tme_all_cellMarker)
load("input/TARGET_all_cellMarker.Rdata")
tme_all=rbind(tme_all,as.data.frame(TARGET_all_cellMarker))
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
tme_all_cellMarker=tme_all[,grep("_cellMarker", colnames(tme_all), value = TRUE)]
load("input/pancancer_WNT_score.Rdata")
load("input/TARGET_WNT_score.Rdata")
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
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tme_all_fs=tme_all_cellMarker[rownames(group_fs),]
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
cancer_cor$tme=sub("_cellMarker$", "",cancer_cor$tme)
colnames(tme_all_fs)=sub("_cellMarker$", "",colnames(tme_all_fs))
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
# TCGA pan-cancer T-cell status analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(TCellSI)
setwd("input")
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_TCSS=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_TCSS=TCSS_Calculate(pancancer_exp_fs)
  im_TCSS=as.data.frame(t(im_TCSS))
  tme_all_TCSS=rbind(tme_all_TCSS,im_TCSS)
}
save(tme_all_TCSS,file="tme_all_TCSS.Rdata")
# ------------------------------------------------------------------------------
# TARGET pan-cancer T-cell status analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(TCellSI)
setwd("input")
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_all_TCSS=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_TCSS=TCSS_Calculate(pancancer_exp_fs)
  im_TCSS=as.data.frame(t(im_TCSS))
  TARGET_all_TCSS=rbind(TARGET_all_TCSS,im_TCSS)
}
save(TARGET_all_TCSS,file="TARGET_all_TCSS.Rdata")
# ------------------------------------------------------------------------------
# WNT-T cell status correlation analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/tme_all_TCSS.Rdata")
load("input/TARGET_all_TCSS.Rdata")
T_status=tme_all_TCSS
T_status=rbind(T_status,TARGET_all_TCSS)
load("input/pancancer_WNT_score.Rdata")
load("input/TARGET_WNT_score.Rdata")
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
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tide_all_fs=T_status[intersect(rownames(group_fs),rownames(T_status)),]
  pancancer_wnt_fs=pancancer_wnt[rownames(tide_all_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tide_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tide_all_fs[,j],method="pearson")
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
p=cancer_cor[,c(1,2,4)]
library(reshape2)
p=dcast(p, cancer ~ tme, value.var = "p.value")
rownames(p)=p[,1]
p=p[,-1]
m=cancer_cor[,c(1,2,3)]
library(reshape2)
m=dcast(m, cancer ~ tme, value.var = "cor")
rownames(m)=m[,1]
m=m[,-1]
tmp=cancer_cor[,c(1,2,5)]
library(reshape2)
tmp=dcast(tmp, cancer ~ tme, value.var = "pstar")
rownames(tmp)=tmp[,1]
tmp=tmp[,-1]
Data_GSVA = data.frame(t(m))
Data_GSVA$Term = rownames(Data_GSVA) 
library(tidyverse)
dft = Data_GSVA %>% pivot_longer(-Term)
Data_GSVA_P = data.frame(t(tmp))
Data_GSVA_P$Term = rownames(Data_GSVA_P)
dft_P = Data_GSVA_P %>% pivot_longer(-Term)
range(dft$value)
Data_GSVA_P = data.frame(t(p))
Data_GSVA_P$Term = rownames(Data_GSVA_P) 
dft_P = Data_GSVA_P %>% pivot_longer(-Term)
Plotdata = dft
colnames(Plotdata)[3] = "Correlation"
Plotdata$Pvalue = dft_P$value
Plotdata$Pvalue[Plotdata$Pvalue < 10^(-10)] = 10^(-10)
Plotdata$Pvalue = -log10(Plotdata$Pvalue)
Plotdata$Term=factor(Plotdata$Term,levels=c("Quiescence","Regulating","Proliferation","Helper","Cytotoxicity","Progenitor_exhaustion","Terminal_exhaustion","Senescence"))
Dotplot = ggplot(data = Plotdata, aes(x = name, y = Term, size = Pvalue)) +  
  geom_point(shape = 21, aes(fill = Correlation), position = position_dodge(0)) +  
  scale_size_continuous(range = c(0, 6),   
                        name = expression(-log[10](p)),
                        limits = c(0, 10),   
                        breaks = seq(0, 10, 2),   
                        labels = seq(0, 10, 2))  +  
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",   
                       name = "Correlation",  
                       limits = c(-0.85, 0.85), midpoint = 0, breaks = seq(-0.8, 0.8, 0.4)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_minimal() +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),  
        panel.grid = element_blank(),  
        plot.title = element_text(hjust = 0.5, size = 10, colour = 'black'),
        axis.ticks.y = element_blank(),  
        axis.title = element_blank(),  
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 16),
        axis.text.y = element_text(colour = 'black', size = 16))
Dotplot2 = Dotplot + theme(legend.position = "bottom",  
                            legend.text = element_text(color = "black", size = 16),
                            legend.title = element_text(size = 16, color = "black")) +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5,  
                               barwidth = 15, barheight = 1.5, ticks = TRUE))  
Dotplot2 
# ------------------------------------------------------------------------------
# TCGA immuno-infiltration analysis (ESTMATE)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_estimate=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_estimate=deconvo_tme(eset = pancancer_exp_fs, method = "estimate")
  tme_all_estimate=rbind(tme_all_estimate,im_estimate)
}
save(tme_all_estimate,file="tme_all_estimate.Rdata")
# ------------------------------------------------------------------------------
# TARGET immuno-infiltration analysis (ESTMATE)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_all_estimate=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_estimate=deconvo_tme(eset = pancancer_exp_fs, method = "estimate")
  TARGET_all_estimate=rbind(TARGET_all_estimate,im_estimate)
}
save(TARGET_all_estimate,file="TARGET_all_estimate.Rdata")
# ------------------------------------------------------------------------------
# ESTIMATE correlation analysis of WNT-immunity
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/tme_all_estimate.Rdata")
tme_all=as.data.frame(tme_all_estimate)
load("input/TARGET_all_estimate.Rdata")
tme_all=rbind(as.data.frame(TARGET_all_estimate),tme_all)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("input/pancancer_WNT_score.Rdata")
load("input/TARGET_WNT_score.Rdata")
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
  group_fs=pancancer_group[pancancer_group==fs[i],]
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
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "   ", 
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*  ", 
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "** ", "***")))
cancer_cor$cancer=paste0(cancer_cor$pstar,cancer_cor$cancer)
# ------------------------------------------------------------------------------
# Stromal
# ------------------------------------------------------------------------------
cancer_all=cancer_cor[cancer_cor$tme=="StromalScore_estimate",]
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +  
  geom_bar(stat = 'identity') + 
  ylab('') +  
  xlab("Correlation coefficient") +
  scale_fill_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +  
  theme_bw() +
  scale_x_continuous(breaks = c(-0.5, 0)) +
  theme(  
    panel.grid.major.y = element_blank(),  
    panel.grid.major.x = element_blank(),  
    legend.position = "none",   
    plot.title = element_text(size = 16, colour = 'black'),
    axis.title.x = element_text(size = 16, colour = 'black'),
    axis.title.y = element_text(size = 16, colour = 'black'),
    axis.text.x = element_text(size = 16, colour = 'black'),
    axis.text.y = element_text(size = 16, colour = 'black')
  ) +  
  ggtitle("Stromal score")
# ------------------------------------------------------------------------------
# Immune
# ------------------------------------------------------------------------------
cancer_all=cancer_cor[cancer_cor$tme=="ImmuneScore_estimate",]
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +  
  geom_bar(stat = 'identity') + 
  ylab('') +  
  xlab('Correlation coefficient') +  
  scale_fill_manual(values = c('#90d7ec', "#d3d7d4")) +  
  theme_bw() +
  scale_x_continuous(breaks = c(-0.5, 0)) +
  theme(
    panel.grid.major.y = element_blank(),  
    panel.grid.major.x = element_blank(),  
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  ggtitle("Immune score")
# ------------------------------------------------------------------------------
# Tumor purity
# ------------------------------------------------------------------------------
cancer_all=cancer_cor[cancer_cor$tme=="TumorPurity_estimate",]
cancer_all$cancer = factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
cancer_all$relation = ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +  
  geom_bar(stat = 'identity') +
  ylab('') +  
  xlab('Correlation coefficient') +  
  scale_fill_manual(values = c("#d3d7d4", '#f58f98')) +  
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.5)) +
  theme(
    panel.grid.major.y = element_blank(),  
    panel.grid.major.x = element_blank(),  
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  ggtitle("Tumor purity")
# ------------------------------------------------------------------------------
# Correlation analysis of WNT-Immunomodulator/Chemokine
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
Chemokinegene=read.delim("Chemokine.txt",header=TRUE,sep='\t',check.names=F)
colnames(Chemokinegene)="gene"
Immunomodulatorgene=read.delim("Immunomodulator.txt",header=TRUE,sep='\t',check.names=F)
colnames(Immunomodulatorgene)="gene"
Chemokine_Immunomodulator=rbind(Immunomodulatorgene,Chemokinegene)
Chemokine_Immunomodulator=as.data.frame(unique(Chemokine_Immunomodulator$gene))
Chemokine_Immunomodulator=Chemokine_Immunomodulator[order(Chemokine_Immunomodulator$`unique(Chemokine_Immunomodulator$gene)`),,drop=F]
load("input/pancancer_exp.Rdata")
load("input/TARGET.Rdata")
pancancer_exp=merge(pancancer_exp,TARGET,by="row.names")
rownames(pancancer_exp)=pancancer_exp[,1]
pancancer_exp=pancancer_exp[,-1]
wntgene_exp=pancancer_exp[intersect(Chemokine_Immunomodulator$`unique(Chemokine_Immunomodulator$gene)`,rownames(pancancer_exp)),]
wntgene_exp=as.data.frame(t(wntgene_exp))
load("input/pancancer_WNT_score.Rdata")
load("input/TARGET_WNT_score.Rdata")
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
# Collation of available pan-cancer TIDE data
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
folder_path = "input/TIDE/Results/Tumor_Dysf_Excl_scores"    
file_names = list.files(folder_path, pattern = "TCGA.*\\.txt$", full.names = TRUE)  
data_list = lapply(file_names, read.table, header = TRUE)  
tide_all = do.call(rbind, data_list)
setwd("input")
save(tide_all,file="tide_all.Rdata")
# ------------------------------------------------------------------------------
# Correlation analysis of WNT-tide
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(limma)
load("input/tide_all.Rdata")
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_WNT_Score.Rdata")
pancancer_WNT_score=pancancer_WNT_Score$final_activity_score[rownames(pancancer_group),]
pancancer_WNT_score=as.matrix(pancancer_WNT_score)
rownames(pancancer_WNT_score)=substr(rownames(pancancer_WNT_score), 1, 12)
pancancer_WNT_score=as.data.frame(avereps(pancancer_WNT_score))
pancancer_group$sample=rownames(pancancer_group)
pancancer_group$sample=substr(pancancer_group$sample, 1, 12)
pancancer_group=pancancer_group[!duplicated(pancancer_group$sample), ]
rownames(pancancer_group)=pancancer_group$sample
co_sample=intersect(rownames(pancancer_group),rownames(tide_all))
pancancer_group=pancancer_group[co_sample,]
pancancer_WNT_score=pancancer_WNT_score[co_sample,]
tide_all=tide_all[co_sample,]
pancancer_wnt=pancancer_WNT_score
T_status=tide_all[,c("Dysfunction","Exclusion")]
fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tide_all_fs=T_status[intersect(rownames(group_fs),rownames(T_status)),]
  pancancer_wnt_fs=pancancer_wnt[rownames(tide_all_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tide_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tide_all_fs[,j],method="pearson")
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
p=cancer_cor[,c(1,2,4)]
library(reshape2)
p=dcast(p, cancer ~ tme, value.var = "p.value")
rownames(p)=p[,1]
p=p[,-1]
m=cancer_cor[,c(1,2,3)]
library(reshape2)
m=dcast(m, cancer ~ tme, value.var = "cor")
rownames(m)=m[,1]
m=m[,-1]
tmp=cancer_cor[,c(1,2,5)]
library(reshape2)
tmp=dcast(tmp, cancer ~ tme, value.var = "pstar")
rownames(tmp)=tmp[,1]
tmp=tmp[,-1]
Data_GSVA = data.frame(t(m))
Data_GSVA$Term = rownames(Data_GSVA) 
library(tidyverse)
dft = Data_GSVA %>% pivot_longer(-Term)
Data_GSVA_P = data.frame(t(tmp))
Data_GSVA_P$Term = rownames(Data_GSVA_P)
dft_P = Data_GSVA_P %>% pivot_longer(-Term)
range(dft$value)
Data_GSVA_P = data.frame(t(p))
Data_GSVA_P$Term = rownames(Data_GSVA_P) 
dft_P = Data_GSVA_P %>% pivot_longer(-Term)
Plotdata = dft
colnames(Plotdata)[3] = "Correlation"
Plotdata$Pvalue = dft_P$value
Plotdata$Pvalue[Plotdata$Pvalue < 10^(-10)] = 10^(-10)
Plotdata$Pvalue = -log10(Plotdata$Pvalue)
Plotdata$Term=factor(Plotdata$Term,levels=c("Dysfunction","Exclusion"))
Dotplot = ggplot(data = Plotdata, aes(x = name, y = Term, size = Pvalue)) +  
  geom_point(shape = 21, aes(fill = Correlation), position = position_dodge(0)) +  
  scale_size_continuous(range = c(0, 6),   
                        name = expression(-log[10](p)),
                        limits = c(0, 10),   
                        breaks = seq(0, 10, 2),   
                        labels = seq(0, 10, 2))  +  
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",   
                       name = "Correlation",  
                       limits = c(-0.85, 0.85), midpoint = 0, breaks = seq(-0.8, 0.8, 0.4)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_minimal() +  
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),  
        panel.grid = element_blank(),  
        plot.title = element_text(hjust = 0.5, size = 10, colour = 'black'),
        axis.ticks.y = element_blank(),  
        axis.title = element_blank(),  
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 16),
        axis.text.y = element_text(colour = 'black', size = 16))
Dotplot2 = Dotplot + theme(legend.position = "bottom",  
                            legend.text = element_text(color = "black", size = 16),
                            legend.title = element_text(size = 16, color = "black")) +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5,  
                               barwidth = 15, barheight = 1.5, ticks = TRUE))  
Dotplot2 
# ------------------------------------------------------------------------------
# Calculation of IMvigor's WNT scores
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/IMvigor.Rdata")
IMvigor=log2(IMvigor+1)
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=IMvigor,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
IMvigor_WNT_score=WNT_Score
save(IMvigor_WNT_score,file="IMvigor_WNT_score.Rdata")
# ------------------------------------------------------------------------------
# PDCD1/CD274 correlations
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor.Rdata")
IMvigor_PD=IMvigor[c("PDCD1","CD274"),]
IMvigor_PD=as.data.frame(t(IMvigor_PD))
load("input/IMvigor_WNT_score.Rdata")
wnt_pd=merge(IMvigor_WNT_score,IMvigor_PD,by="row.names")
library(ggplot2)
library(ggExtra)
library(ggpubr)
p1=ggscatter(wnt_pd, x = "PDCD1", y = "activity_score",
             add = "reg.line",
             conf.int = TRUE, 
             add.params = list(color = "#d71345", fill = "#BBBDBE"),
             color = "#90d7ec", size = 2) +
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
  ggtitle("IMvigor") +
  xlab("PDCD1 expression") +
  ylab("Wnt/β-catenin pathway activity score")
p1
p2 = ggscatter(wnt_pd, x = "CD274", y = "activity_score",
                add = "reg.line",
                conf.int = TRUE, 
                add.params = list(color = "#d71345", fill = "#BBBDBE"),
                color = "#90d7ec", size = 2) +
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
  ggtitle("IMvigor") +
  xlab("CD274 expression") +
  ylab("Wnt/β-catenin pathway activity score")
p2
# ------------------------------------------------------------------------------
# Prognostic analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor_WNT_score.Rdata")
load("input/IMvigor210CoreBiologies.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score$final_activity_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]
wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`=="desert",]
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os","censOS")]
wnt_IMvigor_prognosis=na.omit(wnt_IMvigor_prognosis)
quartiles = quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os, censOS) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os, censOS)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os, censOS)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
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
  data = wnt_IMvigor_prognosis,  
  xlim = c(0, 24),
  size = 0.5,  
  break.time.by = 12,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = paste0("OS(n=",nrow(wnt_IMvigor_prognosis),",desert)"),
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
# Comparison of desert multi-indicator prognosis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
load("input/IMvigor.Rdata")
IMvigor=as.data.frame(t(IMvigor))
wnt_IMvigor=merge(IMvigor,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ] 
wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`=="desert",]
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
for(i in 1:length(gene)){
  wnt_IMvigor_prognosis=wnt_IMvigor[,c(gene[i],"os","censOS")]
  colnames(wnt_IMvigor_prognosis)=c("gene","os","censOS")
  wnt_IMvigor_prognosis=na.omit(wnt_IMvigor_prognosis)
  quartiles <- quantile(wnt_IMvigor_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)  
  wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$gene, breaks = quartiles, include.lowest = TRUE,   
                                    labels = c("Low","Medium","High"))  
  wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(os, censOS) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ] 
wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`=="desert",]
model=c("activity_score","pos_score","neg_score","FMOne mutation burden per MB","Neoantigen burden per MB")
model_p=c()
for(i in 1:length(model)){
  wnt_IMvigor_prognosis=wnt_IMvigor[,c(model[i],"os","censOS")]
  colnames(wnt_IMvigor_prognosis)=c("model","os","censOS")
  wnt_IMvigor_prognosis=na.omit(wnt_IMvigor_prognosis)
  quartiles = quantile(wnt_IMvigor_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)  
  wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$model, breaks = quartiles, include.lowest = TRUE,   
                                    labels = c("Low","Medium","High"))  
  wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(os, censOS) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score",
                    "FMOne mutation burden per MB","Neoantigen burden per MB")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(t(IMvigor)),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
wnt_IMvigor=merge(im_ssGSEA,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ] 
wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`=="desert",]
signature=c("CD8+ T effector","Pan-F-TBRS","Antigen processing machinery","Immune checkpoint")
signature_p=c()
for(i in 1:length(signature)){
  wnt_IMvigor_prognosis=wnt_IMvigor[,c(signature[i],"os","censOS")]
  colnames(wnt_IMvigor_prognosis)=c("signature","os","censOS")
  wnt_IMvigor_prognosis=na.omit(wnt_IMvigor_prognosis)
  quartiles = quantile(wnt_IMvigor_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)  
  wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$signature, breaks = quartiles, include.lowest = TRUE,   
                                    labels = c("Low","Medium","High"))  
  wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(os, censOS) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  signature_p=c(signature_p,p.val)
}
signature_p=as.data.frame(signature_p)
signature_p$signature=signature
colnames(signature_p)=c("value","ID")
p_all=rbind(gene_p,model_p)
p_all=rbind(p_all,signature_p)
p_all$value=-log10(p_all$value)
library(ggplot2)
y_intercept_value = -log10(0.05)  
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) +   
  labs(title = "",  
       subtitle = "",  
       x = "Immunotherapy indicators",  
       y = expression(-log[10](italic(p)))) +  
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),  
    legend.title = element_blank(),  
    legend.position = "none",  
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")
  ) +  
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +  
  ggtitle("IMvigor(desert,OS)") +  
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),   
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)
# ------------------------------------------------------------------------------
# Immunotherapy 1
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","Best Confirmed Overall Response")]
wnt_IMvigor=na.omit(wnt_IMvigor)
table(wnt_IMvigor$'Best Confirmed Overall Response')
wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Best Confirmed Overall Response`!="NE",]
wnt_IMvigor$'Best Confirmed Overall Response'=factor(wnt_IMvigor$'Best Confirmed Overall Response',
                                                     levels=c("CR","PR","SD","PD"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("CR","PR","SD","PD"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggboxplot(wnt_IMvigor,   
              x='Best Confirmed Overall Response',   
              y="activity_score",   
              fill = 'Best Confirmed Overall Response',   
              width = 0.5,  
              xlab = "",  
              palette = c('#d71345','#f58220','#90d7ec','#009ad6')) +  
  ylab(label = "Wnt/β-catenin pathway activity score") +  
  guides(fill="none") +  
  scale_x_discrete(labels = c("CR\n(n=25)",   
                              "PR\n(n=43)",  
                              "SD\n(n=63)",  
                              "PD\n(n=167)")) +  
  theme(axis.text = element_text(size=10, color = "black"),  
        axis.title = element_text(size=16, color = "black"),  
        plot.title = element_text(size=16, color = "black"),  
        legend.text = element_text(size=16, color = "black"),  
        legend.title = element_text(size=16, color = "black")) +  
  stat_compare_means(comparisons=my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.signif",  
                     size = 4) +  
  ggtitle("IMvigor") 
p
# ------------------------------------------------------------------------------
# Immunotherapy 2
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","binaryResponse")]
wnt_IMvigor=na.omit(wnt_IMvigor)
table(wnt_IMvigor$binaryResponse)
wnt_IMvigor$binaryResponse=factor(wnt_IMvigor$binaryResponse,
                                  levels=c("CR/PR","SD/PD"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("CR/PR","SD/PD"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggboxplot(wnt_IMvigor,   
              x='binaryResponse',   
              y="activity_score",   
              fill = 'binaryResponse',   
              width = 0.5,  
              xlab = "",  
              palette = c('#f58f98','#90d7ec')) +  
  ylab(label = "Wnt/β-catenin pathway activity score") +  
  guides(fill=FALSE) +  
  scale_x_discrete(labels = c("CR/PR\n(n=68)",   
                              "SD/PD\n(n=230)")) +  
  theme(axis.text = element_text(size=10, color = "black"),  
        axis.title = element_text(size=16, color = "black"),  
        plot.title = element_text(size=16, color = "black"),  
        legend.text = element_text(size=16, color = "black"),  
        legend.title = element_text(size=16, color = "black")) +  
  stat_compare_means(comparisons=my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.signif",  
                     size = 4) +  
  ggtitle("IMvigor") 
p
# ------------------------------------------------------------------------------
# Immunotherapy 3
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","IC Level")]
wnt_IMvigor=na.omit(wnt_IMvigor)
table(wnt_IMvigor$"IC Level")
wnt_IMvigor$"IC Level"=factor(wnt_IMvigor$"IC Level",
                              levels=c("IC0","IC1","IC2+"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("IC0","IC1","IC2+"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggboxplot(wnt_IMvigor,   
              x="IC Level",   
              y="activity_score",   
              fill = "IC Level",   
              width = 0.5,  
              xlab = "",  
              palette = c('#cde6c7','#45b97c',"#005831")) +  
  ylab(label = "Wnt/β-catenin pathway activity score") +  
  guides(fill=FALSE) +  
  scale_x_discrete(labels = c("IC0\n(n=97)",   
                              "IC1\n(n=132)",  
                              "IC2+\n(n=118)")) +  
  theme(axis.text = element_text(size=16, color = "black"),  
        axis.title = element_text(size=16, color = "black"),  
        plot.title = element_text(size=16, color = "black"),  
        legend.text = element_text(size=16, color = "black"),  
        legend.title = element_text(size=16, color = "black")) +  
  stat_compare_means(comparisons=my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.signif",  
                     size = 4) +  
  ggtitle("IMvigor")   
p
# ------------------------------------------------------------------------------
# Immunotherapy 4
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","TC Level")]
wnt_IMvigor=na.omit(wnt_IMvigor)
table(wnt_IMvigor$"TC Level")
wnt_IMvigor$"TC Level"=factor(wnt_IMvigor$"TC Level",
                              levels=c("TC0","TC1","TC2+"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("TC0","TC1","TC2+"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggboxplot(wnt_IMvigor,   
              x="TC Level",   
              y="activity_score",   
              fill = "TC Level",   
              width = 0.5,  
              xlab = "",  
              palette = c('#9b95c9','#6950a1',"#8552a1")) +  
  ylab(label = "Wnt/β-catenin pathway activity score") +  
  guides(fill=FALSE) +  
  scale_x_discrete(labels = c("TC0\n(n=275)",   
                              "TC1\n(n=22)",  
                              "TC2+\n(n=50)")) +  
  theme(axis.text = element_text(size=16, color="black"),  
        axis.title = element_text(size=16, color="black"),  
        plot.title = element_text(size=16, color="black"),  
        legend.text = element_text(size=16, color="black"),  
        legend.title = element_text(size=16, color="black")) +  
  stat_compare_means(comparisons=my.comparisons,  
                     method = "wilcox.test",  
                     label = "p.signif",  
                     size = 4) +  
  ggtitle("IMvigor")
p
# ------------------------------------------------------------------------------
# Immunotherapy 5
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","Immune phenotype")]
wnt_IMvigor=na.omit(wnt_IMvigor)
table(wnt_IMvigor$`Immune phenotype`)
wnt_IMvigor$`Immune phenotype`=factor(wnt_IMvigor$`Immune phenotype`,
                                      levels=c("desert","excluded","inflamed"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("desert","excluded","inflamed"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggboxplot(wnt_IMvigor,   
              x="Immune phenotype",   
              y="activity_score",   
              fill = "Immune phenotype",   
              width = 0.5,  
              xlab = "",  
              palette = c('#9b95c9','#dec674',"#84bf96")) +  
  ylab(label = "Wnt/β-catenin pathway activity score") +  
  guides(fill=FALSE) +  
  scale_x_discrete(labels = c("Desert\n(n=76)",   
                              "Excluded\n(n=134)",   
                              "Inflamed\n(n=74)")) +  
  theme(axis.text = element_text(size=16, color="black"),
        axis.title = element_text(size=16, color="black"),  
        plot.title = element_text(size=16, color="black"),  
        legend.text = element_text(size=16, color="black"),  
        legend.title = element_text(size=16, color="black")) +  
  stat_compare_means(comparisons=my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.signif",  
                     size = 4) +  
  ggtitle("IMvigor") 
p
# ------------------------------------------------------------------------------
# Immunotherapy for different immune phenotypes
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/IMvigor210CoreBiologies.Rdata")
load("input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","binaryResponse","Immune phenotype")]
colnames(wnt_IMvigor)=c("activity_score","Response","Immune phenotype")
wnt_IMvigor=na.omit(wnt_IMvigor)
wnt_IMvigor$`Immune phenotype`=gsub("^(\\w)(.*)", "\\U\\1\\L\\2", wnt_IMvigor$`Immune phenotype`, perl = TRUE)
p = ggplot(wnt_IMvigor, aes(`Immune phenotype`, activity_score, fill = Response, color = Response)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +  
  theme_bw() +  
  theme(  
    legend.position = "right",  
    panel.grid.major.x = element_line(linewidth = 1),  
    panel.grid.major.y = element_line(linewidth = 1),  
    panel.grid.minor.x = element_line(linewidth = 0.8),  
    panel.grid.minor.y = element_line(linewidth = 0.8),  
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
    axis.line.y = element_line(colour = 'black', linewidth = 1),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#f58f98", "#90d7ec")) +  
  stat_compare_means(aes(group = Response, label = ..p.signif..), method = "wilcox.test") +  
  ggtitle("IMvigor") 
p
# ------------------------------------------------------------------------------
# Immunotherapy GSE91061
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
exp=read.csv(file = "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv")
Cytolytic=read.delim(file = "GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt")
load("input/IMvigor210CoreBiologies.Rdata")
rm(exprSet)
rm(phenoData)
annoData1=annoData[,c(1,2)]
exp=merge(annoData1,exp,by.x="entrez_id",by.y="X")
exp=exp[exp$symbol!='',]
rownames(exp)=exp$symbol
exp=exp[,-c(1,2)]
g_l=annoData[,c("entrez_id","length")]
colnames(g_l)=c("gene_id","length")
library(org.Hs.eg.db)
g2s=toTable(org.Hs.egSYMBOL)
genelength =  merge(g_l,g2s,by='gene_id')
head(genelength)
genelength=genelength[,c(3,2)]
colnames(genelength) = c("gene","length")
counts = exp[rownames(exp) %in% genelength$gene,]
count = as.matrix(counts)
genelength[match(rownames(counts),
                 genelength$gene),"length"]
eff_length = genelength[match(rownames(counts),genelength$gene),"length"]/1000
x = count / eff_length
mat_tpm = t( t(x) / colSums(x) ) * 1e6 
rownames(mat_tpm) = rownames(count) 
mat_tpm[1:4,1:4]
table(rowSums(mat_tpm==0)== dim(mat_tpm)[2])
keep.exprs = mat_tpm > 0
head(keep.exprs)
table(rowSums(keep.exprs))
keep = rowSums(keep.exprs) > 0
summary(keep)
mat_tpm.filtered = mat_tpm[keep,]
rm(keep.exprs)
rm(eff_length)
mat_tpm.filtered[1:4,1:4]
GSE91061=as.data.frame(log2(mat_tpm.filtered+1))
GSE91061 = GSE91061[rowSums(GSE91061)>0,]
library(GSEABase)
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=GSE91061,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
GSE91061_WNT_score=WNT_Score$final_activity_score
gene=as.data.frame(t(GSE91061[c("PDCD1","CD274"),]))
WNT=merge(GSE91061_WNT_score,gene,by="row.names")
Cytolytic$V1=gsub("-", ".", Cytolytic$V1)
WNT=merge(WNT,Cytolytic,by.x="Row.names",by.y="V1")
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
p = ggscatter(WNT, x = "V2", y = "activity_score",  
              add = "reg.line", 
              conf.int = TRUE,   
              add.params = list(color="#d71345", fill = "#BBBDBE"),  
              color = "#90d7ec", size = 2) +
  stat_cor(method = "pearson",  
           label.x = 0,   
           label.y = -2,  
           size = 6) +  
  theme_classic() +  
  theme(axis.title.x = element_text(size = 16, color = "black"),   
        axis.title.y = element_text(size = 16, color = "black"),  
        axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),  
        plot.title = element_text(size = 16, color = "black")) +
  ggtitle("GSE91061") +  
  xlab("Cytolytic score") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  theme(legend.position = "none") 
p
p = ggscatter(WNT, x = "CD274", y = "activity_score",  
              add = "reg.line", 
              conf.int = TRUE,   
              add.params = list(color = "#d71345", fill = "#BBBDBE"),  
              color = "#90d7ec", size = 2) +
  stat_cor(method = "pearson",  
           label.x = 0,   
           label.y = -2,  
           size = 6) +  
  theme_classic() +  
  theme(axis.title.x = element_text(size = 16, color = "black"),   
        axis.title.y = element_text(size = 16, color = "black"),  
        axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),  
        plot.title = element_text(size = 16, color = "black")) + 
  ggtitle("GSE91061") +  
  xlab("CD274 expression") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  theme(legend.position = "none")  
p
p = ggscatter(WNT, x = "PDCD1", y = "activity_score",  
              add = "reg.line",
              conf.int = TRUE,   
              add.params = list(color = "#d71345", fill = "#BBBDBE"),  
              color = "#90d7ec", size = 2) +
  stat_cor(method = "pearson",  
           label.x = 0,   
           label.y = -2,  
           size = 6) +  
  theme_classic() +  
  theme(axis.title.x = element_text(size = 16, color = "black"),   
        axis.title.y = element_text(size = 16, color = "black"),  
        axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),  
        plot.title = element_text(size = 16, color = "black")) +
  ggtitle("GSE91061") +  
  xlab("PDCD1 expression") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  theme(legend.position = "none")
p
# ------------------------------------------------------------------------------
# Collection of stemness gene sets
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(openxlsx)
setwd("input")
stem1=read.xlsx("tumor_stemness_PNAS_2019_109_genes.xlsx", sheet=1)
stem2=read.xlsx("tumor_stemness_PNAS_2019_109_genes.xlsx", sheet=2)
stem3=read.xlsx("13287_2022_2913_MOESM2_ESM.xlsx", sheet=1)
stem1 = lapply(stem1, function(x) na.omit(x)) 
stem2 = lapply(stem2, function(x) na.omit(x)) 
stem3 = lapply(stem3, function(x) na.omit(x)) 
stem=append(stem1, stem2)
stem=append(stem, stem3)
write.gmt=function (gs.list, file) {
  if (length(grep(".gmt$", file)) == 0) 
    stop("write.gmt(): file should be .gmt")
  gs.names=names(gs.list)
  gs.desc=sapply(gs.list, function(x) x[1])
  gs.lines=unlist(lapply(gs.list, function(x) {
    paste(x[-1], collapse = "\t")
  }))
  gs.lines=paste(gs.names, gs.desc, gs.lines, sep = "\t")
  writeLines(gs.lines, con = file)
}
write.gmt(stem,file="input/stemness.gmt")
# ------------------------------------------------------------------------------
# TCGA pan-cancer stemness analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
library(GSEABase)
stemness=getGmt("stemness.gmt")
stemness_gmt=lapply(stemness, function(x){ x@geneIds })
names(stemness_gmt)=names(stemness)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_stemness=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_stemness=calculate_sig_score(eset = pancancer_exp_fs,signature = stemness_gmt,method="ssgsea",mini_gene_count=1)
  tme_all_stemness=rbind(tme_all_stemness,im_stemness)
}
save(tme_all_stemness,file="tme_all_stemness.Rdata")
# ------------------------------------------------------------------------------
# TARGET pan-cancer stemness analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(IOBR)
library(GSEABase)
stemness=getGmt("stemness.gmt")
stemness_gmt=lapply(stemness, function(x){ x@geneIds })
names(stemness_gmt)=names(stemness)
load("input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_all_stemness=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_stemness=calculate_sig_score(eset = pancancer_exp_fs,signature = stemness_gmt,method="ssgsea",mini_gene_count=1)
  TARGET_all_stemness=rbind(TARGET_all_stemness,im_stemness)
}
save(TARGET_all_stemness,file="TARGET_all_stemness.Rdata")
# ------------------------------------------------------------------------------
# Correlation analysis of stemness
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/tme_all_stemness.Rdata")
load("input/TARGET_all_stemness.Rdata")
tme_all=as.data.frame(tme_all_stemness)
tme_all=rbind(tme_all,as.data.frame(TARGET_all_stemness))
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("input/pancancer_WNT_score.Rdata")
load("input/TARGET_WNT_score.Rdata")
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
# Joint CELL article on stemness
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(readxl)
setwd("input")
RNA <- read_excel("1-s2.0-S0092867418303581-mmc1.xlsx", sheet = "StemnessScores_RNAexp")[,-c(2,3)]
load("input/pancancer_cluster.Rdata")
RNA$TCGAlong.id=substr(RNA$TCGAlong.id,1,15)
colnames(RNA)
WNT=merge(RNA,pancancer_cluster,by.x="TCGAlong.id",by.y="Tag")
WNT$project="TCGA"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = `EREG-mRNAsi`, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("EREG-mRNAsi") +  
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
rm(list=ls())
gc()
library(readxl)
setwd("input")
DNA=read_excel("1-s2.0-S0092867418303581-mmc1.xlsx", sheet = "StemnessScores_DNAmeth")[,-c(2,3)]
load("input/pancancer_cluster.Rdata")
DNA$TCGAlong.id=substr(DNA$TCGAlong.id,1,15)
colnames(DNA)
WNT=merge(DNA,pancancer_cluster,by.x="TCGAlong.id",by.y="Tag")
WNT$project="TCGA"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("panWSS1","panWSS2","panWSS3","panWSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("panWSS1","panWSS2","panWSS3","panWSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p=ggplot(WNT, aes(x = cluster, y = `ENHsi`, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("ENHsi") +  
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
# Multiple immunotherapy data sets and calculations
# ------------------------------------------------------------------------------
rm(list = ls())
gc()
setwd("input")
PRJNA482620_GBM=readRDS("input/GBM-PRJNA482620.Response.Rds")
GSE93157_LUSC=readRDS("input/LUSC-GSE93157.Response.Rds")
GSE78220_Melanoma=readRDS("input/Melanoma-GSE78220.Response.Rds")
GSE93157_Melanoma=readRDS("input/Melanoma-GSE93157.Response.Rds")
GSE96619_Melanoma=readRDS("input/Melanoma-GSE96619.Response.Rds")
phs000452_Melanoma=readRDS("input/Melanoma-phs000452.Response.Rds")
GSE93157_NSCLC=readRDS("input/nonsqNSCLC-GSE93157.Response.Rds")
GSE126044_NSCLC=readRDS("input/NSCLC_GSE126044.Response.Rds")
GSE135222_NSCLC=readRDS("input/NSCLC_GSE135222.Response.Rds")
PRJEB25780_STAD=readRDS("input/STAD-PRJEB25780.Response.Rds")
datasets=c("PRJNA482620_GBM","GSE78220_Melanoma","GSE96619_Melanoma","phs000452_Melanoma","GSE126044_NSCLC","GSE135222_NSCLC","PRJEB25780_STAD")
plots=list() 
for (i in 1:length(datasets)){
  exp=as.data.frame(get(datasets[i]))
  rownames(exp)=exp[,1]
  exp=exp[,-1]
  library(dplyr)
  exp = exp %>%  
    filter(rowSums(is.na(.)) <= ncol(.)/2)
  if(any(is.na(exp))){
    library(mice)
    exp=mice(exp)
    exp=complete(exp)
    exp=as.data.frame(exp)
  }
  exp=log2(exp+1)
  library(GSEABase)
  geneSets=getGmt("wpags_wpigs.gmt")
  source("functions/Calculate_bioactivity_scores.R")
  WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                         expression_profile=exp,
                                         foundation="relative_ssGSEA",
                                         activation_geneset=NA,
                                         inhibition_geneset=NA,
                                         geneSets_gmt=geneSets,
                                         min.sz=1,
                                         max.sz=10000,
                                         geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                             rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  WNT_score=WNT_Score$final_activity_score
  gene=as.data.frame(t(exp[c("PDCD1","CD274"),]))
  WNT_gene=merge(WNT_score,gene,by="row.names")
  library(ggplot2)
  library(ggExtra)
  library(ggpubr)
  library(gridExtra)
  p1 = ggscatter(WNT_gene[, c("PDCD1", "activity_score")],   
                 x = "PDCD1",   
                 y = "activity_score",  
                 add = "reg.line",
                 conf.int = TRUE,   
                 add.params = list(color = "#d71345", fill = "#BBBDBE"),  
                 color = "#90d7ec",   
                 size = 2) + 
    stat_cor(method = "pearson", size = 6,
             label.x = 0,   
             label.y = -2) +  
    theme_classic() +  
    theme(axis.title.x = element_text(size = 16, color = "black"),   
          axis.title.y = element_text(size = 16, color = "black"),  
          axis.text.x = element_text(size = 16, color = "black"),  
          axis.text.y = element_text(size = 16, color = "black"),  
          plot.title = element_text(size = 16, color = "black")) +
    ggtitle(paste0(datasets[i])) +  
    xlab("PDCD1 expression") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    theme(legend.position = "none")
  plots[[length(plots) + 1]] = p1
  p2 = ggscatter(WNT_gene[, c("CD274", "activity_score")],   
                 x = "CD274",   
                 y = "activity_score",  
                 add = "reg.line", 
                 conf.int = TRUE,   
                 add.params = list(color = "#d71345", fill = "#BBBDBE"),  
                 color = "#90d7ec",   
                 size = 2) + 
    stat_cor(method = "pearson", size = 6,
             label.x = 0,   
             label.y = -2) +  
    theme_classic() +  
    theme(axis.title.x = element_text(size = 16, color = "black"),   
          axis.title.y = element_text(size = 16, color = "black"),  
          axis.text.x = element_text(size = 16, color = "black"),  
          axis.text.y = element_text(size = 16, color = "black"),  
          plot.title = element_text(size = 16, color = "black")) + 
    ggtitle(paste0(datasets[i])) +  
    xlab("CD274 expression") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    theme(legend.position = "none") 
  p2 
  plots[[length(plots) +1]]=p2
}
grid.arrange(grobs = plots, ncol =4)
# ------------------------------------------------------------------------------
# Timer algorithm
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(IOBR)
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_timer=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_timer=deconvo_tme(eset = pancancer_exp_fs, method = "timer", group_list = rep(tolower(fs[i]),dim(pancancer_exp_fs)[2]))
  tme_all_timer=rbind(tme_all_timer,im_timer)
}
save(tme_all_timer,file="tme_all_timer.Rdata")
# ------------------------------------------------------------------------------
# Visualisation of Timer algorithm results
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/tme_all_timer.Rdata")
tme_all=as.data.frame(tme_all_timer)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("input/pancancer_WNT_score.Rdata")
load("input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
pancancer_group=pancancer_group[rownames(tme_all),]
fs=unique(pancancer_group$Project)
all_fs_all=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group==fs[i],]
  tme_all_fs=tme_all[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_WNT_score[rownames(group_fs),]
  all_fs=merge(pancancer_wnt_fs,tme_all_fs,by="row.names")
  #B(2,5),T4(2,6),T8(2,7),Neu(2,8),macro(2,9),DC(1,10),This can be replaced by
  all_fs=all_fs[,c(2,9)]
  library(dplyr)
  all_fs = all_fs %>%  
    mutate(Group = ntile(Macrophage_TIMER,3))
  all_fs$Group = factor(all_fs$Group, levels =1:3, labels = c("Low", "Medium", "High"))
  all_fs$project=fs[i]
  all_fs_all=rbind(all_fs_all,all_fs)
}
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
all_fs_all$Group=factor(all_fs_all$Group,levels = c("Low","Medium","High"))
p = ggplot(all_fs_all, aes(project, activity_score, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +
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
    axis.text.y = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black", hjust = 0)
  ) +
  scale_fill_manual(values = c("#90d7ec", "#77787b", "#f58f98")) +
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +
  #ggtitle(expression(paste("CD8"^"+", " T cell infiltration")))
  ggtitle("Macrophage infiltration")#This can be replaced by
p
# ------------------------------------------------------------------------------
# Preparation of lasso regression algorithm
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(IOBR)
setwd("input")
load("input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("input/IMvigor.Rdata")
IMvigor=log2(IMvigor+1)
im_xcell=deconvo_tme(eset = IMvigor, method = "xcell")
im_xcell=as.data.frame(im_xcell)
rownames(im_xcell)=im_xcell[,1]
im_xcell=im_xcell[,-1]
colnames(im_xcell)=gsub("_xCell$", "", colnames(im_xcell))
im_xcell=as.data.frame(scale(im_xcell))
phenoData=na.omit(phenoData[,6,drop=F])
xcell_pheno=merge(phenoData,im_xcell,by="row.names")
rownames(xcell_pheno)=xcell_pheno[,1]
xcell_pheno=xcell_pheno[,-1]
colnames(xcell_pheno)[1]=c("phenotype")
table(xcell_pheno$phenotype)
xcell_pheno$phenotype=ifelse(xcell_pheno$phenotype=="inflamed","hot","cold")
table(xcell_pheno$phenotype)
save(xcell_pheno,file="xcell_pheno.Rdata")
# ------------------------------------------------------------------------------
# lasso regression construction
# ------------------------------------------------------------------------------
rm(list = ls())
gc()
setwd("input")
load("input/xcell_pheno.Rdata")
library(glmnet)
library(caret)
library(pROC) 
library(ggplot2)
library(vip)
library(dplyr)
library(purrr)
library(ggthemes)
xcell_pheno$phenotype = ifelse(xcell_pheno$phenotype == "hot", 1, 0)
xcell_pheno$phenotype = as.factor(xcell_pheno$phenotype)
processed_data = xcell_pheno
colnames(processed_data)[1] = "target"
x = as.matrix(processed_data[, -1])
y = processed_data$target
processed_data$target_factor = factor(ifelse(processed_data$target == 1, "hot", "cold"),
                                       levels = c("hot", "cold"))
set.seed(123)
train_control = trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  selectionFunction = "best",
  savePredictions = TRUE,
  verboseIter = FALSE
)
lasso_model = train(
  x = x,
  y = processed_data$target_factor,
  method = "glmnet",
  family = "binomial",
  trControl = train_control,
  tuneGrid = expand.grid(
    alpha = 1,
    lambda = 10^seq(-3, 0, length = 100)
  ),
  metric = "ROC"
)
print(lasso_model$bestTune)
tune_results = lasso_model$results
ggplot(tune_results, aes(x = lambda, y = ROC)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = lasso_model$bestTune$lambda, 
             linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Lasso model tuning",
       subtitle = paste("Optimal λ:", round(lasso_model$bestTune$lambda, 5)),
       x = "λ (log scale)",
       y = "AUC (Cross-validation)") +
  theme_minimal(base_size = 14)
coef_data = as.matrix(coef(lasso_model$finalModel, s = lasso_model$finalModel$lambda))
coef_df = data.frame(
  Lambda = rep(lasso_model$finalModel$lambda, each = nrow(coef_data)),
  Variable = rep(rownames(coef_data), times = ncol(coef_data)),
  Value = as.vector(coef_data)
)
coef_df = coef_df[coef_df$Variable != "(Intercept)", ]
ggplot(coef_df, aes(x = log(Lambda), y = Value, color = Variable)) +
  geom_line(linewidth = 0.2) +
  geom_vline(xintercept = log(lasso_model$bestTune$lambda), 
             linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_reverse() +
  labs(title = "Lasso coefficient paths",
       x = "log(λ)",
       y = "Coefficient value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
final_coef = coef(lasso_model$finalModel, s = lasso_model$bestTune$lambda)
selected_vars = rownames(final_coef)[which(final_coef != 0)]
selected_vars = selected_vars[selected_vars != "(Intercept)"]
selected_vars
coef_df_selected = coef_df[coef_df$Variable %in% selected_vars, ]
ggplot(coef_df_selected, aes(x = log(Lambda), y = Value, color = Variable)) +
  geom_line(linewidth = 0.5) +
  geom_vline(
    xintercept = log(lasso_model$bestTune$lambda),
    linetype = "dashed", color = "red", linewidth = 0.5
  ) +
  scale_x_reverse() +
  labs(
    title = "",
    x = "log(λ)",
    y = "Coefficient Value"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
vip(lasso_model, num_features = 15, geom = "point", 
    aesthetics = list(size = 4, color = "steelblue")) +
  labs(title = "Top 15 important variables",
       y = "Importance") +
  theme_minimal(base_size = 14)
importance_df = vi(lasso_model)
top_15_vars = head(importance_df$Variable, 15)
top_15_vars
set.seed(123)
best_lambda = lasso_model$bestTune$lambda
fixed_lambda_grid = expand.grid(alpha = 1, lambda = best_lambda)
refit_control = trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE,
  verboseIter = FALSE
)
refit_model = train(
  x = x,
  y = processed_data$target_factor,
  method = "glmnet",
  family = "binomial",
  trControl = refit_control,
  tuneGrid = fixed_lambda_grid,
  metric = "ROC"
)
if (!is.null(refit_model$pred)) {
  cv_predictions = refit_model$pred %>% 
    arrange(rowIndex) %>%
    mutate(obs = factor(obs, levels = c("cold", "hot")))
  roc_results = cv_predictions %>%
    group_by(Resample) %>%
    summarise(
      roc_obj = list(suppressMessages(roc(response = obs, predictor = hot))),
      auc = as.numeric(roc(response = obs, predictor = hot)$auc),
      .groups = 'drop'
    )
  roc_data = map2_dfr(
    roc_results$roc_obj,
    roc_results$Resample,
    ~ data.frame(
      FPR = 1 - .x$specificities,
      Sensitivity = .x$sensitivities,
      Resample = .y,
      AUC = sprintf("AUC: %.3f", roc_results$auc[roc_results$Resample == .y]),
      stringsAsFactors = FALSE
    )
  )
  mean_roc = suppressMessages(roc(
    response = cv_predictions$obs,
    predictor = cv_predictions$hot
  ))
  roc_plot = ggplot() +
    geom_line(data = roc_data, 
              aes(x = FPR, y = Sensitivity, group = Resample),
              color = "grey70", 
              alpha = 0.4,
              linewidth = 0.6) +
    geom_line(data = data.frame(
      FPR = 1 - mean_roc$specificities,
      Sensitivity = mean_roc$sensitivities),
      aes(x = FPR, y = Sensitivity),
      color = "#009ad6", 
      linewidth = 1.5) +
    geom_abline(slope = 1, intercept = 0, 
                linetype = "dashed", 
                color = "grey40") +
    annotate("text", 
             x = 0.7, y = 0.25, 
             label = sprintf("Mean AUC = %.3f", 
                             mean_roc$auc),
             size = 5,
             fontface = "bold") +
    labs(x = "1-Specificity",
         y = "Sensitivity",
         title = "ROC curves",
         subtitle = paste("Using optimal λ =", round(best_lambda, 5))) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.minor = element_blank())
  print(roc_plot)
  cat("\nCross-validated AUC statistics:\n")
  print(summary(roc_results$auc))
  cat("SD:", sd(roc_results$auc), "\n")
} else {
  message("Tip: Please set savePredictions = TRUE in the trainControl first!")
}
coefs = coef(lasso_model$finalModel, s = lasso_model$bestTune$lambda)
selected_vars = rownames(coefs)[which(coefs != 0)][-1]
cat("=== Selected Variables (n =", length(selected_vars), ") ===\n")
print(selected_vars)
save(lasso_model, file = "lasso_model.Rdata")
# ------------------------------------------------------------------------------
# GSE19423
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(dplyr)
library(IOBR)
library(stringr)
setwd("input")
exp=read.delim("GSE19423_GPL6102_1_13_1_Matrix.matrix.txt")
rownames(exp)=exp[,1]
exp=exp[,-1]
exp = exp[rowSums(exp)>0,]
im_timer=deconvo_tme(eset = exp, method = "xcell")
im_timer=as.data.frame(im_timer)
rownames(im_timer)=im_timer[,1]
im_timer=im_timer[,-1]
colnames(im_timer) = gsub("_xCell$", "", colnames(im_timer))
im_timer = as.data.frame(scale(im_timer))
load("input/lasso_model.Rdata")
im_timer$subtype = predict(object = lasso_model, 
                            newdata = as.matrix(im_timer))
table(im_timer$subtype)
im_timer=im_timer[im_timer$subtype=="cold",]
exp=exp[,rownames(im_timer)]
library(GSEABase)
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=exp,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
data=read.delim("GSE19423_series_matrix.txt")[c(1,19,21),]
pro=as.data.frame(t(data))
pro=pro[-1,]
colnames(pro)=c("id","status","time")
pro$status=ifelse(pro$status=="overall survival: death",0,1)
pro$time=as.numeric(str_extract(pro$time, "(?<=: ).*$"))
WNT_score_pro=merge(pro,WNT_score,by.x="id",by.y="row.names")
WNT_score_pro=WNT_score_pro[,c("activity_score","status","time")]
WNT_score_pro=na.omit(WNT_score_pro)
quartiles=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)  
WNT_score_pro$level = cut(WNT_score_pro$activity_score, breaks = quartiles, include.lowest = TRUE,   
                          labels = c("Low","Medium","High"))  
WNT_score_pro$level=factor(WNT_score_pro$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(time, status)~ level,data= WNT_score_pro,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(time, status)~ level,data= WNT_score_pro,p.adjust.method = "none")
mycol=c('#009ad6',"#BBBDBE",'#d71345')
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
  data = WNT_score_pro,  
  xlim = c(0, 135),
  size = 0.5,  
  break.time.by = 12,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = "OS(n=42,GSE19423)",
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