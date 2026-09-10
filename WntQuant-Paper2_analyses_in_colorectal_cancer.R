# ==============================================================================
# WntQuant-Paper2 -- Colorectal cancer analyses
# ------------------------------------------------------------------------------
# Purpose:    Large-scale validation of the WntQuant-Paper score in colorectal cancer, including stratification by CMS molecular subtypes and association with molecular and clinical features.
# Inputs:     Processed CRC expression matrices with CMS subtype annotations.
# Outputs:    Wnt activity scores and subtype-level statistical results.
# NOTE:       Part of the code involves the core interests of the laboratory and is not open source at this time.
# ==============================================================================


rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("exp_data")
vector1=c("GSE44076","GSE44861","GSE21510",
          "GSE68468","GSE37178","GSE18105",
          "GSE21815","GSE17537","GSE29621",
          "GSE38832","GSE72970","GSE33113",
          "UCSC.CRC_CMS2")
vector2=c("GSE39582","GSE2109","GSE13067",
          "GSE13294","GSE14333","GSE17536",
          "GSE20916","GSE23878","GSE35896",
          "GSE37892","KFSYSCC","syn26720761",
          "syn2623706")
vector=vector2#Replacement of vector1 and vector2 is implemented here
for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
p_data_ee_all=data.frame()
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=1000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  WNT=scores$final_activity_score
  group=get(paste0(vector[i],"_CMS_G"))
  WNT_group=merge(WNT,group,by.x="row.names",by.y="Tag")
  #####################CMS1 vs. Other
  WNT_group_CMS1=WNT_group
  WNT_group_CMS1$group=ifelse(WNT_group_CMS1$group=="CMS1","CMS1","Other")
  WNT_group_CMS1$group=factor(WNT_group_CMS1$group,levels=c("CMS1","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS1_single=WNT_group_CMS1[,c(ee[j],"group")]
    if(ee[j]=="activity_score"){
      P.Value=wilcox.test(activity_score~group, WNT_group_CMS1_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="CMS1",]$activity_score)-
        mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="Other",]$activity_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=wilcox.test(pos_score~group, WNT_group_CMS1_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="CMS1",]$pos_score)-
        mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="Other",]$pos_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=wilcox.test(neg_score~group, WNT_group_CMS1_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="CMS1",]$neg_score)-
        mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="Other",]$neg_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="CMS1"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  #####################CMS2 vs. Other
  WNT_group_CMS2=WNT_group
  WNT_group_CMS2$group=ifelse(WNT_group_CMS2$group=="CMS2","CMS2","Other")
  WNT_group_CMS2$group=factor(WNT_group_CMS2$group,levels=c("CMS2","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS2_single=WNT_group_CMS2[,c(ee[j],"group")]
    if(ee[j]=="activity_score"){
      P.Value=wilcox.test(activity_score~group, WNT_group_CMS2_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="CMS2",]$activity_score)-
        mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="Other",]$activity_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=wilcox.test(pos_score~group, WNT_group_CMS2_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="CMS2",]$pos_score)-
        mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="Other",]$pos_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=wilcox.test(neg_score~group, WNT_group_CMS2_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="CMS2",]$neg_score)-
        mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="Other",]$neg_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="CMS2"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  #####################CMS3 vs. Other
  WNT_group_CMS3=WNT_group
  WNT_group_CMS3$group=ifelse(WNT_group_CMS3$group=="CMS3","CMS3","Other")
  WNT_group_CMS3$group=factor(WNT_group_CMS3$group,levels=c("CMS3","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS3_single=WNT_group_CMS3[,c(ee[j],"group")]
    if(ee[j]=="activity_score"){
      P.Value=wilcox.test(activity_score~group, WNT_group_CMS3_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="CMS3",]$activity_score)-
        mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="Other",]$activity_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=wilcox.test(pos_score~group, WNT_group_CMS3_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="CMS3",]$pos_score)-
        mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="Other",]$pos_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=wilcox.test(neg_score~group, WNT_group_CMS3_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="CMS3",]$neg_score)-
        mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="Other",]$neg_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="CMS3"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  #####################CMS4 vs. Other
  WNT_group_CMS4=WNT_group
  WNT_group_CMS4$group=ifelse(WNT_group_CMS4$group=="CMS4","CMS4","Other")
  WNT_group_CMS4$group=factor(WNT_group_CMS4$group,levels=c("CMS4","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS4_single=WNT_group_CMS4[,c(ee[j],"group")]
    if(ee[j]=="activity_score"){
      P.Value=wilcox.test(activity_score~group, WNT_group_CMS4_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="CMS4",]$activity_score)-
        mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="Other",]$activity_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=wilcox.test(pos_score~group, WNT_group_CMS4_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="CMS4",]$pos_score)-
        mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="Other",]$pos_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=wilcox.test(neg_score~group, WNT_group_CMS4_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="CMS4",]$neg_score)-
        mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="Other",]$neg_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="CMS4"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  #####################NOLBL vs. Other
  WNT_group_NOLBL=WNT_group
  WNT_group_NOLBL$group=ifelse(WNT_group_NOLBL$group=="NOLBL","NOLBL","Other")
  WNT_group_NOLBL$group=factor(WNT_group_NOLBL$group,levels=c("NOLBL","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_NOLBL_single=WNT_group_NOLBL[,c(ee[j],"group")]
    if(ee[j]=="activity_score"){
      P.Value=wilcox.test(activity_score~group, WNT_group_NOLBL_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="NOLBL",]$activity_score)-
        mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="Other",]$activity_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=wilcox.test(pos_score~group, WNT_group_NOLBL_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="NOLBL",]$pos_score)-
        mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="Other",]$pos_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=wilcox.test(neg_score~group, WNT_group_NOLBL_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="NOLBL",]$neg_score)-
        mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="Other",]$neg_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="NOLBL"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
}
geneSets=getGmt("wnt_usual_using.gmt")
p_data_ee_all_extra=data.frame()
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=gsva(expr=as.matrix(exp), geneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  scores=as.data.frame(t(scores))
  scores=as.data.frame(scale(scores))
  group=get(paste0(vector[i],"_CMS_G"))
  WNT_group=merge(scores,group,by.x="row.names",by.y="Tag")
  #####################CMS1 vs. Other
  WNT_group_CMS1=WNT_group
  WNT_group_CMS1$group=ifelse(WNT_group_CMS1$group=="CMS1","CMS1","Other")
  WNT_group_CMS1$group=factor(WNT_group_CMS1$group,levels=c("CMS1","Other"))
  ee=colnames(WNT_group_CMS1)[2:11]
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS1_single=WNT_group_CMS1[,c(ee[j],"group")]
    colnames(WNT_group_CMS1_single)=c("activity_score","group")
    P.Value=wilcox.test(activity_score~group, WNT_group_CMS1_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="CMS1",]$activity_score)-
      mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="Other",]$activity_score)
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS1"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  #####################CMS2 vs. Other
  WNT_group_CMS2=WNT_group
  WNT_group_CMS2$group=ifelse(WNT_group_CMS2$group=="CMS2","CMS2","Other")
  WNT_group_CMS2$group=factor(WNT_group_CMS2$group,levels=c("CMS2","Other"))
  ee=colnames(WNT_group_CMS2)[2:11]
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS2_single=WNT_group_CMS2[,c(ee[j],"group")]
    colnames(WNT_group_CMS2_single)=c("activity_score","group")
    P.Value=wilcox.test(activity_score~group, WNT_group_CMS2_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="CMS2",]$activity_score)-
      mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="Other",]$activity_score)
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS2"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  #####################CMS3 vs. Other
  WNT_group_CMS3=WNT_group
  WNT_group_CMS3$group=ifelse(WNT_group_CMS3$group=="CMS3","CMS3","Other")
  WNT_group_CMS3$group=factor(WNT_group_CMS3$group,levels=c("CMS3","Other"))
  ee=colnames(WNT_group_CMS3)[2:11]
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS3_single=WNT_group_CMS3[,c(ee[j],"group")]
    colnames(WNT_group_CMS3_single)=c("activity_score","group")
    P.Value=wilcox.test(activity_score~group, WNT_group_CMS3_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="CMS3",]$activity_score)-
      mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="Other",]$activity_score)
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS3"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  #####################CMS4 vs. Other
  WNT_group_CMS4=WNT_group
  WNT_group_CMS4$group=ifelse(WNT_group_CMS4$group=="CMS4","CMS4","Other")
  WNT_group_CMS4$group=factor(WNT_group_CMS4$group,levels=c("CMS4","Other"))
  ee=colnames(WNT_group_CMS4)[2:11]
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS4_single=WNT_group_CMS4[,c(ee[j],"group")]
    colnames(WNT_group_CMS4_single)=c("activity_score","group")
    P.Value=wilcox.test(activity_score~group, WNT_group_CMS4_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="CMS4",]$activity_score)-
      mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="Other",]$activity_score)
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS4"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  #####################NOLBL vs. Other
  WNT_group_NOLBL=WNT_group
  WNT_group_NOLBL$group=ifelse(WNT_group_NOLBL$group=="NOLBL","NOLBL","Other")
  WNT_group_NOLBL$group=factor(WNT_group_NOLBL$group,levels=c("NOLBL","Other"))
  ee=colnames(WNT_group_NOLBL)[2:11]
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_NOLBL_single=WNT_group_NOLBL[,c(ee[j],"group")]
    colnames(WNT_group_NOLBL_single)=c("activity_score","group")
    P.Value=wilcox.test(activity_score~group, WNT_group_NOLBL_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="NOLBL",]$activity_score)-
      mean(WNT_group_NOLBL_single[WNT_group_NOLBL_single$group=="Other",]$activity_score)
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="NOLBL"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
}

p_data_all=rbind(p_data_ee_all,p_data_ee_all_extra)
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","Wnt/β-catenin pathway activity score",
                          ifelse(p_data_all$geneset=="pos_score","Wnt/β-catenin pathway activation score",
                                 ifelse(p_data_all$geneset=="neg_score","Wnt/β-catenin pathway inhibition score",
                                        p_data_all$geneset)))
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "", 
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*", 
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))
p_data_all$geneset=factor(p_data_all$geneset,levels = c("Wnt/β-catenin pathway activity score",
                                                        "Wnt/β-catenin pathway activation score",
                                                        "Wnt/β-catenin pathway inhibition score",unique(p_data_ee_all_extra$geneset)))
p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
library(ggplot2)
library(dplyr)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)
ggplot(p_data_all, aes(dataset, geneset)) + 
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), col = "black", size = 6) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    text = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 16, color = "black")
  ) +
  labs(fill = paste0("Mean Difference"))
# ------------------------------------------------------------------------------
# WNT score calculation in syn2623706
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/crc_data_all.Rdata")
load("input/info_cms.Rdata")
load("input/info_icms.Rdata")
group=info_icms
co_sample=intersect(colnames(crc_data_all),group$Tag)
crc_data_all=crc_data_all[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(crc_data_all)
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
WNT_Score_group=merge(WNT_Score$final_activity_score,group,by.x="row.names",by.y="Tag",all=T)
table(WNT_Score_group$group)
if(group$group[1]%in%c("CMS1","CMS2","CMS3","CMS4","NOLBL")){
  WNT_Score_group$group=factor(WNT_Score_group$group,levels=c("CMS1","CMS2","CMS3","CMS4","NOLBL"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("CMS1","CMS2","CMS3","CMS4","NOLBL"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p = ggplot(WNT_Score_group, aes(x = group, y = activity_score, fill = group)) +
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) + 
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) + 
    xlab("") +
    ylab("Wnt/β-catenin pathway activity score") +
    guides(fill = "none") +
    scale_x_discrete(labels = c("CMS1\n(n=457)", 
                                "CMS2\n(n=1110)",
                                "CMS3\n(n=409)",
                                "CMS4\n(n=770)",
                                "NOLBL\n(n=486)")) +
    theme_bw() +
    theme(
      text = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black")
    ) +
    stat_compare_means(comparisons = my.comparisons, 
                       method = "wilcox.test",
                       label = "p.format",
                       size = 4) +
    scale_fill_manual(values = c('#f47920','#0073C2',"#ef5b9c",'#45b97c',"#a1a3a6")) +
    ggtitle("syn2623706")
  
  p
}else{
  WNT_Score_group$group=factor(WNT_Score_group$group,levels=c("iCMS2","iCMS3","indeterminate"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("iCMS2","iCMS3","indeterminate"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p = ggplot(WNT_Score_group, aes(x = group, y = activity_score, fill = group)) +
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) + 
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) + 
    xlab("") +
    ylab("Wnt/β-catenin pathway activity score") +
    guides(fill = "none") +
    scale_x_discrete(labels = c("iCMS2\n(n=1530)", 
                                "iCMS3\n(n=1319)",
                                "indeterminate\n(n=381)")) +
    theme_bw() +
    theme(
      text = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black")
    ) +
    stat_compare_means(comparisons = my.comparisons, 
                       method = "wilcox.test",
                       label = "p.format",
                       size = 4) +
    scale_fill_manual(values = c('#8552a1','#e0861a',"#a1a3a6")) +
    ggtitle("syn2623706")
  p
}
data = WNT_Score$final_activity_score[, 1, drop = F]  
colnames(data) = "p_score"  
p1 = ggplot(data, aes(x = p_score)) +  
  geom_density(fill = "#d71345", alpha = 0.5) +  
  geom_histogram(aes(y = ..density..), color = "transparent", fill = "transparent", bins = 30,  
                 breaks = seq(-4.3, 3.6, by = 1)) +  
  labs(title = "Distribution of Wnt/β-catenin pathway activity score",  
       x = "",  
       y = "Density") +  
  theme_minimal() +  
  theme(plot.title = element_text(size = 16, color = "black"),  
        axis.title = element_text(size = 16, color = "black"),  
        axis.text = element_text(size = 16, color = "black"),  
        legend.text = element_text(size = 16, color = "black"))

data = WNT_Score$final_activity_score[, 2, drop = F]  
colnames(data) = "p_score"  
p2 = ggplot(data, aes(x = p_score)) +  
  geom_density(fill = "#f58f98", alpha = 0.5) +  
  geom_histogram(aes(y = ..density..), color = "transparent", fill = "transparent", bins = 30,  
                 breaks = seq(-2.3, 1.9, by = 1)) +  
  labs(title = "Distribution of Wnt/β-catenin pathway activation score",  
       x = "",  
       y = "Density") +  
  theme_minimal() +  
  theme(plot.title = element_text(size = 16, color = "black"),  
        axis.title = element_text(size = 16, color = "black"),  
        axis.text = element_text(size = 16, color = "black"),  
        legend.text = element_text(size = 16, color = "black"))
data = WNT_Score$final_activity_score[, 3, drop = F]  
colnames(data) = "p_score"  
p3 = ggplot(data, aes(x = p_score)) +  
  geom_density(fill = "#90d7ec", alpha = 0.5) +  
  geom_histogram(aes(y = ..density..), color = "transparent", fill = "transparent", bins = 30,  
                 breaks = seq(-2.3, 2.6, by = 1)) +  
  labs(title = "Distribution of Wnt/β-catenin pathway inhibition score",  
       x = "",  
       y = "Density") +  
  theme_minimal() +  
  theme(plot.title = element_text(size = 16, color = "black"),  
        axis.title = element_text(size = 16, color = "black"),  
        axis.text = element_text(size = 16, color = "black"),  
        legend.text = element_text(size = 16, color = "black"))
combined_plot = p1 / p2 / p3  
print(combined_plot)
# ------------------------------------------------------------------------------
# CMS/iCMS combined calculations
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/crc_data_all.Rdata")
load("input/info_cms.Rdata")
load("input/info_icms.Rdata")
info_cms_icms=merge(info_cms,info_icms,by="Tag")
write.csv(info_cms_icms,file = "info_cms_icms.csv")
info_cms_icms$group=paste0(info_cms_icms$group.x,"_",info_cms_icms$group.y)
info_cms_icms=info_cms_icms[,c("Tag","group")]
group=info_cms_icms
co_sample=intersect(colnames(crc_data_all),group$Tag)
crc_data_all=crc_data_all[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(crc_data_all)
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
WNT_Score_group=merge(WNT_Score$final_activity_score,group,by.x="row.names",by.y="Tag",all=T)
table(WNT_Score_group$group)
fs=unique(WNT_Score_group$group)
p_data_all=data.frame()
for(i in 1:length(fs)){
  WNT_Score_group_fsi=WNT_Score_group[WNT_Score_group$group==fs[i],]
  for(j in 1:length(fs)){
    WNT_Score_group_fsj=WNT_Score_group[WNT_Score_group$group==fs[j],]
    P.Value=wilcox.test(WNT_Score_group_fsi$activity_score, WNT_Score_group_fsj$activity_score, alternative = "two.sided")$p.value
    FC=mean(WNT_Score_group_fsi$activity_score)-mean(WNT_Score_group_fsj$activity_score)
    p_data=data.frame(dataset1=fs[i],dataset2=fs[j],p=P.Value,FC=FC)
    p_data_all=rbind(p_data_all,p_data)
  }
}
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "", 
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*", 
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))
p_data_all$dataset1 = paste0(p_data_all$dataset1,"(","G1",")")  
p_data_all$dataset2 = paste0(p_data_all$dataset2,"(","G2",")")  
library(ggplot2)  
library(dplyr)  
library(RColorBrewer)  
p1 = ggplot(p_data_all, aes(dataset1, dataset2)) +   
  geom_tile(aes(fill = FC), color = "black") +  
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +  
  geom_text(aes(label = pstar), color = "black", size = 3) +  
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16), 
        axis.text.y = element_text(size = 16, color = "black"),  
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black")) +  
  labs(fill = paste0("Mean Difference"))
WNT_Score_group$group = factor(WNT_Score_group$group, levels = c("CMS1_iCMS2", "CMS1_iCMS3", "CMS1_indeterminate",  
                                                                 "CMS2_iCMS2", "CMS2_iCMS3", "CMS2_indeterminate",  
                                                                 "CMS3_iCMS2", "CMS3_iCMS3", "CMS3_indeterminate",  
                                                                 "CMS4_iCMS2", "CMS4_iCMS3", "CMS4_indeterminate",  
                                                                 "NOLBL_iCMS2", "NOLBL_iCMS3", "NOLBL_indeterminate"))
median_scores = aggregate(activity_score ~ group, data = WNT_Score_group, FUN = median)  
sorted_sites = median_scores$group[order(median_scores$activity_score, decreasing = TRUE)]  
colors = colorRampPalette(c("#f58f98", "white", "#90d7ec"))(length(sorted_sites))  
WNT_Score_group$group = factor(WNT_Score_group$group, levels = sorted_sites)
p2 = ggplot(WNT_Score_group, aes(x = group, y = activity_score, fill = group)) +  
  geom_boxplot(color = "black") +  
  theme_bw() +  
  theme(axis.text.x = element_text(angle = 80, hjust = 1, color = "black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        legend.position = 'none') +  
  scale_fill_manual(values = colors) +  
  labs(y = "Wnt/β-catenin pathway activity score", x = "") + 
  ggtitle("syn2623706") +  
  theme(plot.title = element_text(size = 16, color = "black"))
library(patchwork)  
combined_plot = p2 + p1  
combined_plot
# ------------------------------------------------------------------------------
# Comparison with iCMS2/iCMS3 markers
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("exp_data")
vector1=c("GSE44076","GSE44861","GSE21510",
          "GSE68468","GSE37178","GSE18105",
          "GSE21815","GSE17537","GSE29621",
          "GSE38832","GSE72970","GSE33113",
          "UCSC.CRC_CMS2")
vector2=c("GSE39582","GSE2109","GSE13067",
          "GSE13294","GSE14333","GSE17536",
          "GSE20916","GSE23878","GSE35896",
          "GSE37892","KFSYSCC","syn26720761",
          "syn2623706")
vector=vector2
for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
icms2_icms3=getGmt("iCMS2_iCMS3.gmt")
p_data_ee_all=data.frame()
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=1000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  WNT=scores$final_activity_score
  icms=gsva(expr=as.matrix(exp), icms2_icms3, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  icms=as.data.frame(t(icms))
  icms=as.data.frame(scale(icms))
  icms_wnt=merge(WNT,icms,by="row.names")
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  ##############################iCMS2_UP
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Up")]
    if(ee[j]=="activity_score"){
      P.Value=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="iCMS2_Up"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  p_data_ee=data.frame()
  ##############################iCMS2_DN
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Down")]
    if(ee[j]=="activity_score"){
      P.Value=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="iCMS2_Down"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  p_data_ee=data.frame()
  ##############################iCMS3_UP
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Up")]
    if(ee[j]=="activity_score"){
      P.Value=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="iCMS3_Up"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
  p_data_ee=data.frame()
  ##############################iCMS3_DN
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Down")]
    if(ee[j]=="activity_score"){
      P.Value=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$activity_score, WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="pos_score"){
      P.Value=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$pos_score, WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }else if(ee[j]=="neg_score"){
      P.Value=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$p.value
      Correlation=cor.test(WNT_group_single$neg_score, WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$estimate
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
      p_data_ee=rbind(p_data_ee,p_data)
    }
  }
  p_data_ee$contrast="iCMS3_Down"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
}
p_data_ee_all_extra=data.frame()
geneSets=getGmt("wnt_usual_using.gmt")
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  exp=get(paste0(vector[i],"_CMS"))
  scores=gsva(expr=as.matrix(exp), geneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  scores=as.data.frame(t(scores))
  WNT=as.data.frame(scale(scores))
  icms=gsva(expr=as.matrix(exp), icms2_icms3, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  icms=as.data.frame(t(icms))
  icms=as.data.frame(scale(icms))
  icms_wnt=merge(WNT,icms,by="row.names")
  ee=colnames(icms_wnt)[2:11]
  p_data_ee=data.frame()
  ##############################iCMS2_UP
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Up")]
    P.Value=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$p.value
    Correlation=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Up,alternative="two.sided",method="pearson")$estimate
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS2_Up"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  p_data_ee=data.frame()
  ##############################iCMS2_DN
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Down")]
    P.Value=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$p.value
    Correlation=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Down,alternative="two.sided",method="pearson")$estimate
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS2_Down"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  p_data_ee=data.frame()
  ##############################iCMS3_UP
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Up")]
    P.Value=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$p.value
    Correlation=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Up,alternative="two.sided",method="pearson")$estimate
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS3_Up"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  p_data_ee=data.frame()
  ##############################iCMS3_DN
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Down")]
    P.Value=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$p.value
    Correlation=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Down,alternative="two.sided",method="pearson")$estimate
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,Correlation=Correlation)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS3_Down"
  p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
}
p_data_all=rbind(p_data_ee_all,p_data_ee_all_extra)
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","Wnt pathway activity score",
                          ifelse(p_data_all$geneset=="pos_score","Wnt pathway activation score",
                                 ifelse(p_data_all$geneset=="neg_score","Wnt pathway inhibition score",
                                        p_data_all$geneset)))
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "", 
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*", 
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))
p_data_all$geneset=factor(p_data_all$geneset,levels = c("Wnt pathway activity score",
                                                        "Wnt pathway activation score",
                                                        "Wnt pathway inhibition score",unique(p_data_ee_all_extra$geneset)))
p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
library(ggplot2)
library(dplyr)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)
ggplot(p_data_all, aes(dataset, geneset)) + 
  geom_tile(aes(fill = Correlation), color = "black") +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black")
  ) +
  labs(fill = paste0("Correlation"))
# ------------------------------------------------------------------------------
# Comparison of expression levels of characterised genes in syn2623706
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(reshape2)
setwd("input")
load("input/crc_data_all.Rdata")
load("input/info_cms.Rdata")
load("input/info_icms.Rdata")
group=info_cms#Substitution can be done here, the replacement variables are info_cms and info_icms.
co_sample=intersect(colnames(crc_data_all),group$Tag)
crc_data_all=crc_data_all[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(crc_data_all)
exp=as.data.frame(t(exp[intersect(c("AXIN2","ZNRF3","RNF43","LGR5","TCF7"),rownames(exp)),]))
biomarker_group=merge(exp,group,by.x="row.names",by.y="Tag",all=T)
biomarker_group=biomarker_group[,-1]
biomarker_group=melt(biomarker_group)
colnames(biomarker_group)=c("group","gene","value")
table(biomarker_group$group)
if(group$group[1]%in%c("CMS1","CMS2","CMS3","CMS4","NOLBL")){
  biomarker_group$group=factor(biomarker_group$group,levels=c("CMS1","CMS2","CMS3","CMS4","NOLBL"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  p = ggplot(biomarker_group, aes(gene, value, fill = group, color = group)) +  
    geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
    labs(x = " ", y = "Expression") +  
    theme_grey() +  
    theme(  
      panel.grid.major.x = element_line(linewidth = 1),  
      panel.grid.major.y = element_line(linewidth = 1),  
      panel.grid.minor.x = element_line(linewidth = 0.8),  
      panel.grid.minor.y = element_line(linewidth = 0.8),  
      axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
      axis.line.y = element_line(colour = 'black', linewidth = 1),  
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 16),  
      axis.text.y = element_text(color = "black", size = 16), 
      axis.title.y = element_text(angle = 90, size = 16, color = "black"),  
      axis.title.x = element_text(size = 16, color = "black"),
      legend.text = element_text(color = "black", size = 16), 
      legend.title = element_text(color = "black", size = 16)
    ) +  
    scale_fill_manual(values = c('#f47920', '#0073C2', "#ef5b9c", '#45b97c', "#a1a3a6")) +  
    stat_compare_means(aes(group = group, label = ..p.signif..)) +  
    ggtitle("syn2623706") 
  p
}else{
  biomarker_group$group=factor(biomarker_group$group,levels=c("iCMS2","iCMS3","indeterminate"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  p = ggplot(biomarker_group, aes(gene, value, fill = group, color = group)) +  
    geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
    labs(x = " ", y = "Expression") +  
    theme_grey() +  
    theme(  
      panel.grid.major.x = element_line(linewidth = 1),  
      panel.grid.major.y = element_line(linewidth = 1),  
      panel.grid.minor.x = element_line(linewidth = 0.8),  
      panel.grid.minor.y = element_line(linewidth = 0.8),  
      axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
      axis.line.y = element_line(colour = 'black', linewidth = 1),  
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16), 
      axis.text.y = element_text(colour = "black", size = 16), 
      axis.title.x = element_text(size = 16, colour = "black"),
      axis.title.y = element_text(angle = 90, size = 16, colour = "black"), 
      legend.text = element_text(colour = "black", size = 16),  
      legend.title = element_text(colour = "black", size = 16) 
    ) +  
    scale_fill_manual(values = c('#8552a1', '#e0861a', "#a1a3a6")) +  
    stat_compare_means(aes(group = group, label = ..p.signif..)) +  
    ggtitle("syn2623706") 
  p
}
# ------------------------------------------------------------------------------
# Comparison of original signatures in syn2623706
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("input")
load("input/crc_data_all.Rdata")
load("input/info_cms.Rdata")
load("input/info_icms.Rdata")
group=info_cms#Substitution can be done here, the replacement variables are info_cms and info_icms.
co_sample=intersect(colnames(crc_data_all),group$Tag)
crc_data_all=crc_data_all[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(crc_data_all)
GeneSets=getGmt("NM_genesets.gmt")
geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
geneset=as.data.frame(t(geneset))
geneset=as.data.frame(scale(geneset))
geneset_wnt=geneset[,c("WNT_FLIER"),drop=F]
WNT=merge(geneset_wnt,group,by.x="row.names",by.y="Tag")
table(WNT$group)
if(group$group[1]%in%c("CMS1","CMS2","CMS3","CMS4","NOLBL")){
  WNT$group=factor(WNT$group,levels=c("CMS1","CMS2","CMS3","CMS4","NOLBL"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("CMS1","CMS2","CMS3","CMS4","NOLBL"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p <- ggplot(WNT, aes(x = group, y = WNT_FLIER, fill = group)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("ssGSEA score of FLIER_Wnt") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c("CMS1\n(n=457)",   
                                "CMS2\n(n=1110)",  
                                "CMS3\n(n=409)",  
                                "CMS4\n(n=770)",  
                                "NOLBL\n(n=486)")) +  
    theme_bw() +  
    theme(  
      text = element_text(color = "black", size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16)
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c('#f47920', '#0073C2', "#ef5b9c", '#45b97c', "#a1a3a6")) +  
    ggtitle("syn2623706")  
  p  
}else{
  WNT$group=factor(WNT$group,levels=c("iCMS2","iCMS3","indeterminate"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("iCMS2","iCMS3","indeterminate"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p <- ggplot(WNT, aes(x = group, y = WNT_FLIER, fill = group)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("ssGSEA score of FLIER_Wnt") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c("iCMS2\n(n=1530)",   
                                "iCMS3\n(n=1319)",  
                                "indeterminate\n(n=381)")) +  
    theme_bw() +  
    theme(  
      text = element_text(color = "black", size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16)
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c('#8552a1', '#e0861a', "#a1a3a6")) +  
    ggtitle("syn2623706")  
  p  
}
# ------------------------------------------------------------------------------
# Number of clusters for syn2623706 (K-means typing)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/crc_data_all.Rdata")
exp=as.data.frame(crc_data_all)
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
save(WNT_Score,file="WNT_Score.Rdata")
genesets=as.matrix(WNT_Score$single_score_matrix)
library(mlbench)
library(factoextra)
library(ggplot2)
library (cluster)
library(fpc)
set.seed(123)
res=get_clust_tendency(genesets,300)
res$hopkins_stat
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
km_result$cluster=ifelse(km_result$cluster==3,1,ifelse(km_result$cluster==2,4,ifelse(km_result$cluster==1,2,3)))
dd=cbind(genesets, cluster = km_result$cluster)
fviz_cluster(km_result, data = genesets,  
             palette = c("#d71345", "#45b97c", "#145b7d", "#ffd400"),  
             ellipse.type = "euclid",  
             star.plot = TRUE,   
             repel = TRUE,  
             geom = "point") +   
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
cluster$cluster=ifelse(cluster$cluster=="cluster1","WSS1",
                       ifelse(cluster$cluster=="cluster2","WSS2",
                              ifelse(cluster$cluster=="cluster3","WSS3","WSS4")))
WNT=merge(WNT_Score$final_activity_score,cluster,by.x="row.names",by.y="Tag")
WNT$project="CMS"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("WSS1","WSS2","WSS3","WSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p <- ggplot(WNT, aes(x = cluster, y = activity_score, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("WSS1\n(n=931)",   
                              "WSS2\n(n=865)",  
                              "WSS3\n(n=905)",  
                              "WSS4\n(n=531)")) +  
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
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  ggtitle("syn2623706")
print(p)
save(cluster,file="cluster.Rdata")
# ------------------------------------------------------------------------------
# Percentage of each subtype in different activities of WNT and prognostic analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("input/cluster.Rdata")
load("input/info_cms.Rdata")
load("input/info_icms.Rdata")
cluster_group=merge(cluster,info_icms,by="Tag")
cluster_group$level=cluster_group$cluster
cluster1=cluster_group[cluster_group$level=="WSS1",]
cluster2=cluster_group[cluster_group$level=="WSS2",]
cluster3=cluster_group[cluster_group$level=="WSS3",]
cluster4=cluster_group[cluster_group$level=="WSS4",]
cluster1=cluster1 %>% group_by(group) %>% tally() %>% 
  mutate(prop = round((n/sum(n)) * 100 , digits = 2))
cluster2=cluster2 %>% group_by(group) %>% tally() %>% 
  mutate(prop = round((n/sum(n)) * 100 , digits = 2))
cluster3=cluster3 %>% group_by(group) %>% tally() %>% 
  mutate(prop = round((n/sum(n)) * 100 , digits = 2))
cluster4=cluster4 %>% group_by(group) %>% tally() %>% 
  mutate(prop = round((n/sum(n)) * 100 , digits = 2))
cluster1$cut="WSS1"
cluster2$cut="WSS2"
cluster3$cut="WSS3"
cluster4$cut="WSS4"
HML=rbind(cluster1,cluster2)
HML=rbind(HML,cluster3)
HML=rbind(HML,cluster4)
colnames(HML)[1]="Subtype"
HML$cut=factor(HML$cut,level=c("WSS1","WSS2","WSS3","WSS4"))
if(cluster_group$group[1]%in%c("CMS1","CMS2","CMS3","CMS4","NOLBL")){
  mycol=c('#f47920','#0073C2',"#ef5b9c",'#45b97c',"#a1a3a6")
}else{
  mycol=c('#8552a1','#e0861a',"#a1a3a6")
}
ggplot(HML, aes(x = 3,
                y = prop,
                fill = Subtype)) +
  geom_col(width = 1.5,
           color = 'black', alpha = 0.8) +
  facet_grid(. ~ cut) +
  coord_polar(theta = "y") +
  xlim(c(0.2, 3.8)) +
  scale_fill_manual(values = mycol) +
  theme_bw() +
  theme(  
    strip.text.x = element_text(size = 16, color = 'black'),
    legend.title = element_text(size = 16, color = 'black'),
    legend.text = element_text(size = 16, color = 'black'),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.title = element_text(size = 16, color = 'black')
  ) +   
  labs(fill = "Subtype", y = "", x = "")  
# ------------------------------------------------------------------------------
# OS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("input/cluster.Rdata")
setwd("input")
clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c("sample","osMo","osStat")]

TCGACRC_clinical_merged=read.table(file = "TCGACRC_clinical-merged.tsv",header = T,sep = "\t",check.names = F)
TCGACRC_clinical_merged=TCGACRC_clinical_merged[,c(1,11,12)]
colnames(TCGACRC_clinical_merged)=c("sample","osMo","osStat")

gse17536_survival=read.table(file = "gse17536_survival.tsv",header = T,sep = "\t",check.names = F)
gse17536_survival=gse17536_survival[,c(1,3,2)]
colnames(gse17536_survival)=c("sample","osMo","osStat")

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[c(1,22,23),]))
GSE39582_series_matrix=GSE39582_series_matrix[-1,]
GSE39582_series_matrix=GSE39582_series_matrix[,c(1,3,2)]
colnames(GSE39582_series_matrix)=c("sample","osMo","osStat")
GSE39582_series_matrix$osMo=sub(".*: ", "",GSE39582_series_matrix$osMo)
GSE39582_series_matrix$osStat=sub(".*: ", "",GSE39582_series_matrix$osStat)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$osMo!="N/A",]
GSE39582_series_matrix$osMo=as.numeric(GSE39582_series_matrix$osMo)
GSE39582_series_matrix$osStat=as.numeric(GSE39582_series_matrix$osStat)

all_data=rbind(clinical_molecular_public_all,TCGACRC_clinical_merged)
all_data=rbind(all_data,gse17536_survival)
all_data=rbind(all_data,GSE39582_series_matrix)
unique_ids=unique(all_data$sample)
all_data=all_data[match(unique_ids, all_data$sample), ]
OS_data=all_data
save(OS_data,file="OS_data.Rdata")
survival_data=OS_data
survival_data=na.omit(survival_data)
survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(survival)
library(survminer)
fitd = survdiff(Surv(osMo, osStat) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(osMo, osStat)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(osMo, osStat)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#45b97c","#145b7d","#ffd400")
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
  xlim = c(0, 204),
  size = 0.5,  
  break.time.by = 12,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = "OS(n=1230)",
  ylim = c(0.2, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0.2, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) +
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# RFS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("input/cluster.Rdata")
setwd("input")
clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c("sample","rfsMo","rfsStat")]

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,10,11)]
colnames(FRENCH_clinical)=c("sample","rfsMo","rfsStat")

gse14333_dfs_times=read.table(file = "gse14333_dfs_times.tsv",header = T,sep = "\t",check.names = F)
gse14333_dfs_times$sample=rownames(gse14333_dfs_times)
gse14333_dfs_times$sample=gsub("\\.CEL$", "",gse14333_dfs_times$sample)
gse14333_dfs_times=gse14333_dfs_times[,c(3,1,2)]
colnames(gse14333_dfs_times)=c("sample","rfsMo","rfsStat")

gse17536_survival=read.table(file = "gse17536_survival.tsv",header = T,sep = "\t",check.names = F)
gse17536_survival=gse17536_survival[,c(1,7,6)]
colnames(gse17536_survival)=c("sample","rfsMo","rfsStat")

gse37892_dfs_times=read.table(file = "gse37892_dfs_times.tsv",header = T,sep = "\t",check.names = F)
gse37892_dfs_times$sample=rownames(gse37892_dfs_times)
gse37892_dfs_times=gse37892_dfs_times[,c(3,2,1)]
colnames(gse37892_dfs_times)=c("sample","rfsMo","rfsStat")

GSE33113_series_matrix=read.delim(file = "GSE33113_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE33113_series_matrix=GSE33113_series_matrix[1:34,]
GSE33113_series_matrix=GSE33113_series_matrix[c(1,13,14),]
GSE33113_series_matrix=as.data.frame(t(GSE33113_series_matrix[,-1]))
GSE33113_series_matrix=GSE33113_series_matrix[GSE33113_series_matrix$'13'!="",]
GSE33113_series_matrix=GSE33113_series_matrix[,c(1,3,2)]
colnames(GSE33113_series_matrix)=c("sample","rfsMo","rfsStat")
GSE33113_series_matrix$rfsMo=sub(".*: ", "",GSE33113_series_matrix$rfsMo)
GSE33113_series_matrix$rfsStat=sub(".*: ", "",GSE33113_series_matrix$rfsStat)
GSE33113_series_matrix$rfsStat=ifelse(GSE33113_series_matrix$rfsStat=="yes",1,0)
GSE33113_series_matrix$rfsMo=as.numeric(GSE33113_series_matrix$rfsMo)/30

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[c(1,20,21),]))
GSE39582_series_matrix=GSE39582_series_matrix[-1,]
GSE39582_series_matrix=GSE39582_series_matrix[,c(1,3,2)]
colnames(GSE39582_series_matrix)=c("sample","rfsMo","rfsStat")
GSE39582_series_matrix$rfsMo=sub(".*: ", "",GSE39582_series_matrix$rfsMo)
GSE39582_series_matrix$rfsStat=sub(".*: ", "",GSE39582_series_matrix$rfsStat)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$rfsMo!="N/A",]
GSE39582_series_matrix$rfsMo=as.numeric(GSE39582_series_matrix$rfsMo)
GSE39582_series_matrix$rfsStat=as.numeric(GSE39582_series_matrix$rfsStat)

all_data=rbind(clinical_molecular_public_all,FRENCH_clinical)
all_data=rbind(all_data,gse14333_dfs_times)
all_data=rbind(all_data,gse17536_survival)
all_data=rbind(all_data,gse37892_dfs_times)
all_data=rbind(all_data,GSE33113_series_matrix)
all_data=rbind(all_data,GSE39582_series_matrix)
unique_ids=unique(all_data$sample)
all_data=all_data[match(unique_ids, all_data$sample), ]
RFS_data=all_data
save(RFS_data,file="RFS_data.Rdata")
survival_data=RFS_data
survival_data=na.omit(survival_data)
survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(survival)
library(survminer)
fitd = survdiff(Surv(rfsMo, rfsStat) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(rfsMo, rfsStat)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(rfsMo, rfsStat)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#45b97c","#145b7d","#ffd400")
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
  xlim = c(0, 204),
  size = 0.5,  
  break.time.by = 12,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = "RFS(n=912)",
  ylim = c(0.38, 1)
)  
p.lab = paste0("log-rank test p",  
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))  
p$plot = p$plot +   
  annotate("text", x = 0, y = 0.4, hjust = 0, vjust = 0,   
           fontface = "italic", label = p.lab, size =6) +
  theme(text = element_text(size = 16))
print(p)
# ------------------------------------------------------------------------------
# SAR
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("input/cluster.Rdata")
load("input/OS_data.Rdata")
load("input/RFS_data.Rdata")
survial=merge(OS_data,RFS_data,by="sample")
survial=na.omit(survial)
survial=survial[survial$rfsStat==1,]
survial=survial[,-c(4,5)]
survival_data=survial
survival_data=na.omit(survival_data)
survival_data_cms=merge(survival_data,cluster,by.x="sample",by.y="Tag")
survival_data_cms$level=survival_data_cms$cluster
survival_data_cms$level=factor(survival_data_cms$level,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(survival)
library(survminer)
fitd = survdiff(Surv(osMo, osStat) ~ level,data= survival_data_cms,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(osMo, osStat)~ level,data= survival_data_cms,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(osMo, osStat)~ level,data= survival_data_cms,p.adjust.method = "none")
mycol=c("#d71345","#45b97c","#145b7d","#ffd400")
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
  xlim = c(0, 180),
  size = 0.5,  
  break.time.by = 12,
  legend.title = "",
  xlab = "Time (months)",
  ylab = "Probability",
  risk.table.y.text = FALSE,
  tables.height = 0.25,
  title = "SAR(n=202)",
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
# Expression of marker genes in new typing
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(reshape2)
setwd("input")
load("input/crc_data_all.Rdata")
load("input/cluster.Rdata")
group=cluster
colnames(group)=c("Tag","group")
co_sample=intersect(colnames(crc_data_all),group$Tag)
crc_data_all=crc_data_all[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(crc_data_all)
exp=as.data.frame(t(exp[intersect(c("AXIN2","ZNRF3","RNF43","LGR5","TCF7"),rownames(exp)),]))
biomarker_group=merge(exp,group,by.x="row.names",by.y="Tag",all=T)
biomarker_group=biomarker_group[,-1]
biomarker_group=melt(biomarker_group)
colnames(biomarker_group)=c("group","gene","value")
table(biomarker_group$group)
biomarker_group$group=factor(biomarker_group$group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(biomarker_group, aes(gene, value, fill = group, color = group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Expression") +  
  theme_grey() +  
  theme(  
    panel.grid.major.x = element_line(linewidth = 1),  
    panel.grid.major.y = element_line(linewidth = 1),  
    panel.grid.minor.x = element_line(linewidth = 0.8),  
    panel.grid.minor.y = element_line(linewidth = 0.8),  
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
    axis.line.y = element_line(colour = 'black', linewidth = 1),  
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16),  
    axis.text.y = element_text(color = "black", size = 16),   
    axis.title.y = element_text(angle = 90, size = 16),   
    axis.title.x = element_text(size = 16),   
    plot.title = element_text(size = 16, colour = "black")   
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = group, label = ..p.signif..), size = 5, color = "black") +
  ggtitle("syn2623706")   
p  
# ------------------------------------------------------------------------------
# MSI-H/MSS/MSI-L
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c(1,11)]
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$msi=ifelse(clinical_molecular_public_all$msi=="msi","MSI-H","MSI-L/MSS")

TCGACRC_clinical_merged=read.table(file = "TCGACRC_clinical-merged.tsv",header = T,sep = "\t",check.names = F)
TCGACRC_clinical_merged=TCGACRC_clinical_merged[,c(1,14)]
TCGACRC_clinical_merged=TCGACRC_clinical_merged[TCGACRC_clinical_merged$microsatelite%in%c("MSS","MSI-H","MSI-L"),]
colnames(TCGACRC_clinical_merged)=c("sample","msi")
TCGACRC_clinical_merged$msi=ifelse(TCGACRC_clinical_merged$msi=="MSI-H","MSI-H","MSI-L/MSS")

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,13)]
FRENCH_clinical=na.omit(FRENCH_clinical)
colnames(FRENCH_clinical)=c("sample","msi")
FRENCH_clinical$msi=ifelse(FRENCH_clinical$msi=="pMMR","MSI-L/MSS","MSI-H")

gse13067_annotation=read.table(file = "gse13067_annotation.tsv",header = T,sep = "\t",check.names = F)
colnames(gse13067_annotation)=c("sample","msi")
gse13067_annotation$msi=ifelse(gse13067_annotation$msi=="msi","MSI-H","MSI-L/MSS")

gse13294_annotation=read.table(file = "gse13294_annotation.tsv",header = T,sep = "\t",check.names = F)
colnames(gse13294_annotation)=c("sample","msi")
gse13294_annotation$msi=ifelse(gse13294_annotation$msi=="msi","MSI-H","MSI-L/MSS")

gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,12)]
colnames(gse35896_annotation)=c("sample","msi")
gse35896_annotation$msi=sub(".*: ", "",gse35896_annotation$msi)
gse35896_annotation=gse35896_annotation[gse35896_annotation$msi%in%c("MSI","MSS"),]
gse35896_annotation$msi=ifelse(gse35896_annotation$msi=="MSI","MSI-H","MSI-L/MSS")

auxiliary_coad=read.table(file = "auxiliary_coad-20140108.txt",header = T,sep = "\t",check.names = F)
auxiliary_coad=auxiliary_coad[,c(1,3)]
colnames(auxiliary_coad)=c("sample","msi")
auxiliary_coad$msi=ifelse(auxiliary_coad$msi=="MSI-H","MSI-H","MSI-L/MSS")

auxiliary_read=read.table(file = "auxiliary_read-20140108.txt",header = T,sep = "\t",check.names = F)
auxiliary_read=auxiliary_read[,c(1,3)]
colnames(auxiliary_read)=c("sample","msi")
auxiliary_read$msi=ifelse(auxiliary_read$msi=="MSI-H","MSI-H","MSI-L/MSS")

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,24),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","msi")
GSE39582_series_matrix$msi=sub(".*: ", "",GSE39582_series_matrix$msi)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$msi%in%c("pMMR","dMMR"),]
GSE39582_series_matrix$msi=ifelse(GSE39582_series_matrix$msi=="pMMR","MSI-L/MSS","MSI-H")

clinic=read.csv(file = "SG-BULK_patient_clinical_information.csv")
clinic$patient_id=paste0("X",clinic$patient_id)
clinic=clinic[,c(1,17)]
clinic=na.omit(clinic)
colnames(clinic)=c("sample","msi")
clinic$msi=ifelse(clinic$msi=="MSI","MSI-H","MSI-L/MSS")

MSI_all=rbind(clinical_molecular_public_all,TCGACRC_clinical_merged)
MSI_all=rbind(MSI_all,FRENCH_clinical)
MSI_all=rbind(MSI_all,gse13067_annotation)
MSI_all=rbind(MSI_all,gse13294_annotation)
MSI_all=rbind(MSI_all,gse35896_annotation)
MSI_all=rbind(MSI_all,auxiliary_coad)
MSI_all=rbind(MSI_all,auxiliary_read)
MSI_all=rbind(MSI_all,GSE39582_series_matrix)
MSI_all=rbind(MSI_all,clinic)
table(MSI_all$msi)
unique_ids=unique(MSI_all$sample)
MSI_all=MSI_all[match(unique_ids, MSI_all$sample), ]

MSI=MSI_all
MSI=merge(cluster,MSI,by.x="Tag",by.y="sample")
MSI=MSI[,c("cluster","msi")]
MSI=na.omit(MSI)
MSI$status="Microsatellite_status"
WNT_cluster1=MSI[MSI$cluster=="WSS1",]
WNT_cluster2=MSI[MSI$cluster=="WSS2",]
WNT_cluster3=MSI[MSI$cluster=="WSS3",]
WNT_cluster4=MSI[MSI$cluster=="WSS4",]
MSI=MSI[,-3]
MSI=as.data.frame(table(MSI))
contingency_table=acast(MSI, cluster ~ msi, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1=ggplot(WNT_cluster1, aes(x = "", fill = msi)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MSI-H" = "#d71345", "MSI-L/MSS" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "syn2623706(n=1345,Microsatellite Status)", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",  
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )

p2=ggplot(WNT_cluster2, aes(x = status, fill = msi)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MSI-H" = "#d71345", "MSI-L/MSS" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    legend.position = "none",  
    text = element_text(color = "black", size = 16)  
  )  

p3=ggplot(WNT_cluster3, aes(x = status, fill = msi)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MSI-H" = "#d71345", "MSI-L/MSS" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    legend.position = "none",  
    text = element_text(color = "black", size = 16)  
  )  
p4=ggplot(WNT_cluster4, aes(x = status, fill = msi)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MSI-H" = "#d71345", "MSI-L/MSS" = "#145b7d")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    legend.position = "right",
    text = element_text(color = "black", size = 16)
  )  
final_plot=p1 | p2 | p3 | p4  
final_plot 
# ------------------------------------------------------------------------------
# MUTATION_COUNT
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
load("input/mutation_count.Rdata")
pancancer_df$sampleId=substr(pancancer_df$sampleId,1,12)
WNT=merge(pancancer_df,cluster,by.x="sampleId",by.y="Tag")
WNT$MUTATION_COUNT=as.numeric(WNT$MUTATION_COUNT)
WNT$project="CMS"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("WSS1","WSS2","WSS3","WSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p <- ggplot(WNT, aes(x = cluster, y = log2(MUTATION_COUNT), fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab(expression("Mutation count (log"[2] * "-transformation)")) +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("WSS1\n(n=69)",   
                              "WSS2\n(n=61)",  
                              "WSS3\n(n=52)",  
                              "WSS4\n(n=28)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  ggtitle("syn2623706") 
p
# ------------------------------------------------------------------------------
# TMB
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("input/crc_data_all.Rdata")
load("input/info_cms.Rdata")
load("input/info_icms.Rdata")
group=info_icms
co_sample=intersect(colnames(crc_data_all),group$Tag)
crc_data_all=crc_data_all[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(crc_data_all)
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
load("input/cluster.Rdata")
setwd("input")
load("input/TMB/TMB.Rdata")
TMB$Tag=substr(row.names(TMB),1,12)
TMB=TMB[,-2]
TMB=aggregate(. ~ Tag, data = TMB, FUN = mean)
WNT=merge(TMB,cluster,by.x="Tag",by.y="Tag")
WNT$total_perMB=as.numeric(WNT$total_perMB)
WNT$project="CMS"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("WSS1","WSS2","WSS3","WSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p <- ggplot(WNT, aes(x = cluster, y = log2(total_perMB), fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab(expression("Tumor mutation burden (log"[2] * "-transformation)")) +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("WSS1\n(n=148)",   
                              "WSS2\n(n=137)",  
                              "WSS3\n(n=127)",  
                              "WSS4\n(n=71)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  ggtitle("syn2623706") 
p
WNT=merge(WNT,WNT_Score$final_activity_score,by.x="Tag",by.y="row.names")
library(ggplot2)
library(ggpubr)
colnames(WNT)[5]="WNT/β-catenin pathway activity score"
colnames(WNT)[2]="Tumor mutation burden"
library(GGally) 
p = ggpairs(  
  WNT,  
  columns = c(2, 5),  
  title = "syn2623706",  
  axisLabels = "show",  
  columnLabels = c("Tumor mutation burden", "Wnt/β-catenin pathway activity score"),  
  aes(colour = cluster, alpha = 0.9),  
  upper = list(continuous = wrap("cor", size = 8, display_grid = TRUE)),  
  lower = list(continuous = wrap("smooth", alpha = 0.3, size = 1.5)),  
  diag = list(continuous = "barDiag")  
) +  
  scale_color_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +
  scale_fill_brewer(palette = 1) +    # Customize the diagonal-histogram color  
  theme_bw() +  
  theme(  
    strip.text = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) 
p
# ------------------------------------------------------------------------------
# FGA
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
load("input/FGA.Rdata")
pancancer_df$sampleId=substr(pancancer_df$sampleId,1,12)
WNT=merge(pancancer_df,cluster,by.x="sampleId",by.y="Tag")
WNT$FRACTION_GENOME_ALTERED=as.numeric(WNT$FRACTION_GENOME_ALTERED)
WNT$project="CMS"
table(WNT$cluster)
WNT$cluster=factor(WNT$cluster,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("WSS1","WSS2","WSS3","WSS4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT, aes(x = cluster, y = FRACTION_GENOME_ALTERED, fill = cluster)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Friction genome altered") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("WSS1\n(n=169)",   
                              "WSS2\n(n=157)",  
                              "WSS3\n(n=141)",  
                              "WSS4\n(n=82)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  ggtitle("syn2623706")
p
# ------------------------------------------------------------------------------
# age
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c("sample","age")]
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$age=ifelse(clinical_molecular_public_all$age<=40,"age≤40",
                                         ifelse(clinical_molecular_public_all$age<=60,"40<age≤60","age>60"))

TCGACRC_clinical_merged=read.table(file = "TCGACRC_clinical-merged.tsv",header = T,sep = "\t",check.names = F)
TCGACRC_clinical_merged=TCGACRC_clinical_merged[,c(1,2)]
colnames(TCGACRC_clinical_merged)=c("sample","age")
TCGACRC_clinical_merged$age=ifelse(TCGACRC_clinical_merged$age<=40,"age≤40",
                                   ifelse(TCGACRC_clinical_merged$age<=60,"40<age≤60","age>60"))

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,3)]
colnames(FRENCH_clinical)=c("sample","age")
FRENCH_clinical$age=ifelse(FRENCH_clinical$age<=40,"age≤40",
                           ifelse(FRENCH_clinical$age<=60,"40<age≤60","age>60"))

gse17536_annotation=read.table(file = "gse17536_annotation.tsv",header = T,sep = "\t",check.names = F)
gse17536_annotation=gse17536_annotation[,c(1,2)]
gse17536_annotation$age=ifelse(gse17536_annotation$age<=40,"age≤40",
                               ifelse(gse17536_annotation$age<=60,"40<age≤60","age>60"))

gse2109_annoatation=read.table(file = "gse2109_annoatation.tsv",header = T,sep = "\t",check.names = F)
gse2109_annoatation=gse2109_annoatation[,c(1,32)]
colnames(gse2109_annoatation)=c("sample","age")
gse2109_annoatation$age=ifelse(gse2109_annoatation$age%in%c("30-39","30-40"),"age≤40",
                               ifelse(gse2109_annoatation$age%in%c("40-50","50-59","50-60"),"40<age≤60","age>60"))

gse37892_annotation=read.table(file = "gse37892_annotation.tsv",header = T,sep = "\t",check.names = F)
gse37892_annotation=gse37892_annotation[,c(1,2)]
gse37892_annotation$age=ifelse(gse37892_annotation$age<=40,"age≤40",
                               ifelse(gse37892_annotation$age<=60,"40<age≤60","age>60"))

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,12),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","age")
GSE39582_series_matrix$age=sub(".*: ", "",GSE39582_series_matrix$age)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$age!="N/A",]
GSE39582_series_matrix$age=as.numeric(GSE39582_series_matrix$age)
GSE39582_series_matrix$age=ifelse(GSE39582_series_matrix$age<=40,"age≤40",
                                  ifelse(GSE39582_series_matrix$age<=60,"40<age≤60","age>60"))

GSE33113_series_matrix=read.delim(file = "GSE33113_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE33113_series_matrix=GSE33113_series_matrix[1:34,]
GSE33113_series_matrix=GSE33113_series_matrix[c(1,11),]
GSE33113_series_matrix=as.data.frame(t(GSE33113_series_matrix[,-1]))
colnames(GSE33113_series_matrix)=c("sample","age")
GSE33113_series_matrix$age=sub(".*: ", "",GSE33113_series_matrix$age)
GSE33113_series_matrix=GSE33113_series_matrix[!GSE33113_series_matrix$age%in%c("","col026","col015","col01","col006","col004"),]
GSE33113_series_matrix$age=gsub(",", ".", GSE33113_series_matrix$age)  
GSE33113_series_matrix$age=as.numeric(GSE33113_series_matrix$age)
GSE33113_series_matrix$age=ifelse(GSE33113_series_matrix$age<=40,"age≤40",
                                  ifelse(GSE33113_series_matrix$age<=60,"40<age≤60","age>60"))
age_all=rbind(clinical_molecular_public_all,TCGACRC_clinical_merged)
age_all=rbind(age_all,FRENCH_clinical)
age_all=rbind(age_all,gse17536_annotation)
age_all=rbind(age_all,gse2109_annoatation)
age_all=rbind(age_all,gse37892_annotation)
age_all=rbind(age_all,GSE39582_series_matrix)
age_all=rbind(age_all,GSE33113_series_matrix)
table(age_all$age)
unique_ids=unique(age_all$sample)
age_all=age_all[match(unique_ids, age_all$sample), ]
age=age_all
age=merge(cluster,age,by.x="Tag",by.y="sample")
age=na.omit(age)
age$age=factor(age$age,levels=c("age≤40","40<age≤60","age>60"))
age$status="Age"
WNT_cluster1=age[age$cluster=="WSS1",]
WNT_cluster2=age[age$cluster=="WSS2",]
WNT_cluster3=age[age$cluster=="WSS3",]
WNT_cluster4=age[age$cluster=="WSS4",]
age=age[,-c(1,4)]
age=as.data.frame(table(age))
library(reshape2)
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
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=1702, Age)") +
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
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
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
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
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
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# gender
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c(1,4)]
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$gender=ifelse(clinical_molecular_public_all$gender=="male","Male","Female")

TCGACRC_clinical_merged=read.table(file = "TCGACRC_clinical-merged.tsv",header = T,sep = "\t",check.names = F)
TCGACRC_clinical_merged=TCGACRC_clinical_merged[,c(1,3)]
colnames(TCGACRC_clinical_merged)=c("sample","gender")
TCGACRC_clinical_merged$gender=ifelse(TCGACRC_clinical_merged$gender=="male","Male","Female")

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,4)]
colnames(FRENCH_clinical)=c("sample","gender")
FRENCH_clinical$gender=ifelse(FRENCH_clinical$gender=="M","Male","Female")

gse14333_annotation=read.table(file = "gse14333_annotation.tsv",header = T,sep = "\t",check.names = F)
gse14333_annotation=gse14333_annotation[,c(1,5)]
gse14333_annotation$gender=ifelse(gse14333_annotation$gender=="male","Male","Female")

gse20916_annotation=read.table(file = "gse20916_annotation.tsv",header = T,sep = "\t",check.names = F)
gse20916_annotation=gse20916_annotation[,c(1,3)]
gse20916_annotation$gender=ifelse(gse20916_annotation$gender=="male","Male","Female")

gse17536_annotation=read.table(file = "gse17536_annotation.tsv",header = T,sep = "\t",check.names = F)
gse17536_annotation=gse17536_annotation[,c(1,5)]
gse17536_annotation$gender=ifelse(gse17536_annotation$gender=="male","Male","Female")

gse2109_annoatation=read.table(file = "gse2109_annoatation.tsv",header = T,sep = "\t",check.names = F)
gse2109_annoatation=gse2109_annoatation[,c(1,17)]
gse2109_annoatation=na.omit(gse2109_annoatation)
colnames(gse2109_annoatation)=c("sample","gender")

gse23878_annotation=read.table(file = "gse23878_annotation.tsv",header = T,sep = "\t",check.names = F)
gse23878_annotation=gse23878_annotation[,c(1,2)]
gse23878_annotation$gender=ifelse(gse23878_annotation$gender=="male","Male","Female")

gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,13)]
colnames(gse35896_annotation)=c("sample","gender")
gse35896_annotation$gender=gsub(".*: ", "",gse35896_annotation$gender)

gse37892_annotation=read.table(file = "gse37892_annotation.tsv",header = T,sep = "\t",check.names = F)
gse37892_annotation=gse37892_annotation[,c(1,3)]
gse37892_annotation$gender=ifelse(gse37892_annotation$gender=="male","Male","Female")

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,11),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","gender")
GSE39582_series_matrix$gender=gsub(".*: ", "",GSE39582_series_matrix$gender)

GSE33113_series_matrix=read.delim(file = "GSE33113_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE33113_series_matrix=GSE33113_series_matrix[1:34,]
GSE33113_series_matrix=GSE33113_series_matrix[c(1,12),]
GSE33113_series_matrix=as.data.frame(t(GSE33113_series_matrix[,-1]))
colnames(GSE33113_series_matrix)=c("sample","gender")
GSE33113_series_matrix=GSE33113_series_matrix[GSE33113_series_matrix$gender!="",]
GSE33113_series_matrix$gender=gsub(".*: ", "",GSE33113_series_matrix$gender)
GSE33113_series_matrix$gender=ifelse(GSE33113_series_matrix$gender=="m","Male","Female")

gender_all=rbind(clinical_molecular_public_all,TCGACRC_clinical_merged)
gender_all=rbind(gender_all,FRENCH_clinical)
gender_all=rbind(gender_all,gse14333_annotation)
gender_all=rbind(gender_all,gse20916_annotation)
gender_all=rbind(gender_all,gse17536_annotation)
gender_all=rbind(gender_all,gse2109_annoatation)
gender_all=rbind(gender_all,gse23878_annotation)
gender_all=rbind(gender_all,gse35896_annotation)
gender_all=rbind(gender_all,gse37892_annotation)
gender_all=rbind(gender_all,GSE39582_series_matrix)
gender_all=rbind(gender_all,GSE33113_series_matrix)
table(gender_all$gender)
unique_ids=unique(gender_all$sample)
gender_all=gender_all[match(unique_ids, gender_all$sample), ]
gender=gender_all
gender=merge(cluster,gender,by.x="Tag",by.y="sample")
gender=gender[,c("cluster","gender")]
gender=na.omit(gender)
gender$status="gender"
WNT_cluster1=gender[gender$cluster=="WSS1",]
WNT_cluster2=gender[gender$cluster=="WSS2",]
WNT_cluster3=gender[gender$cluster=="WSS3",]
WNT_cluster4=gender[gender$cluster=="WSS4",]
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
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=2012, Gender)") +
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
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
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
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
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
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot)  
# ------------------------------------------------------------------------------
# APC
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(reshape2)
load("input/cluster.Rdata")
setwd("input")
gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,8)]
colnames(gse35896_annotation)=c("sample","APC")
gse35896_annotation$APC=gsub(".*: ", "",gse35896_annotation$APC)
gse35896_annotation$APC=ifelse(gse35896_annotation$APC=="Y","MT","WT")

COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("KRAS","BRAF","TP53","APC"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),APC=COADREAD_mc3_gene_level$APC)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

APC_all=rbind(gse35896_annotation,COADREAD_mc3_gene_level)
table(APC_all$APC)
unique_ids=unique(APC_all$sample)
APC_all=APC_all[match(unique_ids, APC_all$sample), ]
APC=APC_all
APC=merge(cluster,APC,by.x="Tag",by.y="sample")
APC=APC[,c("cluster","APC")]
APC=na.omit(APC)
APC$status="APC"
WNT_cluster1=APC[APC$cluster=="WSS1",]
WNT_cluster2=APC[APC$cluster=="WSS2",]
WNT_cluster3=APC[APC$cluster=="WSS3",]
WNT_cluster4=APC[APC$cluster=="WSS4",]
APC=APC[,-3]
APC=as.data.frame(table(APC))
contingency_table=acast(APC, cluster ~ APC, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = APC)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=388, APC)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = APC)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = APC)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = APC)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# TP53
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,9)]
colnames(gse35896_annotation)=c("sample","TP53")
gse35896_annotation$TP53=gsub(".*: ", "",gse35896_annotation$TP53)
gse35896_annotation$TP53=ifelse(gse35896_annotation$TP53=="Y","MT","WT")

COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("KRAS","BRAF","TP53","APC"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),TP53=COADREAD_mc3_gene_level$TP53)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,20)]
colnames(FRENCH_clinical)=c("sample","TP53")
FRENCH_clinical=FRENCH_clinical[FRENCH_clinical$TP53!="ND",]
FRENCH_clinical$TP53=ifelse(FRENCH_clinical$TP53=="M","MT","WT")

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,27),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","TP53")
GSE39582_series_matrix$TP53=gsub(".*: ", "",GSE39582_series_matrix$TP53)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$TP53!="N/A",]
GSE39582_series_matrix$TP53=ifelse(GSE39582_series_matrix$TP53=="M","MT","WT")

TP53_all=rbind(gse35896_annotation,COADREAD_mc3_gene_level)
TP53_all=rbind(TP53_all,FRENCH_clinical)
TP53_all=rbind(TP53_all,GSE39582_series_matrix)
table(TP53_all$TP53)
unique_ids=unique(TP53_all$sample)
TP53_all=TP53_all[match(unique_ids, TP53_all$sample), ]
TP53=TP53_all
TP53=merge(cluster,TP53,by.x="Tag",by.y="sample")
TP53=TP53[,c("cluster","TP53")]
TP53=na.omit(TP53)
TP53$status="TP53"
WNT_cluster1=TP53[TP53$cluster=="WSS1",]
WNT_cluster2=TP53[TP53$cluster=="WSS2",]
WNT_cluster3=TP53[TP53$cluster=="WSS3",]
WNT_cluster4=TP53[TP53$cluster=="WSS4",]
TP53=TP53[,-3]
TP53=as.data.frame(table(TP53))
contingency_table=acast(TP53, cluster ~ TP53, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = TP53)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=738, TP53)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = TP53)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = TP53)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = TP53)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# BRAF
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,7)]
colnames(gse35896_annotation)=c("sample","BRAF")
gse35896_annotation$BRAF=gsub(".*: ", "",gse35896_annotation$BRAF)
gse35896_annotation$BRAF=ifelse(gse35896_annotation$BRAF=="Y","MT","WT")

clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c(1,14)]
colnames(clinical_molecular_public_all)=c("sample","BRAF")
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$BRAF=ifelse(clinical_molecular_public_all$BRAF==1,"MT","WT")

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,18)]
colnames(FRENCH_clinical)=c("sample","BRAF")
FRENCH_clinical=FRENCH_clinical[FRENCH_clinical$BRAF!="ND",]
FRENCH_clinical$BRAF=ifelse(FRENCH_clinical$BRAF=="M","MT","WT")

COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("KRAS","BRAF","TP53","APC"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),BRAF=COADREAD_mc3_gene_level$BRAF)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,35),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","BRAF")
GSE39582_series_matrix$BRAF=gsub(".*: ", "",GSE39582_series_matrix$BRAF)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$BRAF!="N/A",]
GSE39582_series_matrix$BRAF=ifelse(GSE39582_series_matrix$BRAF=="M","MT","WT")

BRAF_all=rbind(gse35896_annotation,clinical_molecular_public_all)
BRAF_all=rbind(BRAF_all,FRENCH_clinical)
BRAF_all=rbind(BRAF_all,COADREAD_mc3_gene_level)
BRAF_all=rbind(BRAF_all,GSE39582_series_matrix)
table(BRAF_all$BRAF)
unique_ids=unique(BRAF_all$sample)
BRAF_all=BRAF_all[match(unique_ids, BRAF_all$sample), ]
BRAF=BRAF_all
BRAF=merge(cluster,BRAF,by.x="Tag",by.y="sample")
BRAF=BRAF[,c("cluster","BRAF")]
BRAF=na.omit(BRAF)
BRAF$status="BRAF"
WNT_cluster1=BRAF[BRAF$cluster=="WSS1",]
WNT_cluster2=BRAF[BRAF$cluster=="WSS2",]
WNT_cluster3=BRAF[BRAF$cluster=="WSS3",]
WNT_cluster4=BRAF[BRAF$cluster=="WSS4",]
BRAF=BRAF[,-3]
BRAF=as.data.frame(table(BRAF))
contingency_table=acast(BRAF, cluster ~ BRAF, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = BRAF)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=976, BRAF)") + 
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none", 
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = BRAF)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = BRAF)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = BRAF)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# KRAS
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,6)]
colnames(gse35896_annotation)=c("sample","KRAS")
gse35896_annotation$KRAS=gsub(".*: ", "",gse35896_annotation$KRAS)
gse35896_annotation$KRAS=ifelse(gse35896_annotation$KRAS=="Y","MT","WT")

clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c(1,13)]
colnames(clinical_molecular_public_all)=c("sample","KRAS")
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$KRAS=ifelse(clinical_molecular_public_all$KRAS==1,"MT","WT")

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,16)]
colnames(FRENCH_clinical)=c("sample","KRAS")
FRENCH_clinical=FRENCH_clinical[FRENCH_clinical$KRAS!="ND",]
FRENCH_clinical$KRAS=ifelse(FRENCH_clinical$KRAS=="M","MT","WT")

COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("KRAS","BRAF","TP53","APC"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),KRAS=COADREAD_mc3_gene_level$KRAS)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,31),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","KRAS")
GSE39582_series_matrix$KRAS=gsub(".*: ", "",GSE39582_series_matrix$KRAS)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$KRAS!="N/A",]
GSE39582_series_matrix$KRAS=ifelse(GSE39582_series_matrix$KRAS=="M","MT","WT")

KRAS_all=rbind(gse35896_annotation,clinical_molecular_public_all)
KRAS_all=rbind(KRAS_all,FRENCH_clinical)
KRAS_all=rbind(KRAS_all,COADREAD_mc3_gene_level)
KRAS_all=rbind(KRAS_all,GSE39582_series_matrix)
table(KRAS_all$KRAS)
unique_ids=unique(KRAS_all$sample)
KRAS_all=KRAS_all[match(unique_ids, KRAS_all$sample), ]

KRAS=KRAS_all
KRAS=merge(cluster,KRAS,by.x="Tag",by.y="sample")
KRAS=KRAS[,c("cluster","KRAS")]
KRAS=na.omit(KRAS)
KRAS$status="KRAS"
WNT_cluster1=KRAS[KRAS$cluster=="WSS1",]
WNT_cluster2=KRAS[KRAS$cluster=="WSS2",]
WNT_cluster3=KRAS[KRAS$cluster=="WSS3",]
WNT_cluster4=KRAS[KRAS$cluster=="WSS4",]
KRAS=KRAS[,-3]
KRAS=as.data.frame(table(KRAS))
contingency_table=acast(KRAS, cluster ~ KRAS, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = KRAS)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=1009, KRAS)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none", 
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = KRAS)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none", 
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = KRAS)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = KRAS)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# PIK3CA
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")

gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,10)]
colnames(gse35896_annotation)=c("sample","PIK3CA")
gse35896_annotation$PIK3CA=ifelse(gse35896_annotation$PIK3CA=="pik3ca.mutation: Y","MT","WT")

COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("KRAS","BRAF","TP53","APC","PIK3CA"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),PIK3CA=COADREAD_mc3_gene_level$PIK3CA)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

PIK3CA_all=gse35896_annotation
PIK3CA_all=rbind(PIK3CA_all,COADREAD_mc3_gene_level)
table(PIK3CA_all$PIK3CA)
unique_ids=unique(PIK3CA_all$sample)
PIK3CA_all=PIK3CA_all[match(unique_ids, PIK3CA_all$sample), ]

PIK3CA=PIK3CA_all
PIK3CA=merge(cluster,PIK3CA,by.x="Tag",by.y="sample")
PIK3CA=PIK3CA[,c("cluster","PIK3CA")]
PIK3CA=na.omit(PIK3CA)
PIK3CA$status="PIK3CA"
WNT_cluster1=PIK3CA[PIK3CA$cluster=="WSS1",]
WNT_cluster2=PIK3CA[PIK3CA$cluster=="WSS2",]
WNT_cluster3=PIK3CA[PIK3CA$cluster=="WSS3",]
WNT_cluster4=PIK3CA[PIK3CA$cluster=="WSS4",]
PIK3CA=PIK3CA[,-3]
PIK3CA=as.data.frame(table(PIK3CA))
contingency_table=acast(PIK3CA, cluster ~ PIK3CA, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = PIK3CA)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=388, PIK3CA)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = PIK3CA)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = PIK3CA)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = PIK3CA)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# PTEN
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
gse35896_annotation=read.table(file = "gse35896_annotation.tsv",header = T,sep = "\t",check.names = F)
gse35896_annotation=gse35896_annotation[,c(1,11)]
colnames(gse35896_annotation)=c("sample","PTEN")
gse35896_annotation$PTEN=ifelse(gse35896_annotation$PTEN=="pten.mutation: Y","MT","WT")

COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("KRAS","BRAF","TP53","APC","PTEN"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),PTEN=COADREAD_mc3_gene_level$PTEN)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

PTEN_all=rbind(COADREAD_mc3_gene_level,gse35896_annotation)
table(PTEN_all$PIK3CA)
unique_ids=unique(PTEN_all$sample)
PTEN_all=PTEN_all[match(unique_ids, PTEN_all$sample), ]

PTEN=PTEN_all
PTEN=merge(cluster,PTEN,by.x="Tag",by.y="sample")
PTEN=PTEN[,c("cluster","PTEN")]
PTEN=na.omit(PTEN)
PTEN$status="PTEN"
WNT_cluster1=PTEN[PTEN$cluster=="WSS1",]
WNT_cluster2=PTEN[PTEN$cluster=="WSS2",]
WNT_cluster3=PTEN[PTEN$cluster=="WSS3",]
WNT_cluster4=PTEN[PTEN$cluster=="WSS4",]
PTEN=PTEN[,-3]
PTEN=as.data.frame(table(PTEN))
contingency_table=acast(PTEN, cluster ~ PTEN, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = PTEN)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=388, PTEN)") + 
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = PTEN)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = PTEN)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = PTEN)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# RNF43
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("RNF43","BRAF","TP53","APC","PTEN"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),RNF43=COADREAD_mc3_gene_level$RNF43)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

RNF43_all=COADREAD_mc3_gene_level
table(RNF43_all$RNF43)
unique_ids=unique(RNF43_all$sample)
RNF43_all=RNF43_all[match(unique_ids, RNF43_all$sample), ]

RNF43=RNF43_all
RNF43=merge(cluster,RNF43,by.x="Tag",by.y="sample")
RNF43=RNF43[,c("cluster","RNF43")]
RNF43=na.omit(RNF43)
RNF43$status="RNF43"
WNT_cluster1=RNF43[RNF43$cluster=="WSS1",]
WNT_cluster2=RNF43[RNF43$cluster=="WSS2",]
WNT_cluster3=RNF43[RNF43$cluster=="WSS3",]
WNT_cluster4=RNF43[RNF43$cluster=="WSS4",]
RNF43=RNF43[,-3]
RNF43=as.data.frame(table(RNF43))
contingency_table=acast(RNF43, cluster ~ RNF43, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = RNF43)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=326, RNF43)") + 
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = RNF43)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = RNF43)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = RNF43)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# ZNRF3
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
COADREAD_mc3_gene_level=read.delim(file = "COADREAD_mc3_gene_level.txt",header = T,sep = "\t",check.names = F)
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[COADREAD_mc3_gene_level$xena_sample%in%c("RNF43","ZNRF3","TP53","APC","PTEN"),]
rownames(COADREAD_mc3_gene_level)=COADREAD_mc3_gene_level[,1]
COADREAD_mc3_gene_level=COADREAD_mc3_gene_level[,-1]
COADREAD_mc3_gene_level=t(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=ifelse(COADREAD_mc3_gene_level==1,"MT","WT")
COADREAD_mc3_gene_level=as.data.frame(COADREAD_mc3_gene_level)
COADREAD_mc3_gene_level=data.frame(sample=rownames(COADREAD_mc3_gene_level),ZNRF3=COADREAD_mc3_gene_level$ZNRF3)
COADREAD_mc3_gene_level$sample=substr(COADREAD_mc3_gene_level$sample,1,12)

ZNRF3_all=COADREAD_mc3_gene_level
table(ZNRF3_all$ZNRF3)
unique_ids=unique(ZNRF3_all$sample)
ZNRF3_all=ZNRF3_all[match(unique_ids, ZNRF3_all$sample), ]

ZNRF3=ZNRF3_all
ZNRF3=merge(cluster,ZNRF3,by.x="Tag",by.y="sample")
ZNRF3=ZNRF3[,c("cluster","ZNRF3")]
ZNRF3=na.omit(ZNRF3)
ZNRF3$status="ZNRF3"
WNT_cluster1=ZNRF3[ZNRF3$cluster=="WSS1",]
WNT_cluster2=ZNRF3[ZNRF3$cluster=="WSS2",]
WNT_cluster3=ZNRF3[ZNRF3$cluster=="WSS3",]
WNT_cluster4=ZNRF3[ZNRF3$cluster=="WSS4",]
ZNRF3=ZNRF3[,-3]
ZNRF3=as.data.frame(table(ZNRF3))
contingency_table=acast(ZNRF3, cluster ~ ZNRF3, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = ZNRF3)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=326, ZNRF3)") + 
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )  
p2 = ggplot(WNT_cluster2, aes(x = status, fill = ZNRF3)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p3 = ggplot(WNT_cluster3, aes(x = status, fill = ZNRF3)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )  
p4 = ggplot(WNT_cluster4, aes(x = status, fill = ZNRF3)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("MT" = "#C24976", "WT" = "#469393")) +  
  theme_minimal() +   
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  ) 
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# stage
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")

clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c("sample","stage")]
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$stage=paste0("Stage",clinical_molecular_public_all$stage)

TCGACRC_clinical_merged=read.table(file = "TCGACRC_clinical-merged.tsv",header = T,sep = "\t",check.names = F)
TCGACRC_clinical_merged=TCGACRC_clinical_merged[,c(1,4)]
TCGACRC_clinical_merged=na.omit(TCGACRC_clinical_merged)
TCGACRC_clinical_merged$stage=ifelse(TCGACRC_clinical_merged$stage%in%c("Stage I","Stage IA"),1,
                                     ifelse(TCGACRC_clinical_merged$stage%in%c("Stage II","Stage IIA","Stage IIB","Stage IIC"),2,
                                            ifelse(TCGACRC_clinical_merged$stage%in%c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC"),3,
                                                   4)))
colnames(TCGACRC_clinical_merged)=c("sample","stage")
TCGACRC_clinical_merged$stage=paste0("Stage",TCGACRC_clinical_merged$stage)

FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,5)]
colnames(FRENCH_clinical)=c("sample","stage")
FRENCH_clinical$stage=paste0("Stage",FRENCH_clinical$stage)

gse17536_annotation=read.table(file = "gse17536_annotation.tsv",header = T,sep = "\t",check.names = F)
gse17536_annotation=gse17536_annotation[,c(1,3)]
colnames(gse17536_annotation)=c("sample","stage")
gse17536_annotation$stage=paste0("Stage",gse17536_annotation$stage)

gse2109_annoatation=read.table(file = "gse2109_annoatation.tsv",header = T,sep = "\t",check.names = F)
gse2109_annoatation=gse2109_annoatation[,c(1,28)]
gse2109_annoatation=na.omit(gse2109_annoatation)
gse2109_annoatation=gse2109_annoatation[gse2109_annoatation$`Pathological Stage`!=c("Unknown"),]
table(gse2109_annoatation$`Pathological Stage`)
gse2109_annoatation$`Pathological Stage` <- gsub("[A-Z]", "", gse2109_annoatation$`Pathological Stage`)  
colnames(gse2109_annoatation)=c("sample","stage")
gse2109_annoatation$stage=paste0("Stage",gse2109_annoatation$stage)

gse37892_annotation=read.table(file = "gse37892_annotation.tsv",header = T,sep = "\t",check.names = F)
gse37892_annotation=gse37892_annotation[,c(1,5)]
colnames(gse37892_annotation)=c("sample","stage")
gse37892_annotation$stage=paste0("Stage",gse37892_annotation$stage)

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,13),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","stage")
GSE39582_series_matrix$stage=gsub(".*: ", "",GSE39582_series_matrix$stage)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$stage!="N/A",]
GSE39582_series_matrix$stage=paste0("Stage",GSE39582_series_matrix$stage)

GSE33113_series_matrix=read.delim(file = "GSE33113_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE33113_series_matrix=GSE33113_series_matrix[1:34,]
GSE33113_series_matrix=GSE33113_series_matrix[c(1,9),]
GSE33113_series_matrix=as.data.frame(t(GSE33113_series_matrix[,-1]))
colnames(GSE33113_series_matrix)=c("sample","stage")
GSE33113_series_matrix$stage="Stage2"

stage_all=rbind(clinical_molecular_public_all,TCGACRC_clinical_merged)
stage_all=rbind(stage_all,FRENCH_clinical)
stage_all=rbind(stage_all,gse17536_annotation)
stage_all=rbind(stage_all,gse2109_annoatation)
stage_all=rbind(stage_all,gse37892_annotation)
stage_all=rbind(stage_all,GSE39582_series_matrix)
stage_all=rbind(stage_all,GSE33113_series_matrix)
table(stage_all$stage)
unique_ids=unique(stage_all$sample)
stage_all=stage_all[match(unique_ids, stage_all$sample), ]

stage=stage_all
stage=merge(cluster,stage,by.x="Tag",by.y="sample")
stage=na.omit(stage)
stage$status="stage"
WNT_cluster1=stage[stage$cluster=="WSS1",]
WNT_cluster2=stage[stage$cluster=="WSS2",]
WNT_cluster3=stage[stage$cluster=="WSS3",]
WNT_cluster4=stage[stage$cluster=="WSS4",]
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
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=1845, Stage)") +
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
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
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
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
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
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )  
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# grade
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")

clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c("sample","grade")]
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
clinical_molecular_public_all$grade=paste0("Grade",clinical_molecular_public_all$grade)

gse17536_annotation=read.table(file = "gse17536_annotation.tsv",header = T,sep = "\t",check.names = F)
gse17536_annotation=gse17536_annotation[,c(1,6)]
gse17536_annotation$grade=ifelse(gse17536_annotation$grade=="1 - Well differentiated (WD)",1,
                                 ifelse(gse17536_annotation$grade=="2 - Moderately differentiated (MD)",2,3))
gse17536_annotation$grade=paste0("Grade",gse17536_annotation$grade)

grade_all=rbind(clinical_molecular_public_all,gse17536_annotation)

grade=grade_all
grade=merge(cluster,grade,by.x="Tag",by.y="sample")
grade=na.omit(grade)
grade$status="grade"
WNT_cluster1=grade[grade$cluster=="WSS1",]
WNT_cluster2=grade[grade$cluster=="WSS2",]
WNT_cluster3=grade[grade$cluster=="WSS3",]
WNT_cluster4=grade[grade$cluster=="WSS4",]
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
  scale_fill_manual(values = c("Grade1" = "#afb4db", "Grade2" = "#9b95c9", "Grade3" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=174, Grade)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  ) 
p2 = ggplot(WNT_cluster2, aes(x = status, fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Grade1" = "#afb4db", "Grade2" = "#9b95c9", "Grade3" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  ) 
p3 = ggplot(WNT_cluster3, aes(x = status, fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Grade1" = "#afb4db", "Grade2" = "#9b95c9", "Grade3" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = grade)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Grade1" = "#afb4db", "Grade2" = "#9b95c9", "Grade3" = "#6950a1")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# CIMP cluster
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
clinical_molecular_public_all=read.table(file = "clinical_molecular_public_all.txt",header = T,sep = "\t",check.names = F)
clinical_molecular_public_all=clinical_molecular_public_all[,c("sample","cimp")]
clinical_molecular_public_all=na.omit(clinical_molecular_public_all)
cimp=clinical_molecular_public_all
cimp=merge(cluster,cimp,by.x="Tag",by.y="sample")
cimp=na.omit(cimp)
cimp$status="cimp"
WNT_cluster1=cimp[cimp$cluster=="WSS1",]
WNT_cluster2=cimp[cimp$cluster=="WSS2",]
WNT_cluster3=cimp[cimp$cluster=="WSS3",]
WNT_cluster4=cimp[cimp$cluster=="WSS4",]
cimp=cimp[,-c(1,4)]
cimp=as.data.frame(table(cimp))
contingency_table=acast(cimp, cluster ~ cimp, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = cimp)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CIMP.High" = "#a7324a", "CIMP.Low" = "#f58f98", "CIMP.Neg" = "#a1a3a6")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=957, CIMP cluster)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = cimp)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CIMP.High" = "#a7324a", "CIMP.Low" = "#f58f98", "CIMP.Neg" = "#a1a3a6")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = cimp)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CIMP.High" = "#a7324a", "CIMP.Low" = "#f58f98", "CIMP.Neg" = "#a1a3a6")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = cimp)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("CIMP.High" = "#a7324a", "CIMP.Low" = "#f58f98", "CIMP.Neg" = "#a1a3a6")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# site
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
TCGACRC_clinical_merged=read.table(file = "TCGACRC_clinical-merged.tsv",header = T,sep = "\t",check.names = F)
TCGACRC_clinical_merged=TCGACRC_clinical_merged[,c(1,8)]
TCGACRC_clinical_merged=na.omit(TCGACRC_clinical_merged)
table(TCGACRC_clinical_merged$tumorLocation)
colnames(TCGACRC_clinical_merged)=c("sample","site")
TCGACRC_clinical_merged$site=ifelse(TCGACRC_clinical_merged$site%in%c("Cecum","Ascending Colon","Hepatic Flexure","Transverse Colon"),"Right Colon",
                                    ifelse(TCGACRC_clinical_merged$site%in%c("Splenic Flexure","Descending Colon","Sigmoid Colon"),"Left Colon",
                                           ifelse(TCGACRC_clinical_merged$site=="Rectum","Rectum","Rectosigmoid Junction")))
FRENCH_clinical=read.table(file = "FRENCH_clinical.tsv",header = T,sep = "\t",check.names = F)
FRENCH_clinical=FRENCH_clinical[,c(1,9)]
table(FRENCH_clinical$tumorLocation)
colnames(FRENCH_clinical)=c("sample","site")
FRENCH_clinical$site=ifelse(FRENCH_clinical$site=="distal","Left Colon","Right Colon")

gse14333_annotation=read.table(file = "gse14333_annotation.tsv",header = T,sep = "\t",check.names = F)
gse14333_annotation=gse14333_annotation[,c(1,6)]
colnames(gse14333_annotation)=c("sample","site")
gse14333_annotation=na.omit(gse14333_annotation)
table(gse14333_annotation$site)
gse14333_annotation=gse14333_annotation[gse14333_annotation$site!="colon",]
gse14333_annotation$site=ifelse(gse14333_annotation$site=="right","Right Colon",
                                ifelse(gse14333_annotation$site=="left","Left Colon","Rectum"))

gse37892_annotation=read.table(file = "gse37892_annotation.tsv",header = T,sep = "\t",check.names = F)
gse37892_annotation=gse37892_annotation[,c(1,4)]
gse37892_annotation=na.omit(gse37892_annotation)
colnames(gse37892_annotation)=c("sample","site")
gse37892_annotation$site=ifelse(gse37892_annotation$site=="right","Right Colon",
                                ifelse(gse37892_annotation$site=="left","Left Colon","Rectum"))

GSE42284_series_matrix=read.delim(file = "GSE42284_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE42284_series_matrix=GSE42284_series_matrix[1:48,]
GSE42284_series_matrix_a=GSE42284_series_matrix[c(1,14),]
GSE42284_series_matrix_a=as.data.frame(t(GSE42284_series_matrix_a[,-1]))
colnames(GSE42284_series_matrix_a)=c("sample","site")
table(GSE42284_series_matrix_a$site)
GSE42284_series_matrix_a=GSE42284_series_matrix_a[GSE42284_series_matrix_a$site!="braf: wildtype / other",]
GSE42284_series_matrix_a$site=ifelse(GSE42284_series_matrix_a$site=="location: left colon","Left Colon",
                                     ifelse(GSE42284_series_matrix_a$site=="location: right colon","Right Colon","Rectum"))

GSE39582_series_matrix=read.delim(file = "GSE39582_series_matrix.txt",header = T,sep = "\t",check.names = F)
GSE39582_series_matrix=GSE39582_series_matrix[c(1,17),]
GSE39582_series_matrix=as.data.frame(t(GSE39582_series_matrix[,-1]))
colnames(GSE39582_series_matrix)=c("sample","site")
table(GSE39582_series_matrix$site)
GSE39582_series_matrix=GSE39582_series_matrix[GSE39582_series_matrix$site!="tumor.location: N/A",]
GSE39582_series_matrix$site=ifelse(GSE39582_series_matrix$site=="tumor.location: distal","Left Colon",
                                   ifelse(GSE39582_series_matrix$site=="tumor.location: proximal","Right Colon","Rectum"))

site_all=rbind(TCGACRC_clinical_merged,FRENCH_clinical)
site_all=rbind(site_all,gse14333_annotation)
site_all=rbind(site_all,gse37892_annotation)
site_all=rbind(site_all,GSE42284_series_matrix_a)
site_all=rbind(site_all,GSE39582_series_matrix)
table(site_all$KRAS)
unique_ids=unique(site_all$sample)
site_all=site_all[match(unique_ids, site_all$sample), ]

site=site_all
site=merge(cluster,site,by.x="Tag",by.y="sample")
site=na.omit(site)
site$status="site"
site$site=factor(site$site,level=c("Left Colon","Right Colon","Rectum","Rectosigmoid Junction"))
WNT_cluster1=site[site$cluster=="WSS1",]
WNT_cluster2=site[site$cluster=="WSS2",]
WNT_cluster3=site[site$cluster=="WSS3",]
WNT_cluster4=site[site$cluster=="WSS4",]
site=site[,-c(1,4)]
site=as.data.frame(table(site))
contingency_table=reshape2::acast(site, cluster ~ site, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = site)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Left Colon" = "#7fafc6", "Right Colon" = "#f7c67d",   
                               "Rectum" = "#b4d9b9", "Rectosigmoid Junction" = "#9a82bb")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=1388, Site)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = site)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Left Colon" = "#7fafc6", "Right Colon" = "#f7c67d",   
                               "Rectum" = "#b4d9b9", "Rectosigmoid Junction" = "#9a82bb")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = site)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Left Colon" = "#7fafc6", "Right Colon" = "#f7c67d",   
                               "Rectum" = "#b4d9b9", "Rectosigmoid Junction" = "#9a82bb")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = site)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Left Colon" = "#7fafc6", "Right Colon" = "#f7c67d",   
                               "Rectum" = "#b4d9b9", "Rectosigmoid Junction" = "#9a82bb")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
library(patchwork)
final_plot = p1 | p2 | p3 | p4  
print(final_plot) 
# ------------------------------------------------------------------------------
# alcohol&tobacco
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
load("input/cluster.Rdata")
setwd("input")
gse2109_annoatation=read.table(file = "gse2109_annoatation.tsv",header = T,sep = "\t",check.names = F)
gse2109_annoatation=gse2109_annoatation[,c(1,2,50)]
gse2109_annoatation$bad_habit=paste0(gse2109_annoatation$`Alcohol Consumption?`,gse2109_annoatation$`Tobacco Use`)
gse2109_annoatation$bad_habit=ifelse(gse2109_annoatation$bad_habit=="YesYes","Alcohol&Tobacco",
                                     ifelse(gse2109_annoatation$bad_habit=="YesNo","Alcohol",
                                            ifelse(gse2109_annoatation$bad_habit=="NoYes","Tobacco","None")))
bad_habit=merge(cluster,gse2109_annoatation,by.x="Tag",by.y="sample")
bad_habit=na.omit(bad_habit)
bad_habit$status="bad_habit"
bad_habit$bad_habit=factor(bad_habit$bad_habit,level=c("Alcohol&Tobacco","Alcohol","Tobacco","None"))
WNT_cluster1=bad_habit[bad_habit$cluster=="WSS1",]
WNT_cluster2=bad_habit[bad_habit$cluster=="WSS2",]
WNT_cluster3=bad_habit[bad_habit$cluster=="WSS3",]
WNT_cluster4=bad_habit[bad_habit$cluster=="WSS4",]
bad_habit=bad_habit[,c(2,5)]
bad_habit=as.data.frame(table(bad_habit))
contingency_table=reshape2::acast(bad_habit, cluster ~ bad_habit, value.var = "Freq", fill = 0) 
chi_square_result=chisq.test(contingency_table)
chi_square_result=chi_square_result$p.value
p.lab=paste0("chi-square test p",  
             ifelse(chi_square_result < 0.001, " < 0.001",
                    paste0(" = ", round(chi_square_result, 3)))) 
library(ggplot2)
p1 = ggplot(WNT_cluster1, aes(x = "", fill = bad_habit)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Alcohol&Tobacco" = "#ed1941", "Alcohol" = "#ba8448",   
                               "Tobacco" = "#f15a22", "None" = "#afdfe4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS1)",subtitle = paste0(p.lab)) +  
  ggtitle("syn2623706(n=286, Bad Habit)") +
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16),
    plot.subtitle = element_text(size = 16, face = "italic")
  )
p2 = ggplot(WNT_cluster2, aes(x = status, fill = bad_habit)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Alcohol&Tobacco" = "#ed1941", "Alcohol" = "#ba8448",   
                               "Tobacco" = "#f15a22", "None" = "#afdfe4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS2)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p3 = ggplot(WNT_cluster3, aes(x = status, fill = bad_habit)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Alcohol&Tobacco" = "#ed1941", "Alcohol" = "#ba8448",   
                               "Tobacco" = "#f15a22", "None" = "#afdfe4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS3)") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(color = "black", size = 16)
  )
p4 = ggplot(WNT_cluster4, aes(x = status, fill = bad_habit)) +  
  geom_bar(position = "fill", color = "black") +  
  scale_fill_manual(values = c("Alcohol&Tobacco" = "#ed1941", "Alcohol" = "#ba8448",   
                               "Tobacco" = "#f15a22", "None" = "#afdfe4")) +  
  theme_minimal() +  
  labs(title = "", x = "", y = "Proportion(WSS4)", fill = "") +  
  theme(  
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(color = "black", size = 16)
  )
final_plot = p1 | p2 | p3 | p4  
print(final_plot)
# ------------------------------------------------------------------------------
# drug
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
load("input/cluster.Rdata")
load("input/crc_data_all.Rdata")
exp=as.data.frame(crc_data_all)
drug=getGmt("CRC_drug.gmt")
drug_cluster=gsva(expr=as.matrix(exp), drug, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
drug_cluster=as.data.frame(t(drug_cluster))
drug_cluster=as.data.frame(scale(drug_cluster))
drug_cluster_wnt=merge(cluster,drug_cluster,by.x="Tag",by="row.names")
colnames(drug_cluster_wnt)[3:length(colnames(drug_cluster_wnt))]=toupper(colnames(drug_cluster_wnt)[3:length(colnames(drug_cluster_wnt))])
library(ggplot2)
library(ggpubr)
library(patchwork) 
library(tidyr)
library(dplyr)
data_long=drug_cluster_wnt %>%
  pivot_longer(
    cols = 3:32,
    names_to = "variable",
    values_to = "value"
  )
p = ggplot(data_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    color = "black"
  ) +
  stat_compare_means(
    aes(label = ..p.signif..),
    label.y = sapply(levels(data_long$variable), function(var) {
      min(data_long$value[data_long$variable == var], na.rm = TRUE) * 0.95
    }),
    label.x = 2.5,
    size = 6,
    vjust = 0.5
  ) +
  facet_wrap(~variable, nrow = 5, ncol = 6, scales = "free_y") +
  theme_bw() +
  labs(x = "", y = "ssGSEA score") +
  theme(
    legend.position = "none",
    text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(
      size = 16, 
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 8, color = "black", face = "bold")
  ) +
  scale_fill_manual(values = c("WSS1"="#d71345","WSS2"= "#45b97c", "WSS3"="#145b7d","WSS4"= "#ffd400"))
print(p)
# ------------------------------------------------------------------------------
# Large-scale CRC drug
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("exp_data")
vector=c("GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE21815","GSE17537","GSE29621",
         "GSE38832","GSE72970","GSE33113",
         "UCSC.CRC_CMS2","GSE39582","GSE2109","GSE13067",
         "GSE13294","GSE14333","GSE17536",
         "GSE20916","GSE23878","GSE35896",
         "GSE37892","KFSYSCC","syn26720761",
         "syn2623706")
for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
cancer_cor_all=data.frame()
for (i in 1:length(vector)){
  rt=get(paste0(vector[i],"_CMS"))
  WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                         expression_profile=rt,
                                         foundation="relative_ssGSEA",
                                         activation_geneset=NA,
                                         inhibition_geneset=NA,
                                         geneSets_gmt=geneSets,
                                         min.sz=1,
                                         max.sz=10000,
                                         geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                             rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  library(GSVA)
  library(GSEABase)
  library(reshape)
  setwd("input")
  GeneSets=getGmt("CRC_drug.gmt")
  geneset=gsva(expr=as.matrix(rt), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  geneset=as.data.frame(scale(geneset))
  fs=unique(colnames(geneset))
  cancer_cor=data.frame()
  for(j in 1:length(fs)){
    cor_test=cor.test(WNT_Score$final_activity_score$activity_score,geneset[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[j],tme="activity_score",cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_cor=rbind(cancer_cor,cancer_tme)
  }
  cancer_cor$cor[is.na(cancer_cor$cor)]=0
  cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
  cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "", 
                          ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*", 
                                 ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
  cancer_cor$dataset=vector[i]
  cancer_cor_all=rbind(cancer_cor_all,cancer_cor)
}
library(ggplot2)
library(dplyr)
cancer_cor_all$dataset=ifelse(cancer_cor_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",cancer_cor_all$dataset)
cancer_cor_all$cancer=toupper(cancer_cor_all$cancer)
cancer_cor_all$cancer=factor(cancer_cor_all$cancer,levels = unique(cancer_cor_all$cancer))
ggplot(cancer_cor_all, aes(dataset, cancer)) + 
  geom_tile(aes(fill = cor))+
  scale_fill_gradient2(high = "#f58f98",mid = "white",low = "#90d7ec")+
  geom_text(aes(label=pstar),col ="black",size =3)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))
# ------------------------------------------------------------------------------
# molecular characterisation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
library(reshape)
setwd("input")
load("input/cluster.Rdata")
load("input/crc_data_all.Rdata")
GeneSets=getGmt("NM_genesets.gmt")
geneset=gsva(expr=as.matrix(crc_data_all), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
geneset=as.data.frame(t(geneset))
geneset=as.data.frame(scale(geneset))
geneset=geneset[,-c(4,5)]
name=as.data.frame(colnames(geneset))
# ------------------------------------------------------------------------------
# Infiltration ESTIMATE
# ------------------------------------------------------------------------------
geneset_zdk=geneset[,c("IMMUNE_ESTIMATE","STROMAL_ESTIMATE")]
colnames(geneset_zdk)=c("Immune infiltration","Stromal infiltration")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"), 
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Infiltration ESTIMATE")
p 
# ------------------------------------------------------------------------------
# Signatures
# ------------------------------------------------------------------------------
geneset_zdk=geneset[,c("SERRATED_UP","EMT_CORE_GENES","MYC_TARGETS_ZELLER","TGFB_CORE_GENES","CSC_BATLLE",
                       "MATRIX_REMODEL_REACTOME","WOUND_RESPONSE_GO_BP","CRYPT_BASE","CRYPT_TOP","EPITH_LOBODA","MESENCH_LOBODA")]
colnames(geneset_zdk)=c("Serrated lesion","EMT activation","MYC targets","TGF-β activation","Cancer stem cell",
                        "Matrix remodeling","Wound response","Crypt base","Crypt top","Epithelial","Mesenchymal")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Signatures") 
p 
# ------------------------------------------------------------------------------
# Pathways
# ------------------------------------------------------------------------------
geneset_zdk=geneset[,c("FGF_ACT_CP","NOTCH_REACTOME","NFKB_BIOCARTA",
                       "MAPK_KEGG","PI3K_ACT_REACTOME","SRC_ACT_BIOCARTA","JAK_STAT_KEGG","CASPASE_BIOCARTA",
                       "PROTEASOME_KEGG","CELL_CYCLE_BIOCARTA","SHH_KEGG",
                       "INTEGRIN_BETA3_CP","VEGF_VEGFR_REACTOME","TRANSLATION_RIBOS_REACTOME")]
colnames(geneset_zdk)=c("FGF","NOTCH","NF-κB","MAPK","PI3K","SRC","JAK-STAT","Caspases","Proteosome",
                        "Cell cycle","SHH","Integrin-β3","VEGF,VEGFR","Translation ribosome")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Pathways")   
p 
# ------------------------------------------------------------------------------
# Immune
# ------------------------------------------------------------------------------
geneset_zdk=geneset[,c("IMMUNE_RESP_GO_BP","IMMUNE_MDSC_ALBELDA",
                       "IMMUNE_NKC_BREAST","IMMUNE_TH1_GALON","IMMUNE_TH17_GOUNARI",
                       "IMMUNE_THF_BREAST","IMMUNE_TREG_GALON","PD1_REACTOME",
                       "COMPLEMENT_COAG_KEGG")]
colnames(geneset_zdk)=c("Immune response","MDSC",
                        "NK cell infiltration","T(H)1 infiltration","T(H)17 activation",
                        "TF(H) infiltration","T(reg) activation","PD1 activation",
                        "Complement activation")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"), 
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Immune")  
p
# ------------------------------------------------------------------------------
# Metabolism
# ------------------------------------------------------------------------------
geneset_zdk=geneset[,c("ARACHNOID_METAB_KEGG","FATTY_ACID_METAB_KEGG","FRUTOSE_MANNOSE_METAB_KEGG",
                       "GALACTOSE_METAB_KEGG","GLUTATHIONE_KEGG","GLYCEROPHOSPHOLIPID_METAB_KEGG",
                       "LINOLEIC_METAB_KEGG","LYSOPHOSPHOLIPID_PID","NUCLEOTIDE_METAB_REACTOME",
                       "PENTOSE_GLUC_METAB_KEGG","STARCH_SUCROSE_METAB_KEGG","TYROSINE_METAB_KEGG",
                       "AMINO_SUGAR_NUCLEO_METAB_KEGG","ALANINE_ASPARTATE_GLUTAMATE_KEGG")]
colnames(geneset_zdk)=c("Arachnoid"," Fatty acid","Frutose Mannose",
                        "Galactose","Glutamine","Glycerophospholipid",
                        "Linoleic","Lysophospholipid","Nucleotide",
                        "Pentose","Starch sucrose","Tyrosine",
                        "Sugar,aa,nucleotide","Alainine,aspartate,glutamate")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"), 
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Metabolism") 
p
# ------------------------------------------------------------------------------
# WNT score calculation in syn26720761
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
load("exp_data/syn26720761_CMS.Rdata")
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=syn26720761_CMS,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
clinic=read.csv(file = "SG-BULK_patient_clinical_information.csv")
clinic$patient_id=paste0("X",clinic$patient_id)
WNT_Score_group1=merge(WNT_Score$final_activity_score,clinic,by.x="row.names",by.y="patient_id",all=T)
picture=c("group3")#The parameters ‘iCMS’, ‘CMS’, ‘group3’, ‘group5’ can be replaced.
if(picture=="CMS"){
  WNT_Score_group=WNT_Score_group1[,c("activity_score","CMS")]
  WNT_Score_group[is.na(WNT_Score_group)]="NOLBL"
  table(WNT_Score_group$CMS)
  WNT_Score_group$CMS=factor(WNT_Score_group$CMS,levels=c("CMS1","CMS2","CMS3","CMS4","NOLBL"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("CMS1","CMS2","CMS3","CMS4","NOLBL"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p = ggplot(WNT_Score_group, aes(x = CMS, y = activity_score, fill = CMS)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +  
    xlab("") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c("CMS1\n(n=21)",   
                                "CMS2\n(n=56)",  
                                "CMS3\n(n=24)",  
                                "CMS4\n(n=27)",  
                                "NOLBL\n(n=34)")) +  
    theme_bw() +  
    theme(  
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black"),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, color = "black")
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c('#f47920','#0073C2',"#ef5b9c",'#45b97c',"#a1a3a6")) +  
    ggtitle("syn26720761") 
  p
}else if(picture=="iCMS"){
  WNT_Score_group=WNT_Score_group1[,c("activity_score","iCMS")]
  table(WNT_Score_group$iCMS)
  WNT_Score_group$iCMS=factor(WNT_Score_group$iCMS,levels=c("iCMS2","iCMS3","indeterminate"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("iCMS2","iCMS3","indeterminate"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p = ggplot(WNT_Score_group, aes(x = iCMS, y = activity_score, fill = iCMS)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c("iCMS2\n(n=79)",   
                                "iCMS3\n(n=71)",  
                                "indeterminate\n(n=12)")) +  
    theme_bw() +  
    theme(  
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black"),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, color = "black")
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c('#8552a1','#e0861a',"#a1a3a6")) +  
    ggtitle("syn26720761") 
  p
}else if(picture=="group3"){
  WNT_Score_group=WNT_Score_group1[,c("activity_score","group3")]
  WNT_Score_group[is.na(WNT_Score_group)]="indeterminate"
  table(WNT_Score_group$group3)
  WNT_Score_group$group3=factor(WNT_Score_group$group3,levels=c("iCMS2_MSS","iCMS3_MSI","iCMS3_MSS","indeterminate"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("iCMS2_MSS","iCMS3_MSI","iCMS3_MSS","indeterminate"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p = ggplot(WNT_Score_group, aes(x = group3, y = activity_score, fill = group3)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c("iCMS2_MSS\n(n=75)",   
                                "iCMS3_MSI\n(n=30)",  
                                "iCMS3_MSS\n(n=41)",  
                                "indeterminate\n(n=16)")) +  
    theme_bw() +  
    theme(  
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black"),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, color = "black")
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c('#9b95c9','#fcaf17',"#c37e00","#a1a3a6")) +  
    ggtitle("syn26720761") 
  p
}else{
  WNT_Score_group=WNT_Score_group1[,c("activity_score","group5")]
  WNT_Score_group[is.na(WNT_Score_group)]="indeterminate"
  table(WNT_Score_group$group5)
  WNT_Score_group$group5=factor(WNT_Score_group$group5,levels=c("iCMS2_fibrotic","iCMS2_MSS","iCMS3_fibrotic","iCMS3_MSI","iCMS3_MSS","indeterminate"))
  library(ggplot2)
  library(ggsignif)
  library(ggpubr)
  library(RColorBrewer)
  combinations=combn(c("iCMS2_fibrotic","iCMS2_MSS","iCMS3_fibrotic","iCMS3_MSI","iCMS3_MSS","indeterminate"), 2)
  my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  p = ggplot(WNT_Score_group, aes(x = group5, y = activity_score, fill = group5)) +  
    geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
    geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
    xlab("") +  
    ylab("Wnt/β-catenin pathway activity score") +  
    guides(fill = "none") +  
    scale_x_discrete(labels = c("iCMS2_fibrotic\n(n=11)",   
                                "iCMS2_MSS\n(n=51)",  
                                "iCMS3_fibrotic\n(n=10)",  
                                "iCMS3_MSI\n(n=25)",  
                                "iCMS3_MSS\n(n=18)",  
                                "indeterminate\n(n=47)")) +  
    theme_bw() +  
    theme(  
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 16, color = "black"),
      plot.title = element_text(size = 16, color = "black"),
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 16, color = "black")
    ) +  
    stat_compare_means(comparisons = my.comparisons,   
                       method = "wilcox.test",  
                       label = "p.format",  
                       size = 4) +  
    scale_fill_manual(values = c('#8552a1','#9b95c9','#e0861a','#c37e00',"#fcaf17","#a1a3a6")) +  
    ggtitle("syn26720761")
  p
}
# ------------------------------------------------------------------------------
# syn26720761 Multi-clinical Analysis
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
load("exp_data/syn26720761_CMS.Rdata")
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=syn26720761_CMS,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
clinic=read.csv(file = "SG-BULK_patient_clinical_information.csv")
clinic$patient_id=paste0("X",clinic$patient_id)
WNT_Score_group=merge(WNT_Score$final_activity_score,clinic,by.x="row.names",by.y="patient_id",all=T)
WNT_Score_group0=WNT_Score_group[,c(2,5)]
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
p = ggscatter(WNT_Score_group0, x = "TMB", y = "activity_score",  
              add = "reg.line", 
              conf.int = TRUE,   
              add.params = list(color="#d71345", fill = "#BBBDBE"),  
              color = "#90d7ec", size = 3 ) +
  stat_cor(method = "pearson",  
           label.x = 0,   
           label.y = -4,  
           size = 6) +  
  theme_classic() +  
  theme(axis.title.x = element_text(size = 16, color = "black"),   
        axis.title.y = element_text(size = 16, color = "black"),  
        axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black")) +  
  ggtitle("syn26720761") +  
  xlab("TMB") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  theme(legend.position = "none")
p
WNT_Score_group0=WNT_Score_group[,c(2,6:18)]
WNT_Score_group0[,2:14]=ifelse(WNT_Score_group0[,2:14]=="wt","WT","MT")
WNT_Score_group0=reshape2::melt(WNT_Score_group0,id.vars = c("activity_score"))
colnames(WNT_Score_group0)=c("activity_score","gene","Group")
WNT_Score_group0$Group=factor(WNT_Score_group0$Group,levels=c("MT","WT"))
p = ggplot(WNT_Score_group0, aes(gene, activity_score, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +  # activity, activation, inhibition  
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
    axis.title.y = element_text(angle = 90, size = 16, colour = "black")
  ) +  
  scale_fill_manual(values = c("#C24976", "#469393")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") +  
  ggtitle("syn26720761", subtitle = "") +
  theme(plot.title = element_text(size = 16, colour = "black"))
p
WNT_Score_group0=WNT_Score_group[,c(2,20)]
WNT_Score_group0=na.omit(WNT_Score_group0)
table(WNT_Score_group0$MSI.Status)
colnames(WNT_Score_group0)=c("activity_score","msi")
WNT_Score_group0$msi=factor(WNT_Score_group0$msi,levels=c("MSI","MSS"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("MSI","MSS"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = msi, y = activity_score, fill = msi)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("MSI\n(n=33)",   
                              "MSS\n(n=126)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#d71345', '#145b7d')) +  
  ggtitle("syn26720761") 
p
WNT_Score_group0=WNT_Score_group[,c(2,21)]
WNT_Score_group0=na.omit(WNT_Score_group0)
table(WNT_Score_group0$CRIS)
colnames(WNT_Score_group0)=c("activity_score","CRIS")
WNT_Score_group0$msi=factor(WNT_Score_group0$CRIS,levels=c("CRIS-A","CRIS-B","CRIS-C","CRIS-D","CRIS-E"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("CRIS-A","CRIS-B","CRIS-C","CRIS-D","CRIS-E"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = CRIS, y = activity_score, fill = CRIS)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("CRIS-A\n(n=47)",   
                              "CRIS-B\n(n=23)",  
                              "CRIS-C\n(n=42)",  
                              "CRIS-D\n(n=23)",  
                              "CRIS-E\n(n=24)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#145b7d', '#8552a1', "#fcaf17", "#009ad6", "#007d65")) +  
  ggtitle("syn26720761")  
p
WNT_Score_group0=WNT_Score_group[,c(2,22)]
WNT_Score_group0=na.omit(WNT_Score_group0)
table(WNT_Score_group0$Gender)
colnames(WNT_Score_group0)=c("activity_score","Gender")
WNT_Score_group0$Gender=factor(WNT_Score_group0$Gender,levels=c("Male","Female"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("Male","Female"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = Gender, y = activity_score, fill = Gender)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("Male\n(n=65)",   
                              "Female\n(n=94)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#6EB9C3', '#C98B88')) +  
  ggtitle("syn26720761") 
p
WNT_Score_group0=WNT_Score_group[,c(2,23)]
WNT_Score_group0=na.omit(WNT_Score_group0)
colnames(WNT_Score_group0)=c("activity_score","age")
WNT_Score_group0$age=ifelse(WNT_Score_group0$age<=40,"age≤40",
                            ifelse(WNT_Score_group0$age<=60,"40<age≤60","age>60"))
WNT_Score_group0$age=factor(WNT_Score_group0$age,levels=c("age≤40","40<age≤60","age>60"))
table(WNT_Score_group0$age)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("age≤40","40<age≤60","age>60"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = age, y = activity_score, fill = age)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("age≤40\n(n=6)",   
                              "40<age≤60\n(n=49)",  
                              "age>60\n(n=104)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#90d7ec', '#009ad6', "#145b7d")) +  
  ggtitle("syn26720761") 
p
WNT_Score_group0=WNT_Score_group[,c(2,24)]
WNT_Score_group0=na.omit(WNT_Score_group0)
colnames(WNT_Score_group0)=c("activity_score","site")
table(WNT_Score_group0$site)
WNT_Score_group0=WNT_Score_group0[WNT_Score_group0$site!="Colon",]
WNT_Score_group0$site=ifelse(WNT_Score_group0$site%in%c("Cecum","Ascending colon","Hepatic flexure","Transverse colon"),"Right Colon",
                             ifelse(WNT_Score_group0$site%in%c("Splenic flexure","Descending colon","Sigmoid colon"),"Left Colon",
                                    ifelse(WNT_Score_group0$site=="Rectum","Rectum","Rectosigmoid Junction")))
WNT_Score_group0$site=factor(WNT_Score_group0$site,levels=c("Right Colon","Left Colon","Rectum","Rectosigmoid Junction"))
table(WNT_Score_group0$site)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("Right Colon","Left Colon","Rectum","Rectosigmoid Junction"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = site, y = activity_score, fill = site)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("Right Colon\n(n=51)",   
                              "Left Colon\n(n=65)",  
                              "Rectum\n(n=30)",  
                              "Rectosigmoid Junction\n(n=12)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#7fafc6', '#f7c67d', "#b4d9b9", "#9a82bb")) +  
  ggtitle("syn26720761")  
p
WNT_Score_group0=WNT_Score_group[,c(2,26)]
WNT_Score_group0=na.omit(WNT_Score_group0)
colnames(WNT_Score_group0)=c("activity_score","grade")
table(WNT_Score_group0$grade)
WNT_Score_group0=WNT_Score_group0[!WNT_Score_group0$grade%in%c("Unknown","Mucinous and Neuroendocrine carcinoma",
                                                               "Mucinous (colloid)","Mucinous","Adenosquamous carcinoma","2 to 3"),]
WNT_Score_group0$grade=ifelse(WNT_Score_group0$grade=="2 (mucinous)",2,WNT_Score_group0$grade)
WNT_Score_group0$grade=paste0("Grade",WNT_Score_group0$grade)
WNT_Score_group0$grade=factor(WNT_Score_group0$grade,levels=c("Grade1","Grade2","Grade3"))
table(WNT_Score_group0$grade)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("Grade1","Grade2","Grade3"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = grade, y = activity_score, fill = grade)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("Grade 1\n(n=8)",   
                              "Grade 2\n(n=129)",  
                              "Grade 3\n(n=10)")) +  
  theme_bw() +   
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#afb4db', '#9b95c9', "#6950a1")) +  
  ggtitle("syn26720761") 
p
WNT_Score_group0=WNT_Score_group[,c(2,28)]
WNT_Score_group0=na.omit(WNT_Score_group0)
colnames(WNT_Score_group0)=c("activity_score","stage")
table(WNT_Score_group0$stage)
WNT_Score_group0$stage=ifelse(WNT_Score_group0$stage%in%c("II","IIA","IIB","IIC"),2,
                              ifelse(WNT_Score_group0$stage%in%c("IIIA","IIIB","IIIC"),3,
                                     ifelse(WNT_Score_group0$stage%in%c("IV","IVA","IVB"),4,1)))
WNT_Score_group0$stage=paste0("Stage",WNT_Score_group0$stage)
WNT_Score_group0$stage=factor(WNT_Score_group0$stage,levels=c("Stage1","Stage2","Stage3","Stage4"))
table(WNT_Score_group0$stage)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("Stage1","Stage2","Stage3","Stage4"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group0, aes(x = stage, y = activity_score, fill = stage)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("Stage 1\n(n=12)",   
                              "Stage 2\n(n=48)",  
                              "Stage 3\n(n=50)",  
                              "Stage 4\n(n=49)")) +
  theme_bw() +   
  theme(  
    axis.text = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    axis.title.y = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 16, colour = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#abc88b', '#84bf96', "#45b97c", "#005831")) +  
  ggtitle("syn26720761")
p
# ------------------------------------------------------------------------------
# syn26720761 molecular characterisation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
load("exp_data/syn26720761_CMS.Rdata")
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                       expression_profile=syn26720761_CMS,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))

library(GSVA)
library(GSEABase)
library(reshape)
setwd("input")
GeneSets=getGmt("NM_genesets.gmt")
geneset=gsva(expr=as.matrix(syn26720761_CMS), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
geneset=as.data.frame(t(geneset))
geneset=as.data.frame(scale(geneset))
geneset=geneset[,-c(4,5)]
name=as.data.frame(colnames(geneset))
geneset_zdk1=geneset[,c("IMMUNE_ESTIMATE","STROMAL_ESTIMATE")]
colnames(geneset_zdk1)=c("Immune infiltration","Stromal infiltration")
geneset_zdk2=geneset[,c("SERRATED_UP","EMT_CORE_GENES","MYC_TARGETS_ZELLER","TGFB_CORE_GENES","CSC_BATLLE",
                        "MATRIX_REMODEL_REACTOME","WOUND_RESPONSE_GO_BP","CRYPT_BASE","CRYPT_TOP","EPITH_LOBODA","MESENCH_LOBODA")]
colnames(geneset_zdk2)=c("Serrated lesion","EMT activation","MYC targets","TGF-β activation","Cancer stem cell",
                         "Matrix remodeling","Wound response","Crypt base","Crypt top","Epithelial","Mesenchymal")
geneset_zdk3=geneset[,c("FGF_ACT_CP","NOTCH_REACTOME","NFKB_BIOCARTA",
                        "MAPK_KEGG","PI3K_ACT_REACTOME","SRC_ACT_BIOCARTA","JAK_STAT_KEGG","CASPASE_BIOCARTA",
                        "PROTEASOME_KEGG","CELL_CYCLE_BIOCARTA","SHH_KEGG",
                        "INTEGRIN_BETA3_CP","VEGF_VEGFR_REACTOME","TRANSLATION_RIBOS_REACTOME")]
colnames(geneset_zdk3)=c("FGF","NOTCH","NF-κB","MAPK","PI3K","SRC","JAK-STAT","Caspases","Proteosome",
                         "Cell cycle","SHH","Integrin-β3","VEGF,VEGFR","Translation ribosome")
geneset_zdk4=geneset[,c("IMMUNE_RESP_GO_BP","IMMUNE_MDSC_ALBELDA",
                        "IMMUNE_NKC_BREAST","IMMUNE_TH1_GALON","IMMUNE_TH17_GOUNARI",
                        "IMMUNE_THF_BREAST","IMMUNE_TREG_GALON","PD1_REACTOME",
                        "COMPLEMENT_COAG_KEGG")]
colnames(geneset_zdk4)=c("Immune response","MDSC",
                         "NK cell infiltration","T(H)1 infiltration","T(H)17 activation",
                         "TF(H) infiltration","T(reg) activation","PD1 activation",
                         "Complement activation")
geneset_zdk5=geneset[,c("ARACHNOID_METAB_KEGG","FATTY_ACID_METAB_KEGG","FRUTOSE_MANNOSE_METAB_KEGG",
                        "GALACTOSE_METAB_KEGG","GLUTATHIONE_KEGG","GLYCEROPHOSPHOLIPID_METAB_KEGG",
                        "LINOLEIC_METAB_KEGG","LYSOPHOSPHOLIPID_PID","NUCLEOTIDE_METAB_REACTOME",
                        "PENTOSE_GLUC_METAB_KEGG","STARCH_SUCROSE_METAB_KEGG","TYROSINE_METAB_KEGG",
                        "AMINO_SUGAR_NUCLEO_METAB_KEGG","ALANINE_ASPARTATE_GLUTAMATE_KEGG")]
colnames(geneset_zdk5)=c("Arachnoid"," Fatty acid","Frutose Mannose",
                         "Galactose","Glutamine","Glycerophospholipid",
                         "Linoleic","Lysophospholipid","Nucleotide",
                         "Pentose","Starch sucrose","Tyrosine",
                         "Sugar,aa,nucleotide","Alainine,aspartate,glutamate")
geneset_zdk=cbind(geneset_zdk1,geneset_zdk2)
geneset_zdk=cbind(geneset_zdk,geneset_zdk3)
geneset_zdk=cbind(geneset_zdk,geneset_zdk4)
geneset_zdk=cbind(geneset_zdk,geneset_zdk5)
fs=unique(colnames(geneset_zdk))
cancer_cor=data.frame()
for(i in 1:length(fs)){
  cor_test=cor.test(WNT_Score$final_activity_score$activity_score,geneset_zdk[,i],method="pearson")
  cancer_tme=data.frame(cancer=fs[i],tme="activity_score",cor=cor_test$estimate,p.value=cor_test$p.value)
  cancer_cor=rbind(cancer_cor,cancer_tme)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "   ", 
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*  ", 
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "** ", "***")))
cancer_cor$cancer=paste0(cancer_cor$pstar,cancer_cor$cancer)
cancer_cor$cancer <- factor(cancer_cor$cancer,levels = rev(cancer_cor$cancer))
library(ggplot2)
cancer_cor$relation <- ifelse(cancer_cor$cor > 0&cancer_cor$p.value<0.05,'pos',ifelse(cancer_cor$cor < 0&cancer_cor$p.value<0.05,"neg","none"))
ggplot(cancer_cor, aes(x=reorder(cancer, cor), y=cor, fill=relation)) +  
  geom_bar(stat='identity') +
  xlab('Correlation coefficient') +  
  ylab('WNT/β-catenin pathway activity score') +
  scale_fill_manual(values=c('#145b7d', "#d3d7d4", '#d71345')) +  
  theme_bw() +  
  theme(  
    panel.grid.major.y = element_blank(),  
    panel.grid.major.x = element_blank(),  
    legend.position = "none",  
    axis.text.x = element_text(angle = 75, hjust = 1)
  )  +
  ggtitle("syn26720761")
# ------------------------------------------------------------------------------
# All CRC datasets with molecular profiling
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("exp_data")
vector=c("GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE21815","GSE17537","GSE29621",
         "GSE38832","GSE72970","GSE33113",
         "UCSC.CRC_CMS2","GSE39582","GSE2109","GSE13067",
         "GSE13294","GSE14333","GSE17536",
         "GSE20916","GSE23878","GSE35896",
         "GSE37892","KFSYSCC","syn26720761",
         "syn2623706")
for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}
setwd("input")
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
cancer_cor_all=data.frame()
for (i in 1:length(vector)){
  rt=get(paste0(vector[i],"_CMS"))
  WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                         expression_profile=rt,
                                         foundation="relative_ssGSEA",
                                         activation_geneset=NA,
                                         inhibition_geneset=NA,
                                         geneSets_gmt=geneSets,
                                         min.sz=1,
                                         max.sz=10000,
                                         geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                             rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  library(GSVA)
  library(GSEABase)
  library(reshape)
  setwd("input")
  GeneSets=getGmt("NM_genesets.gmt")
  geneset=gsva(expr=as.matrix(rt), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  geneset=as.data.frame(scale(geneset))
  geneset=geneset[,-c(4,5)]
  name=as.data.frame(colnames(geneset))
  geneset_zdk1=geneset[,c("IMMUNE_ESTIMATE","STROMAL_ESTIMATE")]
  colnames(geneset_zdk1)=c("Immune infiltration","Stromal infiltration")
  geneset_zdk2=geneset[,c("SERRATED_UP","EMT_CORE_GENES","MYC_TARGETS_ZELLER","TGFB_CORE_GENES","CSC_BATLLE",
                          "MATRIX_REMODEL_REACTOME","WOUND_RESPONSE_GO_BP","CRYPT_BASE","CRYPT_TOP","EPITH_LOBODA","MESENCH_LOBODA")]
  colnames(geneset_zdk2)=c("Serrated lesion","EMT activation","MYC targets","TGF-β activation","Cancer stem cell",
                           "Matrix remodeling","Wound response","Crypt base","Crypt top","Epithelial","Mesenchymal")
  geneset_zdk3=geneset[,c("FGF_ACT_CP","NOTCH_REACTOME","NFKB_BIOCARTA",
                          "MAPK_KEGG","PI3K_ACT_REACTOME","SRC_ACT_BIOCARTA","JAK_STAT_KEGG","CASPASE_BIOCARTA",
                          "PROTEASOME_KEGG","CELL_CYCLE_BIOCARTA","SHH_KEGG",
                          "INTEGRIN_BETA3_CP","VEGF_VEGFR_REACTOME","TRANSLATION_RIBOS_REACTOME")]
  colnames(geneset_zdk3)=c("FGF","NOTCH","NF-κB","MAPK","PI3K","SRC","JAK-STAT","Caspases","Proteosome",
                           "Cell cycle","SHH","Integrin-β3","VEGF,VEGFR","Translation ribosome")
  geneset_zdk4=geneset[,c("IMMUNE_RESP_GO_BP","IMMUNE_MDSC_ALBELDA",
                          "IMMUNE_NKC_BREAST","IMMUNE_TH1_GALON","IMMUNE_TH17_GOUNARI",
                          "IMMUNE_THF_BREAST","IMMUNE_TREG_GALON","PD1_REACTOME",
                          "COMPLEMENT_COAG_KEGG")]
  colnames(geneset_zdk4)=c("Immune response","MDSC",
                           "NK cell infiltration","T(H)1 infiltration","T(H)17 activation",
                           "TF(H) infiltration","T(reg) activation","PD1 activation",
                           "Complement activation")
  geneset_zdk5=geneset[,c("ARACHNOID_METAB_KEGG","FATTY_ACID_METAB_KEGG","FRUTOSE_MANNOSE_METAB_KEGG",
                          "GALACTOSE_METAB_KEGG","GLUTATHIONE_KEGG","GLYCEROPHOSPHOLIPID_METAB_KEGG",
                          "LINOLEIC_METAB_KEGG","LYSOPHOSPHOLIPID_PID","NUCLEOTIDE_METAB_REACTOME",
                          "PENTOSE_GLUC_METAB_KEGG","STARCH_SUCROSE_METAB_KEGG","TYROSINE_METAB_KEGG",
                          "AMINO_SUGAR_NUCLEO_METAB_KEGG","ALANINE_ASPARTATE_GLUTAMATE_KEGG")]
  colnames(geneset_zdk5)=c("Arachnoid"," Fatty acid","Frutose Mannose",
                           "Galactose","Glutamine","Glycerophospholipid",
                           "Linoleic","Lysophospholipid","Nucleotide",
                           "Pentose","Starch sucrose","Tyrosine",
                           "Sugar,aa,nucleotide","Alainine,aspartate,glutamate")
  geneset_zdk=cbind(geneset_zdk1,geneset_zdk2)
  geneset_zdk=cbind(geneset_zdk,geneset_zdk3)
  geneset_zdk=cbind(geneset_zdk,geneset_zdk4)
  geneset_zdk=cbind(geneset_zdk,geneset_zdk5)
  fs=unique(colnames(geneset_zdk))
  cancer_cor=data.frame()
  for(j in 1:length(fs)){
    cor_test=cor.test(WNT_Score$final_activity_score$activity_score,geneset_zdk[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[j],tme="activity_score",cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_cor=rbind(cancer_cor,cancer_tme)
  }
  cancer_cor$cor[is.na(cancer_cor$cor)]=0
  cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
  cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "", 
                          ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*", 
                                 ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
  cancer_cor$dataset=vector[i]
  cancer_cor_all=rbind(cancer_cor_all,cancer_cor)
}
library(ggplot2)
library(dplyr)
cancer_cor_all$dataset=ifelse(cancer_cor_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",cancer_cor_all$dataset)
cancer_cor_all$cancer=factor(cancer_cor_all$cancer,levels = cancer_cor$cancer)
ggplot(cancer_cor_all, aes(dataset, cancer)) + 
  geom_tile(aes(fill = cor))+
  scale_fill_gradient2(high = "#f58f98",mid = "white",low = "#90d7ec")+
  geom_text(aes(label=pstar),col ="black",size =3)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))
# ------------------------------------------------------------------------------
# Data validation for robustness
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
load("exp_data/UCSC.CRC_CMS2_CMS.Rdata")
load("exp_data/UCSC.CRC_CMS2_CMS_G.Rdata")
group=UCSC.CRC_CMS2_CMS_G
co_sample=intersect(colnames(UCSC.CRC_CMS2_CMS),group$Tag)
UCSC.CRC_CMS2_CMS=UCSC.CRC_CMS2_CMS[,co_sample]
group=group[group$Tag%in%co_sample,]
exp=as.data.frame(UCSC.CRC_CMS2_CMS)
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
WNT_Score_group=merge(WNT_Score$final_activity_score,group,by.x="row.names",by.y="Tag",all=T)
table(WNT_Score_group$group)
WNT_Score_group$group=factor(WNT_Score_group$group,levels=c("CMS1","CMS2","CMS3","CMS4","NOLBL"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
combinations=combn(c("CMS1","CMS2","CMS3","CMS4","NOLBL"), 2)
my.comparisons=split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
p = ggplot(WNT_Score_group, aes(x = group, y = activity_score, fill = group)) +  
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +   
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +   
  xlab("") +  
  ylab("Wnt/β-catenin pathway activity score") +  
  guides(fill = "none") +  
  scale_x_discrete(labels = c("CMS1\n(n=39)",   
                              "CMS2\n(n=81)",  
                              "CMS3\n(n=61)",  
                              "CMS4\n(n=117)",  
                              "NOLBL\n(n=85)")) +  
  theme_bw() +  
  theme(  
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),  
    plot.title = element_text(size = 16, color = "black")
  ) +  
  stat_compare_means(comparisons = my.comparisons,   
                     method = "wilcox.test",  
                     label = "p.format",  
                     size = 4) +  
  scale_fill_manual(values = c('#f47920','#0073C2',"#ef5b9c",'#45b97c',"#a1a3a6")) +  
  ggtitle("UCSC.CRC")  
p
# ------------------------------------------------------------------------------
# robust test 1: randomly remove N dataset(s)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
vector_all=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
             "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
             "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
             "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
             "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q],".Rdata"))
  load(paste0(vector_all[q],"_G.Rdata"))
}
vector_marker = unique(sub("_.*", "", vector_all))
delete_number=c(1,2,3,4,5,6,7,8,9,10)
robustness_delete_all_50=list()
for (delete in 1:length(delete_number)){
  print(delete)
  robustness=data.frame()
  for (robust in 1:50){
    set.seed(Sys.time() + sample(1:10000,1))
    print(robust)
    random_element = sample(vector_marker, delete_number[delete])
    if(delete_number[delete]!=0){
      vector = vector_all[!grepl(paste0(random_element, collapse = "|"), vector_all)]
    }else{
      vector = vector_all
    }
    group_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
    source("functions/Get_denovo_genesets.R")
    up.denovo_genesets=Get_denovo_genesets(file_paths="input",
                                           expression_accession_vector=vector,
                                           group_HL=group_HL,
                                           gene_difference_method=c("limma"),
                                           alternative=c("greater"),
                                           p_combine_method=c("geometric_mean"),
                                           threshold_or_rank=c("rank"),
                                           top_genes=150,
                                           gene_Pfilter=NA,
                                           gene_FCfilter=NA,
                                           export_file=F)
    dn.denovo_genesets=Get_denovo_genesets(file_paths="input",
                                           expression_accession_vector=vector,
                                           group_HL=group_HL,
                                           gene_difference_method=c("limma"),
                                           alternative=c("less"),
                                           p_combine_method=c("geometric_mean"),
                                           threshold_or_rank=c("rank"),
                                           top_genes=150,
                                           gene_Pfilter=NA,
                                           gene_FCfilter=NA,
                                           export_file=F)
    WNT_denovo=cbind(up.denovo_genesets,dn.denovo_genesets)
    WNT_denovo = rbind(colnames(WNT_denovo), WNT_denovo)
    colnames(WNT_denovo) = paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")  
    WNT_knowledge = read.delim("WNT_knowledge.txt",header=TRUE,sep='\t',check.names=F)
    if(length(rownames(WNT_knowledge))>=length(rownames(WNT_denovo))){
      empty = data.frame(matrix("", nrow = (length(rownames(WNT_knowledge))-length(rownames(WNT_denovo))), ncol = length(unique(gsub("_.*", "",vector)))*2))
      colnames(empty)=colnames(WNT_denovo)
      WNT_denovo = rbind(WNT_denovo,empty)
      all_genesets=cbind(WNT_denovo,WNT_knowledge)
    }else{
      empty = data.frame(matrix("", nrow = (length(rownames(WNT_denovo))-length(rownames(WNT_knowledge))), ncol = length(colnames(WNT_knowledge))))
      colnames(empty)=colnames(WNT_knowledge)
      WNT_knowledge = rbind(WNT_knowledge,empty)
      all_genesets=cbind(WNT_denovo,WNT_knowledge)
    }
    ###########################################Add enrichr's gene set
    files=c("WNT_enrichr1.gmt","WNT_enrichr2.gmt","WNT_enrichr3.gmt","WNT_enrichr4.gmt")
    combined_list=list()
    for (file in files) {  
      gmt_lines=readLines(file)
      gene_sets=strsplit(gmt_lines, "\t")
      names(gene_sets)=sapply(gene_sets, function(x) paste0(x[1],"(",x[2],")","_enrichr"))
      for(i in 1:length(names(gene_sets))){
        gene_sets[[i]][1]=names(gene_sets)[i]
        gene_sets[[i]]=gene_sets[[i]][-2]
      }
      combined_list=append(combined_list, gene_sets)
    }  
    combined_list_names=as.data.frame(table(names(combined_list)))
    combined_list=combined_list[combined_list_names$Var1]
    for (i in 1:length(combined_list_names$Var1)) {  
      var_name=paste("WPRGS", i+152, "_knowledge", sep = "")  
      names(combined_list)[i]=var_name  
    }
    all_genesets=apply(all_genesets,2,function(x) x[nchar(x)>=1])
    all_genesets=append(all_genesets, combined_list)
    max_length=max(sapply(all_genesets, length))
    all_wnt_genesets=lapply(all_genesets, function(x) {
      if (length(x) < max_length) {
        c(x, rep("", max_length - length(x)))
      } else {
        x
      }
    })
    all_wnt_genesets=as.data.frame(all_wnt_genesets)
    save(all_wnt_genesets,file="all_wnt_genesets.Rdata")
    load("input/all_wnt_genesets.Rdata")
    all_genesets=all_wnt_genesets[-1,]
    UP.genesets=all_genesets
    DN.genesets=all_genesets
    source("functions/Precise_cleaning_system.R")
    group_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
    UP.clean_genesets_rank=Precise_cleaning_system(file_paths="input",
                                                   expression_accession_vector=vector,
                                                   group_HL=group_HL,
                                                   gene_difference_method=c("limma"),
                                                   alternative=c("greater"),
                                                   p_combine_method=c("geometric_mean"),
                                                   using_FC=F,
                                                   na_ratio=0.3,
                                                   using_KNN=F,
                                                   statistics=c("sumz"),
                                                   purpose=c("cleaned"),
                                                   threshold_type=c("rank"),
                                                   rank_threshold=150,
                                                   quantile_threshold=NA,
                                                   activation_geneset=UP.genesets,
                                                   inhibition_geneset=NA,
                                                   geneSets_gmt=NA,
                                                   min.sz=1,
                                                   max.sz=10000,
                                                   export_file=F)
    DN.clean_genesets_rank=Precise_cleaning_system(file_paths="input",
                                                   expression_accession_vector=vector,
                                                   group_HL=group_HL,
                                                   gene_difference_method=c("limma"),
                                                   alternative=c("less"),
                                                   p_combine_method=c("geometric_mean"),
                                                   using_FC=F,
                                                   na_ratio=0.3,
                                                   using_KNN=F,
                                                   statistics=c("sumz"),
                                                   purpose=c("cleaned"),
                                                   threshold_type=c("rank"),
                                                   rank_threshold=150,
                                                   quantile_threshold=NA,
                                                   activation_geneset=NA,
                                                   inhibition_geneset=DN.genesets,
                                                   geneSets_gmt=NA,
                                                   min.sz=1,
                                                   max.sz=10000,
                                                   export_file=F)
    UP.clean_genesets=UP.clean_genesets_rank$update_activation_geneset
    DN.clean_genesets=DN.clean_genesets_rank$update_inhibition_geneset
    source("functions/Joint_genesets.R")
    joint_genesets=Joint_genesets(file_paths="input",
                                  activation_geneset=UP.clean_genesets,
                                  inhibition_geneset=DN.clean_genesets,
                                  delete_GN=T,
                                  delete_time="after",
                                  integration_method=c("Jaccard"),
                                  min_GN=5,
                                  de_redundant_basis="igraph_component",
                                  similarity=0.5,
                                  export_file=F)
    activation_geneset=joint_genesets$joint_activation_geneset
    inhibition_geneset=joint_genesets$joint_inhibition_geneset
    if(nrow(activation_geneset)>nrow(inhibition_geneset)){
      empty=data.frame(matrix("", nrow=(length(rownames(activation_geneset))-length(rownames(inhibition_geneset))), ncol=length(colnames(inhibition_geneset))))
      colnames(empty)=colnames(inhibition_geneset)
      inhibition_geneset=rbind(inhibition_geneset,empty)
      rm(empty)
    }else if(nrow(activation_geneset)<nrow(inhibition_geneset)){
      empty=data.frame(matrix("", nrow=(length(rownames(inhibition_geneset))-length(rownames(activation_geneset))), ncol=length(colnames(activation_geneset))))
      colnames(empty)=colnames(activation_geneset)
      activation_geneset=rbind(activation_geneset,empty)
      rm(empty)
    }else{
      print("Same number of rows")
    }
    wpags_wpigs=cbind(activation_geneset,inhibition_geneset)
    wpags_wpigs=rbind(colnames(wpags_wpigs),wpags_wpigs)
    colnames(wpags_wpigs)=c(paste0("wpags", 1:length(colnames(activation_geneset))), 
                            paste0("wpigs", 1:length(colnames(inhibition_geneset))))
    wpags_wpigs=apply(wpags_wpigs,2,function(x) x[nchar(x)>=1])
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
    write.gmt(wpags_wpigs,file="input/wpags_wpigs_r.gmt")
    load("exp_data/UCSC.CRC_CMS2_CMS.Rdata")
    load("exp_data/UCSC.CRC_CMS2_CMS_G.Rdata")
    crc_data_all=UCSC.CRC_CMS2_CMS
    info_cms=UCSC.CRC_CMS2_CMS_G
    source("functions/Calculate_bioactivity_scores.R")
    geneSets=getGmt("wpags_wpigs_r.gmt")
    WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                           expression_profile=crc_data_all,
                                           foundation=c("relative_ssGSEA"),
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",length(colnames(activation_geneset))),
                                                               rep("neg",length(colnames(inhibition_geneset)))))
    WNT_Score_group=merge(WNT_Score$final_activity_score,info_cms,by.x="row.names",by.y="Tag",all=T)
    WNT_Score_group=na.omit(WNT_Score_group)
    WNT_Score_group=WNT_Score_group[WNT_Score_group$group%in%c("CMS1","CMS2","CMS3","CMS4"),]
    CMS1_score=median(WNT_Score_group[WNT_Score_group$group=="CMS1",]$activity_score)
    CMS1_score
    CMS2_score=median(WNT_Score_group[WNT_Score_group$group=="CMS2",]$activity_score)
    CMS2_score
    CMS3_score=median(WNT_Score_group[WNT_Score_group$group=="CMS3",]$activity_score)
    CMS3_score
    CMS4_score=median(WNT_Score_group[WNT_Score_group$group=="CMS4",]$activity_score)
    CMS4_score
    #CMS2>>CMS4/CMS3>>CMS1
    CMS23_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS2","CMS3"),]
    CMS23_score$group=factor(CMS23_score$group,levels = c("CMS2","CMS3"))
    CMS24_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS2","CMS4"),]
    CMS24_score$group=factor(CMS24_score$group,levels = c("CMS2","CMS4"))
    CMS34_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS3","CMS4"),]
    CMS34_score$group=factor(CMS34_score$group,levels = c("CMS3","CMS4"))
    CMS31_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS3","CMS1"),]
    CMS31_score$group=factor(CMS31_score$group,levels = c("CMS3","CMS1"))
    CMS41_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS4","CMS1"),]
    CMS41_score$group=factor(CMS41_score$group,levels = c("CMS4","CMS1"))
    CMS23_p=wilcox.test(activity_score~group, CMS23_score,alternative="greater")$p.value
    CMS24_p=wilcox.test(activity_score~group, CMS24_score,alternative="greater")$p.value
    CMS34_p=wilcox.test(activity_score~group, CMS34_score,alternative="two.sided")$p.value
    CMS31_p=wilcox.test(activity_score~group, CMS31_score,alternative="greater")$p.value
    CMS41_p=wilcox.test(activity_score~group, CMS41_score,alternative="greater")$p.value
    CMS_result=as.data.frame(c(CMS23_p,CMS24_p,CMS31_p,CMS41_p,CMS34_p))
    CMS_result=as.data.frame(t(CMS_result))
    colnames(CMS_result)=c("CMS23_p","CMS24_p","CMS31_p","CMS41_p","CMS34_p")
    rownames(CMS_result)=paste0("cycle_",robust)
    robustness=rbind(robustness,CMS_result)
  }
  robustness_delete_all_50[[delete]]=robustness
  names(robustness_delete_all_50)[delete]=paste0("delete_",delete)
}
save(robustness_delete_all_50,file="robustness_delete_all_50.Rdata")
# ------------------------------------------------------------------------------
# visualisation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("robustness_delete_all_50.Rdata")
data_correct_all=data.frame()
for(i in 1:length(names(robustness_delete_all_50))){
  data=robustness_delete_all_50[[i]]
  data$correctness=ifelse(data$CMS23_p<0.05&data$CMS24_p<0.05&data$CMS31_p<0.05&data$CMS41_p<0.05,"T","F")
  T_number=data[data$correct=="T",]
  data_correct=data.frame(Deletion=i,Correct=length(T_number$correct)/50)
  data_correct_all=rbind(data_correct_all,data_correct)
}
data_correct_all
data_correct_all$Deletion=as.character(data_correct_all$Deletion)
data_correct_all$Deletion=factor(data_correct_all$Deletion,levels = c("1","2","3","4","5",
                                                                      "6","7","8","9","10"))
library(ggplot2)
ggplot(data_correct_all, aes(x = Deletion, y = Correct, fill = Correct)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) +   
  geom_text(aes(label = Correct), vjust = -0.5, size = 5, color = "black") +  
  labs(title = "Robustness - Deletion Operation",  
       subtitle = "UCSC.CRC",  
       x = "Delete the number of training datasets",  
       y = "Correctness") +  
  theme_minimal(base_size = 16) + 
  theme(  
    plot.title = element_text(hjust = 0, size = 16, color = "black"),
    plot.subtitle = element_text(hjust = 0, size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),  
    legend.title = element_blank(),  
    legend.position = "none"  
  ) +  
  scale_fill_gradient(low = "#7bbfea", high = "#009ad6")  
# ------------------------------------------------------------------------------
# robust test 2: randomly swap N dataset(s)
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
vector_all=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
             "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
             "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
             "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
             "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q],".Rdata"))
  load(paste0(vector_all[q],"_G.Rdata"))
}
vector_marker = unique(sub("_.*", "", vector_all))
exchange_number=c(1,2,3,4,5,6,7,8,9,10)
robustness_exchange_all_50=list()
for (exchange in 1:length(exchange_number)){
  print(exchange)
  robustness=data.frame()
  for (robust in 1:50){
    print(robust)
    set.seed(Sys.time() + sample(1:10000,1))
    random_element = sample(vector_marker, exchange_number[exchange])
    vector = vector_all
    group_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
    HL = which(grepl(paste0(random_element, collapse = "|"), group_HL$Accession))
    temp = group_HL$Group1_Status[HL]
    group_HL$Group1_Status[HL] = group_HL$Group2_Status[HL]
    group_HL$Group2_Status[HL] = temp
    rm(temp)
    source("functions/Get_denovo_genesets.R")
    up.denovo_genesets=Get_denovo_genesets(file_paths="input",
                                           expression_accession_vector=vector,
                                           group_HL=group_HL,
                                           gene_difference_method=c("limma"),
                                           alternative=c("greater"),
                                           p_combine_method=c("geometric_mean"),
                                           threshold_or_rank=c("rank"),
                                           top_genes=150,
                                           gene_Pfilter=NA,
                                           gene_FCfilter=NA,
                                           export_file=F)
    dn.denovo_genesets=Get_denovo_genesets(file_paths="input",
                                           expression_accession_vector=vector,
                                           group_HL=group_HL,
                                           gene_difference_method=c("limma"),
                                           alternative=c("less"),
                                           p_combine_method=c("geometric_mean"),
                                           threshold_or_rank=c("rank"),
                                           top_genes=150,
                                           gene_Pfilter=NA,
                                           gene_FCfilter=NA,
                                           export_file=F)
    WNT_denovo=cbind(up.denovo_genesets,dn.denovo_genesets)
    WNT_denovo = rbind(colnames(WNT_denovo), WNT_denovo)
    colnames(WNT_denovo) = paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")  
    WNT_knowledge = read.delim("WNT_knowledge.txt",header=TRUE,sep='\t',check.names=F)
    if(length(rownames(WNT_knowledge))>=length(rownames(WNT_denovo))){
      empty = data.frame(matrix("", nrow = (length(rownames(WNT_knowledge))-length(rownames(WNT_denovo))), ncol = length(unique(gsub("_.*", "",vector)))*2))
      colnames(empty)=colnames(WNT_denovo)
      WNT_denovo = rbind(WNT_denovo,empty)
      all_genesets=cbind(WNT_denovo,WNT_knowledge)
    }else{
      empty = data.frame(matrix("", nrow = (length(rownames(WNT_denovo))-length(rownames(WNT_knowledge))), ncol = length(colnames(WNT_knowledge))))
      colnames(empty)=colnames(WNT_knowledge)
      WNT_knowledge = rbind(WNT_knowledge,empty)
      all_genesets=cbind(WNT_denovo,WNT_knowledge)
    }
    ###########################################Add enrichr's gene set
    files=c("WNT_enrichr1.gmt","WNT_enrichr2.gmt","WNT_enrichr3.gmt","WNT_enrichr4.gmt")
    combined_list=list()
    for (file in files) {  
      gmt_lines=readLines(file)
      gene_sets=strsplit(gmt_lines, "\t")
      names(gene_sets)=sapply(gene_sets, function(x) paste0(x[1],"(",x[2],")","_enrichr"))
      for(i in 1:length(names(gene_sets))){
        gene_sets[[i]][1]=names(gene_sets)[i]
        gene_sets[[i]]=gene_sets[[i]][-2]
      }
      combined_list=append(combined_list, gene_sets)
    }  
    combined_list_names=as.data.frame(table(names(combined_list)))
    combined_list=combined_list[combined_list_names$Var1]
    for (i in 1:length(combined_list_names$Var1)) {  
      var_name=paste("WPRGS", i+152, "_knowledge", sep = "")  
      names(combined_list)[i]=var_name  
    }
    all_genesets=apply(all_genesets,2,function(x) x[nchar(x)>=1])
    all_genesets=append(all_genesets, combined_list)
    max_length=max(sapply(all_genesets, length))
    all_wnt_genesets=lapply(all_genesets, function(x) {
      if (length(x) < max_length) {
        c(x, rep("", max_length - length(x)))
      } else {
        x
      }
    })
    all_wnt_genesets=as.data.frame(all_wnt_genesets)
    save(all_wnt_genesets,file="all_wnt_genesets.Rdata")
    load("input/all_wnt_genesets.Rdata")
    all_genesets=all_wnt_genesets[-1,]
    UP.genesets=all_genesets
    DN.genesets=all_genesets
    source("functions/Precise_cleaning_system.R")
    UP.clean_genesets_rank=Precise_cleaning_system(file_paths="input",
                                                   expression_accession_vector=vector,
                                                   group_HL=group_HL,
                                                   gene_difference_method=c("limma"),
                                                   alternative=c("greater"),
                                                   p_combine_method=c("geometric_mean"),
                                                   using_FC=F,
                                                   na_ratio=0.3,
                                                   using_KNN=F,
                                                   statistics=c("sumz"),
                                                   purpose=c("cleaned"),
                                                   threshold_type=c("rank"),
                                                   rank_threshold=150,
                                                   quantile_threshold=NA,
                                                   activation_geneset=UP.genesets,
                                                   inhibition_geneset=NA,
                                                   geneSets_gmt=NA,
                                                   min.sz=1,
                                                   max.sz=10000,
                                                   export_file=F)
    DN.clean_genesets_rank=Precise_cleaning_system(file_paths="input",
                                                   expression_accession_vector=vector,
                                                   group_HL=group_HL,
                                                   gene_difference_method=c("limma"),
                                                   alternative=c("less"),
                                                   p_combine_method=c("geometric_mean"),
                                                   using_FC=F,
                                                   na_ratio=0.3,
                                                   using_KNN=F,
                                                   statistics=c("sumz"),
                                                   purpose=c("cleaned"),
                                                   threshold_type=c("rank"),
                                                   rank_threshold=150,
                                                   quantile_threshold=NA,
                                                   activation_geneset=NA,
                                                   inhibition_geneset=DN.genesets,
                                                   geneSets_gmt=NA,
                                                   min.sz=1,
                                                   max.sz=10000,
                                                   export_file=F)
    UP.clean_genesets=UP.clean_genesets_rank$update_activation_geneset
    DN.clean_genesets=DN.clean_genesets_rank$update_inhibition_geneset
    source("functions/Joint_genesets.R")
    joint_genesets=Joint_genesets(file_paths="input",
                                  activation_geneset=UP.clean_genesets,
                                  inhibition_geneset=DN.clean_genesets,
                                  delete_GN=T,
                                  delete_time="after",
                                  integration_method=c("Jaccard"),
                                  min_GN=5,
                                  de_redundant_basis="igraph_component",
                                  similarity=0.5,
                                  export_file=F)
    activation_geneset=joint_genesets$joint_activation_geneset
    inhibition_geneset=joint_genesets$joint_inhibition_geneset
    if(nrow(activation_geneset)>nrow(inhibition_geneset)){
      empty=data.frame(matrix("", nrow=(length(rownames(activation_geneset))-length(rownames(inhibition_geneset))), ncol=length(colnames(inhibition_geneset))))
      colnames(empty)=colnames(inhibition_geneset)
      inhibition_geneset=rbind(inhibition_geneset,empty)
      rm(empty)
    }else if(nrow(activation_geneset)<nrow(inhibition_geneset)){
      empty=data.frame(matrix("", nrow=(length(rownames(inhibition_geneset))-length(rownames(activation_geneset))), ncol=length(colnames(activation_geneset))))
      colnames(empty)=colnames(activation_geneset)
      activation_geneset=rbind(activation_geneset,empty)
      rm(empty)
    }else{
      print("Same number of rows")
    }
    wpags_wpigs=cbind(activation_geneset,inhibition_geneset)
    wpags_wpigs=rbind(colnames(wpags_wpigs),wpags_wpigs)
    colnames(wpags_wpigs)=c(paste0("wpags", 1:length(colnames(activation_geneset))), 
                            paste0("wpigs", 1:length(colnames(inhibition_geneset))))
    wpags_wpigs=apply(wpags_wpigs,2,function(x) x[nchar(x)>=1])
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
    write.gmt(wpags_wpigs,file="input/wpags_wpigs_r.gmt")
    load("exp_data/UCSC.CRC_CMS2_CMS.Rdata")
    load("exp_data/UCSC.CRC_CMS2_CMS_G.Rdata")
    crc_data_all=UCSC.CRC_CMS2_CMS
    info_cms=UCSC.CRC_CMS2_CMS_G
    source("functions/Calculate_bioactivity_scores.R")
    geneSets=getGmt("wpags_wpigs_r.gmt")
    WNT_Score=Calculate_bioactivity_scores(file_paths="input",
                                           expression_profile=crc_data_all,
                                           foundation=c("relative_ssGSEA"),
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",length(colnames(activation_geneset))),
                                                               rep("neg",length(colnames(inhibition_geneset)))))
    WNT_Score_group=merge(WNT_Score$final_activity_score,info_cms,by.x="row.names",by.y="Tag",all=T)
    WNT_Score_group=na.omit(WNT_Score_group)
    WNT_Score_group=WNT_Score_group[WNT_Score_group$group%in%c("CMS1","CMS2","CMS3","CMS4"),]
    CMS1_score=median(WNT_Score_group[WNT_Score_group$group=="CMS1",]$activity_score)
    CMS1_score
    CMS2_score=median(WNT_Score_group[WNT_Score_group$group=="CMS2",]$activity_score)
    CMS2_score
    CMS3_score=median(WNT_Score_group[WNT_Score_group$group=="CMS3",]$activity_score)
    CMS3_score
    CMS4_score=median(WNT_Score_group[WNT_Score_group$group=="CMS4",]$activity_score)
    CMS4_score
    #CMS2>>CMS4/CMS3>>CMS1
    CMS23_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS2","CMS3"),]
    CMS23_score$group=factor(CMS23_score$group,levels = c("CMS2","CMS3"))
    CMS24_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS2","CMS4"),]
    CMS24_score$group=factor(CMS24_score$group,levels = c("CMS2","CMS4"))
    CMS34_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS3","CMS4"),]
    CMS34_score$group=factor(CMS34_score$group,levels = c("CMS3","CMS4"))
    CMS31_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS3","CMS1"),]
    CMS31_score$group=factor(CMS31_score$group,levels = c("CMS3","CMS1"))
    CMS41_score=WNT_Score_group[WNT_Score_group$group%in%c("CMS4","CMS1"),]
    CMS41_score$group=factor(CMS41_score$group,levels = c("CMS4","CMS1"))
    CMS23_p=wilcox.test(activity_score~group, CMS23_score,alternative="greater")$p.value
    CMS24_p=wilcox.test(activity_score~group, CMS24_score,alternative="greater")$p.value
    CMS34_p=wilcox.test(activity_score~group, CMS34_score,alternative="two.sided")$p.value
    CMS31_p=wilcox.test(activity_score~group, CMS31_score,alternative="greater")$p.value
    CMS41_p=wilcox.test(activity_score~group, CMS41_score,alternative="greater")$p.value
    CMS_result=as.data.frame(c(CMS23_p,CMS24_p,CMS31_p,CMS41_p,CMS34_p))
    CMS_result=as.data.frame(t(CMS_result))
    colnames(CMS_result)=c("CMS23_p","CMS24_p","CMS31_p","CMS41_p","CMS34_p")
    rownames(CMS_result)=paste0("cycle_",robust)
    robustness=rbind(robustness,CMS_result)
  }
  robustness_exchange_all_50[[exchange]]=robustness
  names(robustness_exchange_all_50)[exchange]=paste0("exchange_",exchange)
}
save(robustness_exchange_all_50,file="robustness_exchange_all_50.Rdata")
# ------------------------------------------------------------------------------
# visualisation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
load("robustness_exchange_all_50.Rdata")
data_correct_all=data.frame()
for(i in 1:length(names(robustness_exchange_all_50))){
  data=robustness_exchange_all_50[[i]]
  data$correctness=ifelse(data$CMS23_p<0.05&data$CMS24_p<0.05&data$CMS31_p<0.05&data$CMS41_p<0.05,"T","F")
  T_number=data[data$correct=="T",]
  data_correct=data.frame(Deletion=i,Correct=length(T_number$correct)/50)
  data_correct_all=rbind(data_correct_all,data_correct)
}
data_correct_all
data_correct_all$Deletion=as.character(data_correct_all$Deletion)
data_correct_all$Deletion=factor(data_correct_all$Deletion,levels = c("1","2","3","4","5",
                                                                      "6","7","8","9","10"))
library(ggplot2)
ggplot(data_correct_all, aes(x = Deletion, y = Correct, fill = Correct)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) +   
  geom_text(aes(label = Correct), vjust = -0.5, size = 5, color = "black") +  
  labs(title = "Robustness - Error Introduction",  
       subtitle = "UCSC.CRC",  
       x = "Number of erroneous training datasets",  
       y = "Correctness") +  
  theme_minimal(base_size = 16) +
  theme(  
    plot.title = element_text(hjust = 0, size = 16, color = "black"),
    plot.subtitle = element_text(hjust = 0, size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black"), 
    axis.title = element_text(size = 16, color = "black"),
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),  
    legend.title = element_blank(),  
    legend.position = "none"  
  ) +  
  scale_fill_gradient(low = "#7bbfea", high = "#009ad6") 