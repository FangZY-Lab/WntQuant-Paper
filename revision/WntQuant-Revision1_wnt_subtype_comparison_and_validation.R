# ==============================================================================
# WntQuant revision -- Signature comparison and subtype validation
# ------------------------------------------------------------------------------
# Purpose:    Compare the Wnt signature with published signatures in the
#             training, validation and large CRC cohorts, and characterise
#             the resulting subtypes.
# Inputs:     Processed CRC expression matrices, Wnt activity scores and
#             study metadata under input/.
# Outputs:    Sankey diagram, accuracy heatmaps, fGSEA direction results,
#             subtype proportion and prognosis figures.
# ==============================================================================

################################################################################sanky plot####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
sanky=read.delim("sankey.txt")
head(sanky)
if (!require("ggsankey")) install.packages("ggsankey")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(ggsankey)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
if (!require("ggalluvial")) install.packages("ggalluvial")
library(ggalluvial)

sanky_alluvial = sanky %>%
  count(GEO_Accession, type, sample_number, data_separation, genetics, test_cell, data_type)

nature_base_colors <- c("#8FC1E3", "#FFC8A2", "#96D496", "#F4A6A6", "#CAB3E6",
                              "#C8A89B", "#F1C6DE", "#BFBFBF", "#DEDE8D", "#8ED8E8",
                              "#D6E6F5", "#FFE0C2", "#CDF0C6", "#FFD4D3", "#E2D8EB",
                              "#D8C6BD", "#FADBEC", "#E3E3E3", "#EDEDC6", "#C7ECF2")
color_palette <- colorRampPalette(nature_base_colors)(length(unique(sanky_alluvial$GEO_Accession)))
p_alluvial = ggplot(sanky_alluvial,
                    aes(axis1 = type,
                        axis2 = GEO_Accession,
                        axis3 = sample_number,
                        axis4 = data_separation,
                        axis5 = genetics,
                        axis6 = test_cell,
                        axis7 = data_type,
                        y = n)) +
  geom_alluvium(aes(fill = GEO_Accession), alpha = 0.9,
                width = 1/4, knot.pos = 0.5,segments = 20) +
  geom_stratum(width = 1/3, alpha = 0.5, color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 5, color = "black") +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(limits = c("GEO_Accession", "type", "sample_number", "data_separation",
                              "genetics", "test_cell", "data_type"),
                   expand = c(0.1, 0.1)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white",color = NA),
        panel.background = element_rect(fill = "white",color = NA),
        plot.margin = margin(5, 8, 5, 8)) +  # Restore reasonable margins
  labs(x = NULL, y = NULL)
print(p_alluvial)
setwd("C:/Users/user/Desktop")
ggsave("sankey_plot.pdf",
       plot = p_alluvial,
       width = 22,
       height = 18,
       dpi = 300,
       bg = "white")
################################################################################Comparison with other signatures in the training datasets####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
vector=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
         "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
         "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
         "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
         "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}
for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  WNT_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
  rownames(WNT_HL)=WNT_HL$Accession
  WNT_HL=WNT_HL[,-1]
  WNT_HL_selected=WNT_HL[vector[i],,drop=F]
  group$group=ifelse(group$group==WNT_HL_selected$Group1,WNT_HL_selected$Group1_Status,WNT_HL_selected$Group2_Status)
  group$Tag <- paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}
all_group=all_group
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=1000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}
all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("C:/Users/user/Desktop/input")
  GeneSets=getGmt("wnt_usual_using.gmt")
  geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  usual=as.data.frame(scale(geneset))
  usual$ID=paste0(vector[cycle])
  all_usual=rbind(all_usual,usual)
}
all_data=merge(all_scores[,c(1,2,3),drop=F],all_usual,by="row.names")
all_data=merge(all_data,all_group,by.x="Row.names",by.y="Tag")
fs=unique(all_data$ID)
p_data_all=data.frame()
for(i in 1:length(fs)){
  all_data_fs=all_data[all_data$ID==fs[i],]
  gene_set=colnames(all_data_fs)[2:19]
  for(j in 1:length(gene_set)){
    all_data_fs_fs=all_data_fs[,c(gene_set[j],"group")]
    all_data_fs_fs$group=factor(all_data_fs_fs$group,levels=c("H","L"))
    colnames(all_data_fs_fs)=c("pathway","group")
    P.Value=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")$p.value
    FC=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group H"]]-
      t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group L"]]
    p_data=data.frame(dataset=fs[i],geneset=gene_set[j],p=P.Value,FC=FC)
    p_data_all=rbind(p_data_all,p_data)
  }
}
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
gene_set=sub("activity_score", "WntQuant", gene_set)
gene_set=sub("pos_score", "WntQuant_activation", gene_set)
gene_set=sub("neg_score", "WntQuant_inhibition", gene_set)
p_data_all$geneset=sub("activity_score", "WntQuant", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "WntQuant_activation", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "WntQuant_inhibition", p_data_all$geneset)
p_data_all$geneset=factor(p_data_all$geneset, levels = gene_set)
p = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-3, 3),
    oob = scales::squish
  ) +
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
    legend.title = element_text(size = 16, color = "black", angle = 270)
  ) +
  labs(fill = "Mean difference") +
  guides(fill = guide_colorbar(
    title.position = "right",
    label.position = "left"
  ))
setwd("C:/Users/user/Desktop/")
ggsave("training_benchmark.pdf",
       plot = p,
       width = 12,
       height = 7,
       device = "pdf")
################################################################################Comparison with other signatures in the validation datasets####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
vector=c("GSE102208","GSE114059_a","GSE114059_b",
         "GSE114059_c","GSE156082","GSE158578_a",
         "GSE158578_b","GSE188465","GSE188466",
         "GSE21815","GSE17537","GSE72970",
         "GSE188467","GSE29435","GSE29621",
         "GSE38832","GSE39904","GSE57760",
         "GSE64337_a","GSE64337_b","GSE94714",
         "GSE69687","GSE237130","UCSC.CRC_APC",
         "UCSC.CRC_TP53","UCSC.LIHC_CTNNB1","UCSC.LIHC_AXIN1",
         "UCSC.CRC_TumorNormal","GSE222023","GSE215947",
         "GSE178260","GSE135679","GSE60564",
         "GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE56896","GSE17385","GSE63072_a",
         "GSE63072_b","GSE61687","UCSC.CRC_CMS2",
         "GSE33113","GSE39582","GSE2109",
         "GSE13067","GSE13294","GSE14333",
         "GSE17536","GSE20916","GSE23878",
         "GSE35896","GSE37892","KFSYSCC",
         "syn26720761","syn2623706")
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}
for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  WNT_HL=read.delim("WNT_HL_V.txt",header=TRUE,sep='\t',check.names=F)
  rownames(WNT_HL)=WNT_HL$Accession
  WNT_HL=WNT_HL[,-1]
  WNT_HL_selected=WNT_HL[vector[i],,drop=F]
  group$group=ifelse(group$group==WNT_HL_selected$Group1,WNT_HL_selected$Group1_Status,WNT_HL_selected$Group2_Status)
  group$Tag <- paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}
all_group=all_group
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=10000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}
all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("C:/Users/user/Desktop/input")
  GeneSets=getGmt("wnt_usual_using.gmt")
  geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  usual=as.data.frame(scale(geneset))
  usual$ID=paste0(vector[cycle])
  all_usual=rbind(all_usual,usual)
}
all_data=merge(all_scores[,c(1,2,3),drop=F],all_usual,by="row.names")
all_data=merge(all_data,all_group,by.x="Row.names",by.y="Tag")
fs=unique(all_data$ID)
p_data_all=data.frame()
for(i in 1:length(fs)){
  all_data_fs=all_data[all_data$ID==fs[i],]
  gene_set=colnames(all_data_fs)[2:19]
  for(j in 1:length(gene_set)){
    all_data_fs_fs=all_data_fs[,c(gene_set[j],"group")]
    all_data_fs_fs$group=factor(all_data_fs_fs$group,levels=c("H","L"))
    colnames(all_data_fs_fs)=c("pathway","group")
    P.Value=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")$p.value
    FC=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group H"]]-
      t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group L"]]
    p_data=data.frame(dataset=fs[i],geneset=gene_set[j],p=P.Value,FC=FC)
    p_data_all=rbind(p_data_all,p_data)
  }
}
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
gene_set=sub("activity_score", "WntQuant", gene_set)
gene_set=sub("pos_score", "WntQuant_activation", gene_set)
gene_set=sub("neg_score", "WntQuant_inhibition", gene_set)
p_data_all$geneset=sub("activity_score", "WntQuant", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "WntQuant_activation", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "WntQuant_inhibition", p_data_all$geneset)
p_data_all$geneset=factor(p_data_all$geneset, levels = gene_set)
p = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-3, 3),
    oob = scales::squish
  ) +
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
    legend.title = element_text(size = 16, color = "black", angle = 270)
  ) +
  labs(fill = "Mean difference") +
  guides(fill = guide_colorbar(
    title.position = "right",
    label.position = "left"
  ))
setwd("C:/Users/user/Desktop/")
ggsave("testing_benchmark.pdf",
       plot = p,
       width = 25,
       height = 8,
       device = "pdf")
################################################################################fGSEA training direction#####
rm(list=ls())
gc()
library(GSEABase)
library(limma)
library(dplyr)
setwd("C:/Users/user/Desktop/input")
vector=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
         "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
         "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
         "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
         "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}
group_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  group=as.data.frame(group)
  Group_HL=group_HL
  rownames(Group_HL)=Group_HL$Accession
  Group_HL=Group_HL[,-1]
  group_HL_selected=Group_HL[paste0(vector[i]),,drop=F]
  group$group=ifelse(group$group==group_HL_selected$Group1,group_HL_selected$Group1_Status,group_HL_selected$Group2_Status)
  assign(paste0(vector[i],"_G"),group)
}
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
for (i in 1:length(vector)){
  result=Get_limma(gene=as.data.frame(get(vector[i])),
                   group=as.data.frame(get(paste0(vector[i],"_G"))),
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
cairo_pdf("GSEA_heatmap_training.pdf", width = 5.5, height = 5, family = "Arial")
draw(Heatmap(fgsea_all_transposed,
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
))
dev.off()
################################################################################fGSEA validation direction#####
rm(list=ls())
gc()
library(GSEABase)
library(limma)
library(dplyr)
setwd("C:/Users/user/Desktop/input")
vector=c("GSE102208","GSE114059_a","GSE114059_b",
         "GSE114059_c","GSE156082","GSE158578_a",
         "GSE158578_b","GSE188465","GSE188466",
         "GSE21815","GSE17537","GSE72970",
         "GSE188467","GSE29435","GSE29621",
         "GSE38832","GSE57760",
         "GSE64337_a","GSE64337_b","GSE94714",
         "GSE69687","GSE237130","UCSC.CRC_APC",
         "UCSC.CRC_TP53","UCSC.LIHC_CTNNB1","UCSC.LIHC_AXIN1",
         "UCSC.CRC_TumorNormal","GSE222023","GSE215947",
         "GSE178260","GSE135679","GSE60564",
         "GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE56896","GSE17385","GSE63072_a",
         "GSE63072_b","GSE61687","UCSC.CRC_CMS2",
         "GSE33113","GSE39582","GSE2109",
         "GSE13067","GSE13294","GSE14333",
         "GSE17536","GSE20916","GSE23878",
         "GSE35896","GSE37892","KFSYSCC",
         "syn26720761","syn2623706")
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}
group_HL=read.delim("WNT_HL_V.txt",header=TRUE,sep='\t',check.names=F)
for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  group=as.data.frame(group)
  Group_HL=group_HL
  rownames(Group_HL)=Group_HL$Accession
  Group_HL=Group_HL[,-1]
  group_HL_selected=Group_HL[paste0(vector[i]),,drop=F]
  group$group=ifelse(group$group==group_HL_selected$Group1,group_HL_selected$Group1_Status,group_HL_selected$Group2_Status)
  assign(paste0(vector[i],"_G"),group)
}
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
for (i in 1:length(vector)){
  result=Get_limma(gene=as.data.frame(get(vector[i])),
                   group=as.data.frame(get(paste0(vector[i],"_G"))),
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
                 nproc = 1)
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
setwd("C:/Users/user/Desktop/")
cairo_pdf("GSEA_heatmap_testing.pdf", width = 11.5, height = 6, family = "Arial")
draw(Heatmap(fgsea_all_transposed,
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
))
dev.off()
################################################################################Validation in large CRC cohorts (CMS)####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("C:/Users/user/Desktop/input/exp_data")
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
setwd("C:/Users/user/Desktop/input")
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
p_data_ee_all=data.frame()
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
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
  WNT_group=WNT_group[WNT_group$group!="NOLBL",]
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
  WNT_group=WNT_group[WNT_group$group!="NOLBL",]
  #####################CMS1 vs. Other
  WNT_group_CMS1=WNT_group
  WNT_group_CMS1$group=ifelse(WNT_group_CMS1$group=="CMS1","CMS1","Other")
  WNT_group_CMS1$group=factor(WNT_group_CMS1$group,levels=c("CMS1","Other"))
  ee=colnames(WNT_group_CMS1)[2:16]
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
  ee=colnames(WNT_group_CMS2)[2:16]
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
  ee=colnames(WNT_group_CMS3)[2:16]
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
  ee=colnames(WNT_group_CMS4)[2:16]
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
}

p_data_all=rbind(p_data_ee_all,p_data_ee_all_extra)
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","WntQuant",
                          ifelse(p_data_all$geneset=="pos_score","WntQuant_activation",
                                 ifelse(p_data_all$geneset=="neg_score","WntQuant_inhibition",
                                        p_data_all$geneset)))
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))
p_data_all$geneset=factor(p_data_all$geneset,levels = c("WntQuant",
                                                        "WntQuant_activation",
                                                        "WntQuant_inhibition",unique(p_data_ee_all_extra$geneset)))
p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
library(ggplot2)
library(dplyr)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)
p = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-1.5, 1.5),
    oob = scales::squish
  ) +
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
    legend.title = element_text(size = 16, color = "black", angle = 270)
  ) +
  labs(fill = "Mean difference") +
  guides(fill = guide_colorbar(
    title.position = "right",
    label.position = "left"
  ))
setwd("C:/Users/user/Desktop/")
ggsave("CRC_benchmark2.pdf",
       plot = p,
       width = 20,
       height = 6,
       device = "pdf")
################################################################################Validation in large CRC cohorts (iCMS)####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("C:/Users/user/Desktop/input/exp_data")
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
vector=vector1
for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}
setwd("C:/Users/user/Desktop/input")
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
icms2_icms3=getGmt("iCMS2_iCMS3.gmt")
p_data_ee_all=data.frame()
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
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
  ee=colnames(icms_wnt)[2:16]
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
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","WntQuant",
                          ifelse(p_data_all$geneset=="pos_score","WntQuant_activation",
                                 ifelse(p_data_all$geneset=="neg_score","WntQuant_inhibition",
                                        p_data_all$geneset)))
p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))
p_data_all$geneset=factor(p_data_all$geneset,levels = c("WntQuant",
                                                        "WntQuant_activation",
                                                        "WntQuant_inhibition",unique(p_data_ee_all_extra$geneset)))
p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
library(ggplot2)
library(dplyr)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)
p = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = Correlation), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-0.5, 0.5),
    oob = scales::squish
  ) +
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
    legend.title = element_text(size = 16, color = "black", angle = 270)
  ) +
  labs(fill = "Correlation") +
  guides(fill = guide_colorbar(
    title.position = "right",
    label.position = "left"
  ))
setwd("C:/Users/user/Desktop/")
ggsave("CRC_benchmark1_cor.pdf",
       plot = p,
       width = 21,
       height = 6.5,
       device = "pdf")
################################################################################Number of clusters for syn2623706 (k-means)####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/crc_data_all.Rdata")
exp=as.data.frame(crc_data_all)
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
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
#elbow method
fviz_nbclust(genesets, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  ggtitle("") +
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
  ggtitle("") +
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
p = ggplot(WNT, aes(x = cluster, y = activity_score, fill = cluster)) +
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("syn2623706_plot.pdf",
       plot = p,
       width = 4,
       height = 6,
       device = "pdf")
################################################################################Marker gene expression in the new subtypes####
rm(list=ls())
gc()
library(reshape2)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/crc_data_all.Rdata")
load("C:/Users/user/Desktop/input/cluster.Rdata")
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
  stat_compare_means(aes(group = group, label = ..p.signif..), size = 5, color = "black") + # Set size for comparison labels
  ggtitle("syn2623706")
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("syn2623706_markers.pdf",
       plot = p,
       width = 4,
       height = 6,
       device = "pdf")
################################################################################Subtype proportions and prognosis within Wnt activity strata####
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/cluster.Rdata")
load("C:/Users/user/Desktop/input/info_cms.Rdata")
load("C:/Users/user/Desktop/input/info_icms.Rdata")
cluster_group=merge(cluster,info_cms,by="Tag")
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
p=ggplot(HML, aes(x = 3,
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
    plot.title = element_text(size = 16, color = 'black'),
    legend.position = "top"
  ) +
  labs(fill = "", y = "", x = "")
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("CMS1_2_3_4_circle.pdf",
       plot = p,
       width = 8,
       height = 8,
       device = "pdf")
##################################################OS
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/cluster.Rdata")
setwd("C:/Users/user/Desktop/input")
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
p = ggsurvplot(
  fit = fit,
  data = survival_data_cms,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 1.2,
  xlim = c(0, 204),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("WSS1","WSS2","WSS3","WSS4"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "OS(n=1230)",
  ylim = c(0.2, 1),
  ggtheme = theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ),
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3
)
p.lab = paste0("log-rank test p",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))
p$plot = p$plot +
  annotate("text", x = 10, y = 0.25,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
##################################################RFS
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/cluster.Rdata")
setwd("C:/Users/user/Desktop/input")
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
p = ggsurvplot(
  fit = fit,
  data = survival_data_cms,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 1.2,
  xlim = c(0, 204),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("WSS1","WSS2","WSS3","WSS4"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "RFS(n=912)",
  ylim = c(0.38, 1),
  ggtheme = theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ),
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3
)
p.lab = paste0("log-rank test p",
               ifelse(p.val < 0.001, " < 0.001",
                      paste0(" = ", round(p.val, 3))))
p$plot = p$plot +
  annotate("text", x = 10, y = 0.42,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
###################################################SAR
rm(list=ls())
gc()
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/cluster.Rdata")
load("C:/Users/user/Desktop/input/OS_data.Rdata")
load("C:/Users/user/Desktop/input/RFS_data.Rdata")
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
p = ggsurvplot(
  fit = fit,
  data = survival_data_cms,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 1.2,
  xlim = c(0, 180),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("WSS1","WSS2","WSS3","WSS4"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "SAR(n=202)",
  ylim = c(0, 1),
  ggtheme = theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ),
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3
)
p.lab = paste0("log-rank test p",
               ifelse(p.val < 0.001, " < 0.001",
                      paste0(" = ", round(p.val, 3))))
p$plot = p$plot +
  annotate("text", x = 10, y = 0.04,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
################################################################################Name conversion####
rm(list=ls())
gc()
library(GSEABase)
GeneSets=getGmt("NM_genesets.gmt")
geneSets_gmt_self=lapply(GeneSets, function(x){ x@geneIds })
names(geneSets_gmt_self)=names(GeneSets)
geneSets_gmt=geneSets_gmt_self[c("IMMUNE_ESTIMATE","STROMAL_ESTIMATE",
                                 "SERRATED_UP","EMT_CORE_GENES","MYC_TARGETS_ZELLER","TGFB_CORE_GENES","CSC_BATLLE",
                                 "MATRIX_REMODEL_REACTOME","WOUND_RESPONSE_GO_BP","CRYPT_BASE","CRYPT_TOP","EPITH_LOBODA","MESENCH_LOBODA",
                                 "FGF_ACT_CP","NOTCH_REACTOME","NFKB_BIOCARTA",
                                 "MAPK_KEGG","PI3K_ACT_REACTOME","SRC_ACT_BIOCARTA","JAK_STAT_KEGG","CASPASE_BIOCARTA",
                                 "PROTEASOME_KEGG","CELL_CYCLE_BIOCARTA","SHH_KEGG",
                                 "INTEGRIN_BETA3_CP","VEGF_VEGFR_REACTOME","TRANSLATION_RIBOS_REACTOME",
                                 "IMMUNE_RESP_GO_BP","IMMUNE_MDSC_ALBELDA",
                                 "IMMUNE_NKC_BREAST","IMMUNE_TH1_GALON","IMMUNE_TH17_GOUNARI",
                                 "IMMUNE_THF_BREAST","IMMUNE_TREG_GALON","PD1_REACTOME",
                                 "COMPLEMENT_COAG_KEGG","ARACHNOID_METAB_KEGG","FATTY_ACID_METAB_KEGG","FRUTOSE_MANNOSE_METAB_KEGG",
                                 "GALACTOSE_METAB_KEGG","GLUTATHIONE_KEGG","GLYCEROPHOSPHOLIPID_METAB_KEGG",
                                 "LINOLEIC_METAB_KEGG","LYSOPHOSPHOLIPID_PID","NUCLEOTIDE_METAB_REACTOME",
                                 "PENTOSE_GLUC_METAB_KEGG","STARCH_SUCROSE_METAB_KEGG","TYROSINE_METAB_KEGG",
                                 "AMINO_SUGAR_NUCLEO_METAB_KEGG","ALANINE_ASPARTATE_GLUTAMATE_KEGG")]
for(i in 1:length(geneSets_gmt)) {
  geneSets_gmt[[i]] <- c(names(geneSets_gmt)[i], geneSets_gmt[[i]])
}
names(geneSets_gmt) = c("Immune infiltration", "Stromal infiltration",
                        "Serrated lesion", "EMT activation", "MYC targets", "TGF-β activation", "Cancer stem cell",
                        "Matrix remodeling", "Wound response", "Crypt base", "Crypt top", "Epithelial", "Mesenchymal",
                        "FGF", "NOTCH", "NF-κB", "MAPK", "PI3K", "SRC", "JAK-STAT", "Caspases", "Proteasome",
                        "Cell cycle", "SHH", "Integrin-β3", "VEGF,VEGFR", "Translation ribosome",
                        "Immune response", "MDSC",
                        "NK cell infiltration", "T helper 1 infiltration", "T helper 17 activation",
                        "T follicular helper infiltration", "Regulatory T cell activation", "PD1 activation",
                        "Complement activation",
                        "Arachidonic", "Fatty acid", "Fructose Mannose",
                        "Galactose", "Glutamine", "Glycerophospholipid",
                        "Linoleic", "Lysophospholipid", "Nucleotide",
                        "Pentose", "Starch sucrose", "Tyrosine",
                        "Sugar,AA,nucleotide", "Alanine,aspartate,glutamate")
write.gmt=function (gs.list, file) {
  if (length(grep(".gmt$", file)) == 0)
    stop("write.gmt(): file should be .gmt")
  gs.names=names(gs.list)
  gs.desc=sapply(gs.list, function(x) x[1]) # Use the first element of each list as gs.desc
  gs.lines=unlist(lapply(gs.list, function(x) {
    paste(x[-1], collapse = "\t") # Exclude the first element
  }))
  gs.lines=paste(gs.names, gs.desc, gs.lines, sep = "\t")
  writeLines(gs.lines, con = file)
}
write.gmt(geneSets_gmt,file="C:/Users/user/Desktop/input/all_NM_genesets.gmt")
################################################################################Molecular feature analysis####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
library(reshape)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/cluster.Rdata")
load("C:/Users/user/Desktop/input/crc_data_all.Rdata")
GeneSets=getGmt("all_NM_genesets.gmt")
geneset=gsva(expr=as.matrix(crc_data_all), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
geneset=as.data.frame(t(geneset))
geneset=as.data.frame(scale(geneset))
col=colnames(geneset)
data_p_all=data.frame()
for(i in 1:length(col)){
  geneset_col=geneset[,col[i],drop=F]
  colnames(geneset_col)="sig"
  geneset_col_cluster=merge(geneset_col,cluster,by.x="row.names",by.y="Tag")
  library(dplyr)
  median_summary = geneset_col_cluster %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      median_sig = median(sig, na.rm = TRUE),
      iqr_sig = IQR(sig, na.rm = TRUE)
    ) %>%
    arrange(cluster)
  print(median_summary)
  median_summary=as.data.frame(median_summary)
  kruskal_test = kruskal.test(sig ~ cluster, data = geneset_col_cluster)
  data_p=data.frame(sig=col[i],
                    WSS1=median_summary[1,3],
                    WSS2=median_summary[2,3],
                    WSS3=median_summary[3,3],
                    WSS4=median_summary[4,3],
                    p.value=kruskal_test$p.value)
  data_p_all=rbind(data_p_all,data_p)
}
library(ComplexHeatmap)
library(circlize)

categories <- list(
  "Infiltration ESTIMATE" = c("Immune infiltration", "Stromal infiltration"),
  "Signatures" = c("Serrated lesion", "EMT activation", "MYC targets", "TGF-β activation",
                   "Cancer stem cell", "Matrix remodeling", "Wound response", "Crypt base",
                   "Crypt top", "Epithelial", "Mesenchymal"),
  "Pathways" = c("FGF", "NOTCH", "NF-κB", "MAPK", "PI3K", "SRC", "JAK-STAT", "Caspases",
                 "Proteasome", "Cell cycle", "SHH", "Integrin-β3", "VEGF,VEGFR", "Translation ribosome"),
  "Immune" = c("Immune response", "MDSC", "NK cell infiltration", "T helper 1 infiltration",
               "T helper 17 activation", "T follicular helper infiltration", "Regulatory T cell activation",
               "PD1 activation", "Complement activation"),
  "Metabolism" = c("Arachidonic", "Fatty acid", "Fructose Mannose", "Galactose", "Glutamine",
                   "Glycerophospholipid", "Linoleic", "Lysophospholipid", "Nucleotide", "Pentose",
                   "Starch sucrose", "Tyrosine", "Sugar,AA,nucleotide", "Alanine,aspartate,glutamate")
)

category_colors <- c(
  "Infiltration ESTIMATE" = "#6baed6",
  "Signatures" = "#fdae61",
  "Pathways" = "#7bc87c",
  "Immune" = "#f391a9",
  "Metabolism" = "#c5b0d5"
)

sig_order <- c("Immune infiltration", "Stromal infiltration",
               "Serrated lesion", "EMT activation", "MYC targets", "TGF-β activation", "Cancer stem cell",
               "Matrix remodeling", "Wound response", "Crypt base", "Crypt top", "Epithelial", "Mesenchymal",
               "FGF", "NOTCH", "NF-κB", "MAPK", "PI3K", "SRC", "JAK-STAT", "Caspases", "Proteasome",
               "Cell cycle", "SHH", "Integrin-β3", "VEGF,VEGFR", "Translation ribosome",
               "Immune response", "MDSC", "NK cell infiltration", "T helper 1 infiltration", "T helper 17 activation",
               "T follicular helper infiltration", "Regulatory T cell activation", "PD1 activation", "Complement activation",
               "Arachidonic", "Fatty acid", "Fructose Mannose", "Galactose", "Glutamine", "Glycerophospholipid",
               "Linoleic", "Lysophospholipid", "Nucleotide", "Pentose", "Starch sucrose", "Tyrosine",
               "Sugar,AA,nucleotide", "Alanine,aspartate,glutamate")

data_p_all_ordered <- data_p_all[match(sig_order, data_p_all$sig), ]
heatmap_data <- as.matrix(data_p_all_ordered[, c("WSS1", "WSS2", "WSS3", "WSS4")])
heatmap_data_t <- t(heatmap_data)
colnames(heatmap_data_t) <- data_p_all_ordered$sig
p_value_data <- -log10(data_p_all_ordered$p.value)

category_vec <- character(ncol(heatmap_data_t))
for(cat_name in names(categories)) {
  category_vec[colnames(heatmap_data_t) %in% categories[[cat_name]]] <- cat_name
}

col_ha <- HeatmapAnnotation(
  Category = category_vec,
  col = list(Category = category_colors),
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16, col = "black"),  # Size 16, black
  annotation_legend_param = list(
    title = "Category",
    title_gp = gpar(fontsize = 16, col = "black"),  # Size 16, black
    labels_gp = gpar(fontsize = 16, col = "black")   # Size 16, black
  )
)

pval_ha <- HeatmapAnnotation(
  `-log[10](italic(p))` = anno_barplot(
    p_value_data,
    gp = gpar(fill = "#7d5886", col = NA),
    axis_param = list(gp = gpar(fontsize = 16, col = "black")),  # Axis text size 16, black
    height = unit(2, "cm"),
    ylim = c(0, max(p_value_data) * 1.1)
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "right",  # Move to the right
  annotation_name_gp = gpar(fontsize = 16, col = "black")  # Size 16, black
)

ht <- Heatmap(heatmap_data_t,
              name = "Median",
              col = colorRamp2(c(-1, 0, 1), c("#90d7ec", "white", "#f58f98")),
              top_annotation = col_ha,
              rect_gp = gpar(col = "black", lwd = 1),
              bottom_annotation = pval_ha,
              column_names_side = "bottom",
              column_names_rot = 60,
              column_names_gp = gpar(fontsize = 16, col = "black"),  # Size 16, black
              row_names_gp = gpar(fontsize = 16, col = "black"),     # Size 16, black
              row_names_rot = 0,
              column_title = "",
              row_title = "",
              show_column_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(
                title = "Median",
                title_gp = gpar(fontsize = 16, col = "black"),  # Size 16, black
                labels_gp = gpar(fontsize = 16, col = "black")   # Size 16, black
              ))

draw(ht,
     padding = unit(c(2, 2, 2, 2), "cm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right")

decorate_annotation("-log[10](italic(p))", {
  grid.lines(x = unit(c(0, 1), "npc"),
             y = unit(rep(-log10(1e-90), 2), "native"),
             gp = gpar(col = "red", lwd = 2, lty = 2))
  grid.text("p = 1×e-90",
            x = unit(1, "npc") + unit(2, "mm"),
            y = unit(-log10(1e-90), "native"),
            just = c("left", "center"),
            gp = gpar(col = "red", fontsize = 16))  # Font size 16
})
################################################################################TCGA pan-cancer classification matrix####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
genesets=as.matrix(pancancer_WNT_Score$single_score_matrix)
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
genesets=genesets[rownames(pancancer_group),]
library(mlbench)
library(factoextra)
library(ggplot2)
library (cluster)
library(fpc)
set.seed(123)
#0.8216041
#res=get_clust_tendency(genesets,300)
#res$hopkins_stat
#elbow method
fviz_nbclust(genesets, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  ggtitle("Optimal number of clusters") +
  labs(x = "Number of clusters", y = "Total within sum of square") +
  theme_minimal() +  # Apply minimal theme
  theme(
    text = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black")
  )
km_result = kmeans(genesets, centers=4,nstart = 1000,algorithm="Hartigan-Wong",iter.max=100)
# Judgement
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
  ggtitle("") +
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("pancancer_plot.pdf",
       plot = p,
       width = 6,
       height = 6,
       device = "pdf")
################################################################################Cumulative bar chart for TCGA####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
WNT=merge(pancancer_cluster,pancancer_group,by.x="Tag",by.y="row.names")
df_counts = as.data.frame(table(WNT$Project, WNT$cluster))
colnames(df_counts) = c("project", "cluster", "count")
df_percent = df_counts %>%
  group_by(project) %>%
  mutate(percent = count / sum(count) * 100)
p=ggplot(df_percent, aes(x = project, y = percent, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("#f7acbc", "#fedcbd", "#90d7ec", "#afb4db")) +
  labs(x = "", y = "Percentage (%)", fill = "Cluster") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("pancancer_percentage.pdf",
       plot = p,
       width = 10,
       height = 4,
       device = "pdf")
################################################################################OS####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
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
p = ggsurvplot(
  fit = fit,
  data = survival_data_cms,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 1.2,
  xlim = c(0, 376),
  break.time.by = 48,
  legend.title = "",
  legend.labs = c("panWSS1","panWSS2","panWSS3","panWSS4"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS(n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1),
  ggtheme = theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ),
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3
)
p.lab = paste0("log-rank test p",
               ifelse(p.val < 0.001, " < 0.001",
                      paste0(" = ", round(p.val, 3))))
p$plot = p$plot +
  annotate("text", x = 10, y = 0.05,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
################################################################################DSS####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
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
p = ggsurvplot(
  fit = fit,
  data = survival_data_cms,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 1.2,
  xlim = c(0, 376),
  break.time.by = 48,
  legend.title = "",
  legend.labs = c("panWSS1","panWSS2","panWSS3","panWSS4"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("DSS(n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1),
  ggtheme = theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ),
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3
)
p.lab = paste0("log-rank test p",
               ifelse(p.val < 0.001, " < 0.001",
                      paste0(" = ", round(p.val, 3))))
p$plot = p$plot +
  annotate("text", x = 10, y = 0.05,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
################################################################################PFI####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
library(dplyr)
library(ggplot2)
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
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
p = ggsurvplot(
  fit = fit,
  data = survival_data_cms,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 1.2,
  xlim = c(0, 376),
  break.time.by = 48,
  legend.title = "",
  legend.labs = c("panWSS1","panWSS2","panWSS3","panWSS4"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("PFI(n=",length(survival_data_cms$sample),")"),
  ylim = c(0, 1),
  ggtheme = theme_minimal() +
    theme(
      text = element_text(size = 16),
      legend.position = "top",
      panel.grid.minor = element_blank()
    ),
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3
)
p.lab = paste0("log-rank test p",
               ifelse(p.val < 0.001, " < 0.001",
                      paste0(" = ", round(p.val, 3))))
p$plot = p$plot +
  annotate("text", x = 10, y = 0.05,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
################################################################################Biological characterisation of the subtypes at the mRNA level####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
library(reshape)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
load("C:/Users/user/Desktop/input/hallmarker_activity.Rdata")
hallmarker_activity=as.data.frame(hallmarker_activity)
rownames(hallmarker_activity)=hallmarker_activity[,1]
hallmarker_activity=hallmarker_activity[,-1]
hallmarker_activity=as.data.frame(scale(hallmarker_activity))
hallmarker_activity = hallmarker_activity[, !names(hallmarker_activity) %in% "WNT_β_catenin_signaling"]
colnames(hallmarker_activity) = gsub("_", " ", colnames(hallmarker_activity))
col=colnames(hallmarker_activity)
data_p_all=data.frame()
for(i in 1:length(col)){
  geneset_col=hallmarker_activity[,col[i],drop=F]
  colnames(geneset_col)="sig"
  geneset_col_cluster=merge(geneset_col,pancancer_cluster,by.x="row.names",by.y="Tag")
  library(dplyr)
  median_summary = geneset_col_cluster %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      median_sig = median(sig, na.rm = TRUE),
      iqr_sig = IQR(sig, na.rm = TRUE)
    ) %>%
    arrange(cluster)
  print(median_summary)
  median_summary=as.data.frame(median_summary)
  kruskal_test = kruskal.test(sig ~ cluster, data = geneset_col_cluster)
  data_p=data.frame(sig=col[i],
                    panWSS1=median_summary[1,3],
                    panWSS2=median_summary[2,3],
                    panWSS3=median_summary[3,3],
                    panWSS4=median_summary[4,3],
                    p.value=kruskal_test$p.value)
  data_p_all=rbind(data_p_all,data_p)
}
replacements = c(
  "IL6 JAK STAT3 signaling" = "IL6/JAK/STAT3 signaling",
  "IFN α response" = "IFN-α response",
  "IFN γ response" = "IFN-γ response",
  "IL2 STAT5 signaling" = "IL-2/STAT5 signaling",
  "TNF α signaling via NFkB" = "TNF-α signaling via NF-κB",
  "PI3K AKT mTOR signaling" = "PI3K/AKT/mTOR signaling",
  "TGF β signaling" = "TGF-β signaling"
)
for (i in seq_along(replacements)) {
  data_p_all$sig = gsub(names(replacements)[i], replacements[i], data_p_all$sig)
}
data_p_all$p.value[data_p_all$p.value == 0] = 1e-323
library(ComplexHeatmap)
library(circlize)

category_vec <- c(
  "Adipogenesis" = "Development",
  "Allograft rejection" = "Immune",
  "Androgen response" = "Signaling",
  "Angiogenesis" = "Development",
  "Apical junction" = "Cellular component",
  "Apical surface" = "Cellular component",
  "Apoptosis" = "Pathway",
  "Bile acid metabolism" = "Metabolic",
  "Cholesterol homeostasis" = "Metabolic",
  "Coagulation" = "Immune",
  "Complement" = "Immune",
  "DNA repair" = "DNA damage",
  "E2F targets" = "Proliferation",
  "Epithelial mesenchymal transition" = "Development",
  "Estrogen response early" = "Signaling",
  "Estrogen response late" = "Signaling",
  "Fatty acid metabolism" = "Metabolic",
  "G2M checkpoint" = "Proliferation",
  "Glycolysis" = "Metabolic",
  "Hedgehog signaling" = "Signaling",
  "Heme metabolism" = "Metabolic",
  "Hypoxia" = "Pathway",
  "IL-2/STAT5 signaling" = "Signaling",
  "IL6/JAK/STAT3 signaling" = "Immune",
  "Inflammatory response" = "Immune",
  "IFN-α response" = "Immune",
  "IFN-γ response" = "Immune",
  "KRAS signaling down" = "Signaling",
  "KRAS signaling up" = "Signaling",
  "Mitosis spindle" = "Proliferation",
  "mTORC1 signaling" = "Signaling",
  "MYC targets version1" = "Proliferation",
  "MYC targets version2" = "Proliferation",
  "Myogenesis" = "Development",
  "NOTCH signaling" = "Signaling",
  "Oxidative phosphorylation" = "Metabolic",
  "P53 pathway" = "Proliferation",
  "Pancreatic β cells" = "Development",
  "Peroxisome" = "Cellular component",
  "PI3K/AKT/mTOR signaling" = "Signaling",
  "Protein secretion" = "Pathway",
  "Reactive oxygen species pathway" = "Pathway",
  "Spermatogenesis" = "Development",
  "TGF-β signaling" = "Signaling",
  "TNF-α signaling via NF-κB" = "Signaling",
  "Unfolded protein response" = "Pathway",
  "UV response down" = "DNA damage",
  "UV response up" = "DNA damage",
  "Xenobiotic metabolism" = "Metabolic"
)

category_colors <- c(
  "Cellular component" = "#66c2a5",
  "Development" = "#fc8d62",
  "DNA damage" = "#8da0cb",
  "Immune" = "#e78ac3",
  "Metabolic" = "#a6d854",
  "Pathway" = "#ffd92f",
  "Proliferation" = "#e5c494",
  "Signaling" = "#b3b3b3"
)

# Add a category column
data_p_all$category <- category_vec[data_p_all$sig]

# Sort by category and p-value
data_p_all_ordered <- data_p_all[order(data_p_all$category, -log10(data_p_all$p.value)), ]

# Prepare the heatmap data, laid out horizontally (columns are gene sets, rows are samples)
heatmap_data <- as.matrix(data_p_all_ordered[, c("panWSS1", "panWSS2", "panWSS3", "panWSS4")])
heatmap_data_t <- t(heatmap_data)
colnames(heatmap_data_t) <- data_p_all_ordered$sig
p_value_data <- -log10(data_p_all_ordered$p.value)

# Build the column annotation (categories)
col_ha <- HeatmapAnnotation(
  Category = data_p_all_ordered$category,
  col = list(Category = category_colors),
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16, col = "black"),
  annotation_legend_param = list(
    title = "Category",
    title_gp = gpar(fontsize = 16, col = "black"),
    labels_gp = gpar(fontsize = 16, col = "black")
  )
)

# Build the column annotation (p-values)
pval_ha <- HeatmapAnnotation(
  `-log[10](p)` = anno_barplot(
    p_value_data,
    gp = gpar(fill = "#7d5886", col = NA),
    axis_param = list(gp = gpar(fontsize = 16, col = "black")),
    height = unit(2, "cm"),
    ylim = c(0, max(p_value_data) * 1.1)
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 16, col = "black")
)

# Create the heatmap, laid out horizontally
ht <- Heatmap(heatmap_data_t,
              name = "Median",
              col = colorRamp2(c(-1, 0, 1), c("#90d7ec", "white", "#f58f98")),
              top_annotation = col_ha,
              bottom_annotation = pval_ha,
              rect_gp = gpar(col = "black", lwd = 1),
              column_names_side = "bottom",
              column_names_rot = 60,
              column_names_gp = gpar(fontsize = 16, col = "black"),
              row_names_gp = gpar(fontsize = 16, col = "black"),
              row_names_rot = 0,
              column_title = "",
              row_title = "",
              show_column_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(
                title = "Median",
                title_gp = gpar(fontsize = 16, col = "black"),
                labels_gp = gpar(fontsize = 16, col = "black")
              ))

# Draw the heatmap
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.5,
                      width = 0.9, height = 0.9,
                      just = "center"))

draw(ht,
     padding = unit(c(2, 2, 2, 2), "cm"),
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     newpage = FALSE)

# Add the p-value reference line
decorate_annotation("-log[10](p)", {
  grid.lines(x = unit(c(0, 1), "npc"),
             y = unit(rep(-log10(1e-90), 2), "native"),
             gp = gpar(col = "red", lwd = 2, lty = 2))
  grid.text("p = 1×e-90",
            x = unit(1, "npc") + unit(2, "mm"),
            y = unit(-log10(1e-90), "native"),
            just = c("left", "center"),
            gp = gpar(col = "red", fontsize = 16))
})
################################################################################TCGA-StromalScore_estimate####
rm(list=ls())
gc()
library(reshape2)
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
group=pancancer_cluster
colnames(group)=c("Tag","group")
load("C:/Users/user/Desktop/input/tme_all_estimate.Rdata")
tme_all=as.data.frame(tme_all_estimate)
library(ggplot2)
cor_test <- cor.test(tme_all$StromalScore_estimate,
                     tme_all$ImmuneScore_estimate)
ggplot(tme_all, aes(x = StromalScore_estimate, y = ImmuneScore_estimate)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  labs(title = "TCGA: StromalScore vs ImmuneScore",
       subtitle = paste0("Correlation: r = ", round(cor_test$estimate, 3),
                         ", p = ", format.pval(cor_test$p.value, digits = 3)),
       x = "Stromal Score Estimate",
       y = "Immune Score Estimate") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "darkred"))
biomarker_group=merge(tme_all,group,by.x="ID",by.y="Tag",all=T)
print(colnames(biomarker_group))
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("StromalScore_estimate","group")]
colnames(biomarker_group)=c("StromalScore_estimate","cluster")
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
p=ggplot(WNT, aes(x = cluster, y = StromalScore_estimate, fill = cluster)) +
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  xlab("") +
  ylab("Stromal score") +
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_SS.pdf",
       plot = p,
       width = 5.5,
       height = 6,
       device = "pdf")
################################################################################TCGA-ImmuneScore_estimate####
rm(list=ls())
gc()
library(reshape2)
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
group=pancancer_cluster
colnames(group)=c("Tag","group")
load("C:/Users/user/Desktop/input/tme_all_estimate.Rdata")
tme_all=as.data.frame(tme_all_estimate)
library(ggplot2)
cor_test <- cor.test(tme_all$StromalScore_estimate,
                     tme_all$ImmuneScore_estimate)
ggplot(tme_all, aes(x = StromalScore_estimate, y = ImmuneScore_estimate)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  labs(title = "TCGA: StromalScore vs ImmuneScore",
       subtitle = paste0("Correlation: r = ", round(cor_test$estimate, 3),
                         ", p = ", format.pval(cor_test$p.value, digits = 3)),
       x = "Stromal Score Estimate",
       y = "Immune Score Estimate") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "darkred"))
biomarker_group=merge(tme_all,group,by.x="ID",by.y="Tag",all=T)
print(colnames(biomarker_group))
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("ImmuneScore_estimate","group")]
colnames(biomarker_group)=c("ImmuneScore_estimate","cluster")
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
p=ggplot(WNT, aes(x = cluster, y = ImmuneScore_estimate, fill = cluster)) +
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  xlab("") +
  ylab("Immune score") +
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_IS.pdf",
       plot = p,
       width = 5.5,
       height = 6,
       device = "pdf")
################################################################################TCGA-ESTIMATEScore_estimate####
rm(list=ls())
gc()
library(reshape2)
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
group=pancancer_cluster
colnames(group)=c("Tag","group")
load("C:/Users/user/Desktop/input/tme_all_estimate.Rdata")
tme_all=as.data.frame(tme_all_estimate)
library(ggplot2)
cor_test <- cor.test(tme_all$StromalScore_estimate,
                     tme_all$ImmuneScore_estimate)
ggplot(tme_all, aes(x = StromalScore_estimate, y = ImmuneScore_estimate)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  labs(title = "TCGA: StromalScore vs ImmuneScore",
       subtitle = paste0("Correlation: r = ", round(cor_test$estimate, 3),
                         ", p = ", format.pval(cor_test$p.value, digits = 3)),
       x = "Stromal Score Estimate",
       y = "Immune Score Estimate") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "darkred"))
biomarker_group=merge(tme_all,group,by.x="ID",by.y="Tag",all=T)
print(colnames(biomarker_group))
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("ESTIMATEScore_estimate","group")]
colnames(biomarker_group)=c("ESTIMATEScore_estimate","cluster")
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
p=ggplot(WNT, aes(x = cluster, y = ESTIMATEScore_estimate, fill = cluster)) +
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  xlab("") +
  ylab("ESTIMATE score") +
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_ES.pdf",
       plot = p,
       width = 5.5,
       height = 6,
       device = "pdf")
################################################################################TCGA-TumorPurity_estimate####
rm(list=ls())
gc()
library(reshape2)
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
group=pancancer_cluster
colnames(group)=c("Tag","group")
load("C:/Users/user/Desktop/input/tme_all_estimate.Rdata")
tme_all=as.data.frame(tme_all_estimate)
library(ggplot2)
cor_test <- cor.test(tme_all$StromalScore_estimate,
                     tme_all$ImmuneScore_estimate)
ggplot(tme_all, aes(x = StromalScore_estimate, y = ImmuneScore_estimate)) +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "lightpink") +
  labs(title = "TCGA: StromalScore vs ImmuneScore",
       subtitle = paste0("Correlation: r = ", round(cor_test$estimate, 3),
                         ", p = ", format.pval(cor_test$p.value, digits = 3)),
       x = "Stromal Score Estimate",
       y = "Immune Score Estimate") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "darkred"))
biomarker_group=merge(tme_all,group,by.x="ID",by.y="Tag",all=T)
print(colnames(biomarker_group))
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("TumorPurity_estimate","group")]
colnames(biomarker_group)=c("TumorPurity_estimate","cluster")
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
p=ggplot(WNT, aes(x = cluster, y = TumorPurity_estimate, fill = cluster)) +
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  xlab("") +
  ylab("Tumor purity") +
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
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_TP.pdf",
       plot = p,
       width = 5.5,
       height = 6,
       device = "pdf")
################################################################################Comparison across TCGA immune subtypes####
rm(list = ls())
gc()

setwd("C:/Users/user/Desktop/input")

# Read the data
Immune_subtype <- read.table(
  "Subtype_Immune_Model_Based.txt",
  header = TRUE,
  sep = "\t",
  fill = FALSE
)

load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")

# Data merging and processing
immune_exp <- merge(
  pancancer_WNT_Score$final_activity_score,
  Immune_subtype,
  by.x = "row.names",
  by.y = "sample"
)

immune_exp <- immune_exp[, c("activity_score", "Subtype_Immune_Model_Based")]

table(immune_exp$Subtype_Immune_Model_Based)

immune_exp$Subtype_Immune_Model_Based <- factor(
  immune_exp$Subtype_Immune_Model_Based,
  levels = c(
    "Wound Healing (Immune C1)",
    "IFN-gamma Dominant (Immune C2)",
    "Inflammatory (Immune C3)",
    "Lymphocyte Depleted (Immune C4)",
    "Immunologically Quiet (Immune C5)",
    "TGF-beta Dominant (Immune C6)"
  )
)

# Load packages
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)

# Build the comparison pairs
combinations <- combn(
  c(
    "Wound Healing (Immune C1)",
    "IFN-gamma Dominant (Immune C2)",
    "Inflammatory (Immune C3)",
    "Lymphocyte Depleted (Immune C4)",
    "Immunologically Quiet (Immune C5)",
    "TGF-beta Dominant (Immune C6)"
  ),
  2
)

my.comparisons <- split(
  combinations,
  rep(1:ncol(combinations), each = nrow(combinations))
)

# Method 1: vertical layout with coord_flip()
p <- ggboxplot(
  immune_exp,
  x = "Subtype_Immune_Model_Based",
  y = "activity_score",
  fill = "Subtype_Immune_Model_Based",
  xlab = "",
  ylab = "Wnt/β-catenin pathway activity score",
  palette = c('#f58f98', '#b2d235', "#1d953f", "#90d7ec", "#009ad6", "#ef5b9c")
) +
  guides(fill = "none") +
  scale_x_discrete(
    labels = c(
      "Wound healing",
      "IFN-γ dominant",
      "Inflammatory",
      "Lymphocyte depleted",
      "Immunologically quiet",
      "TGF-β dominant"
    )
  ) +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  stat_compare_means(
    comparisons = my.comparisons,
    method = "wilcox.test",
    label = "p.signif",
    size = 4
  ) +
  ggtitle("")
print(p)
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_C123456.pdf",
       plot = p,
       width = 4,
       height = 7,
       device = "pdf")
################################################################################MSI####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
WNT_Score=pancancer_WNT_Score$final_activity_score
WNT_Score=WNT_Score[as.numeric(substr(rownames(WNT_Score), 14, 15)) <= 10, ]
WNT_Score=as.matrix(WNT_Score)
rownames(WNT_Score)=substr(rownames(WNT_Score), 1, 12)
library(limma)
WNT_Score=as.data.frame(avereps(WNT_Score))
#BiocManager::install("BiocOncoTK")
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
cancer_all = cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05,
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$relation <- ifelse(cancer_all$cor > 0 & cancer_all$p.value < 0.05, 'pos',
                              ifelse(cancer_all$cor < 0 & cancer_all$p.value < 0.05, "neg", "none"))
cancer_all$cancer_clean <- gsub("\\*", "", cancer_all$cancer)
cancer_all$cancer_clean <- factor(cancer_all$cancer_clean, levels = rev(cancer_all$cancer_clean))
p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer_clean, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer_clean), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) +
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  geom_text(aes(label = Signature,
                x = ifelse(cor > 0, -max(abs(cor)) * 0.1, max(abs(cor)) * 0.1),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("MSI") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_MSI.pdf",
       plot = p,
       width = 3.2,
       height = 12,
       device = "pdf")
################################################################################TMB####
rm(list=ls())
gc()
#BiocManager::install("PoisonAlien/TCGAmutations")
library(TCGAmutations)
tcga_available()
rm(list=ls())
setwd("C:/Users/user/Desktop/input/TMB")
fs <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
        "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV",
        "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
        "UCEC","UCS","UVM")
lapply(fs,function(x){
  library(TCGAmutations)
  print(x)
  maf <- TCGAmutations::tcga_load(study = x)
  save(maf, file = paste0("C:/Users/user/Desktop/input/TMB/TCGA_",x,"_MC3_maf.Rdata"))
})
rm(list=ls())
setwd("C:/Users/user/Desktop/input/TMB")
fs <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
        "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV",
        "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
        "UCEC","UCS","UVM")
TCGAtmb <- lapply(fs,function(x){
  library(maftools)
  load(paste0("C:/Users/user/Desktop/input/TMB/TCGA_",x,"_MC3_maf.Rdata"))
  TMB = data.frame(tmb(maf = maf))
  TMB$Tumor_Sample_Barcode <- as.character(TMB$Tumor_Sample_Barcode)
  TMB$Group = x
  print(x)
  return(TMB)
})
TCGAtmb <- do.call(rbind,TCGAtmb)
head(TCGAtmb)
TCGAtmb$Tumor_Sample_Barcode=substr(TCGAtmb$Tumor_Sample_Barcode, 1, 15)
TMB=TCGAtmb[,c(1,3,5)]
rownames(TMB)=TMB[,1]
TMB=TMB[,-1]
save(TMB,file="TMB.Rdata")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
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
cancer_all = cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05,
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$relation <- ifelse(cancer_all$cor > 0 & cancer_all$p.value < 0.05, 'pos',
                              ifelse(cancer_all$cor < 0 & cancer_all$p.value < 0.05, "neg", "none"))
cancer_all$cancer_clean <- gsub("\\*", "", cancer_all$cancer)
cancer_all$cancer_clean <- factor(cancer_all$cancer_clean, levels = rev(cancer_all$cancer_clean))
p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer_clean, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer_clean), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) +
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  geom_text(aes(label = Signature,
                x = ifelse(cor > 0, -max(abs(cor)) * 0.1, max(abs(cor)) * 0.1),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("TMB")+
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.15)),
    breaks = c(-0.15, 0, 0.15),  # Show only these three ticks
    labels = c("-0.15", "0", "0.15")  # Corresponding labels
  )
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_TMB.pdf",
       plot = p,
       width = 3.2,
       height = 12,
       device = "pdf")
################################################################################FGA####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
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
WNT_Score_FGA=merge(WNT_Score,pancancer_df,by.x="row.names",by.y = "sample")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","CRC","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
WNT_Score_FGA$friction_genome_altered=as.numeric(WNT_Score_FGA$friction_genome_altered)
cancer_all=data.frame()
for(i in 1:length(fs)){
  WNT_Score_fs=WNT_Score_FGA[WNT_Score_FGA$Group==fs[i],]
  cor_test=cor.test(WNT_Score_fs$activity_score,WNT_Score_fs$friction_genome_altered,method="pearson")
  cancer=data.frame(cancer=fs[i],cor=cor_test$estimate,p.value=cor_test$p.value)
  cancer_all=rbind(cancer_all,cancer)
}
library(ggsci)
library(ggplot2)
Color = pal_npg(alpha = 0.7)(6)
cancer_all = cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05,
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$relation <- ifelse(cancer_all$cor > 0 & cancer_all$p.value < 0.05, 'pos',
                              ifelse(cancer_all$cor < 0 & cancer_all$p.value < 0.05, "neg", "none"))
cancer_all$cancer_clean <- gsub("\\*", "", cancer_all$cancer)
cancer_all$cancer_clean <- factor(cancer_all$cancer_clean, levels = rev(cancer_all$cancer_clean))
p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer_clean, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer_clean), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) +
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  geom_text(aes(label = Signature,
                x = ifelse(cor > 0, -max(abs(cor)) * 0.1, max(abs(cor)) * 0.1),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("FGA") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_FGA.pdf",
       plot = p,
       width = 3.2,
       height = 12,
       device = "pdf")
################################################################################Mutation count####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
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
cancer_all = cancer_all[order(cancer_all$cor),]
cancer_all$Signature = ifelse(cancer_all$p.value < 0.05,
                              ifelse(cancer_all$p.value < 0.01,
                                     ifelse(cancer_all$p.value < 0.001,'***','**'),'*'),'')
cancer_all$relation <- ifelse(cancer_all$cor > 0 & cancer_all$p.value < 0.05, 'pos',
                              ifelse(cancer_all$cor < 0 & cancer_all$p.value < 0.05, "neg", "none"))
cancer_all$cancer_clean <- gsub("\\*", "", cancer_all$cancer)
cancer_all$cancer_clean <- factor(cancer_all$cancer_clean, levels = rev(cancer_all$cancer_clean))
p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer_clean, cor))) +
  ylab('') +
  xlab('Correlation coefficient') +
  geom_segment(aes(yend = cancer_clean), xend = 0, colour = 'grey50') +
  geom_point(size = 3.8, aes(colour = relation)) +
  scale_colour_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  geom_text(aes(label = Signature,
                x = ifelse(cor > 0, -max(abs(cor)) * 0.1, max(abs(cor)) * 0.1),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +
  ggtitle("Mutation count") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15)))
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("TCGA_MC.pdf",
       plot = p,
       width = 3.2,
       height = 12,
       device = "pdf")
################################################################################WNT-immune cell marker correlation analysis####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/tme_all_cellMarker.Rdata")
tme_all=as.data.frame(tme_all_cellMarker)
load("C:/Users/user/Desktop/input/TARGET_all_cellMarker.Rdata")
tme_all=rbind(tme_all,as.data.frame(TARGET_all_cellMarker))
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
tme_all_cellMarker=tme_all[,grep("_cellMarker", colnames(tme_all), value = TRUE)]
load("C:/Users/user/Desktop/input/pancancer_WNT_score.Rdata")
pancancer_wnt=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
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
cancer_cor$tme = recode(cancer_cor$tme,
                         "Activated_B_cell" = "Activated B cell",
                         "Activated_CD4_T_cell" = "Activated CD4+ T cell",
                         "Activated_CD8_T_cell" = "Activated CD8+ T cell",
                         "Activated_dendritic_cell" = "Activated dendritic cell",
                         "CD56bright_natural_killer_cell" = "CD56 bright natural killer cell",
                         "CD56dim_natural_killer_cell" = "CD56 dim natural killer cell",
                         "Central_memory_CD4_T_cell" = "Central memory CD4+ T cell",
                         "Central_memory_CD8_T_cell" = "Central memory CD8+ T cell",
                         "Effector_memeory_CD4_T_cell" = "Effector memory CD4+ T cell",
                         "Effector_memeory_CD8_T_cell" = "Effector memory CD8+ T cell",
                         "Gamma_delta_T_cell" = "γδ T cell",
                         "Immature__B_cell" = "Immature B cell",
                         "Immature_dendritic_cell" = "Immature dendritic cell",
                         "Mast_cell" = "Mast cell",
                         "Memory_B_cell" = "Memory B cell",
                         "Natural_killer_cell" = "Natural killer cell",
                         "Natural_killer_T_cell" = "Natural killer T cell",
                         "Plasmacytoid_dendritic_cell" = "Plasmacytoid dendritic cell",
                         "Regulatory_T_cell" = "Regulatory T cell",
                         "T_follicular_helper_cell" = "T follicular helper cell",
                         "Type_1_T_helper_cell" = "Type 1 T helper cell",
                         "Type_17_T_helper_cell" = "Type 17 T helper cell",
                         "Type_2_T_helper_cell" = "Type 2 T helper cell"
)
p=ggplot(cancer_cor, aes(cancer, tme)) +
  geom_tile(aes(fill = cor)) +
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 6) + # Text colour black
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black", face = "plain"), # Black, not bold
    axis.text.y = element_text(size = 16, color = "black", face = "plain"), # Black, not bold
    legend.text = element_text(size = 16, color = "black", face = "plain"), # Black, not bold
    legend.title = element_text(size = 16, color = "black", face = "plain") # Black, not bold
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation"))
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("immune_cellmarker.pdf",
       plot = p,
       width = 10,
       height = 10,
       device = "pdf")
################################################################################WNT-T cell state correlation analysis####
rm(list=ls())
gc()

# Load data
load("C:/Users/user/Desktop/input/tme_all_TCSS.Rdata")
load("C:/Users/user/Desktop/input/pancancer_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")

# Data preprocessing
T_status <- tme_all_TCSS
pancancer_wnt <- pancancer_WNT_Score$final_activity_score
pancancer_group <- pancancer_group[pancancer_group$Group == "Tumor",]

# Compute the correlation
cancer_cor <- data.frame()
for(i in unique(pancancer_group$Project)) {
  samples <- rownames(pancancer_group[pancancer_group$Project == i,])
  common_samples <- intersect(samples, rownames(T_status))

  if(length(common_samples) > 0) {
    tide_fs <- T_status[common_samples, ]
    wnt_fs <- pancancer_wnt[common_samples, "activity_score", drop = FALSE]

    for(j in colnames(tide_fs)) {
      cor_test <- cor.test(wnt_fs$activity_score, tide_fs[, j], method = "pearson")
      cancer_cor <- rbind(cancer_cor, data.frame(
        cancer = i, tme = j,
        cor = cor_test$estimate,
        p.value = cor_test$p.value
      ))
    }
  }
}

# Handle missing values
cancer_cor$cor[is.na(cancer_cor$cor)] <- 0
cancer_cor$p.value[is.na(cancer_cor$p.value)] <- 1

# Add the significance marks
cancer_cor$pstar <- ifelse(cancer_cor$p.value > 0.05, "",
                           ifelse(cancer_cor$p.value > 0.01, "*",
                                  ifelse(cancer_cor$p.value > 0.001, "**", "***")))

# Use reshape2::dcast
library(reshape2)

# Load reshape2 explicitly
m <- reshape2::dcast(cancer_cor, cancer ~ tme, value.var = "cor")
tmp <- reshape2::dcast(cancer_cor, cancer ~ tme, value.var = "pstar")
p <- reshape2::dcast(cancer_cor, cancer ~ tme, value.var = "p.value")

rownames(m) <- m[,1]; m <- m[,-1]
rownames(tmp) <- tmp[,1]; tmp <- tmp[,-1]
rownames(p) <- p[,1]; p <- p[,-1]

# Prepare the plotting data
library(tidyverse)

Plotdata <- data.frame(t(m)) %>%
  mutate(Term = rownames(.)) %>%
  pivot_longer(-Term, names_to = "name", values_to = "Correlation") %>%
  mutate(
    Pvalue = map2_dbl(Term, name, ~ {
      p_val <- t(p)[.x, .y]
      if(is.na(p_val)) 1 else p_val
    }),
    Pvalue = -log10(Pvalue),
    Pvalue = ifelse(Pvalue > 10, 10, Pvalue),
    Term = factor(Term, levels = c("Quiescence","Regulating","Proliferation","Helper","Cytotoxicity",
                                   "Progenitor_exhaustion","Terminal_exhaustion","Senescence"))
  )
Plotdata$Term = gsub("_", " ", Plotdata$Term)
# Plot
Dotplot2 <- ggplot(Plotdata, aes(x = name, y = Term, size = Pvalue)) +
  geom_point(shape = 21, aes(fill = Correlation)) +
  scale_size_continuous(range = c(0, 6), name = expression(-log[10](p)),
                        limits = c(0, 10), breaks = seq(0, 10, 2)) +
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",
                       name = "Correlation", limits = c(-0.85, 0.85), midpoint = 0) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 16),
        axis.text.y = element_text(colour = 'black', size = 16),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 16),
        legend.title = element_text(size = 16, color = "black")) +
  guides(fill = guide_colorbar(title.position = "left", title.hjust = 0.5,
                               barwidth = 15, barheight = 1.5))

Dotplot2
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("immune_T_status.pdf",
       plot = Dotplot2,
       width = 12,
       height = 4,
       device = "pdf")
################################################################################WNT-immune ESTIMATE correlation analysis####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/tme_all_estimate.Rdata")
tme_all=as.data.frame(tme_all_estimate)
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("C:/Users/user/Desktop/input/pancancer_WNT_score.Rdata")
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
#####################################Stromal
cancer_all=cancer_cor[cancer_cor$tme=="StromalScore_estimate",]
cancer_all$cancer <- factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
cancer_all$relation <- ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +
  geom_bar(stat = 'identity') +
  ylab('') +
  xlab("Correlation coefficient") +
  scale_fill_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  # Add significance asterisks
  geom_text(aes(label = pstar,
                x = ifelse(cor > 0, -0.02, 0.02),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  scale_x_continuous(breaks = c(-0.5, 0),
                     limits = c(-0.81, 0.2)) +
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
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("cor_SS.pdf",
       plot = p,
       width = 5,
       height = 10,
       device = "pdf")
#####################################Immune
cancer_all=cancer_cor[cancer_cor$tme=="ImmuneScore_estimate",]
cancer_all$cancer <- factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
library(ggplot2)
cancer_all$relation <- ifelse(cancer_all$cor > 0&cancer_all$p.value<0.05,'pos',ifelse(cancer_all$cor < 0&cancer_all$p.value<0.05,"neg","none"))
p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +
  geom_bar(stat = 'identity') +
  ylab('') +
  xlab("Correlation coefficient") +
  scale_fill_manual(values = c('#90d7ec', "#d3d7d4", '#f58f98')) +
  # Add significance asterisks
  geom_text(aes(label = pstar,
                x = ifelse(cor > 0, -0.02, 0.02),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  scale_x_continuous(breaks = c(-0.5, 0),
                     limits = c(-0.8, 0.1)) +
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
  ggtitle("Immune score")
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("cor_IS.pdf",
       plot = p,
       width = 5,
       height = 10,
       device = "pdf")
#####################################Tumor purity
cancer_all=cancer_cor[cancer_cor$tme=="TumorPurity_estimate",]
cancer_all$cancer <- factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
cancer_all$relation <- factor(ifelse(cancer_all$cor > 0 & cancer_all$p.value < 0.05, 'pos',
                                     ifelse(cancer_all$cor < 0 & cancer_all$p.value < 0.05, "neg", "none")),
                              levels = c("pos", "none", "neg"))  # Specify the level order explicitly

p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +
  geom_bar(stat = 'identity') +
  ylab('') +
  xlab("Correlation coefficient") +
  scale_fill_manual(values = c('neg' = '#90d7ec',  # Specify the colour mapping explicitly
                               'none' = "#d3d7d4",
                               'pos' = '#f58f98')) +
  geom_text(aes(label = pstar,
                x = ifelse(cor > 0, -0.02, 0.02),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 0.5),
                     limits = c(-0.1, 0.8)) +
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
  ggtitle("Tumor purity")
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("cor_TP.pdf",
       plot = p,
       width = 5,
       height = 10,
       device = "pdf")
#####################################ESTIMATE score
cancer_all=cancer_cor[cancer_cor$tme=="ESTIMATEScore_estimate",]
cancer_all$cancer <- factor(cancer_all$cancer,levels = rev(cancer_all$cancer))
cancer_all$relation <- factor(ifelse(cancer_all$cor > 0 & cancer_all$p.value < 0.05, 'pos',
                                     ifelse(cancer_all$cor < 0 & cancer_all$p.value < 0.05, "neg", "none")),
                              levels = c("pos", "none", "neg"))

p=ggplot(cancer_all, aes(x = cor, y = reorder(cancer, cor), fill = relation)) +
  geom_bar(stat = 'identity') +
  ylab('') +
  xlab("Correlation coefficient") +
  scale_fill_manual(values = c('neg' = '#90d7ec',
                               'none' = "#d3d7d4",
                               'pos' = '#f58f98')) +
  geom_text(aes(label = pstar,
                x = ifelse(cor > 0, -0.02, 0.02),
                hjust = ifelse(cor > 0, 1, 0)),
            size = 5, color = "black", fontface = "bold") +
  theme_bw() +
  scale_x_continuous(breaks = c(-0.5,0),
                     limits = c(-0.81, 0.1)) +
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
  ggtitle("ESTIMATE score")
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("cor_ES.pdf",
       plot = p,
       width = 5,
       height = 10,
       device = "pdf")
################################################################################Comparison with other signatures in the training datasets (with accuracy)####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
vector=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
         "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
         "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
         "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
         "GSE73906_f","GSE123807","GSE89124")

for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}

for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  WNT_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
  rownames(WNT_HL)=WNT_HL$Accession
  WNT_HL=WNT_HL[,-1]
  WNT_HL_selected=WNT_HL[vector[i],,drop=F]
  group$group=ifelse(group$group==WNT_HL_selected$Group1,WNT_HL_selected$Group1_Status,WNT_HL_selected$Group2_Status)
  group$Tag <- paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}

all_group=all_group
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=1000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}

all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("C:/Users/user/Desktop/input")
  GeneSets=getGmt("wnt_usual_using.gmt")
  geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  usual=as.data.frame(scale(geneset))
  usual$ID=paste0(vector[cycle])
  all_usual=rbind(all_usual,usual)
}

all_data=merge(all_scores[,c(1,2,3),drop=F],all_usual,by="row.names")
all_data=merge(all_data,all_group,by.x="Row.names",by.y="Tag")
fs=unique(all_data$ID)
p_data_all=data.frame()

for(i in 1:length(fs)){
  all_data_fs=all_data[all_data$ID==fs[i],]
  gene_set=colnames(all_data_fs)[2:19]
  for(j in 1:length(gene_set)){
    all_data_fs_fs=all_data_fs[,c(gene_set[j],"group")]
    all_data_fs_fs$group=factor(all_data_fs_fs$group,levels=c("H","L"))
    colnames(all_data_fs_fs)=c("pathway","group")
    P.Value=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")$p.value
    FC=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group H"]]-
      t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group L"]]
    p_data=data.frame(dataset=fs[i],geneset=gene_set[j],p=P.Value,FC=FC)
    p_data_all=rbind(p_data_all,p_data)
  }
}

p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
library(patchwork)

# Rename labels
gene_set=sub("activity_score", "WntQuant", gene_set)
gene_set=sub("pos_score", "WntQuant_activation", gene_set)
gene_set=sub("neg_score", "WntQuant_inhibition", gene_set)
p_data_all$geneset=sub("activity_score", "WntQuant", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "WntQuant_activation", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "WntQuant_inhibition", p_data_all$geneset)

# 1. Compute the accuracy of each signature (gene set)
p_data_all <- p_data_all %>%
  mutate(is_correct = case_when(
    geneset %in% c("FEVR_Wnt_Down", "WntQuant_inhibition") ~ (FC < 0 & p < 0.05),
    TRUE ~ (FC > 0 & p < 0.05)
  ))

accuracy_df <- p_data_all %>%
  group_by(geneset) %>%
  summarise(accuracy = mean(is_correct) * 100) # Compute the percentage accuracy

# 2. Sort the signatures by accuracy, high to low
geneset_order <- accuracy_df %>%
  arrange(accuracy) %>%
  pull(geneset)

p_data_all$geneset = factor(p_data_all$geneset, levels = geneset_order)
accuracy_df$geneset = factor(accuracy_df$geneset, levels = geneset_order)

# 3. Sort datasets alphabetically by name
sorted_datasets <- sort(unique(p_data_all$dataset))
p_data_all$dataset = factor(p_data_all$dataset, levels = sorted_datasets)

# 4. Draw the heatmap (datasets alphabetically on x, gene sets on y)
p1 = ggplot(p_data_all, aes(x = dataset, y = geneset)) +
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  geom_text(aes(label = pstar), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13.5, color = "black"),
    axis.text.y = element_text(size = 13.5, color = "black"),
    text = element_text(size = 13.5, color = "black"),
    legend.position = "bottom"
  ) +
  labs(fill = "Mean difference")

# 5. Draw the accuracy bar chart on the right (shares the gene sets with the heatmap, so they align perfectly)
p2 = ggplot(accuracy_df, aes(x = accuracy, y = geneset)) +
  geom_col(fill = "#9b95c9", color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", accuracy)), hjust = -0.2, size = 4, color = "black") +
  scale_x_continuous(limits = c(0, 118)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 13.5, color = "black"),
    text = element_text(size = 13.5, color = "black"),
    panel.grid.major.y = element_blank()
  ) +
  labs(x = "Accuracy (%)")

# 6. Combine with patchwork, placing the bars to the right of the heatmap
p = p1 + p2 + plot_layout(widths = c(4, 1))

# Save the figure
setwd("C:/Users/user/Desktop/")
ggsave("training_benchmark.pdf",
       plot = p,
       width = 12,
       height = 7,
       device = "pdf")
################################################################################Comparison with other signatures in the validation datasets (with accuracy)####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
vector=c("GSE102208","GSE114059_a","GSE114059_b",
         "GSE114059_c","GSE156082","GSE158578_a",
         "GSE158578_b","GSE188465","GSE188466",
         "GSE21815","GSE17537","GSE72970",
         "GSE188467","GSE29435","GSE29621",
         "GSE38832","GSE57760",
         "GSE64337_a","GSE64337_b","GSE94714",
         "GSE69687","GSE237130","UCSC.CRC_APC",
         "UCSC.CRC_TP53","UCSC.LIHC_CTNNB1","UCSC.LIHC_AXIN1",
         "UCSC.CRC_TumorNormal","GSE222023","GSE215947",
         "GSE178260","GSE135679","GSE60564",
         "GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE56896","GSE17385","GSE63072_a",
         "GSE63072_b","GSE61687","UCSC.CRC_CMS2",
         "GSE33113","GSE39582","GSE2109",
         "GSE13067","GSE13294","GSE14333",
         "GSE17536","GSE20916","GSE23878",
         "GSE35896","GSE37892","KFSYSCC",
         "syn26720761","syn2623706")

for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}

for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  WNT_HL=read.delim("WNT_HL_V.txt",header=TRUE,sep='\t',check.names=F)
  rownames(WNT_HL)=WNT_HL$Accession
  WNT_HL=WNT_HL[,-1]
  WNT_HL_selected=WNT_HL[vector[i],,drop=F]
  group$group=ifelse(group$group==WNT_HL_selected$Group1,WNT_HL_selected$Group1_Status,WNT_HL_selected$Group2_Status)
  group$Tag <- paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}

all_group=all_group
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=10000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}

all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("C:/Users/user/Desktop/input")
  GeneSets=getGmt("wnt_usual_using.gmt")
  geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  usual=as.data.frame(scale(geneset))
  usual$ID=paste0(vector[cycle])
  all_usual=rbind(all_usual,usual)
}

all_data=merge(all_scores[,c(1,2,3),drop=F],all_usual,by="row.names")
all_data=merge(all_data,all_group,by.x="Row.names",by.y="Tag")
fs=unique(all_data$ID)
p_data_all=data.frame()

for(i in 1:length(fs)){
  all_data_fs=all_data[all_data$ID==fs[i],]
  gene_set=colnames(all_data_fs)[2:19]
  for(j in 1:length(gene_set)){
    all_data_fs_fs=all_data_fs[,c(gene_set[j],"group")]
    all_data_fs_fs$group=factor(all_data_fs_fs$group,levels=c("H","L"))
    colnames(all_data_fs_fs)=c("pathway","group")
    P.Value=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")$p.value
    FC=t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group H"]]-
      t.test(pathway~group, all_data_fs_fs,alternative="two.sided")[["estimate"]][["mean in group L"]]
    p_data=data.frame(dataset=fs[i],geneset=gene_set[j],p=P.Value,FC=FC)
    p_data_all=rbind(p_data_all,p_data)
  }
}

p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
library(patchwork)

# Rename labels
gene_set=sub("activity_score", "WntQuant", gene_set)
gene_set=sub("pos_score", "WntQuant_activation", gene_set)
gene_set=sub("neg_score", "WntQuant_inhibition", gene_set)
p_data_all$geneset=sub("activity_score", "WntQuant", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "WntQuant_activation", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "WntQuant_inhibition", p_data_all$geneset)

# 1. Compute the accuracy of each signature (gene set)
p_data_all <- p_data_all %>%
  mutate(is_correct = case_when(
    geneset %in% c("FEVR_Wnt_Down", "WntQuant_inhibition") ~ (FC < 0 & p < 0.05),
    TRUE ~ (FC > 0 & p < 0.05)
  ))

accuracy_df <- p_data_all %>%
  group_by(geneset) %>%
  summarise(accuracy = mean(is_correct) * 100)

# 2. Sort the signatures by accuracy, high to low (y axis)
geneset_order <- accuracy_df %>%
  arrange(accuracy) %>%
  pull(geneset)

p_data_all$geneset = factor(p_data_all$geneset, levels = geneset_order)
accuracy_df$geneset = factor(accuracy_df$geneset, levels = geneset_order)

# 3. Sort datasets alphabetically by name
sorted_datasets <- sort(unique(p_data_all$dataset))
p_data_all$dataset = factor(p_data_all$dataset, levels = sorted_datasets)

# 4. Draw the heatmap (datasets alphabetically on x, gene sets by accuracy on y)
p1 = ggplot(p_data_all, aes(x = dataset, y = geneset)) +
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-3, 3),
    oob = scales::squish
  ) +
  geom_text(aes(label = pstar), col = "black", size = 4.5) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13.5, color = "black"),
    axis.text.y = element_text(size = 13.5, color = "black"),
    text = element_text(size = 13.5, color = "black"),
    legend.position = "bottom"
  ) +
  labs(fill = "Mean difference")

# 5. Draw the accuracy bar chart on the right
p2 = ggplot(accuracy_df, aes(x = accuracy, y = geneset)) +
  geom_col(fill = "#9b95c9", color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", accuracy)), hjust = -0.2, size = 4, color = "black") +
  scale_x_continuous(limits = c(0, 118)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 13.5, color = "black"),
    text = element_text(size = 13.5, color = "black"),
    panel.grid.major.y = element_blank()
  ) +
  labs(x = "Accuracy (%)")

# 6. Combine with patchwork, placing the bars to the right of the heatmap
p = p1 + p2 + plot_layout(widths = c(5, 1))

# Save the figure
setwd("C:/Users/user/Desktop/")
ggsave("testing_benchmark.pdf",
       plot = p,
       width = 25,
       height = 8,
       device = "pdf")
################################################################################Validation in large CRC cohorts (CMS, with correct counts)####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("C:/Users/user/Desktop/input/exp_data")

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
vector=vector1

for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}

setwd("C:/Users/user/Desktop/input")
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
p_data_ee_all=data.frame()

for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
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
  WNT_group=WNT_group[WNT_group$group!="NOLBL",]

  # CMS1 vs. Other
  WNT_group_CMS1=WNT_group
  WNT_group_CMS1$group=ifelse(WNT_group_CMS1$group=="CMS1","CMS1","Other")
  WNT_group_CMS1$group=factor(WNT_group_CMS1$group,levels=c("CMS1","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS1_single=WNT_group_CMS1[,c(ee[j],"group")]
    P.Value=wilcox.test(get(ee[j])~group, WNT_group_CMS1_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="CMS1",1])-
      mean(WNT_group_CMS1_single[WNT_group_CMS1_single$group=="Other",1])
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS1"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # CMS2 vs. Other
  WNT_group_CMS2=WNT_group
  WNT_group_CMS2$group=ifelse(WNT_group_CMS2$group=="CMS2","CMS2","Other")
  WNT_group_CMS2$group=factor(WNT_group_CMS2$group,levels=c("CMS2","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS2_single=WNT_group_CMS2[,c(ee[j],"group")]
    P.Value=wilcox.test(get(ee[j])~group, WNT_group_CMS2_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="CMS2",1])-
      mean(WNT_group_CMS2_single[WNT_group_CMS2_single$group=="Other",1])
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS2"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # CMS3 vs. Other
  WNT_group_CMS3=WNT_group
  WNT_group_CMS3$group=ifelse(WNT_group_CMS3$group=="CMS3","CMS3","Other")
  WNT_group_CMS3$group=factor(WNT_group_CMS3$group,levels=c("CMS3","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS3_single=WNT_group_CMS3[,c(ee[j],"group")]
    P.Value=wilcox.test(get(ee[j])~group, WNT_group_CMS3_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="CMS3",1])-
      mean(WNT_group_CMS3_single[WNT_group_CMS3_single$group=="Other",1])
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS3"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # CMS4 vs. Other
  WNT_group_CMS4=WNT_group
  WNT_group_CMS4$group=ifelse(WNT_group_CMS4$group=="CMS4","CMS4","Other")
  WNT_group_CMS4$group=factor(WNT_group_CMS4$group,levels=c("CMS4","Other"))
  ee=c("activity_score","pos_score","neg_score")
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_CMS4_single=WNT_group_CMS4[,c(ee[j],"group")]
    P.Value=wilcox.test(get(ee[j])~group, WNT_group_CMS4_single,alternative="two.sided")$p.value
    FC=mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="CMS4",1])-
      mean(WNT_group_CMS4_single[WNT_group_CMS4_single$group=="Other",1])
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="CMS4"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
}

geneSets=getGmt("wnt_usual_using.gmt")
p_data_ee_all_extra=data.frame()

for(i in 1:length(vector)){
  scores=gsva(expr=as.matrix(get(paste0(vector[i],"_CMS"))), geneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  scores=as.data.frame(t(scores))
  scores=as.data.frame(scale(scores))
  group=get(paste0(vector[i],"_CMS_G"))
  WNT_group=merge(scores,group,by.x="row.names",by.y="Tag")
  WNT_group=WNT_group[WNT_group$group!="NOLBL",]

  for(c_name in c("CMS1", "CMS2", "CMS3", "CMS4")) {
    WNT_group_c=WNT_group
    WNT_group_c$group=ifelse(WNT_group_c$group==c_name, c_name, "Other")
    WNT_group_c$group=factor(WNT_group_c$group, levels=c(c_name, "Other"))
    ee=colnames(WNT_group_c)[2:16]
    p_data_ee=data.frame()
    for(j in 1:length(ee)){
      WNT_group_single=WNT_group_c[,c(ee[j],"group")]
      colnames(WNT_group_single)=c("activity_score","group")
      P.Value=wilcox.test(activity_score~group, WNT_group_single,alternative="two.sided")$p.value
      FC=mean(WNT_group_single[WNT_group_single$group==c_name,]$activity_score)-
        mean(WNT_group_single[WNT_group_single$group=="Other",]$activity_score)
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=P.Value,FC=FC)
      p_data_ee=rbind(p_data_ee,p_data)
    }
    p_data_ee$contrast=c_name
    p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  }
}

p_data_all=rbind(p_data_ee_all,p_data_ee_all_extra)
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","WntQuant",
                          ifelse(p_data_all$geneset=="pos_score","WntQuant_activation",
                                 ifelse(p_data_all$geneset=="neg_score","WntQuant_inhibition",
                                        p_data_all$geneset)))

p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)

library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Updated correctness logic:
# - FEVR_Wnt_Down and WntQuant_inhibition: FC < 0 & p < 0.05 in CMS2, or FC > 0 & p < 0.05 in CMS1
# - Other signatures: FC > 0 & p < 0.05 in CMS2, or FC < 0 & p < 0.05 in CMS1
p_data_all <- p_data_all %>%
  mutate(is_correct = case_when(
    geneset %in% c("FEVR_Wnt_Down", "WntQuant_inhibition") ~ (
      (contrast == "CMS2" & FC < 0 & p < 0.05) |
        (contrast == "CMS1" & FC > 0 & p < 0.05)
    ),
    TRUE ~ (
      (contrast == "CMS2" & FC > 0 & p < 0.05) |
        (contrast == "CMS1" & FC < 0 & p < 0.05)
    )
  ))

# 2. Count the correct calls per signature and sort from most to fewest
correct_summary <- p_data_all %>%
  group_by(geneset, contrast) %>%
  summarise(correct_count = sum(is_correct), .groups = "drop")

total_correct <- correct_summary %>%
  group_by(geneset) %>%
  summarise(total = sum(correct_count), .groups = "drop") %>%
  arrange(total) # Sort by total ascending (ggplot2 draws bottom-up, so ascending puts the largest totals on top)

# Convert geneset to a factor ordered by the total number of correct calls
ordered_genesets <- total_correct$geneset
p_data_all$geneset = factor(p_data_all$geneset, levels = ordered_genesets)
correct_summary$geneset = factor(correct_summary$geneset, levels = ordered_genesets)

# 3. Draw the heatmap (p1)
p1 = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = FC), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-1.5, 1.5),
    oob = scales::squish
  ) +
  geom_text(aes(label = pstar), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  ) +
  labs(fill = "Mean difference")

# 4. Draw the stacked bar chart on the right (p2)
p2 = ggplot(correct_summary, aes(x = correct_count, y = geneset, fill = contrast)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(data = total_correct, aes(x = total, y = geneset, label = total),
            inherit.aes = FALSE, hjust = -0.2, size = 4.5, color = "black") +
  scale_x_continuous(limits = c(0, max(total_correct$total) * 1.25)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 14, color = "black"),
    text = element_text(size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Correct Count", fill = "Contrast")

# 5. Assemble and save
p = p1 + p2 + plot_layout(widths = c(4, 1.2))

setwd("C:/Users/user/Desktop/")
ggsave("CRC_benchmark1.pdf",
       plot = p,
       width = 24,
       height = 7,
       device = "pdf")
################################################################################Validation in large CRC cohorts (iCMS)####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("C:/Users/user/Desktop/input/exp_data")
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
vector=vector1
for (q in 1:length(vector)) {
  load(paste0(vector[q],"_CMS.Rdata"))
  load(paste0(vector[q],"_CMS_G.Rdata"))
}
setwd("C:/Users/user/Desktop/input")
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
icms2_icms3=getGmt("iCMS2_iCMS3.gmt")
p_data_ee_all=data.frame()

for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
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

  # iCMS2_Up
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Up")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Up, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS2_Up"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # iCMS2_Down
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Down")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Down, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS2_Down"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # iCMS3_Up
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Up")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Up, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS3_Up"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # iCMS3_Down
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Down")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Down, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS3_Down"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
}

p_data_ee_all_extra=data.frame()
geneSets=getGmt("wnt_usual_using.gmt")
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=gsva(expr=as.matrix(exp), geneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  scores=as.data.frame(t(scores))
  WNT=as.data.frame(scale(scores))
  icms=gsva(expr=as.matrix(exp), icms2_icms3, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  icms=as.data.frame(t(icms))
  icms=as.data.frame(scale(icms))
  icms_wnt=merge(WNT,icms,by="row.names")
  ee=colnames(icms_wnt)[2:16]

  for(c_name in c("iCMS2_Up", "iCMS2_Down", "iCMS3_Up", "iCMS3_Down")) {
    p_data_ee=data.frame()
    for(j in 1:length(ee)){
      WNT_group_single=icms_wnt[,c(ee[j], c_name)]
      res=cor.test(WNT_group_single[,1], WNT_group_single[,2], alternative="two.sided", method="pearson")
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
      p_data_ee=rbind(p_data_ee,p_data)
    }
    p_data_ee$contrast=c_name
    p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  }
}

p_data_all=rbind(p_data_ee_all,p_data_ee_all_extra)
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","WntQuant",
                          ifelse(p_data_all$geneset=="pos_score","WntQuant_activation",
                                 ifelse(p_data_all$geneset=="neg_score","WntQuant_inhibition",
                                        p_data_all$geneset)))

p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)

library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Updated iCMS correctness logic:
# - FEVR_Wnt_Down & WntQuant_inhibition：
#   negatively correlated with iCMS2_Up / iCMS3_Down (r < 0 & p < 0.05)
#   positively correlated with iCMS2_Down / iCMS3_Up (r > 0 & p < 0.05)
# - Other conventional signatures:
#   positively correlated with iCMS2_Up / iCMS3_Down (r > 0 & p < 0.05)
#   negatively correlated with iCMS2_Down / iCMS3_Up (r < 0 & p < 0.05)
p_data_all <- p_data_all %>%
  mutate(is_correct = case_when(
    geneset %in% c("FEVR_Wnt_Down", "WntQuant_inhibition") ~ (
      (contrast %in% c("iCMS2_Up", "iCMS3_Down") & Correlation < 0 & p < 0.05) |
        (contrast %in% c("iCMS2_Down", "iCMS3_Up") & Correlation > 0 & p < 0.05)
    ),
    TRUE ~ (
      (contrast %in% c("iCMS2_Up", "iCMS3_Down") & Correlation > 0 & p < 0.05) |
        (contrast %in% c("iCMS2_Down", "iCMS3_Up") & Correlation < 0 & p < 0.05)
    )
  ))

# 2. Count the correct calls per signature and sort from most to fewest
correct_summary <- p_data_all %>%
  group_by(geneset, contrast) %>%
  summarise(correct_count = sum(is_correct), .groups = "drop")

total_correct <- correct_summary %>%
  group_by(geneset) %>%
  summarise(total = sum(correct_count), .groups = "drop") %>%
  arrange(total) # Sort ascending so the largest totals end up on top

ordered_genesets <- total_correct$geneset
p_data_all$geneset = factor(p_data_all$geneset, levels = ordered_genesets)
correct_summary$geneset = factor(correct_summary$geneset, levels = ordered_genesets)

# 3. Draw the heatmap (p1)
p1 = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = Correlation), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-0.5, 0.5),
    oob = scales::squish
  ) +
  geom_text(aes(label = pstar), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  ) +
  labs(fill = "Correlation") +
  guides(fill = guide_colorbar(
    title.position = "right",
    label.position = "left"
  ))

# 4. Draw the stacked bar chart on the right (p2)
p2 = ggplot(correct_summary, aes(x = correct_count, y = geneset, fill = contrast)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(data = total_correct, aes(x = total, y = geneset, label = total),
            inherit.aes = FALSE, hjust = -0.2, size = 4.5, color = "black") +
  scale_x_continuous(limits = c(0, max(total_correct$total, na.rm = TRUE) * 1.3)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 14, color = "black"),
    text = element_text(size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Correct Count", fill = "Contrast")

# 5. Assemble and save
p = p1 + p2 + plot_layout(widths = c(4, 1.2))

setwd("C:/Users/user/Desktop/")
ggsave("CRC_benchmark1_cor.pdf",
       plot = p,
       width = 24,
       height = 7,
       device = "pdf")
################################################################################Validation in large CRC cohorts (iCMS)####
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
setwd("C:/Users/user/Desktop/input/exp_data")
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
setwd("C:/Users/user/Desktop/input")
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
icms2_icms3=getGmt("iCMS2_iCMS3.gmt")
p_data_ee_all=data.frame()

for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
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

  # iCMS2_Up
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Up")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Up, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS2_Up"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # iCMS2_Down
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS2_Down")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS2_Down, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS2_Down"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # iCMS3_Up
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Up")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Up, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS3_Up"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)

  # iCMS3_Down
  p_data_ee=data.frame()
  for(j in 1:length(ee)){
    WNT_group_single=icms_wnt[,c(ee[j],"iCMS3_Down")]
    res=cor.test(WNT_group_single[,1], WNT_group_single$iCMS3_Down, alternative="two.sided", method="pearson")
    p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
    p_data_ee=rbind(p_data_ee,p_data)
  }
  p_data_ee$contrast="iCMS3_Down"
  p_data_ee_all=rbind(p_data_ee_all,p_data_ee)
}

p_data_ee_all_extra=data.frame()
geneSets=getGmt("wnt_usual_using.gmt")
for(i in 1:length(vector)){
  exp=get(paste0(vector[i],"_CMS"))
  scores=gsva(expr=as.matrix(exp), geneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  scores=as.data.frame(t(scores))
  WNT=as.data.frame(scale(scores))
  icms=gsva(expr=as.matrix(exp), icms2_icms3, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  icms=as.data.frame(t(icms))
  icms=as.data.frame(scale(icms))
  icms_wnt=merge(WNT,icms,by="row.names")
  ee=colnames(icms_wnt)[2:16]

  for(c_name in c("iCMS2_Up", "iCMS2_Down", "iCMS3_Up", "iCMS3_Down")) {
    p_data_ee=data.frame()
    for(j in 1:length(ee)){
      WNT_group_single=icms_wnt[,c(ee[j], c_name)]
      res=cor.test(WNT_group_single[,1], WNT_group_single[,2], alternative="two.sided", method="pearson")
      p_data=data.frame(dataset=vector[i],geneset=ee[j],p=res$p.value,Correlation=res$estimate)
      p_data_ee=rbind(p_data_ee,p_data)
    }
    p_data_ee$contrast=c_name
    p_data_ee_all_extra=rbind(p_data_ee_all_extra,p_data_ee)
  }
}

p_data_all=rbind(p_data_ee_all,p_data_ee_all_extra)
p_data_all$geneset=ifelse(p_data_all$geneset=="activity_score","WntQuant",
                          ifelse(p_data_all$geneset=="pos_score","WntQuant_activation",
                                 ifelse(p_data_all$geneset=="neg_score","WntQuant_inhibition",
                                        p_data_all$geneset)))

p_data_all$pstar=ifelse(p_data_all$p > 0.05, "",
                        ifelse(p_data_all$p <= 0.05 & p_data_all$p > 0.01, "*",
                               ifelse(p_data_all$p <= 0.01 & p_data_all$p > 0.001, "**", "***")))

p_data_all$dataset=ifelse(p_data_all$dataset=="UCSC.CRC_CMS2","UCSC.CRC",p_data_all$dataset)
p_data_all$dataset=paste0(p_data_all$dataset,"-",p_data_all$contrast)

library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Updated iCMS correctness logic:
# - FEVR_Wnt_Down & WntQuant_inhibition：
#   negatively correlated with iCMS2_Up / iCMS3_Down (r < 0 & p < 0.05)
#   positively correlated with iCMS2_Down / iCMS3_Up (r > 0 & p < 0.05)
# - Other conventional signatures:
#   positively correlated with iCMS2_Up / iCMS3_Down (r > 0 & p < 0.05)
#   negatively correlated with iCMS2_Down / iCMS3_Up (r < 0 & p < 0.05)
p_data_all <- p_data_all %>%
  mutate(is_correct = case_when(
    geneset %in% c("FEVR_Wnt_Down", "WntQuant_inhibition") ~ (
      (contrast %in% c("iCMS2_Up", "iCMS3_Down") & Correlation < 0 & p < 0.05) |
        (contrast %in% c("iCMS2_Down", "iCMS3_Up") & Correlation > 0 & p < 0.05)
    ),
    TRUE ~ (
      (contrast %in% c("iCMS2_Up", "iCMS3_Down") & Correlation > 0 & p < 0.05) |
        (contrast %in% c("iCMS2_Down", "iCMS3_Up") & Correlation < 0 & p < 0.05)
    )
  ))

# 2. Count the correct calls per signature and sort from most to fewest
correct_summary <- p_data_all %>%
  group_by(geneset, contrast) %>%
  summarise(correct_count = sum(is_correct), .groups = "drop")

total_correct <- correct_summary %>%
  group_by(geneset) %>%
  summarise(total = sum(correct_count), .groups = "drop") %>%
  arrange(total) # Sort ascending so the largest totals end up on top

ordered_genesets <- total_correct$geneset
p_data_all$geneset = factor(p_data_all$geneset, levels = ordered_genesets)
correct_summary$geneset = factor(correct_summary$geneset, levels = ordered_genesets)

# 3. Draw the heatmap (p1)
p1 = ggplot(p_data_all, aes(dataset, geneset)) +
  geom_tile(aes(fill = Correlation), color = "black") +
  scale_fill_gradient2(
    high = "#f58f98",
    mid = "white",
    low = "#90d7ec",
    limits = c(-0.5, 0.5),
    oob = scales::squish
  ) +
  geom_text(aes(label = pstar), col = "black", size = 5) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  ) +
  labs(fill = "Correlation") +
  guides(fill = guide_colorbar(
    title.position = "right",
    label.position = "left"
  ))

# 4. Draw the stacked bar chart on the right (p2)
p2 = ggplot(correct_summary, aes(x = correct_count, y = geneset, fill = contrast)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(data = total_correct, aes(x = total, y = geneset, label = total),
            inherit.aes = FALSE, hjust = -0.2, size = 4.5, color = "black") +
  scale_x_continuous(limits = c(0, max(total_correct$total, na.rm = TRUE) * 1.3)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 14, color = "black"),
    text = element_text(size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "Correct Count", fill = "Contrast")

# 5. Assemble and save
p = p1 + p2 + plot_layout(widths = c(4, 1.2))

setwd("C:/Users/user/Desktop/")
ggsave("CRC_benchmark2_cor.pdf",
       plot = p,
       width = 24,
       height = 7,
       device = "pdf")
################################################################################Summary statistics (training datasets)####
rm(list=ls())
gc()
library(GSEABase)
library(writexl) # Used to export the Excel sheet
setwd("C:/Users/user/Desktop/input")

vector=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
         "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
         "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
         "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
         "GSE73906_f","GSE123807","GSE89124")

for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}

for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  WNT_HL=read.delim("WNT_HL.txt",header=TRUE,sep='\t',check.names=F)
  rownames(WNT_HL)=WNT_HL$Accession
  WNT_HL=WNT_HL[,-1]
  WNT_HL_selected=WNT_HL[vector[i],,drop=F]
  group$group=ifelse(group$group==WNT_HL_selected$Group1,WNT_HL_selected$Group1_Status,WNT_HL_selected$Group2_Status)
  group$Tag <- paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}

all_group=all_group
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=1000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}

all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("C:/Users/user/Desktop/input")
  GeneSets=getGmt("wnt_usual_using.gmt")
  geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  usual=as.data.frame(scale(geneset))
  usual$ID=paste0(vector[cycle])
  all_usual=rbind(all_usual,usual)
}

all_data=merge(all_scores[,c(1,2,3),drop=F],all_usual,by="row.names")
all_data=merge(all_data,all_group,by.x="Row.names",by.y="Tag")
fs=unique(all_data$ID)

p_data_all=data.frame()

for(i in 1:length(fs)){
  all_data_fs=all_data[all_data$ID==fs[i],]
  gene_set=colnames(all_data_fs)[2:19]
  for(j in 1:length(gene_set)){
    all_data_fs_fs=all_data_fs[,c(gene_set[j],"group")]
    all_data_fs_fs$group=factor(all_data_fs_fs$group,levels=c("H","L"))
    colnames(all_data_fs_fs)=c("pathway","group")

    # Run the t-test
    tt = t.test(pathway ~ group, all_data_fs_fs, alternative = "two.sided")

    # Extract the detailed metrics
    mean_H = tt$estimate["mean in group H"]
    mean_L = tt$estimate["mean in group L"]
    FC = mean_H - mean_L                  # Mean difference (effect estimate)
    t_stat = tt$statistic                 # t statistic
    df_val = tt$parameter                 # Degrees of freedom
    p_val = tt$p.value                    # Raw p-value
    ci_lower = tt$conf.int[1]             # 95% confidence interval lower bound
    ci_upper = tt$conf.int[2]             # 95% confidence interval upper bound

    p_data=data.frame(
      dataset = fs[i],
      geneset = gene_set[j],
      mean_H = mean_H,
      mean_L = mean_L,
      mean_diff_FC = FC,
      t_statistic = t_stat,
      df = df_val,
      p_value = p_val,
      CI_95_lower = ci_lower,
      CI_95_upper = ci_upper
    )
    p_data_all=rbind(p_data_all, p_data)
  }
}

# Add FDR correction (Benjamini-Hochberg)
p_data_all$p_adj_FDR = p.adjust(p_data_all$p_value, method = "BH")

# Rename labels
p_data_all$geneset=sub("activity_score", "WntQuant", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "WntQuant_activation", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "WntQuant_inhibition", p_data_all$geneset)

# Export as an Excel sheet
setwd("C:/Users/user/Desktop/")
write_xlsx(p_data_all, path = "benchmark_detailed_statistics.xlsx")
cat("Excel workbook exported to the desktop: benchmark_detailed_statistics.xlsx\n")
################################################################################Summary statistics (validation datasets)####
rm(list=ls())
gc()
library(GSEABase)
library(writexl) # Used to export the Excel sheet
setwd("C:/Users/user/Desktop/input")

vector=c("GSE102208","GSE114059_a","GSE114059_b",
         "GSE114059_c","GSE156082","GSE158578_a",
         "GSE158578_b","GSE188465","GSE188466",
         "GSE21815","GSE17537","GSE72970",
         "GSE188467","GSE29435","GSE29621",
         "GSE38832","GSE57760",
         "GSE64337_a","GSE64337_b","GSE94714",
         "GSE69687","GSE237130","UCSC.CRC_APC",
         "UCSC.CRC_TP53","UCSC.LIHC_CTNNB1","UCSC.LIHC_AXIN1",
         "UCSC.CRC_TumorNormal","GSE222023","GSE215947",
         "GSE178260","GSE135679","GSE60564",
         "GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE56896","GSE17385","GSE63072_a",
         "GSE63072_b","GSE61687","UCSC.CRC_CMS2",
         "GSE33113","GSE39582","GSE2109",
         "GSE13067","GSE13294","GSE14333",
         "GSE17536","GSE20916","GSE23878",
         "GSE35896","GSE37892","KFSYSCC",
         "syn26720761","syn2623706")

for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}

for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  WNT_HL=read.delim("WNT_HL_V.txt",header=TRUE,sep='\t',check.names=F)
  rownames(WNT_HL)=WNT_HL$Accession
  WNT_HL=WNT_HL[,-1]
  WNT_HL_selected=WNT_HL[vector[i],,drop=F]
  group$group=ifelse(group$group==WNT_HL_selected$Group1,WNT_HL_selected$Group1_Status,WNT_HL_selected$Group2_Status)
  group$Tag <- paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}

all_group=all_group
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                      expression_profile=exp,
                                      foundation="relative_ssGSEA",
                                      activation_geneset=NA,
                                      inhibition_geneset=NA,
                                      geneSets_gmt=geneSets,
                                      min.sz=1,
                                      max.sz=10000,
                                      geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                          rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}

all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("C:/Users/user/Desktop/input")
  GeneSets=getGmt("wnt_usual_using.gmt")
  geneset=gsva(expr=as.matrix(exp), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
  geneset=as.data.frame(t(geneset))
  usual=as.data.frame(scale(geneset))
  usual$ID=paste0(vector[cycle])
  all_usual=rbind(all_usual,usual)
}

all_data=merge(all_scores[,c(1,2,3),drop=F],all_usual,by="row.names")
all_data=merge(all_data,all_group,by.x="Row.names",by.y="Tag")
fs=unique(all_data$ID)

p_data_all=data.frame()

for(i in 1:length(fs)){
  all_data_fs=all_data[all_data$ID==fs[i],]
  gene_set=colnames(all_data_fs)[2:19]
  for(j in 1:length(gene_set)){
    all_data_fs_fs=all_data_fs[,c(gene_set[j],"group")]
    all_data_fs_fs$group=factor(all_data_fs_fs$group,levels=c("H","L"))
    colnames(all_data_fs_fs)=c("pathway","group")

    # Run the t-test
    tt = t.test(pathway ~ group, all_data_fs_fs, alternative = "two.sided")

    # Extract the detailed metrics
    mean_H = tt$estimate["mean in group H"]
    mean_L = tt$estimate["mean in group L"]
    FC = mean_H - mean_L                  # Mean difference (effect estimate)
    t_stat = tt$statistic                 # t statistic
    df_val = tt$parameter                 # Degrees of freedom
    p_val = tt$p.value                    # Raw p-value
    ci_lower = tt$conf.int[1]             # 95% confidence interval lower bound
    ci_upper = tt$conf.int[2]             # 95% confidence interval upper bound

    p_data=data.frame(
      dataset = fs[i],
      geneset = gene_set[j],
      mean_H = mean_H,
      mean_L = mean_L,
      mean_diff_FC = FC,
      t_statistic = t_stat,
      df = df_val,
      p_value = p_val,
      CI_95_lower = ci_lower,
      CI_95_upper = ci_upper
    )
    p_data_all=rbind(p_data_all, p_data)
  }
}

# Add FDR correction (Benjamini-Hochberg)
p_data_all$p_adj_FDR = p.adjust(p_data_all$p_value, method = "BH")

# Rename labels
p_data_all$geneset=sub("activity_score", "WntQuant", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "WntQuant_activation", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "WntQuant_inhibition", p_data_all$geneset)

# Export as an Excel sheet
setwd("C:/Users/user/Desktop/")
write_xlsx(p_data_all, path = "testing_benchmark_detailed_statistics.xlsx")
cat("Excel workbook exported to the desktop: testing_benchmark_detailed_statistics.xlsx\n")
################################################################################Multivariable Cox regression for OS####
rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(survival)
library(broom)
library(tidyr)
library(forestploter)

setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/cluster.Rdata")

clinical_molecular_public_all <- read.table(
  file = "clinical_molecular_public_all.txt",
  header = TRUE, sep = "\t", check.names = FALSE
)

data <- merge(cluster, clinical_molecular_public_all,
              by.x = "Tag", by.y = "sample")

data_analysis <- data[, c("Tag","cluster","dataset","age","gender","stage",
                          "pt","pn","pm","tnm","msi","cimp","kras_mut",
                          "braf_mut","cms_label","osMo","osStat")]

data_analysis <- data_analysis %>%
  mutate(
    # Reclassify tnm
    tnm = case_when(
      tnm == "I" ~ "I",
      tnm %in% c("IIA", "IIB/IIC") ~ "II",
      tnm %in% c("IIIA", "IIIB", "IIIC") ~ "III",
      tnm == "IV" ~ "IV",
      TRUE ~ NA_character_
    ),
    tnm = factor(tnm, levels = c("I", "II", "III", "IV")),

    # Merge clusters: WSS1, WSS2 and WSS3 become Other; WSS4 is kept
    cluster = case_when(
      cluster %in% c("WSS1", "WSS2", "WSS3") ~ "Other",
      cluster == "WSS4" ~ "WSS4",
      TRUE ~ NA_character_
    ),
    cluster = factor(cluster, levels = c("Other", "WSS4")),

    # Merge cms_label: CMS1, CMS2, CMS3 and NOLBL become Other; CMS4 is kept
    cms_label = case_when(
      cms_label %in% c("CMS1", "CMS2", "CMS3", "NOLBL") ~ "Other",
      cms_label == "CMS4" ~ "CMS4",
      TRUE ~ NA_character_
    ),
    cms_label = factor(cms_label, levels = c("Other", "CMS4")),

    # Normalise the subgroup names (gender, msi)
    gender = case_when(
      tolower(gender) == "female" ~ "Female",
      tolower(gender) == "male" ~ "Male",
      TRUE ~ gender
    ),
    msi = case_when(
      tolower(msi) == "msi" ~ "MSI",
      tolower(msi) == "mss" ~ "MSS",
      TRUE ~ msi
    )
  )

# Lock the survival data (no extra n_missing row filtering)
dat_os_filtered <- data_analysis %>%
  filter(!is.na(osMo), !is.na(osStat), osMo > 0)

cat("\nSamples analysed:", nrow(dat_os_filtered), "\n")

dat_os_filtered <- dat_os_filtered %>%
  mutate(
    gender    = factor(gender),
    msi       = factor(msi),
    cimp      = factor(cimp),
    kras_mut  = factor(kras_mut),
    braf_mut  = factor(braf_mut),
    osMo      = as.numeric(osMo),
    osStat    = as.numeric(as.character(osStat))
  ) %>%
  mutate(
    cluster   = relevel(cluster, ref = "WSS4"),
    cms_label = relevel(cms_label, ref = "CMS4"),
    cimp      = relevel(cimp, ref = "CIMP.Neg")
  )

formula_full <- Surv(osMo, osStat) ~ cluster + age + gender +
  tnm + msi + cimp + kras_mut + braf_mut + cms_label

fit_full <- coxph(formula_full, data = dat_os_filtered)
cat("\n===== Main model summary =====\n")
print(summary(fit_full))

res_tidy <- tidy(fit_full, exponentiate = TRUE, conf.int = TRUE)

factors_in_model <- c("cluster", "gender", "tnm", "msi", "cimp", "kras_mut", "braf_mut", "cms_label")
factors_in_model <- intersect(factors_in_model, names(dat_os_filtered))

# Display-name dictionary for the SCI subtypes
var_labels <- c(
  cluster   = "WSS labels",
  gender    = "Gender",
  tnm       = "TNM stages",
  msi       = "Microsatellite",
  cimp      = "Methylation",
  kras_mut  = "KRAS mutation",
  braf_mut  = "BRAF mutation",
  cms_label = "CMS labels"
)

plot_rows <- list()

for (v in factors_in_model) {
  lvls <- levels(dat_os_filtered[[v]])
  display_var_name <- ifelse(v %in% names(var_labels), var_labels[v], v)

  # Main variable heading row
  plot_rows[[paste0(v, "_header")]] <- data.frame(
    col1 = display_var_name,
    col2 = "",
    col3 = "",
    col4 = "",
    col5 = paste(rep(" ", 25), collapse = " "),
    col6 = "",
    HR   = NA,
    lower = NA,
    upper = NA,
    stringsAsFactors = FALSE
  )

  # Rows for each subtype / level
  for (i in seq_along(lvls)) {
    l <- lvls[i]
    n_val <- sum(dat_os_filtered[[v]] == l, na.rm = TRUE)

    if (i == 1) {
      # Reference group
      plot_rows[[paste0(v, "_", l)]] <- data.frame(
        col1 = "",
        col2 = paste0("   ", l),
        col3 = as.character(n_val),
        col4 = "Reference",
        col5 = paste(rep(" ", 25), collapse = " "),
        col6 = "",
        HR   = 1.0,
        lower = 1.0,
        upper = 1.0,
        stringsAsFactors = FALSE
      )
    } else {
      # Non-reference group
      term_name <- paste0(v, l)
      sub_res <- res_tidy %>% filter(term == term_name)

      hr_str <- if(nrow(sub_res) > 0) sprintf("%.2f (%.2f-%.2f)", sub_res$estimate, sub_res$conf.low, sub_res$conf.high) else ""
      p_str  <- if(nrow(sub_res) > 0) {
        if(sub_res$p.value < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value)
      } else ""

      plot_rows[[paste0(v, "_", l)]] <- data.frame(
        col1 = "",
        col2 = paste0("   ", l),
        col3 = as.character(n_val),
        col4 = hr_str,
        col5 = paste(rep(" ", 25), collapse = " "),
        col6 = p_str,
        HR   = if(nrow(sub_res) > 0) sub_res$estimate else NA,
        lower = if(nrow(sub_res) > 0) sub_res$conf.low else NA,
        upper = if(nrow(sub_res) > 0) sub_res$conf.high else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Handle the continuous variable (age)
continuous_vars <- intersect(c("age"), names(dat_os_filtered))
for (cv in continuous_vars) {
  sub_res <- res_tidy %>% filter(term == cv)
  n_val <- sum(!is.na(dat_os_filtered[[cv]]))
  display_cv_name <- ifelse(cv == "age", "Age", cv)

  hr_str <- if(nrow(sub_res) > 0) {
    sprintf("%.2f (%.2f-%.2f)", sub_res$estimate, sub_res$conf.low, sub_res$conf.high)
  } else {
    ""
  }

  p_str <- if(nrow(sub_res) > 0) {
    if(sub_res$p.value < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value)
  } else {
    ""
  }

  plot_rows[[cv]] <- data.frame(
    col1 = display_cv_name,
    col2 = "",
    col3 = as.character(n_val),
    col4 = hr_str,
    col5 = paste(rep(" ", 25), collapse = " "),
    col6 = p_str,
    HR   = if(nrow(sub_res) > 0) sub_res$estimate else NA,
    lower = if(nrow(sub_res) > 0) sub_res$conf.low else NA,
    upper = if(nrow(sub_res) > 0) sub_res$conf.high else NA,
    stringsAsFactors = FALSE
  )
}

# Merge all rows
df_forest_input <- bind_rows(plot_rows)

# Set the header
colnames(df_forest_input)[1:6] <- c(
  "Variable",
  "Subgroup",
  "N",
  "Hazard ratio (95% CI)",
  "",
  "p-value"
)

pdf("C:/Users/user/Desktop/Multivariable_Cox_Forestplot.pdf", width = 11, height = 7.5)

p_final <- forest(
  df_forest_input[, 1:6],
  est = df_forest_input$HR,
  lower = df_forest_input$lower,
  upper = df_forest_input$upper,
  size = 0.4,
  ci_column = 5,
  ref_line = 1,
  xlim = c(0, 15),
  ticks_at = c(0, 1, 5, 10, 15),
  xlab = "Hazard ratio (95% CI)"
)

print(p_final)
dev.off()

cat("\nMissing-value row filtering disabled; updated forest plot saved to the desktop: C:/Users/user/Desktop/Multivariable_Cox_Forestplot.pdf\n")

################################################################################Multivariable Cox regression for RFS####
rm(list = ls())
gc()

library(dplyr)
library(ggplot2)
library(survival)
library(broom)
library(tidyr)
library(forestploter)

setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/cluster.Rdata")

clinical_molecular_public_all <- read.table(
  file = "clinical_molecular_public_all.txt",
  header = TRUE, sep = "\t", check.names = FALSE
)

data <- merge(cluster, clinical_molecular_public_all,
              by.x = "Tag", by.y = "sample")

# Set the survival endpoint to the RFS fields (rfsMo, rfsStat)
data_analysis <- data[, c("Tag","cluster","dataset","age","gender","stage",
                          "pt","pn","pm","tnm","msi","cimp","kras_mut",
                          "braf_mut","cms_label","rfsMo","rfsStat")]

data_analysis <- data_analysis %>%
  mutate(
    # Reclassify tnm
    tnm = case_when(
      tnm == "I" ~ "I",
      tnm %in% c("IIA", "IIB/IIC") ~ "II",
      tnm %in% c("IIIA", "IIIB", "IIIC") ~ "III",
      tnm == "IV" ~ "IV",
      TRUE ~ NA_character_
    ),
    tnm = factor(tnm, levels = c("I", "II", "III", "IV")),

    # Merge clusters: WSS1, WSS2 and WSS3 become Other; WSS4 is kept
    cluster = case_when(
      cluster %in% c("WSS1", "WSS2", "WSS3") ~ "Other",
      cluster == "WSS4" ~ "WSS4",
      TRUE ~ NA_character_
    ),
    cluster = factor(cluster, levels = c("Other", "WSS4")),

    # Merge cms_label: CMS1, CMS2, CMS3 and NOLBL become Other; CMS4 is kept
    cms_label = case_when(
      cms_label %in% c("CMS1", "CMS2", "CMS3", "NOLBL") ~ "Other",
      cms_label == "CMS4" ~ "CMS4",
      TRUE ~ NA_character_
    ),
    cms_label = factor(cms_label, levels = c("Other", "CMS4")),

    # Normalise the subgroup names (gender, msi)
    gender = case_when(
      tolower(gender) == "female" ~ "Female",
      tolower(gender) == "male" ~ "Male",
      TRUE ~ gender
    ),
    msi = case_when(
      tolower(msi) == "msi" ~ "MSI",
      tolower(msi) == "mss" ~ "MSS",
      TRUE ~ msi
    )
  )

# Lock the RFS survival data (no extra n_missing row filtering)
dat_rfs_filtered <- data_analysis %>%
  filter(!is.na(rfsMo), !is.na(rfsStat), rfsMo > 0)

cat("\nRFS samples analysed:", nrow(dat_rfs_filtered), "\n")

dat_rfs_filtered <- dat_rfs_filtered %>%
  mutate(
    gender    = factor(gender),
    msi       = factor(msi),
    cimp      = factor(cimp),
    kras_mut  = factor(kras_mut),
    braf_mut  = factor(braf_mut),
    rfsMo     = as.numeric(rfsMo),
    rfsStat   = as.numeric(as.character(rfsStat))
  ) %>%
  mutate(
    cluster   = relevel(cluster, ref = "WSS4"),     # cluster uses WSS4 as reference
    cms_label = relevel(cms_label, ref = "CMS4"),   # cms_label uses CMS4 as reference
    cimp      = relevel(cimp, ref = "CIMP.Neg")     # cimp uses CIMP.Neg as reference
  )

formula_full <- Surv(rfsMo, rfsStat) ~ cluster + age + gender +
  tnm + msi + cimp + kras_mut + braf_mut + cms_label

fit_full <- coxph(formula_full, data = dat_rfs_filtered)
cat("\n===== Main model summary (RFS) =====\n")
print(summary(fit_full))

res_tidy <- tidy(fit_full, exponentiate = TRUE, conf.int = TRUE)

factors_in_model <- c("cluster", "gender", "tnm", "msi", "cimp", "kras_mut", "braf_mut", "cms_label")
factors_in_model <- intersect(factors_in_model, names(dat_rfs_filtered))

# Display-name dictionary for the SCI subtypes
var_labels <- c(
  cluster   = "WSS labels",
  gender    = "Gender",
  tnm       = "TNM stages",
  msi       = "Microsatellite",
  cimp      = "Methylation",
  kras_mut  = "KRAS mutation",
  braf_mut  = "BRAF mutation",
  cms_label = "CMS labels"
)

plot_rows <- list()

for (v in factors_in_model) {
  lvls <- levels(dat_rfs_filtered[[v]])
  display_var_name <- ifelse(v %in% names(var_labels), var_labels[v], v)

  # Main variable heading row
  plot_rows[[paste0(v, "_header")]] <- data.frame(
    col1 = display_var_name,
    col2 = "",
    col3 = "",
    col4 = "",
    col5 = paste(rep(" ", 25), collapse = " "),
    col6 = "",
    HR   = NA,
    lower = NA,
    upper = NA,
    stringsAsFactors = FALSE
  )

  # Rows for each subtype / level
  for (i in seq_along(lvls)) {
    l <- lvls[i]
    n_val <- sum(dat_rfs_filtered[[v]] == l, na.rm = TRUE)

    if (i == 1) {
      # Reference group
      plot_rows[[paste0(v, "_", l)]] <- data.frame(
        col1 = "",
        col2 = paste0("   ", l),
        col3 = as.character(n_val),
        col4 = "Reference",
        col5 = paste(rep(" ", 25), collapse = " "),
        col6 = "",
        HR   = 1.0,
        lower = 1.0,
        upper = 1.0,
        stringsAsFactors = FALSE
      )
    } else {
      # Non-reference group
      term_name <- paste0(v, l)
      sub_res <- res_tidy %>% filter(term == term_name)

      hr_str <- if(nrow(sub_res) > 0) sprintf("%.2f (%.2f-%.2f)", sub_res$estimate, sub_res$conf.low, sub_res$conf.high) else ""
      p_str  <- if(nrow(sub_res) > 0) {
        if(sub_res$p.value < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value)
      } else ""

      plot_rows[[paste0(v, "_", l)]] <- data.frame(
        col1 = "",
        col2 = paste0("   ", l),
        col3 = as.character(n_val),
        col4 = hr_str,
        col5 = paste(rep(" ", 25), collapse = " "),
        col6 = p_str,
        HR   = if(nrow(sub_res) > 0) sub_res$estimate else NA,
        lower = if(nrow(sub_res) > 0) sub_res$conf.low else NA,
        upper = if(nrow(sub_res) > 0) sub_res$conf.high else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Handle the continuous variable (age)
continuous_vars <- intersect(c("age"), names(dat_rfs_filtered))
for (cv in continuous_vars) {
  sub_res <- res_tidy %>% filter(term == cv)
  n_val <- sum(!is.na(dat_rfs_filtered[[cv]]))
  display_cv_name <- ifelse(cv == "age", "Age", cv)

  hr_str <- if(nrow(sub_res) > 0) {
    sprintf("%.2f (%.2f-%.2f)", sub_res$estimate, sub_res$conf.low, sub_res$conf.high)
  } else {
    ""
  }

  p_str <- if(nrow(sub_res) > 0) {
    if(sub_res$p.value < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value)
  } else {
    ""
  }

  plot_rows[[cv]] <- data.frame(
    col1 = display_cv_name,
    col2 = "",
    col3 = as.character(n_val),
    col4 = hr_str,
    col5 = paste(rep(" ", 25), collapse = " "),
    col6 = p_str,
    HR   = if(nrow(sub_res) > 0) sub_res$estimate else NA,
    lower = if(nrow(sub_res) > 0) sub_res$conf.low else NA,
    upper = if(nrow(sub_res) > 0) sub_res$conf.high else NA,
    stringsAsFactors = FALSE
  )
}

# Merge all rows
df_forest_input <- bind_rows(plot_rows)

# Set the header
colnames(df_forest_input)[1:6] <- c(
  "Variable",
  "Subgroup",
  "N",
  "Hazard ratio (95% CI)",
  "",
  "p-value"
)

pdf("C:/Users/user/Desktop/Multivariable_Cox_Forestplot_RFS.pdf", width = 11, height = 7.5)

p_final <- forest(
  df_forest_input[, 1:6],
  est = df_forest_input$HR,
  lower = df_forest_input$lower,
  upper = df_forest_input$upper,
  size = 0.4,
  ci_column = 5,
  ref_line = 1,
  xlim = c(0, 15),
  ticks_at = c(0, 1, 5, 10, 15),
  xlab = "Hazard ratio (95% CI)"
)

print(p_final)
dev.off()

cat("\nRFS forest plot (WSS4/Other, x-axis 0-15) saved to the desktop: C:/Users/user/Desktop/Multivariable_Cox_Forestplot_RFS.pdf\n")
################################################################################Organising the TCGA files####
rm(list = ls())
gc()

library(dplyr)
library(reshape2)

setwd("C:/Users/user/Desktop/input")

load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")

# Core check: no NA or empty cluster labels in pancancer_cluster
pancancer_cluster = pancancer_cluster[!is.na(pancancer_cluster$cluster) & pancancer_cluster$cluster != "", ]

survival_data = read.delim("Survival_SupplementalTable_S1_20171025_xena_sp",
                            check.names = FALSE)
mut_data = read.delim("mc3.v0.2.8.PUBLIC.nonsilentGene.xena",
                       check.names = FALSE)

# 1. Age handling
age = survival_data[, c(1, 4)]
colnames(age) = c("Tag", "age")
age = na.omit(age)
age$age = ifelse(age$age <= 40, "age<=40",
                  ifelse(age$age <= 60, "40<age<=60", "age>60"))
age$age = factor(age$age, levels = c("age<=40", "40<age<=60", "age>60"))

# 2. Sex handling
gender = survival_data[, c(1, 5)]
colnames(gender) = c("Tag", "gender")
gender$gender = ifelse(gender$gender == "MALE", "Male", "Female")

# 3. Race handling
race = survival_data[, c(1, 6)]
colnames(race) = c("Tag", "race")
race$race[race$race %in% c("", "[Not Evaluated]", "[Unknown]")] = "Other/Unknown"
race$race = ifelse(race$race == "WHITE", "White",
                    ifelse(race$race == "BLACK OR AFRICAN AMERICAN", "Black",
                           ifelse(race$race == "ASIAN", "Asian", "Other/Unknown")))

# 4. Stage handling: prefer the pathologic stage (column 7), fill gaps from the clinical stage (column 8), map to Stage 0-4
stage = survival_data[, c(1, 7, 8)]
colnames(stage) = c("Tag", "stage_path", "stage_clin")

# Strip whitespace
stage$stage_path = trimws(as.character(stage$stage_path))
stage$stage_clin = trimws(as.character(stage$stage_clin))

# Values that cannot be mapped become NA
na_set = c("", "[Discrepancy]", "[Unknown]", "Stage X", "X", "Unknown", "I/II NOS")
stage$stage_path[stage$stage_path %in% na_set] = NA
stage$stage_clin[stage$stage_clin %in% na_set] = NA

# Normalise any spelling to Stage 0/1/2/3/4, or NA
normalize_stage = function(x) {
  ifelse(
    x %in% c("0", "Stage 0", "IS", "Stage IS"), "Stage 0",
    ifelse(
      x %in% c("I", "Stage I", "Stage IA", "Stage IA1", "Stage IA2",
               "Stage IB", "Stage IB1", "Stage IB2", "Stage IC"), "Stage 1",
      ifelse(
        x %in% c("II", "IIa", "IIb", "Stage II", "Stage IIA", "Stage IIA1",
                 "Stage IIA2", "Stage IIB", "Stage IIC"), "Stage 2",
        ifelse(
          x %in% c("III", "Stage III", "Stage IIIA", "Stage IIIB",
                   "Stage IIIC", "Stage IIIC1", "Stage IIIC2"), "Stage 3",
          ifelse(
            x %in% c("IV", "Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"), "Stage 4",
            NA
          )
        )
      )
    )
  )
}

# Pathologic stage wins; clinical stage fills the gaps
stage$stage_path_n = normalize_stage(stage$stage_path)
stage$stage_clin_n = normalize_stage(stage$stage_clin)
stage$stage = ifelse(!is.na(stage$stage_path_n), stage$stage_path_n, stage$stage_clin_n)

# Keep Tag and stage only, then convert to a factor
stage = stage[, c("Tag", "stage")]
stage$stage = factor(stage$stage, levels = c("Stage 0", "Stage 1", "Stage 2", "Stage 3", "Stage 4"))

# Check the result
table(stage$stage, useNA = "ifany")

# 5. Tumour status and response (grade removed)
tumor_status = survival_data[, c(1, 15)]
colnames(tumor_status) = c("Tag", "tumor_status")
tumor_status$tumor_status[tumor_status$tumor_status %in% c("", "[Discrepancy]")] = NA
tumor_status$tumor_status = ifelse(tumor_status$tumor_status == "TUMOR FREE",
                                    "Tumor free", "With tumor")

outcome = survival_data[, c(1, 23)]
colnames(outcome) = c("Tag", "outcome")
outcome$outcome[outcome$outcome %in% c("", "[Discrepancy]", "[Unknown]", "[Not Evaluated]")] = NA
outcome$outcome = ifelse(
  outcome$outcome %in% c("Complete Remission/Response",
                         "No Measureable Tumor or Tumor Markers",
                         "Normalization of Tumor Markers, but Residual Tumor Mass",
                         "Partial Remission/Response"),
  "CR/PR", "SD/PD")

# 6. Mutation data (fixes the earlier matrix-structure crash)
gene_list = c("APC","AXIN1","AXIN2","RNF43",
               "DKK1","DKK3","TCF7L2","GSK3B","CTNNB1")

mut_sub = mut_data[mut_data[[1]] %in% gene_list, ]
rownames(mut_sub) = mut_sub[[1]]
mut_sub = mut_sub[, -1, drop = FALSE]

mut_mat = as.matrix(t(mut_sub))
mut_mat[mut_mat == 1] = "MT"
mut_mat[mut_mat != "MT"] = "WT"

mut_t = as.data.frame(mut_mat, stringsAsFactors = FALSE)
mut_t$Tag = rownames(mut_t)

# 7. Survival data cleaning
os_tbl = survival_data[, c(1, 26, 27)]
colnames(os_tbl) = c("Tag", "OS", "OS.time")
os_tbl = os_tbl[!is.na(os_tbl$OS) & !is.na(os_tbl$OS.time) & os_tbl$OS.time > 0, ]

dss_tbl = survival_data[, c(1, 28, 29)]
colnames(dss_tbl) = c("Tag", "DSS", "DSS.time")
dss_tbl = dss_tbl[!is.na(dss_tbl$DSS) & !is.na(dss_tbl$DSS.time) & dss_tbl$DSS.time > 0, ]

dfi_tbl = survival_data[, c(1, 30, 31)]
colnames(dfi_tbl) = c("Tag", "DFI", "DFI.time")
dfi_tbl = dfi_tbl[!is.na(dfi_tbl$DFI) & !is.na(dfi_tbl$DFI.time) & dfi_tbl$DFI.time > 0, ]

pfi_tbl = survival_data[, c(1, 32, 33)]
colnames(pfi_tbl) = c("Tag", "PFI", "PFI.time")
pfi_tbl = pfi_tbl[!is.na(pfi_tbl$PFI) & !is.na(pfi_tbl$PFI.time) & pfi_tbl$PFI.time > 0, ]

# 8. Build the final data frames with merge()
build_data = function(surv_tbl) {
  d = pancancer_cluster
  d = merge(d, age,          by = "Tag", all = TRUE)
  d = merge(d, gender,       by = "Tag", all = TRUE)
  d = merge(d, race,         by = "Tag", all = TRUE)
  d = merge(d, stage,        by = "Tag", all = TRUE)
  d = merge(d, tumor_status, by = "Tag", all = TRUE)
  d = merge(d, outcome,      by = "Tag", all = TRUE)
  d = merge(d, mut_t,        by = "Tag", all = TRUE)
  d = merge(d, surv_tbl,     by = "Tag")
  return(d)
}

data_OS = build_data(os_tbl)
data_DSS = build_data(dss_tbl)
data_DFI = build_data(dfi_tbl)
data_PFI = build_data(pfi_tbl)

# 9. Inspect the data, including the NA check on the cluster variable
check_data = function(d, name, evt, time) {
  cat("\n=====", name, "=====\n")
  cat("Dimensions:", dim(d), "\n")
  cat("Missing cluster:", sum(is.na(d$cluster)), "\n")
  cat("Missing event:", sum(is.na(d[[evt]])), "\n")
  cat("Missing time:", sum(is.na(d[[time]])), "\n")
  cat("Event distribution:\n"); print(table(d[[evt]], useNA = "ifany"))
}

check_data(data_OS,  "data_OS",  "OS",  "OS.time")
check_data(data_DSS, "data_DSS", "DSS", "DSS.time")
check_data(data_DFI, "data_DFI", "DFI", "DFI.time")
check_data(data_PFI, "data_PFI", "PFI", "PFI.time")

# 10. Save
save(data_OS,  file = "data_OS.Rdata")
save(data_DSS, file = "data_DSS.Rdata")
save(data_DFI, file = "data_DFI.Rdata")
save(data_PFI, file = "data_PFI.Rdata")

################################################################################Multivariable visualisation of OS####
rm(list = ls())
gc()
library(survival)
library(dplyr)
library(broom)
library(forestploter)
setwd("C:/Users/user/Desktop/input")
load("data_OS.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
data_OS = merge(data_OS, pancancer_group, by.x = "Tag", by.y = "row.names")

# Keep tumour samples only: the numeric suffix of Tag must be <= 10
data_OS = data_OS %>%
  mutate(sample_code = as.numeric(sub(".*-", "", Tag))) %>%
  filter(sample_code <= 10) %>%
  select(-sample_code)

# ---------- Cancer type: 12 organ-system groups
data_OS = data_OS %>%
  mutate(
    cancer_types = case_when(
      is.na(Project) ~ "Unknown",
      Project %in% c("OV", "UCEC", "CESC", "BRCA") ~ "Gynecologic",
      Project %in% c("BLCA", "PRAD", "TGCT", "KIRC", "KICH", "KIRP") ~ "Urologic",
      Project %in% c("THCA", "ACC") ~ "Endocrine",
      Project %in% c("ESCA", "STAD", "COAD", "READ") ~ "Core_gastrointestinal",
      Project %in% c("LIHC", "PAAD", "CHOL") ~ "Developmental_gastrointestinal",
      Project %in% c("HNSC") ~ "Head_and_neck",
      Project %in% c("LUAD", "LUSC", "MESO") ~ "Thoracic",
      Project %in% c("GBM", "LGG") ~ "Central_nervous_system",
      Project %in% c("SARC", "UCS") ~ "Soft_tissue",
      Project %in% c("PCPG") ~ "Neural_crest_derived",
      Project %in% c("SKCM", "UVM") ~ "Melanocytic",
      Project %in% c("LAML", "DLBC", "THYM") ~ "Hematologic_and_lymphatic",
      TRUE ~ "Unknown"
    )
  )
data_OS$cancer_types = as.character(data_OS$cancer_types)
data_OS$cancer_types[is.na(data_OS$cancer_types) | data_OS$cancer_types == ""] = "Unknown"
data_OS$cancer_types = factor(data_OS$cancer_types,
                               levels = c("Gynecologic", "Urologic", "Endocrine",
                                          "Core_gastrointestinal", "Developmental_gastrointestinal",
                                          "Head_and_neck", "Thoracic", "Central_nervous_system",
                                          "Soft_tissue", "Neural_crest_derived", "Melanocytic",
                                          "Hematologic_and_lymphatic", "Unknown"))

# ---------- Stage: unmappable values and Stage 0 become Unknown
data_OS$stage = trimws(as.character(data_OS$stage))
data_OS$stage[data_OS$stage %in% c("", "NA", "[Discrepancy]", "[Unknown]",
                                   "Stage X", "X", "Unknown", "I/II NOS",
                                   "Stage 0")] = "Unknown"
data_OS$stage[is.na(data_OS$stage)] = "Unknown"
data_OS$stage = factor(data_OS$stage,
                        levels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Unknown"))

# ---------- WSS subtype
data_OS$level = as.character(data_OS$cluster)
data_OS$level[is.na(data_OS$level) | data_OS$level == ""] = "Unknown"
data_OS$level = factor(data_OS$level,
                        levels = c("panWSS4", "panWSS1", "panWSS2", "panWSS3", "Unknown"))

# ---------- Demographics and clinical variables: missing becomes Unknown
data_OS$age = as.character(data_OS$age)
data_OS$age[is.na(data_OS$age) | data_OS$age == ""] = "Unknown"
data_OS$age = factor(data_OS$age,
                      levels = c("age<=40", "40<age<=60", "age>60", "Unknown"))

data_OS$gender = as.character(data_OS$gender)
data_OS$gender[is.na(data_OS$gender) | data_OS$gender == ""] = "Unknown"
data_OS$gender = factor(data_OS$gender, levels = c("Female", "Male", "Unknown"))

if ("race" %in% colnames(data_OS)) {
  data_OS$race = as.character(data_OS$race)
  data_OS$race[data_OS$race %in% c("", "Other/Unknown")] = "Unknown"
  data_OS$race[is.na(data_OS$race)] = "Unknown"
  data_OS$race = factor(data_OS$race, levels = c("White", "Black", "Asian", "Unknown"))
}

if ("outcome" %in% colnames(data_OS)) {
  data_OS$outcome = as.character(data_OS$outcome)
  data_OS$outcome[data_OS$outcome %in% c("", "NA")] = "Unknown"
  data_OS$outcome[is.na(data_OS$outcome)] = "Unknown"
  data_OS$outcome = factor(data_OS$outcome,
                            levels = c("CR/PR", "SD/PD", "Unknown"))
}

# ---------- Mutation variables: APC only
gene_list = c("APC")
for (g in gene_list) {
  if (g %in% colnames(data_OS)) {
    data_OS[[g]] = as.character(data_OS[[g]])
    data_OS[[g]][data_OS[[g]] %in% c("", "NA")] = "Unknown"
    data_OS[[g]][is.na(data_OS[[g]])] = "Unknown"
    data_OS[[g]] = factor(data_OS[[g]], levels = c("WT", "MT", "Unknown"))
  }
}

# ---------- Covariates and the Cox model (tumor_status and the other genes removed)
all_covariates = c("level", "cancer_types", "stage",
                    "age", "gender", "race", "outcome", gene_list)
all_covariates = all_covariates[all_covariates %in% colnames(data_OS)]

# Drop factor levels that have no samples
for (v in all_covariates) {
  if (is.factor(data_OS[[v]])) data_OS[[v]] = droplevels(data_OS[[v]])
}

cox_formula = as.formula(paste("Surv(OS.time, OS) ~", paste(all_covariates, collapse = " + ")))
fit_cox = coxph(cox_formula, data = data_OS)

res_tidy = tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)

# ---------- Display-name mapping
var_labels = c(
  "level"        = "panWSS labels",
  "cancer_types" = "Cancer types",
  "stage"        = "Stage",
  "gender"       = "Gender",
  "age"          = "Age",
  "race"         = "Race",
  "outcome"      = "Outcome"
)
for (g in gene_list) {
  var_labels[g] = paste0(g, " mutation")
}

# Readable labels for the 12 cancer-type levels
cancer_level_labels = c(
  "Gynecologic"                      = "Gynecologic",
  "Urologic"                         = "Urologic",
  "Endocrine"                        = "Endocrine",
  "Core_gastrointestinal"            = "Core gastrointestinal",
  "Developmental_gastrointestinal"   = "Developmental gastrointestinal",
  "Head_and_neck"                    = "Head and neck",
  "Thoracic"                         = "Thoracic",
  "Central_nervous_system"           = "Central nervous system",
  "Soft_tissue"                      = "Soft tissue",
  "Neural_crest_derived"             = "Neural-crest-derived",
  "Melanocytic"                      = "Melanocytic",
  "Hematologic_and_lymphatic"        = "Hematologic and lymphatic",
  "Unknown"                          = "Unknown"
)

factors_in_model = all_covariates
plot_rows = list()

# ---------- Build the forest-plot table rows
for (v in factors_in_model) {
  if (!v %in% names(data_OS)) next
  lvls = levels(data_OS[[v]])
  if (is.null(lvls)) next

  display_var_name = ifelse(v %in% names(var_labels), var_labels[v], v)

  plot_rows[[paste0(v, "_header")]] = data.frame(
    col1 = display_var_name, col2 = "", col3 = "", col4 = "",
    col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
    HR = NA, lower = NA, upper = NA, stringsAsFactors = FALSE
  )

  for (i in seq_along(lvls)) {
    l = lvls[i]
    n_val = sum(data_OS[[v]] == l, na.rm = TRUE)

    if (n_val == 0) next

    l_display = if (v == "cancer_types" && l %in% names(cancer_level_labels)) {
      cancer_level_labels[[l]]
    } else {
      l
    }

    if (i == 1) {
      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = "Reference", col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
        HR = 1.0, lower = 1.0, upper = 1.0, stringsAsFactors = FALSE
      )
    } else {
      term_name = paste0(v, l)
      sub_res = res_tidy %>% filter(term == term_name)

      hr_str = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1])) {
        sprintf("%.2f (%.2f–%.2f)", sub_res$estimate[1], sub_res$conf.low[1], sub_res$conf.high[1])
      } else ""

      p_str = if (nrow(sub_res) > 0 && !is.na(sub_res$p.value[1])) {
        if (sub_res$p.value[1] < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value[1])
      } else ""

      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = hr_str, col5 = paste(rep(" ", 25), collapse = " "), col6 = p_str,
        HR = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1]))  sub_res$estimate[1]  else NA,
        lower = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.low[1]))  sub_res$conf.low[1]  else NA,
        upper = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.high[1])) sub_res$conf.high[1] else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_forest_input = bind_rows(plot_rows)
colnames(df_forest_input)[1:6] = c("Variable", "Subgroup", "N",
                                    "Hazard ratio (95% CI)", "", "p-value")

# ---------- Render and export (HR axis max 7)
pdf("C:/Users/user/Desktop/Multivariate_Cox_OS_APC_Forestplot.pdf", width = 20, height = 24)

p_final = forest(
  df_forest_input[, 1:6],
  est = df_forest_input$HR,
  lower = df_forest_input$lower,
  upper = df_forest_input$upper,
  size = 0.5,
  ci_column = 5,
  ref_line = 1,
  xlim = c(0, 7),
  ticks_at = c(0, 1, 2, 3, 4, 5, 6, 7),
  xlab = "Hazard ratio (95% CI)"
)

print(p_final)
dev.off()

print(p_final)
cat("\nMultivariable OS Cox forest plot (APC only, HR axis max 7) saved to the desktop.\n")

################################################################################Multivariable visualisation of DSS####
rm(list = ls())
gc()

library(survival)
library(dplyr)
library(broom)
library(forestploter)

setwd("C:/Users/user/Desktop/input")
load("data_DSS.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")

data_DSS = merge(data_DSS, pancancer_group, by.x = "Tag", by.y = "row.names")

# Keep tumour samples only: drop samples whose numeric suffix is > 10
data_DSS = data_DSS %>%
  mutate(sample_code = as.numeric(sub(".*-", "", Tag))) %>%
  filter(sample_code <= 10) %>%
  select(-sample_code)

# ---------- Cancer type: 12 organ-system groups
data_DSS = data_DSS %>%
  mutate(
    cancer_types = case_when(
      is.na(Project) ~ "Unknown",
      Project %in% c("OV", "UCEC", "CESC", "BRCA") ~ "Gynecologic",
      Project %in% c("BLCA", "PRAD", "TGCT", "KIRC", "KICH", "KIRP") ~ "Urologic",
      Project %in% c("THCA", "ACC") ~ "Endocrine",
      Project %in% c("ESCA", "STAD", "COAD", "READ") ~ "Core_gastrointestinal",
      Project %in% c("LIHC", "PAAD", "CHOL") ~ "Developmental_gastrointestinal",
      Project %in% c("HNSC") ~ "Head_and_neck",
      Project %in% c("LUAD", "LUSC", "MESO") ~ "Thoracic",
      Project %in% c("GBM", "LGG") ~ "Central_nervous_system",
      Project %in% c("SARC", "UCS") ~ "Soft_tissue",
      Project %in% c("PCPG") ~ "Neural_crest_derived",
      Project %in% c("SKCM", "UVM") ~ "Melanocytic",
      Project %in% c("LAML", "DLBC", "THYM") ~ "Hematologic_and_lymphatic",
      TRUE ~ "Unknown"
    )
  )
data_DSS$cancer_types = as.character(data_DSS$cancer_types)
data_DSS$cancer_types[is.na(data_DSS$cancer_types) | data_DSS$cancer_types == ""] = "Unknown"
data_DSS$cancer_types = factor(data_DSS$cancer_types,
                                levels = c("Gynecologic", "Urologic", "Endocrine",
                                           "Core_gastrointestinal", "Developmental_gastrointestinal",
                                           "Head_and_neck", "Thoracic", "Central_nervous_system",
                                           "Soft_tissue", "Neural_crest_derived", "Melanocytic",
                                           "Hematologic_and_lymphatic", "Unknown"))
table(data_DSS$cancer_types, useNA = "ifany")

# ---------- Stage: unmappable values and Stage 0 become Unknown
data_DSS$stage = trimws(as.character(data_DSS$stage))
data_DSS$stage[data_DSS$stage %in% c("", "NA", "[Discrepancy]", "[Unknown]",
                                     "Stage X", "X", "Unknown", "I/II NOS",
                                     "Stage 0")] = "Unknown"
data_DSS$stage[is.na(data_DSS$stage)] = "Unknown"
data_DSS$stage = factor(data_DSS$stage,
                         levels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Unknown"))
table(data_DSS$stage, useNA = "ifany")

# ---------- WSS subtype
data_DSS$level = as.character(data_DSS$cluster)
data_DSS$level[is.na(data_DSS$level) | data_DSS$level == ""] = "Unknown"
data_DSS$level = factor(data_DSS$level,
                         levels = c("panWSS4", "panWSS1", "panWSS2", "panWSS3", "Unknown"))

# ---------- Demographics and clinical variables: missing becomes Unknown
data_DSS$age = as.character(data_DSS$age)
data_DSS$age[is.na(data_DSS$age) | data_DSS$age == ""] = "Unknown"
data_DSS$age = factor(data_DSS$age,
                       levels = c("age<=40", "40<age<=60", "age>60", "Unknown"))

data_DSS$gender = as.character(data_DSS$gender)
data_DSS$gender[is.na(data_DSS$gender) | data_DSS$gender == ""] = "Unknown"
data_DSS$gender = factor(data_DSS$gender, levels = c("Female", "Male", "Unknown"))

if ("race" %in% colnames(data_DSS)) {
  data_DSS$race = as.character(data_DSS$race)
  data_DSS$race[data_DSS$race %in% c("", "Other/Unknown")] = "Unknown"
  data_DSS$race[is.na(data_DSS$race)] = "Unknown"
  data_DSS$race = factor(data_DSS$race, levels = c("White", "Black", "Asian", "Unknown"))
}

if ("outcome" %in% colnames(data_DSS)) {
  data_DSS$outcome = as.character(data_DSS$outcome)
  data_DSS$outcome[data_DSS$outcome %in% c("", "NA")] = "Unknown"
  data_DSS$outcome[is.na(data_DSS$outcome)] = "Unknown"
  data_DSS$outcome = factor(data_DSS$outcome,
                             levels = c("CR/PR", "SD/PD", "Unknown"))
}

# ---------- Mutation variables: APC only
gene_list = c("APC")
for (g in gene_list) {
  if (g %in% colnames(data_DSS)) {
    data_DSS[[g]] = as.character(data_DSS[[g]])
    data_DSS[[g]][data_DSS[[g]] %in% c("", "NA")] = "Unknown"
    data_DSS[[g]][is.na(data_DSS[[g]])] = "Unknown"
    data_DSS[[g]] = factor(data_DSS[[g]], levels = c("WT", "MT", "Unknown"))
  }
}

# ---------- Covariates and the Cox model (tumor_status and the other genes removed)
all_covariates = c("level", "cancer_types", "stage",
                    "age", "gender", "race", "outcome", gene_list)
all_covariates = all_covariates[all_covariates %in% colnames(data_DSS)]

# Drop factor levels that have no samples
for (v in all_covariates) {
  if (is.factor(data_DSS[[v]])) data_DSS[[v]] = droplevels(data_DSS[[v]])
}

cox_formula = as.formula(paste("Surv(DSS.time, DSS) ~", paste(all_covariates, collapse = " + ")))
fit_cox = coxph(cox_formula, data = data_DSS)

res_tidy = tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)

# ---------- Display-name mapping
var_labels = c(
  "level"        = "panWSS labels",
  "cancer_types" = "Cancer types",
  "stage"        = "Stage",
  "gender"       = "Gender",
  "age"          = "Age",
  "race"         = "Race",
  "outcome"      = "Outcome"
)
for (g in gene_list) {
  var_labels[g] = paste0(g, " mutation")
}

# Readable labels for the 12 cancer-type levels
cancer_level_labels = c(
  "Gynecologic"                      = "Gynecologic",
  "Urologic"                         = "Urologic",
  "Endocrine"                        = "Endocrine",
  "Core_gastrointestinal"            = "Core gastrointestinal",
  "Developmental_gastrointestinal"   = "Developmental gastrointestinal",
  "Head_and_neck"                    = "Head and neck",
  "Thoracic"                         = "Thoracic",
  "Central_nervous_system"           = "Central nervous system",
  "Soft_tissue"                      = "Soft tissue",
  "Neural_crest_derived"             = "Neural-crest-derived",
  "Melanocytic"                      = "Melanocytic",
  "Hematologic_and_lymphatic"        = "Hematologic and lymphatic",
  "Unknown"                          = "Unknown"
)

factors_in_model = all_covariates
plot_rows = list()

# ---------- Build the forest-plot table rows
for (v in factors_in_model) {
  if (!v %in% names(data_DSS)) next
  lvls = levels(data_DSS[[v]])
  if (is.null(lvls)) next

  display_var_name = ifelse(v %in% names(var_labels), var_labels[v], v)

  plot_rows[[paste0(v, "_header")]] = data.frame(
    col1 = display_var_name, col2 = "", col3 = "", col4 = "",
    col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
    HR = NA, lower = NA, upper = NA, stringsAsFactors = FALSE
  )

  for (i in seq_along(lvls)) {
    l = lvls[i]
    n_val = sum(data_DSS[[v]] == l, na.rm = TRUE)

    if (n_val == 0) next

    l_display = if (v == "cancer_types" && l %in% names(cancer_level_labels)) {
      cancer_level_labels[[l]]
    } else {
      l
    }

    if (i == 1) {
      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = "Reference", col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
        HR = 1.0, lower = 1.0, upper = 1.0, stringsAsFactors = FALSE
      )
    } else {
      term_name = paste0(v, l)
      sub_res = res_tidy %>% filter(term == term_name)

      hr_str = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1])) {
        sprintf("%.2f (%.2f–%.2f)", sub_res$estimate[1], sub_res$conf.low[1], sub_res$conf.high[1])
      } else ""

      p_str = if (nrow(sub_res) > 0 && !is.na(sub_res$p.value[1])) {
        if (sub_res$p.value[1] < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value[1])
      } else ""

      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = hr_str, col5 = paste(rep(" ", 25), collapse = " "), col6 = p_str,
        HR = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1]))  sub_res$estimate[1]  else NA,
        lower = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.low[1]))  sub_res$conf.low[1]  else NA,
        upper = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.high[1])) sub_res$conf.high[1] else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_forest_input = bind_rows(plot_rows)
colnames(df_forest_input)[1:6] = c("Variable", "Subgroup", "N",
                                    "Hazard ratio (95% CI)", "", "p-value")

# ---------- Render and export the DSS forest plot (APC only, HR axis max 11)
pdf("C:/Users/user/Desktop/Multivariate_Cox_DSS_APC_0_11_Forestplot.pdf", width = 20, height = 24)

p_final = forest(
  df_forest_input[, 1:6],
  est = df_forest_input$HR,
  lower = df_forest_input$lower,
  upper = df_forest_input$upper,
  size = 0.5,
  ci_column = 5,
  ref_line = 1,
  xlim = c(0, 11),
  ticks_at = c(0, 2, 4, 6, 8, 10, 11),
  xlab = "Hazard ratio (95% CI)"
)

print(p_final)
dev.off()

print(p_final)
cat("\nMultivariable DSS Cox forest plot (APC only, HR axis max 11) saved to the desktop.\n")

################################################################################Multivariable visualisation of DFI####
rm(list = ls())
gc()

library(survival)
library(dplyr)
library(broom)
library(forestploter)

setwd("C:/Users/user/Desktop/input")
load("data_DFI.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")

data_DFI = merge(data_DFI, pancancer_group, by.x = "Tag", by.y = "row.names")

# Keep tumour samples only: drop samples whose numeric suffix is > 10
data_DFI = data_DFI %>%
  mutate(sample_code = as.numeric(sub(".*-", "", Tag))) %>%
  filter(sample_code <= 10) %>%
  select(-sample_code)

# ---------- Cancer type: 12 organ-system groups
data_DFI = data_DFI %>%
  mutate(
    cancer_types = case_when(
      is.na(Project) ~ "Unknown",
      Project %in% c("OV", "UCEC", "CESC", "BRCA") ~ "Gynecologic",
      Project %in% c("BLCA", "PRAD", "TGCT", "KIRC", "KICH", "KIRP") ~ "Urologic",
      Project %in% c("THCA", "ACC") ~ "Endocrine",
      Project %in% c("ESCA", "STAD", "COAD", "READ") ~ "Core_gastrointestinal",
      Project %in% c("LIHC", "PAAD", "CHOL") ~ "Developmental_gastrointestinal",
      Project %in% c("HNSC") ~ "Head_and_neck",
      Project %in% c("LUAD", "LUSC", "MESO") ~ "Thoracic",
      Project %in% c("GBM", "LGG") ~ "Central_nervous_system",
      Project %in% c("SARC", "UCS") ~ "Soft_tissue",
      Project %in% c("PCPG") ~ "Neural_crest_derived",
      Project %in% c("SKCM", "UVM") ~ "Melanocytic",
      Project %in% c("LAML", "DLBC", "THYM") ~ "Hematologic_and_lymphatic",
      TRUE ~ "Unknown"
    )
  )
data_DFI$cancer_types = as.character(data_DFI$cancer_types)
data_DFI$cancer_types[is.na(data_DFI$cancer_types) | data_DFI$cancer_types == ""] = "Unknown"
data_DFI$cancer_types = factor(data_DFI$cancer_types,
                                levels = c("Gynecologic", "Urologic", "Endocrine",
                                           "Core_gastrointestinal", "Developmental_gastrointestinal",
                                           "Head_and_neck", "Thoracic", "Central_nervous_system",
                                           "Soft_tissue", "Neural_crest_derived", "Melanocytic",
                                           "Hematologic_and_lymphatic", "Unknown"))
table(data_DFI$cancer_types, useNA = "ifany")

# ---------- Stage: unmappable values and Stage 0 become Unknown
data_DFI$stage = trimws(as.character(data_DFI$stage))
data_DFI$stage[data_DFI$stage %in% c("", "NA", "[Discrepancy]", "[Unknown]",
                                     "Stage X", "X", "Unknown", "I/II NOS",
                                     "Stage 0")] = "Unknown"
data_DFI$stage[is.na(data_DFI$stage)] = "Unknown"
data_DFI$stage = factor(data_DFI$stage,
                         levels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Unknown"))
table(data_DFI$stage, useNA = "ifany")

# ---------- WSS subtype
data_DFI$level = as.character(data_DFI$cluster)
data_DFI$level[is.na(data_DFI$level) | data_DFI$level == ""] = "Unknown"
data_DFI$level = factor(data_DFI$level,
                         levels = c("panWSS4", "panWSS1", "panWSS2", "panWSS3", "Unknown"))

# ---------- Demographics and clinical variables: missing becomes Unknown
data_DFI$age = as.character(data_DFI$age)
data_DFI$age[is.na(data_DFI$age) | data_DFI$age == ""] = "Unknown"
data_DFI$age = factor(data_DFI$age,
                       levels = c("age<=40", "40<age<=60", "age>60", "Unknown"))

data_DFI$gender = as.character(data_DFI$gender)
data_DFI$gender[is.na(data_DFI$gender) | data_DFI$gender == ""] = "Unknown"
data_DFI$gender = factor(data_DFI$gender, levels = c("Female", "Male", "Unknown"))

if ("race" %in% colnames(data_DFI)) {
  data_DFI$race = as.character(data_DFI$race)
  data_DFI$race[data_DFI$race %in% c("", "Other/Unknown")] = "Unknown"
  data_DFI$race[is.na(data_DFI$race)] = "Unknown"
  data_DFI$race = factor(data_DFI$race, levels = c("White", "Black", "Asian", "Unknown"))
}

if ("outcome" %in% colnames(data_DFI)) {
  data_DFI$outcome = as.character(data_DFI$outcome)
  data_DFI$outcome[data_DFI$outcome %in% c("", "NA")] = "Unknown"
  data_DFI$outcome[is.na(data_DFI$outcome)] = "Unknown"
  data_DFI$outcome = factor(data_DFI$outcome,
                             levels = c("CR/PR", "SD/PD", "Unknown"))
}

# ---------- Mutation variables: APC only
gene_list = c("APC")
for (g in gene_list) {
  if (g %in% colnames(data_DFI)) {
    data_DFI[[g]] = as.character(data_DFI[[g]])
    data_DFI[[g]][data_DFI[[g]] %in% c("", "NA")] = "Unknown"
    data_DFI[[g]][is.na(data_DFI[[g]])] = "Unknown"
    data_DFI[[g]] = factor(data_DFI[[g]], levels = c("WT", "MT", "Unknown"))
  }
}

# ---------- Covariates and the Cox model (tumor_status and the other genes removed)
all_covariates = c("level", "cancer_types", "stage",
                    "age", "gender", "race", "outcome", gene_list)
all_covariates = all_covariates[all_covariates %in% colnames(data_DFI)]

# Drop factor levels that have no samples
for (v in all_covariates) {
  if (is.factor(data_DFI[[v]])) data_DFI[[v]] = droplevels(data_DFI[[v]])
}

cox_formula = as.formula(paste("Surv(DFI.time, DFI) ~", paste(all_covariates, collapse = " + ")))
fit_cox = coxph(cox_formula, data = data_DFI)

res_tidy = tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)

# ---------- Display-name mapping
var_labels = c(
  "level"        = "panWSS labels",
  "cancer_types" = "Cancer types",
  "stage"        = "Stage",
  "gender"       = "Gender",
  "age"          = "Age",
  "race"         = "Race",
  "outcome"      = "Outcome"
)
for (g in gene_list) {
  var_labels[g] = paste0(g, " mutation")
}

# Readable labels for the 12 cancer-type levels
cancer_level_labels = c(
  "Gynecologic"                      = "Gynecologic",
  "Urologic"                         = "Urologic",
  "Endocrine"                        = "Endocrine",
  "Core_gastrointestinal"            = "Core gastrointestinal",
  "Developmental_gastrointestinal"   = "Developmental gastrointestinal",
  "Head_and_neck"                    = "Head and neck",
  "Thoracic"                         = "Thoracic",
  "Central_nervous_system"           = "Central nervous system",
  "Soft_tissue"                      = "Soft tissue",
  "Neural_crest_derived"             = "Neural-crest-derived",
  "Melanocytic"                      = "Melanocytic",
  "Hematologic_and_lymphatic"        = "Hematologic and lymphatic",
  "Unknown"                          = "Unknown"
)

factors_in_model = all_covariates
plot_rows = list()

# ---------- Build the forest-plot table rows
for (v in factors_in_model) {
  if (!v %in% names(data_DFI)) next
  lvls = levels(data_DFI[[v]])
  if (is.null(lvls)) next

  display_var_name = ifelse(v %in% names(var_labels), var_labels[v], v)

  plot_rows[[paste0(v, "_header")]] = data.frame(
    col1 = display_var_name, col2 = "", col3 = "", col4 = "",
    col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
    HR = NA, lower = NA, upper = NA, stringsAsFactors = FALSE
  )

  for (i in seq_along(lvls)) {
    l = lvls[i]
    n_val = sum(data_DFI[[v]] == l, na.rm = TRUE)

    if (n_val == 0) next

    l_display = if (v == "cancer_types" && l %in% names(cancer_level_labels)) {
      cancer_level_labels[[l]]
    } else {
      l
    }

    if (i == 1) {
      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = "Reference", col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
        HR = 1.0, lower = 1.0, upper = 1.0, stringsAsFactors = FALSE
      )
    } else {
      term_name = paste0(v, l)
      sub_res = res_tidy %>% filter(term == term_name)

      hr_str = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1])) {
        sprintf("%.2f (%.2f–%.2f)", sub_res$estimate[1], sub_res$conf.low[1], sub_res$conf.high[1])
      } else ""

      p_str = if (nrow(sub_res) > 0 && !is.na(sub_res$p.value[1])) {
        if (sub_res$p.value[1] < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value[1])
      } else ""

      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = hr_str, col5 = paste(rep(" ", 25), collapse = " "), col6 = p_str,
        HR = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1]))  sub_res$estimate[1]  else NA,
        lower = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.low[1]))  sub_res$conf.low[1]  else NA,
        upper = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.high[1])) sub_res$conf.high[1] else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_forest_input = bind_rows(plot_rows)
colnames(df_forest_input)[1:6] = c("Variable", "Subgroup", "N",
                                    "Hazard ratio (95% CI)", "", "p-value")

# ---------- Render and export the DFI forest plot (APC only, HR axis max 8)
pdf("C:/Users/user/Desktop/Multivariate_Cox_DFI_APC_0_8_Forestplot.pdf", width = 20, height = 24)

p_final = forest(
  df_forest_input[, 1:6],
  est = df_forest_input$HR,
  lower = df_forest_input$lower,
  upper = df_forest_input$upper,
  size = 0.5,
  ci_column = 5,
  ref_line = 1,
  xlim = c(0, 8),
  ticks_at = c(0, 2, 4, 6, 8),
  xlab = "Hazard ratio (95% CI)"
)

print(p_final)
dev.off()

print(p_final)
cat("\nMultivariable DFI Cox forest plot (APC only, HR axis max 8) saved to the desktop.\n")

################################################################################Multivariable visualisation of PFI####
rm(list = ls())
gc()

library(survival)
library(dplyr)
library(broom)
library(forestploter)

setwd("C:/Users/user/Desktop/input")
load("data_PFI.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")

data_PFI = merge(data_PFI, pancancer_group, by.x = "Tag", by.y = "row.names")

# Keep tumour samples only: drop samples whose numeric suffix is > 10
data_PFI = data_PFI %>%
  mutate(sample_code = as.numeric(sub(".*-", "", Tag))) %>%
  filter(sample_code <= 10) %>%
  select(-sample_code)

# ---------- Cancer type: 12 organ-system groups
data_PFI = data_PFI %>%
  mutate(
    cancer_types = case_when(
      is.na(Project) ~ "Unknown",
      Project %in% c("OV", "UCEC", "CESC", "BRCA") ~ "Gynecologic",
      Project %in% c("BLCA", "PRAD", "TGCT", "KIRC", "KICH", "KIRP") ~ "Urologic",
      Project %in% c("THCA", "ACC") ~ "Endocrine",
      Project %in% c("ESCA", "STAD", "COAD", "READ") ~ "Core_gastrointestinal",
      Project %in% c("LIHC", "PAAD", "CHOL") ~ "Developmental_gastrointestinal",
      Project %in% c("HNSC") ~ "Head_and_neck",
      Project %in% c("LUAD", "LUSC", "MESO") ~ "Thoracic",
      Project %in% c("GBM", "LGG") ~ "Central_nervous_system",
      Project %in% c("SARC", "UCS") ~ "Soft_tissue",
      Project %in% c("PCPG") ~ "Neural_crest_derived",
      Project %in% c("SKCM", "UVM") ~ "Melanocytic",
      Project %in% c("LAML", "DLBC", "THYM") ~ "Hematologic_and_lymphatic",
      TRUE ~ "Unknown"
    )
  )
data_PFI$cancer_types = as.character(data_PFI$cancer_types)
data_PFI$cancer_types[is.na(data_PFI$cancer_types) | data_PFI$cancer_types == ""] = "Unknown"
data_PFI$cancer_types = factor(data_PFI$cancer_types,
                                levels = c("Gynecologic", "Urologic", "Endocrine",
                                           "Core_gastrointestinal", "Developmental_gastrointestinal",
                                           "Head_and_neck", "Thoracic", "Central_nervous_system",
                                           "Soft_tissue", "Neural_crest_derived", "Melanocytic",
                                           "Hematologic_and_lymphatic", "Unknown"))
table(data_PFI$cancer_types, useNA = "ifany")

# ---------- Stage: unmappable values and Stage 0 become Unknown
data_PFI$stage = trimws(as.character(data_PFI$stage))
data_PFI$stage[data_PFI$stage %in% c("", "NA", "[Discrepancy]", "[Unknown]",
                                     "Stage X", "X", "Unknown", "I/II NOS",
                                     "Stage 0")] = "Unknown"
data_PFI$stage[is.na(data_PFI$stage)] = "Unknown"
data_PFI$stage = factor(data_PFI$stage,
                         levels = c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Unknown"))
table(data_PFI$stage, useNA = "ifany")

# ---------- WSS subtype, with panWSS4 as the reference
data_PFI$level = as.character(data_PFI$cluster)
data_PFI$level[is.na(data_PFI$level) | data_PFI$level == ""] = "Unknown"
data_PFI$level = factor(data_PFI$level,
                         levels = c("panWSS4", "panWSS1", "panWSS2", "panWSS3", "Unknown"))

# ---------- Demographics and clinical variables: missing becomes Unknown
data_PFI$age = as.character(data_PFI$age)
data_PFI$age[is.na(data_PFI$age) | data_PFI$age == ""] = "Unknown"
data_PFI$age = factor(data_PFI$age,
                       levels = c("age<=40", "40<age<=60", "age>60", "Unknown"))

data_PFI$gender = as.character(data_PFI$gender)
data_PFI$gender[is.na(data_PFI$gender) | data_PFI$gender == ""] = "Unknown"
data_PFI$gender = factor(data_PFI$gender, levels = c("Female", "Male", "Unknown"))

if ("race" %in% colnames(data_PFI)) {
  data_PFI$race = as.character(data_PFI$race)
  data_PFI$race[data_PFI$race %in% c("", "Other/Unknown")] = "Unknown"
  data_PFI$race[is.na(data_PFI$race)] = "Unknown"
  data_PFI$race = factor(data_PFI$race, levels = c("White", "Black", "Asian", "Unknown"))
}

if ("outcome" %in% colnames(data_PFI)) {
  data_PFI$outcome = as.character(data_PFI$outcome)
  data_PFI$outcome[data_PFI$outcome %in% c("", "NA")] = "Unknown"
  data_PFI$outcome[is.na(data_PFI$outcome)] = "Unknown"
  data_PFI$outcome = factor(data_PFI$outcome,
                             levels = c("CR/PR", "SD/PD", "Unknown"))
}

# ---------- Mutation variables: APC only
gene_list = c("APC")
for (g in gene_list) {
  if (g %in% colnames(data_PFI)) {
    data_PFI[[g]] = as.character(data_PFI[[g]])
    data_PFI[[g]][data_PFI[[g]] %in% c("", "NA")] = "Unknown"
    data_PFI[[g]][is.na(data_PFI[[g]])] = "Unknown"
    data_PFI[[g]] = factor(data_PFI[[g]], levels = c("WT", "MT", "Unknown"))
  }
}

# ---------- Covariates and the Cox model (tumor_status and the other genes removed)
all_covariates = c("level", "cancer_types", "stage",
                    "age", "gender", "race", "outcome", gene_list)
all_covariates = all_covariates[all_covariates %in% colnames(data_PFI)]

# Drop factor levels that have no samples
for (v in all_covariates) {
  if (is.factor(data_PFI[[v]])) data_PFI[[v]] = droplevels(data_PFI[[v]])
}

cox_formula = as.formula(paste("Surv(PFI.time, PFI) ~", paste(all_covariates, collapse = " + ")))
fit_cox = coxph(cox_formula, data = data_PFI)

res_tidy = tidy(fit_cox, exponentiate = TRUE, conf.int = TRUE)

# ---------- Display-name mapping
var_labels = c(
  "level"        = "panWSS labels",
  "cancer_types" = "Cancer types",
  "stage"        = "Stage",
  "gender"       = "Gender",
  "age"          = "Age",
  "race"         = "Race",
  "outcome"      = "Outcome"
)
for (g in gene_list) {
  var_labels[g] = paste0(g, " mutation")
}

# Readable labels for the 12 cancer-type levels
cancer_level_labels = c(
  "Gynecologic"                      = "Gynecologic",
  "Urologic"                         = "Urologic",
  "Endocrine"                        = "Endocrine",
  "Core_gastrointestinal"            = "Core gastrointestinal",
  "Developmental_gastrointestinal"   = "Developmental gastrointestinal",
  "Head_and_neck"                    = "Head and neck",
  "Thoracic"                         = "Thoracic",
  "Central_nervous_system"           = "Central nervous system",
  "Soft_tissue"                      = "Soft tissue",
  "Neural_crest_derived"             = "Neural-crest-derived",
  "Melanocytic"                      = "Melanocytic",
  "Hematologic_and_lymphatic"        = "Hematologic and lymphatic",
  "Unknown"                          = "Unknown"
)

factors_in_model = all_covariates
plot_rows = list()

# ---------- Build the forest-plot table rows
for (v in factors_in_model) {
  if (!v %in% names(data_PFI)) next
  lvls = levels(data_PFI[[v]])
  if (is.null(lvls)) next

  display_var_name = ifelse(v %in% names(var_labels), var_labels[v], v)

  plot_rows[[paste0(v, "_header")]] = data.frame(
    col1 = display_var_name, col2 = "", col3 = "", col4 = "",
    col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
    HR = NA, lower = NA, upper = NA, stringsAsFactors = FALSE
  )

  for (i in seq_along(lvls)) {
    l = lvls[i]
    n_val = sum(data_PFI[[v]] == l, na.rm = TRUE)

    if (n_val == 0) next

    l_display = if (v == "cancer_types" && l %in% names(cancer_level_labels)) {
      cancer_level_labels[[l]]
    } else {
      l
    }

    if (i == 1) {
      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = "Reference", col5 = paste(rep(" ", 25), collapse = " "), col6 = "",
        HR = 1.0, lower = 1.0, upper = 1.0, stringsAsFactors = FALSE
      )
    } else {
      term_name = paste0(v, l)
      sub_res = res_tidy %>% filter(term == term_name)

      hr_str = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1])) {
        sprintf("%.2f (%.2f–%.2f)", sub_res$estimate[1], sub_res$conf.low[1], sub_res$conf.high[1])
      } else ""

      p_str = if (nrow(sub_res) > 0 && !is.na(sub_res$p.value[1])) {
        if (sub_res$p.value[1] < 0.001) "<0.001" else sprintf("%.3f", sub_res$p.value[1])
      } else ""

      plot_rows[[paste0(v, "_", l)]] = data.frame(
        col1 = "", col2 = paste0("   ", l_display), col3 = as.character(n_val),
        col4 = hr_str, col5 = paste(rep(" ", 25), collapse = " "), col6 = p_str,
        HR = if (nrow(sub_res) > 0 && !is.na(sub_res$estimate[1]))  sub_res$estimate[1]  else NA,
        lower = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.low[1]))  sub_res$conf.low[1]  else NA,
        upper = if (nrow(sub_res) > 0 && !is.na(sub_res$conf.high[1])) sub_res$conf.high[1] else NA,
        stringsAsFactors = FALSE
      )
    }
  }
}

df_forest_input = bind_rows(plot_rows)
colnames(df_forest_input)[1:6] = c("Variable", "Subgroup", "N",
                                    "Hazard ratio (95% CI)", "", "p-value")

# ---------- Render and export the PFI forest plot (APC only, HR axis max 7)
pdf("C:/Users/user/Desktop/Multivariate_Cox_PFI_APC_0_7_Forestplot.pdf", width = 20, height = 24)

p_final = forest(
  df_forest_input[, 1:6],
  est = df_forest_input$HR,
  lower = df_forest_input$lower,
  upper = df_forest_input$upper,
  size = 0.5,
  ci_column = 5,
  ref_line = 1,
  xlim = c(0, 7),
  ticks_at = c(0, 1, 2, 3, 4, 5, 6, 7),
  xlab = "Hazard ratio (95% CI)"
)

print(p_final)
dev.off()

print(p_final)
cat("\nMultivariable PFI Cox forest plot (APC only, HR axis max 7) saved to the desktop.\n")
