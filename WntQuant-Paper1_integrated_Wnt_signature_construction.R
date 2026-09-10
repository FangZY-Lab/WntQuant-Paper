# ==============================================================================
# WntQuant-Paper1 -- Integrated Wnt signature construction
# ------------------------------------------------------------------------------
# Purpose:    Construct the integrated, direction-aware Wnt pathway signature (WntQuant-Paper score) from curated activation (WPAGS) and inhibition (WPIGS) gene sets, with PCA-based training-data quality control.
# Inputs:     Expression matrices and group annotations for 24 training cohorts.
# Outputs:    Directional gene-set GMT files and per-sample Wnt activity scores.
# NOTE:       Part of the code involves the core interests of the laboratory and is not open source at this time.
# ==============================================================================


rm(list=ls())
gc()
setwd("input")
list=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
       "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
       "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
       "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
       "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(list)) {
  load(paste0(list[q],".Rdata"))
  load(paste0(list[q],"_G.Rdata"))
}
library(ggplot2)
library(plyr)
for (i in 1:length(list)) {
  exp=get(list[i])
  pca_data=as.data.frame(exp)
  group=get(paste0(list[i],"_G"))
  Group=as.data.frame(group)
  pca_data = scale(pca_data)
  pca_data = t(scale(t(pca_data),center=TRUE,scale=F))
  pca = prcomp(t(pca_data),center=FALSE,scale.=FALSE)
  PCA_sum = summary(pca)
  PC12 = pca$x[,1:2]
  pc = PCA_sum$importance[2,]*100
  PC12 = as.data.frame(PC12)
  PC12$samples = row.names(PC12)
  colnames(Group) = c("samples","group")
  df = merge(PC12,Group,by="samples")
  color=c("#f58f98","#90d7ec")
  p1 = ggplot(data = df, aes(x = PC1, y = PC2, color = group, shape = group)) +  
    ggtitle(paste0(list[i])) +  
    theme_bw() +  
    geom_point(size = 3) +  
    theme(panel.grid = element_blank()) +  
    geom_vline(xintercept = 0, linetype = 'dotdash', size = 1, color = 'grey') +  
    geom_hline(yintercept = 0, linetype = 'dotdash', size = 1, color = 'grey') +  
    labs(x = paste0("PC1 (", pc[1], "%)"), y = paste0("PC2 (", pc[2], "%)")) +  
    scale_color_manual(values = color) +  
    theme(axis.title.x = element_text(size = 16,color = "black"),  
          axis.title.y = element_text(size = 16, angle = 90,color = "black"),  
          axis.text.y = element_text(size = 16,color = "black"),  
          axis.text.x = element_text(size = 16,color = "black"),  
          plot.title = element_text(size = 20,color = "black"),  
          legend.title = element_text(size = 14,color = "black"),  
          legend.text = element_text(size = 12,color = "black")) 
  library(plyr)
  cluster_border = ddply(df, 'group', function(df) df[chull(df[[2]], df[[3]]), ])
  p=p1 + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = c("#f58f98","#90d7ec"))
  p
  ggsave(paste0(list[i],"_PCA",".pdf"), p, width = 5, height = 3)
}
# ------------------------------------------------------------------------------
# WNT gene set preparation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
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
# ------------------------------------------------------------------------------
# Add enrichr's gene set
# ------------------------------------------------------------------------------
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
write.gmt(all_genesets,file="input/all_genesets.gmt")
# ------------------------------------------------------------------------------
# Cleaning the gene set
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
vector=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
         "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
         "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
         "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
         "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}
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
up.p.matrix=UP.clean_genesets_rank$cleaning_system
save(up.p.matrix,file="up.p.matrix.Rdata")
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
dn.p.matrix=DN.clean_genesets_rank$cleaning_system
save(dn.p.matrix,file="dn.p.matrix.Rdata")
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
colnames(wpags_wpigs)=c(paste0("WPAGS", 1:length(colnames(activation_geneset))), 
                        paste0("WPIGS", 1:length(colnames(inhibition_geneset))))
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
write.gmt(wpags_wpigs,file="input/wpags_wpigs.gmt")
save(joint_genesets,file="joint_genesets.Rdata")
# ------------------------------------------------------------------------------
# Mining of significant genes
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(dplyr)
setwd("input")
load("input/up.p.matrix.Rdata")
load("input/dn.p.matrix.Rdata")
up.p.matrix=up.p.matrix[1:150,length(colnames(up.p.matrix)),drop=F]
up.p.matrix$gene=rownames(up.p.matrix)
dn.p.matrix=dn.p.matrix[1:150,length(colnames(dn.p.matrix)),drop=F]
dn.p.matrix$gene=rownames(dn.p.matrix)
library(ggplot2)
ggplot(up.p.matrix, aes(x = reorder(gene, -p_score), y = p_score, fill = p_score)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) + 
  labs(
    title = "",
    subtitle = "",
    x = "Activation genes (top 150)",  
    y = expression(-log[10](italic(p)[UP]) * "-combined")
  ) +
  theme_minimal(base_size = 16) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, color = "black"), 
    axis.text = element_text(size = 16, color = "black"),  
    axis.title = element_text(size = 30, color = "black"),  
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 0, hjust = 1, color = "black")
  ) +  
  scale_fill_gradient(low = "#fffffb", high = "#f58f98") +
  coord_flip()
ggplot(dn.p.matrix, aes(x = reorder(gene, -p_score), y = p_score, fill = p_score)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) + 
  labs(
    title = "",
    subtitle = "",
    x = "Inhibition genes (top 150)",  
    y = expression(-log[10](italic(p)[DN]) * "-combined")
  ) +
  theme_minimal(base_size = 16) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, color = "black"), 
    axis.text = element_text(size = 16, color = "black"),  
    axis.title = element_text(size = 30, color = "black"),  
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 0, hjust = 1, color = "black")
  ) +  
  scale_fill_gradient(low = "#fffffb", high = "#90d7ec") +
  coord_flip()
# ------------------------------------------------------------------------------
# fGSEA training data
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
library(limma)
library(dplyr)
setwd("input")
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
# ------------------------------------------------------------------------------
# fGSEA validation data
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
library(limma)
library(dplyr)
setwd("input")
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
# ------------------------------------------------------------------------------
# Training dataset validation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
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
  group$Tag=paste0(group$Tag,"_",vector[i])
  if(i==1){
    all_group=group
  }else{
    all_group=rbind(all_group,group)
  }
}
all_group=all_group
save(all_group,file="all_group.Rdata")
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="input",
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
load("input/all_group.Rdata")
drawed_data=merge(all_scores,all_group,by.x="row.names",by.y="Tag")
colnames(drawed_data)[6]="Group"
library(ggplot2)
library(ggpubr)
p <- ggplot(drawed_data, aes(fill=Group, y=activity_score, x=ID)) +
  geom_bar(position=position_dodge(),
           stat="summary",
           width=0.7,
           colour = "black",
           linewidth=0.5) +
  theme_bw() +
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +
  theme(legend.position = "top",
        panel.grid.major.x = element_line(linewidth = 0.2),
        panel.grid.major.y = element_line(linewidth = 0.2),
        panel.grid.minor.x = element_line(linewidth = 0.2),
        panel.grid.minor.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(colour = 'black', linewidth = 0.2),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_text(colour = "black", size = 16),
        axis.title.y = element_text(colour = "black", size = 16),
        legend.text = element_text(colour = "black", size = 16),
        legend.title = element_text(colour = "black", size = 16),
        plot.title = element_text(colour = "black", size = 16)) +
  stat_summary(fun.data = 'mean_se', 
               geom = "errorbar", 
               colour = "black",
               width = 0.3,
               position=position_dodge(0.7)) +
  scale_fill_manual(values = c("#f58f98","#90d7ec")) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method="t.test") +
  ggtitle("Training datasets")
p
# ------------------------------------------------------------------------------
# Testing dataset validation
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
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
all_group_v=all_group
save(all_group_v,file="all_group_v.Rdata")
all_scores=data.frame()
geneSets=getGmt("wpags_wpigs.gmt")
source("functions/Calculate_bioactivity_scores.R")
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="input",
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
load("input/all_group_v.Rdata")
drawed_data=merge(all_scores,all_group_v,by.x="row.names",by.y="Tag")
colnames(drawed_data)[6]="Group"
library(ggplot2)
library(ggpubr)
p <- ggplot(drawed_data, aes(fill=Group, y=activity_score, x=ID)) +
  geom_bar(position=position_dodge(),
           stat="summary",
           width=0.7,
           colour = "black",
           linewidth=0.5) +
  theme_bw() +
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +
  theme(legend.position = "top",
        panel.grid.major.x = element_line(linewidth = 0.2),
        panel.grid.major.y = element_line(linewidth = 0.2),
        panel.grid.minor.x = element_line(linewidth = 0.2),
        panel.grid.minor.y = element_line(linewidth = 0.2),
        axis.line.x = element_line(colour = 'black', linewidth = 0.2),
        axis.line.y = element_line(colour = 'black', linewidth = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        axis.title.x = element_text(colour = "black", size = 16),
        axis.title.y = element_text(colour = "black", size = 16),
        legend.text = element_text(colour = "black", size = 16),
        legend.title = element_text(colour = "black", size = 16),
        plot.title = element_text(colour = "black", size = 16)) +
  stat_summary(fun.data = 'mean_se', 
               geom = "errorbar", 
               colour = "black",
               width = 0.3,
               position=position_dodge(0.7)) +
  scale_fill_manual(values = c("#f58f98","#90d7ec")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  stat_compare_means(aes(group = Group, label = ..p.signif..), method="t.test") +
  ggtitle("Testing datasets")+
  coord_flip() 
p
# ------------------------------------------------------------------------------
# PPI network
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
setwd("input")
library(clusterProfiler)
activation_inhibition_geneset=read.gmt("wpags_wpigs.gmt")
colnames(activation_inhibition_geneset)=c("group","gene")
activation_inhibition_geneset=activation_inhibition_geneset[,c("gene","group")]
activation_inhibition_geneset$group = substr(activation_inhibition_geneset$group, 1, 5)  
activation_inhibition_geneset = unique(activation_inhibition_geneset)
library(tidyverse)
library(STRINGdb)
string_db=STRINGdb$new( version="12.0",
                        species=9606,
                        score_threshold=600,
                        input_directory="input")
data_frame=data.frame(gene = activation_inhibition_geneset$gene)
dat_map=string_db$map(my_data_frame=data_frame, 
                      my_data_frame_id_col_names="gene",
                      removeUnmappedRows=TRUE)
hits=dat_map$STRING_id 
dat_link=string_db$get_interactions(hits)
dat_link$from=dat_map[match(dat_link$from,dat_map$STRING_id),'gene']
dat_link$to=dat_map[match(dat_link$to,dat_map$STRING_id),'gene'] 
colnames(dat_link)=c('node1','node2','combined_score')
dat_link=dat_link %>% distinct(node1, node2, .keep_all = T)
dat_link=dat_link[dat_link$node1%in%activation_inhibition_geneset$gene,]
dat_link=dat_link[dat_link$node2%in%activation_inhibition_geneset$gene,]
write.csv(dat_link,'string_link.csv',row.names = F,quote = F)
write.csv(activation_inhibition_geneset,'group_PPI.csv',row.names = F,quote = F)
# ------------------------------------------------------------------------------
# Comparison of other signatures in the training datasets
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
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
source("functions/Calculate_bioactivity_scores.R")
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
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
  scores=scores$final_activity_score
  scores$ID=paste0(vector[cycle])
  all_scores=rbind(all_scores,scores)
}
all_usual=data.frame()
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  setwd("input")
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
  gene_set=colnames(all_data_fs)[2:14]
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
gene_set=sub("activity_score", "Wnt/β-catenin pathway activity score", gene_set)
gene_set=sub("pos_score", "Wnt/β-catenin pathway activation score", gene_set)
gene_set=sub("neg_score", "Wnt/β-catenin pathway inhibition score", gene_set)
p_data_all$geneset=sub("activity_score", "Wnt/β-catenin pathway activity score", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "Wnt/β-catenin pathway activation score", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "Wnt/β-catenin pathway inhibition score", p_data_all$geneset)
p_data_all$geneset=factor(p_data_all$geneset, levels = gene_set)
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
# Comparison of other signatures in the testing datasets
# ------------------------------------------------------------------------------
rm(list=ls())
gc()
library(GSEABase)
setwd("input")
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
source("functions/Calculate_bioactivity_scores.R")
for(cycle in 1:length(vector)){
  exp=get(vector[cycle])
  colnames(exp)=paste0(colnames(exp),"_",vector[cycle])
  scores=Calculate_bioactivity_scores(file_paths="input",
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
  setwd("input")
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
  gene_set=colnames(all_data_fs)[2:14]
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
gene_set=sub("activity_score", "Wnt/β-catenin pathway activity score", gene_set)
gene_set=sub("pos_score", "Wnt/β-catenin pathway activation score", gene_set)
gene_set=sub("neg_score", "Wnt/β-catenin pathway inhibition score", gene_set)
p_data_all$geneset=sub("activity_score", "Wnt/β-catenin pathway activity score", p_data_all$geneset)
p_data_all$geneset=sub("pos_score", "Wnt/β-catenin pathway activation score", p_data_all$geneset)
p_data_all$geneset=sub("neg_score", "Wnt/β-catenin pathway inhibition score", p_data_all$geneset)
p_data_all$geneset=factor(p_data_all$geneset, levels = gene_set)
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