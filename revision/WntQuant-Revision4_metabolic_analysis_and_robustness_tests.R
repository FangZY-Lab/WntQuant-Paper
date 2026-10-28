# ==============================================================================
# WntQuant revision -- Metabolic programmes and robustness tests
# ------------------------------------------------------------------------------
# Purpose:    Link Wnt activity to TCGA metabolic programmes and quantify how
#             robust the gene sets are to gene and cohort dropout and to
#             added noise.
# Inputs:     TCGA metabolism and pan-cancer cluster objects under input/.
# Outputs:    Correlation panels, robustness curves and the redundancy-method
#             comparison.
# ==============================================================================

################################################################################
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/TCGA_metabolism.Rdata")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
library(limma)
library(dplyr)
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
#panWSS1
group_panWSS1=pancancer_cluster
group_panWSS1$cluster=ifelse(group_panWSS1$cluster=="panWSS1","H","L")
colnames(group_panWSS1)=c("Tag","group")
result=Get_limma(gene=TCGA_metabolism,
                 group=group_panWSS1,
                 alternative_limma="greater")

library(ggplot2)
library(ggrepel)
min_p_value_row = result[which.min(result$P.Value), ]
p1=ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value < 0.05,
                                ifelse(logFC > 0, "Up", "Down"),
                                "None")),
             alpha = 0.8, size = 3.8) +

  geom_point(data = min_p_value_row, color = "black", size = 4, shape = 1, stroke = 1.5) +
  geom_label_repel(data = min_p_value_row,
                   aes(label = gene_id),
                   box.padding = 0.8,
                   point.padding = 0.5,
                   min.segment.length = 0,
                   force = 20,
                   color = "black",
                   fill = "white",
                   alpha = 0.8,
                   size = 5,
                   fontface = "bold") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", alpha = 1) +
  scale_color_manual(values = c("Up" = "#d71345",
                                "Down" = "#9999FF",
                                "None" = "grey80")) +
  theme_classic() +
  labs(title = "panWSS1 vs. Others",
       x = expression(log[2]~"("*italic(FC)*")"),
       y = expression(-log[10]~"("*italic(p)*")"),
       color = "") +
  theme(plot.title = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 20, 10, 20))
#panWSS2
group_panWSS2=pancancer_cluster
group_panWSS2$cluster=ifelse(group_panWSS2$cluster=="panWSS2","H","L")
colnames(group_panWSS2)=c("Tag","group")
result=Get_limma(gene=TCGA_metabolism,
                 group=group_panWSS2,
                 alternative_limma="greater")

library(ggplot2)
library(ggrepel)
min_p_value_row = result[which.min(result$P.Value), ]
p2=ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value < 0.05,
                                ifelse(logFC > 0, "Up", "Down"),
                                "None")),
             alpha = 0.8, size = 3.8) +

  geom_point(data = min_p_value_row, color = "black", size = 4, shape = 1, stroke = 1.5) +
  geom_label_repel(data = min_p_value_row,
                   aes(label = gene_id),
                   box.padding = 0.8,
                   point.padding = 0.5,
                   min.segment.length = 0,
                   force = 20,
                   color = "black",
                   fill = "white",
                   alpha = 0.8,
                   size = 5,
                   fontface = "bold") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", alpha = 1) +
  scale_color_manual(values = c("Up" = "#f47920",
                                "Down" = "#9999FF",
                                "None" = "grey80")) +
  theme_classic() +
  labs(title = "panWSS2 vs. Others",
       x = expression(log[2]~"("*italic(FC)*")"),
       y = expression(-log[10]~"("*italic(p)*")"),
       color = "") +
  theme(plot.title = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 20, 10, 20))
#panWSS3
group_panWSS3=pancancer_cluster
group_panWSS3$cluster=ifelse(group_panWSS3$cluster=="panWSS3","H","L")
colnames(group_panWSS3)=c("Tag","group")
result=Get_limma(gene=TCGA_metabolism,
                 group=group_panWSS3,
                 alternative_limma="greater")

library(ggplot2)
library(ggrepel)
min_p_value_row = result[which.min(result$P.Value), ]
p3=ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value < 0.05,
                                ifelse(logFC > 0, "Up", "Down"),
                                "None")),
             alpha = 0.8, size = 3.8) +

  geom_point(data = min_p_value_row, color = "black", size = 4, shape = 1, stroke = 1.5) +
  geom_label_repel(data = min_p_value_row,
                   aes(label = gene_id),
                   box.padding = 0.8,
                   point.padding = 0.5,
                   min.segment.length = 0,
                   force = 20,
                   color = "black",
                   fill = "white",
                   alpha = 0.8,
                   size = 5,
                   fontface = "bold") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", alpha = 1) +
  scale_color_manual(values = c("Up" = "#145b7d",
                                "Down" = "#9999FF",
                                "None" = "grey80")) +
  theme_classic() +
  labs(title = "panWSS3 vs. Others",
       x = expression(log[2]~"("*italic(FC)*")"),
       y = expression(-log[10]~"("*italic(p)*")"),
       color = "") +
  theme(plot.title = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 20, 10, 20))
#panWSS4
group_panWSS4=pancancer_cluster
group_panWSS4$cluster=ifelse(group_panWSS4$cluster=="panWSS4","H","L")
colnames(group_panWSS4)=c("Tag","group")
result=Get_limma(gene=TCGA_metabolism,
                 group=group_panWSS4,
                 alternative_limma="greater")

library(ggplot2)
library(ggrepel)
min_p_value_row = result[which.min(result$P.Value), ]
p4=ggplot(result, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(P.Value < 0.05,
                                ifelse(logFC > 0, "Up", "Down"),
                                "None")),
             alpha = 0.8, size = 3.8) +

  geom_point(data = min_p_value_row, color = "black", size = 4, shape = 1, stroke = 1.5) +
  geom_label_repel(data = min_p_value_row,
                   aes(label = gene_id),
                   box.padding = 0.8,
                   point.padding = 0.5,
                   min.segment.length = 0,
                   force = 20,
                   color = "black",
                   fill = "white",
                   alpha = 0.8,
                   size = 5,
                   fontface = "bold") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", alpha = 1) +
  scale_color_manual(values = c("Up" = "#6950a1",
                                "Down" = "#9999FF",
                                "None" = "grey80")) +
  theme_classic() +
  labs(title = "panWSS4 vs. Others",
       x = expression(log[2]~"("*italic(FC)*")"),
       y = expression(-log[10]~"("*italic(p)*")"),
       color = "") +
  theme(plot.title = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "none",
        plot.margin = margin(10, 20, 10, 20))
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("panWSS1_meta.pdf",
       plot = p1,
       width = 5.5,
       height = 5,
       device = "pdf")
ggsave("panWSS2_meta.pdf",
       plot = p2,
       width = 5.5,
       height = 5,
       device = "pdf")
ggsave("panWSS3_meta.pdf",
       plot = p3,
       width = 5.5,
       height =5 ,
       device = "pdf")
ggsave("panWSS4_meta.pdf",
       plot = p4,
       width = 5.5,
       height = 5,
       device = "pdf")
################################################################################TCGA-Xenobiotics metabolism####
rm(list=ls())
gc()
library(reshape2)
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_cluster.Rdata")
group=pancancer_cluster
colnames(group)=c("Tag","group")
load("C:/Users/user/Desktop/input/TCGA_metabolism.Rdata")
tme_all=as.data.frame(t(TCGA_metabolism))
biomarker_group=merge(tme_all,group,by.x="row.names",by.y="Tag",all=T)
print(colnames(biomarker_group))
biomarker_group=biomarker_group[,-1]
biomarker_group=biomarker_group[,c("Xenobiotics metabolism","group")]
colnames(biomarker_group)=c("Xenobiotics metabolism","cluster")
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
p=ggplot(WNT, aes(x = cluster, y = `Xenobiotics metabolism`, fill = cluster)) +
  geom_violin(width = 0.5, scale = "area", trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) +
  xlab("") +
  ylab("Xenobiotics metabolism") +
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
ggsave("Xenobiotics_boxplot.pdf",
       plot = p,
       width = 5.5,
       height = 6,
       device = "pdf")
################################################################################TCGA-Xenobiotics metabolism(cor)####
rm(list=ls())
gc()
library(reshape2)
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_WNT_Score.Rdata")
activity_score=pancancer_WNT_Score$final_activity_score
load("C:/Users/user/Desktop/input/TCGA_metabolism.Rdata")
tme_all=as.data.frame(t(TCGA_metabolism))
biomarker_group=merge(tme_all,activity_score,by="row.names")
p <- ggscatter(biomarker_group, x = "Xenobiotics metabolism", y = "activity_score",
               add = "reg.line",
               conf.int = TRUE,
               add.params = list(color = "#d71345", fill = "#BBBDBE"),
               color = "#90d7ec", size = 1) +
  stat_cor(method = "pearson",
           label.x = 0,
           label.y = -3,
           size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 10, color = "black"),
    legend.position = "none"
  ) +
  ggtitle("METAFlux") +
  xlab("Xenobiotics metabolism") +
  ylab("Wnt/β-catenin pathway activity score")
p
setwd("C:/Users/user/Desktop/input_wnt/pictures")
ggsave("Xenobiotics_corrplot.pdf",
       plot = p,
       width = 5,
       height = 5,
       device = "pdf")
################################################################################Robustness test (randomly remove 30% of gene sets)####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")

# Load all datasets
vector_all <- c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
                "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
                "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
                "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
                "GSE73906_f","GSE123807","GSE89124")

for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q], ".Rdata"))
  load(paste0(vector_all[q], "_G.Rdata"))
}

# Get the unique study identifier (strip the suffix)
vector_marker <- unique(sub("_.*", "", vector_all))

# Set the number of repeats and the deletion proportion
n_repeats <- 30
delete_ratio <- 0.3

# Store all results
robustness_results <- list()

for (repeat_idx in 1:n_repeats) {
  set.seed(Sys.time() + sample(1:10000, 1))
  print(paste("Starting round ", repeat_idx, ""))

  # Randomly drop 30% of the studies
  n_to_delete <- round(length(vector_marker) * delete_ratio)
  random_studies <- sample(vector_marker, n_to_delete)

  # Build the vector of datasets to keep (drop every dataset belonging to the excluded study)
  vector_remaining <- vector_all[!grepl(paste0(random_studies, collapse = "|"), vector_all)]

  print(paste("Removed studies:", paste(random_studies, collapse = ", ")))
  print(paste("Datasets remaining:", length(vector_remaining)))

  # Read the group information
  group_HL <- read.delim("WNT_HL.txt", header = TRUE, sep = '\t', check.names = F)

  # Generate the de novo gene sets
  source("C:/Users/user/Desktop/functions_to_package/Get_denovo_genesets.R")

  up.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  dn.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  # Merge the de novo gene sets
  WNT_denovo <- cbind(up.denovo_genesets, dn.denovo_genesets)
  WNT_denovo <- rbind(colnames(WNT_denovo), WNT_denovo)
  colnames(WNT_denovo) <- paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")

  # Read the known gene sets
  WNT_knowledge <- read.delim("WNT_knowledge.txt", header = TRUE, sep = '\t', check.names = F)

  # Adjust the row count to match
  if (nrow(WNT_knowledge) >= nrow(WNT_denovo)) {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_knowledge) - nrow(WNT_denovo)),
                               ncol = ncol(WNT_denovo)))
    colnames(empty) <- colnames(WNT_denovo)
    WNT_denovo <- rbind(WNT_denovo, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  } else {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_denovo) - nrow(WNT_knowledge)),
                               ncol = ncol(WNT_knowledge)))
    colnames(empty) <- colnames(WNT_knowledge)
    WNT_knowledge <- rbind(WNT_knowledge, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  }

  # Add the enrichr gene sets
  files <- c("WNT_enrichr1.gmt", "WNT_enrichr2.gmt", "WNT_enrichr3.gmt", "WNT_enrichr4.gmt")
  combined_list <- list()

  for (file in files) {
    gmt_lines <- readLines(file)
    gene_sets <- strsplit(gmt_lines, "\t")
    names(gene_sets) <- sapply(gene_sets, function(x) paste0(x[1], "(", x[2], ")", "_enrichr"))
    for (i in 1:length(names(gene_sets))) {
      gene_sets[[i]][1] <- names(gene_sets)[i]
      gene_sets[[i]] <- gene_sets[[i]][-2]
    }
    combined_list <- append(combined_list, gene_sets)
  }

  combined_list_names <- as.data.frame(table(names(combined_list)))
  combined_list <- combined_list[combined_list_names$Var1]

  for (i in 1:length(combined_list_names$Var1)) {
    var_name <- paste("WPRGS", i + 152, "_knowledge", sep = "")
    names(combined_list)[i] <- var_name
  }

  # Process the gene set data
  all_genesets <- apply(all_genesets, 2, function(x) x[nchar(x) >= 1])
  all_genesets <- append(all_genesets, combined_list)

  max_length <- max(sapply(all_genesets, length))
  all_wnt_genesets <- lapply(all_genesets, function(x) {
    if (length(x) < max_length) {
      c(x, rep("", max_length - length(x)))
    } else {
      x
    }
  })

  all_wnt_genesets <- as.data.frame(all_wnt_genesets)
  save(all_wnt_genesets, file = "all_wnt_genesets.Rdata")

  # Load and process the gene sets
  load("C:/Users/user/Desktop/input/all_wnt_genesets.Rdata")
  all_genesets <- all_wnt_genesets[-1, ]
  UP.genesets <- all_genesets
  DN.genesets <- all_genesets

  # Clean the environment
  source("C:/Users/user/Desktop/functions_to_package/Precise_cleaning_system.R")

  UP.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = UP.genesets,
    inhibition_geneset = NA,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  DN.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = NA,
    inhibition_geneset = DN.genesets,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  UP.clean_genesets <- UP.clean_genesets_rank$update_activation_geneset
  DN.clean_genesets <- DN.clean_genesets_rank$update_inhibition_geneset

  # Combined gene set
  source("C:/Users/user/Desktop/functions_to_package/Joint_genesets.R")

  joint_genesets <- Joint_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    activation_geneset = UP.clean_genesets,
    inhibition_geneset = DN.clean_genesets,
    delete_GN = T,
    delete_time = "after",
    integration_method = c("Jaccard"),
    min_GN = 5,
    de_redundant_basis = "igraph_component",
    similarity = 0.5,
    export_file = F
  )

  activation_geneset <- joint_genesets$joint_activation_geneset
  inhibition_geneset <- joint_genesets$joint_inhibition_geneset

  # Pad rows so the lengths match
  if (nrow(activation_geneset) > nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(activation_geneset) - nrow(inhibition_geneset)),
                               ncol = ncol(inhibition_geneset)))
    colnames(empty) <- colnames(inhibition_geneset)
    inhibition_geneset <- rbind(inhibition_geneset, empty)
  } else if (nrow(activation_geneset) < nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(inhibition_geneset) - nrow(activation_geneset)),
                               ncol = ncol(activation_geneset)))
    colnames(empty) <- colnames(activation_geneset)
    activation_geneset <- rbind(activation_geneset, empty)
  }

  # Merge the activated and inhibited gene sets
  wpags_wpigs <- cbind(activation_geneset, inhibition_geneset)
  wpags_wpigs <- rbind(colnames(wpags_wpigs), wpags_wpigs)
  colnames(wpags_wpigs) <- c(paste0("wpags", 1:ncol(activation_geneset)),
                             paste0("wpigs", 1:ncol(inhibition_geneset)))

  wpags_wpigs <- apply(wpags_wpigs, 2, function(x) x[nchar(x) >= 1])

  # Write the GMT file
  write.gmt <- function(gs.list, file) {
    if (length(grep(".gmt$", file)) == 0)
      stop("write.gmt(): file should be .gmt")
    gs.names <- names(gs.list)
    gs.desc <- sapply(gs.list, function(x) x[1])
    gs.lines <- unlist(lapply(gs.list, function(x) {
      paste(x[-1], collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
  }

  write.gmt(wpags_wpigs, file = "C:/Users/user/Desktop/input/wpags_wpigs_r.gmt")

  # Load the test data and compute the activity score
  load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")

  crc_data_all <- UCSC.CRC_CMS2_CMS

  source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

  geneSets <- getGmt("wpags_wpigs_r.gmt")

  WNT_Score <- Calculate_bioactivity_scores(
    file_paths = "C:/Users/user/Desktop/input",
    expression_profile = crc_data_all,
    foundation = c("relative_ssGSEA"),
    activation_geneset = NA,
    inhibition_geneset = NA,
    geneSets_gmt = geneSets,
    min.sz = 1,
    max.sz = 10000,
    geneset_direction = c(rep("pos", ncol(activation_geneset)),
                          rep("neg", ncol(inhibition_geneset)))
  )

  # Save results
  CMS_result <- WNT_Score$final_activity_score
  CMS_result$rep <- paste0("repeat_", repeat_idx)

  robustness_results[[repeat_idx]] <- CMS_result
  names(robustness_results)[repeat_idx] <- paste0("repeat_", repeat_idx)

  print(paste("Finished round ", repeat_idx, ""))
}

# Save all results
save(robustness_results, file = "robustness_30percent_30repeats.Rdata")
print("Robustness test complete")
################################################################################Robustness test (randomly delete 30% of genes per dataset)####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")

# Load all datasets
vector_all <- c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
                "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
                "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
                "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
                "GSE73906_f","GSE123807","GSE89124")

# Load all data first
all_data <- list()
all_data_G <- list()

for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q], ".Rdata"))
  load(paste0(vector_all[q], "_G.Rdata"))
  all_data[[vector_all[q]]] <- get(vector_all[q])
  all_data_G[[vector_all[q]]] <- get(paste0(vector_all[q], "_G"))
}

# Set the number of repeats and the deletion proportion
n_repeats <- 30
delete_ratio <- 0.3

# Store all results
robustness_results <- list()

for (repeat_idx in 1:n_repeats) {
  set.seed(Sys.time() + sample(1:10000, 1))
  print(paste("Starting round ", repeat_idx, ""))

  # Randomly delete 30% of the genes in each dataset
  modified_data <- list()

  for (dataset_name in vector_all) {
    original_data <- all_data[[dataset_name]]
    original_genes <- rownames(original_data)

    # Randomly select 30% of the genes for deletion
    n_genes_to_delete <- round(length(original_genes) * delete_ratio)
    genes_to_delete <- sample(original_genes, n_genes_to_delete)

    # Build the modified data with the selected genes removed
    modified_data[[dataset_name]] <- original_data[!rownames(original_data) %in% genes_to_delete, ]

    print(paste("dataset", dataset_name, "removed", n_genes_to_delete, " genes"))
  }

  # Store the modified data in the global environment for later use
  for (dataset_name in vector_all) {
    assign(dataset_name, modified_data[[dataset_name]])
    assign(paste0(dataset_name, "_G"), all_data_G[[dataset_name]])
  }

  # Read the group information
  group_HL <- read.delim("WNT_HL.txt", header = TRUE, sep = '\t', check.names = F)

  # Generate the de novo gene sets
  source("C:/Users/user/Desktop/functions_to_package/Get_denovo_genesets.R")

  up.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  dn.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  # Merge the de novo gene sets
  WNT_denovo <- cbind(up.denovo_genesets, dn.denovo_genesets)
  WNT_denovo <- rbind(colnames(WNT_denovo), WNT_denovo)
  colnames(WNT_denovo) <- paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")

  # Read the known gene sets
  WNT_knowledge <- read.delim("WNT_knowledge.txt", header = TRUE, sep = '\t', check.names = F)

  # Adjust the row count to match
  if (nrow(WNT_knowledge) >= nrow(WNT_denovo)) {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_knowledge) - nrow(WNT_denovo)),
                               ncol = ncol(WNT_denovo)))
    colnames(empty) <- colnames(WNT_denovo)
    WNT_denovo <- rbind(WNT_denovo, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  } else {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_denovo) - nrow(WNT_knowledge)),
                               ncol = ncol(WNT_knowledge)))
    colnames(empty) <- colnames(WNT_knowledge)
    WNT_knowledge <- rbind(WNT_knowledge, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  }

  # Add the enrichr gene sets
  files <- c("WNT_enrichr1.gmt", "WNT_enrichr2.gmt", "WNT_enrichr3.gmt", "WNT_enrichr4.gmt")
  combined_list <- list()

  for (file in files) {
    gmt_lines <- readLines(file)
    gene_sets <- strsplit(gmt_lines, "\t")
    names(gene_sets) <- sapply(gene_sets, function(x) paste0(x[1], "(", x[2], ")", "_enrichr"))
    for (i in 1:length(names(gene_sets))) {
      gene_sets[[i]][1] <- names(gene_sets)[i]
      gene_sets[[i]] <- gene_sets[[i]][-2]
    }
    combined_list <- append(combined_list, gene_sets)
  }

  combined_list_names <- as.data.frame(table(names(combined_list)))
  combined_list <- combined_list[combined_list_names$Var1]

  for (i in 1:length(combined_list_names$Var1)) {
    var_name <- paste("WPRGS", i + 152, "_knowledge", sep = "")
    names(combined_list)[i] <- var_name
  }

  # Process the gene set data
  all_genesets <- apply(all_genesets, 2, function(x) x[nchar(x) >= 1])
  all_genesets <- append(all_genesets, combined_list)

  max_length <- max(sapply(all_genesets, length))
  all_wnt_genesets <- lapply(all_genesets, function(x) {
    if (length(x) < max_length) {
      c(x, rep("", max_length - length(x)))
    } else {
      x
    }
  })

  all_wnt_genesets <- as.data.frame(all_wnt_genesets)
  save(all_wnt_genesets, file = "all_wnt_genesets.Rdata")

  # Load and process the gene sets
  load("C:/Users/user/Desktop/input/all_wnt_genesets.Rdata")
  all_genesets <- all_wnt_genesets[-1, ]
  UP.genesets <- all_genesets
  DN.genesets <- all_genesets

  # Clean the environment
  source("C:/Users/user/Desktop/functions_to_package/Precise_cleaning_system.R")

  UP.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = UP.genesets,
    inhibition_geneset = NA,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  DN.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = NA,
    inhibition_geneset = DN.genesets,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  UP.clean_genesets <- UP.clean_genesets_rank$update_activation_geneset
  DN.clean_genesets <- DN.clean_genesets_rank$update_inhibition_geneset

  # Combined gene set
  source("C:/Users/user/Desktop/functions_to_package/Joint_genesets.R")

  joint_genesets <- Joint_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    activation_geneset = UP.clean_genesets,
    inhibition_geneset = DN.clean_genesets,
    delete_GN = T,
    delete_time = "after",
    integration_method = c("Jaccard"),
    min_GN = 5,
    de_redundant_basis = "igraph_component",
    similarity = 0.5,
    export_file = F
  )

  activation_geneset <- joint_genesets$joint_activation_geneset
  inhibition_geneset <- joint_genesets$joint_inhibition_geneset

  # Pad rows so the lengths match
  if (nrow(activation_geneset) > nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(activation_geneset) - nrow(inhibition_geneset)),
                               ncol = ncol(inhibition_geneset)))
    colnames(empty) <- colnames(inhibition_geneset)
    inhibition_geneset <- rbind(inhibition_geneset, empty)
  } else if (nrow(activation_geneset) < nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(inhibition_geneset) - nrow(activation_geneset)),
                               ncol = ncol(activation_geneset)))
    colnames(empty) <- colnames(activation_geneset)
    activation_geneset <- rbind(activation_geneset, empty)
  }

  # Merge the activated and inhibited gene sets
  wpags_wpigs <- cbind(activation_geneset, inhibition_geneset)
  wpags_wpigs <- rbind(colnames(wpags_wpigs), wpags_wpigs)
  colnames(wpags_wpigs) <- c(paste0("wpags", 1:ncol(activation_geneset)),
                             paste0("wpigs", 1:ncol(inhibition_geneset)))

  wpags_wpigs <- apply(wpags_wpigs, 2, function(x) x[nchar(x) >= 1])

  # Write the GMT file
  write.gmt <- function(gs.list, file) {
    if (length(grep(".gmt$", file)) == 0)
      stop("write.gmt(): file should be .gmt")
    gs.names <- names(gs.list)
    gs.desc <- sapply(gs.list, function(x) x[1])
    gs.lines <- unlist(lapply(gs.list, function(x) {
      paste(x[-1], collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
  }

  write.gmt(wpags_wpigs, file = "C:/Users/user/Desktop/input/wpags_wpigs_r.gmt")

  # Load the test data and compute the activity score
  load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")

  crc_data_all <- UCSC.CRC_CMS2_CMS

  source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

  geneSets <- getGmt("wpags_wpigs_r.gmt")

  WNT_Score <- Calculate_bioactivity_scores(
    file_paths = "C:/Users/user/Desktop/input",
    expression_profile = crc_data_all,
    foundation = c("relative_ssGSEA"),
    activation_geneset = NA,
    inhibition_geneset = NA,
    geneSets_gmt = geneSets,
    min.sz = 1,
    max.sz = 10000,
    geneset_direction = c(rep("pos", ncol(activation_geneset)),
                          rep("neg", ncol(inhibition_geneset)))
  )

  # Save results
  CMS_result <- WNT_Score$final_activity_score
  CMS_result$rep <- paste0("repeat_", repeat_idx)

  robustness_results[[repeat_idx]] <- CMS_result
  names(robustness_results)[repeat_idx] <- paste0("repeat_", repeat_idx)

  print(paste("Finished round ", repeat_idx, ""))

  # Clean the environment for the next replicate
  rm(list = setdiff(ls(), c("all_data", "all_data_G", "vector_all", "n_repeats",
                            "delete_ratio", "robustness_results", "repeat_idx")))
  gc()
}

# Save all results
save(robustness_results, file = "robustness_30percent_genes_30repeats.Rdata")
print("Gene-level robustness test complete")
################################################################################Robustness test (add Gaussian noise to each dataset)####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")

# Load all datasets
vector_all <- c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
                "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
                "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
                "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
                "GSE73906_f","GSE123807","GSE89124")

# Load all data first
all_data <- list()
all_data_G <- list()

for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q], ".Rdata"))
  load(paste0(vector_all[q], "_G.Rdata"))
  all_data[[vector_all[q]]] <- get(vector_all[q])
  all_data_G[[vector_all[q]]] <- get(paste0(vector_all[q], "_G"))
}

# Set the number of repeats and the noise parameters
n_repeats <- 30
noise_sd_ratio <- 0.1  # Noise SD as a ratio of the data SD

# Store all results
robustness_results <- list()

for (repeat_idx in 1:n_repeats) {
  set.seed(Sys.time() + sample(1:10000, 1))
  print(paste("Starting round ", repeat_idx, ""))

  # Add Gaussian noise to each dataset
  noisy_data <- list()

  for (dataset_name in vector_all) {
    original_data <- all_data[[dataset_name]]

    # Use the SD of the data as the noise baseline
    data_sd <- sd(as.matrix(original_data), na.rm = TRUE)
    noise_sd <- data_sd * noise_sd_ratio

    # Generate random Gaussian noise with the same dimensions as the data
    noise_matrix <- matrix(rnorm(nrow(original_data) * ncol(original_data),
                                 mean = 0, sd = noise_sd),
                           nrow = nrow(original_data), ncol = ncol(original_data))

    # Add the noise to the original data
    noisy_data[[dataset_name]] <- as.matrix(original_data) + noise_matrix

    # Keep the row and column names
    rownames(noisy_data[[dataset_name]]) <- rownames(original_data)
    colnames(noisy_data[[dataset_name]]) <- colnames(original_data)

    # Convert to a data frame
    noisy_data[[dataset_name]] <- as.data.frame(noisy_data[[dataset_name]])

    print(paste("dataset", dataset_name, "Added Gaussian noise with SD ", round(noise_sd, 4), ""))
  }

  # Store the noisy data in the global environment for later use
  for (dataset_name in vector_all) {
    assign(dataset_name, noisy_data[[dataset_name]])
    assign(paste0(dataset_name, "_G"), all_data_G[[dataset_name]])
  }

  # Read the group information
  group_HL <- read.delim("WNT_HL.txt", header = TRUE, sep = '\t', check.names = F)

  # Generate the de novo gene sets
  source("C:/Users/user/Desktop/functions_to_package/Get_denovo_genesets.R")

  up.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  dn.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  # Merge the de novo gene sets
  WNT_denovo <- cbind(up.denovo_genesets, dn.denovo_genesets)
  WNT_denovo <- rbind(colnames(WNT_denovo), WNT_denovo)
  colnames(WNT_denovo) <- paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")

  # Read the known gene sets
  WNT_knowledge <- read.delim("WNT_knowledge.txt", header = TRUE, sep = '\t', check.names = F)

  # Adjust the row count to match
  if (nrow(WNT_knowledge) >= nrow(WNT_denovo)) {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_knowledge) - nrow(WNT_denovo)),
                               ncol = ncol(WNT_denovo)))
    colnames(empty) <- colnames(WNT_denovo)
    WNT_denovo <- rbind(WNT_denovo, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  } else {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_denovo) - nrow(WNT_knowledge)),
                               ncol = ncol(WNT_knowledge)))
    colnames(empty) <- colnames(WNT_knowledge)
    WNT_knowledge <- rbind(WNT_knowledge, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  }

  # Add the enrichr gene sets
  files <- c("WNT_enrichr1.gmt", "WNT_enrichr2.gmt", "WNT_enrichr3.gmt", "WNT_enrichr4.gmt")
  combined_list <- list()

  for (file in files) {
    gmt_lines <- readLines(file)
    gene_sets <- strsplit(gmt_lines, "\t")
    names(gene_sets) <- sapply(gene_sets, function(x) paste0(x[1], "(", x[2], ")", "_enrichr"))
    for (i in 1:length(names(gene_sets))) {
      gene_sets[[i]][1] <- names(gene_sets)[i]
      gene_sets[[i]] <- gene_sets[[i]][-2]
    }
    combined_list <- append(combined_list, gene_sets)
  }

  combined_list_names <- as.data.frame(table(names(combined_list)))
  combined_list <- combined_list[combined_list_names$Var1]

  for (i in 1:length(combined_list_names$Var1)) {
    var_name <- paste("WPRGS", i + 152, "_knowledge", sep = "")
    names(combined_list)[i] <- var_name
  }

  # Process the gene set data
  all_genesets <- apply(all_genesets, 2, function(x) x[nchar(x) >= 1])
  all_genesets <- append(all_genesets, combined_list)

  max_length <- max(sapply(all_genesets, length))
  all_wnt_genesets <- lapply(all_genesets, function(x) {
    if (length(x) < max_length) {
      c(x, rep("", max_length - length(x)))
    } else {
      x
    }
  })

  all_wnt_genesets <- as.data.frame(all_wnt_genesets)
  save(all_wnt_genesets, file = "all_wnt_genesets.Rdata")

  # Load and process the gene sets
  load("C:/Users/user/Desktop/input/all_wnt_genesets.Rdata")
  all_genesets <- all_wnt_genesets[-1, ]
  UP.genesets <- all_genesets
  DN.genesets <- all_genesets

  # Clean the environment
  source("C:/Users/user/Desktop/functions_to_package/Precise_cleaning_system.R")

  UP.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = UP.genesets,
    inhibition_geneset = NA,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  DN.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_all,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = NA,
    inhibition_geneset = DN.genesets,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  UP.clean_genesets <- UP.clean_genesets_rank$update_activation_geneset
  DN.clean_genesets <- DN.clean_genesets_rank$update_inhibition_geneset

  # Combined gene set
  source("C:/Users/user/Desktop/functions_to_package/Joint_genesets.R")

  joint_genesets <- Joint_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    activation_geneset = UP.clean_genesets,
    inhibition_geneset = DN.clean_genesets,
    delete_GN = T,
    delete_time = "after",
    integration_method = c("Jaccard"),
    min_GN = 5,
    de_redundant_basis = "igraph_component",
    similarity = 0.5,
    export_file = F
  )

  activation_geneset <- joint_genesets$joint_activation_geneset
  inhibition_geneset <- joint_genesets$joint_inhibition_geneset

  # Pad rows so the lengths match
  if (nrow(activation_geneset) > nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(activation_geneset) - nrow(inhibition_geneset)),
                               ncol = ncol(inhibition_geneset)))
    colnames(empty) <- colnames(inhibition_geneset)
    inhibition_geneset <- rbind(inhibition_geneset, empty)
  } else if (nrow(activation_geneset) < nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(inhibition_geneset) - nrow(activation_geneset)),
                               ncol = ncol(activation_geneset)))
    colnames(empty) <- colnames(activation_geneset)
    activation_geneset <- rbind(activation_geneset, empty)
  }

  # Merge the activated and inhibited gene sets
  wpags_wpigs <- cbind(activation_geneset, inhibition_geneset)
  wpags_wpigs <- rbind(colnames(wpags_wpigs), wpags_wpigs)
  colnames(wpags_wpigs) <- c(paste0("wpags", 1:ncol(activation_geneset)),
                             paste0("wpigs", 1:ncol(inhibition_geneset)))

  wpags_wpigs <- apply(wpags_wpigs, 2, function(x) x[nchar(x) >= 1])

  # Write the GMT file
  write.gmt <- function(gs.list, file) {
    if (length(grep(".gmt$", file)) == 0)
      stop("write.gmt(): file should be .gmt")
    gs.names <- names(gs.list)
    gs.desc <- sapply(gs.list, function(x) x[1])
    gs.lines <- unlist(lapply(gs.list, function(x) {
      paste(x[-1], collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
  }

  write.gmt(wpags_wpigs, file = "C:/Users/user/Desktop/input/wpags_wpigs_r.gmt")

  # Load the test data and compute the activity score
  load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")

  crc_data_all <- UCSC.CRC_CMS2_CMS

  source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

  geneSets <- getGmt("wpags_wpigs_r.gmt")

  WNT_Score <- Calculate_bioactivity_scores(
    file_paths = "C:/Users/user/Desktop/input",
    expression_profile = crc_data_all,
    foundation = c("relative_ssGSEA"),
    activation_geneset = NA,
    inhibition_geneset = NA,
    geneSets_gmt = geneSets,
    min.sz = 1,
    max.sz = 10000,
    geneset_direction = c(rep("pos", ncol(activation_geneset)),
                          rep("neg", ncol(inhibition_geneset)))
  )

  # Save results
  CMS_result <- WNT_Score$final_activity_score
  CMS_result$rep <- paste0("repeat_", repeat_idx)
  CMS_result$noise_sd_ratio <- noise_sd_ratio

  robustness_results[[repeat_idx]] <- CMS_result
  names(robustness_results)[repeat_idx] <- paste0("repeat_", repeat_idx)

  print(paste("Finished round ", repeat_idx, ""))

  # Clean the environment for the next replicate
  rm(list = setdiff(ls(), c("all_data", "all_data_G", "vector_all", "n_repeats",
                            "noise_sd_ratio", "robustness_results", "repeat_idx")))
  gc()
}

# Save all results
save(robustness_results, file = "robustness_gaussian_noise_30repeats.Rdata")
print("Gaussian-noise robustness test complete")

################################################################################Robustness visualisation 1####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/robustness_30percent_30repeats.Rdata")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")
crc_data_all=UCSC.CRC_CMS2_CMS
geneSets <- getGmt("wpags_wpigs.gmt")
WNT_Score <- Calculate_bioactivity_scores(
  file_paths = "C:/Users/user/Desktop/input",
  expression_profile = crc_data_all,
  foundation = c("relative_ssGSEA"),
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                      rep("neg",sum(grepl("WPIGS", names(geneSets)))))
)
WNT_Score <- WNT_Score$final_activity_score
cor.data_all=data.frame()
for(i in 1:length(names(robustness_results))){
  score=robustness_results[[i]]
  cor.data=cor.test(score$activity_score,WNT_Score$activity_score)
  cor.data0=data.frame(rep=i,cor=cor.data$estimate)
  cor.data_all=rbind(cor.data_all,cor.data0)
}
# Load the required packages
library(ggplot2)
library(dplyr)
library(patchwork) # Helper to add statistics to the plot elegantly

# Compute the key summary statistics
n_samples <- nrow(cor.data_all)
mean_value <- mean(cor.data_all$cor)
median_value <- median(cor.data_all$cor)

# Create the labels for the statistics
stats_label <- paste(
  sprintf("n: %d", n_samples),
  sprintf("Mean: %.4f", mean_value),
  sprintf("Median: %.4f", median_value),
  sep = "\n" # Separate with line breaks
)

# Draw the density distribution with ggplot2
p <- ggplot(cor.data_all, aes(x = cor)) +
  # Draw the density curve with fill
  geom_density(fill = "#69b3a2",
               alpha = 0.6,
               color = "#69b3a2",
               linewidth = 0.8) +
  # Add a vertical line at the mean
  geom_vline(xintercept = mean_value,
             color = "red",
             linetype = "dashed",
             linewidth = 0.8) +
  # Add a vertical line at the median
  geom_vline(xintercept = median_value,
             color = "blue",
             linetype = "dashed",
             linewidth = 0.8) +
  # Set axis labels and title
  labs(
    title = "Randomly delete 30% of the training dataset",
    subtitle = "",
    x = "Pearson correlation coefficient",
    y = "Density"
  ) +
  # Use the classic theme
  theme_classic() +
  # Custom theme: all text at size 15
  theme(
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),  # Axis tick values
    legend.text = element_text(size = 15), # Legend text
    legend.title = element_text(size = 15) # Legend title
  ) +
  # Set a suitable x-range, adjusted dynamically to the data
  scale_x_continuous(
    limits = c(min(cor.data_all$cor) - 0.005,
               max(cor.data_all$cor) + 0.005),
    breaks = seq(0.93, 0.99, by = 0.01)
  ) +
  # Add the statistics label to the plot; black, size 15, not bold
  annotate("text",
           x = min(cor.data_all$cor) - 0.005,
           y = Inf,
           label = stats_label,
           hjust = 0,
           vjust = 1.2,
           size = 5, # 15 pt font corresponds to size = 5
           color = "black",
           fontface = "plain") # Explicitly not bold

# Display the figure
print(p)
################################################################################Robustness visualisation 2####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/robustness_30percent_genes_30repeats.Rdata")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")
crc_data_all=UCSC.CRC_CMS2_CMS
geneSets <- getGmt("wpags_wpigs.gmt")
WNT_Score <- Calculate_bioactivity_scores(
  file_paths = "C:/Users/user/Desktop/input",
  expression_profile = crc_data_all,
  foundation = c("relative_ssGSEA"),
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                      rep("neg",sum(grepl("WPIGS", names(geneSets)))))
)
WNT_Score <- WNT_Score$final_activity_score
cor.data_all=data.frame()
for(i in 1:length(names(robustness_results))){
  score=robustness_results[[i]]
  cor.data=cor.test(score$activity_score,WNT_Score$activity_score)
  cor.data0=data.frame(rep=i,cor=cor.data$estimate)
  cor.data_all=rbind(cor.data_all,cor.data0)
}
# Load the required packages
library(ggplot2)
library(dplyr)
library(patchwork) # Helper to add statistics to the plot elegantly

# Compute the key summary statistics
n_samples <- nrow(cor.data_all)
mean_value <- mean(cor.data_all$cor)
median_value <- median(cor.data_all$cor)

# Create the labels for the statistics
stats_label <- paste(
  sprintf("n: %d", n_samples),
  sprintf("Mean: %.4f", mean_value),
  sprintf("Median: %.4f", median_value),
  sep = "\n" # Separate with line breaks
)

# Draw the density distribution with ggplot2
p <- ggplot(cor.data_all, aes(x = cor)) +
  # Draw the density curve with fill
  geom_density(fill = "#69b3a2",
               alpha = 0.6,
               color = "#69b3a2",
               linewidth = 0.8) +
  # Add a vertical line at the mean
  geom_vline(xintercept = mean_value,
             color = "red",
             linetype = "dashed",
             linewidth = 0.8) +
  # Add a vertical line at the median
  geom_vline(xintercept = median_value,
             color = "blue",
             linetype = "dashed",
             linewidth = 0.8) +
  # Set axis labels and title
  labs(
    title = "Randomly delete 30% of the genes",
    subtitle = "",
    x = "Pearson correlation coefficient",
    y = "Density"
  ) +
  # Use the classic theme
  theme_classic() +
  # Custom theme: all text at size 15
  theme(
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),  # Axis tick values
    legend.text = element_text(size = 15), # Legend text
    legend.title = element_text(size = 15) # Legend title
  ) +
  # Set a suitable x-range, adjusted dynamically to the data
  scale_x_continuous(
    limits = c(min(cor.data_all$cor) - 0.005,
               max(cor.data_all$cor) + 0.005),
    breaks = seq(0.60, 0.92, by = 0.05)
  ) +
  # Add the statistics label to the plot; black, size 15, not bold
  annotate("text",
           x = min(cor.data_all$cor) - 0.005,
           y = Inf,
           label = stats_label,
           hjust = 0,
           vjust = 1.2,
           size = 5, # 15 pt font corresponds to size = 5
           color = "black",
           fontface = "plain") # Explicitly not bold

# Display the figure
print(p)
################################################################################Robustness visualisation 3####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/robustness_gaussian_noise_30repeats.Rdata")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")
crc_data_all=UCSC.CRC_CMS2_CMS
geneSets <- getGmt("wpags_wpigs.gmt")
WNT_Score <- Calculate_bioactivity_scores(
  file_paths = "C:/Users/user/Desktop/input",
  expression_profile = crc_data_all,
  foundation = c("relative_ssGSEA"),
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                      rep("neg",sum(grepl("WPIGS", names(geneSets)))))
)
WNT_Score <- WNT_Score$final_activity_score
cor.data_all=data.frame()
for(i in 1:length(names(robustness_results))){
  score=robustness_results[[i]]
  cor.data=cor.test(score$activity_score,WNT_Score$activity_score)
  cor.data0=data.frame(rep=i,cor=cor.data$estimate)
  cor.data_all=rbind(cor.data_all,cor.data0)
}
# Load the required packages
library(ggplot2)
library(dplyr)
library(patchwork) # Helper to add statistics to the plot elegantly

# Compute the key summary statistics
n_samples <- nrow(cor.data_all)
mean_value <- mean(cor.data_all$cor)
median_value <- median(cor.data_all$cor)

# Create the labels for the statistics
stats_label <- paste(
  sprintf("n: %d", n_samples),
  sprintf("Mean: %.4f", mean_value),
  sprintf("Median: %.4f", median_value),
  sep = "\n" # Separate with line breaks
)

# Draw the density distribution with ggplot2
p <- ggplot(cor.data_all, aes(x = cor)) +
  # Draw the density curve with fill
  geom_density(fill = "#69b3a2",
               alpha = 0.6,
               color = "#69b3a2",
               linewidth = 0.8) +
  # Add a vertical line at the mean
  geom_vline(xintercept = mean_value,
             color = "red",
             linetype = "dashed",
             linewidth = 0.8) +
  # Add a vertical line at the median
  geom_vline(xintercept = median_value,
             color = "blue",
             linetype = "dashed",
             linewidth = 0.8) +
  # Set axis labels and title
  labs(
    title = "Introduce random Gaussian noise",
    subtitle = "",
    x = "Pearson correlation coefficient",
    y = "Density"
  ) +
  # Use the classic theme
  theme_classic() +
  # Custom theme: all text at size 15
  theme(
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),  # Axis tick values
    legend.text = element_text(size = 15), # Legend text
    legend.title = element_text(size = 15) # Legend title
  ) +
  # Set a suitable x-range, adjusted dynamically to the data
  scale_x_continuous(
    limits = c(min(cor.data_all$cor) - 0.005,
               max(cor.data_all$cor) + 0.005),
    breaks = seq(0.97, 0.99, by = 0.005)
  ) +
  # Add the statistics label to the plot; black, size 15, not bold
  annotate("text",
           x = min(cor.data_all$cor) - 0.005,
           y = Inf,
           label = stats_label,
           hjust = 0,
           vjust = 1.2,
           size = 5, # 15 pt font corresponds to size = 5
           color = "black",
           fontface = "plain") # Explicitly not bold

# Display the figure
print(p)
################################################################################Revised sensitivity calculation####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")

# Load all datasets
vector_all <- c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
                "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
                "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
                "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
                "GSE73906_f","GSE123807","GSE89124")

for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q], ".Rdata"))
  load(paste0(vector_all[q], "_G.Rdata"))
}

# Use the full dataset
vector_remaining <- vector_all

# Set the similarity range (0 to 1, step 0.1)
similarity_values <- seq(0, 1, by = 0.1)

# Store all results and the gene-set statistics
similarity_results <- list()
geneset_count_summary <- data.frame(
  similarity = numeric(),
  n_activation = integer(),
  n_inhibition = integer(),
  total_genesets = integer()
)

for (sim in similarity_values) {
  print(paste("========================================"))
  print(paste("Computing similarity =", sim))

  # Read the group information
  group_HL <- read.delim("WNT_HL.txt", header = TRUE, sep = '\t', check.names = F)

  # Generate the de novo gene sets
  source("C:/Users/user/Desktop/functions_to_package/Get_denovo_genesets.R")

  up.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  dn.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  # Merge the de novo gene sets
  WNT_denovo <- cbind(up.denovo_genesets, dn.denovo_genesets)
  WNT_denovo <- rbind(colnames(WNT_denovo), WNT_denovo)
  colnames(WNT_denovo) <- paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")

  # Read the known gene sets
  WNT_knowledge <- read.delim("WNT_knowledge.txt", header = TRUE, sep = '\t', check.names = F)

  # Adjust the row count to match
  if (nrow(WNT_knowledge) >= nrow(WNT_denovo)) {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_knowledge) - nrow(WNT_denovo)),
                               ncol = ncol(WNT_denovo)))
    colnames(empty) <- colnames(WNT_denovo)
    WNT_denovo <- rbind(WNT_denovo, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  } else {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_denovo) - nrow(WNT_knowledge)),
                               ncol = ncol(WNT_knowledge)))
    colnames(empty) <- colnames(WNT_knowledge)
    WNT_knowledge <- rbind(WNT_knowledge, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  }

  # Add the enrichr gene sets
  files <- c("WNT_enrichr1.gmt", "WNT_enrichr2.gmt", "WNT_enrichr3.gmt", "WNT_enrichr4.gmt")
  combined_list <- list()

  for (file in files) {
    gmt_lines <- readLines(file)
    gene_sets <- strsplit(gmt_lines, "\t")
    names(gene_sets) <- sapply(gene_sets, function(x) paste0(x[1], "(", x[2], ")", "_enrichr"))
    for (i in 1:length(names(gene_sets))) {
      gene_sets[[i]][1] <- names(gene_sets)[i]
      gene_sets[[i]] <- gene_sets[[i]][-2]
    }
    combined_list <- append(combined_list, gene_sets)
  }

  combined_list_names <- as.data.frame(table(names(combined_list)))
  combined_list <- combined_list[combined_list_names$Var1]

  for (i in 1:length(combined_list_names$Var1)) {
    var_name <- paste("WPRGS", i + 152, "_knowledge", sep = "")
    names(combined_list)[i] <- var_name
  }

  # Process the gene set data
  all_genesets <- apply(all_genesets, 2, function(x) x[nchar(x) >= 1])
  all_genesets <- append(all_genesets, combined_list)

  max_length <- max(sapply(all_genesets, length))
  all_wnt_genesets <- lapply(all_genesets, function(x) {
    if (length(x) < max_length) {
      c(x, rep("", max_length - length(x)))
    } else {
      x
    }
  })

  all_wnt_genesets <- as.data.frame(all_wnt_genesets)
  save(all_wnt_genesets, file = "all_wnt_genesets.Rdata")

  # Load and process the gene sets
  load("C:/Users/user/Desktop/input/all_wnt_genesets.Rdata")
  all_genesets <- all_wnt_genesets[-1, ]
  UP.genesets <- all_genesets
  DN.genesets <- all_genesets

  # Clean the environment
  source("C:/Users/user/Desktop/functions_to_package/Precise_cleaning_system.R")

  UP.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = UP.genesets,
    inhibition_geneset = NA,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  DN.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = NA,
    inhibition_geneset = DN.genesets,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  UP.clean_genesets <- UP.clean_genesets_rank$update_activation_geneset
  DN.clean_genesets <- DN.clean_genesets_rank$update_inhibition_geneset

  # Combined gene set
  source("C:/Users/user/Desktop/functions_to_package/Joint_genesets.R")

  joint_genesets <- Joint_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    activation_geneset = UP.clean_genesets,
    inhibition_geneset = DN.clean_genesets,
    delete_GN = T,
    delete_time = "after",
    integration_method = c("Jaccard"),
    min_GN = 5,
    de_redundant_basis = "igraph_component",
    similarity = sim,  # Use the current similarity value
    export_file = F
  )

  activation_geneset <- joint_genesets$joint_activation_geneset
  inhibition_geneset <- joint_genesets$joint_inhibition_geneset

  # Count the gene sets generated so far
  n_act <- ncol(activation_geneset)
  n_inh <- ncol(inhibition_geneset)
  n_total <- n_act + n_inh

  print(paste("[Gene sets] activated:", n_act, "| inhibited:", n_inh, "| total:", n_total))

  # Record in the summary table
  geneset_count_summary <- rbind(geneset_count_summary, data.frame(
    similarity = sim,
    n_activation = n_act,
    n_inhibition = n_inh,
    total_genesets = n_total
  ))

  # Pad rows so the lengths match
  if (nrow(activation_geneset) > nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(activation_geneset) - nrow(inhibition_geneset)),
                               ncol = ncol(inhibition_geneset)))
    colnames(empty) <- colnames(inhibition_geneset)
    inhibition_geneset <- rbind(inhibition_geneset, empty)
  } else if (nrow(activation_geneset) < nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(inhibition_geneset) - nrow(activation_geneset)),
                               ncol = ncol(activation_geneset)))
    colnames(empty) <- colnames(activation_geneset)
    activation_geneset <- rbind(activation_geneset, empty)
  }

  # Merge the activated and inhibited gene sets
  wpags_wpigs <- cbind(activation_geneset, inhibition_geneset)
  wpags_wpigs <- rbind(colnames(wpags_wpigs), wpags_wpigs)
  colnames(wpags_wpigs) <- c(paste0("wpags", 1:ncol(activation_geneset)),
                             paste0("wpigs", 1:ncol(inhibition_geneset)))

  wpags_wpigs <- apply(wpags_wpigs, 2, function(x) x[nchar(x) >= 1])

  # Write the GMT file
  write.gmt <- function(gs.list, file) {
    if (length(grep(".gmt$", file)) == 0)
      stop("write.gmt(): file should be .gmt")
    gs.names <- names(gs.list)
    gs.desc <- sapply(gs.list, function(x) x[1])
    gs.lines <- unlist(lapply(gs.list, function(x) {
      paste(x[-1], collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
  }

  write.gmt(wpags_wpigs, file = "C:/Users/user/Desktop/input/wpags_wpigs_r.gmt")

  # Load the test data and compute the activity score
  load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")

  crc_data_all <- UCSC.CRC_CMS2_CMS

  source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

  geneSets <- getGmt("wpags_wpigs_r.gmt")

  WNT_Score <- Calculate_bioactivity_scores(
    file_paths = "C:/Users/user/Desktop/input",
    expression_profile = crc_data_all,
    foundation = c("relative_ssGSEA"),
    activation_geneset = NA,
    inhibition_geneset = NA,
    geneSets_gmt = geneSets,
    min.sz = 1,
    max.sz = 10000,
    geneset_direction = c(rep("pos", ncol(activation_geneset)),
                          rep("neg", ncol(inhibition_geneset)))
  )

  # Save the results for the current similarity (scores plus the gene-set statistics for that setting)
  CMS_result <- WNT_Score$final_activity_score
  CMS_result$similarity <- sim

  similarity_results[[paste0("sim_", sim)]] <- list(
    activity_score = CMS_result,
    n_activation = n_act,
    n_inhibition = n_inh,
    total_genesets = n_total
  )

  print(paste("Done similarity =", sim))
}

# Save all results and the summary table of gene-set counts
save(similarity_results, geneset_count_summary, file = "similarity_sensitivity_results.Rdata")

print("========================================")
print("Similarity sensitivity analysis and gene-set tally complete")
print("Gene-set count summary:")
print(geneset_count_summary)


library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Load the results saved earlier
load("similarity_sensitivity_results.Rdata")

# 2. Data conversion (long format)
df_long <- geneset_count_summary %>%
  select(similarity, n_activation, n_inhibition) %>%
  pivot_longer(cols = c(n_activation, n_inhibition),
               names_to = "Geneset_Type",
               values_to = "Count")

# Convert to plain English labels
df_long$Geneset_Type <- factor(df_long$Geneset_Type,
                               levels = c("n_activation", "n_inhibition"),
                               labels = c("Activation", "Inhibition"))

# 3. Draw the stacked bar chart
p <- ggplot(df_long, aes(x = as.factor(similarity), y = Count, fill = Geneset_Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6, alpha = 0.85) +
  # Add a total label above each bar (black, 14 pt, not bold)
  geom_text(data = geneset_count_summary,
            aes(x = as.factor(similarity), y = total_genesets, label = total_genesets),
            inherit.aes = FALSE, vjust = -0.5, size = 5, color = "black", fontface = "plain") +
  # Custom colour scheme
  scale_fill_manual(values = c("Activation" = "#E64B35FF",
                               "Inhibition" = "#4DBBD5FF")) +
  # Remove the main and sub titles
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Similarity",
    y = "Number of Gene Sets",
    fill = NULL
  ) +
  # Theme: 14 pt, black, not bold
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 14, color = "black", face = "plain"),
    axis.text = element_text(size = 14, color = "black", face = "plain"),
    axis.title = element_text(size = 14, color = "black", face = "plain"),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  ) +
  # Leave headroom so the total-count labels at the top are not clipped
  ylim(0, max(geneset_count_summary$total_genesets) * 1.15)
p

library(ggplot2)
library(reshape2)

# 1. Load results
load("similarity_sensitivity_results.Rdata")

# 2. Extract final scores into a matrix
sims <- seq(0, 1, by = 0.1)
score_list <- list()

for (sim in sims) {
  key <- paste0("sim_", sim)
  df <- similarity_results[[key]]$activity_score
  num_col <- which(sapply(df, is.numeric))[1]
  score_list[[paste0("sim_", sim)]] <- df[, num_col]
}

score_matrix <- do.call(cbind, score_list)
colnames(score_matrix) <- as.character(sims)

# 3. Compute correlation matrix and p-value matrix
cor_matrix <- cor(score_matrix, use = "complete.obs", method = "pearson")

p_matrix <- matrix(NA, nrow = ncol(score_matrix), ncol = ncol(score_matrix))
rownames(p_matrix) <- colnames(p_matrix) <- colnames(score_matrix)

for(i in 1:ncol(score_matrix)) {
  for(j in 1:ncol(score_matrix)) {
    test <- cor.test(score_matrix[,i], score_matrix[,j], method = "pearson")
    p_matrix[i, j] <- test$p.value
  }
}

# 4. Melt matrices for ggplot2
cor_melt <- melt(cor_matrix)
colnames(cor_melt) <- c("Var1", "Var2", "Correlation")

p_melt <- melt(p_matrix)
colnames(p_melt) <- c("Var1", "Var2", "Pvalue")

plot_data <- merge(cor_melt, p_melt, by = c("Var1", "Var2"))

# Keep the axes strictly aligned and in the same order
sim_labels <- as.character(sims)
plot_data$Var1 <- factor(plot_data$Var1, levels = sim_labels)
plot_data$Var2 <- factor(plot_data$Var2, levels = sim_labels)

# Add significance stars
plot_data$Sig <- ifelse(plot_data$Pvalue < 0.001, "***",
                        ifelse(plot_data$Pvalue < 0.01, "**",
                               ifelse(plot_data$Pvalue < 0.05, "*", "")))

# Combine correlation coefficient R and significance stars into label
plot_data$Label <- ifelse(plot_data$Sig == "",
                          sprintf("%.2f", plot_data$Correlation),
                          paste0(sprintf("%.2f", plot_data$Correlation), "\n", plot_data$Sig))

# 5. Plot Modern Bubble Correlation Matrix with lighter colors
p <- ggplot(plot_data, aes(x = Var1, y = Var2)) +
  # Light grid at the bottom to aid alignment
  geom_point(color = "gray95", size = 6) +
  # Show the correlation with dot size and colour
  geom_point(aes(size = abs(Correlation), color = Correlation), alpha = 0.85) +
  scale_size(range = c(2, 10), guide = "none") +
  # Softer, lower-saturation blue and red palette
  scale_color_gradient2(low = "#A9CCE3", mid = "white", high = "#F5B7B1", midpoint = 0,
                        limit = c(-1, 1), name = "Pearson R") +
  # Add the values and significance asterisks
  geom_text(aes(label = Label), size = 3.5, color = "black", fontface = "plain") +
  labs(title = NULL, subtitle = NULL, x = "Similarity", y = "Similarity") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 14, color = "black", face = "plain"),
    axis.text = element_text(size = 14, color = "black", face = "plain"),
    axis.title = element_text(size = 14, color = "black", face = "plain"),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.title = element_text(size = 14, color = "black", face = "plain"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  ) +
  coord_fixed()

p
################################################################################Comparison of redundancy-reduction methods####
rm(list=ls())
gc()
library(GSEABase)
setwd("C:/Users/user/Desktop/input")

# Load all datasets
vector_all <- c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
                "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
                "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
                "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
                "GSE73906_f","GSE123807","GSE89124")

for (q in 1:length(vector_all)) {
  load(paste0(vector_all[q], ".Rdata"))
  load(paste0(vector_all[q], "_G.Rdata"))
}

# Use the full dataset
vector_remaining <- vector_all

# Fix the similarity parameter at 0.5
fixed_similarity <- 0.5

# Set the 7 redundancy-reduction methods
de_redundant_methods <- c(
  "igraph_component", "agglomeration_ward.D",
  "agglomeration_ward.D2", "agglomeration_single", "agglomeration_complete",
  "agglomeration_average", "agglomeration_mcquitty"
)

# Store the results and gene-set statistics of every method
de_redundant_results <- list()
geneset_count_summary_methods <- data.frame(
  method = character(),
  n_activation = integer(),
  n_inhibition = integer(),
  total_genesets = integer(),
  stringsAsFactors = FALSE
)

# Loop over the 7 redundancy-reduction methods
for (method_idx in 1:length(de_redundant_methods)) {
  current_method <- de_redundant_methods[method_idx]
  print(paste("################################################"))
  print(paste("Running redundancy method [", method_idx, "/ 7]:", current_method, "| fixed similarity =", fixed_similarity))
  print(paste("################################################"))

  # Read the group information
  group_HL <- read.delim("WNT_HL.txt", header = TRUE, sep = '\t', check.names = F)

  # Generate the de novo gene sets
  source("C:/Users/user/Desktop/functions_to_package/Get_denovo_genesets.R")

  up.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  dn.denovo_genesets <- Get_denovo_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    threshold_or_rank = c("rank"),
    top_genes = 150,
    gene_Pfilter = NA,
    gene_FCfilter = NA,
    export_file = F
  )

  # Merge the de novo gene sets
  WNT_denovo <- cbind(up.denovo_genesets, dn.denovo_genesets)
  WNT_denovo <- rbind(colnames(WNT_denovo), WNT_denovo)
  colnames(WNT_denovo) <- paste0("WPRGS", 1:ncol(WNT_denovo), "_denovo")

  # Read the known gene sets
  WNT_knowledge <- read.delim("WNT_knowledge.txt", header = TRUE, sep = '\t', check.names = F)

  # Adjust the row count to match
  if (nrow(WNT_knowledge) >= nrow(WNT_denovo)) {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_knowledge) - nrow(WNT_denovo)),
                               ncol = ncol(WNT_denovo)))
    colnames(empty) <- colnames(WNT_denovo)
    WNT_denovo <- rbind(WNT_denovo, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  } else {
    empty <- data.frame(matrix("", nrow = (nrow(WNT_denovo) - nrow(WNT_knowledge)),
                               ncol = ncol(WNT_knowledge)))
    colnames(empty) <- colnames(WNT_knowledge)
    WNT_knowledge <- rbind(WNT_knowledge, empty)
    all_genesets <- cbind(WNT_denovo, WNT_knowledge)
  }

  # Add the enrichr gene sets
  files <- c("WNT_enrichr1.gmt", "WNT_enrichr2.gmt", "WNT_enrichr3.gmt", "WNT_enrichr4.gmt")
  combined_list <- list()

  for (file in files) {
    gmt_lines <- readLines(file)
    gene_sets <- strsplit(gmt_lines, "\t")
    names(gene_sets) <- sapply(gene_sets, function(x) paste0(x[1], "(", x[2], ")", "_enrichr"))
    for (i in 1:length(names(gene_sets))) {
      gene_sets[[i]][1] <- names(gene_sets)[i]
      gene_sets[[i]] <- gene_sets[[i]][-2]
    }
    combined_list <- append(combined_list, gene_sets)
  }

  combined_list_names <- as.data.frame(table(names(combined_list)))
  combined_list <- combined_list[combined_list_names$Var1]

  for (i in 1:length(combined_list_names$Var1)) {
    var_name <- paste("WPRGS", i + 152, "_knowledge", sep = "")
    names(combined_list)[i] <- var_name
  }

  # Process the gene set data
  all_genesets <- apply(all_genesets, 2, function(x) x[nchar(x) >= 1])
  all_genesets <- append(all_genesets, combined_list)

  max_length <- max(sapply(all_genesets, length))
  all_wnt_genesets <- lapply(all_genesets, function(x) {
    if (length(x) < max_length) {
      c(x, rep("", max_length - length(x)))
    } else {
      x
    }
  })

  all_wnt_genesets <- as.data.frame(all_wnt_genesets)
  save(all_wnt_genesets, file = "all_wnt_genesets.Rdata")

  # Load and process the gene sets
  load("C:/Users/user/Desktop/input/all_wnt_genesets.Rdata")
  all_genesets <- all_wnt_genesets[-1, ]
  UP.genesets <- all_genesets
  DN.genesets <- all_genesets

  # Clean the environment
  source("C:/Users/user/Desktop/functions_to_package/Precise_cleaning_system.R")

  UP.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("greater"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = UP.genesets,
    inhibition_geneset = NA,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  DN.clean_genesets_rank <- Precise_cleaning_system(
    file_paths = "C:/Users/user/Desktop/input",
    expression_accession_vector = vector_remaining,
    group_HL = group_HL,
    gene_difference_method = c("limma"),
    alternative = c("less"),
    p_combine_method = c("geometric_mean"),
    using_FC = F,
    na_ratio = 0.3,
    using_KNN = F,
    statistics = c("sumz"),
    purpose = c("cleaned"),
    threshold_type = c("rank"),
    rank_threshold = 150,
    quantile_threshold = NA,
    activation_geneset = NA,
    inhibition_geneset = DN.genesets,
    geneSets_gmt = NA,
    min.sz = 1,
    max.sz = 10000,
    export_file = F
  )

  UP.clean_genesets <- UP.clean_genesets_rank$update_activation_geneset
  DN.clean_genesets <- DN.clean_genesets_rank$update_inhibition_geneset

  # Combined gene set (current de_redundant_basis with similarity fixed at 0.5)
  source("C:/Users/user/Desktop/functions_to_package/Joint_genesets.R")

  joint_genesets <- Joint_genesets(
    file_paths = "C:/Users/user/Desktop/input",
    activation_geneset = UP.clean_genesets,
    inhibition_geneset = DN.clean_genesets,
    delete_GN = T,
    delete_time = "after",
    integration_method = c("Jaccard"),
    min_GN = 5,
    de_redundant_basis = current_method,
    similarity = fixed_similarity,
    export_file = F
  )

  activation_geneset <- joint_genesets$joint_activation_geneset
  inhibition_geneset <- joint_genesets$joint_inhibition_geneset

  # Count the gene sets generated so far
  n_act <- ncol(activation_geneset)
  n_inh <- ncol(inhibition_geneset)
  n_total <- n_act + n_inh

  print(paste("[Gene sets] activated:", n_act, "| inhibited:", n_inh, "| total:", n_total))

  # Record in the summary table
  geneset_count_summary_methods <- rbind(geneset_count_summary_methods, data.frame(
    method = current_method,
    n_activation = n_act,
    n_inhibition = n_inh,
    total_genesets = n_total
  ))

  # Pad rows so the lengths match
  if (nrow(activation_geneset) > nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(activation_geneset) - nrow(inhibition_geneset)),
                               ncol = ncol(inhibition_geneset)))
    colnames(empty) <- colnames(inhibition_geneset)
    inhibition_geneset <- rbind(inhibition_geneset, empty)
  } else if (nrow(activation_geneset) < nrow(inhibition_geneset)) {
    empty <- data.frame(matrix("", nrow = (nrow(inhibition_geneset) - nrow(activation_geneset)),
                               ncol = ncol(activation_geneset)))
    colnames(empty) <- colnames(activation_geneset)
    activation_geneset <- rbind(activation_geneset, empty)
  }

  # Merge the activated and inhibited gene sets
  wpags_wpigs <- cbind(activation_geneset, inhibition_geneset)
  wpags_wpigs <- rbind(colnames(wpags_wpigs), wpags_wpigs)
  colnames(wpags_wpigs) <- c(paste0("wpags", 1:ncol(activation_geneset)),
                             paste0("wpigs", 1:ncol(inhibition_geneset)))

  wpags_wpigs <- apply(wpags_wpigs, 2, function(x) x[nchar(x) >= 1])

  # Write the GMT file
  write.gmt <- function(gs.list, file) {
    if (length(grep(".gmt$", file)) == 0)
      stop("write.gmt(): file should be .gmt")
    gs.names <- names(gs.list)
    gs.desc <- sapply(gs.list, function(x) x[1])
    gs.lines <- unlist(lapply(gs.list, function(x) {
      paste(x[-1], collapse = "\t")
    }))
    gs.lines <- paste(gs.names, gs.desc, gs.lines, sep = "\t")
    writeLines(gs.lines, con = file)
  }

  write.gmt(wpags_wpigs, file = "C:/Users/user/Desktop/input/wpags_wpigs_r.gmt")

  # Load the test data and compute the activity score
  load("C:/Users/user/Desktop/input/exp_data/UCSC.CRC_CMS2_CMS.Rdata")

  crc_data_all <- UCSC.CRC_CMS2_CMS

  source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")

  geneSets <- getGmt("wpags_wpigs_r.gmt")

  WNT_Score <- Calculate_bioactivity_scores(
    file_paths = "C:/Users/user/Desktop/input",
    expression_profile = crc_data_all,
    foundation = c("relative_ssGSEA"),
    activation_geneset = NA,
    inhibition_geneset = NA,
    geneSets_gmt = geneSets,
    min.sz = 1,
    max.sz = 10000,
    geneset_direction = c(rep("pos", ncol(activation_geneset)),
                          rep("neg", ncol(inhibition_geneset)))
  )

  CMS_result <- WNT_Score$final_activity_score
  CMS_result$method <- current_method

  de_redundant_results[[current_method]] <- list(
    activity_score = CMS_result,
    n_activation = n_act,
    n_inhibition = n_inh,
    total_genesets = n_total
  )

  print(paste("Done method:", current_method))
}

# Save the final output of every method together with the gene-set summary
save(de_redundant_results, geneset_count_summary_methods, file = "de_redundant_fixed_sim0.5_results.Rdata")

print("========================================")
print("All 7 methods finished at similarity = 0.5")
print("Gene-set summary:")
print(geneset_count_summary_methods)
library(ggplot2)
library(dplyr)
library(tidyr)


# If the script just ran, geneset_count_summary_methods already exists;
# Uncomment the next line when reloading the data:
# load("de_redundant_fixed_sim0.5_results.Rdata")

# 1. Reshape to long format, drop agglomeration_single, and strip the "agglomeration_" prefix from the method names
df_summary_clean <- geneset_count_summary_methods %>%
  filter(method != "agglomeration_single") %>% # <--- Remove agglomeration_single here
  mutate(method = gsub("agglomeration_", "", method))

df_long <- df_summary_clean %>%
  pivot_longer(cols = c(n_activation, n_inhibition),
               names_to = "type",
               values_to = "count") %>%
  mutate(type = case_when(
    type == "n_activation" ~ "Activation",
    type == "n_inhibition" ~ "Inhibition"
  ))

# 2. Count the total per method, used for sorting
df_total <- df_summary_clean %>%
  mutate(method = factor(method, levels = method[order(total_genesets)]))

df_long$method <- factor(df_long$method, levels = levels(df_total$method))

# 3. Draw the stacked bar chart: no title, no subtitle, plain black, 14 pt
p_stacked <- ggplot(df_long, aes(x = method, y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  # Segment labels inside the bars (14 pt, not bold)
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            size = 4.9, color = "white", fontface = "plain") +
  # Total label at the top of each bar (black, not bold, 14 pt)
  geom_text(data = df_total,
            aes(x = method, y = total_genesets, label = paste0("Total: ", total_genesets), fill = NULL),
            vjust = -0.6, size = 4.9, color = "black", fontface = "plain") +
  scale_fill_manual(values = c("Activation" = "#E64B35FF", "Inhibition" = "#4DBBD5FF")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Method",
    y = "Number of Gene Sets",
    fill = "Gene Set Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(color = "black", size = 14, face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 14, face = "plain"),
    axis.text.y = element_text(color = "black", size = 14, face = "plain"),
    axis.title = element_text(color = "black", size = 14, face = "plain"),
    legend.text = element_text(color = "black", size = 14, face = "plain"),
    legend.title = element_text(color = "black", size = 14, face = "plain"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(df_summary_clean$total_genesets) * 1.15)

print(p_stacked)


library(ggplot2)
library(dplyr)

# 1. Extract the activity score vector of each method (excluding agglomeration_single)
methods_to_compare <- setdiff(names(de_redundant_results), "agglomeration_single")
score_list <- lapply(methods_to_compare, function(m) de_redundant_results[[m]]$activity_score$activity_score)
score_df <- as.data.frame(score_list)
colnames(score_df) <- methods_to_compare

# 2. Pearson correlation of the other methods against igraph_component as baseline
baseline_scores <- score_df$igraph_component
cor_results <- sapply(colnames(score_df), function(m) {
  if (m == "igraph_component") return(1.0) # baseline correlates with itself at 1
  cor(score_df[[m]], baseline_scores, method = "pearson", use = "complete.obs")
})

# 3. Prepare the plotting data and strip the agglomeration_ prefix
cor_df <- data.frame(
  Method = names(cor_results),
  Correlation = cor_results
)
# Strip the agglomeration_ prefix for display
cor_df$Method <- gsub("agglomeration_", "", cor_df$Method)

# Sort by correlation
cor_df <- cor_df %>% arrange(Correlation)
cor_df$Method <- factor(cor_df$Method, levels = cor_df$Method)

# 4. Draw the horizontal bar chart (14 pt, plain black, regular weight, no title, no bold, x-range 0.98-1)
p_cor <- ggplot(cor_df, aes(x = Method, y = Correlation)) +
  geom_col(fill = "gray50", width = 0.6) +
  coord_flip(ylim = c(0.99, 1)) + # <--- Edit here to restrict the displayed range to 0.98-1
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Method",
    y = "Pearson correlation with baseline"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(size = 14, color = "black", face = "plain"),
    axis.text = element_text(size = 14, color = "black", face = "plain"),
    axis.title = element_text(size = 14, color = "black", face = "plain"),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

print(p_cor)
