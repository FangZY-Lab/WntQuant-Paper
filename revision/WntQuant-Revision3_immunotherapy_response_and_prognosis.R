# ==============================================================================
# WntQuant revision -- Immunotherapy response and prognosis
# ------------------------------------------------------------------------------
# Purpose:    Relate Wnt activity to response (R vs NR) and to OS/PFS across
#             immunotherapy cohorts, including ESTIMATE-adjusted analyses.
# Inputs:     ICI clinical and expression tables under input/Wnt
#             immunotherapy/.
# Outputs:    Response comparisons, Cox forest plots and multi-cohort summary
#             statistics.
# ==============================================================================

################################################################################Prognostic analysis####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")

ici_precog_info=read.csv("ici_precog_all_patient_clinical_data_FOR_DOWNLOAD.csv")


datasets=unique(ici_precog_info$Dataset)
datasets=unique(c(datasets,"Zhao_NatMed_2019","Kim_NatMed_2018","GSE67501","Braun_NatMed_2020",
                  "GSE93157","Gide_Cancercell_2019","GSE115821","GSE106128","GSE100797",
                  "GSE93157","Müller_JCI_2018","GSE100797","GSE135222","GSE19423","GSE100797",
                  "GSE78220","IMVigor210"))
vector = datasets[datasets != "Gide_Cell_2019"]

# Load all data first
for (q in 1:length(vector)) {
  if(file.exists(paste0(vector[q],".Rdata")) && file.exists(paste0(vector[q],"_info.Rdata"))) {
    load(paste0(vector[q],".Rdata"))
    load(paste0(vector[q],"_info.Rdata"))
  } else {
    cat(paste0("Warning: ", vector[q], " data file not found\n"))
  }
}

# Check for valid datasets
valid_vector <- character()
valid_vector_CD274 <- character()

for (q in 1:length(vector)) {
  dataset_name <- vector[q]
  data_obj_name <- dataset_name
  info_obj_name <- paste0(dataset_name, "_info")

  if(exists(info_obj_name) && exists(data_obj_name)) {
    info_obj <- get(info_obj_name)
    data_obj <- get(data_obj_name)

    # Check that the os.time and os.status columns are present
    has_os_time <- "os.time" %in% colnames(info_obj)
    has_os_status <- "os.status" %in% colnames(info_obj)

    if(has_os_time && has_os_status) {
      # Check that neither column is entirely NA
      os_time_valid <- !all(is.na(info_obj$os.time))
      os_status_valid <- !all(is.na(info_obj$os.status))

      if(os_time_valid && os_status_valid) {
        valid_vector <- c(valid_vector, dataset_name)
        cat(paste0("✓ ", dataset_name, ": valid (os.time and os.status present and not all NA)\n"))

        # Require more than 10 samples and the presence of CD274
        sample_size <- nrow(info_obj)

        # Check whether CD274 is present
        has_cd274 <- "CD274" %in% rownames(data_obj)

        if(sample_size > 10 && has_cd274) {
          valid_vector_CD274 <- c(valid_vector_CD274, dataset_name)
          cat(paste0("  OK meets the CD274 criterion: n = ", sample_size, ", CD274 present = ", has_cd274, "\n"))
        } else {
          cat(paste0("  FAILS the CD274 criterion: n = ", sample_size, ", CD274 present = ", has_cd274, "\n"))
        }
      } else {
        cat(paste0("✗ ", dataset_name, ": columns present but all NA\n"))
      }
    } else {
      missing_cols <- c()
      if(!has_os_time) missing_cols <- c(missing_cols, "os.time")
      if(!has_os_status) missing_cols <- c(missing_cols, "os.status")
      cat(paste0("✗ ", dataset_name, ": missing columns - ", paste(missing_cols, collapse = ", "), "\n"))
    }
  } else {
    if(!exists(info_obj_name)) cat(paste0("✗ ", dataset_name, ": no clinical object\n"))
    if(!exists(data_obj_name)) cat(paste0("✗ ", dataset_name, ": no expression object\n"))
  }
}

cat("\nAll valid datasets (", length(valid_vector), "):", paste(valid_vector, collapse = ", "), "\n")
cat("\nDatasets with CD274 and n > 10 (", length(valid_vector_CD274), "):", paste(valid_vector_CD274, collapse = ", "), "\n")

# Save results
save(valid_vector, valid_vector_CD274, file = "filtered_datasets.RData")

# Print detailed statistics
if(length(valid_vector_CD274) > 0) {
  cat("\nDetailed statistics:\n")
  for(dataset in valid_vector_CD274) {
    info_obj <- get(paste0(dataset, "_info"))
    data_obj <- get(dataset)
    sample_size <- nrow(info_obj)
    cat(paste0("- ", dataset, ": n = ", sample_size, ", CD274 expression range = ",
               paste(round(range(data_obj["CD274", ], na.rm = TRUE), 3), collapse = " ~ "), "\n"))
  }
}
quantile_CD274_p.data_all=data.frame()
valid_vector_CD274=valid_vector_CD274[19:length(valid_vector_CD274)]
for(i in 1:length(valid_vector_CD274)){
  exp=get(valid_vector_CD274[i])
  prognosis=get(paste0(valid_vector_CD274[i],"_info"))
  immunotherapy_scores=exp[c("PDCD1","CD274"),]
  immunotherapy_scores=as.data.frame(t(immunotherapy_scores))
  immunotherapy_scores$ID=rownames(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$CD274),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.3, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$CD274,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$CD274 <= threshold,]
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    WNT_score_pro=merge(prognosis,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","os.status","os.time")]
    WNT_score_pro=na.omit(WNT_score_pro)
    colnames(WNT_score_pro)=c("activity_score","status","time")
    WNT_score_pro=WNT_score_pro[,c("activity_score","status","time")]
    quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
    quartiles2=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/2), na.rm = TRUE)
    WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
    WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
    WNT_score_pro$level2 = cut(WNT_score_pro$activity_score, breaks = quartiles2, include.lowest = TRUE,
                               labels = c("Low","High"))
    WNT_score_pro$level2=factor(WNT_score_pro$level2,levels=c("Low","High"))
    library(survival)
    library(survminer)
    fitd3 = survdiff(Surv(time, status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
    p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
    print(p.val3)
    library(survival)
    library(survminer)
    fitd2 = survdiff(Surv(time, status) ~ level2,data= WNT_score_pro,na.action = na.exclude)
    p.val2 = 1 - pchisq(fitd2$chisq, length(fitd2$n) - 1)
    print(p.val2)
    p.data=data.frame(p2=p.val2,p3=p.val3,percent=percent_vector[j])
    p.data_all=rbind(p.data_all,p.data)
  }
  p.data_all$datasets=valid_vector_CD274[i]
  quantile_CD274_p.data_all=rbind(quantile_CD274_p.data_all,p.data_all)
}
save(quantile_CD274_p.data_all,file="quantile_CD274_p.data_all.Rdata")
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("quantile_CD274_p.data_all.Rdata")
# quantile_CD274_p.data_all = quantile_CD274_p.data_all[quantile_CD274_p.data_all$p3 < 0.05, ]
# table(quantile_CD274_p.data_all$datasets)
quantile_CD274_p.data_all=quantile_CD274_p.data_all[,-1]
matrix_data = dcast(quantile_CD274_p.data_all, datasets ~ percent, value.var = "p3")
rownames(matrix_data) = matrix_data$datasets
matrix_data$datasets = NULL
matrix_data = as.matrix(matrix_data)
binary_matrix = ifelse(matrix_data < 0.05, 1, 0)
library(pheatmap)
pheatmap(binary_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "red"),
         legend_breaks = c(0, 1),
         legend_labels = c("p3 ≥ 0.05", "p3 < 0.05"),
         main = "Significance Heatmap (p3 < 0.05)")
################################################################################IMvigor-OS_all####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
IMvigor_info=data.frame(tag=rownames(phenoData),
                        cancer="BLCA",
                        responder=phenoData$binaryResponse,
                        os.time=phenoData$os,
                        os.status=phenoData$censOS,
                        pfs.time=NA,
                        pfs.status=NA)
save(IMvigor_info,file="IMvigor_info.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score$final_activity_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os","censOS")]
wnt_IMvigor_prognosis=na.omit(wnt_IMvigor_prognosis)
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
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
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
################################################################################IMvigor-OS####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score$final_activity_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]
wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`%in%c("desert","excluded"),]
# wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`%in%c("excluded"),]
# wnt_IMvigor=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`%in%c("desert"),]
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os","censOS")]
wnt_IMvigor_prognosis=na.omit(wnt_IMvigor_prognosis)
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
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
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
################################################################################Prognostic analysis: high Wnt activity predicts a better outcome####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score$final_activity_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]
wnt_IMvigor_hot=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`%in%c("inflamed"),]
wnt_IMvigor_hot_prognosis=wnt_IMvigor_hot[,c("activity_score","os","censOS")]
wnt_IMvigor_cold=wnt_IMvigor[wnt_IMvigor$`Immune phenotype`%in%c("desert","excluded"),]
wnt_IMvigor_cold_prognosis=wnt_IMvigor_cold[,c("activity_score","os","censOS")]
wnt_IMvigor_cold_prognosis=na.omit(wnt_IMvigor_cold_prognosis)
quartiles <- quantile(wnt_IMvigor_cold_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_cold_prognosis$level = cut(wnt_IMvigor_cold_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                       labels = c("Cold_low","Cold_medium","Cold_high"))
wnt_IMvigor_hot_prognosis$level="Hot"
wnt_IMvigor_prognosis=rbind(wnt_IMvigor_hot_prognosis,wnt_IMvigor_cold_prognosis)
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Hot","Cold_low","Cold_medium","Cold_high"))

library(survival)
library(survminer)
fitd = survdiff(Surv(os, censOS) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os, censOS)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os, censOS)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#ffd400','#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)

# Debug: inspect the structure of the ps object
print("Structure of the ps object:")
print(str(ps))
print("p-value matrix:")
print(ps$p.value)
print("Row names:")
print(rownames(ps$p.value))
print("Column names:")
print(colnames(ps$p.value))

# Corrected p-value extraction
# Method 1: inspect the whole matrix directly
hot_vs_cold_high_p <- NA

# Check whether these combinations exist in the matrix
if ("Hot" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  hot_vs_cold_high_p <- ps$p.value["Hot", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Hot" %in% colnames(ps$p.value)) {
  hot_vs_cold_high_p <- ps$p.value["Cold_high", "Hot"]
} else {
  # If none of the above works, compute the Hot vs Cold_high p-value by hand
  temp_data <- wnt_IMvigor_prognosis[wnt_IMvigor_prognosis$level %in% c("Hot", "Cold_high"), ]
  temp_data$level <- factor(temp_data$level)
  temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
  hot_vs_cold_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)
}

# Format the p-value labels
hot_vs_cold_high_label <- ifelse(is.na(hot_vs_cold_high_p), "NA",
                                 ifelse(hot_vs_cold_high_p < 0.001, "p < 0.001",
                                        ifelse(hot_vs_cold_high_p < 0.01, paste0("p = ", round(hot_vs_cold_high_p, 3)),
                                               ifelse(hot_vs_cold_high_p < 0.05, paste0("p = ", round(hot_vs_cold_high_p, 3)), "ns"))))

p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Hot","Cold_low","Cold_medium","Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS","(Versus hot tumors,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic") +
  # Add the Hot vs Cold_high significance annotation
  annotate("text", x = 0, y = 0.08,
           label = paste0("Hot vs. Cold_high: ", hot_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic",color="red")
print(p)
################################################################################Immunotherapy response####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score$final_activity_score,phenoData,by="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","binaryResponse","Immune phenotype")]
colnames(wnt_IMvigor)=c("activity_score","Response","Immune phenotype")
wnt_IMvigor=na.omit(wnt_IMvigor)
wnt_IMvigor$`Immune phenotype`=gsub("^(\\w)(.*)", "\\U\\1\\L\\2", wnt_IMvigor$`Immune phenotype`, perl = TRUE)
p = ggplot(wnt_IMvigor, aes(`Immune phenotype`, activity_score, fill = Response, color = Response)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +
  theme_bw() +
  theme(
    legend.position = "top",
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
################################################################################Analysis of R vs NR####
rm(list=ls())
gc()

# Set the working directory and load packages
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
library(reshape2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(GSEABase)
library(dplyr)

# Read the data and build the dataset list
ici_precog_info <- read.csv("ici_precog_all_patient_clinical_data_FOR_DOWNLOAD.csv")
datasets <- unique(ici_precog_info$Dataset)
datasets <- unique(c(datasets, "Zhao_NatMed_2019", "Kim_NatMed_2018", "GSE67501", "Braun_NatMed_2020",
                     "GSE93157", "Gide_CancerCell_2019", "GSE115821", "GSE106128", "GSE100797",
                     "GSE93157", "Müller_JCI_2018", "GSE100797", "GSE135222", "GSE19423", "GSE100797",
                     "GSE78220", "IMVigor210"))
vector <- datasets[datasets != "Gide_Cell_2019"]

# Check for valid datasets
vector_responder <- c()

for (q in 1:length(vector)) {
  data_file <- paste0(vector[q], ".Rdata")
  info_file <- paste0(vector[q], "_info.Rdata")

  if(file.exists(data_file) && file.exists(info_file)) {
    load(data_file)
    load(info_file)

    info_var_name <- paste0(vector[q], "_info")

    if(exists(info_var_name)) {
      info_data <- get(info_var_name)

      if("responder" %in% colnames(info_data)) {
        responder_values <- info_data$responder
        if(!all(is.na(responder_values)) && length(responder_values) > 0) {
          vector_responder <- c(vector_responder, vector[q])
          cat(paste0("OK dataset ", vector[q], " has a valid responder column\n"))
        } else {
          cat(paste0("WARNING: ", vector[q], " has a responder column that is entirely NA or empty\n"))
        }
      } else {
        cat(paste0("WARNING: ", vector[q], " has no responder column\n"))
      }
    } else {
      cat(paste0("WARNING: ", vector[q], " failed to load its info table\n"))
    }
  } else {
    cat(paste0("WARNING: ", vector[q], " data file not found\n"))
  }
}

cat("\nDatasets with a valid responder column:\n")
print(vector_responder)
cat("Number of valid datasets:", length(vector_responder), "\n")

# Stop if no valid dataset is available
if(length(vector_responder) == 0) {
  stop("No dataset has a valid responder column")
}

# Load the required functions and gene sets
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
geneSets <- getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")

# Compute the Wnt activity score
all_WNT_scores <- list()

for(i in 1:length(vector_responder)) {
  dataset_name <- vector_responder[i]
  cat("Processing dataset:", dataset_name, "\n")

  tryCatch({
    # Get the expression data and clinical information
    exp_data <- get(dataset_name)
    info_data <- get(paste0(dataset_name, "_info"))

    # Compute the Wnt activity score
    WNT_Score <- Calculate_bioactivity_scores(
      file_paths = "C:/Users/user/Desktop/input",
      expression_profile = exp_data,
      foundation = "relative_ssGSEA",
      activation_geneset = NA,
      inhibition_geneset = NA,
      geneSets_gmt = geneSets,
      min.sz = 2,
      max.sz = 10000,
      geneset_direction = c(
        rep("pos", sum(grepl("WPAGS", names(geneSets)))),
        rep("neg", sum(grepl("WPIGS", names(geneSets))))
      )
    )

    WNT_score <- WNT_Score$final_activity_score

    # Merge the response information
    responder_info <- info_data[, c("tag", "responder")]
    WNT_score_pro <- merge(responder_info, WNT_score, by.x = "tag", by.y = "row.names")

    # Data cleaning and formatting
    WNT_score_pro$responder <- factor(WNT_score_pro$responder, levels = c("TRUE", "FALSE"))
    WNT_score_pro <- WNT_score_pro[, c("responder", "activity_score")]
    WNT_score_pro$dataset <- dataset_name

    # Remove NA values
    WNT_score_pro <- na.omit(WNT_score_pro)

    # Save results
    all_WNT_scores[[dataset_name]] <- WNT_score_pro
    cat("OK processed dataset:", dataset_name, "\n")

  }, error = function(e) {
    cat("FAILED on dataset", dataset_name, ": error:", e$message, "\n")
  })
}

# Merge all data
combined_data <- do.call(rbind, all_WNT_scores)
combined_data <- na.omit(combined_data)

# Sample size
sample_sizes <- data.frame(
  dataset = unique(combined_data$dataset)
)

sample_sizes$total_n <- sapply(sample_sizes$dataset, function(ds) {
  nrow(combined_data[combined_data$dataset == ds, ])
})

sample_sizes$responder_n <- sapply(sample_sizes$dataset, function(ds) {
  sum(combined_data$dataset == ds & combined_data$responder == "TRUE")
})

sample_sizes$non_responder_n <- sapply(sample_sizes$dataset, function(ds) {
  sum(combined_data$dataset == ds & combined_data$responder == "FALSE")
})

sample_sizes$label <- paste0("n=", sample_sizes$total_n,
                             " (R:", sample_sizes$responder_n,
                             "/NR:", sample_sizes$non_responder_n, ")")

# Build new x-axis labels that include the sample size
new_labels <- setNames(
  paste0(sample_sizes$dataset, "\n(n=", sample_sizes$total_n, ")"),
  sample_sizes$dataset
)
new_labels <- setNames(
  paste0(sample_sizes$dataset, "(n=", sample_sizes$total_n, ")"),
  sample_sizes$dataset
)
# Create the figure
p <- ggplot(combined_data, aes(x = dataset, y = activity_score, fill = responder)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.fill = 'black',
    outlier.size = 0.25,
    color = "black",
    alpha = 0.8
  ) +
  labs(
    x = "",
    y = "WNT Activity Score",
    title = "WNT Activity Score by Responder Status Across Datasets",
    fill = "Response Status"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      colour = "black",
      size = 9  # Reduce the font size slightly to fit more content
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.title.x = element_text(size = 12, margin = margin(t = 25)),  # Increase the bottom margin
    plot.title = element_text(
      size = 14,
      colour = "black",
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    legend.position = "top",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
    labels = c("TRUE" = "Responder", "FALSE" = "Non-responder")
  ) +
  scale_x_discrete(labels = new_labels) +  # Use the new labels that include the sample size
  # Add the statistical test
  stat_compare_means(
    aes(group = responder, label = ..p.signif..),
    method = "wilcox.test",
    size = 4,
    color = "black",
    vjust = -0.5
  )

# Display the figure
print(p)
################################################################################Analysis restricted to datasets with more than 90 samples (immune response data)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020","CheckMate_Meta","Gide_CancerCell_2019","Goswami_2025","GSE159067",
         "GSE94873","IMmotion150","IMVigor210","Kinget_2024","Liu_NatMed_2019","phs000452",
         "Ravi_NatGen_2023","Riaz_Nivolumab_2017")
filtered_vector <- c()
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  exp=get(vector[q])
  if("CD274" %in% rownames(exp)) {
    filtered_vector = c(filtered_vector, vector[q])
    print(paste0(vector[q], ": ", length(colnames(exp)), " samples"))
    load(paste0(vector[q],"_info.Rdata"))
  }
}
vector = filtered_vector
wilcox_p.data_all=data.frame()
for(i in 1:length(vector)){
  exp=get(vector[i])
  immunotherapy_scores=exp[c("PDCD1","CD274"),]
  immunotherapy_scores=as.data.frame(t(immunotherapy_scores))
  immunotherapy_scores$ID=rownames(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$CD274),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.1, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$CD274,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$CD274 <= threshold,]
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    responder=get(paste0(vector[i],"_info"))
    WNT_score_pro=merge(responder,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","responder")]
    WNT_score_pro=na.omit(WNT_score_pro)
    WNT_score_pro$responder=factor(WNT_score_pro$responder,levels = c("TRUE","FALSE"))
    if(length(unique(WNT_score_pro$responder)) >= 2 &&
       all(table(WNT_score_pro$responder) > 3)){
      p_two=wilcox.test(activity_score ~ responder, data = WNT_score_pro)
      p_two=p_two$p.value
      p_greater=wilcox.test(activity_score ~ responder, data = WNT_score_pro,alternative = c("greater"))
      p_greater=p_greater$p.value
      p_less=wilcox.test(activity_score ~ responder, data = WNT_score_pro,alternative = c("less"))
      p_less=p_less$p.value
      p.data=data.frame(p_two=p_two,p_greater=p_greater,p_less=p_less,percent=percent_vector[j])
      p.data_all=rbind(p.data_all,p.data)
      print(paste0(p_two,"/",p_greater,"/",p_less))
    }else{
      p.data=data.frame(p_two=1,p_greater=1,p_less=1,percent=percent_vector[j])
      p.data_all=rbind(p.data_all,p.data)
      print(paste0("1","/","1","/","1"))
    }
  }
  p.data_all$datasets=vector[i]
  wilcox_p.data_all=rbind(wilcox_p.data_all,p.data_all)
}
save(wilcox_p.data_all,file="wilcox_p.data_all.Rdata")
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("wilcox_p.data_all.Rdata")
wilcox_p.data_all$pstar=ifelse(wilcox_p.data_all$p_two > 0.05, "",
                                        ifelse(wilcox_p.data_all$p_two <= 0.05 & wilcox_p.data_all$p_two > 0.01, "*",
                                               ifelse(wilcox_p.data_all$p_two <= 0.01 & wilcox_p.data_all$p_two > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
wilcox_p.data_all$percent=as.factor(wilcox_p.data_all$percent)
ggplot(wilcox_p.data_all, aes(percent, datasets)) +
  geom_tile(aes(fill = -log10(p_two)), color = "black") +  # Add a black box
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 3) +  # 10 pt font, roughly 3 mm
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "-log10(P)")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )+ggtitle("CD274_response")
################################################################################Analysis restricted to datasets with more than 90 samples (OS immune prognosis data)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020","CheckMate_Meta","Gide_CancerCell_2019","Goswami_2025","GSE159067",
         "GSE94873","IMmotion150","IMVigor210","Kinget_2024","Liu_NatMed_2019","phs000452",
         "Ravi_NatGen_2023","Riaz_Nivolumab_2017")
filtered_vector <- c()
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  exp <- get(vector[q])

  # Check whether CD274 is present
  if("CD274" %in% rownames(exp)) {
    # Load the matching info file
    load(paste0(vector[q],"_info.Rdata"))

    # Get the name of the info data frame (assumed to be vector[q] plus an _info suffix)
    info_name <- paste0(vector[q], "_info")
    info_data <- get(info_name)

    # Check that os.time and os.status exist and are not entirely NA
    if(all(c("os.time", "os.status") %in% colnames(info_data))) {
      # Check that neither column is entirely NA
      if(!all(is.na(info_data$os.time)) && !all(is.na(info_data$os.status))) {
        filtered_vector <- c(filtered_vector, vector[q])
        print(paste0(vector[q], ": ", length(colnames(exp)), " samples"))
      } else {
        print(paste0(vector[q], ": CD274 present but os.time or os.status all NA"))
      }
    } else {
      print(paste0(vector[q], ": CD274 present but os.time or os.status missing"))
    }
  }
}
vector = filtered_vector
os_p.data_all=data.frame()
for(i in 1:length(vector)){
  exp=get(vector[i])
  prognosis=get(paste0(vector[i],"_info"))
  immunotherapy_scores=exp[c("PDCD1","CD274"),]
  immunotherapy_scores=as.data.frame(t(immunotherapy_scores))
  immunotherapy_scores$ID=rownames(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$CD274),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.1, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$CD274,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$CD274 <= threshold,]
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    WNT_score_pro=merge(prognosis,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","os.status","os.time")]
    WNT_score_pro=na.omit(WNT_score_pro)
    colnames(WNT_score_pro)=c("activity_score","status","time")
    WNT_score_pro=WNT_score_pro[,c("activity_score","status","time")]
    quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
    quartiles2=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/2), na.rm = TRUE)
    WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
    WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
    WNT_score_pro$level2 = cut(WNT_score_pro$activity_score, breaks = quartiles2, include.lowest = TRUE,
                               labels = c("Low","High"))
    WNT_score_pro$level2=factor(WNT_score_pro$level2,levels=c("Low","High"))
    library(survival)
    library(survminer)
    fitd3 = survdiff(Surv(time, status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
    p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
    print(p.val3)
    library(survival)
    library(survminer)
    fitd2 = survdiff(Surv(time, status) ~ level2,data= WNT_score_pro,na.action = na.exclude)
    p.val2 = 1 - pchisq(fitd2$chisq, length(fitd2$n) - 1)
    print(p.val2)
    p.data=data.frame(p2=p.val2,p3=p.val3,percent=percent_vector[j])
    p.data_all=rbind(p.data_all,p.data)
  }
  p.data_all$datasets=vector[i]
  os_p.data_all=rbind(os_p.data_all,p.data_all)
}
save(os_p.data_all,file="os_p.data_all.Rdata")
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("os_p.data_all.Rdata")
os_p.data_all$pstar=ifelse(os_p.data_all$p3 > 0.05, "",
                                    ifelse(os_p.data_all$p3 <= 0.05 & os_p.data_all$p3 > 0.01, "*",
                                           ifelse(os_p.data_all$p3 <= 0.01 & os_p.data_all$p3 > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
os_p.data_all$percent=as.factor(os_p.data_all$percent)
ggplot(os_p.data_all, aes(percent, datasets)) +
  geom_tile(aes(fill = -log10(p3)), color = "black") +  # Add a black box
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 3) +  # 10 pt font, roughly 3 mm
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "-log10(P)")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )+
  ggtitle("CD274_OS_p3")
################################################################################Analysis restricted to datasets with more than 90 samples (PFS immune prognosis data)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020","CheckMate_Meta","Gide_CancerCell_2019","Goswami_2025","GSE159067",
         "GSE94873","IMmotion150","IMVigor210","Kinget_2024","Liu_NatMed_2019","phs000452",
         "Ravi_NatGen_2023","Riaz_Nivolumab_2017")
filtered_vector <- c()
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  exp <- get(vector[q])

  # Check whether CD274 is present
  if("CD274" %in% rownames(exp)) {
    # Load the matching info file
    load(paste0(vector[q],"_info.Rdata"))

    # Get the name of the info data frame (assumed to be vector[q] plus an _info suffix)
    info_name <- paste0(vector[q], "_info")
    info_data <- get(info_name)

    # Check that pfs.time and pfs.status exist and are not entirely NA
    if(all(c("pfs.time", "pfs.status") %in% colnames(info_data))) {
      # Check that neither column is entirely NA
      if(!all(is.na(info_data$pfs.time)) && !all(is.na(info_data$pfs.status))) {
        filtered_vector <- c(filtered_vector, vector[q])
        print(paste0(vector[q], ": ", length(colnames(exp)), " samples"))
      } else {
        print(paste0(vector[q], ": CD274 present but pfs.time or pfs.status all NA"))
      }
    } else {
      print(paste0(vector[q], ": CD274 present but pfs.time or pfs.status missing"))
    }
  }
}
vector = filtered_vector
pfs_p.data_all=data.frame()
for(i in 1:length(vector)){
  exp=get(vector[i])
  prognosis=get(paste0(vector[i],"_info"))
  immunotherapy_scores=exp[c("PDCD1","CD274"),]
  immunotherapy_scores=as.data.frame(t(immunotherapy_scores))
  immunotherapy_scores$ID=rownames(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$CD274),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.1, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$CD274,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$CD274 <= threshold,]
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    WNT_score_pro=merge(prognosis,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","pfs.status","pfs.time")]
    WNT_score_pro=na.omit(WNT_score_pro)
    colnames(WNT_score_pro)=c("activity_score","status","time")
    WNT_score_pro=WNT_score_pro[,c("activity_score","status","time")]
    quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
    quartiles2=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/2), na.rm = TRUE)
    WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
    WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
    WNT_score_pro$level2 = cut(WNT_score_pro$activity_score, breaks = quartiles2, include.lowest = TRUE,
                               labels = c("Low","High"))
    WNT_score_pro$level2=factor(WNT_score_pro$level2,levels=c("Low","High"))
    library(survival)
    library(survminer)
    fitd3 = survdiff(Surv(time, status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
    p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
    print(p.val3)
    library(survival)
    library(survminer)
    fitd2 = survdiff(Surv(time, status) ~ level2,data= WNT_score_pro,na.action = na.exclude)
    p.val2 = 1 - pchisq(fitd2$chisq, length(fitd2$n) - 1)
    print(p.val2)
    p.data=data.frame(p2=p.val2,p3=p.val3,percent=percent_vector[j])
    p.data_all=rbind(p.data_all,p.data)
  }
  p.data_all$datasets=vector[i]
  pfs_p.data_all=rbind(pfs_p.data_all,p.data_all)
}
save(pfs_p.data_all,file="pfs_p.data_all.Rdata")
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_p.data_all.Rdata")
pfs_p.data_all$pstar=ifelse(pfs_p.data_all$p2 > 0.05, "",
                                     ifelse(pfs_p.data_all$p2 <= 0.05 & pfs_p.data_all$p2 > 0.01, "*",
                                            ifelse(pfs_p.data_all$p2 <= 0.01 & pfs_p.data_all$p2 > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
pfs_p.data_all$percent=as.factor(pfs_p.data_all$percent)
ggplot(pfs_p.data_all, aes(percent, datasets)) +
  geom_tile(aes(fill = -log10(p2)), color = "black") +  # Add a black box
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 3) +  # 10 pt font, roughly 3 mm
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "-log10(P)")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )+
  ggtitle("CD274_PFS_p2")
################################################################################Analysis restricted to datasets with more than 90 samples (immune response data, ESTIMATE)####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020","CheckMate_Meta","Gide_CancerCell_2019","Goswami_2025","GSE159067",
         "GSE94873","IMmotion150","IMVigor210","Kinget_2024","Liu_NatMed_2019","phs000452",
         "Ravi_NatGen_2023","Riaz_Nivolumab_2017")

filtered_vector <- c()
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  exp <- get(vector[q])
  load(paste0(vector[q],"_info.Rdata"))
  info_name <- paste0(vector[q], "_info")
  info_data <- get(info_name)

  # Run only two checks
  column_exists <- "responder" %in% colnames(info_data)
  data_valid <- !all(is.na(info_data$responder))

  if(column_exists && data_valid) {
    filtered_vector <- c(filtered_vector, vector[q])
    print(paste0(vector[q], ": passed the check, ", length(colnames(exp)), " samples"))
  }
}
estimate_wilcox_p.data_all=data.frame()
for(i in 1:length(vector)){
  exp=get(vector[i])
  immunotherapy_scores=deconvo_tme(eset = exp, method = "estimate")
  immunotherapy_scores=as.data.frame(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$ImmuneScore_estimate),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.1, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$ImmuneScore_estimate,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$ImmuneScore_estimate <= threshold,]
    if(vector[i]%in%c("GSE159067","Goswami_2025")){
      estimate_cal$ID = gsub("-", ".", gsub("^X", "", estimate_cal$ID))
    }
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    responder=get(paste0(vector[i],"_info"))
    WNT_score_pro=merge(responder,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","responder")]
    WNT_score_pro=na.omit(WNT_score_pro)
    WNT_score_pro$responder=factor(WNT_score_pro$responder,levels = c("TRUE","FALSE"))
    if(length(unique(WNT_score_pro$responder)) >= 2 &&
       all(table(WNT_score_pro$responder) > 3)){
      p_two=wilcox.test(activity_score ~ responder, data = WNT_score_pro)
      p_two=p_two$p.value
      p_greater=wilcox.test(activity_score ~ responder, data = WNT_score_pro,alternative = c("greater"))
      p_greater=p_greater$p.value
      p_less=wilcox.test(activity_score ~ responder, data = WNT_score_pro,alternative = c("less"))
      p_less=p_less$p.value
      p.data=data.frame(p_two=p_two,p_greater=p_greater,p_less=p_less,percent=percent_vector[j])
      p.data_all=rbind(p.data_all,p.data)
      print(paste0(p_two,"/",p_greater,"/",p_less))
    }else{
      p.data=data.frame(p_two=1,p_greater=1,p_less=1,percent=percent_vector[j])
      p.data_all=rbind(p.data_all,p.data)
      print(paste0("1","/","1","/","1"))
    }
  }
  p.data_all$datasets=vector[i]
  estimate_wilcox_p.data_all=rbind(estimate_wilcox_p.data_all,p.data_all)
}
save(estimate_wilcox_p.data_all,file="estimate_wilcox_p.data_all.Rdata")
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("estimate_wilcox_p.data_all.Rdata")
estimate_wilcox_p.data_all$pstar=ifelse(estimate_wilcox_p.data_all$p_two > 0.05, "",
                        ifelse(estimate_wilcox_p.data_all$p_two <= 0.05 & estimate_wilcox_p.data_all$p_two > 0.01, "*",
                               ifelse(estimate_wilcox_p.data_all$p_two <= 0.01 & estimate_wilcox_p.data_all$p_two > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
estimate_wilcox_p.data_all$percent=as.factor(estimate_wilcox_p.data_all$percent)
ggplot(estimate_wilcox_p.data_all, aes(percent, datasets)) +
  geom_tile(aes(fill = -log10(p_two)), color = "black") +  # Add a black box
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 3) +  # 10 pt font, roughly 3 mm
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "-log10(P)")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )+ggtitle("ImmuneScore_response")
################################################################################Analysis restricted to datasets with more than 90 samples (OS immune prognosis data, ESTIMATE)####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020","CheckMate_Meta","Gide_CancerCell_2019","Goswami_2025","GSE159067",
         "GSE94873","IMmotion150","IMVigor210","Kinget_2024","Liu_NatMed_2019","phs000452",
         "Ravi_NatGen_2023","Riaz_Nivolumab_2017")
filtered_vector <- c()
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  exp <- get(vector[q])
  load(paste0(vector[q],"_info.Rdata"))
  info_name <- paste0(vector[q], "_info")
  info_data <- get(info_name)

  # Run only two checks
  columns_exist <- all(c("os.time", "os.status") %in% colnames(info_data))
  data_valid <- !all(is.na(info_data$os.time)) && !all(is.na(info_data$os.status))

  if(columns_exist && data_valid) {
    filtered_vector <- c(filtered_vector, vector[q])
    print(paste0(vector[q], ": passed the check, ", length(colnames(exp)), " samples"))
  }
}
vector = filtered_vector
estimate_os_p.data_all=data.frame()
for(i in 1:length(vector)){
  exp=get(vector[i])
  prognosis=get(paste0(vector[i],"_info"))
  immunotherapy_scores=deconvo_tme(eset = exp, method = "estimate")
  immunotherapy_scores=as.data.frame(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$ImmuneScore_estimate),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.1, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$ImmuneScore_estimate,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$ImmuneScore_estimate <= threshold,]
    if(vector[i]%in%c("GSE159067","Goswami_2025")){
      estimate_cal$ID = gsub("-", ".", gsub("^X", "", estimate_cal$ID))
    }
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    WNT_score_pro=merge(prognosis,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","os.status","os.time")]
    WNT_score_pro=na.omit(WNT_score_pro)
    colnames(WNT_score_pro)=c("activity_score","status","time")
    WNT_score_pro=WNT_score_pro[,c("activity_score","status","time")]
    quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
    quartiles2=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/2), na.rm = TRUE)
    WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
    WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
    WNT_score_pro$level2 = cut(WNT_score_pro$activity_score, breaks = quartiles2, include.lowest = TRUE,
                               labels = c("Low","High"))
    WNT_score_pro$level2=factor(WNT_score_pro$level2,levels=c("Low","High"))
    library(survival)
    library(survminer)
    fitd3 = survdiff(Surv(time, status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
    p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
    print(p.val3)
    library(survival)
    library(survminer)
    fitd2 = survdiff(Surv(time, status) ~ level2,data= WNT_score_pro,na.action = na.exclude)
    p.val2 = 1 - pchisq(fitd2$chisq, length(fitd2$n) - 1)
    print(p.val2)
    p.data=data.frame(p2=p.val2,p3=p.val3,percent=percent_vector[j])
    p.data_all=rbind(p.data_all,p.data)
  }
  p.data_all$datasets=vector[i]
  estimate_os_p.data_all=rbind(estimate_os_p.data_all,p.data_all)
}
save(estimate_os_p.data_all,file="estimate_os_p.data_all.Rdata")
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("estimate_os_p.data_all.Rdata")
estimate_os_p.data_all$pstar=ifelse(estimate_os_p.data_all$p2 > 0.05, "",
                                        ifelse(estimate_os_p.data_all$p2 <= 0.05 & estimate_os_p.data_all$p2 > 0.01, "*",
                                               ifelse(estimate_os_p.data_all$p2 <= 0.01 & estimate_os_p.data_all$p2 > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
estimate_os_p.data_all$percent=as.factor(estimate_os_p.data_all$percent)
ggplot(estimate_os_p.data_all, aes(percent, datasets)) +
  geom_tile(aes(fill = -log10(p2)), color = "black") +  # Add a black box
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 3) +  # 10 pt font, roughly 3 mm
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "-log10(P)")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )+ggtitle("ImmuneScore_OS_p2")
################################################################################Analysis restricted to datasets with more than 90 samples (PFS immune prognosis data, ESTIMATE)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020","CheckMate_Meta","Gide_CancerCell_2019","Goswami_2025","GSE159067",
         "GSE94873","IMmotion150","IMVigor210","Kinget_2024","Liu_NatMed_2019","phs000452",
         "Ravi_NatGen_2023","Riaz_Nivolumab_2017")
filtered_vector <- c()
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  exp <- get(vector[q])
  load(paste0(vector[q],"_info.Rdata"))
  info_name <- paste0(vector[q], "_info")
  info_data <- get(info_name)

  # Run only two checks
  columns_exist <- all(c("pfs.time", "pfs.status") %in% colnames(info_data))
  data_valid <- !all(is.na(info_data$pfs.time)) && !all(is.na(info_data$pfs.status))

  if(columns_exist && data_valid) {
    filtered_vector <- c(filtered_vector, vector[q])
    print(paste0(vector[q], ": passed the check, ", length(colnames(exp)), " samples"))
  }
}
vector = filtered_vector
estimate_pfs_p.data_all=data.frame()
for(i in 1:length(vector)){
  exp=get(vector[i])
  prognosis=get(paste0(vector[i],"_info"))
  immunotherapy_scores=deconvo_tme(eset = exp, method = "estimate")
  immunotherapy_scores=as.data.frame(immunotherapy_scores)
  immunotherapy_scores=immunotherapy_scores[order(immunotherapy_scores$ImmuneScore_estimate),]
  p.data_all=data.frame()
  percent_vector = seq(from = 0.1, to = 1, by = 0.01)
  for(j in 1:length(percent_vector)){
    print(percent_vector[j])
    threshold = quantile(immunotherapy_scores$ImmuneScore_estimate,
                         probs = percent_vector[j],
                         na.rm = TRUE)
    print(threshold)
    estimate_cal = immunotherapy_scores[immunotherapy_scores$ImmuneScore_estimate <= threshold,]
    if(vector[i]%in%c("GSE159067","Goswami_2025")){
      estimate_cal$ID = gsub("-", ".", gsub("^X", "", estimate_cal$ID))
    }
    exp_cal=exp[,intersect(estimate_cal$ID,colnames(exp))]
    library(GSEABase)
    geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
    source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
    WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                           expression_profile=exp_cal,
                                           foundation="relative_ssGSEA",
                                           activation_geneset=NA,
                                           inhibition_geneset=NA,
                                           geneSets_gmt=geneSets,
                                           min.sz=1,
                                           max.sz=10000,
                                           geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                               rep("neg",sum(grepl("WPIGS", names(geneSets))))))
    WNT_score=WNT_Score$final_activity_score
    WNT_score_pro=merge(prognosis,WNT_score,by.x="tag",by.y="row.names")
    WNT_score_pro=WNT_score_pro[,c("activity_score","pfs.status","pfs.time")]
    WNT_score_pro=na.omit(WNT_score_pro)
    colnames(WNT_score_pro)=c("activity_score","status","time")
    WNT_score_pro=WNT_score_pro[,c("activity_score","status","time")]
    quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
    quartiles2=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/2), na.rm = TRUE)
    WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
    WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
    WNT_score_pro$level2 = cut(WNT_score_pro$activity_score, breaks = quartiles2, include.lowest = TRUE,
                               labels = c("Low","High"))
    WNT_score_pro$level2=factor(WNT_score_pro$level2,levels=c("Low","High"))
    library(survival)
    library(survminer)
    fitd3 = survdiff(Surv(time, status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
    p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
    print(p.val3)
    library(survival)
    library(survminer)
    fitd2 = survdiff(Surv(time, status) ~ level2,data= WNT_score_pro,na.action = na.exclude)
    p.val2 = 1 - pchisq(fitd2$chisq, length(fitd2$n) - 1)
    print(p.val2)
    p.data=data.frame(p2=p.val2,p3=p.val3,percent=percent_vector[j])
    p.data_all=rbind(p.data_all,p.data)
  }
  p.data_all$datasets=vector[i]
  estimate_pfs_p.data_all=rbind(estimate_pfs_p.data_all,p.data_all)
}
save(estimate_pfs_p.data_all,file="estimate_pfs_p.data_all.Rdata")
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("estimate_pfs_p.data_all.Rdata")
estimate_pfs_p.data_all$pstar=ifelse(estimate_pfs_p.data_all$p3 > 0.05, "",
                                    ifelse(estimate_pfs_p.data_all$p3 <= 0.05 & estimate_pfs_p.data_all$p3 > 0.01, "*",
                                           ifelse(estimate_pfs_p.data_all$p3 <= 0.01 & estimate_pfs_p.data_all$p3 > 0.001, "**", "***")))

library(ggplot2)
library(dplyr)
estimate_pfs_p.data_all$percent=as.factor(estimate_pfs_p.data_all$percent)
ggplot(estimate_pfs_p.data_all, aes(percent, datasets)) +
  geom_tile(aes(fill = -log10(p3)), color = "black") +  # Add a black box
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +
  geom_text(aes(label = pstar), color = "black", size = 3) +  # 10 pt font, roughly 3 mm
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "-log10(P)")) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )+ggtitle("ImmuneScore_PFS_p3")
################################################################################Supervised clustering: lasso####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("C:/Users/user/Desktop/input/IMvigor.Rdata")
IMvigor=log2(IMvigor+1)
# IMvigor_G=data.frame(Tag=rownames(phenoData),group=phenoData$`Immune phenotype`)
# IMvigor_G=na.omit(IMvigor_G)
# IMvigor_G$group=ifelse(IMvigor_G$group=="inflamed","hot","cold")
# IMvigor=IMvigor[,IMvigor_G$Tag]
# save(IMvigor,file="IMvigor.Rdata")
# save(IMvigor_G,file="IMvigor_G.Rdata")
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
################################################################################Building the lasso regression####
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/xcell_pheno.Rdata")

library(glmnet)
library(caret)
library(pROC)
library(ggplot2)
library(vip)
library(dplyr)
library(purrr)
library(ggthemes)
library(scales)
library(RColorBrewer)

theme_set(theme_minimal(base_size = 14))
my_colors <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#3B1F2B")

# Data preprocessing
xcell_pheno$phenotype <- ifelse(xcell_pheno$phenotype == "hot", 1, 0)
xcell_pheno$phenotype <- as.factor(xcell_pheno$phenotype)
processed_data <- xcell_pheno
colnames(processed_data)[1] <- "target"

# Prepare matrix format
x <- as.matrix(processed_data[, -1])
y <- processed_data$target

# Create factor response variable (required by caret)
processed_data$target_factor <- factor(ifelse(processed_data$target == 1, "hot", "cold"),
                                       levels = c("hot", "cold"))

# Set up cross-validation
set.seed(123)
train_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  selectionFunction = "best",
  savePredictions = TRUE,
  verboseIter = FALSE
)

# Train LASSO model
cat("Training LASSO model...\n")
lasso_model <- train(
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

# View best parameters
cat("\n=== Best Tuning Results ===\n")
print(lasso_model$bestTune)

# 1. Tuning results plot
tune_results <- lasso_model$results
p1 <- ggplot(tune_results, aes(x = lambda, y = ROC)) +
  geom_point(size = 3, color = my_colors[1], alpha = 0.7) +
  geom_line(color = my_colors[1], linewidth = 1.2, alpha = 0.8) +
  geom_vline(xintercept = lasso_model$bestTune$lambda,
             linetype = "dashed", color = my_colors[3], linewidth = 1.2) +
  geom_point(data = tune_results[tune_results$lambda == lasso_model$bestTune$lambda, ],
             aes(x = lambda, y = ROC),
             color = my_colors[3], size = 5, shape = 18) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "LASSO model tuning results",
       subtitle = paste("Optimal λ =", round(lasso_model$bestTune$lambda, 5),
                        "| AUC =", round(max(tune_results$ROC), 4)),
       x = "Regularization parameter (λ)",
       y = "Cross-validated AUC") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

print(p1)

# 2. Coefficient paths plot
coef_data <- as.matrix(coef(lasso_model$finalModel, s = lasso_model$finalModel$lambda))
coef_df <- data.frame(
  Lambda = rep(lasso_model$finalModel$lambda, each = nrow(coef_data)),
  Variable = rep(rownames(coef_data), times = ncol(coef_data)),
  Value = as.vector(coef_data)
)
coef_df <- coef_df[coef_df$Variable != "(Intercept)", ]

# Full coefficient paths
p2_full <- ggplot(coef_df, aes(x = log(Lambda), y = Value, group = Variable)) +
  geom_line(aes(color = abs(Value)), alpha = 0.8, linewidth = 0.8) +
  scale_color_gradientn(colors = c("#4575B4", "#91BFDB", "#E0F3F8", "#FEE090", "#FC8D59", "#D73027"),
                        name = "|Coefficient|") +
  geom_vline(xintercept = log(lasso_model$bestTune$lambda),
             linetype = "dashed", color = "#E31A1C", linewidth = 1.2) +
  scale_x_reverse() +
  labs(title = "LASSO coefficient paths",
       x = expression(log[10](λ)),
       y = "Coefficient value") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

print(p2_full)

# Selected variables coefficient paths
final_coef <- coef(lasso_model$finalModel, s = lasso_model$bestTune$lambda)
selected_vars <- rownames(final_coef)[which(final_coef != 0)]
selected_vars <- selected_vars[selected_vars != "(Intercept)"]

cat("\nNumber of selected variables:", length(selected_vars), "\n")

if(length(selected_vars) > 0) {
  coef_df_selected <- coef_df[coef_df$Variable %in% selected_vars, ]

  # Choose appropriate color scheme based on number of variables
  if(length(selected_vars) <= 8) {
    color_scheme <- brewer.pal(max(3, length(selected_vars)), "Set2")
  } else if(length(selected_vars) <= 12) {
    color_scheme <- brewer.pal(length(selected_vars), "Set3")
  } else {
    color_scheme <- colorRampPalette(brewer.pal(8, "Set2"))(length(selected_vars))
  }

  p2_selected <- ggplot(coef_df_selected, aes(x = log(Lambda), y = Value, color = Variable)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_vline(xintercept = log(lasso_model$bestTune$lambda),
               linetype = "dashed", color = "#E31A1C", linewidth = 1.2, alpha = 0.8) +
    scale_x_reverse() +
    scale_color_manual(values = color_scheme) +
    labs(title = "LASSO coefficient paths",
         subtitle = paste(length(selected_vars), "selected variables"),
         x = expression(log[10](λ)),
         y = "Coefficient value",
         color = "Variables") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )

  # Simplify legend if too many variables
  if(length(selected_vars) > 10) {
    p2_selected <- p2_selected +
      theme(legend.position = "none") +
      labs(title = "LASSO coefficient paths",
           subtitle = paste(length(selected_vars), "selected variables (legend omitted)"))
  }

  print(p2_selected)
} else {
  cat("No variables selected (null model)\n")
}

# 3. Variable importance plot
importance_data <- vi(lasso_model)
top_n <- min(8, nrow(importance_data))

p3 <- importance_data %>%
  slice_max(Importance, n = top_n) %>%
  mutate(Variable = factor(Variable, levels = Variable[order(Importance)])) %>%
  ggplot(aes(x = Importance, y = Variable)) +
  geom_segment(aes(xend = 0, yend = Variable),
               color = "gray70", linewidth = 1) +
  geom_point(aes(color = Importance), size = 4, alpha = 0.8) +
  scale_color_gradient(low = my_colors[1], high = my_colors[3]) +
  labs(title = paste("Top", top_n, "most important variables"),
       x = "Variable importance",
       y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 16),
    legend.position = "none"
  )

print(p3)

# 4. ROC curve - with individual CV curves (fixed version)
create_roc_plot <- function() {
  # Use cross-validation predictions
  cv_predictions <- lasso_model$pred %>%
    filter(lambda == lasso_model$bestTune$lambda)

  # Calculate mean AUC across resamples
  cv_aucs <- cv_predictions %>%
    group_by(Resample) %>%
    summarise(
      auc = as.numeric(auc(roc(obs, hot, levels = c("cold", "hot"), quiet = TRUE))),
      .groups = 'drop'
    )

  mean_auc <- mean(cv_aucs$auc)

  # Create individual ROC curves for each resample with error handling
  roc_curves <- list()
  valid_resamples <- c()

  for (resample in unique(cv_predictions$Resample)) {
    resample_data <- cv_predictions %>% filter(Resample == resample)

    # Check if there are enough observations
    if (nrow(resample_data) > 1 && length(unique(resample_data$obs)) > 1) {
      tryCatch({
        roc_obj <- roc(resample_data$obs, resample_data$hot,
                       levels = c("cold", "hot"), quiet = TRUE)

        if (length(roc_obj$sensitivities) > 0 && length(roc_obj$specificities) > 0) {
          roc_curves[[resample]] <- data.frame(
            Sensitivity = roc_obj$sensitivities,
            Specificity = roc_obj$specificities,
            Resample = resample
          )
          valid_resamples <- c(valid_resamples, resample)
        }
      }, error = function(e) {
        message("Skipping resample ", resample, " due to error: ", e$message)
      })
    }
  }

  # Combine all valid ROC curves
  if (length(roc_curves) > 0) {
    roc_curves_df <- bind_rows(roc_curves)
  } else {
    roc_curves_df <- data.frame(Sensitivity = numeric(), Specificity = numeric(), Resample = character())
  }

  # Create mean ROC curve using all data
  roc_obj <- roc(response = cv_predictions$obs,
                 predictor = cv_predictions$hot,
                 levels = c("cold", "hot"),
                 quiet = TRUE)

  p4 <- ggplot() +
    # Individual CV curves in light gray (if available)
    {if (nrow(roc_curves_df) > 0)
      geom_line(data = roc_curves_df,
                aes(x = 1 - Specificity, y = Sensitivity, group = Resample),
                color = "gray80", alpha = 0.4, linewidth = 0.5)
    } +
    # Mean ROC curve
    geom_line(data = data.frame(
      Sensitivity = roc_obj$sensitivities,
      Specificity = roc_obj$specificities
    ), aes(x = 1 - Specificity, y = Sensitivity),
    color = my_colors[1], linewidth = 1.5, alpha = 0.9) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "gray50", alpha = 0.7) +
    annotate("text", x = 0.6, y = 0.2,
             label = sprintf("Mean AUC = %.4f", mean_auc),
             size = 5, fontface = "bold", color = my_colors[3]) +
    labs(title = "ROC curve",
         x = "1 - Specificity",
         y = "Sensitivity") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      aspect.ratio = 1
    )

  cat("\n=== Cross-validation Performance ===\n")
  cat("Mean AUC:", round(mean_auc, 4), "\n")
  cat("SD:", round(sd(cv_aucs$auc), 4), "\n")
  cat("Range:", round(min(cv_aucs$auc), 4), "-", round(max(cv_aucs$auc), 4), "\n")
  cat("Number of valid ROC curves:", length(valid_resamples), "/", length(unique(cv_predictions$Resample)), "\n")
  cat("Number of selected variables:", length(selected_vars), "\n")

  print(p4)
  return(p4)
}

# Plot ROC curve
p4 <- create_roc_plot()

# 5. Coefficients of selected variables
cat("\n=== Selected Variables (n =", length(selected_vars), ") ===\n")
if(length(selected_vars) > 0) {
  selected_coefs <- final_coef[selected_vars, ]
  coef_summary <- data.frame(
    Variable = selected_vars,
    Coefficient = round(selected_coefs, 4),
    Abs_Coefficient = abs(round(selected_coefs, 4))
  ) %>% arrange(desc(Abs_Coefficient))

  print(coef_summary)

  # Create coefficient bar plot
  p5 <- coef_summary %>%
    mutate(Variable = factor(Variable, levels = Variable[order(Abs_Coefficient)])) %>%
    ggplot(aes(x = Coefficient, y = Variable, fill = Coefficient > 0)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = c(my_colors[2], my_colors[1])) +
    labs(title = "LASSO Model Coefficients",
         subtitle = paste("Selected variables (n =", length(selected_vars), ")"),
         x = "Coefficient Value",
         y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "none",
      panel.grid.major.y = element_blank()
    )

  print(p5)
} else {
  cat("No variables selected (null model)\n")
}

# Save model and plots
save(lasso_model, file = "lasso_model.Rdata")
################################################################################Analysis of R vs NR with lasso####
rm(list=ls())
gc()

# Set the working directory and load packages
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
library(reshape2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(GSEABase)
library(dplyr)

# Read the data and build the dataset list
ici_precog_info <- read.csv("ici_precog_all_patient_clinical_data_FOR_DOWNLOAD.csv")
datasets <- unique(ici_precog_info$Dataset)
datasets <- unique(c(datasets, "Zhao_NatMed_2019", "Kim_NatMed_2018", "GSE67501", "Braun_NatMed_2020",
                     "GSE93157", "Gide_CancerCell_2019", "GSE115821", "GSE106128", "GSE100797",
                     "GSE93157", "Müller_JCI_2018", "GSE100797", "GSE135222", "GSE19423", "GSE100797",
                     "GSE78220", "IMVigor210"))
vector <- datasets[datasets != "Gide_Cell_2019"]

# Check for valid datasets
vector_responder <- c()

for (q in 1:length(vector)) {
  data_file <- paste0(vector[q], ".Rdata")
  info_file <- paste0(vector[q], "_info.Rdata")

  if(file.exists(data_file) && file.exists(info_file)) {
    load(data_file)
    load(info_file)

    info_var_name <- paste0(vector[q], "_info")

    if(exists(info_var_name)) {
      info_data <- get(info_var_name)

      if("responder" %in% colnames(info_data)) {
        responder_values <- info_data$responder
        if(!all(is.na(responder_values)) && length(responder_values) > 0) {
          vector_responder <- c(vector_responder, vector[q])
          cat(paste0("OK dataset ", vector[q], " has a valid responder column\n"))
        } else {
          cat(paste0("WARNING: ", vector[q], " has a responder column that is entirely NA or empty\n"))
        }
      } else {
        cat(paste0("WARNING: ", vector[q], " has no responder column\n"))
      }
    } else {
      cat(paste0("WARNING: ", vector[q], " failed to load its info table\n"))
    }
  } else {
    cat(paste0("WARNING: ", vector[q], " data file not found\n"))
  }
}

cat("\nDatasets with a valid responder column:\n")
print(vector_responder)
cat("Number of valid datasets:", length(vector_responder), "\n")

# Stop if no valid dataset is available
if(length(vector_responder) == 0) {
  stop("No dataset has a valid responder column")
}

# Load the required functions and gene sets
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
geneSets <- getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
t_test_result_all=data.frame()
vector_responder=vector_responder[-c(17,46,55)]
for(i in 1:length(vector_responder)){
  exp=get(vector_responder[i])
  if(length(rownames(exp))<800){
    next
  }
  im_timer <- deconvo_tme(eset = exp, method = "xcell")
  im_timer <- as.data.frame(im_timer)
  rownames(im_timer) <- im_timer[,1]
  im_timer <- im_timer[,-1]
  colnames(im_timer) <- gsub("_xCell$", "", colnames(im_timer))
  im_timer <- as.data.frame(scale(im_timer))
  load("C:/Users/user/Desktop/input/lasso_model.Rdata")
  im_timer$subtype = predict(object = lasso_model,
                              newdata = as.matrix(im_timer))
  table(im_timer$subtype)
  im_timer=im_timer[im_timer$subtype=="cold",]
  exp_responder=exp[,rownames(im_timer)]
  WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                         expression_profile=exp_responder,
                                         foundation="relative_ssGSEA",
                                         activation_geneset=NA,
                                         inhibition_geneset=NA,
                                         geneSets_gmt=geneSets,
                                         min.sz=1,
                                         max.sz=10000,
                                         geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                             rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  WNT_score=WNT_Score$final_activity_score
  pro=get(paste0(vector_responder[i],"_info"))
  WNT_score_pro=merge(pro,WNT_score,by.x="tag",by.y="row.names")
  WNT_score_pro=WNT_score_pro[,c("tag","cancer","responder","activity_score")]
  WNT_score_pro$responder=factor(WNT_score_pro$responder,levels = c("TRUE","FALSE"))
  if(length(unique(WNT_score_pro$responder))!=1){
    t_test_result = t.test(activity_score ~ responder, data = WNT_score_pro)
    t_test_result=data.frame(datasets=vector_responder[i],cancer=unique(WNT_score_pro$cancer),p.value=t_test_result$p.value)
    t_test_result_all=rbind(t_test_result_all,t_test_result)
  }
}
save(t_test_result_all,file="t_test_result_all.Rdata")
################################################################################OS analysis and lasso analysis####
rm(list=ls())
gc()

# Set the working directory and load packages
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
library(reshape2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(GSEABase)
library(dplyr)

# Read the data and build the dataset list
ici_precog_info <- read.csv("ici_precog_all_patient_clinical_data_FOR_DOWNLOAD.csv")
datasets <- unique(ici_precog_info$Dataset)
datasets <- unique(c(datasets, "Zhao_NatMed_2019", "Kim_NatMed_2018", "GSE67501", "Braun_NatMed_2020",
                     "GSE93157", "Gide_CancerCell_2019", "GSE115821", "GSE106128", "GSE100797",
                     "GSE93157", "Müller_JCI_2018", "GSE100797", "GSE135222", "GSE19423", "GSE100797",
                     "GSE78220", "IMVigor210"))
vector <- datasets[datasets != "Gide_Cell_2019"]

# Check for valid survival data (os.time and os.status)
vector_survival <- c()

for (q in 1:length(vector)) {
  data_file <- paste0(vector[q], ".Rdata")
  info_file <- paste0(vector[q], "_info.Rdata")

  if(file.exists(data_file) && file.exists(info_file)) {
    load(data_file)
    load(info_file)

    info_var_name <- paste0(vector[q], "_info")

    if(exists(info_var_name)) {
      info_data <- get(info_var_name)

      # Check that os.time and os.status are both present
      has_os_time <- "os.time" %in% colnames(info_data)
      has_os_status <- "os.status" %in% colnames(info_data)

      if(has_os_time && has_os_status) {
        os_time_values <- info_data$os.time
        os_status_values <- info_data$os.status

        # Check that both columns contain valid data
        if(!all(is.na(os_time_values)) && length(os_time_values) > 0 &&
           !all(is.na(os_status_values)) && length(os_status_values) > 0) {
          vector_survival <- c(vector_survival, vector[q])
          cat(paste0("OK dataset ", vector[q], " has valid os.time and os.status columns\n"))
        } else {
          cat(paste0("WARNING: ", vector[q], " has os.time or os.status entirely NA or empty\n"))
        }
      } else {
        missing_cols <- c()
        if(!has_os_time) missing_cols <- c(missing_cols, "os.time")
        if(!has_os_status) missing_cols <- c(missing_cols, "os.status")
        cat(paste0("WARNING: ", vector[q], " missing columns: ", paste(missing_cols, collapse = ", "), "\n"))
      }
    } else {
      cat(paste0("WARNING: ", vector[q], " failed to load its info table\n"))
    }
  } else {
    cat(paste0("WARNING: ", vector[q], " data file not found\n"))
  }
}

cat("\nDatasets with valid os.time and os.status columns:\n")
print(vector_survival)
cat("Number of valid datasets:", length(vector_survival), "\n")
library(IOBR)
# Stop if no valid dataset is available
if(length(vector_survival) == 0) {
  stop("No dataset has valid os.time and os.status columns")
}

# Load the required functions and gene sets
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
geneSets <- getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")

os_test_result_all=data.frame()
vector_survival=vector_survival[-22]
for(i in 1:length(vector_survival)){
  exp=get(vector_survival[i])
  if(length(rownames(exp))<800){
    next
  }
  im_timer <- deconvo_tme(eset = exp, method = "xcell")
  im_timer <- as.data.frame(im_timer)
  rownames(im_timer) <- im_timer[,1]
  im_timer <- im_timer[,-1]
  colnames(im_timer) <- gsub("_xCell$", "", colnames(im_timer))
  im_timer <- as.data.frame(scale(im_timer))
  load("C:/Users/user/Desktop/input/lasso_model.Rdata")
  im_timer$subtype = predict(object = lasso_model,
                             newdata = as.matrix(im_timer))
  table(im_timer$subtype)
  im_timer=im_timer[im_timer$subtype=="cold",]
  exp_responder=exp[,rownames(im_timer)]
  WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                         expression_profile=exp_responder,
                                         foundation="relative_ssGSEA",
                                         activation_geneset=NA,
                                         inhibition_geneset=NA,
                                         geneSets_gmt=geneSets,
                                         min.sz=1,
                                         max.sz=10000,
                                         geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                             rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  WNT_score=WNT_Score$final_activity_score
  pro=get(paste0(vector_survival[i],"_info"))
  WNT_score_pro=merge(pro,WNT_score,by.x="tag",by.y="row.names")
  WNT_score_pro=WNT_score_pro[,c("tag","cancer","os.time","os.status","activity_score")]
  quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                             labels = c("Low","Medium","High"))
  WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd3 = survdiff(Surv(os.time, os.status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
  p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
  os_test_result=data.frame(p3=p.val3,datasets=vector_survival[i])
  os_test_result_all=rbind(os_test_result_all,os_test_result)
}
save(os_test_result_all,file="os_test_result_all.Rdata")
################################################################################PFS analysis and lasso analysis####
rm(list=ls())
gc()

# Set the working directory and load packages
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
library(reshape2)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(GSEABase)
library(dplyr)

# Read the data and build the dataset list
ici_precog_info <- read.csv("ici_precog_all_patient_clinical_data_FOR_DOWNLOAD.csv")
datasets <- unique(ici_precog_info$Dataset)
datasets <- unique(c(datasets, "Zhao_NatMed_2019", "Kim_NatMed_2018", "GSE67501", "Braun_NatMed_2020",
                     "GSE93157", "Gide_CancerCell_2019", "GSE115821", "GSE106128", "GSE100797",
                     "GSE93157", "Müller_JCI_2018", "GSE100797", "GSE135222", "GSE19423", "GSE100797",
                     "GSE78220", "IMVigor210"))
vector <- datasets[datasets != "Gide_Cell_2019"]

# Check for valid survival data (pfs.time and pfs.status)
vector_survival <- c()

for (q in 1:length(vector)) {
  data_file <- paste0(vector[q], ".Rdata")
  info_file <- paste0(vector[q], "_info.Rdata")

  if(file.exists(data_file) && file.exists(info_file)) {
    load(data_file)
    load(info_file)

    info_var_name <- paste0(vector[q], "_info")

    if(exists(info_var_name)) {
      info_data <- get(info_var_name)

      # Check that pfs.time and pfs.status are both present
      has_pfs_time <- "pfs.time" %in% colnames(info_data)
      has_pfs_status <- "pfs.status" %in% colnames(info_data)

      if(has_pfs_time && has_pfs_status) {
        pfs_time_values <- info_data$pfs.time
        pfs_status_values <- info_data$pfs.status

        # Check that both columns contain valid data
        if(!all(is.na(pfs_time_values)) && length(pfs_time_values) > 0 &&
           !all(is.na(pfs_status_values)) && length(pfs_status_values) > 0) {
          vector_survival <- c(vector_survival, vector[q])
          cat(paste0("OK dataset ", vector[q], " has valid pfs.time and pfs.status columns\n"))
        } else {
          cat(paste0("WARNING: ", vector[q], " has pfs.time or pfs.status entirely NA or empty\n"))
        }
      } else {
        missing_cols <- c()
        if(!has_pfs_time) missing_cols <- c(missing_cols, "pfs.time")
        if(!has_pfs_status) missing_cols <- c(missing_cols, "pfs.status")
        cat(paste0("WARNING: ", vector[q], " missing columns: ", paste(missing_cols, collapse = ", "), "\n"))
      }
    } else {
      cat(paste0("WARNING: ", vector[q], " failed to load its info table\n"))
    }
  } else {
    cat(paste0("WARNING: ", vector[q], " data file not found\n"))
  }
}

cat("\nDatasets with valid pfs.time and pfs.status columns:\n")
print(vector_survival)
cat("Number of valid datasets:", length(vector_survival), "\n")
library(IOBR)
# Stop if no valid dataset is available
if(length(vector_survival) == 0) {
  stop("No dataset has valid os.time and os.status columns")
}

# Load the required functions and gene sets
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
geneSets <- getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")

pfs_test_result_all=data.frame()
vector_survival=vector_survival[-19]
for(i in 1:length(vector_survival)){
  exp=get(vector_survival[i])
  if(length(rownames(exp))<800){
    next
  }
  im_timer <- deconvo_tme(eset = exp, method = "xcell")
  im_timer <- as.data.frame(im_timer)
  rownames(im_timer) <- im_timer[,1]
  im_timer <- im_timer[,-1]
  colnames(im_timer) <- gsub("_xCell$", "", colnames(im_timer))
  im_timer <- as.data.frame(scale(im_timer))
  load("C:/Users/user/Desktop/input/lasso_model.Rdata")
  im_timer$subtype = predict(object = lasso_model,
                             newdata = as.matrix(im_timer))
  table(im_timer$subtype)
  im_timer=im_timer[im_timer$subtype=="cold",]
  exp_responder=exp[,rownames(im_timer)]
  WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                         expression_profile=exp_responder,
                                         foundation="relative_ssGSEA",
                                         activation_geneset=NA,
                                         inhibition_geneset=NA,
                                         geneSets_gmt=geneSets,
                                         min.sz=1,
                                         max.sz=10000,
                                         geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                             rep("neg",sum(grepl("WPIGS", names(geneSets))))))
  WNT_score=WNT_Score$final_activity_score
  pro=get(paste0(vector_survival[i],"_info"))
  WNT_score_pro=merge(pro,WNT_score,by.x="tag",by.y="row.names")
  WNT_score_pro=WNT_score_pro[,c("tag","cancer","pfs.time","pfs.status","activity_score")]
  quartiles3=quantile(WNT_score_pro$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro$level3 = cut(WNT_score_pro$activity_score, breaks = quartiles3, include.lowest = TRUE,
                             labels = c("Low","Medium","High"))
  WNT_score_pro$level3=factor(WNT_score_pro$level3,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd3 = survdiff(Surv(pfs.time, pfs.status) ~ level3,data= WNT_score_pro,na.action = na.exclude)
  p.val3 = 1 - pchisq(fitd3$chisq, length(fitd3$n) - 1)
  pfs_test_result=data.frame(p3=p.val3,datasets=vector_survival[i])
  pfs_test_result_all=rbind(pfs_test_result_all,pfs_test_result)
}
save(pfs_test_result_all,file="pfs_test_result_all.Rdata")
################################################################################Integrated prognostic analysis across datasets (grouped)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015","GSE106128","Braun_NatMed_2020","GSE253260","Zhao_NatMed_2019","Tsimberidou_JExpMed_2024")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}
library(IOBR)
# Dai_Nature_2015=Dai_2025
# Dai_Nature_2015_info=Dai_2025_info
# save(Dai_Nature_2015,file="Dai_Nature_2015.Rdata")
# save(Dai_Nature_2015_info,file="Dai_Nature_2015_info.Rdata")
# Tsimberidou_JExpMed_2024=AMADEUS_2024
# Tsimberidou_JExpMed_2024_info=AMADEUS_2024_info
# save(Tsimberidou_JExpMed_2024,file="Tsimberidou_JExpMed_2024.Rdata")
# save(Tsimberidou_JExpMed_2024_info,file="Tsimberidou_JExpMed_2024_info.Rdata")
for(i in 1:length(vector)){
  exp=get(vector[i])
  im_timer = deconvo_tme(eset = exp, method = "xcell")
  im_timer = as.data.frame(im_timer)
  rownames(im_timer) = im_timer[,1]
  im_timer = im_timer[,-1]
  colnames(im_timer) = gsub("_xCell$", "", colnames(im_timer))
  im_timer = as.data.frame(scale(im_timer))
  load("C:/Users/user/Desktop/input/lasso_model.Rdata")
  im_timer$subtype = predict(object = lasso_model,
                             newdata = as.matrix(im_timer))
  table(im_timer$subtype)
  group=data.frame(Tag=rownames(im_timer),group=im_timer$subtype)
  assign(paste0(vector[i],"_G"),group)
  save(list=paste0(vector[i],"_G"),file=paste0(vector[i],"_G.Rdata"))
}
################################################################################Integrated prognostic analysis across datasets (differential analysis)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015","GSE106128","Braun_NatMed_2020","GSE253260","Zhao_NatMed_2019","Tsimberidou_JExpMed_2024","IMvigor")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
WNT_score_group_all=data.frame()
for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  print(table(group$group))
}

for(i in 1:length(vector)){
  exp=get(vector[i])
  library(GSEABase)
  geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
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
  WNT_score=WNT_Score$final_activity_score
  group=get(paste0(vector[i],"_G"))
  WNT_score_group=merge(WNT_score,group,by.x="row.names",by.y="Tag")
  WNT_score_group$datasets=vector[i]
  WNT_score_group_all=rbind(WNT_score_group_all,WNT_score_group)
}
library(ggplot2)
library(ggpubr)
p <- ggplot(WNT_score_group_all, aes(fill=group, y=activity_score, x=datasets)) +
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
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black", size = 16), # x-axis text
        axis.text.y = element_text(colour = "black", size = 16), # y-axis text
        axis.title.x = element_text(colour = "black", size = 16), # x-axis title
        axis.title.y = element_text(colour = "black", size = 16), # y-axis title
        legend.text = element_text(colour = "black", size = 16), # Legend text
        legend.title = element_text(colour = "black", size = 16), # Legend title
        plot.title = element_text(colour = "black", size = 16)) + # Figure title
  stat_summary(fun.data = 'mean_se',
               geom = "errorbar",
               colour = "black",
               width = 0.3,
               position=position_dodge(0.7)) +
  scale_fill_manual(values = c("#f58f98","#90d7ec")) +
  theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
  stat_compare_means(aes(group = group, label = ..p.signif..), method="t.test") +
  ggtitle("Training datasets")
p
################################################################################Single dataset (Dai_Nature_2015)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Dai_Nature_2015_G=Dai_Nature_2015_G[Dai_Nature_2015_G$group=="cold",]
Dai_Nature_2015=Dai_Nature_2015[,colnames(Dai_Nature_2015)%in%Dai_Nature_2015_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Dai_Nature_2015,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Dai_Nature_2015_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
wnt_IMvigor_prognosis$os.time=wnt_IMvigor_prognosis$os.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 72),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Dai_Nature_2015",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Dai_Nature_2015_G
group=group[group$group=="cold",]
exp=Dai_Nature_2015
prognosis=Dai_Nature_2015_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                             labels = c("Low","Medium","High"))
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("Dai_Nature_2015(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Dai_Nature_2015_os_p=p_all
save(Dai_Nature_2015_os_p,file="Dai_Nature_2015_os_p.Rdata")
################################################################################Single dataset (GSE106128)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE106128")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
GSE106128_G=GSE106128_G[GSE106128_G$group=="cold",]
GSE106128=GSE106128[,colnames(GSE106128)%in%GSE106128_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=GSE106128,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,GSE106128_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
wnt_IMvigor_prognosis$os.time=wnt_IMvigor_prognosis$os.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 70),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="GSE106128",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE106128")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=GSE106128_G
group=group[group$group=="cold",]
exp=GSE106128
prognosis=GSE106128_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                             labels = c("Low","Medium","High"))
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("GSE106128(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
GSE106128_os_p=p_all
save(GSE106128_os_p,file="GSE106128_os_p.Rdata")
################################################################################Single dataset (Braun_NatMed_2020)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Braun_NatMed_2020_G=Braun_NatMed_2020_G[Braun_NatMed_2020_G$group=="cold",]
Braun_NatMed_2020=Braun_NatMed_2020[,colnames(Braun_NatMed_2020)%in%Braun_NatMed_2020_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Braun_NatMed_2020,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Braun_NatMed_2020_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
wnt_IMvigor_prognosis$os.time=wnt_IMvigor_prognosis$os.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 74),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Braun_NatMed_2020",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Braun_NatMed_2020_G
group=group[group$group=="cold",]
exp=Braun_NatMed_2020
prognosis=Braun_NatMed_2020_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                             labels = c("Low","Medium","High"))
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("Braun_NatMed_2020(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Braun_NatMed_2020_os_p=p_all
save(Braun_NatMed_2020_os_p,file="Braun_NatMed_2020_os_p.Rdata")
################################################################################Single dataset (GSE253260)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
GSE253260_G=GSE253260_G[GSE253260_G$group=="cold",]
GSE253260=GSE253260[,colnames(GSE253260)%in%GSE253260_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=GSE253260,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,GSE253260_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
wnt_IMvigor_prognosis$os.time=wnt_IMvigor_prognosis$os.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 74),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="GSE253260",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=GSE253260_G
group=group[group$group=="cold",]
exp=GSE253260
prognosis=GSE253260_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  if(length(unique(quartiles)) < 4){
    n <- nrow(data_prognosis)
    ranks <- rank(data_prognosis$gene, ties.method = "first")
    group_size <- ceiling(n / 3)
    data_prognosis$level <- "Medium"
    data_prognosis$level[ranks <= group_size] <- "Low"
    data_prognosis$level[ranks > 2 * group_size] <- "High"
    data_prognosis$level <- factor(data_prognosis$level, levels = c("Low","Medium","High"))
  }else{
    data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  }
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("GSE253260(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
GSE253260_os_p=p_all
save(GSE253260_os_p,file="GSE253260_os_p.Rdata")
################################################################################Single dataset (Zhao_NatMed_2019)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Zhao_NatMed_2019")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Zhao_NatMed_2019_G=Zhao_NatMed_2019_G[Zhao_NatMed_2019_G$group=="cold",]
Zhao_NatMed_2019=Zhao_NatMed_2019[,colnames(Zhao_NatMed_2019)%in%Zhao_NatMed_2019_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Zhao_NatMed_2019,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Zhao_NatMed_2019_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
wnt_IMvigor_prognosis$os.time=wnt_IMvigor_prognosis$os.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 57),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Zhao_NatMed_2019",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Zhao_NatMed_2019")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Zhao_NatMed_2019_G
group=group[group$group=="cold",]
exp=Zhao_NatMed_2019
prognosis=Zhao_NatMed_2019_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  if(length(unique(quartiles)) < 4){
    n <- nrow(data_prognosis)
    ranks <- rank(data_prognosis$gene, ties.method = "first")
    group_size <- ceiling(n / 3)
    data_prognosis$level <- "Medium"
    data_prognosis$level[ranks <= group_size] <- "Low"
    data_prognosis$level[ranks > 2 * group_size] <- "High"
    data_prognosis$level <- factor(data_prognosis$level, levels = c("Low","Medium","High"))
  }else{
    data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  }
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("Zhao_NatMed_2019(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Zhao_NatMed_2019_os_p=p_all
save(Zhao_NatMed_2019_os_p,file="Zhao_NatMed_2019_os_p.Rdata")
################################################################################Single dataset (Tsimberidou_JExpMed_2024)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Tsimberidou_JExpMed_2024")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Tsimberidou_JExpMed_2024_G=Tsimberidou_JExpMed_2024_G[Tsimberidou_JExpMed_2024_G$group=="cold",]
Tsimberidou_JExpMed_2024=Tsimberidou_JExpMed_2024[,colnames(Tsimberidou_JExpMed_2024)%in%Tsimberidou_JExpMed_2024_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Tsimberidou_JExpMed_2024,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Tsimberidou_JExpMed_2024_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
wnt_IMvigor_prognosis$os.time=wnt_IMvigor_prognosis$os.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 44),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Tsimberidou_JExpMed_2024",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Tsimberidou_JExpMed_2024")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Tsimberidou_JExpMed_2024_G
group=group[group$group=="cold",]
exp=Tsimberidou_JExpMed_2024
prognosis=Tsimberidou_JExpMed_2024_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  if(length(unique(quartiles)) < 4){
    n <- nrow(data_prognosis)
    ranks <- rank(data_prognosis$gene, ties.method = "first")
    group_size <- ceiling(n / 3)
    data_prognosis$level <- "Medium"
    data_prognosis$level[ranks <= group_size] <- "Low"
    data_prognosis$level[ranks > 2 * group_size] <- "High"
    data_prognosis$level <- factor(data_prognosis$level, levels = c("Low","Medium","High"))
  }else{
    data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  }
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("Tsimberidou_JExpMed_2024(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Tsimberidou_JExpMed_2024_os_p=p_all
save(Tsimberidou_JExpMed_2024_os_p,file="Tsimberidou_JExpMed_2024_os_p.Rdata")
################################################################################Single dataset (IMvigor)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("IMvigor")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
IMvigor_G=IMvigor_G[IMvigor_G$group=="cold",]
IMvigor=IMvigor[,colnames(IMvigor)%in%IMvigor_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=IMvigor,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,IMvigor_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","os.time","os.status")]
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(os.time, os.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(os.time, os.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 44),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="IMvigor",
  subtitle = paste0("OS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("IMvigor")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=IMvigor_G
group=group[group$group=="cold",]
exp=IMvigor
prognosis=IMvigor_info[,c("tag","os.time","os.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  if(length(unique(quartiles)) < 4){
    n <- nrow(data_prognosis)
    ranks <- rank(data_prognosis$gene, ties.method = "first")
    group_size <- ceiling(n / 3)
    data_prognosis$level <- "Medium"
    data_prognosis$level[ranks <= group_size] <- "Low"
    data_prognosis$level[ranks > 2 * group_size] <- "High"
    data_prognosis$level <- factor(data_prognosis$level, levels = c("Low","Medium","High"))
  }else{
    data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  }
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("IMvigor(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
IMvigor_os_p=p_all
save(IMvigor_os_p,file="IMvigor_os_p.Rdata")
################################################################################Single dataset (Dai_Nature_2015)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Dai_Nature_2015_G=Dai_Nature_2015_G[Dai_Nature_2015_G$group=="cold",]
Dai_Nature_2015=Dai_Nature_2015[,colnames(Dai_Nature_2015)%in%Dai_Nature_2015_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Dai_Nature_2015,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Dai_Nature_2015_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","pfs.time","pfs.status")]
wnt_IMvigor_prognosis$pfs.time=wnt_IMvigor_prognosis$pfs.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(pfs.time, pfs.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(pfs.time, pfs.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(pfs.time, pfs.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 72),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Dai_Nature_2015",
  subtitle = paste0("PFS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Dai_Nature_2015_G
group=group[group$group=="cold",]
exp=Dai_Nature_2015
prognosis=Dai_Nature_2015_info[,c("tag","pfs.time","pfs.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                             labels = c("Low","Medium","High"))
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("Dai_Nature_2015(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Dai_Nature_2015_pfs_p=p_all
save(Dai_Nature_2015_pfs_p,file="Dai_Nature_2015_pfs_p.Rdata")
################################################################################Single dataset (GSE253260)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
GSE253260_G=GSE253260_G[GSE253260_G$group=="cold",]
GSE253260=GSE253260[,colnames(GSE253260)%in%GSE253260_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=GSE253260,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,GSE253260_info,by.x="row.names",by.y="tag")
wnt_IMvigor_prognosis=wnt_IMvigor[,c("activity_score","pfs.time","pfs.status")]
wnt_IMvigor_prognosis$pfs.time=wnt_IMvigor_prognosis$pfs.time/30
quartiles <- quantile(wnt_IMvigor_prognosis$activity_score, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
wnt_IMvigor_prognosis$level = cut(wnt_IMvigor_prognosis$activity_score, breaks = quartiles, include.lowest = TRUE,
                                  labels = c("Low","Medium","High"))
wnt_IMvigor_prognosis$level=factor(wnt_IMvigor_prognosis$level,levels=c("Low","Medium","High"))
library(survival)
library(survminer)
fitd = survdiff(Surv(pfs.time, pfs.status) ~ level,data= wnt_IMvigor_prognosis,na.action = na.exclude)
p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit = survfit(Surv(pfs.time, pfs.status)~ level,data= wnt_IMvigor_prognosis,type= "kaplan-meier",error= "greenwood",conf.type = "plain",na.action = na.exclude)
ps = pairwise_survdiff(Surv(pfs.time, pfs.status)~ level,data= wnt_IMvigor_prognosis,p.adjust.method = "none")
mycol=c('#009ad6',"grey50",'#d71345')
names(fit$strata)=gsub("group=", "", names(fit$strata))
library(RColorBrewer)
library(tibble)
library(ggpp)
p = ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = F,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 69),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="GSE253260",
  subtitle = paste0("PFS","(Immune cold,n=",nrow(wnt_IMvigor_prognosis),")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")
print(p)
rm(list = ls())
gc()
library(GSVA)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=GSE253260_G
group=group[group$group=="cold",]
exp=GSE253260
prognosis=GSE253260_info[,c("tag","pfs.time","pfs.status")]
colnames(prognosis)=c("id","time","status")
data=merge(as.data.frame(t(exp)),prognosis,by.x="row.names",by.y="id")
gene=c("AXIN2","ZNRF3","RNF43","LGR5","TCF7","PDCD1","CD274","TGFB1")
gene_p=c()
data=data[data$Row.names%in%group$Tag,]
for(i in 1:length(gene)){
  data_prognosis=data[,c(gene[i],"time","status")]
  colnames(data_prognosis)=c("gene","time","status")
  data_prognosis=na.omit(data_prognosis)
  quartiles <- quantile(data_prognosis$gene, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  if(length(unique(quartiles)) < 4){
    n <- nrow(data_prognosis)
    ranks <- rank(data_prognosis$gene, ties.method = "first")
    group_size <- ceiling(n / 3)
    data_prognosis$level <- "Medium"
    data_prognosis$level[ranks <= group_size] <- "Low"
    data_prognosis$level[ranks > 2 * group_size] <- "High"
    data_prognosis$level <- factor(data_prognosis$level, levels = c("Low","Medium","High"))
  }else{
    data_prognosis$level = cut(data_prognosis$gene, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  }
  data_prognosis$level=factor(data_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= data_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  gene_p=c(gene_p,p.val)
}
gene_p=as.data.frame(gene_p)
gene_p$signature=gene
colnames(gene_p)=c("value","ID")
exp_cal=exp[,group$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=exp_cal,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
WNT_score_pro=merge(prognosis,WNT_score,by.x="id",by.y="row.names")
model=c("activity_score","pos_score","neg_score")
model_p=c()
for(i in 1:length(model)){
  WNT_score_pro_prognosis=WNT_score_pro[,c(model[i],"time","status")]
  colnames(WNT_score_pro_prognosis)=c("model","time","status")
  WNT_score_pro_prognosis=na.omit(WNT_score_pro_prognosis)
  quartiles <- quantile(WNT_score_pro_prognosis$model, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  WNT_score_pro_prognosis$level = cut(WNT_score_pro_prognosis$model, breaks = quartiles, include.lowest = TRUE,
                                      labels = c("Low","Medium","High"))
  WNT_score_pro_prognosis$level=factor(WNT_score_pro_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= WNT_score_pro_prognosis,na.action = na.exclude)
  p.val = 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  p.val
  model_p=c(model_p,p.val)
}
model_p=as.data.frame(model_p)
model_p$signature=model
model_p$signature=c("Wnt/β-catenin pathway activity score",
                    "Wnt/β-catenin pathway activation score",
                    "Wnt/β-catenin pathway inhibition score")
colnames(model_p)=c("value","ID")
library(IOBR)
library(GSEABase)
immunotherapy=getGmt("immunotherapy.gmt")
immunotherapy_gmt=lapply(immunotherapy, function(x){ x@geneIds })
names(immunotherapy_gmt)=names(immunotherapy)
im_ssGSEA=calculate_sig_score(eset = as.data.frame(exp_cal),signature = immunotherapy_gmt,method="ssgsea",mini_gene_count=1)
im_ssGSEA=as.data.frame(im_ssGSEA)
rownames(im_ssGSEA)=im_ssGSEA[,1]
im_ssGSEA=im_ssGSEA[,-1]
im_ssGSEA_prognosis=merge(im_ssGSEA,prognosis,by.x="row.names",by.y="id")
signature=c(
  "CD8+_T_effector",
  "Pan_F_TBRS",
  "Antigen_processing_machinery",
  "Immune_checkpoint",
  "T_cell_inflamed_GEP",
  "IFNG_signature",
  "T_cell_accumulation",
  "T_cell_exhaustion",
  "T_cell_regulatory",
  "ICB_resistance",
  "MDSC",
  "TAM",
  "CAF"
)
signature_p=c()
for(i in 1:length(signature)){
  ssGSEA_prognosis=im_ssGSEA_prognosis[,c(signature[i],"time","status")]
  colnames(ssGSEA_prognosis)=c("signature","time","status")
  ssGSEA_prognosis=na.omit(ssGSEA_prognosis)
  quartiles <- quantile(ssGSEA_prognosis$signature, probs = seq(0, 1, by = 1/3), na.rm = TRUE)
  ssGSEA_prognosis$level = cut(ssGSEA_prognosis$signature, breaks = quartiles, include.lowest = TRUE,
                               labels = c("Low","Medium","High"))
  ssGSEA_prognosis$level=factor(ssGSEA_prognosis$level,levels=c("Low","Medium","High"))
  library(survival)
  library(survminer)
  fitd = survdiff(Surv(time, status) ~ level,data= ssGSEA_prognosis,na.action = na.exclude)
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
y_intercept_value <- -log10(0.05)
p_all$ID = gsub("_", " ", p_all$ID)
p_all$ID = ifelse(p_all$ID=="Pan F TBRS","Pan-F-TBRS",p_all$ID)
ggplot(p_all, aes(x = reorder(ID, -value), y = value, fill = value)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(title = "",
       subtitle = "",
       x = "Immunotherapy indicators",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 16) +  # Set base size to 16
  theme(
    plot.title = element_text(hjust = 0, size = 16, color = "black"), # Set title size to 16, colour black
    axis.text = element_text(size = 16, color = "black"),  # Set axis text size to 16, colour black
    axis.title = element_text(size = 16, color = "black"),  # Set axis title size to 16, colour black
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 80, hjust = 1, color = "black")  # Set x-axis text size to 16, colour black
  ) +
  scale_fill_gradient(low = "#90d7ec", high = "#f58f98") +
  ggtitle("GSE253260(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
GSE253260_pfs_p=p_all
save(GSE253260_pfs_p,file="GSE253260_pfs_p.Rdata")
################################################################################Heatmap visualisation of everything####
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/input")
vector <- c("Dai_Nature_2015_os", "GSE106128_os", "Braun_NatMed_2020_os", "GSE253260_os", "Zhao_NatMed_2019_os","Tsimberidou_JExpMed_2024_os",
            "Dai_Nature_2015_pfs","GSE253260_pfs","IMvigor_os")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/", vector[l], "_p.Rdata"))
  exp=get(paste0(vector[l], "_p"))
  colnames(exp)=c(vector[l],"ID")
  assign(paste0(vector[l], "_p"),exp)
}
exp=merge(Dai_Nature_2015_os_p,GSE106128_os_p,by="ID")
exp=merge(exp,Braun_NatMed_2020_os_p,by="ID")
exp=merge(exp,GSE253260_os_p,by="ID")
exp=merge(exp,Zhao_NatMed_2019_os_p,by="ID")
exp=merge(exp,Tsimberidou_JExpMed_2024_os_p,by="ID")
exp=merge(exp,Dai_Nature_2015_pfs_p,by="ID")
exp=merge(exp,GSE253260_pfs_p,by="ID")
exp=merge(exp,IMvigor_os_p,by="ID")
rownames(exp)=exp[,1]
exp=exp[,-1]
exp$sum=rowMeans(exp)

# Sort by sum
exp_sorted <- exp[order(exp$sum, decreasing = TRUE), ]
annotation_data <- data.frame(Sum = exp_sorted$sum)
rownames(annotation_data) <- rownames(exp_sorted)
heatmap_data <- exp_sorted[, -which(colnames(exp_sorted) == "sum")]

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tibble)

# 1. Create the feature-name column and sort it
heatmap_df <- as.data.frame(heatmap_data) %>%
  tibble::rownames_to_column("feature") %>%
  arrange(desc(exp_sorted$sum)) %>%  # Sort by the sum
  mutate(feature = factor(feature, levels = feature))  # Keep the sort order

# 2. Convert to long format for the heatmap
heatmap_long <- heatmap_df %>%
  pivot_longer(-feature, names_to = "dataset", values_to = "value")

# 3. Build the data for the significance asterisks
significance_df <- heatmap_long %>%
  mutate(star = ifelse(value > -log10(0.001), "***",ifelse(value > -log10(0.01), "**",ifelse(value > -log10(0.05), "*",""))))

# 4. Draw the heatmap
heatmap_plot <- ggplot(heatmap_long, aes(x = dataset, y = feature, fill = value)) +
  geom_tile(color = "#a1a3a6", size = 0.5) +
  geom_text(data = significance_df, aes(label = star),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",
                       midpoint = median(heatmap_long$value)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "", fill = expression(-log[10](italic(p) ))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "bottom"
  )

# 5. Prepare the bar chart data
bar_data <- data.frame(
  feature = heatmap_df$feature,
  importance = exp_sorted$sum[match(heatmap_df$feature, rownames(exp_sorted))]
)

# 6. Draw the bar chart
barplot_plot <- ggplot(bar_data, aes(x = feature, y = importance)) +
  geom_col(fill = "#ffce7b", width = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "Mean") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.position = "none"
  )

# 7. Combine the figures
combined_plot <- heatmap_plot + barplot_plot +
  plot_layout(widths = c(3, 1))  # Adjust the width ratio

print(combined_plot)
################################################################################Processing all cohorts in a loop####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("C:/Users/user/Desktop/input/IMvigor.Rdata")
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
library(limma)
limma.one.sided = function(fit, lower = FALSE){
  se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total=fit$df.prior + fit$df.residual
  pt(fit$t, df = df.total, lower.tail = lower)
}
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(hot-cold, levels=design)
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
exp_im=as.data.frame(t(xcell_pheno[,2:length(xcell_pheno)]))
group_im=data.frame(Tag=rownames(xcell_pheno),group=xcell_pheno$phenotype)
resultup=Get_limma(gene=exp_im,
                   group=group_im,
                   alternative_limma="greater")
resultup=resultup[order(resultup$P.Value), ]
resultup=resultup[resultup$P.Value<0.05,]
resultdn=Get_limma(gene=exp_im,
                   group=group_im,
                   alternative_limma="less")
resultdn=resultdn[order(resultdn$P.Value), ]
resultdn=resultdn[resultdn$P.Value<0.05,]
save(resultup,file="resultup.Rdata")
save(resultdn,file="resultdn.Rdata")
exp_im=exp_im[c(resultup$gene_id,resultdn$gene_id),]
exp_im=as.data.frame(t(exp_im))
load("C:/Users/user/Desktop/input/hot_cancer_sig.Rdata")
load("C:/Users/user/Desktop/input/cold_cancer_sig.Rdata")
library(GSVA)
scores=ssgseaParam(exprData=as.matrix(IMvigor),geneSets=list("Hot_cancer_signature" = hot_cancer_sig$gene_id,
                                                                 "Cold_cancer_signature" = cold_cancer_sig$gene_id),minSize=1,maxSize=10000)
scores=gsva(scores)
scores=as.data.frame(scale(t(scores)))
exp_im=merge(exp_im,scores,by="row.names")
rownames(exp_im)=exp_im[,1]
exp_im=exp_im[,-1]
exp_im_pheno <- merge(xcell_pheno[, "phenotype", drop = FALSE], exp_im, by = "row.names")
rownames(exp_im_pheno) <- exp_im_pheno[, 1]
exp_im_pheno <- exp_im_pheno[, -1]
library(pheatmap)
library(RColorBrewer)
phenotype_annotation <- data.frame(Phenotype = exp_im_pheno$phenotype)
rownames(phenotype_annotation) <- rownames(exp_im_pheno)
expression_matrix <- as.matrix(exp_im_pheno[, -which(colnames(exp_im_pheno) == "phenotype")])
phenotype_order <- order(exp_im_pheno$phenotype)
expression_matrix_ordered <- expression_matrix[phenotype_order, ]
phenotype_annotation_ordered <- phenotype_annotation[phenotype_order, , drop = FALSE]
signature_features <- c("Hot_cancer_signature", "Cold_cancer_signature")
celltype_features <- setdiff(colnames(expression_matrix_ordered), signature_features)
up_cells <- intersect(celltype_features, resultup$gene_id)
down_cells <- intersect(celltype_features, resultdn$gene_id)
ordered_features <- c("Hot_cancer_signature", up_cells, "Cold_cancer_signature", down_cells)
expression_matrix_ordered_final <- expression_matrix_ordered[, ordered_features]
feature_annotation <- data.frame(
  row.names = ordered_features,
  FeatureType = ifelse(ordered_features %in% signature_features, "Signature", "CellType"),
  Regulation = ifelse(ordered_features %in% c("Hot_cancer_signature", up_cells), "Up", "Down")
)
phenotype_colors <- list(
  Phenotype = c(cold = "#90d7ec", hot = "#f58f98"),
  Regulation = c(Up = "#e74c3c", Down = "#3498db"),
  FeatureType = c(Signature = "#f39c12", CellType = "#9b59b6")
)
data_range <- quantile(expression_matrix_ordered_final, probs = c(0.05, 0.95), na.rm = TRUE)
color_limit <- max(abs(data_range))
breaks <- seq(-color_limit, color_limit, length.out = 100)
gaps_row_positions <- c(
  1 + length(up_cells)
)
pheatmap(t(expression_matrix_ordered_final),
         annotation_col = phenotype_annotation_ordered,
         annotation_row = feature_annotation,
         annotation_colors = phenotype_colors,
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         breaks = breaks,
         main = "Immune Cell Infiltration and Cancer Signatures\nby Phenotype",
         border_color = NA,
         gaps_col = sum(phenotype_annotation_ordered$Phenotype == "cold"),
         gaps_row = gaps_row_positions)
################################################################################Pan-cancer loop heatmap####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)

for(i in 1:length(fs)){
  pancancer_group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  exp=pancancer_exp[,rownames(pancancer_group_fs)]
  im_xcell=deconvo_tme(eset = exp, method = "xcell")
  im_xcell=as.data.frame(im_xcell)
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  colnames(im_xcell)=gsub("_xCell$", "", colnames(im_xcell))
  im_xcell=as.data.frame(scale(im_xcell))
  load("C:/Users/user/Desktop/input/hot_cancer_sig.Rdata")
  load("C:/Users/user/Desktop/input/cold_cancer_sig.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=list("Hot_cancer_signature" = hot_cancer_sig$gene_id,
                                                           "Cold_cancer_signature" = cold_cancer_sig$gene_id),minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_xcell=merge(im_xcell,scores,by="row.names")
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  load("resultup.Rdata")
  load("resultdn.Rdata")
  im_xcell=im_xcell[,c(resultup$gene_id,resultdn$gene_id,"Hot_cancer_signature","Cold_cancer_signature")]

  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)

  expression_matrix <- as.matrix(im_xcell)

  # Order the features to match exp_im
  signature_features <- c("Hot_cancer_signature", "Cold_cancer_signature")
  celltype_features <- setdiff(colnames(expression_matrix), signature_features)

  # Order the cell types by regulatory direction (consistent with exp_im)
  up_cells <- intersect(celltype_features, resultup$gene_id)
  down_cells <- intersect(celltype_features, resultdn$gene_id)

  # Final order: up-regulated features first (Hot signature + up-regulated cell types), then down-regulated features (Cold signature + down-regulated cell types)
  ordered_features <- c("Hot_cancer_signature", up_cells, "Cold_cancer_signature", down_cells)
  expression_matrix_ordered <- expression_matrix[, ordered_features]

  # Build the feature annotation (consistent with exp_im)
  feature_annotation <- data.frame(
    row.names = ordered_features,
    FeatureType = ifelse(ordered_features %in% signature_features, "Signature", "CellType"),
    Regulation = ifelse(ordered_features %in% c("Hot_cancer_signature", up_cells), "Up", "Down")
  )

  # Set the colours (consistent with exp_im)
  feature_colors <- list(
    FeatureType = c(Signature = "#f39c12", CellType = "#9b59b6"),
    Regulation = c(Up = "#e74c3c", Down = "#3498db")
  )

  # Compute the colour range
  data_range <- quantile(expression_matrix_ordered, probs = c(0.05, 0.95), na.rm = TRUE)
  color_limit <- max(abs(data_range))
  breaks <- seq(-color_limit, color_limit, length.out = 100)

  # Compute the gap position (between the Hot and Cold signatures)
  gaps_row_positions <- c(
    1 + length(up_cells)  # Between the Hot signature and the up-regulated cells, before the Cold signature
  )

  # Create the heatmap (same style as exp_im)
  p <- pheatmap(t(expression_matrix_ordered),
                annotation_row = feature_annotation,
                annotation_colors = feature_colors,
                scale = "row",
                cluster_rows = FALSE,    # Do not cluster the rows; use the manual order
                cluster_cols = TRUE,     # Cluster the samples
                show_rownames = TRUE,
                show_colnames = FALSE,
                fontsize_row = 8,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                breaks = breaks,
                main = paste("Immune Cell Infiltration and Cancer Signatures\n", fs[i]),
                border_color = NA,
                gaps_row = gaps_row_positions,
                silent = TRUE)
  output_dir <- "C:/Users/user/Desktop/input/heatmaps"
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  filename <- paste0(output_dir, "/", fs[i], "_consistent_immune_heatmap.png")
  ggsave(filename, p$gtable, width = 12, height = 8, dpi = 300,bg = "white")
}
################################################################################Pan-cancer loop heatmap (with predicted labels)####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)

for(i in 1:length(fs)){
  pancancer_group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  exp=pancancer_exp[,rownames(pancancer_group_fs)]
  im_xcell=deconvo_tme(eset = exp, method = "xcell")
  im_xcell=as.data.frame(im_xcell)
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  colnames(im_xcell)=gsub("_xCell$", "", colnames(im_xcell))
  im_xcell=as.data.frame(scale(im_xcell))
  load("C:/Users/user/Desktop/input/hot_cancer_sig.Rdata")
  load("C:/Users/user/Desktop/input/cold_cancer_sig.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=list("Hot_signature" = hot_cancer_sig$gene_id,
                                                           "Cold_signature" = cold_cancer_sig$gene_id),minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_xcell=merge(im_xcell,scores,by="row.names")
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  load("resultup.Rdata")
  load("resultdn.Rdata")
  im_xcell=im_xcell[,c(resultup$gene_id,resultdn$gene_id,"Hot_signature","Cold_signature")]

  load("lasso_model.Rdata")
  im_xcell_phenotype=data.frame(Tag=rownames(im_xcell),subtype=predict(object = lasso_model,
                                                                       newdata = as.matrix(im_xcell)))
  rownames(im_xcell_phenotype) <- im_xcell_phenotype$Tag

  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)

  expression_matrix <- as.matrix(im_xcell)
  signature_features <- c("Hot_signature", "Cold_signature")
  celltype_features <- setdiff(colnames(expression_matrix), signature_features)
  up_cells <- intersect(celltype_features, resultup$gene_id)
  down_cells <- intersect(celltype_features, resultdn$gene_id)
  ordered_features <- c("Hot_signature", up_cells, "Cold_signature", down_cells)
  expression_matrix_ordered <- expression_matrix[, ordered_features]
  sample_annotation <- data.frame(
    row.names = rownames(expression_matrix_ordered),
    Phenotype = im_xcell_phenotype[rownames(expression_matrix_ordered), "subtype"]
  )
  phenotype_order <- order(sample_annotation$Phenotype,decreasing = T)
  expression_matrix_ordered <- expression_matrix_ordered[phenotype_order, ]
  sample_annotation_ordered <- sample_annotation[phenotype_order, , drop = FALSE]
  cold_count <- sum(sample_annotation_ordered$Phenotype == "cold")
  feature_annotation <- data.frame(
    row.names = ordered_features,
    FeatureType = ifelse(ordered_features %in% signature_features, "Signature", "CellType"),
    Regulation = ifelse(ordered_features %in% c("Hot_signature", up_cells), "Up", "Down")
  )
  annotation_colors <- list(
    Phenotype = c(cold = "#90d7ec", hot = "#f58f98"),
    FeatureType = c(Signature = "#f39c12", CellType = "#9b59b6"),
    Regulation = c(Up = "#e74c3c", Down = "#3498db")
  )
  data_range <- quantile(expression_matrix_ordered, probs = c(0.05, 0.95), na.rm = TRUE)
  color_limit <- max(abs(data_range))
  breaks <- seq(-color_limit, color_limit, length.out = 100)
  gaps_row_positions <- c(
    1 + length(up_cells)
  )
  p <- pheatmap(t(expression_matrix_ordered),
                annotation_col = sample_annotation_ordered,
                annotation_row = feature_annotation,
                annotation_colors = annotation_colors,
                scale = "row",
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                show_colnames = FALSE,
                fontsize_row = 8,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                breaks = breaks,
                main = paste("Immune Cell Infiltration and Cancer Signatures\n", fs[i]),
                border_color = NA,
                gaps_row = gaps_row_positions,
                gaps_col = cold_count,
                silent = TRUE)

  output_dir <- "C:/Users/user/Desktop/input/heatmaps"
  if(!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  filename <- paste0(output_dir, "/","a_", fs[i], "_consistent_immune_heatmap.png")
  ggsave(filename, p$gtable, width = 12, height = 8, dpi = 300,bg = "white")
}
