# ==============================================================================
# WntQuant revision -- Per-cohort immunotherapy analyses
# ------------------------------------------------------------------------------
# Purpose:    Run the Wnt classification and survival workflow cohort by
#             cohort across the public immunotherapy datasets.
# Inputs:     Cohort expression matrices (<Dataset>_G.Rdata) and clinical
#             tables (<Dataset>_info) under input/.
# Outputs:    Per-cohort subtype labels, Hot/Cold groupings, Kaplan-Meier
#             panels and ROC output.
# ==============================================================================

################################################################################IMvigor####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("C:/Users/user/Desktop/input/IMvigor.Rdata")
IMVigor210=log2(IMvigor+1)
save(IMVigor210,file="IMVigor210.Rdata")
IMVigor210_info=data.frame(tag=rownames(phenoData),
                        cancer="BLCA",
                        responder=phenoData$binaryResponse,
                        os.time=phenoData$os,
                        os.status=phenoData$censOS,
                        pfs.time=NA,
                        pfs.status=NA)
IMVigor210_info$responder=ifelse(IMVigor210_info$responder%in%c("CR/PR"),"TRUE",
                              ifelse(IMVigor210_info$responder%in%c("SD/PD"),"FALSE",NA))
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
save(IMVigor210,file="IMVigor210.Rdata")
save(IMVigor210_info,file="IMVigor210_info.Rdata")
################################################################################PRECOG####
rm(list=ls())
gc()
process_expression_matrix <- function(expr_matrix) {
  if (any(expr_matrix < 0, na.rm = TRUE)) {
    cat("Values below 0 detected; applying a minimum-offset correction...\n")
    min_val <- min(expr_matrix, na.rm = TRUE)
    expr_matrix <- expr_matrix - min_val
    cat("Minimum value subtracted", min_val, "\n")
  } else {
    cat("No values below 0 detected\n")
  }
  if (any(expr_matrix > 30, na.rm = TRUE)) {
    cat("Values above 30 detected; applying log2(x+1)...\n")
    expr_matrix <- log2(expr_matrix + 1)
    cat("log2 transformation done\n")
  } else {
    cat("Values within 0-30; left unchanged\n")
  }
  return(expr_matrix)
}
convert_counts_to_cpm = function(expr_matrix) {

  mat <- as.matrix(expr_matrix)

  integer_ratio <- mean(mat == floor(mat), na.rm = TRUE)
  non_negative_ratio <- mean(mat >= 0, na.rm = TRUE)

  if (integer_ratio > 0.8 && non_negative_ratio > 0.95) {
    library(edgeR)
    mat=cpm(mat)
    mat=as.data.frame(mat)
    return(mat)
  } else {
    cat("Not recognised as a counts matrix; left unchanged\n")
    return(mat)
  }
}
setwd("C:/Users/user/Desktop/input")
ici_precog_exp=readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/ici_precog_all_patient_expression_data_FOR_DOWNLOAD.rds")
ici_precog_info=read.csv("ici_precog_all_patient_clinical_data_FOR_DOWNLOAD.csv")
datasets=unique(ici_precog_info$Dataset)
for(i in 1:length(datasets)){
  ici_precog_info_single = ici_precog_info[ici_precog_info$Dataset %in% datasets[i], ]
  match_counts = c()
  match_counts["Patient_ID"] = sum(ici_precog_info_single$Patient_ID %in% rownames(ici_precog_exp))
  match_counts["Sample_ID"] = sum(ici_precog_info_single$Sample_ID %in% rownames(ici_precog_exp))
  match_counts["geo_accession"] = sum(ici_precog_info_single$geo_accession %in% rownames(ici_precog_exp))
  best_match_col = names(which.max(match_counts))
  ici_precog_exp_single = ici_precog_exp[rownames(ici_precog_exp) %in% ici_precog_info_single[[best_match_col]], ]
  cat("Matches per column:\n")
  print(match_counts)
  cat("Best column selected:", best_match_col, "\n")
  cat("Rows after filtering:", nrow(ici_precog_exp_single), "\n")
  ici_precog_exp_single=as.data.frame(t(ici_precog_exp_single))
  has_nan <- any(sapply(ici_precog_exp_single, function(x) any(is.nan(x))))

  if (has_nan) {
    cat("NaN values found; converting to NA...\n")

    mat <- as.matrix(ici_precog_exp_single)
    mat[is.nan(mat)] <- NA
    ici_precog_exp_single <- as.data.frame(mat)

    cat("NaN conversion done\n")
  } else {
    cat("No NaN values found\n")
  }
  ici_precog_exp_single = ici_precog_exp_single[!apply(ici_precog_exp_single, 1, function(x) all(is.na(x))), ]
  if (any(is.na(ici_precog_exp_single))) {
    cat("Missing values found; imputing with KNN...\n")
    cat("Total missing values:", sum(is.na(ici_precog_exp_single)), "\n")
    library(DMwR2)
    ici_precog_exp_single = knnImputation(ici_precog_exp_single, k = 3)
    cat("KNN imputation done\n")
    if (any(is.na(ici_precog_exp_single))) {
      cat("Warning: missing values remain after imputation:", sum(is.na(ici_precog_exp_single)), "\n")
    } else {
      cat("All missing values imputed\n")
    }
  } else {
    cat("No missing values; nothing to impute\n")
  }
  #@@@#
  ici_precog_exp_single = convert_counts_to_cpm(ici_precog_exp_single)
  ici_precog_exp_single = process_expression_matrix(ici_precog_exp_single)
  assign(datasets[i],as.data.frame(ici_precog_exp_single))
  ici_precog_info_single=ici_precog_info_single[,c(best_match_col,
                                                   "TCGA_Tissue",
                                                   "Responder",
                                                   "OS_Days",
                                                   "OS_Censor",
                                                   "PFS_Days",
                                                   "PFS_Censor")]
  colnames(ici_precog_info_single)=c("tag","cancer","responder","os.time","os.status","pfs.time","pfs.status")
  assign(paste0(datasets[i],"_info"),ici_precog_info_single)
}
setwd("C:/Users/user/Desktop/input/ICI")
datasets = datasets[datasets != c("IMVigor210","Zhao_NatMed_2019","Gide_Cell_2019")]
for(dataset in datasets) {
  if(exists(dataset)) {
    save_file <- paste0(dataset, ".Rdata")
    save(list = dataset, file = save_file)
    cat("Saved:", save_file, "\n")
  }
  info_name <- paste0(dataset, "_info")
  if(exists(info_name)) {
    save_info_file <- paste0(info_name, ".Rdata")
    save(list = info_name, file = save_info_file)
    cat("Saved:", save_info_file, "\n")
  }
}
################################################################################PRJNA482620####
rm(list=ls())
gc()
PRJNA482620 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/GBM-PRJNA482620.Response.Rds")
PRJNA482620$GENE_SYMBOL = sub("_.*", "", PRJNA482620$GENE_SYMBOL)
PRJNA482620=PRJNA482620[!grepl("^[0-9]", PRJNA482620$GENE_SYMBOL), ]
PRJNA482620=aggregate(. ~ GENE_SYMBOL, data = PRJNA482620, FUN = mean)
rownames(PRJNA482620)=PRJNA482620[,1]
PRJNA482620=PRJNA482620[,-1]
PRJNA482620=log2(PRJNA482620+1)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("GBM-PRJNA482620.Response.tsv")
Zhao_NatMed_2019=PRJNA482620
Zhao_NatMed_2019_info=data.frame(tag=anno$sample_id,
                           cancer="GBM",
                           responder=anno$response,
                           os.time=anno$overall.survival..days.,
                           os.status=anno$vital.status,
                           pfs.time=NA,
                           pfs.status=NA)
Zhao_NatMed_2019_info$responder=ifelse(Zhao_NatMed_2019_info$responder=="R","TRUE",
                                       ifelse(Zhao_NatMed_2019_info$responder=="N","FALSE",NA))
Zhao_NatMed_2019_info$os.status=ifelse(Zhao_NatMed_2019_info$os.status=="Dead",1,0)
save(Zhao_NatMed_2019,file="Zhao_NatMed_2019.Rdata")
save(Zhao_NatMed_2019_info,file="Zhao_NatMed_2019_info.Rdata")
################################################################################PRJEB25780####
rm(list=ls())
gc()
PRJEB25780 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/STAD-PRJEB25780.Response.Rds")
PRJEB25780$GENE_SYMBOL = sub("_.*", "", PRJEB25780$GENE_SYMBOL)
PRJEB25780=PRJEB25780[!grepl("^[0-9]", PRJEB25780$GENE_SYMBOL), ]
PRJEB25780=aggregate(. ~ GENE_SYMBOL, data = PRJEB25780, FUN = mean)
rownames(PRJEB25780)=PRJEB25780[,1]
PRJEB25780=PRJEB25780[,-1]
PRJEB25780=log2(PRJEB25780+1)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("STAD-PRJEB25780.Response.tsv")
Kim_NatMed_2018=PRJEB25780
Kim_NatMed_2018_info=data.frame(tag=anno$sample_id,
                                 cancer="STAD",
                                 responder=anno$response_NR,
                                 os.time=anno$overall.survival..days.,
                                 os.status=anno$vital.status,
                                 pfs.time=NA,
                                 pfs.status=NA)
Kim_NatMed_2018_info$responder=ifelse(Kim_NatMed_2018_info$responder=="R","TRUE",
                                       ifelse(Kim_NatMed_2018_info$responder=="N","FALSE",NA))
save(Kim_NatMed_2018,file="Kim_NatMed_2018.Rdata")
save(Kim_NatMed_2018_info,file="Kim_NatMed_2018_info.Rdata")
################################################################################GSE67501####
rm(list=ls())
gc()
GSE67501 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/RCC-GSE67501.Response.Rds")
GSE67501 = GSE67501[!grepl("^[0-9]", GSE67501$GENE_SYMBOL), ]
GSE67501 = GSE67501[!apply(GSE67501[, -1], 1, function(x) all(is.na(x))), ]
GSE67501 = aggregate(. ~ GENE_SYMBOL, data = GSE67501, FUN = mean)
rownames(GSE67501)=GSE67501[,1]
GSE67501=GSE67501[,-1]
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("RCC-GSE67501.Response.tsv")
GSE67501_info=data.frame(tag=anno$sample_id,
                         cancer="KIRC",
                         responder=anno$response_NR,
                         os.time=anno$overall.survival..days.,
                         os.status=anno$vital.status,
                         pfs.time=NA,
                         pfs.status=NA)
GSE67501_info$responder=ifelse(GSE67501_info$responder=="R","TRUE",
                                      ifelse(GSE67501_info$responder=="N","FALSE",NA))
save(GSE67501,file="GSE67501.Rdata")
save(GSE67501_info,file="GSE67501_info.Rdata")
################################################################################Braun_2020####
rm(list=ls())
gc()
Braun_NatMed_2020 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/RCC-Braun_2020.Response.Rds")
Braun_NatMed_2020 = Braun_NatMed_2020[!grepl("^[0-9]", Braun_NatMed_2020$GENE_SYMBOL), ]
Braun_NatMed_2020 = Braun_NatMed_2020[!apply(Braun_NatMed_2020[, -1], 1, function(x) all(is.na(x))), ]
Braun_NatMed_2020 = aggregate(. ~ GENE_SYMBOL, data = Braun_NatMed_2020, FUN = mean)
rownames(Braun_NatMed_2020)=Braun_NatMed_2020[,1]
Braun_NatMed_2020=Braun_NatMed_2020[,-1]
Braun_NatMed_2020=log2(Braun_NatMed_2020+1)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("RCC-Braun_2020.Response.tsv")
Braun_NatMed_2020_info=data.frame(tag=anno$sample_id,
                         cancer="KIRC",
                         responder=anno$response_NR,
                         os.time=anno$overall.survival..days.,
                         os.status=anno$vital.status,
                         pfs.time=NA,
                         pfs.status=NA)
Braun_NatMed_2020_info$responder=ifelse(Braun_NatMed_2020_info$responder=="R","TRUE",
                               ifelse(Braun_NatMed_2020_info$responder=="N","FALSE",NA))
Braun_NatMed_2020_info$os.status=ifelse(Braun_NatMed_2020_info$os.status=="Dead",1,0)
save(Braun_NatMed_2020,file="Braun_NatMed_2020.Rdata")
save(Braun_NatMed_2020_info,file="Braun_NatMed_2020_info.Rdata")
################################################################################GSE93157####
rm(list=ls())
gc()
GSE93157 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/nonsqNSCLC-GSE93157.Response.Rds")
GSE93157 = GSE93157[!grepl("^[0-9]", GSE93157$GENE_SYMBOL), ]
GSE93157 = GSE93157[!apply(GSE93157[, -1], 1, function(x) all(is.na(x))), ]
GSE93157 = aggregate(. ~ GENE_SYMBOL, data = GSE93157, FUN = mean)
rownames(GSE93157)=GSE93157[,1]
GSE93157=GSE93157[,-1]
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("nonsqNSCLC-GSE93157.Response.tsv")
GSE93157_info=data.frame(tag=anno$sample_id,
                                  cancer="NSCLC",
                                  responder=anno$response_NR,
                                  os.time=anno$overall.survival..days.,
                                  os.status=anno$vital.status,
                                  pfs.time=NA,
                                  pfs.status=NA)
GSE93157_info$responder=ifelse(GSE93157_info$responder=="R","TRUE",
                                        ifelse(GSE93157_info$responder=="N","FALSE",NA))
GSE93157_info$os.status=ifelse(GSE93157_info$os.status=="Dead",1,0)
save(GSE93157,file="GSE93157.Rdata")
save(GSE93157_info,file="GSE93157_info.Rdata")
################################################################################PRJEB23709####
rm(list=ls())
gc()
PRJEB23709 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/Melanoma-PRJEB23709.Response.Rds")
PRJEB23709$GENE_SYMBOL = sub("_.*", "", PRJEB23709$GENE_SYMBOL)
PRJEB23709=PRJEB23709[!grepl("^[0-9]", PRJEB23709$GENE_SYMBOL), ]
PRJEB23709=aggregate(. ~ GENE_SYMBOL, data = PRJEB23709, FUN = mean)
rownames(PRJEB23709)=PRJEB23709[,1]
PRJEB23709=PRJEB23709[,-1]
PRJEB23709=log2(PRJEB23709+1)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("Melanoma-PRJEB23709.Response.tsv")
Gide_CancerCell_2019=PRJEB23709
Gide_CancerCell_2019_info=data.frame(tag=anno$sample_id,
                                 cancer="SKCM",
                                 responder=anno$response_NR,
                                 os.time=anno$overall.survival..days.,
                                 os.status=anno$vital.status,
                                 pfs.time=NA,
                                 pfs.status=NA)
Gide_CancerCell_2019_info$responder=ifelse(Gide_CancerCell_2019_info$responder=="R","TRUE",
                                       ifelse(Gide_CancerCell_2019_info$responder=="N","FALSE",NA))
Gide_CancerCell_2019_info$os.status=ifelse(Gide_CancerCell_2019_info$os.status=="Dead",1,0)
save(Gide_CancerCell_2019,file="Gide_CancerCell_2019.Rdata")
save(Gide_CancerCell_2019_info,file="Gide_CancerCell_2019_info.Rdata")
################################################################################GSE115821####
rm(list=ls())
gc()
GSE115821 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/Melanoma-GSE115821.Response.Rds")
GSE115821 = GSE115821[!grepl("^[0-9]", GSE115821$GENE_SYMBOL), ]
GSE115821 = GSE115821[!apply(GSE115821[, -1], 1, function(x) all(is.na(x))), ]
GSE115821 = aggregate(. ~ GENE_SYMBOL, data = GSE115821, FUN = mean)
rownames(GSE115821)=GSE115821[,1]
GSE115821=GSE115821[,-1]
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("Melanoma-GSE115821.Response.tsv")
GSE115821_info=data.frame(tag=anno$sample_id,
                         cancer="SKCM",
                         responder=anno$response_NR,
                         os.time=anno$overall.survival..days.,
                         os.status=anno$vital.status,
                         pfs.time=NA,
                         pfs.status=NA)
GSE115821_info$responder=ifelse(GSE115821_info$responder=="R","TRUE",
                               ifelse(GSE115821_info$responder=="N","FALSE",NA))
GSE115821_info$os.status=ifelse(GSE115821_info$os.status=="Dead",1,0)
save(GSE115821,file="GSE115821.Rdata")
save(GSE115821_info,file="GSE115821_info.Rdata")
################################################################################GSE106128####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
GSE106128 = read.delim("GSE106128_series_matrix.txt")
ann=read.delim("GSE106128_family.soft")
ann=ann[,c(1,7)]
ann=ann[ann$ILMN_Gene!="",]
GSE106128=merge(ann,GSE106128,by.x="ID",by.y="ID_REF")
GSE106128=GSE106128[,-1]
GSE106128 = aggregate(. ~ ILMN_Gene, data = GSE106128, FUN = mean)
GSE106128 = GSE106128[!grepl("^[0-9]", GSE106128$ILMN_Gene), ]
rownames(GSE106128)=GSE106128[,1]
GSE106128=GSE106128[,-1]
GSE106128=log2(GSE106128+1)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("Melanoma-GSE106128.Response.tsv")
GSE106128_info=data.frame(tag=anno$sample_id,
                          cancer="SKCM",
                          responder=anno$response_NR,
                          os.time=anno$overall.survival..days.,
                          os.status=anno$vital.status,
                          pfs.time=NA,
                          pfs.status=NA)
GSE106128_info$responder=ifelse(GSE106128_info$responder=="R","TRUE",
                                ifelse(GSE106128_info$responder=="N","FALSE",NA))
GSE106128_info$os.status=ifelse(GSE106128_info$os.status=="Dead",1,0)
save(GSE106128,file="GSE106128.Rdata")
save(GSE106128_info,file="GSE106128_info.Rdata")
################################################################################GSE100797####
rm(list=ls())
gc()
GSE100797 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/Melanoma-GSE100797.Response.Rds")
GSE100797 = GSE100797[!grepl("^[0-9]", GSE100797$GENE_SYMBOL), ]
GSE100797 = GSE100797[!apply(GSE100797[, -1], 1, function(x) all(is.na(x))), ]
GSE100797 = aggregate(. ~ GENE_SYMBOL, data = GSE100797, FUN = mean)
rownames(GSE100797)=GSE100797[,1]
GSE100797=GSE100797[,-1]
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("Melanoma-GSE100797.Response.tsv")
GSE100797_info=data.frame(tag=anno$sample_id,
                          cancer="SKCM",
                          responder=anno$response_NR,
                          os.time=anno$overall.survival..days.,
                          os.status=anno$vital.status,
                          pfs.time=NA,
                          pfs.status=NA)
GSE100797_info$responder=ifelse(GSE100797_info$responder=="R","TRUE",
                                ifelse(GSE100797_info$responder=="N","FALSE",NA))
GSE100797_info$os.status=ifelse(GSE100797_info$os.status=="Dead",1,0)
save(GSE100797,file="GSE100797.Rdata")
save(GSE100797_info,file="GSE100797_info.Rdata")
################################################################################GSE93157####
rm(list=ls())
gc()
GSE93157 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/Melanoma-GSE93157.Response.Rds")
GSE93157 = GSE93157[!grepl("^[0-9]", GSE93157$GENE_SYMBOL), ]
GSE93157 = GSE93157[!apply(GSE93157[, -1], 1, function(x) all(is.na(x))), ]
GSE93157 = aggregate(. ~ GENE_SYMBOL, data = GSE93157, FUN = mean)
rownames(GSE93157)=GSE93157[,1]
GSE93157=GSE93157[,-1]
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("Melanoma-GSE93157.Response.tsv")
GSE93157_info=data.frame(tag=anno$sample_id,
                          cancer="LUSC",
                          responder=anno$response_NR,
                          os.time=anno$overall.survival..days.,
                          os.status=anno$vital.status,
                          pfs.time=NA,
                          pfs.status=NA)
GSE93157_info$responder=ifelse(GSE93157_info$responder=="R","TRUE",
                                ifelse(GSE93157_info$responder=="N","FALSE",NA))
GSE93157_info$os.status=ifelse(GSE93157_info$os.status=="Dead",1,0)
save(GSE93157,file="GSE93157.Rdata")
save(GSE93157_info,file="GSE93157_info.Rdata")
################################################################################E-MTAB-6270####
rm(list=ls())
gc()
Müller_JCI_2018 = readRDS("C:/Users/user/Desktop/input/Wnt_immunotherapy/LGG_E-MTAB-6270.Response.Rds")
Müller_JCI_2018 = Müller_JCI_2018[!grepl("^[0-9]", Müller_JCI_2018$GENE_SYMBOL), ]
Müller_JCI_2018 = Müller_JCI_2018[!apply(Müller_JCI_2018[, -1], 1, function(x) all(is.na(x))), ]
Müller_JCI_2018 = aggregate(. ~ GENE_SYMBOL, data = Müller_JCI_2018, FUN = mean)
rownames(Müller_JCI_2018)=Müller_JCI_2018[,1]
Müller_JCI_2018=Müller_JCI_2018[,-1]
Müller_JCI_2018=log2(Müller_JCI_2018+1)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
anno=read.delim("LGG_E-MTAB-6270.Response.tsv")
Müller_JCI_2018_info=data.frame(tag=anno$sample_id,
                         cancer="LGG",
                         responder=anno$response_NR,
                         os.time=anno$overall.survival..days.,
                         os.status=anno$vital.status,
                         pfs.time=NA,
                         pfs.status=NA)
Müller_JCI_2018_info$responder=ifelse(Müller_JCI_2018_info$responder=="R","TRUE",
                               ifelse(Müller_JCI_2018_info$responder=="N","FALSE",NA))
Müller_JCI_2018_info$os.status=ifelse(Müller_JCI_2018_info$os.status=="Dead",1,0)
save(Müller_JCI_2018,file="Müller_JCI_2018.Rdata")
save(Müller_JCI_2018_info,file="Müller_JCI_2018_info.Rdata")
################################################################################GSE100797####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
exp=read.delim("GSE100797_ProcessedData.txt")
rownames(exp)=exp[,1]
exp=exp[,-1]
GSE100797=exp
data=read.delim("GSE100797_series_matrix.txt")[c(13,14,11,12,18),]
data=as.data.frame(t(data))
data=data[-1,]
colnames(data)=c("os","os.event","pfs","pfs.event","responder")
GSE100797_info=data.frame(tag=rownames(data),
                                cancer="SKCM",
                                responder=data$responder,
                                os.time=data$os,
                                os.status=data$os.event,
                                pfs.time=data$pfs,
                                pfs.status=data$pfs.event)
GSE100797_info$os.time=as.numeric(str_extract(GSE100797_info$os.time, "(?<=: ).*$"))
GSE100797_info$os.status=as.numeric(str_extract(GSE100797_info$os.status, "(?<=: ).*$"))
GSE100797_info$pfs.time=as.numeric(str_extract(GSE100797_info$pfs.time, "(?<=: ).*$"))
GSE100797_info$pfs.status=as.numeric(str_extract(GSE100797_info$pfs.status, "(?<=: ).*$"))
GSE100797_info$responder=ifelse(GSE100797_info$responder%in%c("recist: CR","recist: PR"),"TRUE","FALSE")
save(GSE100797,file="GSE100797.Rdata")
save(GSE100797_info,file="GSE100797_info.Rdata")
################################################################################GSE135222####
rm(list=ls())
gc()
library(IOBR)
library(dplyr)
library(stringr)
setwd("C:/Users/user/Desktop/input")
exp=read.delim("GSE135222_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
ann=read.delim("Human.GRCh38.p13.annot.tsv")[,c(1,2)]
exp=merge(ann,exp,by="GeneID")
exp=exp[,-1]
exp=aggregate(. ~ Symbol, data = exp, FUN = mean)
rownames(exp)=exp[,1]
exp=exp[,-1]
exp=log2(exp+1)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/ICI/GSE135222_info.Rdata")
GSE135222_info1=read.delim("GSE135222_series_matrix.txt")
GSE135222_info1=as.data.frame(t(GSE135222_info1[c(12,13,41),]))
GSE135222_info1=GSE135222_info1[-1,]
colnames(GSE135222_info1)=c("status","time","id")
GSE135222_info1$status=as.numeric(str_extract(GSE135222_info1$status, "(?<=: ).*$"))
GSE135222_info1$time=as.numeric(str_extract(GSE135222_info1$time, "(?<=: ).*$"))
GSE135222_info=merge(GSE135222_info,GSE135222_info1,by.x="tag",by.y="id")
GSE135222_info=data.frame(tag=GSE135222_info$tag,
                          cancer="NSCLC",
                          responder=GSE135222_info$responder,
                          os.time=NA,
                          os.status=NA,
                          pfs.time=GSE135222_info$time,
                          pfs.status=GSE135222_info$status)
GSE135222=exp
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
save(GSE135222,file="GSE135222.Rdata")
save(GSE135222_info,file="GSE135222_info.Rdata")
################################################################################GSE19423####
rm(list=ls())
gc()
library(dplyr)
library(stringr)
setwd("C:/Users/user/Desktop/input")
exp=read.delim("GSE19423_GPL6102_1_13_1_Matrix.matrix.txt")
rownames(exp)=exp[,1]
exp=exp[,-1]
GSE19423=exp
data=read.delim("GSE19423_series_matrix.txt")[c(1,13,19,21),]
data=as.data.frame(t(data))
data=data[-1,]
colnames(data)=c("tag","responder","status","time")
data$status=ifelse(data$status=="overall survival: death",1,0)
data$time=as.numeric(str_extract(data$time, "(?<=: ).*$"))
GSE19423_info=data.frame(tag=data$tag,
                         cancer="SKCM",
                         responder=data$responder,
                         os.time=data$time,
                         os.status=data$status,
                         pfs.time=NA,
                         pfs.status=NA)
GSE19423_info$responder=ifelse(GSE19423_info$responder=="previous mapki: N","FALSE","TRUE")
GSE19423_info$os.status=ifelse(GSE19423_info$os.status=="overall survival: death",1,0)
GSE19423=exp
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
save(GSE19423,file="GSE19423.Rdata")
save(GSE19423_info,file="GSE19423_info.Rdata")
################################################################################GSE100797####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
exp=read.delim("GSE100797_ProcessedData.txt")
rownames(exp)=exp[,1]
exp=exp[,-1]
GSE100797=exp
data=read.delim("GSE100797_series_matrix.txt")[c(13,14,11,12,18),]
data=as.data.frame(t(data))
data=data[-1,]
colnames(data)=c("os","os.event","pfs","pfs.event","responder")
GSE100797_info=data.frame(tag=rownames(data),
                          cancer="SKCM",
                          responder=data$responder,
                          os.time=data$os,
                          os.status=data$os.event,
                          pfs.time=data$pfs,
                          pfs.status=data$pfs.event)
GSE100797_info$os.time=as.numeric(str_extract(GSE100797_info$os.time, "(?<=: ).*$"))
GSE100797_info$os.status=as.numeric(str_extract(GSE100797_info$os.status, "(?<=: ).*$"))
GSE100797_info$pfs.time=as.numeric(str_extract(GSE100797_info$pfs.time, "(?<=: ).*$"))
GSE100797_info$pfs.status=as.numeric(str_extract(GSE100797_info$pfs.status, "(?<=: ).*$"))
GSE100797_info$responder=ifelse(GSE100797_info$responder%in%c("recist: CR","recist: PR"),"TRUE","FALSE")

save(GSE100797,file="GSE100797.Rdata")
save(GSE100797_info,file="GSE100797_info.Rdata")
################################################################################GSE78220####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
library(readxl)
exp=as.data.frame(read_excel("GSE78220_PatientFPKM.xlsx"))
rownames(exp)=exp[,1]
exp=exp[,-1]
exp=log2(exp+1)
data=read.delim("GSE78220_series_matrix.txt")[c(9,15,16,17),]
data=as.data.frame(t(data))
data=data[-1,]
colnames(data)=c("tag","time","status","responder")
GSE78220_info=data.frame(tag=rownames(data),
                         cancer="SKCM",
                         responder=data$responder,
                         os.time=data$time,
                         os.status=data$status,
                         pfs.time=NA,
                         pfs.status=NA)
GSE78220_info$tag=paste0(GSE78220_info$tag,".baseline")
GSE78220_info$responder=ifelse(GSE78220_info$responder=="previous mapki: N","FALSE","TRUE")
GSE78220_info$os.time=as.numeric(str_extract(GSE78220_info$os.time, "(?<=: ).*$"))
GSE78220_info$os.status=ifelse(GSE78220_info$os.status=="vital status: Alive",0,1)
GSE78220=exp
save(GSE78220,file="GSE78220.Rdata")
save(GSE78220_info,file="GSE78220_info.Rdata")
################################################################################Rebuilding the machine-learning model with predicted probabilities####
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
resultup=resultup[!resultup$gene_id%in%c("Sebocytes","GMP","Keratinocytes","Astrocytes"),]
save(resultup,file="resultup.Rdata")
exp_im=exp_im[c(resultup$gene_id),]
exp_im=as.data.frame(t(exp_im))
load("C:/Users/user/Desktop/input/hot_cancer_sig.Rdata")
exp_hot=IMvigor[hot_cancer_sig$gene_id,]
cor_matrix <- cor(t(exp_hot), method = "pearson")
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "ward.D")
modules = cutree(hc, k = 26)
module_list <- list()
for(i in 1:26) {
  module_name <- paste0("Hot_cancer_signature_", i)
  genes_in_module <- names(modules)[modules == i]
  module_list[[module_name]] <- genes_in_module
}
Hot_cancer_signature_list=module_list
exp_for_yao
#save(Hot_cancer_signature_list,file="Hot_cancer_signature_list.Rdata")
library(GSVA)
scores=ssgseaParam(exprData=as.matrix(IMvigor),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
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
signature_features <- grep("Hot_cancer_signature", colnames(expression_matrix_ordered), value = TRUE)
celltype_features <- setdiff(colnames(expression_matrix_ordered), signature_features)
up_cells <- intersect(celltype_features, resultup$gene_id)
ordered_features <- c(sort(signature_features), sort(up_cells))
expression_matrix_ordered_final <- expression_matrix_ordered[, ordered_features]
feature_annotation <- data.frame(
  row.names = ordered_features,
  FeatureType = ifelse(grepl("Hot_cancer_signature", ordered_features), "Signature", "CellType"),
  Regulation = ifelse(grepl("Hot_cancer_signature", ordered_features), "Signature", "Up-regulated")
)
phenotype_colors <- list(
  Phenotype = c(cold = "#3498db", hot = "#e74c3c"),
  Regulation = c(Signature = "#f39c12", "Up-regulated" = "#e74c3c"),
  FeatureType = c(Signature = "#f39c12", CellType = "#2ecc71")
)
data_range <- quantile(expression_matrix_ordered_final, probs = c(0.05, 0.95), na.rm = TRUE)
color_limit <- max(abs(data_range))
breaks <- seq(-color_limit, color_limit, length.out = 100)
cold_samples_count <- sum(phenotype_annotation_ordered$Phenotype == "cold")
signature_count <- length(signature_features)
gaps_row_positions <- signature_count
pheatmap(t(expression_matrix_ordered_final),
         annotation_col = phenotype_annotation_ordered,
         annotation_row = feature_annotation,
         annotation_colors = phenotype_colors,
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         fontsize_col = 8,
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         breaks = breaks,
         main = "Hot Cancer Signatures and Up-regulated Immune Cells\nby Tumor Phenotype",
         border_color = NA,
         gaps_col = cold_samples_count,
         gaps_row = gaps_row_positions,
         angle_col = "45",
         cellwidth = NA,
         cellheight = NA,
         legend = TRUE,
         annotation_legend = TRUE,
         drop_levels = TRUE)
################################################################################Obtaining the hot labels####
rm(list=ls())
gc()
library(IOBR)
library(GSVA)
library(BiocParallel)

# ===== Patch: replace the internal rawEnrichmentAnalysis of IOBR with the new ssGSEA API
raw_old <- getFromNamespace("rawEnrichmentAnalysis", "IOBR")
raw_fixed <- raw_old

body(raw_fixed) <- substitute({

  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))

  expr <- expr[shared.genes, , drop = FALSE]

  if (nrow(expr) < 5000) {
    stop("Not enough genes")
  }

  expr <- apply(expr, 2, rank)

  # Core fix: use ssGSEA rather than GSVA
  params <- GSVA::ssgseaParam(
    exprData = as.matrix(expr),
    geneSets = signatures,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = FALSE
  )

  scores <- GSVA::gsva(
    params,
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )

  # Keep the original IOBR post-processing
  scores <- scores - apply(scores, 1, min)

  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]

  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  if (!is.null(file.name)) {
    write.table(
      scores,
      file = file.name,
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
  }

  scores
})

assignInNamespace(
  "rawEnrichmentAnalysis",
  raw_fixed,
  ns = "IOBR"
)
# ===== End of patch

setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("C:/Users/user/Desktop/input/IMvigor.Rdata")
IMvigor=log2(IMvigor+1)
im_xcell=deconvo_tme(eset = IMvigor, method = "xcell")

# Everything below stays as is
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

# Keep all data for now; do not filter by p-value
volcano_data <- resultup
volcano_data$logP <- -log10(volcano_data$P.Value)

# Set the significance groups
volcano_data$significance <- ifelse(volcano_data$P.Value < 0.05, "Significant", "Not significant")
volcano_data$group <- ifelse(volcano_data$logFC > 0, "Up in hot", "Down in hot")

# Build the combined colour and shape groups
volcano_data$plot_group <- ifelse(
  volcano_data$significance == "Not significant", "Not significant",
  ifelse(volcano_data$gene_id %in% c("Sebocytes", "GMP", "Keratinocytes", "Astrocytes"),
         "Excluded cells",
         volcano_data$group)
)

library(ggplot2)
library(ggrepel)

p <- ggplot(volcano_data, aes(x = logFC, y = logP)) +
  geom_point(data = subset(volcano_data, significance == "Not significant"),
             color = "grey", alpha = 0.5, size = 2) +
  geom_point(data = subset(volcano_data, significance == "Significant"),
             aes(color = plot_group), alpha = 0.7, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", alpha = 0.7, linewidth = 0.8) +
  scale_color_manual(
    values = c("Up in hot" = "#f7acbc",
               "Down in hot" = "blue",
               "Excluded cells" = "#90d7ec",
               "Not significant" = "grey"),
    breaks = c("Up in hot", "Down in hot", "Excluded cells"),
    labels = c("Up in hot", "Down in hot", "Excluded cells")
  ) +
  labs(
    x = expression(log[2]("FC")),
    y = expression(-log[10]("p")),
    color = "Expression"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    legend.position = "top",
    legend.justification = "center",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

p <- p + geom_text_repel(
  data = subset(volcano_data, significance == "Significant"),
  aes(label = gene_id),
  color = "black",
  size = 3,
  box.padding = 0.5,
  max.overlaps = Inf,
  show.legend = FALSE,
  min.segment.length = 0
)

print(p)


################################################################################Data preparation for the random-forest regression####
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
load("resultup.Rdata")
im_xcell=im_xcell[,colnames(im_xcell)%in%resultup$gene_id]
load("Hot_cancer_signature_list.Rdata")
library(openxlsx)

# Find the length of the longest gene set
max_length <- max(sapply(Hot_cancer_signature_list, length))

# Create a data frame and pad the short rows with NA
hot_cancer_wide <- as.data.frame(
  do.call(rbind,
          lapply(Hot_cancer_signature_list, function(x) {
            c(x, rep(NA, max_length - length(x)))
          })
  )
)

# Add the gene-set name as the first column
hot_cancer_wide <- cbind(GeneSet_Name = names(Hot_cancer_signature_list), hot_cancer_wide)

# Set the column names
colnames(hot_cancer_wide) <- c("GeneSet_Name", paste0("Gene_", 1:max_length))

# Export to Excel
write.xlsx(hot_cancer_wide, "Hot_cancer_signatures_wide.xlsx", rowNames = FALSE)
library(GSVA)
scores=ssgseaParam(exprData=as.matrix(IMvigor),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
scores=gsva(scores)
scores=as.data.frame(scale(t(scores)))
im_xcell=merge(im_xcell,scores,by="row.names")
rownames(im_xcell)=im_xcell[,1]
im_xcell=im_xcell[,-1]
phenoData=na.omit(phenoData[,6,drop=F])
xcell_pheno=merge(phenoData,im_xcell,by="row.names")
rownames(xcell_pheno)=xcell_pheno[,1]
xcell_pheno=xcell_pheno[,-1]
colnames(xcell_pheno)[1]=c("phenotype")
table(xcell_pheno$phenotype)
xcell_pheno$phenotype=ifelse(xcell_pheno$phenotype=="inflamed","hot","cold")
table(xcell_pheno$phenotype)
save(xcell_pheno,file="xcell_pheno.Rdata")
################################################################################Building the random-forest model####
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/xcell_pheno.Rdata")
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(vip)
library(dplyr)
library(purrr)
library(ggthemes)
library(scales)
library(RColorBrewer)
library(tidyr)
xcell_pheno$phenotype = ifelse(xcell_pheno$phenotype == "hot", 1, 0)
xcell_pheno$phenotype = as.factor(xcell_pheno$phenotype)
processed_data = xcell_pheno
colnames(processed_data)[1] = "target"
processed_data$target_factor = factor(ifelse(processed_data$target == 1, "hot", "cold"),
                                       levels = c("hot", "cold"))
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
n_vars = ncol(processed_data) - 2
mtry_range = 1:20
cat("Number of features:", n_vars, "\n")
cat("mtry tuning range:", mtry_range, "\n")
cat("Training Random Forest model...\n")
rf_model = train(
  x = processed_data[, -c(1, ncol(processed_data))],
  y = processed_data$target_factor,
  method = "rf",
  trControl = train_control,
  tuneGrid = expand.grid(mtry = mtry_range),
  ntree = 1000,
  importance = TRUE,
  metric = "ROC"
)
save(rf_model,file = "rf_model.Rdata")
######################Plotting the optimal parameter
rm(list=ls())
gc()
load("rf_model.Rdata")
print(rf_model$bestTune)
tune_results <- rf_model$pred
tune_results$pred=ifelse(tune_results$hot>0.5,"hot","cold")
tune_results$pred <- factor(tune_results$pred, levels = c("hot", "cold"))
tune_results$obs <- factor(tune_results$obs, levels = c("hot", "cold"))
fs=unique(tune_results$mtry)
AUC_all=data.frame()
for(i in 1:length(fs)){
  tune_results_fs=tune_results[tune_results$mtry==fs[i],]
  roc_obj = roc(response = tune_results_fs$obs, predictor = tune_results_fs$hot)
  auc_value = auc(roc_obj)
  AUC=data.frame(auc=auc_value,mtry=fs[i])
  AUC_all=rbind(AUC_all,AUC)
}
library(ggplot2)
ggplot(AUC_all, aes(x = mtry, y = auc)) +
  geom_line(color = "gray60", linewidth = 0.8) +
  geom_vline(xintercept = 18, color = "#f7acbc", linetype = "dashed", linewidth = 0.8) +
  geom_point(aes(color = ifelse(mtry == 18, "Optimal parameter", "Other parameters")),
             size = 5) +
  geom_point(data = subset(AUC_all, mtry == 18),
             color = "#f7acbc", size = 5, shape = 1, stroke = 1.5) +
  scale_color_manual(values = c("Optimal parameter" = "#f7acbc", "Other parameters" = "#90d7ec")) +
  scale_x_continuous(breaks = seq(2, 20, by = 2)) +
  labs(title = "Parameter tuning",
       subtitle = "",
       x = "Parameter(mtry)",
       y = "AUC value",
       color = "Parameter type") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, color = "black"),
    plot.subtitle = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )
################################################################################Visualising the model####
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/xcell_pheno.Rdata")
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(vip)
library(dplyr)
library(purrr)
library(ggthemes)
library(scales)
library(RColorBrewer)
library(tidyr)
xcell_pheno$phenotype = ifelse(xcell_pheno$phenotype == "hot", 1, 0)
xcell_pheno$phenotype = as.factor(xcell_pheno$phenotype)
processed_data = xcell_pheno
colnames(processed_data)[1] = "target"
processed_data$target_factor = factor(ifelse(processed_data$target == 1, "hot", "cold"),
                                      levels = c("hot", "cold"))
set.seed(123)
cat("Training Random Forest Model with Optimal Parameters (mtry=4)...\n")
train_control_final <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  selectionFunction = "best",
  savePredictions = TRUE,
  verboseIter = FALSE
)
rf_model_final <- train(
  x = processed_data[, -c(1, ncol(processed_data))],
  y = processed_data$target_factor,
  method = "rf",
  trControl = train_control_final,
  tuneGrid = data.frame(mtry = 18),
  ntree = 1000,
  importance = TRUE,
  metric = "ROC"
)
save(rf_model_final, file = "rf_model_final.Rdata")
cv_predictions <- rf_model_final$pred
cv_predictions$pred=ifelse(cv_predictions$hot>0.5,"hot","cold")
cv_predictions$pred <- factor(cv_predictions$pred, levels = c("hot", "cold"))
cv_predictions$obs <- factor(cv_predictions$obs, levels = c("hot", "cold"))
# Extract the ROC data of each cross-validation fold
roc_data_list <- list()
folds <- unique(cv_predictions$Resample)

for (i in seq_along(folds)) {
  fold_data <- cv_predictions[cv_predictions$Resample == folds[i], ]
  roc_fold <- roc(response = fold_data$obs, predictor = fold_data$hot)
  roc_data_list[[i]] <- data.frame(
    sensitivity = roc_fold$sensitivities,
    specificity = roc_fold$specificities,
    fold = paste0("Fold_", i)
  )
}

roc_data_all <- bind_rows(roc_data_list)

# Compute the mean ROC curve and its confidence interval
roc_obj <- roc(response = cv_predictions$obs, predictor = cv_predictions$hot)
auc_value <- auc(roc_obj)

# Prepare the error-bar data by summarising the sensitivity at several 1-specificity points
specificity_seq <- seq(0, 1, by = 0.01)
sensitivity_matrix <- matrix(NA, nrow = length(folds), ncol = length(specificity_seq))

for (i in seq_along(folds)) {
  fold_data <- cv_predictions[cv_predictions$Resample == folds[i], ]
  roc_fold <- roc(response = fold_data$obs, predictor = fold_data$hot)

  # Interpolate the sensitivity at each specificity
  sensitivity_matrix[i, ] <- approx(
    x = 1 - roc_fold$specificities,
    y = roc_fold$sensitivities,
    xout = specificity_seq,
    method = "linear",
    rule = 2
  )$y
}

# Compute the mean and SD of each point
sensitivity_mean <- apply(sensitivity_matrix, 2, mean, na.rm = TRUE)
sensitivity_sd <- apply(sensitivity_matrix, 2, sd, na.rm = TRUE)

# Build the error-bar data frame
error_ribbon_data <- data.frame(
  specificity = specificity_seq,
  sensitivity_mean = sensitivity_mean,
  sensitivity_lower = pmax(0, sensitivity_mean - sensitivity_sd),
  sensitivity_upper = pmin(1, sensitivity_mean + sensitivity_sd)
)

# Build the ROC plot with error bars
auc_plot <- ggplot() +
  # Draw the mean ROC curve
  geom_line(data = data.frame(
    sensitivity = roc_obj$sensitivities,
    specificity = roc_obj$specificities
  ), aes(x = 1 - specificity, y = sensitivity),
  color = "#2E86AB", size = 1.5) +
  # Add the error band
  geom_ribbon(data = error_ribbon_data,
              aes(x = specificity,
                  ymin = sensitivity_lower,
                  ymax = sensitivity_upper),
              fill = "#2E86AB", alpha = 0.2) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = paste0("ROC curve with variability band"),
       x = "1 - Specificity",
       y = "Sensitivity") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", size = 0.8),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  coord_equal() +
  # Add the AUC annotations
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("Mean AUC = ", round(auc_value, 3)),
           size = 5, color = "#2E86AB", fontface = "bold")

print(auc_plot)

importance_data <- varImp(rf_model_final, scale = TRUE)$importance
if (!"Overall" %in% colnames(importance_data)) {
  # If the column has a different name, use the first numeric column
  numeric_cols <- sapply(importance_data, is.numeric)
  if (any(numeric_cols)) {
    importance_data$Overall <- importance_data[, which(numeric_cols)[1]]
  } else {
    stop("No numeric columns found in importance data")
  }
}

importance_data$Variable <- rownames(importance_data)
importance_data <- importance_data[order(-importance_data$Overall), ]
top_20_importance <- head(importance_data, 20)
name_mapping <- c(
  "pDC" = "Plasmacytoid dendritic cells",
  "CD8+_Tem" = "CD8+ effector memory T cells",
  "CD8+_Tcm" = "CD8+ central memory T cells",
  "Th1_cells" = "Type 1 T helper cells",
  "Th2_cells" = "Type 2 T helper cells",
  "Tgd_cells" = "γδ T cells"
)
top_20_importance$Variable <- rownames(top_20_importance)
top_20_importance$DisplayName <- top_20_importance$Variable
for (old_name in names(name_mapping)) {
  top_20_importance$DisplayName[top_20_importance$DisplayName == old_name] <- name_mapping[old_name]
}
top_20_importance$DisplayName <- gsub("Hot_cancer_signature_", "HCS", top_20_importance$DisplayName)
importance_plot <- ggplot(top_20_importance, aes(x = reorder(DisplayName, Overall), y = Overall)) +
  geom_col(fill = "#2E86AB", alpha = 0.9, width = 0.8) +
  geom_text(aes(label = round(Overall, 1)),
            hjust = -0.2, size = 4.5, color = "black") +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Top 20 variables",
       x = NULL,
       y = "Importance score") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(importance_plot)
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
  load("resultup.Rdata")
  im_xcell=im_xcell[,c(resultup$gene_id)]
  load("Hot_cancer_signature_list.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_xcell=merge(im_xcell,scores,by="row.names")
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  load("rf_model_final.Rdata")
  im_xcell_phenotype=data.frame(Tag=rownames(im_xcell),subtype=predict(object = rf_model_final,
                                                                       newdata = as.matrix(im_xcell),type = "prob"))
  im_xcell_phenotype$subtype=ifelse(im_xcell_phenotype$subtype.hot>0.5,"hot","cold")
  table(im_xcell_phenotype$subtype)
  rownames(im_xcell_phenotype) <- im_xcell_phenotype$Tag

  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)

  expression_matrix <- as.matrix(im_xcell)

  # Get all column names related to Hot_cancer_signature
  signature_features <- grep("^Hot_cancer_signature", colnames(expression_matrix), value = TRUE)
  celltype_features <- setdiff(colnames(expression_matrix), signature_features)
  up_cells <- intersect(celltype_features, resultup$gene_id)

  # Reorder features: Hot_cancer_signature first, then the up-regulated cell types
  ordered_features <- c(signature_features, up_cells)

  expression_matrix_ordered <- expression_matrix[, ordered_features]
  sample_annotation <- data.frame(
    row.names = rownames(expression_matrix_ordered),
    Phenotype = im_xcell_phenotype[rownames(expression_matrix_ordered), "subtype"]
  )
  phenotype_order <- order(sample_annotation$Phenotype)
  expression_matrix_ordered <- expression_matrix_ordered[phenotype_order, ]
  sample_annotation_ordered <- sample_annotation[phenotype_order, , drop = FALSE]
  cold_count <- sum(sample_annotation_ordered$Phenotype == "cold")

  # Simplify the feature annotation, keeping only FeatureType
  feature_annotation <- data.frame(
    row.names = ordered_features,
    FeatureType = ifelse(grepl("^Hot_cancer_signature", ordered_features), "Signature", "CellType")
  )

  # Simplify the annotation colours, keeping only Phenotype and FeatureType
  annotation_colors <- list(
    Phenotype = c(cold = "#90d7ec", hot = "#f58f98"),
    FeatureType = c(Signature = "#f39c12", CellType = "#9b59b6")
  )

  data_range <- quantile(expression_matrix_ordered, probs = c(0.05, 0.95), na.rm = TRUE)
  color_limit <- max(abs(data_range))
  breaks <- seq(-color_limit, color_limit, length.out = 100)

  # Rows are reordered: insert a separator between all Hot_cancer_signature entries and the cell types
  gaps_row_positions <- length(signature_features)  # Insert a separator after all signature features

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

  filename <- paste0(output_dir, "/","z_", fs[i], "_consistent_immune_heatmap.png")
  ggsave(filename, p$gtable, width = 12, height = 8, dpi = 300,bg = "white")
}
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
  load("resultup.Rdata")
  im_xcell=im_xcell[,c(resultup$gene_id)]
  load("Hot_cancer_signature_list.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_xcell=merge(im_xcell,scores,by="row.names")
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  load("rf_model_final.Rdata")
  im_xcell_phenotype=data.frame(Tag=rownames(im_xcell),subtype=predict(object = rf_model_final,
                                                                       newdata = as.matrix(im_xcell),type = "prob"))
  im_xcell_phenotype$subtype=ifelse(im_xcell_phenotype$subtype.hot>0.5,"hot","cold")
  table(im_xcell_phenotype$subtype)
  rownames(im_xcell_phenotype) <- im_xcell_phenotype$Tag

  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)

  expression_matrix <- as.matrix(im_xcell)

  # Get all column names related to Hot_cancer_signature
  signature_features <- grep("^Hot_cancer_signature", colnames(expression_matrix), value = TRUE)
  celltype_features <- setdiff(colnames(expression_matrix), signature_features)
  up_cells <- intersect(celltype_features, resultup$gene_id)

  # Reorder features: Hot_cancer_signature first, then the up-regulated cell types
  ordered_features <- c(signature_features, up_cells)

  expression_matrix_ordered <- expression_matrix[, ordered_features]
  sample_annotation <- data.frame(
    row.names = rownames(expression_matrix_ordered),
    Phenotype = im_xcell_phenotype[rownames(expression_matrix_ordered), "subtype"]
  )
  phenotype_order <- order(sample_annotation$Phenotype)
  expression_matrix_ordered <- expression_matrix_ordered[phenotype_order, ]
  sample_annotation_ordered <- sample_annotation[phenotype_order, , drop = FALSE]
  cold_count <- sum(sample_annotation_ordered$Phenotype == "cold")

  # Simplify the feature annotation, keeping only FeatureType
  feature_annotation <- data.frame(
    row.names = ordered_features,
    FeatureType = ifelse(grepl("^Hot_cancer_signature", ordered_features), "Signature", "CellType")
  )

  # Simplify the annotation colours, keeping only Phenotype and FeatureType
  annotation_colors <- list(
    Phenotype = c(cold = "#90d7ec", hot = "#f58f98"),
    FeatureType = c(Signature = "#f39c12", CellType = "#9b59b6")
  )

  data_range <- quantile(expression_matrix_ordered, probs = c(0.05, 0.95), na.rm = TRUE)
  color_limit <- max(abs(data_range))
  breaks <- seq(-color_limit, color_limit, length.out = 100)

  # Rows are reordered: insert a separator between all Hot_cancer_signature entries and the cell types
  gaps_row_positions <- length(signature_features)  # Insert a separator after all signature features

  p <- pheatmap(t(expression_matrix_ordered),
                annotation_col = sample_annotation_ordered,
                annotation_row = feature_annotation,
                annotation_colors = annotation_colors,
                scale = "row",
                cluster_rows = FALSE,
                cluster_cols = T,
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

  filename <- paste0(output_dir, "/","w_", fs[i], "_consistent_immune_heatmap.png")
  ggsave(filename, p$gtable, width = 12, height = 8, dpi = 300,bg = "white")
}
################################################################################cold/hot proportion statistics across pan-cancer####
rm(list=ls())
gc()
library(IOBR)
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/pancancer_exp.Rdata")
load("C:/Users/user/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
fs=unique(pancancer_group$Project)
cold_hot_all=data.frame()
for(i in 1:length(fs)){
  pancancer_group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  exp=pancancer_exp[,rownames(pancancer_group_fs)]
  im_xcell=deconvo_tme(eset = exp, method = "xcell")
  im_xcell=as.data.frame(im_xcell)
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  colnames(im_xcell)=gsub("_xCell$", "", colnames(im_xcell))
  im_xcell=as.data.frame(scale(im_xcell))
  load("resultup.Rdata")
  im_xcell=im_xcell[,c(resultup$gene_id)]
  load("Hot_cancer_signature_list.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_xcell=merge(im_xcell,scores,by="row.names")
  rownames(im_xcell)=im_xcell[,1]
  im_xcell=im_xcell[,-1]
  load("rf_model_final.Rdata")
  im_xcell_phenotype=data.frame(Tag=rownames(im_xcell),subtype=predict(object = rf_model_final,
                                                                       newdata = as.matrix(im_xcell),type = "prob"))
  im_xcell_phenotype$subtype=ifelse(im_xcell_phenotype$subtype.hot>0.5,"hot","cold")
  table(im_xcell_phenotype$subtype)
  cold_hot=data.frame(Tag=im_xcell_phenotype$Tag,group=im_xcell_phenotype$subtype)
  cold_hot$cancer=fs[i]
  cold_hot_all=rbind(cold_hot_all,cold_hot)
}
save(cold_hot_all,file="cold_hot_all.Rdata")
library(ggplot2)
library(dplyr)

# Calculate proportions of cold and hot for each cancer type
cold_hot_summary <- cold_hot_all %>%
  count(cancer, group) %>%
  group_by(cancer) %>%
  mutate(
    count = n,
    percentage = n / sum(n) * 100
  ) %>%
  ungroup()

# View the summary
print(cold_hot_summary)

# Option 1: Stacked bar plot
ggplot(cold_hot_summary, aes(x = cancer, y = percentage, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Proportion of Cold and Hot Samples by Cancer Type",
       x = "Cancer Type",
       y = "Percentage (%)",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
################################################################################OS analysis and random-forest analysis####
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
  load("C:/Users/user/Desktop/input/resultup.Rdata")
  im_timer=im_timer[,c(resultup$gene_id)]
  load("C:/Users/user/Desktop/input/Hot_cancer_signature_list.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_timer=merge(im_timer,scores,by="row.names")
  rownames(im_timer)=im_timer[,1]
  im_timer=im_timer[,-1]
  load("C:/Users/user/Desktop/input/rf_model_final.Rdata")
  im_timer$subtype = predict(object = rf_model_final,
                             newdata = as.matrix(im_timer),type = "prob")[, "hot"]
  im_timer$subtype=ifelse(im_timer$subtype>0.5,"hot","cold")
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
  load("C:/Users/user/Desktop/input/resultup.Rdata")
  im_timer=im_timer[,c(resultup$gene_id)]
  load("C:/Users/user/Desktop/input/Hot_cancer_signature_list.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_timer=merge(im_timer,scores,by="row.names")
  rownames(im_timer)=im_timer[,1]
  im_timer=im_timer[,-1]
  load("C:/Users/user/Desktop/input/rf_model_final.Rdata")
  im_timer$subtype = predict(object = rf_model_final,
                             newdata = as.matrix(im_timer),type = "prob")[, "hot"]
  im_timer$subtype=ifelse(im_timer$subtype>0.5,"hot","cold")
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
vector=unique(c("GSE253260","Dai_2025","Prins_GBM_2019_RNAseq","GSE135222","Dai_2025",
                "GSE253260","PORTER_2025","Zhao_NatMed_2019","Braun_NatMed_2020","GSE78220"))
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
  load("C:/Users/user/Desktop/input/resultup.Rdata")
  im_timer=im_timer[,c(resultup$gene_id)]
  load("C:/Users/user/Desktop/input/Hot_cancer_signature_list.Rdata")
  library(GSVA)
  scores=ssgseaParam(exprData=as.matrix(exp),geneSets=Hot_cancer_signature_list,minSize=1,maxSize=10000)
  scores=gsva(scores)
  scores=as.data.frame(scale(t(scores)))
  im_timer=merge(im_timer,scores,by="row.names")
  rownames(im_timer)=im_timer[,1]
  im_timer=im_timer[,-1]
  load("C:/Users/user/Desktop/input/rf_model_final.Rdata")
  im_timer$subtype = predict(object = rf_model_final,
                             newdata = as.matrix(im_timer),type = "prob")[, "hot"]
  im_timer$subtype=ifelse(im_timer$subtype>0.5,"hot","cold")
  table(im_timer$subtype)
  group=data.frame(Tag=rownames(im_timer),group=im_timer$subtype)
  assign(paste0(vector[i],"_G"),group)
  save(list=paste0(vector[i],"_G"),file=paste0(vector[i],"_G.Rdata"))
}
################################################################################Group percentages####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")

vector <- unique(c("GSE253260","Dai_2025","Prins_GBM_2019_RNAseq","Dai_2025",
                   "GSE253260","PORTER_2025","Zhao_NatMed_2019","Braun_NatMed_2020","GSE78220","IMvigor"))

# Create an empty data frame to collect the group information of every dataset
all_data <- data.frame()

for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_G.Rdata")

  # Check whether the file exists
  if (file.exists(data_file) & file.exists(info_file)) {
    load(data_file)
    load(info_file)

    # Get the group data
    # Assume the object in the _G file is named dataset_G (e.g. GSE253260_G)
    dataset_name <- vector[q]
    group_var_name <- paste0(dataset_name, "_G")

    if (exists(group_var_name)) {
      group_data <- get(group_var_name)

      # Check whether a group column is present
      if ("group" %in% colnames(group_data)) {
        # Count the hot and cold cases
        count_data <- as.data.frame(table(group_data$group))
        colnames(count_data) <- c("Group", "Count")
        count_data$Dataset <- dataset_name
        count_data$Percentage <- count_data$Count / sum(count_data$Count) * 100

        all_data <- rbind(all_data, count_data)
      }
    }
  }
}

# Create the stacked percentage plot
library(ggplot2)
library(dplyr)
library(scales)

# Make Group a factor with hot first
all_data$Group <- factor(all_data$Group, levels = c("hot", "cold"))
all_data$Dataset <- ifelse(all_data$Dataset=="GSE253260","Nicolle_EBioMedicine_2024",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="Dai_2025","Dai_Nature_2025",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="PORTER_2025","Galsky_ClinCancerRes_2025",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="GSE78220","Hugo_Cell_2016",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="Prins_GBM_2019_RNAseq","Cloughesy_NatMed_2019",all_data$Dataset)

# Compute the cumulative percentage of each bar for label placement
label_data <- all_data %>%
  group_by(Dataset) %>%
  arrange(desc(Group)) %>%
  mutate(
    cumulative_perc = cumsum(Percentage),
    label_pos = cumulative_perc - Percentage / 2
  )

p <- ggplot(all_data, aes(x = reorder(Dataset, -Percentage), y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  # Add the percentage labels with a smaller font
  geom_text(data = label_data,
            aes(x = Dataset, y = label_pos,
                label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
            size = 3, color = "black", fontface = "bold") +  # Change size from 13/3 to 3
  scale_fill_manual(values = c("hot" = "#f7acbc", "cold" = "#90d7ec")) +
  labs(x = "", y = "Percentage (%)", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 13, color = "black"),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(color = "black", size = 13)  # All text black, size 13
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100))

# Display in R
print(p)

################################################################################Single dataset (GSE135222)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE135222")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
GSE135222_G=GSE135222_G[GSE135222_G$group=="cold",]
GSE135222=GSE135222[,colnames(GSE135222)%in%GSE135222_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=GSE135222,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,GSE135222_info,by.x="row.names",by.y="tag")
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
  xlim = c(0, 36),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Jung_NatCommun_2019",
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
vector=c("GSE135222")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=GSE78220_G
group=group[group$group=="cold",]
exp=GSE78220
prognosis=GSE78220_info[,c("tag","pfs.time","pfs.status")]
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
  ggtitle("GSE78220(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Jung_NatCommun_2019_pfs_p=p_all
save(Jung_NatCommun_2019_pfs_p,file="Jung_NatCommun_2019_pfs_p.Rdata")


################################################################################Single dataset (Dai_Nature_2025_os)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_2025")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Dai_2025_G=Dai_2025_G[Dai_2025_G$group=="cold",]
Dai_2025=Dai_2025[,colnames(Dai_2025)%in%Dai_2025_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Dai_2025,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Dai_2025_info,by.x="row.names",by.y="tag")
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
  title="Dai_Nature_2025",
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
vector=c("Dai_2025")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Dai_2025_G
group=group[group$group=="cold",]
exp=Dai_2025
prognosis=Dai_2025_info[,c("tag","os.time","os.status")]
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
  ggtitle("Dai_Nature_2025(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Dai_Nature_2025_os_p=p_all
save(Dai_Nature_2025_os_p,file="Dai_Nature_2025_os_p.Rdata")
################################################################################Single dataset (Nicolle_EBioMedicine_2024_os)####
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
  xlim = c(0, 69),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Nicolle_EBioMedicine_2024",
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
Nicolle_EBioMedicine_2024_os_p=p_all
save(Nicolle_EBioMedicine_2024_os_p,file="Nicolle_EBioMedicine_2024_os_p.Rdata")
################################################################################Single dataset (Nicolle_EBioMedicine_2024_pfs)####
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
  title="Nicolle_EBioMedicine_2024",
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
Nicolle_EBioMedicine_2024_pfs_p=p_all
save(Nicolle_EBioMedicine_2024_pfs_p,file="Nicolle_EBioMedicine_2024_pfs_p.Rdata")
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
################################################################################Single dataset (Hugo_Cell_2016)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE78220")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
GSE78220_G=GSE78220_G[GSE78220_G$group=="cold",]
GSE78220=GSE78220[,colnames(GSE78220)%in%GSE78220_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=GSE78220,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,GSE78220_info,by.x="row.names",by.y="tag")
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
  xlim = c(0, 36),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Hugo_Cell_2016",
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
vector=c("GSE78220")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=GSE78220_G
group=group[group$group=="cold",]
exp=GSE78220
prognosis=GSE78220_info[,c("tag","os.time","os.status")]
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
  ggtitle("GSE78220(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Hugo_Cell_2016_os_p=p_all
save(Hugo_Cell_2016_os_p,file="Hugo_Cell_2016_os_p.Rdata")
################################################################################Single dataset (PORTER_2025)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("PORTER_2025")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
PORTER_2025_G=PORTER_2025_G[PORTER_2025_G$group=="cold",]
PORTER_2025=PORTER_2025[,colnames(PORTER_2025)%in%PORTER_2025_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=PORTER_2025,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,PORTER_2025_info,by.x="row.names",by.y="tag")
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
  xlim = c(0, 36),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Galsky_ClinCancerRes_2025",
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
vector=c("PORTER_2025")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=PORTER_2025_G
group=group[group$group=="cold",]
exp=PORTER_2025
prognosis=PORTER_2025_info[,c("tag","os.time","os.status")]
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
  ggtitle("PORTER_2025(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Galsky_ClinCancerRes_2025_os_p=p_all
save(Galsky_ClinCancerRes_2025_os_p,file="Galsky_ClinCancerRes_2025_os_p.Rdata")
################################################################################Single dataset (Dai_Nature_2025_pfs)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_2025")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Dai_2025_G=Dai_2025_G[Dai_2025_G$group=="cold",]
Dai_2025=Dai_2025[,colnames(Dai_2025)%in%Dai_2025_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Dai_2025,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Dai_2025_info,by.x="row.names",by.y="tag")
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
  title="Dai_Nature_2025",
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
vector=c("Dai_2025")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Dai_2025_G
group=group[group$group=="cold",]
exp=Dai_2025
prognosis=Dai_2025_info[,c("tag","pfs.time","pfs.status")]
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
  ggtitle("Dai_Nature_2025(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Dai_Nature_2025_pfs_p=p_all
save(Dai_Nature_2025_pfs_p,file="Dai_Nature_2025_pfs_p.Rdata")
################################################################################Single dataset (Prins_GBM_2019_RNAseq)####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Prins_GBM_2019_RNAseq")
for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")
  info_file = paste0(vector[q], "_info.Rdata")
  group_file = paste0(vector[q], "_G.Rdata")
  load(data_file)
  load(info_file)
  load(group_file)
}
Prins_GBM_2019_RNAseq_G=Prins_GBM_2019_RNAseq_G[Prins_GBM_2019_RNAseq_G$group=="cold",]
Prins_GBM_2019_RNAseq=Prins_GBM_2019_RNAseq[,colnames(Prins_GBM_2019_RNAseq)%in%Prins_GBM_2019_RNAseq_G$Tag]
library(GSEABase)
geneSets=getGmt("C:/Users/user/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/user/Desktop/functions_to_package/Calculate_bioactivity_scores.R")
WNT_Score=Calculate_bioactivity_scores(file_paths="C:/Users/user/Desktop/input",
                                       expression_profile=Prins_GBM_2019_RNAseq,
                                       foundation="relative_ssGSEA",
                                       activation_geneset=NA,
                                       inhibition_geneset=NA,
                                       geneSets_gmt=geneSets,
                                       min.sz=1,
                                       max.sz=10000,
                                       geneset_direction=c(rep("pos",sum(grepl("WPAGS", names(geneSets)))),
                                                           rep("neg",sum(grepl("WPIGS", names(geneSets))))))
WNT_score=WNT_Score$final_activity_score
wnt_IMvigor=merge(WNT_score,Prins_GBM_2019_RNAseq_info,by.x="row.names",by.y="tag")
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
  xlim = c(0, 15),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Low","Medium","High"),
  xlab = "Time (months)",
  ylab = "Probability",
  title="Cloughesy_NatMed_2019",
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
vector=c("Prins_GBM_2019_RNAseq")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_G.Rdata"))
}
group=Prins_GBM_2019_RNAseq_G
group=group[group$group=="cold",]
exp=Prins_GBM_2019_RNAseq
prognosis=Prins_GBM_2019_RNAseq_info[,c("tag","pfs.time","pfs.status")]
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
  #"Antigen_processing_machinery",
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
  ggtitle("Dai_Nature_2025(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) + # Add horizontal line
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)  # Add horizontal-line labels and adjust the font size
Cloughesy_NatMed_2019_pfs_p=p_all
save(Cloughesy_NatMed_2019_pfs_p,file="Cloughesy_NatMed_2019_pfs_p.Rdata")
################################################################################Heatmap visualisation of everything####
rm(list = ls())
gc()
setwd("C:/Users/user/Desktop/input")
vector <- c("Dai_Nature_2025_os", "Hugo_Cell_2016_os", "Braun_NatMed_2020_os", "Nicolle_EBioMedicine_2024_os", "Zhao_NatMed_2019_os","Cloughesy_NatMed_2019_pfs",
            "Dai_Nature_2025_pfs","Nicolle_EBioMedicine_2024_pfs","IMvigor_os","Galsky_ClinCancerRes_2025_os")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/", vector[l], "_p.Rdata"))
  exp=get(paste0(vector[l], "_p"))
  colnames(exp)=c(vector[l],"ID")
  assign(paste0(vector[l], "_p"),exp)
}
exp=merge(Dai_Nature_2025_os_p,Hugo_Cell_2016_os_p,by="ID")
exp=merge(exp,Braun_NatMed_2020_os_p,by="ID")
exp=merge(exp,Nicolle_EBioMedicine_2024_os_p,by="ID")
exp=merge(exp,Zhao_NatMed_2019_os_p,by="ID")
exp=merge(exp,Cloughesy_NatMed_2019_pfs_p,by="ID")
exp=merge(exp,Dai_Nature_2025_pfs_p,by="ID")
exp=merge(exp,Nicolle_EBioMedicine_2024_pfs_p,by="ID")
exp=merge(exp,IMvigor_os_p,by="ID")
exp=merge(exp,Galsky_ClinCancerRes_2025_os_p,by="ID")
rownames(exp)=exp[,1]
exp=exp[,-1]
exp$sum=rowMeans(exp)

# Sort by sum
exp_sorted <- exp[order(exp$sum), ]
annotation_data <- data.frame(Sum = exp_sorted$sum)
rownames(annotation_data) <- rownames(exp_sorted)
heatmap_data <- exp_sorted[, -which(colnames(exp_sorted) == "sum")]

library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tibble)

# Define the annotation categories
annotation_categories <- list(
  "Wnt features" = c("Wnt/β-catenin pathway activity score", "Wnt/β-catenin pathway activation score", "Wnt/β-catenin pathway inhibition score"),
  "Wnt-related genes" = c("TCF7", "AXIN2", "RNF43", "ZNRF3", "LGR5"),
  "T cell-related" = c("T cell regulatory", "T cell exhaustion", "T cell inflamed GEP", "T cell accumulation", "CD8+ T effector"),
  "Other indicators" = c("TAM", "MDSC", "CAF", "IFNG signature"),
  "Immune checkpoint" = c("ICB resistance", "Immune checkpoint", "CD274", "PDCD1"),
  "TGF-β" = c("Pan-F-TBRS", "TGFB1")
)

# Build the annotation data frame
annotation_df <- data.frame(feature = rownames(heatmap_data))
annotation_df$category <- NA

for (category in names(annotation_categories)) {
  features_in_category <- annotation_categories[[category]]
  annotation_df$category[annotation_df$feature %in% features_in_category] <- category
}

# Sort by category, then by sum within each category
heatmap_df <- as.data.frame(heatmap_data) %>%
  tibble::rownames_to_column("feature") %>%
  left_join(annotation_df, by = "feature") %>%
  # Add the sum values
  mutate(sum = exp_sorted$sum[match(feature, rownames(exp_sorted))]) %>%
  # Sort by category, then by sum descending
  arrange(factor(category, levels = names(annotation_categories)), sum) %>%
  mutate(feature = factor(feature, levels = feature))  # Keep the sort order

# 2. Convert to long format for the heatmap
heatmap_long <- heatmap_df %>%
  pivot_longer(-c(feature, category, sum), names_to = "dataset", values_to = "value")

# 3. Build the data for the significance asterisks
significance_df <- heatmap_long %>%
  mutate(star = ifelse(value > -log10(0.001), "***",
                       ifelse(value > -log10(0.01), "**",
                              ifelse(value > -log10(0.05), "*", ""))))

# 4. Draw the heatmap
heatmap_plot <- ggplot(heatmap_long, aes(x = dataset, y = feature, fill = value)) +
  geom_tile(color = "#a1a3a6", size = 0.5) +
  geom_text(data = significance_df, aes(label = star),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",
                       midpoint = median(heatmap_long$value)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "", fill = expression(-log[10](italic(p)))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "bottom"
  )

# 5. Prepare the bar chart data
bar_data <- data.frame(
  feature = heatmap_df$feature,
  importance = heatmap_df$sum,
  category = heatmap_df$category
)

# 6. Draw the bar chart
barplot_plot <- ggplot(bar_data, aes(x = feature, y = importance, fill = category)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "Mean", fill = "Category") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10)
  ) +
  scale_fill_manual(values = c(
    "Wnt features" = "#1f77b4",
    "Wnt-related genes" = "#ff7f0e",
    "T cell-related" = "#2ca02c",
    "Other indicators" = "#d62728",
    "Immune checkpoint" = "#9467bd",
    "TGF-β" = "#8c564b"
  ))

# 7. Combine the figures
combined_plot <- heatmap_plot + barplot_plot +
  plot_layout(widths = c(3, 1.5))  # Adjust the width ratio to fit the legend

print(combined_plot)
################################################################################IMvigor_TME####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("C:/Users/user/Desktop/input/IMvigor.Rdata")
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
IMvigor=log2(IMvigor+1)
res_ensemble = tme_classifier(eset = IMvigor,method ="ensemble", scale =T)[,c(1,5)]
colnames(res_ensemble)=c("ID","Ensemble")
save(res_ensemble,file="res_ensemble.Rdata")
res_svm = tme_classifier(eset = IMvigor,method ="svm", scale =T)[,c(1,5)]
colnames(res_svm)=c("ID","SVM")
save(res_svm,file="res_svm.Rdata")
res_rf = tme_classifier(eset = IMvigor,method ="rf", scale =T)[,c(1,5)]
colnames(res_rf)=c("ID","RF")
save(res_rf,file="res_rf.Rdata")
res_nnet = tme_classifier(eset = IMvigor,method ="nnet", scale =T)[,c(1,5)]
colnames(res_nnet)=c("ID","NNET")
save(res_nnet,file="res_nnet.Rdata")
res_knn = tme_classifier(eset = IMvigor,method ="knn", scale =T)[,c(1,5)]
colnames(res_knn)=c("ID","KNN")
save(res_knn,file="res_knn.Rdata")
res_dt = tme_classifier(eset = IMvigor,method ="dt", scale =T)[,c(1,5)]
colnames(res_dt)=c("ID","DT")
save(res_dt,file="res_dt.Rdata")
res_xgboost = tme_classifier(eset = IMvigor,method ="xgboost", scale =T)[,c(1,5)]
colnames(res_xgboost)=c("ID","XGBoost")
save(res_xgboost,file="res_xgboost.Rdata")
res_GT=data.frame(ID=rownames(phenoData),GT=phenoData$`Immune phenotype`)
res_all = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                  list(res_ensemble, res_svm, res_rf, res_nnet,
                       res_knn, res_dt, res_xgboost))
res_all = merge(res_all, res_GT, by = "ID", all = TRUE)
head(res_all)
res_all = res_all[!is.na(res_all$GT), ]
res_all$GT = ifelse(res_all$GT %in% c("desert", "excluded"), "Cold", "Hot")
pred_cols = setdiff(colnames(res_all), c("ID", "GT"))
res_all[pred_cols] = lapply(res_all[pred_cols], function(x) {
  ifelse(x == "IA", "Hot", "Cold")
})
head(res_all)
library(ggplot2)
pred_cols <- c("Ensemble", "SVM", "RF", "NNET", "KNN", "DT", "XGBoost")
acc_df <- data.frame(
  Method = pred_cols,
  Accuracy = sapply(pred_cols, function(m) {
    mean(res_all[[m]] == res_all$GT, na.rm = TRUE)
  })
)
acc_df <- acc_df[order(acc_df$Accuracy), ]
acc_df$Method <- factor(acc_df$Method, levels = acc_df$Method)

ggplot(acc_df, aes(x = Method, y = Accuracy, fill = Accuracy)) +
  geom_col(width = 0.7, color = "black") +
  geom_text(aes(label = scales::percent(Accuracy, accuracy = 0.1)),
            hjust = -0.1, size = 4, color = "black", fontface = "plain") +
  scale_fill_gradient2(low = "#A6CEE3",   # Light blue
                       mid = "white",     # White
                       high = "#FB9A99",  # Light red
                       midpoint = median(acc_df$Accuracy)) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, max(acc_df$Accuracy) * 1.15)) +
  labs(x = "", y = "Accuracy", title = "") +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
library(ggplot2)

# Desktop path (Windows)
out_dir <- "C:/Users/user/Desktop"

# 7 methods
pred_cols <- c("Ensemble", "SVM", "RF", "NNET", "KNN", "DT", "XGBoost")

# Keep the factor order consistent
res_all$GT <- factor(res_all$GT, levels = c("Cold", "Hot"))
for (m in pred_cols) {
  res_all[[m]] <- factor(res_all[[m]], levels = c("Cold", "Hot"))
}

# Generate seven PDFs in turn
for (m in pred_cols) {

  # Confusion matrix of the current method
  cm_df <- as.data.frame(table(Predicted = res_all[[m]], Actual = res_all$GT))
  cm_df$Predicted <- factor(cm_df$Predicted, levels = c("Cold", "Hot"))
  cm_df$Actual    <- factor(cm_df$Actual,    levels = c("Cold", "Hot"))

  # Plot
  p <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
    geom_tile(color = "black") +
    geom_text(aes(label = Freq), color = "black", size = 5) +
    scale_fill_gradient2(low = "#A6CEE3", mid = "white", high = "#FB9A99",
                         midpoint = median(cm_df$Freq)) +
    labs(x = "Actual (GT)", y = "Predicted", fill = "Count",
         title = m) +
    theme_bw(base_size = 14) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      plot.title = element_text(size = 14, color = "black", hjust = 0.5),
      legend.position = "right"
    )

  # Save to the desktop
  ggsave(
    filename = file.path(out_dir, paste0("confusion_matrix_", m, ".pdf")),
    plot = p,
    width = 5,
    height = 4.5,
    device = "pdf"
  )
}
################################################################################Immunotherapy response####
rm(list=ls())
gc()
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/res_knn.Rdata")
wnt_IMvigor=merge(IMvigor_WNT_score$final_activity_score,res_knn,by.x="row.names",by.y="ID")
phenoData=phenoData[,"binaryResponse",drop=F]
wnt_IMvigor=merge(wnt_IMvigor,phenoData,by.x="Row.names",by.y="row.names")
wnt_IMvigor=wnt_IMvigor[,c("activity_score","binaryResponse","KNN")]
colnames(wnt_IMvigor)=c("activity_score","Response","KNN")
wnt_IMvigor=na.omit(wnt_IMvigor)
library(ggplot2)
library(ggpubr)
wnt_IMvigor$KNN=ifelse(wnt_IMvigor$KNN=="IA","Hot","Cold")
p = ggplot(wnt_IMvigor, aes(KNN, activity_score, fill = Response, color = Response)) +
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
################################################################################Prognostic analysis: high Wnt predicts a better outcome (KNN)####
rm(list=ls())
gc()

load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
load("C:/Users/user/Desktop/input/res_knn.Rdata")
table(res_knn$KNN)
# Merge the WNT score, phenoData and KNN labels
wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score,
                     phenoData,
                     by = "row.names")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

# IA -> Hot, everything else -> Cold
wnt_IMvigor$KNN <- ifelse(wnt_IMvigor$KNN == "IA", "Hot", "Cold")

# Key fix: keep levels as character, not factor
wnt_IMvigor$level <- as.character(wnt_IMvigor$KNN)
table(wnt_IMvigor$level)

cold_idx <- wnt_IMvigor$KNN == "Cold"

quartiles <- quantile(wnt_IMvigor$activity_score[cold_idx],
                      probs = seq(0, 1, by = 1/3),
                      na.rm = TRUE)

wnt_IMvigor$level[cold_idx] <- as.character(
  cut(wnt_IMvigor$activity_score[cold_idx],
      breaks = quartiles,
      include.lowest = TRUE,
      labels = c("Cold_low", "Cold_medium", "Cold_high"))
)

table(wnt_IMvigor$level, useNA = "ifany")

# Build the data for the four-group survival analysis
wnt_IMvigor_prognosis <- wnt_IMvigor[, c("activity_score", "os", "censOS", "level")]
wnt_IMvigor_prognosis <- na.omit(wnt_IMvigor_prognosis)

wnt_IMvigor_prognosis$level <- factor(
  wnt_IMvigor_prognosis$level,
  levels = c("Hot", "Cold_low", "Cold_medium", "Cold_high")
)

table(wnt_IMvigor_prognosis$level)

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# Overall log-rank test
fitd <- survdiff(Surv(os, censOS) ~ level,
                 data = wnt_IMvigor_prognosis,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# KM fit
fit <- survfit(Surv(os, censOS) ~ level,
               data = wnt_IMvigor_prognosis,
               type = "kaplan-meier",
               error = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# Pairwise comparisons
ps <- pairwise_survdiff(Surv(os, censOS) ~ level,
                        data = wnt_IMvigor_prognosis,
                        p.adjust.method = "none")

# Hot vs Cold_high computed manually
temp_data <- wnt_IMvigor_prognosis[wnt_IMvigor_prognosis$level %in% c("Hot", "Cold_high"), ]
temp_data$level <- droplevels(temp_data$level)

temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
hot_vs_cold_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)

hot_vs_cold_high_label <- ifelse(
  is.na(hot_vs_cold_high_p), "NA",
  ifelse(hot_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(hot_vs_cold_high_p < 0.05,
                paste0("p = ", round(hot_vs_cold_high_p, 3)), "ns"))
)

# Cold_low vs Cold_high computed manually
temp_data2 <- wnt_IMvigor_prognosis[wnt_IMvigor_prognosis$level %in% c("Cold_low", "Cold_high"), ]
temp_data2$level <- droplevels(temp_data2$level)

temp_fit2 <- survdiff(Surv(os, censOS) ~ level, data = temp_data2)
cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

mycol <- c('#ffd400', '#009ad6', "grey50", '#d71345')

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Hot", "Cold_low", "Cold_medium", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS (KNN, n = ", nrow(wnt_IMvigor_prognosis), ")"),
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

p.lab <- paste0("log-rank test p",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))

p$plot <- p$plot +
  annotate("text", x = 0, y = 0.02,
           label = p.lab, size = 6, hjust = 0, fontface = "italic") +
  annotate("text", x = 0, y = 0.10,
           label = paste0("Hot vs. Cold_high: ", hot_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic", color = "red")
print(p)
################################################################################Prognostic analysis: high Wnt predicts a better outcome (RF)####
rm(list=ls())
gc()

load("C:/Users/user/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/user/Desktop/input/IMvigor210CoreBiologies.Rdata")
load("C:/Users/user/Desktop/input/res_rf.Rdata")
table(res_rf$RF)
# Merge the WNT score, phenoData and KNN labels
wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score,
                     phenoData,
                     by = "row.names")

wnt_IMvigor <- merge(wnt_IMvigor, res_rf,
                     by.x = "Row.names", by.y = "ID")

# IA -> Hot, everything else -> Cold
wnt_IMvigor$RF <- ifelse(wnt_IMvigor$RF == "IA", "Hot", "Cold")

# Key fix: keep levels as character, not factor
wnt_IMvigor$level <- as.character(wnt_IMvigor$RF)
table(wnt_IMvigor$level)

cold_idx <- wnt_IMvigor$RF == "Cold"

quartiles <- quantile(wnt_IMvigor$activity_score[cold_idx],
                      probs = seq(0, 1, by = 1/3),
                      na.rm = TRUE)

wnt_IMvigor$level[cold_idx] <- as.character(
  cut(wnt_IMvigor$activity_score[cold_idx],
      breaks = quartiles,
      include.lowest = TRUE,
      labels = c("Cold_low", "Cold_medium", "Cold_high"))
)

table(wnt_IMvigor$level, useNA = "ifany")

# Build the data for the four-group survival analysis
wnt_IMvigor_prognosis <- wnt_IMvigor[, c("activity_score", "os", "censOS", "level")]
wnt_IMvigor_prognosis <- na.omit(wnt_IMvigor_prognosis)

wnt_IMvigor_prognosis$level <- factor(
  wnt_IMvigor_prognosis$level,
  levels = c("Hot", "Cold_low", "Cold_medium", "Cold_high")
)

table(wnt_IMvigor_prognosis$level)

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# Overall log-rank test
fitd <- survdiff(Surv(os, censOS) ~ level,
                 data = wnt_IMvigor_prognosis,
                 na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# KM fit
fit <- survfit(Surv(os, censOS) ~ level,
               data = wnt_IMvigor_prognosis,
               type = "kaplan-meier",
               error = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# Pairwise comparisons
ps <- pairwise_survdiff(Surv(os, censOS) ~ level,
                        data = wnt_IMvigor_prognosis,
                        p.adjust.method = "none")

# Hot vs Cold_high computed manually
temp_data <- wnt_IMvigor_prognosis[wnt_IMvigor_prognosis$level %in% c("Hot", "Cold_high"), ]
temp_data$level <- droplevels(temp_data$level)

temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
hot_vs_cold_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)

hot_vs_cold_high_label <- ifelse(
  is.na(hot_vs_cold_high_p), "NA",
  ifelse(hot_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(hot_vs_cold_high_p < 0.05,
                paste0("p = ", round(hot_vs_cold_high_p, 3)), "ns"))
)

# Cold_low vs Cold_high computed manually
temp_data2 <- wnt_IMvigor_prognosis[wnt_IMvigor_prognosis$level %in% c("Cold_low", "Cold_high"), ]
temp_data2$level <- droplevels(temp_data2$level)

temp_fit2 <- survdiff(Surv(os, censOS) ~ level, data = temp_data2)
cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

mycol <- c('#ffd400', '#009ad6', "grey50", '#d71345')

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Hot", "Cold_low", "Cold_medium", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS (RF, n = ", nrow(wnt_IMvigor_prognosis), ")"),
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

p.lab <- paste0("log-rank test p",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))

p$plot <- p$plot +
  annotate("text", x = 0, y = 0.02,
           label = p.lab, size = 6, hjust = 0, fontface = "italic") +
  annotate("text", x = 0, y = 0.10,
           label = paste0("Hot vs. Cold_high: ", hot_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic", color = "red")
print(p)
################################################################################Group proportions under the KNN algorithm####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")

library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

vector <- unique(c("GSE253260","Dai_2025","Prins_GBM_2019_RNAseq","Dai_2025",
                   "GSE253260","PORTER_2025","Zhao_NatMed_2019","Braun_NatMed_2020",
                   "GSE78220","IMvigor"))

# Create an empty data frame to collect the group information of every dataset
all_data <- data.frame()

for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")

  if (file.exists(data_file)) {
    load(data_file)
    dataset_name <- vector[q]

    if (exists(dataset_name)) {
      eset <- get(dataset_name)

      res_knn <- tryCatch({
        tme_classifier(eset = eset, method = "knn", scale = TRUE)[, c(1, 5)]
      }, error = function(e) {
        message(paste0("Dataset ", dataset_name, " failed: ", e$message))
        return(NULL)
      })

      if (!is.null(res_knn)) {
        colnames(res_knn) <- c("ID", "KNN")

        res_knn$group <- ifelse(res_knn$KNN == "IA", "hot", "cold")

        count_data <- as.data.frame(table(res_knn$group))
        colnames(count_data) <- c("Group", "Count")
        count_data$Dataset <- dataset_name
        count_data$Percentage <- count_data$Count / sum(count_data$Count) * 100

        all_data <- rbind(all_data, count_data)
      }
    }
  }
}

# Create the stacked percentage plot
library(ggplot2)
library(dplyr)
library(scales)

# Make Group a factor with hot first
all_data$Group <- factor(all_data$Group, levels = c("hot", "cold"))

all_data$Dataset <- ifelse(all_data$Dataset=="GSE253260","Nicolle_EBioMedicine_2024",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="Dai_2025","Dai_Nature_2025",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="PORTER_2025","Galsky_ClinCancerRes_2025",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="GSE78220","Hugo_Cell_2016",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="Prins_GBM_2019_RNAseq","Cloughesy_NatMed_2019",all_data$Dataset)

hot_order <- all_data %>%
  filter(Group == "hot") %>%
  arrange(desc(Percentage)) %>%
  pull(Dataset)

all_data$Dataset <- factor(all_data$Dataset, levels = hot_order)

# Compute the cumulative percentage of each bar for label placement
label_data <- all_data %>%
  group_by(Dataset) %>%
  arrange(desc(Group)) %>%
  mutate(
    cumulative_perc = cumsum(Percentage),
    label_pos = cumulative_perc - Percentage / 2
  )

label_data$Dataset <- factor(label_data$Dataset, levels = hot_order)

p <- ggplot(all_data, aes(x = Dataset, y = Percentage, fill = Group)) +   # No longer using reorder
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  geom_text(data = label_data,
            aes(x = Dataset, y = label_pos,
                label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_manual(values = c("hot" = "#f7acbc", "cold" = "#90d7ec")) +
  labs(x = "", y = "Percentage (%)", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 13, color = "black"),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(color = "black", size = 13)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100))

print(p)
################################################################################Group proportions under the KNN algorithm####
rm(list=ls())
gc()
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")

library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

vector <- unique(c("GSE253260","Dai_2025","Prins_GBM_2019_RNAseq","Dai_2025",
                   "GSE253260","PORTER_2025","Zhao_NatMed_2019","Braun_NatMed_2020",
                   "GSE78220","IMvigor"))

# Create an empty data frame to collect the group information of every dataset
all_data <- data.frame()

for (q in 1:length(vector)) {
  data_file = paste0(vector[q], ".Rdata")

  if (file.exists(data_file)) {
    load(data_file)
    dataset_name <- vector[q]

    if (exists(dataset_name)) {
      eset <- get(dataset_name)

      res_knn <- tryCatch({
        tme_classifier(eset = eset, method = "knn", scale = TRUE)[, c(1, 5)]
      }, error = function(e) {
        message(paste0("Dataset ", dataset_name, " failed: ", e$message))
        return(NULL)
      })

      if (!is.null(res_knn)) {
        colnames(res_knn) <- c("ID", "KNN")

        res_knn$group <- ifelse(res_knn$KNN == "IA", "hot", "cold")

        count_data <- as.data.frame(table(res_knn$group))
        colnames(count_data) <- c("Group", "Count")
        count_data$Dataset <- dataset_name
        count_data$Percentage <- count_data$Count / sum(count_data$Count) * 100

        all_data <- rbind(all_data, count_data)
      }
    }
  }
}

# Create the stacked percentage plot
library(ggplot2)
library(dplyr)
library(scales)

# Make Group a factor with hot first
all_data$Group <- factor(all_data$Group, levels = c("hot", "cold"))

all_data$Dataset <- ifelse(all_data$Dataset=="GSE253260","Nicolle_EBioMedicine_2024",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="Dai_2025","Dai_Nature_2025",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="PORTER_2025","Galsky_ClinCancerRes_2025",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="GSE78220","Hugo_Cell_2016",all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset=="Prins_GBM_2019_RNAseq","Cloughesy_NatMed_2019",all_data$Dataset)

hot_order <- all_data %>%
  filter(Group == "hot") %>%
  arrange(desc(Percentage)) %>%
  pull(Dataset)

all_data$Dataset <- factor(all_data$Dataset, levels = hot_order)

# Compute the cumulative percentage of each bar for label placement
label_data <- all_data %>%
  group_by(Dataset) %>%
  arrange(desc(Group)) %>%
  mutate(
    cumulative_perc = cumsum(Percentage),
    label_pos = cumulative_perc - Percentage / 2
  )

label_data$Dataset <- factor(label_data$Dataset, levels = hot_order)

p <- ggplot(all_data, aes(x = Dataset, y = Percentage, fill = Group)) +   # No longer using reorder
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  geom_text(data = label_data,
            aes(x = Dataset, y = label_pos,
                label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_manual(values = c("hot" = "#f7acbc", "cold" = "#90d7ec")) +
  labs(x = "", y = "Percentage (%)", fill = "Group") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    legend.text = element_text(size = 13, color = "black"),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    text = element_text(color = "black", size = 13)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100))

print(p)
################################################################################
################################################################################
# 00_prepare.R
# Data loading + IOBR patch + phenotype annotation + stratified 6-fold + xCell for all samples (once)
################################################################################
rm(list=ls())
gc()

library(IOBR)
library(GSVA)
library(BiocParallel)
library(limma)
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)

# ==================== IOBR patch (applied once globally) ====================
raw_old <- getFromNamespace("rawEnrichmentAnalysis", "IOBR")
raw_fixed <- raw_old

body(raw_fixed) <- substitute({
  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))

  expr <- expr[shared.genes, , drop = FALSE]

  if (nrow(expr) < 5000) {
    stop("Not enough genes")
  }

  expr <- apply(expr, 2, rank)

  params <- GSVA::ssgseaParam(
    exprData = as.matrix(expr),
    geneSets = signatures,
    minSize = 1,
    maxSize = Inf,
    alpha = 0.25,
    normalize = FALSE
  )

  scores <- GSVA::gsva(
    params,
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )

  scores <- scores - apply(scores, 1, min)

  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]

  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]

  if (!is.null(file.name)) {
    write.table(
      scores,
      file = file.name,
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
  }

  scores
})

assignInNamespace(
  "rawEnrichmentAnalysis",
  raw_fixed,
  ns = "IOBR"
)
# ==================== End of IOBR patch ====================

setwd("C:/Users/user/Desktop/input")

# ===== Load the data =
load("IMvigor210CoreBiologies.Rdata")
rm(annoData); rm(exprSet)
load("IMvigor.Rdata")
IMvigor <- log2(IMvigor + 1)

# ===== Phenotype annotation (key fix: unify the column name to phenotype) =
pheno <- na.omit(phenoData[,6,drop=F])
colnames(pheno) <- "phenotype"
pheno$phenotype <- ifelse(pheno$phenotype == "inflamed", "hot", "cold")
pheno$phenotype <- factor(pheno$phenotype, levels = c("hot","cold"))
table(pheno$phenotype)

# Keep the sample order consistent
common_samples <- intersect(rownames(pheno), colnames(IMvigor))
pheno   <- pheno[common_samples, , drop=FALSE]
IMvigor <- IMvigor[, common_samples, drop=FALSE]

# ===== Microenvironment components to keep =
keep_cells <- c(
  "CD8+_Tem", "Macrophages_M1", "Macrophages", "Th2_cells", "pDC", "aDC",
  "CD8+_Tcm", "Th1_cells", "Basophils", "NK_cells", "Tregs", "CD4+_Tem",
  "Monocytes", "ImmuneScore", "CD8+_T-cells", "CD4+_memory_T-cells",
  "Macrophages_M2", "Plasma_cells", "DC", "iDC", "Tgd_cells",
  "MicroenvironmentScore", "CD4+_T-cells", "NKT", "cDC", "B-cells",
  "CD8+_naive_T-cells", "CD4+_naive_T-cells", "CD4+_Tcm",
  "naive_B-cells", "Memory_B-cells", "Class-switched_memory_B-cells",
  "pro_B-cells", "Neutrophils", "Mast_cells", "Eosinophils",
  "Megakaryocytes", "Erythrocytes", "Platelets",
  "Fibroblasts", "Endothelial_cells", "Epithelial_cells",
  "StromaScore"
)

# ===== * Run xCell once for all samples *
cat("Running xCell on ALL samples (only once)...\n")
im_xcell_all <- deconvo_tme(eset = IMvigor, method = "xcell")
im_xcell_all <- as.data.frame(im_xcell_all)
rownames(im_xcell_all) <- im_xcell_all[,1]
im_xcell_all <- im_xcell_all[,-1]
colnames(im_xcell_all) <- gsub("_xCell$", "", colnames(im_xcell_all))

# No scaling here: keep raw values and standardise inside the loop with the training mean/sd
keep_in_data <- intersect(keep_cells, colnames(im_xcell_all))
im_xcell_all <- im_xcell_all[, keep_in_data, drop = FALSE]
cat("Retained cell types:", length(keep_in_data), "\n")

# ===== * Stratified 6-fold (hot/cold ratio preserved per fold) *
set.seed(123)
folds <- createFolds(pheno$phenotype, k = 6, list = TRUE, returnTrain = FALSE)

# Check the hot/cold ratio in each fold
for (i in seq_along(folds)) {
  cat("Fold", i, ":\n")
  print(table(pheno$phenotype[folds[[i]]]))
}

# ===== Save =
save(folds, pheno, IMvigor, keep_cells, keep_in_data, im_xcell_all,
     file = "prepared_data.Rdata")
cat("prepared_data.Rdata saved.\n")
################################################################################
# 01_loop_6folds.R
# 6-fold independent loop: split xCell -> limma -> cluster -> ssGSEA -> RF -> predict -> AUC
# xCell was computed once in 00_prepare.R; only subsetting and score transfer happen here
################################################################################
rm(list=ls())
gc()

setwd("C:/Users/user/Desktop/input")
load("prepared_data.Rdata")

library(IOBR)
library(GSVA)
library(BiocParallel)
library(limma)
library(randomForest)
library(caret)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)

# ==================== IOBR patch (applied again so it survives reloads) ====
# Note: deconvo_tme is not called in this loop any more, but the patch is kept for future extensions
raw_old <- getFromNamespace("rawEnrichmentAnalysis", "IOBR")
raw_fixed <- raw_old

body(raw_fixed) <- substitute({
  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))
  expr <- expr[shared.genes, , drop = FALSE]
  if (nrow(expr) < 5000) stop("Not enough genes")
  expr <- apply(expr, 2, rank)
  params <- GSVA::ssgseaParam(
    exprData = as.matrix(expr),
    geneSets = signatures,
    minSize = 1, maxSize = Inf, alpha = 0.25, normalize = FALSE
  )
  scores <- GSVA::gsva(
    params, verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )
  scores <- scores - apply(scores, 1, min)
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]
  if (!is.null(file.name)) {
    write.table(scores, file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  scores
})

assignInNamespace("rawEnrichmentAnalysis", raw_fixed, ns = "IOBR")
# ==================== End of IOBR patch

# ===== One-sided limma function =
limma.one.sided <- function(fit, lower = FALSE){
  se.coef  <- sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total <- fit$df.prior + fit$df.residual
  pt(fit$t, df = df.total, lower.tail = lower)
}

Get_limma <- function(gene, group, alternative_limma = "greater"){
  design <- model.matrix(~0 + factor(group$group))
  colnames(design) <- levels(factor(group$group))
  rownames(design) <- colnames(gene)
  compare <- makeContrasts(hot - cold, levels = design)
  fit <- lmFit(gene, design)
  fit <- contrasts.fit(fit, compare)
  fit <- eBayes(fit)
  if (alternative_limma == "greater"){
    res <- limma.one.sided(fit, lower = FALSE)
    res <- as.data.frame(res); colnames(res) <- "P.Value"
    res$gene_id <- rownames(res)
    tt <- topTable(fit, coef = 1, number = Inf)
    tt$gene_id <- rownames(tt)
    tt <- tt[, c("gene_id","logFC")]
    res <- merge(res, tt, by = "gene_id")
  } else if (alternative_limma == "less"){
    res <- limma.one.sided(fit, lower = TRUE)
    res <- as.data.frame(res); colnames(res) <- "P.Value"
    res$gene_id <- rownames(res)
    tt <- topTable(fit, coef = 1, number = Inf)
    tt$gene_id <- rownames(tt)
    tt <- tt[, c("gene_id","logFC")]
    res <- merge(res, tt, by = "gene_id")
  } else {
    res <- topTable(fit, coef = 1, number = Inf)
    res <- na.omit(res); res$gene_id <- rownames(res)
  }
  return(res)
}

# ===== Per-fold result storage =
n_folds <- length(folds)          # * now = 6
roc_list <- vector("list", n_folds)
auc_vec  <- rep(NA_real_, n_folds)
skip_folds <- c()

# ===== * Training control: 6 folds x 3 repeats =
train_control <- trainControl(
  method = "repeatedcv", number = 6, repeats = 3,   # * inner 6-fold
  classProbs = TRUE, summaryFunction = twoClassSummary,
  selectionFunction = "best", savePredictions = TRUE,
  verboseIter = TRUE                                 # Enable progress output so the run can be followed
)

# ============ Start the 6-fold loop ============
for (k in 1:n_folds) {
  cat("\n========== Fold", k, "==========\n")

  test_idx  <- folds[[k]]
  train_idx <- setdiff(seq_len(ncol(IMvigor)), test_idx)

  expr_train  <- IMvigor[, train_idx, drop = FALSE]
  expr_test   <- IMvigor[, test_idx,  drop = FALSE]
  pheno_train <- pheno[train_idx, , drop = FALSE]
  pheno_test  <- pheno[test_idx,  , drop = FALSE]

  # ---------- Step 2: subset the precomputed xCell scores and transfer the normalisation
  im_xcell_train_raw <- im_xcell_all[train_idx, , drop = FALSE]
  im_xcell_test_raw  <- im_xcell_all[test_idx,  , drop = FALSE]

  mu  <- apply(im_xcell_train_raw, 2, mean)
  sdv <- apply(im_xcell_train_raw, 2, sd)
  sdv[sdv == 0] <- 1

  im_xcell_train <- as.data.frame(scale(im_xcell_train_raw, center = mu, scale = sdv))
  im_xcell_test  <- as.data.frame(scale(im_xcell_test_raw,  center = mu, scale = sdv))

  cat("Retained cell types:", ncol(im_xcell_train), "\n")

  xcell_pheno_train <- merge(pheno_train, im_xcell_train, by = "row.names")
  rownames(xcell_pheno_train) <- xcell_pheno_train[,1]
  xcell_pheno_train <- xcell_pheno_train[,-1]
  colnames(xcell_pheno_train)[1] <- "phenotype"

  # ---------- Step 3: one-sided limma on the microenvironment components
  exp_im <- as.data.frame(t(xcell_pheno_train[, 2:ncol(xcell_pheno_train)]))
  group_im <- data.frame(Tag = rownames(xcell_pheno_train),
                         group = xcell_pheno_train$phenotype)
  resultup_cell <- Get_limma(gene = exp_im, group = group_im,
                             alternative_limma = "greater")
  resultup_cell <- resultup_cell[order(resultup_cell$P.Value), ]

  sig_cells   <- resultup_cell$gene_id[resultup_cell$P.Value < 0.05]
  n_sig_cells <- length(sig_cells)
  cat("Significant cell types in hot:", n_sig_cells, "\n")

  if (n_sig_cells == 0) {
    cat("Fold", k, ": no significant cell types, skip.\n")
    skip_folds <- c(skip_folds, k)
    next
  }

  # ---------- Step 4: one-sided limma genes, clustered into heat signatures
  exp_gene   <- as.data.frame(expr_train)
  group_gene <- data.frame(Tag = colnames(expr_train),
                           group = pheno_train$phenotype)
  resultup_gene <- Get_limma(gene = exp_gene, group = group_gene,
                             alternative_limma = "greater")

  sig_genes <- resultup_gene$gene_id[
    resultup_gene$P.Value < 0.05 & resultup_gene$logFC >= 0.5
  ]
  cat("Significant genes in hot:", length(sig_genes), "\n")

  if (length(sig_genes) < n_sig_cells) {
    cat("Fold", k, ": not enough significant genes, skip.\n")
    skip_folds <- c(skip_folds, k)
    next
  }

  gene_expr_sig <- expr_train[sig_genes, , drop = FALSE]
  gene_expr_sig <- gene_expr_sig[apply(gene_expr_sig, 1, sd) > 0, , drop = FALSE]
  k_cluster <- min(n_sig_cells, nrow(gene_expr_sig))

  hc <- hclust(dist(gene_expr_sig), method = "complete")
  cluster_assign <- cutree(hc, k = k_cluster)

  Hot_cancer_signature_list <- lapply(seq_len(k_cluster), function(i){
    names(cluster_assign)[cluster_assign == i]
  })
  names(Hot_cancer_signature_list) <- paste0("HotSig_", seq_len(k_cluster))
  cat("Number of hot signatures:", length(Hot_cancer_signature_list), "\n")

  # ---------- Step 5: ssGSEA scoring and merging
  params_train <- ssgseaParam(
    exprData = as.matrix(expr_train),
    geneSets = Hot_cancer_signature_list,
    minSize = 1, maxSize = 10000
  )
  scores_train <- gsva(params_train)
  scores_train <- as.data.frame(t(scores_train))

  params_test <- ssgseaParam(
    exprData = as.matrix(expr_test),
    geneSets = Hot_cancer_signature_list,
    minSize = 1, maxSize = 10000
  )
  scores_test <- gsva(params_test)
  scores_test <- as.data.frame(t(scores_test))

  mu_s  <- apply(scores_train, 2, mean)
  sd_s  <- apply(scores_train, 2, sd)
  sd_s[sd_s == 0] <- 1

  scores_train <- as.data.frame(scale(scores_train, center = mu_s, scale = sd_s))
  scores_test  <- as.data.frame(scale(scores_test,  center = mu_s, scale = sd_s))

  merged_train <- merge(im_xcell_train, scores_train, by = "row.names")
  rownames(merged_train) <- merged_train[,1]
  merged_train <- merged_train[,-1]
  xcell_pheno_train2 <- merge(pheno_train, merged_train, by = "row.names")
  rownames(xcell_pheno_train2) <- xcell_pheno_train2[,1]
  xcell_pheno_train2 <- xcell_pheno_train2[,-1]
  colnames(xcell_pheno_train2)[1] <- "phenotype"

  merged_test <- merge(im_xcell_test, scores_test, by = "row.names")
  rownames(merged_test) <- merged_test[,1]
  merged_test <- merged_test[,-1]
  xcell_pheno_test <- merge(pheno_test, merged_test, by = "row.names")
  rownames(xcell_pheno_test) <- xcell_pheno_test[,1]
  xcell_pheno_test <- xcell_pheno_test[,-1]
  colnames(xcell_pheno_test)[1] <- "phenotype"

  # ---------- Step 6: random-forest tuning (mtry)
  proc_train <- xcell_pheno_train2
  proc_train$target <- ifelse(proc_train$phenotype == "hot", 1, 0)
  proc_train$target_factor <- factor(
    ifelse(proc_train$target == 1, "hot", "cold"),
    levels = c("hot","cold")
  )

  feat_cols <- setdiff(colnames(proc_train),
                       c("phenotype","target","target_factor"))
  n_vars <- length(feat_cols)
  mtry_range <- 1:min(20, n_vars)

  set.seed(123)                                       # Seed 123
  rf_tune <- train(
    x = proc_train[, feat_cols, drop = FALSE],
    y = proc_train$target_factor,
    method = "rf",
    trControl = train_control,
    tuneGrid = expand.grid(mtry = mtry_range),
    ntree = 1000, importance = TRUE, metric = "ROC"
  )
  best_mtry <- rf_tune$bestTune$mtry
  cat("Best mtry:", best_mtry, "\n")

  # ---------- Step 7: train the model with the best parameters
  set.seed(123)                                       # Seed 123
  rf_final <- train(
    x = proc_train[, feat_cols, drop = FALSE],
    y = proc_train$target_factor,
    method = "rf",
    trControl = trainControl(
      method = "repeatedcv", number = 6, repeats = 3, # * inner 6-fold
      classProbs = TRUE, summaryFunction = twoClassSummary,
      selectionFunction = "best", savePredictions = TRUE,
      verboseIter = FALSE
    ),
    tuneGrid = data.frame(mtry = best_mtry),
    ntree = 1000, importance = TRUE, metric = "ROC"
  )

  # ---------- Step 8: predict the test folds and plot the AUC
  proc_test <- xcell_pheno_test
  proc_test$target <- ifelse(proc_test$phenotype == "hot", 1, 0)
  proc_test$target_factor <- factor(
    ifelse(proc_test$target == 1, "hot", "cold"),
    levels = c("hot","cold")
  )
  proc_test <- proc_test[, feat_cols, drop = FALSE]

  pred_prob <- predict(rf_final, newdata = proc_test, type = "prob")

  roc_obj <- roc(response  = pheno_test$phenotype,
                 predictor = pred_prob$hot,
                 levels    = c("cold","hot"))
  auc_val <- as.numeric(auc(roc_obj))
  auc_vec[k] <- auc_val
  roc_list[[k]] <- roc_obj

  cat("Fold", k, "AUC =", round(auc_val, 4), "\n")

  roc_df <- data.frame(
    sensitivity = roc_obj$sensitivities,
    specificity = roc_obj$specificities
  )
  p_fold <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "#2E86AB", size = 1.5) +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
    labs(title = paste0("Fold ", k, " ROC (AUC = ", round(auc_val, 3), ")"),
         x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal(base_size = 14) +
    coord_equal()
  ggsave(paste0("ROC_Fold_", k, ".pdf"), p_fold, width = 5, height = 5)

  save(rf_final, roc_obj, Hot_cancer_signature_list,
       file = paste0("RF_Fold_", k, ".Rdata"))
}

# ============ Summary mean AUC plot ============
mean_auc <- mean(auc_vec, na.rm = TRUE)
cat("\nMean AUC across folds:", round(mean_auc, 4), "\n")
cat("Skipped folds:", ifelse(length(skip_folds) == 0, "none",
                             paste(skip_folds, collapse = ", ")), "\n")

roc_all <- lapply(seq_along(roc_list), function(i){
  if (is.null(roc_list[[i]])) return(NULL)
  data.frame(
    sensitivity = roc_list[[i]]$sensitivities,
    specificity = roc_list[[i]]$specificities,
    fold = paste0("Fold_", i)
  )
})
roc_all <- bind_rows(roc_all)

spec_seq <- seq(0, 1, by = 0.01)
sens_mat <- sapply(roc_list, function(r){
  if (is.null(r)) return(rep(NA, length(spec_seq)))
  approx(x = 1 - r$specificities, y = r$sensitivities,
         xout = spec_seq, method = "linear", rule = 2)$y
})
sens_mean <- rowMeans(sens_mat, na.rm = TRUE)
sens_sd   <- apply(sens_mat, 1, sd, na.rm = TRUE)

mean_roc_df <- data.frame(
  fpr   = spec_seq,
  mean_sens = sens_mean,
  lower = pmax(0, sens_mean - sens_sd),
  upper = pmin(1, sens_mean + sens_sd)
)

p_mean <- ggplot() +
  geom_line(data = roc_all,
            aes(x = 1 - specificity, y = sensitivity, group = fold),
            color = "grey70", alpha = 0.5, size = 0.4) +
  geom_ribbon(data = mean_roc_df,
              aes(x = fpr, ymin = lower, ymax = upper),
              fill = "#2E86AB", alpha = 0.2) +
  geom_line(data = mean_roc_df,
            aes(x = fpr, y = mean_sens),
            color = "#2E86AB", size = 1.5) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray50") +
  labs(title = paste0("Mean ROC across 6 folds (Mean AUC = ",
                      round(mean_auc, 3), ")"),
       x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 14) +
  coord_equal() +
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("Mean AUC = ", round(mean_auc, 3)),
           size = 5, color = "#2E86AB", fontface = "bold")

ggsave("ROC_Mean_6folds.pdf", p_mean, width = 6, height = 6)
################################################################################Dai_Nature_2015####
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier with your ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Dai_Nature_2015, method = "knn", scale = TRUE)

# Assume column 5 holds the classification result; convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
Dai_Nature_2015_G = res_knn[, c(1, ncol(res_knn))]
colnames(Dai_Nature_2015_G) = c("Tag", "group")

# Extract cold samples
Dai_Nature_2015_G = Dai_Nature_2015_G[Dai_Nature_2015_G$group == "Cold",]
Dai_Nature_2015 = Dai_Nature_2015[, colnames(Dai_Nature_2015) %in% Dai_Nature_2015_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Dai_Nature_2015, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("Dai_Nature_2015(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE106128, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
GSE106128_G = res_knn[, c(1, ncol(res_knn))]
colnames(GSE106128_G) = c("Tag", "group")

# Extract cold samples
GSE106128_G = GSE106128_G[GSE106128_G$group == "Cold",]
GSE106128 = GSE106128[, colnames(GSE106128) %in% GSE106128_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE106128")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE106128, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("GSE106128(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

GSE106128_os_p=p_all
save(GSE106128_os_p,file="GSE106128_os_p.Rdata")
################################################################----------------Single dataset (Braun_NatMed_2020)####
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Braun_NatMed_2020, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
Braun_NatMed_2020_G = res_knn[, c(1, ncol(res_knn))]
colnames(Braun_NatMed_2020_G) = c("Tag", "group")

# Extract cold samples
Braun_NatMed_2020_G = Braun_NatMed_2020_G[Braun_NatMed_2020_G$group == "Cold",]
Braun_NatMed_2020 = Braun_NatMed_2020[, colnames(Braun_NatMed_2020) %in% Braun_NatMed_2020_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Braun_NatMed_2020")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Braun_NatMed_2020, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("Braun_NatMed_2020(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
GSE253260_G = res_knn[, c(1, ncol(res_knn))]
colnames(GSE253260_G) = c("Tag", "group")

# Extract cold samples
GSE253260_G = GSE253260_G[GSE253260_G$group == "Cold",]
GSE253260 = GSE253260[, colnames(GSE253260) %in% GSE253260_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("GSE253260(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

GSE253260_os_p=p_all
save(GSE253260_os_p,file="GSE253260_os_p.Rdata")
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
GSE253260_G = res_knn[, c(1, ncol(res_knn))]
colnames(GSE253260_G) = c("Tag", "group")

# Extract cold samples
GSE253260_G = GSE253260_G[GSE253260_G$group == "Cold",]
GSE253260 = GSE253260[, colnames(GSE253260) %in% GSE253260_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("GSE253260(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Zhao_NatMed_2019, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
Zhao_NatMed_2019_G = res_knn[, c(1, ncol(res_knn))]
colnames(Zhao_NatMed_2019_G) = c("Tag", "group")

# Extract cold samples
Zhao_NatMed_2019_G = Zhao_NatMed_2019_G[Zhao_NatMed_2019_G$group == "Cold",]
Zhao_NatMed_2019 = Zhao_NatMed_2019[, colnames(Zhao_NatMed_2019) %in% Zhao_NatMed_2019_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Zhao_NatMed_2019")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Zhao_NatMed_2019, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("Zhao_NatMed_2019(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

Zhao_NatMed_2019_os_p=p_all
save(Zhao_NatMed_2019_os_p,file="Zhao_NatMed_2019_os_p.Rdata")
################################################################----------------Single dataset (Tsimberidou_JExpMed_2024)####
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Tsimberidou_JExpMed_2024, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
Tsimberidou_JExpMed_2024_G = res_knn[, c(1, ncol(res_knn))]
colnames(Tsimberidou_JExpMed_2024_G) = c("Tag", "group")

# Extract cold samples
Tsimberidou_JExpMed_2024_G = Tsimberidou_JExpMed_2024_G[Tsimberidou_JExpMed_2024_G$group == "Cold",]
Tsimberidou_JExpMed_2024 = Tsimberidou_JExpMed_2024[, colnames(Tsimberidou_JExpMed_2024) %in% Tsimberidou_JExpMed_2024_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Tsimberidou_JExpMed_2024")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Tsimberidou_JExpMed_2024, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("Tsimberidou_JExpMed_2024(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

Tsimberidou_JExpMed_2024_os_p=p_all
save(Tsimberidou_JExpMed_2024_os_p,file="Tsimberidou_JExpMed_2024_os_p.Rdata")
################################################################----------------Single dataset (IMvigor)####
################################################################----------------Single dataset (IMvigor)####
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = IMvigor, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
IMvigor_G = res_knn[, c(1, ncol(res_knn))]
colnames(IMvigor_G) = c("Tag", "group")

# Extract cold samples
IMvigor_G = IMvigor_G[IMvigor_G$group == "Cold",]
IMvigor = IMvigor[, colnames(IMvigor) %in% IMvigor_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("IMvigor")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = IMvigor, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("IMvigor(os)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

IMvigor_os_p=p_all
save(IMvigor_os_p,file="IMvigor_os_p.Rdata")
################################################################----------------Single dataset (Dai_Nature_2015 - PFS)####
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Dai_Nature_2015, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
Dai_Nature_2015_G = res_knn[, c(1, ncol(res_knn))]
colnames(Dai_Nature_2015_G) = c("Tag", "group")

# Extract cold samples
Dai_Nature_2015_G = Dai_Nature_2015_G[Dai_Nature_2015_G$group == "Cold",]
Dai_Nature_2015 = Dai_Nature_2015[, colnames(Dai_Nature_2015) %in% Dai_Nature_2015_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("Dai_Nature_2015")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = Dai_Nature_2015, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("Dai_Nature_2015(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

Dai_Nature_2015_pfs_p=p_all
save(Dai_Nature_2015_pfs_p,file="Dai_Nature_2015_pfs_p.Rdata")
################################################################----------------Single dataset (GSE253260 - PFS)####
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
  load(data_file)
  load(info_file)
}

# ==================== Change 1: use TMEclassifier and apply the ifelse conversion ====================
library(TMEclassifier)
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)

# Convert to Hot / Cold with your ifelse rule
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
GSE253260_G = res_knn[, c(1, ncol(res_knn))]
colnames(GSE253260_G) = c("Tag", "group")

# Extract cold samples
GSE253260_G = GSE253260_G[GSE253260_G$group == "Cold",]
GSE253260 = GSE253260[, colnames(GSE253260) %in% GSE253260_G$Tag]
# =========================================================================================

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
library(TMEclassifier)
setwd("C:/Users/user/Desktop/input/Wnt_immunotherapy")
vector=c("GSE253260")
for (l in 1:length(vector)) {
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],".Rdata"))
  load(paste0("C:/Users/user/Desktop/input/Wnt_immunotherapy/",vector[l],"_info.Rdata"))
}

# ==================== Change 2: apply the same ifelse conversion to the second block ====================
source("C:/Users/user/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
res_knn = tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)
res_knn$group = ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

group = res_knn[, c(1, ncol(res_knn))]
colnames(group) = c("Tag", "group")
group = group[group$group == "Cold",]
# =========================================================================================

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
  ggtitle("GSE253260(pfs)") +
  geom_hline(yintercept = y_intercept_value, linetype = "dashed", color = "red", size = 0.6) +
  geom_text(aes(x = Inf, y = y_intercept_value, label = round(y_intercept_value, 6)),
            hjust = 1.2, vjust = -0.5, color = "#f15a22", size = 5)

GSE253260_pfs_p=p_all
save(GSE253260_pfs_p,file="GSE253260_pfs_p.Rdata")
