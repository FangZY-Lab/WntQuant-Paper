# ==============================================================================
# WntQuant revision -- Revised immunotherapy analyses
# ------------------------------------------------------------------------------
# Purpose:    Re-run the immunotherapy response analyses with the revised
#             Hot/Cold assignment and add per-cohort KNN + surv_cutpoint
#             evaluation.
# Inputs:     ICI cohort objects and the revised KNN labels.
# Outputs:    Prognostic models, Kaplan-Meier panels and accuracy summaries.
# ==============================================================================

################################################################################Choosing the best method####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")
rm(annoData)
rm(exprSet)
load("C:/Users/赵定康/Desktop/input/IMvigor.Rdata")
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
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
out_dir <- "C:/Users/赵定康/Desktop"

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
################################################################################Prognostic analysis: high Wnt activity predicts a better outcome####
rm(list=ls())
gc()

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score, phenoData, by = "row.names")
wnt_IMvigor <- wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]

# ========== 2. Split into Hot / Cold
wnt_IMvigor_hot  <- wnt_IMvigor[wnt_IMvigor$`Immune phenotype` %in% c("inflamed"), ]
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$`Immune phenotype` %in% c("desert", "excluded"), ]

wnt_IMvigor_hot_prognosis  <- wnt_IMvigor_hot[,  c("activity_score", "os", "censOS")]
wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os", "censOS")]

# ========== 3. Find the optimal cut point in the Cold group with surv_cutpoint
res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os",
  event     = "censOS",
  variables = c("activity_score")
)

summary(res.cut)
print(res.cut$activity_score$cutpoint)

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low",  "Cold_low",  "Cold_high"
)

# ========== 4. Label the Hot group
wnt_IMvigor_hot_prognosis$level <- "Hot"

# ========== 5. Merge
wnt_IMvigor_prognosis <- rbind(wnt_IMvigor_hot_prognosis, wnt_IMvigor_cold_prognosis)
wnt_IMvigor_prognosis$level <- factor(
  wnt_IMvigor_prognosis$level,
  levels = c("Hot", "Cold_low", "Cold_high")
)

wnt_IMvigor_prognosis <- na.omit(wnt_IMvigor_prognosis)
print(table(wnt_IMvigor_prognosis$level))

# ========== 6. Survival analysis
fitd  <- survdiff(Surv(os, censOS) ~ level, data = wnt_IMvigor_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_prognosis,
  p.adjust.method = "none"
)

# ========== 7. Hot vs Cold_high computed manually
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

# ========== 8. Cold_low vs Cold_high computed manually
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

# ========== 9. Plot
mycol <- c('#ffd400', '#009ad6', '#d71345')
names(fit$strata) <- gsub("level=", "", names(fit$strata))

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
  legend.labs = c("Hot", "Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS (Versus hot tumors, n=", nrow(wnt_IMvigor_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic") +
  annotate("text", x = 0, y = 0.00,
           label = paste0("Hot vs. Cold_high: ", hot_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic", color = "red") +
  annotate("text", x = 0, y = 0.08,
           label = paste0("Cold_low vs. Cold_high: ", cold_low_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic", color = "red")

print(p)
################################################################################Prognostic analysis of the desert group####
rm(list=ls())
gc()

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score, phenoData, by = "row.names")
wnt_IMvigor <- wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]

# ========== 2. Keep only the desert group
wnt_IMvigor_desert <- wnt_IMvigor[wnt_IMvigor$`Immune phenotype` %in% c("desert"), ]

wnt_IMvigor_desert_prognosis <- wnt_IMvigor_desert[, c("activity_score", "os", "censOS")]

cat("Desert 组样本量：", nrow(wnt_IMvigor_desert_prognosis), "\n")

# ========== 3. Find the optimal cut point with surv_cutpoint
res.cut <- surv_cutpoint(
  wnt_IMvigor_desert_prognosis,
  time      = "os",
  event     = "censOS",
  variables = c("activity_score")
)

summary(res.cut)
cat("最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_desert_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Desert_low", "Desert_high"
)

wnt_IMvigor_desert_prognosis$level <- factor(
  wnt_IMvigor_desert_prognosis$level,
  levels = c("Desert_low", "Desert_high")
)

wnt_IMvigor_desert_prognosis <- na.omit(wnt_IMvigor_desert_prognosis)

# Show the number of samples per group
print(table(wnt_IMvigor_desert_prognosis$level))

# ========== 4. Survival analysis
fitd  <- survdiff(Surv(os, censOS) ~ level, data = wnt_IMvigor_desert_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_desert_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_desert_prognosis,
  p.adjust.method = "none"
)

# ========== 5. Extract the p-value of the two-group comparison
desert_low_vs_high_p <- NA

if ("Desert_low" %in% rownames(ps$p.value) && "Desert_high" %in% colnames(ps$p.value)) {
  desert_low_vs_high_p <- ps$p.value["Desert_low", "Desert_high"]
} else if ("Desert_high" %in% rownames(ps$p.value) && "Desert_low" %in% colnames(ps$p.value)) {
  desert_low_vs_high_p <- ps$p.value["Desert_high", "Desert_low"]
} else {
  temp_data <- wnt_IMvigor_desert_prognosis
  temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
  desert_low_vs_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)
}

desert_low_vs_high_label <- ifelse(
  is.na(desert_low_vs_high_p), "NA",
  ifelse(desert_low_vs_high_p < 0.001, "p < 0.001",
         ifelse(desert_low_vs_high_p < 0.05,
                paste0("p = ", round(desert_low_vs_high_p, 3)), "ns"))
)

# ========== 6. Plot (same style as above)
mycol <- c('#009ad6', '#d71345')   # Desert_low blue, Desert_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_desert_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Desert_low", "Desert_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS in desert tumors (n=",
                 nrow(wnt_IMvigor_desert_prognosis), ")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Prognostic analysis of the excluded group####
rm(list=ls())
gc()

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score, phenoData, by = "row.names")
wnt_IMvigor <- wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]

# ========== 2. Keep only the excluded group
wnt_IMvigor_excluded <- wnt_IMvigor[wnt_IMvigor$`Immune phenotype` %in% c("excluded"), ]

wnt_IMvigor_excluded_prognosis <- wnt_IMvigor_excluded[, c("activity_score", "os", "censOS")]

cat("Excluded 组样本量：", nrow(wnt_IMvigor_excluded_prognosis), "\n")

# ========== 3. Find the optimal cut point with surv_cutpoint
res.cut <- surv_cutpoint(
  wnt_IMvigor_excluded_prognosis,
  time      = "os",
  event     = "censOS",
  variables = c("activity_score")
)

summary(res.cut)
cat("最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_excluded_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Excluded_low", "Excluded_high"
)

wnt_IMvigor_excluded_prognosis$level <- factor(
  wnt_IMvigor_excluded_prognosis$level,
  levels = c("Excluded_low", "Excluded_high")
)

wnt_IMvigor_excluded_prognosis <- na.omit(wnt_IMvigor_excluded_prognosis)

# Show the number of samples per group
print(table(wnt_IMvigor_excluded_prognosis$level))

# ========== 4. Survival analysis
fitd  <- survdiff(Surv(os, censOS) ~ level, data = wnt_IMvigor_excluded_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_excluded_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_excluded_prognosis,
  p.adjust.method = "none"
)

# ========== 5. Extract the p-value of the two-group comparison
excluded_low_vs_high_p <- NA

if ("Excluded_low" %in% rownames(ps$p.value) && "Excluded_high" %in% colnames(ps$p.value)) {
  excluded_low_vs_high_p <- ps$p.value["Excluded_low", "Excluded_high"]
} else if ("Excluded_high" %in% rownames(ps$p.value) && "Excluded_low" %in% colnames(ps$p.value)) {
  excluded_low_vs_high_p <- ps$p.value["Excluded_high", "Excluded_low"]
} else {
  temp_data <- wnt_IMvigor_excluded_prognosis
  temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
  excluded_low_vs_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)
}

excluded_low_vs_high_label <- ifelse(
  is.na(excluded_low_vs_high_p), "NA",
  ifelse(excluded_low_vs_high_p < 0.001, "p < 0.001",
         ifelse(excluded_low_vs_high_p < 0.05,
                paste0("p = ", round(excluded_low_vs_high_p, 3)), "ns"))
)

# ========== 6. Plot (same style as above)
mycol <- c('#009ad6', '#d71345')   # Excluded_low blue, Excluded_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_excluded_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Excluded_low", "Excluded_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS in excluded tumors (n=",
                 nrow(wnt_IMvigor_excluded_prognosis), ")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Prognostic analysis of the inflamed group####
rm(list=ls())
gc()

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score, phenoData, by = "row.names")
wnt_IMvigor <- wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]

# ========== 2. Keep only the inflamed group
wnt_IMvigor_inflamed <- wnt_IMvigor[wnt_IMvigor$`Immune phenotype` %in% c("inflamed"), ]

wnt_IMvigor_inflamed_prognosis <- wnt_IMvigor_inflamed[, c("activity_score", "os", "censOS")]

cat("Inflamed 组样本量：", nrow(wnt_IMvigor_inflamed_prognosis), "\n")

# ========== 3. Find the optimal cut point with surv_cutpoint
res.cut <- surv_cutpoint(
  wnt_IMvigor_inflamed_prognosis,
  time      = "os",
  event     = "censOS",
  variables = c("activity_score")
)

summary(res.cut)
cat("最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_inflamed_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Inflamed_low", "Inflamed_high"
)

wnt_IMvigor_inflamed_prognosis$level <- factor(
  wnt_IMvigor_inflamed_prognosis$level,
  levels = c("Inflamed_low", "Inflamed_high")
)

wnt_IMvigor_inflamed_prognosis <- na.omit(wnt_IMvigor_inflamed_prognosis)

# Show the number of samples per group
print(table(wnt_IMvigor_inflamed_prognosis$level))

# ========== 4. Survival analysis
fitd  <- survdiff(Surv(os, censOS) ~ level, data = wnt_IMvigor_inflamed_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_inflamed_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_inflamed_prognosis,
  p.adjust.method = "none"
)

# ========== 5. Extract the p-value of the two-group comparison
inflamed_low_vs_high_p <- NA

if ("Inflamed_low" %in% rownames(ps$p.value) && "Inflamed_high" %in% colnames(ps$p.value)) {
  inflamed_low_vs_high_p <- ps$p.value["Inflamed_low", "Inflamed_high"]
} else if ("Inflamed_high" %in% rownames(ps$p.value) && "Inflamed_low" %in% colnames(ps$p.value)) {
  inflamed_low_vs_high_p <- ps$p.value["Inflamed_high", "Inflamed_low"]
} else {
  temp_data <- wnt_IMvigor_inflamed_prognosis
  temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
  inflamed_low_vs_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)
}

inflamed_low_vs_high_label <- ifelse(
  is.na(inflamed_low_vs_high_p), "NA",
  ifelse(inflamed_low_vs_high_p < 0.001, "p < 0.001",
         ifelse(inflamed_low_vs_high_p < 0.05,
                paste0("p = ", round(inflamed_low_vs_high_p, 3)), "ns"))
)

# ========== 6. Plot (same style as above)
mycol <- c('#009ad6', '#d71345')   # Inflamed_low blue, Inflamed_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_inflamed_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Inflamed_low", "Inflamed_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS in inflamed tumors (n=",
                 nrow(wnt_IMvigor_inflamed_prognosis), ")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Prognostic analysis of all tumours####
rm(list=ls())
gc()

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score, phenoData, by = "row.names")
wnt_IMvigor <- wnt_IMvigor[!is.na(wnt_IMvigor$`Immune phenotype`), ]

# ========== 2. All tumours (no immune-phenotype split)
wnt_IMvigor_all <- wnt_IMvigor

wnt_IMvigor_all_prognosis <- wnt_IMvigor_all[, c("activity_score", "os", "censOS")]

cat("全部肿瘤样本量：", nrow(wnt_IMvigor_all_prognosis), "\n")

# ========== 3. Find the optimal cut point with surv_cutpoint
res.cut <- surv_cutpoint(
  wnt_IMvigor_all_prognosis,
  time      = "os",
  event     = "censOS",
  variables = c("activity_score")
)

summary(res.cut)
cat("最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_all_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "All_low", "All_high"
)

wnt_IMvigor_all_prognosis$level <- factor(
  wnt_IMvigor_all_prognosis$level,
  levels = c("All_low", "All_high")
)

wnt_IMvigor_all_prognosis <- na.omit(wnt_IMvigor_all_prognosis)

# Show the number of samples per group
print(table(wnt_IMvigor_all_prognosis$level))

# ========== 4. Survival analysis
fitd  <- survdiff(Surv(os, censOS) ~ level, data = wnt_IMvigor_all_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_all_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_all_prognosis,
  p.adjust.method = "none"
)

# ========== 5. Extract the p-value of the two-group comparison
all_low_vs_high_p <- NA

if ("All_low" %in% rownames(ps$p.value) && "All_high" %in% colnames(ps$p.value)) {
  all_low_vs_high_p <- ps$p.value["All_low", "All_high"]
} else if ("All_high" %in% rownames(ps$p.value) && "All_low" %in% colnames(ps$p.value)) {
  all_low_vs_high_p <- ps$p.value["All_high", "All_low"]
} else {
  temp_data <- wnt_IMvigor_all_prognosis
  temp_fit <- survdiff(Surv(os, censOS) ~ level, data = temp_data)
  all_low_vs_high_p <- 1 - pchisq(temp_fit$chisq, length(temp_fit$n) - 1)
}

all_low_vs_high_label <- ifelse(
  is.na(all_low_vs_high_p), "NA",
  ifelse(all_low_vs_high_p < 0.001, "p < 0.001",
         ifelse(all_low_vs_high_p < 0.05,
                paste0("p = ", round(all_low_vs_high_p, 3)), "ns"))
)

# ========== 6. Plot (same style as above)
mycol <- c('#009ad6', '#d71345')   # All_low blue, All_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_all_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 25),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("All_low", "All_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS in all tumors (n=",
                 nrow(wnt_IMvigor_all_prognosis), ")"),
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
  annotate("text", x = 0, y = 0,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Dai_Nature_2025_os), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("Dai_2025")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
Dai_2025 <- Dai_2025

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = Dai_2025, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = Dai_2025,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, Dai_2025_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os.time", "os.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$os.time <- wnt_IMvigor_cold_prognosis$os.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os.time",
  event     = "os.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(os.time, os.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(os.time, os.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 72),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Dai_Nature_2025",
  subtitle = paste0("OS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Dai_Nature_2025_pfs), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("Dai_2025")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
Dai_2025 <- Dai_2025

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = Dai_2025, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = Dai_2025,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, Dai_2025_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "pfs.time", "pfs.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$pfs.time <- wnt_IMvigor_cold_prognosis$pfs.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "pfs.time",
  event     = "pfs.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(pfs.time, pfs.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(pfs.time, pfs.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(pfs.time, pfs.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(pfs.time, pfs.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 72),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Dai_Nature_2025",
  subtitle = paste0("PFS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Nicolle_EBioMedicine_2024_os), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("GSE253260")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
GSE253260 <- GSE253260

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = GSE253260,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, GSE253260_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os.time", "os.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$os.time <- wnt_IMvigor_cold_prognosis$os.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os.time",
  event     = "os.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(os.time, os.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(os.time, os.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 69),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Nicolle_EBioMedicine_2024",
  subtitle = paste0("OS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Nicolle_EBioMedicine_2024_pfs), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("GSE253260")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
GSE253260 <- GSE253260

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = GSE253260, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = GSE253260,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, GSE253260_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "pfs.time", "pfs.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$pfs.time <- wnt_IMvigor_cold_prognosis$pfs.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "pfs.time",
  event     = "pfs.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(pfs.time, pfs.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(pfs.time, pfs.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(pfs.time, pfs.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(pfs.time, pfs.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 69),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Nicolle_EBioMedicine_2024",
  subtitle = paste0("PFS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Braun_NatMed_2020_os), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("Braun_NatMed_2020")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
Braun_NatMed_2020 <- Braun_NatMed_2020

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = Braun_NatMed_2020, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = Braun_NatMed_2020,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, Braun_NatMed_2020_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os.time", "os.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$os.time <- wnt_IMvigor_cold_prognosis$os.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os.time",
  event     = "os.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(os.time, os.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(os.time, os.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 74),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Braun_NatMed_2020",
  subtitle = paste0("OS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Zhao_NatMed_2019_os), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("Zhao_NatMed_2019")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
Zhao_NatMed_2019 <- Zhao_NatMed_2019

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = Zhao_NatMed_2019, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = Zhao_NatMed_2019,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, Zhao_NatMed_2019_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os.time", "os.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$os.time <- wnt_IMvigor_cold_prognosis$os.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os.time",
  event     = "os.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(os.time, os.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(os.time, os.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 57),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Zhao_NatMed_2019",
  subtitle = paste0("OS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Hugo_Cell_2016_os), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("GSE78220")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
GSE78220 <- GSE78220

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = GSE78220, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = GSE78220,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, GSE78220_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os.time", "os.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$os.time <- wnt_IMvigor_cold_prognosis$os.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os.time",
  event     = "os.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(os.time, os.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(os.time, os.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 36),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Hugo_Cell_2016",
  subtitle = paste0("OS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Galsky_ClinCancerRes_2025_os), KNN + surv_cutpoint (Cold only)####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("PORTER_2025")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
PORTER_2025 <- PORTER_2025

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = PORTER_2025, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = PORTER_2025,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, PORTER_2025_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os.time", "os.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$os.time <- wnt_IMvigor_cold_prognosis$os.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os.time",
  event     = "os.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(os.time, os.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os.time, os.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(os.time, os.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 36),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Galsky_ClinCancerRes_2025",
  subtitle = paste0("OS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Single dataset (Cloughesy_NatMed_2019_pfs), KNN + surv_cutpoint (Cold only)#####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
load("pfs_test_result_all.Rdata")
load("os_test_result_all.Rdata")
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

vector = c("Prins_GBM_2019_RNAseq")
for (q in 1:length(vector)) {
  data_file  = paste0(vector[q], ".Rdata")
  info_file  = paste0(vector[q], "_info.Rdata")
  load(data_file)
  load(info_file)
}

# ========== 1. All samples, no cold filtering
Prins_GBM_2019_RNAseq <- Prins_GBM_2019_RNAseq

# ========== 2. KNN classification with TMEclassifier
library(TMEclassifier)
source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

res_knn <- tme_classifier(eset = Prins_GBM_2019_RNAseq, method = "knn", scale = TRUE)

# Convert to Hot / Cold with the ifelse rule
res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")
table(res_knn$group)

# ========== 3. Compute the Wnt activity score
library(GSEABase)
geneSets = getGmt("C:/Users/赵定康/Desktop/input/wpags_wpigs.gmt")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

WNT_Score = Calculate_bioactivity_scores(
  file_paths = "C:/Users/赵定康/Desktop/input",
  expression_profile = Prins_GBM_2019_RNAseq,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                        rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
WNT_score = WNT_Score$final_activity_score

# ========== 4. Merge the WNT score, clinical data and KNN labels
wnt_IMvigor <- merge(WNT_score, Prins_GBM_2019_RNAseq_info, by.x = "row.names", by.y = "tag")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

table(wnt_IMvigor$group)

# ========== 5. Keep only the Cold group (Hot is excluded downstream)
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$group == "Cold", ]

wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "pfs.time", "pfs.status")]

# Express time in months
wnt_IMvigor_cold_prognosis$pfs.time <- wnt_IMvigor_cold_prognosis$pfs.time / 30

cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 6. Find the optimal cut point in the Cold group with surv_cutpoint
library(survival)
library(survminer)

res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "pfs.time",
  event     = "pfs.status",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

wnt_IMvigor_cold_prognosis$level <- factor(
  wnt_IMvigor_cold_prognosis$level,
  levels = c("Cold_low", "Cold_high")
)

wnt_IMvigor_cold_prognosis <- na.omit(wnt_IMvigor_cold_prognosis)
print(table(wnt_IMvigor_cold_prognosis$level))

# ========== 7. Survival analysis
fitd  <- survdiff(Surv(pfs.time, pfs.status) ~ level, data = wnt_IMvigor_cold_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(pfs.time, pfs.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(pfs.time, pfs.status) ~ level,
  data = wnt_IMvigor_cold_prognosis,
  p.adjust.method = "none"
)

# ========== 8. Cold_low vs Cold_high computed manually
cold_low_vs_cold_high_p <- NA

if ("Cold_low" %in% rownames(ps$p.value) && "Cold_high" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_low", "Cold_high"]
} else if ("Cold_high" %in% rownames(ps$p.value) && "Cold_low" %in% colnames(ps$p.value)) {
  cold_low_vs_cold_high_p <- ps$p.value["Cold_high", "Cold_low"]
} else {
  temp_data2 <- wnt_IMvigor_cold_prognosis
  temp_fit2 <- survdiff(Surv(pfs.time, pfs.status) ~ level, data = temp_data2)
  cold_low_vs_cold_high_p <- 1 - pchisq(temp_fit2$chisq, length(temp_fit2$n) - 1)
}

cold_low_vs_cold_high_label <- ifelse(
  is.na(cold_low_vs_cold_high_p), "NA",
  ifelse(cold_low_vs_cold_high_p < 0.001, "p < 0.001",
         ifelse(cold_low_vs_cold_high_p < 0.05,
                paste0("p = ", round(cold_low_vs_cold_high_p, 3)), "ns"))
)

# ========== 9. Plot only Cold_low and Cold_high
library(RColorBrewer)
library(tibble)
library(ggpp)

mycol <- c('#009ad6', '#d71345')   # Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  data = wnt_IMvigor_cold_prognosis,
  conf.int = FALSE,
  risk.table = FALSE,
  palette = mycol,
  size = 0.5,
  xlim = c(0, 15),
  break.time.by = 12,
  legend.title = "",
  legend.labs = c("Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = "Cloughesy_NatMed_2019",
  subtitle = paste0("PFS (n=", nrow(wnt_IMvigor_cold_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic")

print(p)
################################################################################Prognostic analysis: high Wnt predicts a better outcome (KNN + surv_cutpoint)####
rm(list=ls())
gc()

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")
load("C:/Users/赵定康/Desktop/input/res_knn.Rdata")
table(res_knn$KNN)

# ========== 1. Merge the WNT score, phenoData and KNN labels
wnt_IMvigor <- merge(IMvigor_WNT_score$final_activity_score,
                     phenoData,
                     by = "row.names")

wnt_IMvigor <- merge(wnt_IMvigor, res_knn,
                     by.x = "Row.names", by.y = "ID")

# IA -> Hot, everything else -> Cold
wnt_IMvigor$KNN <- ifelse(wnt_IMvigor$KNN == "IA", "Hot", "Cold")
table(wnt_IMvigor$KNN)

# ========== 2. Split into Hot / Cold
wnt_IMvigor_hot  <- wnt_IMvigor[wnt_IMvigor$KNN == "Hot", ]
wnt_IMvigor_cold <- wnt_IMvigor[wnt_IMvigor$KNN == "Cold", ]

wnt_IMvigor_hot_prognosis  <- wnt_IMvigor_hot[,  c("activity_score", "os", "censOS")]
wnt_IMvigor_cold_prognosis <- wnt_IMvigor_cold[, c("activity_score", "os", "censOS")]

cat("Hot 组样本量：", nrow(wnt_IMvigor_hot_prognosis), "\n")
cat("Cold 组样本量：", nrow(wnt_IMvigor_cold_prognosis), "\n")

# ========== 3. Find the optimal cut point in the Cold group with surv_cutpoint
res.cut <- surv_cutpoint(
  wnt_IMvigor_cold_prognosis,
  time      = "os",
  event     = "censOS",
  variables = c("activity_score")
)

summary(res.cut)
cat("Cold 组最佳截断点：", res.cut$activity_score$cutpoint, "\n")

res.cat <- surv_categorize(res.cut)
table(res.cat$activity_score)

wnt_IMvigor_cold_prognosis$level <- ifelse(
  res.cat$activity_score == "low", "Cold_low", "Cold_high"
)

# ========== 4. Label the Hot group
wnt_IMvigor_hot_prognosis$level <- "Hot"

# ========== 5. Merge
wnt_IMvigor_prognosis <- rbind(wnt_IMvigor_hot_prognosis, wnt_IMvigor_cold_prognosis)

wnt_IMvigor_prognosis$level <- factor(
  wnt_IMvigor_prognosis$level,
  levels = c("Hot", "Cold_low", "Cold_high")
)

wnt_IMvigor_prognosis <- na.omit(wnt_IMvigor_prognosis)
print(table(wnt_IMvigor_prognosis$level))

# ========== 6. Survival analysis
fitd  <- survdiff(Surv(os, censOS) ~ level, data = wnt_IMvigor_prognosis)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_prognosis,
  type = "kaplan-meier",
  error = "greenwood",
  conf.type = "plain"
)

ps <- pairwise_survdiff(
  Surv(os, censOS) ~ level,
  data = wnt_IMvigor_prognosis,
  p.adjust.method = "none"
)

# ========== 7. Hot vs Cold_high computed manually
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

# ========== 8. Cold_low vs Cold_high computed manually
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

# ========== 9. Plot
mycol <- c('#ffd400', '#009ad6', '#d71345')   # Hot yellow, Cold_low blue, Cold_high red

names(fit$strata) <- gsub("level=", "", names(fit$strata))

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
  legend.labs = c("Hot", "Cold_low", "Cold_high"),
  xlab = "Time (months)",
  ylab = "Probability",
  title = paste0("OS (KNN, n = ", nrow(wnt_IMvigor_prognosis), ")"),
  ylim = c(-0.1, 1),
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
  annotate("text", x = 0, y = -0.08,
           label = p.lab, size = 6, hjust = 0, fontface = "italic") +
  annotate("text", x = 0, y = 0.0,
           label = paste0("Hot vs. Cold_high: ", hot_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic", color = "red") +
  annotate("text", x = 0, y = 0.08,
           label = paste0("Cold_low vs. Cold_high: ", cold_low_vs_cold_high_label),
           size = 6, hjust = 0, fontface = "italic", color = "red")
p
################################################################################cold/hot proportions####
rm(list = ls())
gc()
setwd("C:/Users/赵定康/Desktop/input/Wnt免疫治疗")

library(TMEclassifier)
library(ggplot2)
library(dplyr)
library(scales)

source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")

vector <- unique(c("GSE253260", "Dai_2025", "Prins_GBM_2019_RNAseq",
                   "PORTER_2025", "Zhao_NatMed_2019",
                   "Braun_NatMed_2020", "GSE78220",
                   "IMvigor"))

all_data <- data.frame()

for (q in seq_along(vector)) {

  dataset_name <- vector[q]
  data_file    <- paste0(dataset_name, ".Rdata")

  # Load only the expression matrix, not _G.Rdata
  if (!file.exists(data_file)) {
    cat("跳过（文件不存在）：", dataset_name, "\n")
    next
  }

  if (exists("eset_tmp")) rm(eset_tmp)
  load(data_file)

  if (!exists(dataset_name)) {
    cat("跳过（未找到表达矩阵对象）：", dataset_name, "\n")
    cat("  当前环境对象：", paste(ls(), collapse = ", "), "\n")
    next
  }
  eset_tmp <- get(dataset_name)

  # ========== KNN classification (the only source of grouping)
  res_knn <- tryCatch(
    tme_classifier(eset = eset_tmp, method = "knn", scale = TRUE),
    error = function(e) {
      cat("KNN 失败：", dataset_name, "-", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(res_knn)) next

  # IA -> hot, everything else -> cold
  res_knn$group <- ifelse(res_knn[, 5] == "IA", "Hot", "Cold")

  count_data <- as.data.frame(table(res_knn$group))
  colnames(count_data) <- c("Group", "Count")
  count_data$Dataset    <- dataset_name
  count_data$Percentage <- count_data$Count / sum(count_data$Count) * 100

  all_data <- rbind(all_data, count_data)

  cat("完成：", dataset_name,
      " Hot =", sum(res_knn$group == "Hot"),
      " Cold =", sum(res_knn$group == "Cold"), "\n")
}

# ========== Rename the datasets
all_data$Dataset <- ifelse(all_data$Dataset == "GSE253260", "Nicolle_EBioMedicine_2024", all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset == "Dai_2025", "Dai_Nature_2025", all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset == "PORTER_2025", "Galsky_ClinCancerRes_2025", all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset == "GSE78220", "Hugo_Cell_2016", all_data$Dataset)
all_data$Dataset <- ifelse(all_data$Dataset == "Prins_GBM_2019_RNAseq", "Cloughesy_NatMed_2019", all_data$Dataset)

# ========== Sort by the hot proportion, low to high
all_data$Group <- factor(all_data$Group, levels = c("Hot", "Cold"))

hot_order <- all_data %>%
  filter(Group == "Hot") %>%
  arrange(Percentage) %>%
  pull(Dataset)

# Datasets without a hot group are moved to the front
missing_ds <- setdiff(unique(all_data$Dataset), hot_order)
hot_order  <- c(missing_ds, hot_order)

all_data$Dataset <- factor(all_data$Dataset, levels = hot_order)

# ========== Label positions
label_data <- all_data %>%
  group_by(Dataset) %>%
  arrange(desc(Group)) %>%
  mutate(
    cumulative_perc = cumsum(Percentage),
    label_pos = cumulative_perc - Percentage / 2
  )

# ========== Plot
p <- ggplot(all_data, aes(x = Dataset, y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  geom_text(data = label_data,
            aes(x = Dataset, y = label_pos,
                label = ifelse(Percentage > 5, paste0(round(Percentage, 1), "%"), "")),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_manual(values = c("Hot" = "#f7acbc", "Cold" = "#90d7ec")) +
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
################################################################################benchmark####
rm(list = ls())
gc()

library(GSVA)
library(GSEABase)
library(IOBR)
library(TMEclassifier)

input_dir <- "C:/Users/赵定康/Desktop/input"
wnt_dir   <- "C:/Users/赵定康/Desktop/input/Wnt免疫治疗"

source("C:/Users/赵定康/Desktop/input_wnt/TMEclassifier-master/R/tme_classifier_fix.R")
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

geneSets <- getGmt(file.path(input_dir, "wpags_wpigs.gmt"))

# ---------- Fixed metrics
gene <- c("AXIN2", "ZNRF3", "RNF43", "LGR5", "TCF7", "PDCD1", "CD274", "TGFB1")

signature <- c(
  "CD8+_T_effector", "Pan_F_TBRS", "Antigen_processing_machinery",
  "Immune_checkpoint", "T_cell_inflamed_GEP", "IFNG_signature",
  "T_cell_accumulation", "T_cell_exhaustion", "T_cell_regulatory",
  "ICB_resistance", "MDSC", "TAM", "CAF"
)

immunotherapy <- getGmt(file.path(input_dir, "immunotherapy.gmt"))
immunotherapy_gmt <- lapply(immunotherapy, function(x) x@geneIds)
names(immunotherapy_gmt) <- names(immunotherapy)

wnt_usual <- getGmt(file.path(input_dir, "wnt_usual_using.gmt"))
wnt_usual_gmt <- lapply(wnt_usual, function(x) x@geneIds)
names(wnt_usual_gmt) <- names(wnt_usual)

# ---------- Configuration table (10 entries including IMvigor, OS only)
dataset_config <- data.frame(
  var_name = c("Dai_2025", "Dai_2025",
               "GSE253260", "GSE253260",
               "Braun_NatMed_2020",
               "Zhao_NatMed_2019",
               "GSE78220",
               "PORTER_2025",
               "Prins_GBM_2019_RNAseq",
               "IMvigor"),
  title    = c("Dai_Nature_2025", "Dai_Nature_2025",
               "Nicolle_EBioMedicine_2024", "Nicolle_EBioMedicine_2024",
               "Braun_NatMed_2020",
               "Zhao_NatMed_2019",
               "Hugo_Cell_2016",
               "Galsky_ClinCancerRes_2025",
               "Cloughesy_NatMed_2019",
               "IMvigor210"),
  endpoint = c("os", "pfs", "os", "pfs", "os", "os", "os", "os", "pfs", "os"),
  stringsAsFactors = FALSE
)

# ---------- One data frame per dataset
result_list <- list()

for (i in seq_len(nrow(dataset_config))) {

  v     <- dataset_config$var_name[i]
  title <- dataset_config$title[i]
  ep    <- dataset_config$endpoint[i]

  cat("\n==========", title, ep, "==========\n")

  data_file <- file.path(wnt_dir, paste0(v, ".Rdata"))

  if (!file.exists(data_file)) {
    cat("跳过（文件缺失）：", data_file, "\n")
    next
  }

  if (exists("exp_tmp")) rm(exp_tmp)
  load(data_file)

  if (!exists(v)) {
    cat("跳过（未找到对象）：", v, "\n")
    cat("  当前环境对象：", paste(ls(), collapse = ", "), "\n")
    next
  }
  exp_tmp <- get(v)

  # ---------- KNN classification
  res_knn <- tryCatch(
    tme_classifier(eset = exp_tmp, method = "knn", scale = TRUE),
    error = function(e) {
      cat("KNN 失败：", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(res_knn)) next

  knn_col <- if ("KNN" %in% colnames(res_knn)) "KNN" else colnames(res_knn)[5]
  res_knn$group <- ifelse(res_knn[[knn_col]] == "IA", "Hot", "Cold")

  cold_tags <- res_knn$ID[res_knn$group == "Cold"]
  cold_tags <- intersect(cold_tags, colnames(exp_tmp))

  if (length(cold_tags) < 5) {
    cat("Cold 组样本太少，跳过：", length(cold_tags), "\n")
    next
  }

  cat("Cold n =", length(cold_tags), "\n")

  exp_cal <- exp_tmp[, cold_tags, drop = FALSE]

  # ===== 1. Single-gene expression
  gene_use <- intersect(gene, rownames(exp_cal))
  gene_mat <- as.data.frame(t(exp_cal[gene_use, , drop = FALSE]))
  gene_mat$Sample <- rownames(gene_mat)

  # ===== 2. Wnt score
  WNT_Score <- Calculate_bioactivity_scores(
    file_paths         = input_dir,
    expression_profile = exp_cal,
    foundation         = "relative_ssGSEA",
    activation_geneset = NA,
    inhibition_geneset = NA,
    geneSets_gmt       = geneSets,
    min.sz             = 1,
    max.sz             = 10000,
    geneset_direction  = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                           rep("neg", sum(grepl("WPIGS", names(geneSets)))))
  )
  WNT_score <- WNT_Score$final_activity_score
  WNT_score$Sample <- rownames(WNT_score)

  # ===== 3. Immune signatures
  im_ssGSEA <- calculate_sig_score(
    eset            = as.data.frame(exp_cal),
    signature       = immunotherapy_gmt,
    method          = "ssgsea",
    mini_gene_count = 1
  )
  im_ssGSEA <- as.data.frame(im_ssGSEA)
  rownames(im_ssGSEA) <- im_ssGSEA[, 1]
  im_ssGSEA <- im_ssGSEA[, -1]
  im_ssGSEA$Sample <- rownames(im_ssGSEA)

  # ===== 4. wnt_usual_using.gmt
  wnt_usual_score <- calculate_sig_score(
    eset            = as.data.frame(exp_cal),
    signature       = wnt_usual_gmt,
    method          = "ssgsea",
    mini_gene_count = 1
  )
  wnt_usual_score <- as.data.frame(wnt_usual_score)
  rownames(wnt_usual_score) <- wnt_usual_score[, 1]
  wnt_usual_score <- wnt_usual_score[, -1]
  wnt_usual_score$Sample <- rownames(wnt_usual_score)

  # ===== Merge
  df <- Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE),
               list(gene_mat, WNT_score, im_ssGSEA, wnt_usual_score))

  df <- df[, c("Sample", setdiff(colnames(df), "Sample"))]

  result_list[[paste0(title, "_", ep)]] <- df

  cat("完成，指标列数 =", ncol(df) - 1, "\n")
}

# ---------- Save
save(result_list, file = file.path(input_dir, "result_list.Rdata"))

cat("\n=== 完成 ===\n")
cat("数据集个数：", length(result_list), "\n")
print(names(result_list))
################################################################################Analysis of result_list####
rm(list = ls())
gc()

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(patchwork)

# ---------- 1. Paths and configuration
input_dir <- "C:/Users/赵定康/Desktop/input"
wnt_dir   <- "C:/Users/赵定康/Desktop/input/Wnt免疫治疗"

load(file.path(input_dir, "result_list.Rdata"))

surv_config <- data.frame(
  key        = c("Dai_Nature_2025_os", "Dai_Nature_2025_pfs",
                 "Nicolle_EBioMedicine_2024_os", "Nicolle_EBioMedicine_2024_pfs",
                 "Braun_NatMed_2020_os", "Zhao_NatMed_2019_os",
                 "Hugo_Cell_2016_os", "Galsky_ClinCancerRes_2025_os",
                 "Cloughesy_NatMed_2019_pfs", "IMvigor210_os"),
  var_name   = c("Dai_2025", "Dai_2025", "GSE253260", "GSE253260",
                 "Braun_NatMed_2020", "Zhao_NatMed_2019", "GSE78220",
                 "PORTER_2025", "Prins_GBM_2019_RNAseq", "IMvigor"),
  time_col   = c("os.time", "pfs.time", "os.time", "pfs.time",
                 "os.time", "os.time", "os.time", "os.time", "pfs.time", "os.time"),
  status_col = c("os.status", "pfs.status", "os.status", "pfs.status",
                 "os.status", "os.status", "os.status", "os.status", "pfs.status", "os.status"),
  stringsAsFactors = FALSE
)

# ---------- 2. Load the survival data and merge it with the feature matrix
df_list <- list()
for (i in seq_len(nrow(surv_config))) {
  key  <- surv_config$key[i]
  v    <- surv_config$var_name[i]
  if (!key %in% names(result_list)) next

  info_file <- file.path(wnt_dir, paste0(v, "_info.Rdata"))
  if (!file.exists(info_file)) next

  load(info_file)
  info_tmp <- get(paste0(v, "_info"))
  surv_info <- na.omit(data.frame(
    Sample = info_tmp$tag,
    time   = info_tmp[[surv_config$time_col[i]]] / 30,
    status = as.numeric(info_tmp[[surv_config$status_col[i]]])
  ))

  df_all <- na.omit(merge(result_list[[key]], surv_info, by = "Sample"))
  if (nrow(df_all) > 0) df_list[[key]] <- df_all
}

# ---------- 3. Survival statistics loop (cut point + log-rank + Cox)
p_list  <- hr_list <- list()

for (key in names(df_list)) {
  df_all <- df_list[[key]]
  indicator_cols <- setdiff(colnames(df_all), c("Sample", "time", "status"))

  for (ind in indicator_cols) {
    sub_df <- na.omit(df_all[, c(ind, "time", "status")])
    if (nrow(sub_df) < 10 || length(unique(sub_df[[ind]])) < 3) next
    sub_df[[ind]] <- as.numeric(as.character(sub_df[[ind]]))

    # Optimal cut-off
    cut_res <- tryCatch(surv_cutpoint(sub_df, time = "time", event = "status", variables = ind, minprop = 0.1), error = function(e) NULL)
    if (is.null(cut_res)) next
    cp <- tryCatch(cut_res$cutpoint[ind, "cutpoint"], error = function(e) NULL)
    if (is.null(cp) || is.na(cp)) next

    # Group
    sub_df$level <- factor(ifelse(sub_df[[ind]] > cp, "high", "low"), levels = c("low", "high"))
    if (sum(sub_df$level == "low") < 1 || sum(sub_df$level == "high") < 1) next

    # Log-rank p-value
    fitd <- tryCatch(survdiff(Surv(time, status) ~ level, data = sub_df), error = function(e) NULL)
    if (is.null(fitd)) next
    p_val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

    # Cox model
    cox_fit <- tryCatch(coxph(Surv(time, status) ~ level, data = sub_df), error = function(e) NULL)
    if (is.null(cox_fit)) next
    cox_sum <- summary(cox_fit)
    if (nrow(cox_sum$conf.int) < 1) next

    hr     <- cox_sum$conf.int[1, "exp(coef)"]
    hr_low <- cox_sum$conf.int[1, "lower .95"]
    hr_up  <- cox_sum$conf.int[1, "upper .95"]
    if (any(is.na(c(hr, hr_low, hr_up)))) next

    # Collect the results
    p_list[[length(p_list) + 1]] <- data.frame(Dataset = key, Indicator = ind, p = p_val, stringsAsFactors = FALSE)
    hr_list[[length(hr_list) + 1]] <- data.frame(Dataset = key, Indicator = ind, HR = hr, HR_lower = hr_low, HR_upper = hr_up, stringsAsFactors = FALSE)
  }
}

p_df  <- do.call(rbind, p_list)
hr_df <- do.call(rbind, hr_list)

# ---------- 4. Normalise the names and prepare the plotting data
name_mapping <- c(
  "Pan_F_TBRS"          = "Pan-F-TBRS",
  "ICB_resistance"      = "ICB resistance",
  "Immune_checkpoint"   = "Immune checkpoint",
  "IFNG_signature"      = "IFNG signature",
  "T_cell_inflamed_GEP" = "T cell inflamed GEP",
  "T_cell_exhaustion"   = "T cell exhaustion",
  "T_cell_regulatory"   = "T cell regulatory",
  "T_cell_accumulation" = "T cell accumulation",
  "CD8+_T_effector"     = "CD8+ T effector",
  "activity_score"      = "WntQuant",
  "pos_score"           = "WntQuant_activation",
  "neg_score"           = "WntQuant_inhibition"
)

plot_df <- p_df %>%
  mutate(
    neg_log_p = -log10(p),
    Indicator = recode(Indicator, !!!name_mapping),
    # Rename IMvigor210_os to IMvigor_os
    Dataset = ifelse(Dataset == "IMvigor210_os", "IMvigor_os", Dataset)
  )

mean_df <- plot_df %>%
  group_by(Indicator) %>%
  summarise(mean_val = mean(neg_log_p, na.rm = TRUE)) %>%
  arrange(desc(mean_val))

plot_df$Indicator <- factor(plot_df$Indicator, levels = mean_df$Indicator)
mean_df$Indicator <- factor(mean_df$Indicator, levels = mean_df$Indicator)

p_bar <- ggplot(mean_df, aes(x = Indicator, y = mean_val)) +
  geom_col(fill = "#9b95c9", alpha = 0.85, width = 0.7) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black", face = "plain"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, color = "black", face = "plain"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(y = "")

p_hm <- ggplot(plot_df, aes(x = Indicator, y = Dataset, fill = neg_log_p)) +
  geom_tile(color = "black", linewidth = 0.4) +
  scale_fill_gradient2(
    low = "lightblue",
    mid = "white",
    high = "firebrick3",
    midpoint = 1,
    na.value = "grey90",
    name = "-log10(p)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black", face = "plain"),
    axis.text.y = element_text(size = 14, color = "black", face = "plain"),
    axis.title.x = element_text(size = 14, color = "black", face = "plain"),
    axis.title.y = element_text(size = 14, color = "black", face = "plain"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = "", y = "")

# Assemble the figure
combined_plot <- p_bar / p_hm + plot_layout(heights = c(1, 3))
print(combined_plot)
################################################################################Multivariable Cox regression analysis####
rm(list=ls())
gc()

library(survival)
library(survminer)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

data = merge(IMvigor_WNT_score$final_activity_score, phenoData, by="row.names")
data = data[data$`Immune phenotype` %in% c("excluded","desert"), ]
data = data[, c("os","censOS","activity_score","binaryResponse",
                "FMOne mutation burden per MB","Enrollment IC",
                "Sex","Race","Intravesical BCG administered",
                "Baseline ECOG Score","Tobacco Use History","Met Disease Status",
                "Received platinum","Sample collected pre-platinum",
                "Neoantigen burden per MB","TCGA Subtype")]
colnames(data)

colnames(data) <- c(
  "OS", "CensOS", "WntQuant", "BinaryResponse",
  "FMOne_mutation_burden", "Enrollment_IC",
  "Sex", "Race", "BCG", "ECOG", "Tobacco", "Met_status",
  "Platinum_received", "Pre_platinum_samples",
  "Neoantigen_burden", "TCGA_subtype"
)

# ========== 2. Convert the data types
num_vars <- c("OS", "CensOS", "WntQuant",
              "FMOne_mutation_burden", "Neoantigen_burden")
for (v in num_vars) data[[v]] <- as.numeric(data[[v]])

chr_vars <- c("Race", "Pre_platinum_samples")
for (v in chr_vars) data[[v]] <- as.factor(data[[v]])

data$ECOG <- as.factor(data$ECOG)

data$BinaryResponse    <- factor(data$BinaryResponse, levels = c("SD/PD", "CR/PR"))
data$Sex               <- factor(data$Sex, levels = c("F", "M"))
data$BCG               <- factor(data$BCG, levels = c("N", "Y"))
data$Platinum_received <- factor(data$Platinum_received, levels = c("N", "Y"))
data$Met_status        <- factor(data$Met_status)
data$Tobacco           <- factor(data$Tobacco)
data$TCGA_subtype      <- factor(data$TCGA_subtype)
data$Enrollment_IC     <- factor(data$Enrollment_IC, levels = c("IC0", "IC1", "IC2"))

# ========== 3. Multivariable Cox regression (ECOG excluded, Met_status included)

data_cc <- data[complete.cases(data[, c("OS", "CensOS", "WntQuant")]), ]
cat("完整生存数据样本量：", nrow(data_cc), "\n")

# Drop ECOG, add Met_status back
multi_vars <- c("WntQuant", "BinaryResponse", "Sex",
                "BCG", "Met_status", "Platinum_received",
                "FMOne_mutation_burden", "Neoantigen_burden",
                "TCGA_subtype", "Enrollment_IC")

multi_data <- data_cc[complete.cases(data_cc[, multi_vars]), ]
cat("多因素分析样本量：", nrow(multi_data), "\n")

form_multi <- as.formula(
  paste0("Surv(OS, CensOS) ~ ", paste(multi_vars, collapse = " + "))
)
fit_multi <- coxph(form_multi, data = multi_data)
print(summary(fit_multi))

multi_res <- data.frame(
  Variable = rownames(summary(fit_multi)$coefficients),
  HR       = round(summary(fit_multi)$coefficients[, "exp(coef)"], 3),
  lower95  = round(summary(fit_multi)$conf.int[, "lower .95"], 3),
  upper95  = round(summary(fit_multi)$conf.int[, "upper .95"], 3),
  Pvalue   = round(summary(fit_multi)$coefficients[, "Pr(>|z|)"], 4)
)
print(multi_res)

cat("\n===== WntQuant 结果 =====\n")
print(multi_res[grep("WntQuant", multi_res$Variable), ])

# ========== 5. Forest plot (ggforest)
ggforest(
  fit_multi,
  data = multi_data,
  main = "Multivariable cox regression for OS in cold tumors",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.8,
  refLabel = "Reference",
  noDigits = 3
)


# ========== 6. Interaction visualisation: KM by Wnt high/low within each response stratum

library(survival)
library(survminer)

# ---- 6.1 Split Wnt-high / Wnt-low by the median
multi_data$WntGroup <- factor(
  ifelse(multi_data$WntQuant > median(multi_data$WntQuant, na.rm = TRUE),
         "Wnt High", "Wnt Low"),
  levels = c("Wnt Low", "Wnt High")
)


# ---- 6.2 Stratify by response, two panels in one figure

fit_km <- survfit(Surv(OS, CensOS) ~ WntGroup + BinaryResponse,
                  data = multi_data)

ggsurvplot(
  fit_km,
  data         = multi_data,
  facet.by     = "BinaryResponse",
  pval         = TRUE,
  pval.method  = TRUE,
  risk.table   = TRUE,
  conf.int     = FALSE,
  palette      = c("#2C7FB8", "#D7301F"),
  legend.title = "WntQuant",
  legend.labs  = c("Low", "High"),
  xlab         = "Time (months)",
  ylab         = "Overall survival",
  ggtheme      = theme_bw(base_size = 14),

  # ---- Uniform font: 14 pt, black, not bold
  font.title    = c(14, "plain", "black"),   # Main title (if any)
  font.subtitle = c(14, "plain", "black"),   # Subtitle
  font.caption  = c(14, "plain", "black"),   # Figure caption
  font.x        = c(14, "plain", "black"),   # x-axis title
  font.y        = c(14, "plain", "black"),   # y-axis title
  font.xtickslab = c(14, "plain", "black"),  # x-axis ticks
  font.ytickslab = c(14, "plain", "black"),  # y-axis ticks
  font.legend   = c(14, "plain", "black"),   # Legend text
  font.main     = c(14, "plain", "black"),   # Facet title
  font.facet    = c(14, "plain", "black"),   # Facet title (alternative)

  # ---- risk-table font
  tables.theme = theme_cleantable(base_size = 14),
  risk.table.fontsize = 4.5,                 # Numbers inside the risk table (in pt)

  # ---- P-value font
  pval.size = 5                              # P-value font size (default 5, about 14 pt)
)


# ========== 6. WntQuant x BinaryResponse interaction

library(ggplot2)

# ---- 6.1 Centre the WntQuant scores
multi_data$WntQuant_c <- as.numeric(scale(multi_data$WntQuant,
                                          center = TRUE, scale = FALSE))

# ---- 6.2 Cox model with an interaction term (only to obtain the interaction p-value)
covars_no_int <- setdiff(multi_vars, c("WntQuant", "BinaryResponse"))

form_int <- as.formula(
  paste0("Surv(OS, CensOS) ~ WntQuant_c * BinaryResponse + ",
         paste(covars_no_int, collapse = " + "))
)
fit_int <- coxph(form_int, data = multi_data)

int_term <- grep(":", rownames(summary(fit_int)$coefficients), value = TRUE)
int_p    <- round(summary(fit_int)$coefficients[int_term, "Pr(>|z|)"], 3)
int_p
################################################################################Multivariable Cox regression analysis with KNN####
rm(list=ls())
gc()

library(survival)
library(survminer)

# ========== 1. Load data
load("C:/Users/赵定康/Desktop/input/IMvigor_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/IMvigor210CoreBiologies.Rdata")

data = merge(IMvigor_WNT_score$final_activity_score, phenoData, by="row.names")
load("C:/Users/赵定康/Desktop/input/res_knn.Rdata")
res_knn=res_knn[res_knn$KNN%in%c("IS","IE"),]
data = data[data$Row.names%in%res_knn$ID, ]
data = data[, c("os","censOS","activity_score","binaryResponse",
                "FMOne mutation burden per MB","Enrollment IC",
                "Sex","Race","Intravesical BCG administered",
                "Baseline ECOG Score","Tobacco Use History","Met Disease Status",
                "Received platinum","Sample collected pre-platinum",
                "Neoantigen burden per MB","TCGA Subtype")]

colnames(data) <- c(
  "OS", "CensOS", "WntQuant", "BinaryResponse",
  "FMOne_mutation_burden", "Enrollment_IC",
  "Sex", "Race", "BCG", "ECOG", "Tobacco", "Met_status",
  "Platinum_received", "Pre_platinum_samples",
  "Neoantigen_burden", "TCGA_subtype"
)

# ========== 2. Convert the data types
num_vars <- c("OS", "CensOS", "WntQuant",
              "FMOne_mutation_burden", "Neoantigen_burden")
for (v in num_vars) data[[v]] <- as.numeric(data[[v]])

chr_vars <- c("Race", "Pre_platinum_samples")
for (v in chr_vars) data[[v]] <- as.factor(data[[v]])

data$ECOG <- as.factor(data$ECOG)

data$BinaryResponse    <- factor(data$BinaryResponse, levels = c("SD/PD", "CR/PR"))
data$Sex               <- factor(data$Sex, levels = c("F", "M"))
data$BCG               <- factor(data$BCG, levels = c("N", "Y"))
data$Platinum_received <- factor(data$Platinum_received, levels = c("N", "Y"))
data$Met_status        <- factor(data$Met_status)
data$Tobacco           <- factor(data$Tobacco)
data$TCGA_subtype      <- factor(data$TCGA_subtype)
data$Enrollment_IC     <- factor(data$Enrollment_IC, levels = c("IC0", "IC1", "IC2"))

# ========== 3. Multivariable Cox regression (ECOG and Met_status excluded)

data_cc <- data[complete.cases(data[, c("OS", "CensOS", "WntQuant")]), ]
cat("完整生存数据样本量：", nrow(data_cc), "\n")

# Drop ECOG and Met_status
multi_vars <- c("WntQuant", "BinaryResponse", "Sex",
                "BCG", "Platinum_received",
                "FMOne_mutation_burden", "Neoantigen_burden",
                "TCGA_subtype", "Enrollment_IC")

multi_data <- data_cc[complete.cases(data_cc[, multi_vars]), ]
cat("多因素分析样本量：", nrow(multi_data), "\n")

form_multi <- as.formula(
  paste0("Surv(OS, CensOS) ~ ", paste(multi_vars, collapse = " + "))
)
fit_multi <- coxph(form_multi, data = multi_data)
print(summary(fit_multi))

multi_res <- data.frame(
  Variable = rownames(summary(fit_multi)$coefficients),
  HR       = round(summary(fit_multi)$coefficients[, "exp(coef)"], 3),
  lower95  = round(summary(fit_multi)$conf.int[, "lower .95"], 3),
  upper95  = round(summary(fit_multi)$conf.int[, "upper .95"], 3),
  Pvalue   = round(summary(fit_multi)$coefficients[, "Pr(>|z|)"], 4)
)
print(multi_res)

cat("\n===== WntQuant 结果 =====\n")
print(multi_res[grep("WntQuant", multi_res$Variable), ])

# ========== 5. Forest plot (ggforest)
ggforest(
  fit_multi,
  data = multi_data,
  main = "Multivariable Cox Regression for OS",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.8,
  refLabel = "Reference",
  noDigits = 3
)
