# ==============================================================================
# WntQuant revision -- WSS and panWSS subtype definition
# ------------------------------------------------------------------------------
# Purpose:    Define the transferable WSS and panWSS subtypes, test their
#             internal validity, and project them onto TCGA, TARGET and
#             additional cohorts.
# Inputs:     CRC and pan-cancer expression matrices plus Wnt scores under
#             input/.
# Outputs:    Subtype definition objects (.Rdata), proportion and prognostic
#             figures, TCGA/TARGET comparison.
# ==============================================================================

################################################################################WSS: save a transferable subtype definition####
rm(list=ls()); gc()
library(GSEABase)
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/crc_data_all.Rdata")
load("C:/Users/赵定康/Desktop/input/WNT_Score.Rdata")
WNT_score_sig <- WNT_Score$single_score_matrix
load("C:/Users/赵定康/Desktop/input/cluster.Rdata")

library(factoextra)
genesets <- as.matrix(WNT_score_sig)

# ---- 1. Re-run k-means with a fixed seed
set.seed(123)
km_result <- kmeans(genesets, centers = 4, nstart = 1000,
                    algorithm = "Hartigan-Wong", iter.max = 100)

WSS_centers  <- km_result$centers
WSS_features <- colnames(genesets)

# ---- 2. Map the raw cluster numbers to WSS (from the diagnostic cross-tab)
# raw 1 -> WSS1, raw 2 -> WSS3, raw 3 -> WSS2, raw 4 -> WSS4
raw_to_WSS <- c("1"="WSS1", "2"="WSS3", "3"="WSS2", "4"="WSS4")

# ---- 3. Key step: reorder and rename the centres in WSS order
# WSS_centers rownames are WSS1..WSS4, so no mapping is needed downstream
WSS_centers <- WSS_centers[names(raw_to_WSS), , drop = FALSE]
rownames(WSS_centers) <- raw_to_WSS[rownames(WSS_centers)]
cat("WSS_centers 行名：", paste(rownames(WSS_centers), collapse=", "), "\n")

# ---- 4. Training-set normalisation parameters
train_mean <- colMeans(genesets)
train_sd   <- apply(genesets, 2, sd)
train_sd[train_sd == 0] <- 1e-8

# ---- 5. Nearest-centre assignment (returns rownames(centers), i.e. the WSS labels)
assign_by_centroid <- function(x, centers) {
  d <- matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    d[, k] <- sqrt(rowSums(
      (x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2
    ))
  }
  rownames(centers)[apply(d, 1, which.min)]
}

# ---- 6. Training-set consistency rate (should be 100%)
pred_train <- assign_by_centroid(genesets, WSS_centers)
orig_label <- cluster$cluster[match(rownames(genesets), cluster$Tag)]
cat("训练集自洽率：",
    round(mean(pred_train == orig_label, na.rm = TRUE) * 100, 2), "%\n")
print(table(pred = pred_train, orig = orig_label))

# ---- 7. Save
save(WSS_centers, WSS_features, train_mean, train_sd,
     file = "WSS_subtype_definition.Rdata")
cat("已保存：WSS_subtype_definition.Rdata\n")
cat("中心维度：", dim(WSS_centers), "\n")   # 4 x 63
cat("中心行名：", paste(rownames(WSS_centers), collapse=", "), "\n")
################################################################################Validation in additional datasets####
rm(list=ls()); gc()
library(GSVA); library(GSEABase)

dir_input    <- "C:/Users/赵定康/Desktop/input"
dir_exp_data <- "C:/Users/赵定康/Desktop/input/exp_data"   # <- your actual path

stopifnot(dir.exists(dir_input), dir.exists(dir_exp_data))

setwd(dir_input)
load(file.path(dir_input, "WSS_subtype_definition.Rdata"))
# WSS_centers rownames are now WSS1..WSS4, so raw_to_WSS is no longer needed

geneSets <- getGmt(file.path(dir_input, "wpags_wpigs.gmt"))
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

vector <- c("GSE44076","GSE44861","GSE21510",
            "GSE68468","GSE37178","GSE18105",
            "GSE21815","GSE17537","GSE29621",
            "GSE38832","GSE72970","GSE33113",
            "UCSC.CRC_CMS2","GSE39582","GSE2109","GSE13067",
            "GSE13294","GSE14333","GSE17536",
            "GSE20916","GSE23878","GSE35896",
            "GSE37892","KFSYSCC","syn26720761",
            "syn2623706")

# ---- Nearest-centre assignment (returns rownames(centers) = WSS labels)
assign_by_centroid <- function(x, centers) {
  d <- matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    d[, k] <- sqrt(rowSums(
      (x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2
    ))
  }
  rownames(centers)[apply(d, 1, which.min)]
}

subtype_all <- list(); score_all <- list(); wnt_mat_all <- list()

for (i in seq_along(vector)) {
  q <- vector[i]
  cat("\n=====", i, "/", length(vector), ":", q, "=====\n")

  # ---- Read the expression matrix ----
  f <- file.path(dir_exp_data, paste0(q, "_CMS.Rdata"))
  if (!file.exists(f)) { cat("  文件不存在：", f, "\n"); next }
  ok <- tryCatch({ load(f); TRUE },
                 error = function(e) { cat("  载入失败：", conditionMessage(e), "\n"); FALSE })
  if (!ok) next

  rt <- get(paste0(q, "_CMS"))
  if (!is.matrix(rt) && !is.data.frame(rt)) { cat("  对象类型异常，跳过\n"); next }

  # ---- 1. Wnt scoring (same function, same GMT)
  WNT_Score <- tryCatch(
    Calculate_bioactivity_scores(
      file_paths         = dir_input,
      expression_profile = rt,
      foundation         = "relative_ssGSEA",
      activation_geneset = NA,
      inhibition_geneset = NA,
      geneSets_gmt       = geneSets,
      min.sz             = 1,
      max.sz             = 10000,
      geneset_direction  = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                             rep("neg", sum(grepl("WPIGS", names(geneSets)))))
    ),
    error = function(e) { cat("  打分失败：", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(WNT_Score)) next

  sm <- as.matrix(WNT_Score$single_score_matrix)

  # ---- 2. Align features
  missing <- setdiff(WSS_features, colnames(sm))
  if (length(missing) > 0) {
    cat("  缺失特征数:", length(missing), "（补0）\n")
    for (m in missing) sm <- cbind(sm, setNames(0, m))
  }
  sm <- sm[, WSS_features, drop = FALSE]

  # ---- 3. Standardise with the training-set parameters
  sm_z <- sweep(sm, 2, train_mean, "-")
  sm_z <- sweep(sm_z, 2, train_sd, "/")
  cz   <- sweep(WSS_centers, 2, train_mean, "-")
  cz   <- sweep(cz, 2, train_sd, "/")

  # ---- 4. Assign subtypes
  pred <- assign_by_centroid(sm_z, cz)
  subtype_all[[q]] <- data.frame(Tag = rownames(sm_z), subtype = pred,
                                 dataset = q, stringsAsFactors = FALSE)
  score_all[[q]] <- data.frame(
    Tag = rownames(WNT_Score$final_activity_score),
    activity_score = WNT_Score$final_activity_score$activity_score,
    dataset = q, stringsAsFactors = FALSE
  )
  wnt_mat_all[[q]] <- sm_z

  print(table(pred))
}

subtype_df <- do.call(rbind, subtype_all)
score_df   <- do.call(rbind, score_all)

save(subtype_df, score_df, wnt_mat_all,
     file = file.path(dir_input, "all_cohorts_subtype.Rdata"))

cat("\n完成。subtype_df 行数：", nrow(subtype_df), "\n")
cat("成功队列数：", length(subtype_all), " / ", length(vector), "\n")

# ---- 5. Immediately verify baseline self-consistency
load(file.path(dir_input, "cluster.Rdata"))
base <- subtype_df[subtype_df$dataset == "syn2623706", ]
base <- base[match(cluster$Tag, base$Tag), ]
cat("baseline 自洽率：",
    round(mean(base$subtype == cluster$cluster) * 100, 2), "%\n")
print(table(pred = base$subtype, orig = cluster$cluster))
################################################################################Proportion visualisation####
rm(list=ls()); gc()
load("C:/Users/赵定康/Desktop/input/all_cohorts_subtype.Rdata")
library(dplyr); library(tidyr); library(ggplot2)

prop_df <- subtype_df %>%
  group_by(dataset, subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dataset) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_wide <- prop_df %>%
  select(dataset, subtype, prop) %>%
  pivot_wider(names_from = subtype, values_from = prop, values_fill = 0)

base_vec <- as.numeric(prop_wide[prop_wide$dataset == "syn2623706",
                                 c("WSS1","WSS2","WSS3","WSS4")])
prop_wide$cor_to_baseline <- apply(
  prop_wide[, c("WSS1","WSS2","WSS3","WSS4")], 1,
  function(x) cor(as.numeric(x), base_vec)
)
print(prop_wide[order(-prop_wide$cor_to_baseline), c("dataset","cor_to_baseline")])
prop_df$dataset=ifelse(prop_df$dataset=="UCSC.CRC_CMS2","UCSC.CRC",prop_df$dataset)
ggplot(prop_df, aes(dataset, prop, fill = subtype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("#d71345","#45b97c","#145b7d","#ffd400")) +
  theme_minimal() +
  theme(
    text        = element_text(size = 14, color = "black", face = "plain"),
    axis.text   = element_text(size = 14, color = "black", face = "plain"),
    axis.title  = element_text(size = 14, color = "black", face = "plain"),
    plot.title  = element_text(size = 14, color = "black", face = "plain"),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.title= element_text(size = 14, color = "black", face = "plain"),
    strip.text  = element_text(size = 14, color = "black", face = "plain"),
    axis.text.x = element_text(size = 14, color = "black", face = "plain",
                               angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Proportion",
       title = "Subtype composition across 26 CRC cohorts")

################################################################################Spearman ρ####
rm(list=ls()); gc()
library(GSVA); library(GSEABase); library(dplyr)

dir_input    <- "C:/Users/赵定康/Desktop/input"
dir_exp_data <- "C:/Users/赵定康/Desktop/input/exp_data"
load(file.path(dir_input, "all_cohorts_subtype.Rdata"))
geneSets_NM <- getGmt(file.path(dir_input, "NM_genesets.gmt"))

bio_cols <- c("IMMUNE_ESTIMATE","STROMAL_ESTIMATE",
              "SERRATED_UP","EMT_CORE_GENES","MYC_TARGETS_ZELLER","TGFB_CORE_GENES","CSC_BATLLE",
              "MATRIX_REMODEL_REACTOME","WOUND_RESPONSE_GO_BP","CRYPT_BASE","CRYPT_TOP",
              "EPITH_LOBODA","MESENCH_LOBODA",
              "FGF_ACT_CP","NOTCH_REACTOME","NFKB_BIOCARTA","MAPK_KEGG","PI3K_ACT_REACTOME",
              "SRC_ACT_BIOCARTA","JAK_STAT_KEGG","CASPASE_BIOCARTA","PROTEASOME_KEGG",
              "CELL_CYCLE_BIOCARTA","SHH_KEGG","INTEGRIN_BETA3_CP","VEGF_VEGFR_REACTOME",
              "TRANSLATION_RIBOS_REACTOME",
              "IMMUNE_RESP_GO_BP","IMMUNE_MDSC_ALBELDA","IMMUNE_NKC_BREAST","IMMUNE_TH1_GALON",
              "IMMUNE_TH17_GOUNARI","IMMUNE_THF_BREAST","IMMUNE_TREG_GALON","PD1_REACTOME",
              "COMPLEMENT_COAG_KEGG",
              "ARACHNOID_METAB_KEGG","FATTY_ACID_METAB_KEGG","FRUTOSE_MANNOSE_METAB_KEGG",
              "GALACTOSE_METAB_KEGG","GLUTATHIONE_KEGG","GLYCEROPHOSPHOLIPID_METAB_KEGG",
              "LINOLEIC_METAB_KEGG","LYSOPHOSPHOLIPID_PID","NUCLEOTIDE_METAB_REACTOME",
              "PENTOSE_GLUC_METAB_KEGG","STARCH_SUCROSE_METAB_KEGG","TYROSINE_METAB_KEGG",
              "AMINO_SUGAR_NUCLEO_METAB_KEGG","ALANINE_ASPARTATE_GLUTAMATE_KEGG")

# ---- 1. Compute the subtype x feature mean matrix for the 26 cohorts in bulk
bio_by_subtype_all <- data.frame()
for (q in unique(subtype_df$dataset)) {
  cat(">>", q, "\n")
  f <- file.path(dir_exp_data, paste0(q, "_CMS.Rdata"))
  if (!file.exists(f)) next
  load(f); rt <- get(paste0(q, "_CMS"))

  gs <- tryCatch(gsva(as.matrix(rt), geneSets_NM, kcdf="Gaussian",
                      method="ssgsea", min.sz=1, max.sz=10000),
                 error = function(e) NULL)
  if (is.null(gs)) next

  gs <- as.data.frame(scale(t(gs)))              # Within-cohort z-score
  bio_use <- intersect(bio_cols, colnames(gs))
  gs <- gs[, bio_use, drop = FALSE]
  gs$Tag <- rownames(gs)
  gs <- merge(gs, subtype_df[subtype_df$dataset == q, c("Tag","subtype")],
              by = "Tag")

  bio_by_subtype_all <- rbind(
    bio_by_subtype_all,
    gs %>% group_by(subtype) %>%
      summarise(across(all_of(bio_use), mean), .groups="drop") %>%
      mutate(dataset = q)
  )
}
save(bio_by_subtype_all,
     file = file.path(dir_input, "bio_by_subtype_all.Rdata"))

# ---- 2. The only metric: Spearman rho
make_matrix <- function(dt, ds) {
  m <- dt %>% filter(dataset == ds) %>% select(-dataset) %>%
    tibble::column_to_rownames("subtype") %>% as.matrix()
  m[order(rownames(m)), , drop = FALSE]
}
base_mat <- make_matrix(bio_by_subtype_all, "syn2623706")

res <- data.frame()
for (ds in unique(bio_by_subtype_all$dataset)) {
  if (ds == "syn2623706") next
  m <- make_matrix(bio_by_subtype_all, ds)
  cs <- intersect(rownames(base_mat), rownames(m))
  cc <- intersect(colnames(base_mat), colnames(m))
  b <- base_mat[cs, cc, drop=FALSE]; x <- m[cs, cc, drop=FALSE]
  rho <- cor(as.numeric(b), as.numeric(x), method = "spearman")
  res <- rbind(res, data.frame(dataset = ds, spearman_rho = rho))
}
res <- res[order(-res$spearman_rho), ]
print(res)
cat("\n中位 Spearman ρ =", round(median(res$spearman_rho), 3), "\n")
save(res, file = file.path(dir_input, "pattern_spearman.Rdata"))
library(ggplot2)
ggplot(res, aes(reorder(dataset, -spearman_rho), spearman_rho, fill = spearman_rho)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0.65, linetype = 5, color = "red") +
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",
                       midpoint = median(res$spearman_rho),
                       limits = c(min(res$spearman_rho), max(res$spearman_rho)),
                       name = "Spearman ρ") +
  theme_minimal() +
  theme(
    text        = element_text(size = 14, color = "black", face = "plain"),
    axis.text   = element_text(size = 14, color = "black", face = "plain"),
    axis.title  = element_text(size = 14, color = "black", face = "plain"),
    plot.title  = element_text(size = 14, color = "black", face = "plain"),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.title= element_text(size = 14, color = "black", face = "plain"),
    axis.text.x = element_text(size = 14, color = "black", face = "plain",
                               angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Spearman ρ vs. baseline",
       title = "Biological pattern consistency across 26 CRC cohorts")

################################################################################Internal structural validity: between/within group distance ratio + permutation test####
rm(list=ls()); gc()
load("C:/Users/赵定康/Desktop/input/all_cohorts_subtype.Rdata")
library(dplyr)

# Distance-ratio helper
dist_ratio <- function(D, sub) {
  us <- unique(sub)
  # Within groups
  w <- sapply(us, function(s) {
    ix <- which(sub == s)
    if (length(ix) < 2) return(NA)
    mean(D[ix, ix][lower.tri(D[ix, ix])])
  })
  # Between groups
  b <- c()
  for (i in 1:(length(us) - 1)) for (j in (i + 1):length(us)) {
    i1 <- which(sub == us[i]); i2 <- which(sub == us[j])
    if (length(i1) == 0 || length(i2) == 0) next
    b <- c(b, mean(D[i1, i2]))
  }
  mean(b, na.rm = TRUE) / mean(w, na.rm = TRUE)
}

set.seed(123)
n_perm <- 999
ratio_df <- data.frame()

for (q in names(wnt_mat_all)) {
  cat(">>", q, "\n")
  X <- wnt_mat_all[[q]]              # 63 standardised scores, rows are samples
  sub_q <- subtype_df[subtype_df$dataset == q, c("Tag", "subtype")]
  sub <- sub_q$subtype[match(rownames(X), sub_q$Tag)]

  # Requires at least 2 subtypes with at least 2 samples each
  tb <- table(sub)
  if (length(tb) < 2 || any(tb < 2)) { cat("  样本过少，跳过\n"); next }

  D <- as.matrix(dist(X))
  obs <- dist_ratio(D, sub)

  # Permutation: shuffle the subtype labels and recompute the ratio
  perm <- numeric(n_perm)
  for (b in 1:n_perm) perm[b] <- dist_ratio(D, sample(sub))

  p <- (sum(perm >= obs) + 1) / (n_perm + 1)

  ratio_df <- rbind(ratio_df, data.frame(
    dataset = q, ratio = obs, p_value = p, n = nrow(X),
    n_subtypes = length(tb), min_subtype_n = min(tb)
  ))
}

ratio_df <- ratio_df[order(-ratio_df$ratio), ]
print(ratio_df)
################################################################################Prognostic validation in GSE16125####
rm(list=ls()); gc()
library(GSVA); library(GSEABase)
library(survival); library(survminer)

dir_input    <- "C:/Users/赵定康/Desktop/input"
dir_exp_data <- "C:/Users/赵定康/Desktop/prognosis/CRC_data"
dir_prog     <- "C:/Users/赵定康/Desktop/prognosis"

load(file.path(dir_input, "WSS_subtype_definition.Rdata"))
geneSets <- getGmt(file.path(dir_input, "wpags_wpigs.gmt"))
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

# ---- 1. Load the GSE16125 expression matrix and classify
load(file.path(dir_exp_data, "GSE16125_os_exp.Rdata"))
rt <- GSE16125_os_exp

WNT_Score <- Calculate_bioactivity_scores(
  file_paths         = dir_input,
  expression_profile = rt,
  foundation         = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt       = geneSets,
  min.sz             = 1,
  max.sz             = 10000,
  geneset_direction  = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                         rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
sm <- as.matrix(WNT_Score$single_score_matrix)

missing <- setdiff(WSS_features, colnames(sm))
if (length(missing) > 0) for (m in missing) sm <- cbind(sm, setNames(0, m))
sm <- sm[, WSS_features, drop = FALSE]

sm_z <- sweep(sm, 2, train_mean, "-")
sm_z <- sweep(sm_z, 2, train_sd, "/")
cz   <- sweep(WSS_centers, 2, train_mean, "-")
cz   <- sweep(cz, 2, train_sd, "/")

assign_by_centroid <- function(x, centers) {
  d <- matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    d[, k] <- sqrt(rowSums(
      (x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2
    ))
  }
  rownames(centers)[apply(d, 1, which.min)]
}
pred <- assign_by_centroid(sm_z, cz)
print(table(pred))

# ---- 2. Merge the survival data
load(file.path(dir_prog, "GSE16125_os_survival.Rdata"))
sd <- GSE16125_os_survival
sd$Tag <- rownames(sd)
sub_ds <- data.frame(Tag = rownames(sm_z), subtype = pred, stringsAsFactors = FALSE)
df <- merge(sd, sub_ds, by = "Tag")
cat("合并后样本数：", nrow(df), "\n")

df$grp <- factor(ifelse(df$subtype == "WSS4", "WSS4", "Other"),
                 levels = c("Other","WSS4"))

fit1 <- survfit(Surv(time, status) ~ grp, data = df)
sd1  <- survdiff(Surv(time, status) ~ grp, data = df)
p1   <- 1 - pchisq(sd1$chisq, length(sd1$n) - 1)
cox1 <- coxph(Surv(time, status) ~ grp, data = df)
hr1  <- summary(cox1)$conf.int["grpWSS4", "exp(coef)"]
hr1_l<- summary(cox1)$conf.int["grpWSS4", "lower .95"]
hr1_h<- summary(cox1)$conf.int["grpWSS4", "upper .95"]

cat("\n[WSS4 vs Other]\n")
cat("log-rank p =", signif(p1, 3), "\n")
cat("HR =", round(hr1, 3), "(95% CI", round(hr1_l, 3), "-", round(hr1_h, 3), ")\n")

g1 <- ggsurvplot(
  fit1, data = df,
  pval = TRUE, pval.size = 5, pval.coord = c(10, 0.3),
  conf.int = FALSE, risk.table = FALSE, ncensor.plot = FALSE,
  palette = c("#90d7ec", "#d71345"),
  title = paste0("GSE16125 (OS, n=", nrow(df), ")"),
  xlab = "Time (Months)", ylab = "Probability",
  legend.title = "", legend.labs = c("Other", "WSS4"),
  font.title = 14, font.x = 14, font.y = 14, font.tickslab = 12,
  ggtheme = theme_bw(base_size = 12)
)
print(g1)

df$subtype <- factor(df$subtype, levels = c("WSS1","WSS2","WSS3","WSS4"))

fit2 <- survfit(Surv(time, status) ~ subtype, data = df)
sd2  <- survdiff(Surv(time, status) ~ subtype, data = df)
p2   <- 1 - pchisq(sd2$chisq, length(sd2$n) - 1)

cat("\n[WSS1-4]\n")
cat("log-rank p =", signif(p2, 3), "\n")

g2 <- ggsurvplot(
  fit2, data = df,
  pval = TRUE, pval.size = 5, pval.coord = c(10, 0.3),
  conf.int = FALSE, risk.table = FALSE, ncensor.plot = FALSE,
  palette = c("#d71345","#45b97c","#145b7d","#ffd400"),
  title = paste0("GSE16125 (OS, n=", nrow(df), ")"),
  xlab = "Time (Months)", ylab = "Probability",
  legend.title = "Subtype",
  legend.labs = c("WSS1","WSS2","WSS3","WSS4"),
  font.title = 14, font.x = 14, font.y = 14, font.tickslab = 12,
  ggtheme = theme_bw(base_size = 12)
)
print(g2)
################################################################################panWSS: save a transferable subtype definition####
rm(list=ls()); gc()
library(GSEABase); library(factoextra)

setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/pancancer_WNT_Score.Rdata")
load("C:/Users/赵定康/Desktop/input/pancancer_group.Rdata")

genesets <- as.matrix(pancancer_WNT_Score$single_score_matrix)
pancancer_group <- pancancer_group[pancancer_group$Group == "Tumor", ]
genesets <- genesets[rownames(pancancer_group), ]

# ---- 1. k-means (fixed seed)
set.seed(123)
km_result <- kmeans(genesets, centers = 4, nstart = 1000,
                    algorithm = "Hartigan-Wong", iter.max = 100)

# ---- 2. Median activity per raw cluster, used to reorder the clusters
raw_cluster <- km_result$cluster
acts <- sapply(1:4, function(k) {
  median(pancancer_WNT_Score$final_activity_score[
    rownames(genesets)[raw_cluster == k], "activity_score"])
})
rank_df <- data.frame(raw = 1:4, value = acts)
rank_df <- rank_df[order(-rank_df$value), ]
cat("原始 cluster 按 activity 排序：\n"); print(rank_df)

# ---- 3. Reorder centres by activity from high to low; rename the rows panWSS1-4
panWSS_centers <- km_result$centers
panWSS_centers <- panWSS_centers[rank_df$raw, , drop = FALSE]
rownames(panWSS_centers) <- paste0("panWSS", 1:4)
cat("panWSS_centers 行名：",
    paste(rownames(panWSS_centers), collapse=", "), "\n")

# ---- 4. Training-set parameters
panWSS_features <- colnames(genesets)
train_mean_pan  <- colMeans(genesets)
train_sd_pan    <- apply(genesets, 2, sd)
train_sd_pan[train_sd_pan == 0] <- 1e-8

# ---- 5. Training-set consistency check (should be 100%)
assign_by_centroid <- function(x, centers) {
  d <- matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    d[, k] <- sqrt(rowSums(
      (x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2
    ))
  }
  rownames(centers)[apply(d, 1, which.min)]
}

pred_panWSS <- assign_by_centroid(genesets, panWSS_centers)
orig_panWSS <- paste0("panWSS", match(raw_cluster, rank_df$raw))

cat("训练集自洽率：",
    round(mean(pred_panWSS == orig_panWSS) * 100, 2), "%\n")
print(table(pred = pred_panWSS, orig = orig_panWSS))

# ---- 6. New: per-sample panWSS1-4 string label for TCGA
pancancer_labels <- data.frame(
  Tag     = rownames(genesets),
  subtype = orig_panWSS,           # Already a panWSS1-4 string
  stringsAsFactors = FALSE
)
cat("TCGA 样本 panWSS 标签分布：\n")
print(table(pancancer_labels$subtype))

# ---- 7. Save (includes pancancer_labels; section B does not need rank_df)
save(panWSS_centers, panWSS_features,
     train_mean_pan, train_sd_pan,
     pancancer_labels,
     file = "panWSS_subtype_definition.Rdata")
cat("已保存：panWSS_subtype_definition.Rdata\n")
cat("centers 维度：", dim(panWSS_centers), "\n")
cat("centers 行名：",
    paste(rownames(panWSS_centers), collapse=", "), "\n")
cat("pancancer_labels 行数：", nrow(pancancer_labels), "\n")
################################################################################Organising the TARGET data####
rm(list=ls()); gc()
library(data.table); library(limma); library(sva)

setwd("C:/Users/赵定康/Desktop/input")

files <- c("TARGET-AML.star_tpm.tsv","TARGET-ALL-P1.star_tpm.tsv",
           "TARGET-ALL-P2.star_tpm.tsv","TARGET-ALL-P3.star_tpm.tsv",
           "TARGET-CCSK.star_tpm.tsv","TARGET-NBL.star_tpm.tsv",
           "TARGET-OS.star_tpm.tsv","TARGET-RT.star_tpm.tsv",
           "TARGET-WT.star_tpm.tsv")
proj_names <- c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")

dt_list <- lapply(files, function(f) fread(f, header = TRUE, sep = "\t"))

rt <- Reduce(function(x, y) merge(x, y, by = "Ensembl_ID", all = FALSE), dt_list)
cat("合并后维度：", dim(rt), "\n")

ann <- fread("gencode.v36.annotation.gtf.gene.probemap",
             header = TRUE, sep = "\t")
cat("ann 列名：", paste(colnames(ann), collapse=", "), "\n")

ann_dt <- data.table(
  Ensembl_ID = as.character(ann[[1]]),   # Column 1
  gene       = as.character(ann[[2]])    # Column 2
)

# ---- Drop the trailing "." from all version numbers
ann_dt[, Ensembl_ID := sub("\\..*$", "", Ensembl_ID)]
rt    [, Ensembl_ID := sub("\\..*$", "", Ensembl_ID)]

# ---- Merge
rt_ann <- merge(ann_dt, rt, by = "Ensembl_ID")
cat("合并注释后维度：", dim(rt_ann), "\n")

mat <- as.matrix(rt_ann[, -c("Ensembl_ID","gene"), with = FALSE])
rownames(mat) <- rt_ann$gene
mat <- avereps(mat)
cat("去重后维度：", dim(mat), "\n")

load("C:/Users/赵定康/Desktop/input/pancancer_exp.Rdata")
cat("pancancer_exp 基因数：", nrow(pancancer_exp), "\n")

common_genes <- intersect(rownames(mat), rownames(pancancer_exp))
cat("交集基因数：", length(common_genes), "\n")
mat <- mat[common_genes, ]
cat("裁剪后维度：", dim(mat), "\n")

clin_files <- sub("\\.star_tpm\\.tsv$", ".clinical.tsv", files)
clin_list <- lapply(clin_files, function(f) {
  d <- fread(f, header = TRUE, sep = "\t")
  d <- d[, .(sample, tissue_type.samples)]
  d
})
clin <- rbindlist(clin_list)
setnames(clin, c("Tag","Group"))
clin <- clin[Tag %in% colnames(mat)]
cat("临床样本：", nrow(clin), "\n")


proj_map <- list()
for (i in seq_along(files)) {
  cols_i <- names(fread(files[i], nrows = 0))
  proj_map[[ proj_names[i] ]] <- setdiff(cols_i, "Ensembl_ID")
}
clin[, Project := NA_character_]
for (p in names(proj_map)) {
  clin[Tag %in% proj_map[[p]], Project := p]
}
cat("Project 分布：\n"); print(table(clin$Project, useNA = "ifany"))


clin_tumor <- clin[Group == "Tumor"]
mat_tumor  <- mat[, clin_tumor$Tag]
cat("肿瘤样本数：", ncol(mat_tumor), " 基因数：", nrow(mat_tumor), "\n")
cat("Project 分布（肿瘤）：\n"); print(table(clin_tumor$Project))

mat_log <- log2(mat_tumor + 1)

TARGET_combat <- ComBat(
  dat         = mat_log,
  batch       = clin_tumor$Project,
  par.prior   = TRUE,
  prior.plots = FALSE
)
cat("去批次后维度：", dim(TARGET_combat), "\n")

TARGET_tumor_combat <- as.data.frame(TARGET_combat)
TARGET_G_tumor      <- as.data.frame(clin_tumor)

save(TARGET_tumor_combat, TARGET_G_tumor,
     file = "TARGET_tumor_combat.Rdata")
save(TARGET_tumor_combat, file = "TARGET.Rdata")
save(TARGET_G_tumor,      file = "TARGET_G.Rdata")
cat("已保存\n")
################################################################################TARGET batch-effect check: before vs after correction####
rm(list=ls()); gc()
library(data.table); library(ggplot2); library(patchwork)

dir_input <- "C:/Users/赵定康/Desktop/input"
setwd(dir_input)

files <- c("TARGET-AML.star_tpm.tsv","TARGET-ALL-P1.star_tpm.tsv",
           "TARGET-ALL-P2.star_tpm.tsv","TARGET-ALL-P3.star_tpm.tsv",
           "TARGET-CCSK.star_tpm.tsv","TARGET-NBL.star_tpm.tsv",
           "TARGET-OS.star_tpm.tsv","TARGET-RT.star_tpm.tsv",
           "TARGET-WT.star_tpm.tsv")
proj_names <- c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")

# ---- Read and merge (tumour samples only)
dt_list <- lapply(files, function(f) fread(f, header = TRUE, sep = "\t"))
rt <- Reduce(function(x, y) merge(x, y, by = "Ensembl_ID", all = FALSE), dt_list)

# ---- Annotation
ann <- fread("gencode.v36.annotation.gtf.gene.probemap",
             header = TRUE, sep = "\t")
ann_dt <- data.table(Ensembl_ID = as.character(ann[[1]]),
                     gene       = as.character(ann[[2]]))
ann_dt[, Ensembl_ID := sub("\\..*$", "", Ensembl_ID)]
rt    [, Ensembl_ID := sub("\\..*$", "", Ensembl_ID)]
rt_ann <- merge(ann_dt, rt, by = "Ensembl_ID")

mat <- as.matrix(rt_ann[, -c("Ensembl_ID","gene"), with = FALSE])
rownames(mat) <- rt_ann$gene
mat <- limma::avereps(mat)

# ---- Subset to the pancancer_exp gene set, consistent with the batch-corrected data
load("pancancer_exp.Rdata")
common_genes <- intersect(rownames(mat), rownames(pancancer_exp))
mat <- mat[common_genes, ]
cat("去批次前矩阵：", dim(mat), "\n")

# ---- Read the clinical data, keep tumour samples only, and tag the project
clin_files <- sub("\\.star_tpm\\.tsv$", ".clinical.tsv", files)
clin_list <- lapply(clin_files, function(f) {
  d <- fread(f, header = TRUE, sep = "\t")
  d[, .(sample, tissue_type.samples)]
})
clin <- rbindlist(clin_list)
setnames(clin, c("Tag","Group"))
clin <- clin[Tag %in% colnames(mat) & Group == "Tumor"]

# Tag the project
clin[, Project := NA_character_]
for (i in seq_along(files)) {
  cols_i <- names(fread(files[i], nrows = 0))
  clin[Tag %in% setdiff(cols_i, "Ensembl_ID"), Project := proj_names[i]]
}

mat_tumor <- mat[, clin$Tag]
cat("肿瘤样本：", ncol(mat_tumor), "\n")
print(table(clin$Project))

mat_log <- log2(mat_tumor + 1)

# Take the most variable genes (top 2000)
gvar <- apply(mat_log, 1, var)
gvar <- gvar[!is.na(gvar) & gvar > 0]
top_genes <- names(sort(gvar, decreasing = TRUE))[1:2000]

pca_before <- prcomp(t(mat_log[top_genes, ]), scale. = TRUE)

df_before <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  Project = clin$Project,
  Tag = clin$Tag
)

p_before <- ggplot(df_before, aes(PC1, PC2, color = Project)) +
  geom_point(size = 1, alpha = 0.7) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain"),
    axis.text = element_text(size = 14, color = "black", face = "plain")
  ) +
  labs(title = "Before ComBat",
       subtitle = paste0("PC1 = ", round(summary(pca_before)$importance[2,1]*100, 1),
                         "%, PC2 = ", round(summary(pca_before)$importance[2,2]*100, 1), "%"))

load("TARGET_tumor_combat.Rdata")   # TARGET_tumor_combat, TARGET_G_tumor
mat_after <- as.matrix(TARGET_tumor_combat)
clin_after <- TARGET_G_tumor

# Align the sample order
common <- intersect(colnames(mat_after), clin_after$Tag)
mat_after <- mat_after[, common]
clin_after <- clin_after[match(common, clin_after$Tag), ]

# Use the same top 2000 genes (if they survive batch correction)
top_genes_after <- intersect(top_genes, rownames(mat_after))
pca_after <- prcomp(t(mat_after[top_genes_after, ]), scale. = TRUE)

df_after <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Project = clin_after$Project,
  Tag = clin_after$Tag
)

p_after <- ggplot(df_after, aes(PC1, PC2, color = Project)) +
  geom_point(size = 1, alpha = 0.7) +
  theme_minimal() +
  theme(
    text = element_text(size = 14, color = "black", face = "plain"),
    axis.text = element_text(size = 14, color = "black", face = "plain")
  ) +
  labs(title = "After ComBat",
       subtitle = paste0("PC1 = ", round(summary(pca_after)$importance[2,1]*100, 1),
                         "%, PC2 = ", round(summary(pca_after)$importance[2,2]*100, 1), "%"))

p_before + p_after + plot_layout(guides = "collect")
################################################################################Independent pan-cancer data from TARGET####
rm(list=ls()); gc()
library(GSVA); library(GSEABase)

dir_input <- "C:/Users/赵定康/Desktop/input"

# ---- Load the panWSS definition
load(file.path(dir_input, "panWSS_subtype_definition.Rdata"))
# panWSS_centers, panWSS_features, train_mean_pan, train_sd_pan

geneSets <- getGmt(file.path(dir_input, "wpags_wpigs.gmt"))
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

# ---- Nearest-centre assignment
assign_by_centroid <- function(x, centers) {
  d <- matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    d[, k] <- sqrt(rowSums(
      (x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2
    ))
  }
  rownames(centers)[apply(d, 1, which.min)]
}

# ---- Load the TARGET expression matrix (batch-corrected, tumour samples)
load(file.path(dir_input, "TARGET.Rdata"))   # Object name TARGET_tumor_combat
rt <- TARGET_tumor_combat
cat("TARGET 维度：", dim(rt), "\n")
# ========== 1. Extract the patient ID
get_patient_id <- function(x) {
  x2 <- sub("^TARGET-[0-9]+-", "", x)   # Strip the TARGET-XX- prefix
  parts <- strsplit(x2, "-")
  sapply(parts, function(p) p[1])       # Take the first part as the patient ID
}

sample_ids  <- colnames(rt)
patient_ids <- get_patient_id(sample_ids)

cat("原始样本数：", length(sample_ids), "\n")
cat("唯一病人数：", length(unique(patient_ids)), "\n")

# ========== 2. Average per patient
patient_factor <- factor(patient_ids)

expr_merged <- t(apply(rt, 1, function(x) {
  tapply(x, patient_factor, mean, na.rm = TRUE)
}))
expr_merged <- as.matrix(expr_merged)

rownames(expr_merged) <- rownames(rt)
colnames(expr_merged) <- levels(patient_factor)

expr_merged[is.nan(expr_merged)] <- NA
# ---- 1. Wnt scoring
WNT_Score <- Calculate_bioactivity_scores(
  file_paths         = dir_input,
  expression_profile = as.data.frame(expr_merged),
  foundation         = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt       = geneSets,
  min.sz             = 1,
  max.sz             = 10000,
  geneset_direction  = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                         rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
sm <- as.matrix(WNT_Score$single_score_matrix)
cat("得分矩阵维度：", dim(sm), "\n")

# ---- 2. Align features
missing <- setdiff(panWSS_features, colnames(sm))
cat("缺失特征数：", length(missing), " / ", length(panWSS_features), "\n")
if (length(missing) > 0) {
  print(missing)
  for (m in missing) sm <- cbind(sm, setNames(0, m))
}
sm <- sm[, panWSS_features, drop = FALSE]

# ---- 3. Standardise with the training-set parameters
sm_z <- sweep(sm, 2, train_mean_pan, "-")
sm_z <- sweep(sm_z, 2, train_sd_pan, "/")
cz   <- sweep(panWSS_centers, 2, train_mean_pan, "-")
cz   <- sweep(cz, 2, train_sd_pan, "/")

# ---- 4. Assign the panWSS subtype
pred <- assign_by_centroid(sm_z, cz)
print(table(pred))

# ---- 5. Save
TARGET_panWSS <- data.frame(Tag = rownames(sm_z), subtype = pred,
                            dataset = "TARGET", stringsAsFactors = FALSE)
save(TARGET_panWSS,
     file = file.path(dir_input, "TARGET_panWSS_subtype.Rdata"))
cat("已保存：TARGET_panWSS_subtype.Rdata\n")
cat("样本数：", nrow(TARGET_panWSS), "\n")

# ---- 6. Also save sm_z; it is needed later for PCA and biological scoring
save(sm_z, file = file.path(dir_input, "TARGET_sm_z.Rdata"))
cat("已保存：TARGET_sm_z.Rdata（63维标准化得分）\n")
################################################################################Visualising the panWSS subtype proportions in TARGET####
rm(list=ls()); gc()
load("C:/Users/赵定康/Desktop/input/TARGET_panWSS_subtype.Rdata")
library(dplyr); library(ggplot2)

# ---- Proportions
prop_df <- TARGET_panWSS %>%
  group_by(subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n),
         label = sprintf("%s\n(n=%d, %.1f%%)", subtype, n, prop * 100))

print(prop_df)

# ---- Pie chart
library(ggplot2)

# ---- Compute the label positions (offset plus leader lines)
prop_df <- prop_df[order(prop_df$subtype), ]
prop_df$pos <- cumsum(prop_df$prop) - prop_df$prop / 2   # Sector centre angle

library(dplyr); library(ggplot2)

prop_df <- TARGET_panWSS %>%
  group_by(subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(prop = n / sum(n),
         label = sprintf("%s\n(n=%d, %.1f%%)", subtype, n, prop * 100))

prop_df$subtype <- factor(prop_df$subtype,
                          levels = c("panWSS1","panWSS2","panWSS3","panWSS4"))

ggplot(prop_df, aes(x = "", y = prop, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("#d71345","#f47920","#145b7d","#6950a1")) +
  geom_text(aes(label = sprintf("%.1f%%", prop * 100)),
            position = position_stack(vjust = 0.5),
            size = 5, color = "white") +
  theme_minimal() +
  theme(
    text        = element_text(size = 14, color = "black", face = "plain"),
    axis.text   = element_text(size = 14, color = "black", face = "plain"),
    axis.title  = element_text(size = 14, color = "black", face = "plain"),
    plot.title  = element_text(size = 14, color = "black", face = "plain"),
    legend.text = element_text(size = 14, color = "black", face = "plain"),
    legend.title= element_text(size = 14, color = "black", face = "plain"),
    panel.grid  = element_blank()
  ) +
  labs(x = "", y = "Proportion", fill = "",
       title = paste0("panWSS subtype composition in TARGET (n=",
                      nrow(TARGET_panWSS), ")"))
################################################################################panWSS x Hallmark: TCGA vs TARGET biological pattern comparison####
rm(list=ls()); gc()
library(GSVA); library(GSEABase); library(dplyr); library(tidyr)
library(tibble); library(ggplot2)

dir_input <- "C:/Users/赵定康/Desktop/input"

# ---- 1. Hallmark gene sets, dropping WNT_beta_catenin_signaling
geneSets_H <- getGmt(file.path(dir_input, "geneset_hallmark.gmt"))
rm_name <- grep("WNT.*catenin|catenin.*WNT|WNT_beta",
                names(geneSets_H), ignore.case = TRUE, value = TRUE)
cat("移除：", paste(rm_name, collapse=", "), "\n")
geneSets_H <- geneSets_H[setdiff(names(geneSets_H), rm_name)]
cat("剩余 Hallmark 数：", length(geneSets_H), "\n")

# ---- 2. Load the section A definitions (includes pancancer_labels, excludes rank_df)
load(file.path(dir_input, "panWSS_subtype_definition.Rdata"))
cat("panWSS_centers 行名：",
    paste(rownames(panWSS_centers), collapse=", "), "\n")
cat("pancancer_labels 亚型分布：\n"); print(table(pancancer_labels$subtype))

# ---- 3. Scoring function
compute_profile <- function(expr, labels, geneSets_H) {
  gs <- gsva(as.matrix(expr), geneSets_H, kcdf = "Gaussian",
             method = "ssgsea", min.sz = 1, max.sz = 10000)
  gs <- as.data.frame(scale(t(gs)))
  gs$Tag <- rownames(gs)
  df <- merge(gs, labels, by = "Tag")
  cols <- setdiff(colnames(df), c("Tag","subtype"))
  df %>% group_by(subtype) %>%
    summarise(across(all_of(cols), mean), .groups = "drop")
}

# ---- 4. TCGA baseline
load(file.path(dir_input, "pancancer_group.Rdata"))
load(file.path(dir_input, "pancancer_exp.Rdata"))

tcga_group <- pancancer_group[pancancer_group$Group == "Tumor", ]
expr_tcga  <- pancancer_exp[, rownames(tcga_group)]

labels_tcga <- pancancer_labels
labels_tcga <- labels_tcga[labels_tcga$Tag %in% colnames(expr_tcga), ]
cat("TCGA 匹配样本：", nrow(labels_tcga), "\n")
print(table(labels_tcga$subtype))

profile_tcga <- compute_profile(expr_tcga, labels_tcga, geneSets_H)
profile_tcga$dataset <- "TCGA (baseline)"
cat("profile_tcga 亚型：", paste(profile_tcga$subtype, collapse=", "), "\n")

# ---- 5. TARGET
load(file.path(dir_input, "TARGET.Rdata"))                # TARGET_tumor_combat
load(file.path(dir_input, "TARGET_panWSS_subtype.Rdata")) # TARGET_panWSS

expr_target <- TARGET_tumor_combat

# Key step: replace "." with "-" in TARGET_panWSS$Tag before building labels_target
TARGET_panWSS$Tag <- gsub("\\.", "-", TARGET_panWSS$Tag)

labels_target <- data.frame(Tag     = TARGET_panWSS$Tag,
                            subtype = TARGET_panWSS$subtype,
                            stringsAsFactors = FALSE)

cat("TARGET Tag 交集数：",
    sum(labels_target$Tag %in% colnames(expr_target)), "\n")
labels_target <- labels_target[labels_target$Tag %in% colnames(expr_target), ]
cat("TARGET 匹配样本：", nrow(labels_target), "\n")
print(table(labels_target$subtype))

profile_target <- compute_profile(expr_target, labels_target, geneSets_H)
profile_target$dataset <- "TARGET"
cat("profile_target 亚型：", paste(profile_target$subtype, collapse=", "), "\n")

# ---- 6. Spearman ρ
make_mat <- function(prof) {
  df <- as.data.frame(prof)
  rownames(df) <- df$subtype
  df$subtype <- NULL
  df$dataset <- NULL
  as.matrix(df)
}

b <- make_mat(profile_tcga)
x <- make_mat(profile_target)
cat("TCGA  行名：", paste(rownames(b), collapse=", "), "\n")
cat("TARGET 行名：", paste(rownames(x), collapse=", "), "\n")

cs <- intersect(rownames(b), rownames(x))
cc <- intersect(colnames(b), colnames(x))
cat("共同亚型：", paste(cs, collapse=", "),
    "  共同特征数：", length(cc), "\n")

if (length(cs) >= 2 && length(cc) >= 2) {
  rho <- cor(as.numeric(b[cs, cc]), as.numeric(x[cs, cc]),
             method = "spearman")
  cat("TARGET vs TCGA 的 Spearman ρ =", round(rho, 3), "\n")
} else {
  stop("亚型或特征交集不足，检查 subtype 命名是否统一")
}

# ---- 7. Heatmap
prof_all  <- rbind(profile_tcga, profile_target)
prof_long <- prof_all %>%
  pivot_longer(-c(subtype, dataset), names_to = "feature", values_to = "z")

ord <- prof_long %>% filter(dataset == "TCGA (baseline)") %>%
  group_by(feature) %>% summarise(m = mean(z)) %>%
  arrange(m) %>% pull(feature)
prof_long$feature <- factor(prof_long$feature, levels = ord)

prof_long$subtype <- factor(prof_long$subtype,
                            levels = c("panWSS1","panWSS2","panWSS3","panWSS4"))
prof_long$dataset <- factor(prof_long$dataset,
                            levels = c("TCGA (baseline)","TARGET"))

label_rho <- paste0("Spearman ρ = ", round(rho, 3))

ggplot(prof_long, aes(feature, subtype, fill = z)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ dataset) +
  scale_fill_gradient2(low = "#90d7ec", mid = "white", high = "#f58f98",
                       midpoint = 0, name = "z-score") +
  annotate("text",
           x = length(ord) * 0.98, y = 4.6,
           label = label_rho,
           hjust = 1, vjust = 1, size = 5, color = "black") +
  theme_minimal() +
  theme(
    text        = element_text(size = 14, color = "black", face = "plain"),
    axis.text   = element_text(size = 12, color = "black", face = "plain"),
    axis.text.x = element_text(size = 12, color = "black", face = "plain",
                               angle = 45, hjust = 1),
    strip.text  = element_text(size = 14, color = "black", face = "plain"),
    panel.grid  = element_blank()
  ) +
  labs(x = "", y = "",
       title = "panWSS × Hallmark features: TCGA (baseline) vs TARGET")
################################################################################TARGET prognostic analysis####
rm(list=ls()); gc()
library(GSVA); library(GSEABase)

dir_input <- "C:/Users/赵定康/Desktop/input"

# ---- Load the panWSS definition
load(file.path(dir_input, "TARGET_panWSS_subtype.Rdata"))
load(file.path(dir_input, "TARGET_sur.Rdata"))

library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)

# ---- 0. Aliases and a uniform ID format
TARGET_cluster <- TARGET_sur
TARGET_cluster$sample <- gsub("\\.", "-", TARGET_cluster$sample)
TARGET_panWSS$Tag    <- gsub("\\.", "-", TARGET_panWSS$Tag)

# ---- 1. Merge
survival_data_cms <- merge(TARGET_panWSS,
                           TARGET_cluster,
                           by.x = "Tag", by.y = "sample")

cat("合并后行数:", nrow(survival_data_cms), "\n")
head(survival_data_cms)

# ---- 2. Set the levels
survival_data_cms$level <- survival_data_cms$subtype
survival_data_cms$level <- factor(survival_data_cms$level,
                                  levels = c("panWSS1","panWSS2","panWSS3","panWSS4"))

# ---- 3. Drop missing values
survival_data_cms <- survival_data_cms[!is.na(survival_data_cms$OS) &
                                         !is.na(survival_data_cms$OS.time) &
                                         !is.na(survival_data_cms$level), ]

# Extract the patient ID: strip the -03A / -04A / -09A / -09B / -40A suffixes
survival_data_cms$patient <- sub("-[0-9]+[A-Z]?$", "", survival_data_cms$Tag)

cat("去重前样本数:", nrow(survival_data_cms), "\n")
cat("独立病人数:", length(unique(survival_data_cms$patient)), "\n")

# Prefer primary samples (-03A or -09A), one per patient
survival_data_cms <- survival_data_cms[order(survival_data_cms$Tag), ]
survival_data_cms$is_dx <- grepl("-(03A|09A)$", survival_data_cms$Tag)

# Take the primary samples first
dx <- survival_data_cms[survival_data_cms$is_dx, ]
dx <- dx[!duplicated(dx$patient), ]

# Patients without a primary sample keep their first sample
non_dx <- survival_data_cms[!survival_data_cms$is_dx, ]
non_dx <- non_dx[!non_dx$patient %in% dx$patient, ]
non_dx <- non_dx[!duplicated(non_dx$patient), ]

survival_data_cms <- rbind(dx, non_dx)
survival_data_cms <- survival_data_cms[!is.na(survival_data_cms$level), ]

cat("去重后样本数:", nrow(survival_data_cms), "\n")
print(table(survival_data_cms$level))

# ---- 4. log-rank
fitd <- survdiff(Surv(OS.time, OS) ~ level,
                 data = survival_data_cms, na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# ---- 5. KM
fit <- survfit(Surv(OS.time, OS) ~ level,
               data = survival_data_cms,
               type = "kaplan-meier",
               error = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# ---- 6. Pairwise comparisons
ps <- pairwise_survdiff(Surv(OS.time, OS) ~ level,
                        data = survival_data_cms,
                        p.adjust.method = "BH")
print(ps)

# ---- 7. Plot
mycol <- c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  conf.int = FALSE,
  risk.table = FALSE,
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
  title = paste0("OS(TARGET-Pancancer,n=", nrow(survival_data_cms), ")"),
  ylim = c(0, 1)
)

p.lab <- paste0("log-rank test p",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))

p$plot <- p$plot +
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,
           fontface = "italic", label = p.lab, size = 6) +
  theme(text = element_text(size = 16))

print(p)
################################################################################Independent ICGC pan-cancer data####
rm(list=ls()); gc()
library(GSVA); library(GSEABase)

dir_input <- "C:/Users/赵定康/Desktop/input"

# ---- Load the panWSS definition
load(file.path(dir_input, "panWSS_subtype_definition.Rdata"))
# panWSS_centers, panWSS_features, train_mean_pan, train_sd_pan

geneSets <- getGmt(file.path(dir_input, "wpags_wpigs.gmt"))
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

# ---- Nearest-centre assignment
assign_by_centroid <- function(x, centers) {
  d <- matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    d[, k] <- sqrt(rowSums(
      (x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2
    ))
  }
  rownames(centers)[apply(d, 1, which.min)]
}

library(data.table)
rt <- as.data.frame(fread("C:/Users/赵定康/Desktop/input/donor_exp_seq.all_projects.donor.USonly.xena.tsv.gz"))
rownames(rt)=rt[,1]
rt=rt[,-1]
# ---- 1. Wnt scoring
WNT_Score <- Calculate_bioactivity_scores(
  file_paths         = dir_input,
  expression_profile = as.data.frame(rt),
  foundation         = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt       = geneSets,
  min.sz             = 1,
  max.sz             = 10000,
  geneset_direction  = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                         rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
sm <- as.matrix(WNT_Score$single_score_matrix)
cat("得分矩阵维度：", dim(sm), "\n")

# ---- 2. Align features
missing <- setdiff(panWSS_features, colnames(sm))
cat("缺失特征数：", length(missing), " / ", length(panWSS_features), "\n")
if (length(missing) > 0) {
  print(missing)
  for (m in missing) sm <- cbind(sm, setNames(0, m))
}
sm <- sm[, panWSS_features, drop = FALSE]

# ---- 3. Standardise with the training-set parameters
sm_z <- sweep(sm, 2, train_mean_pan, "-")
sm_z <- sweep(sm_z, 2, train_sd_pan, "/")
cz   <- sweep(panWSS_centers, 2, train_mean_pan, "-")
cz   <- sweep(cz, 2, train_sd_pan, "/")

# ---- 4. Assign the panWSS subtype
pred <- assign_by_centroid(sm_z, cz)
print(table(pred))

# ---- 5. Save
ICGC_panWSS <- data.frame(Tag = rownames(sm_z), subtype = pred,
                            dataset = "ICGC", stringsAsFactors = FALSE)
save(ICGC_panWSS,
     file = file.path(dir_input, "ICGC_panWSS_subtype.Rdata"))
cat("已保存：ICGC_panWSS_subtype.Rdata\n")
cat("样本数：", nrow(ICGC_panWSS), "\n")

# ---- 6. Also save sm_z; it is needed later for PCA and biological scoring
save(sm_z, file = file.path(dir_input, "ICGC_sm_z.Rdata"))
cat("已保存：ICGC_sm_z.Rdata（63维标准化得分）\n")
################################################################################ICGC prognostic analysis####
rm(list=ls()); gc()
library(GSVA); library(GSEABase)

dir_input <- "C:/Users/赵定康/Desktop/input"

# ---- Load the panWSS definition
load(file.path(dir_input, "ICGC_panWSS_subtype.Rdata"))
ICGC_sur=as.data.frame(fread("C:/Users/赵定康/Desktop/input/sp_donor.all_projects.overallSurvival_transfer_specimen.gz"))
ICGC_nor=as.data.frame(fread("C:/Users/赵定康/Desktop/input/sp_donor.all_projects.phenotype_transfer_specimen.gz"))
ICGC_nor <- ICGC_nor[!grepl("TCGA", ICGC_nor$submitted_donor_id), ]
ICGC_sur=data.frame(Tag=ICGC_sur$xena_sample,OS.time=ICGC_sur$OS.time,OS=ICGC_sur$OS)
ICGC_sur=ICGC_sur[ICGC_sur$Tag%in%ICGC_nor$xena_sample,]
ICGC_sur=na.omit(ICGC_sur)
ICGC_sur$OS.time=ICGC_sur$OS.time/30
library(survival)
library(survminer)
library(RColorBrewer)
library(tibble)
library(ggpp)
# ---- 1. Merge
survival_data_cms <- merge(ICGC_panWSS,
                           ICGC_sur,
                           by = "Tag",all=T)

# ---- 2. Set the levels
survival_data_cms$level <- survival_data_cms$subtype
survival_data_cms$level <- factor(survival_data_cms$level,
                                  levels = c("panWSS1","panWSS2","panWSS3","panWSS4"))


# ---- 4. log-rank
fitd <- survdiff(Surv(OS.time, OS) ~ level,
                 data = survival_data_cms, na.action = na.exclude)
p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

# ---- 5. KM
fit <- survfit(Surv(OS.time, OS) ~ level,
               data = survival_data_cms,
               type = "kaplan-meier",
               error = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

# ---- 6. Pairwise comparisons
ps <- pairwise_survdiff(Surv(OS.time, OS) ~ level,
                        data = survival_data_cms,
                        p.adjust.method = "BH")
print(ps)

# ---- 7. Plot
mycol <- c("#d71345","#f47920","#145b7d","#6950a1")
names(fit$strata) <- gsub("level=", "", names(fit$strata))

p <- ggsurvplot(
  fit = fit,
  conf.int = FALSE,
  risk.table = FALSE,
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
  title = paste0("OS(TARGET-Pancancer,n=", nrow(survival_data_cms), ")"),
  ylim = c(0, 1)
)

p.lab <- paste0("log-rank test p",
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ", round(p.val, 3))))

p$plot <- p$plot +
  annotate("text", x = 0, y = 0, hjust = 0, vjust = 0,
           fontface = "italic", label = p.lab, size = 6) +
  theme(text = element_text(size = 16))

print(p)

################################################################################panWSS: save a transferable subtype definition (route A: PCA + RBF soft assignment, MAP for the final subtype)####
################################################################################
## panWSS: save a transferable subtype definition (route A: PCA + RBF soft assignment, MAP for the final subtype)
################################################################################
rm(list=ls()); gc()
library(GSEABase); library(factoextra)

setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/pancancer_WNT_Score.Rdata")
load("C:/Users/赵定康/Desktop/input/pancancer_group.Rdata")

genesets <- as.matrix(pancancer_WNT_Score$single_score_matrix)
pancancer_group <- pancancer_group[pancancer_group$Group == "Tumor", ]
genesets <- genesets[rownames(pancancer_group), ]

# ---- 1. k-means (fixed seed)
set.seed(123)
km_result <- kmeans(genesets, centers = 4, nstart = 1000,
                    algorithm = "Hartigan-Wong", iter.max = 100)

# ---- 2. Median activity per raw cluster, used to reorder the clusters
raw_cluster <- km_result$cluster
acts <- sapply(1:4, function(k) {
  median(pancancer_WNT_Score$final_activity_score[
    rownames(genesets)[raw_cluster == k], "activity_score"])
})
rank_df <- data.frame(raw = 1:4, value = acts)
rank_df <- rank_df[order(-rank_df$value), ]
cat("原始 cluster 按 activity 排序：\n"); print(rank_df)

# ---- 3. Reorder centres by activity from high to low; rename the rows panWSS1-4
panWSS_centers <- km_result$centers
panWSS_centers <- panWSS_centers[rank_df$raw, , drop = FALSE]
rownames(panWSS_centers) <- paste0("panWSS", 1:4)
cat("panWSS_centers 行名：",
    paste(rownames(panWSS_centers), collapse=", "), "\n")

# ---- 4. Training-set parameters (kept for standardisation and alignment)
panWSS_features <- colnames(genesets)
train_mean_pan  <- colMeans(genesets)
train_sd_pan    <- apply(genesets, 2, sd)
train_sd_pan[train_sd_pan == 0] <- 1e-8

# ---- 5. Training-set labels (panWSS1-4 ordered by activity)
orig_panWSS <- paste0("panWSS", match(raw_cluster, rank_df$raw))

# ---- 6. Fit the PCA on the training set and choose the number of components
pca_fit <- prcomp(genesets, center = TRUE, scale. = FALSE)
var_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
n_pc <- which(cumsum(var_explained) >= 0.9)[1]
n_pc <- max(n_pc, 10)                       # At least 10 dimensions
n_pc <- min(n_pc, ncol(pca_fit$rotation))   # Do not exceed the available dimensions
cat("保留 PC 数：", n_pc, "，累计方差：",
    round(sum(var_explained[1:n_pc]) * 100, 2), "%\n")

# Slim PCA object kept for transfer (small file, few dependencies)
pca_fit_slim <- list(
  center   = pca_fit$center,
  rotation = pca_fit$rotation[, 1:n_pc, drop = FALSE],
  sdev     = pca_fit$sdev[1:n_pc]
)

# Project the centres into PC space
centers_pc <- predict(pca_fit, newdata = panWSS_centers)[, 1:n_pc, drop = FALSE]

# Set the softmax temperature from the median within-class distance of the training set
proj_train <- pca_fit$x[, 1:n_pc, drop = FALSE]
within_dist <- sapply(seq_len(nrow(panWSS_centers)), function(k) {
  idx <- orig_panWSS == rownames(panWSS_centers)[k]
  if (sum(idx) < 2) return(NA_real_)
  d <- sqrt(rowSums((proj_train[idx, , drop = FALSE] -
                       matrix(centers_pc[k, ], sum(idx), n_pc, byrow = TRUE))^2))
  median(d)
})
temp_pan <- median(within_dist, na.rm = TRUE)
cat("softmax 温度 temp_pan =", round(temp_pan, 4), "\n")

# ---- 7. Transfer function: PCA projection + RBF soft assignment + MAP subtype
project_pca <- function(x, pca_fit_slim) {
  xc <- sweep(x, 2, pca_fit_slim$center, "-")
  xc %*% pca_fit_slim$rotation
}

assign_soft_pca <- function(x, pca_fit_slim, centers_pc,
                            temp, features, ambiguous_cut = 0.6) {
  x <- x[, features, drop = FALSE]
  x_pc <- project_pca(x, pca_fit_slim)

  # Euclidean distance to each centre
  d <- sapply(seq_len(nrow(centers_pc)), function(k) {
    sqrt(rowSums((x_pc - matrix(centers_pc[k, ], nrow(x_pc),
                                ncol(x_pc), byrow = TRUE))^2))
  })
  colnames(d) <- rownames(centers_pc)

  # RBF soft assignment
  w <- exp(-d^2 / (2 * temp^2))
  prob <- w / rowSums(w)
  colnames(prob) <- rownames(centers_pc)
  max_prob <- apply(prob, 1, max)

  # Vector columns only
  out <- data.frame(
    Tag       = rownames(x),
    subtype   = colnames(prob)[max.col(prob, ties.method = "first")],
    max_prob  = max_prob,
    ambiguous = max_prob < ambiguous_cut,
    stringsAsFactors = FALSE
  )
  # cbind the probability matrix into four separate columns
  cbind(out, as.data.frame(prob))
}

# ---- 8. Training-set consistency check
self_check <- assign_soft_pca(genesets, pca_fit_slim, centers_pc,
                              temp_pan, panWSS_features)
cat("训练集硬标签自洽率：",
    round(mean(self_check$subtype == orig_panWSS) * 100, 2), "%\n")
cat("边界样本比例：",
    round(mean(self_check$ambiguous) * 100, 2), "%\n")
print(table(pred = self_check$subtype, orig = orig_panWSS))

# If the consistency rate is clearly below 100%, scan n_pc for a suitable dimension
if (mean(self_check$subtype == orig_panWSS) < 0.99) {
  cat("\n自洽率偏低，扫描 n_pc：\n")
  for (np in c(10, 20, 30, 50)) {
    np <- min(np, ncol(pca_fit$rotation))
    cp <- predict(pca_fit, newdata = panWSS_centers)[, 1:np, drop = FALSE]
    pt <- pca_fit$x[, 1:np, drop = FALSE]
    wd <- sapply(seq_len(nrow(panWSS_centers)), function(k) {
      idx <- orig_panWSS == rownames(panWSS_centers)[k]
      median(sqrt(rowSums((pt[idx, , drop = FALSE] -
                             matrix(cp[k, ], sum(idx), np, byrow = TRUE))^2)))
    })
    tp <- median(wd, na.rm = TRUE)
    pca_slim_np <- list(center = pca_fit$center,
                        rotation = pca_fit$rotation[, 1:np, drop = FALSE],
                        sdev = pca_fit$sdev[1:np])
    s <- assign_soft_pca(genesets, pca_slim_np, cp, tp, panWSS_features)
    cat("n_pc =", np, " 自洽率 =",
        round(mean(s$subtype == orig_panWSS) * 100, 2),
        "%  边界比例 =", round(mean(s$ambiguous) * 100, 2), "%\n")
  }
}

# ---- 9. panWSS1-4 label per TCGA sample (with probabilities)
pancancer_labels <- self_check   # Already 8 columns: Tag, subtype, max_prob,
# ambiguous, panWSS1..4
cat("TCGA 样本 panWSS 标签分布：\n")
print(table(pancancer_labels$subtype))
cat("边界样本比例：",
    round(mean(pancancer_labels$ambiguous) * 100, 2), "%\n")
cat("pancancer_labels 结构：\n")
str(pancancer_labels)

# ---- 10. Save
save(panWSS_centers, panWSS_features,
     train_mean_pan, train_sd_pan,
     pancancer_labels,
     pca_fit    = pca_fit_slim,
     centers_pc = centers_pc,
     n_pc       = n_pc,
     temp_pan   = temp_pan,
     file = "panWSS_subtype_definition.Rdata")
cat("已保存：panWSS_subtype_definition.Rdata\n")
cat("centers 维度：", dim(panWSS_centers), "\n")
cat("centers 行名：",
    paste(rownames(panWSS_centers), collapse=", "), "\n")
cat("pancancer_labels 行数：", nrow(pancancer_labels), "\n")
cat("保存对象：panWSS_centers, panWSS_features, train_mean_pan, train_sd_pan,\n")
cat("          pancancer_labels, pca_fit(slim), centers_pc, n_pc, temp_pan\n")
################################################################################
## TARGET independent pan-cancer data -- transfer the training panWSS definition (PCA + RBF soft assignment)
################################################################################
rm(list=ls()); gc()
library(GSVA); library(GSEABase)

dir_input <- "C:/Users/赵定康/Desktop/input"

# ---- Load the panWSS definition (route A)
load(file.path(dir_input, "panWSS_subtype_definition.Rdata"))
# Contains: panWSS_centers, panWSS_features, train_mean_pan, train_sd_pan,
#       pancancer_labels, pca_fit(slim), centers_pc, n_pc, temp_pan

geneSets <- getGmt(file.path(dir_input, "wpags_wpigs.gmt"))
source("C:/Users/赵定康/Desktop/要整合成包的函数/Calculate_bioactivity_scores.R")

# ---- Transfer function: PCA projection + RBF soft assignment + MAP subtype
project_pca <- function(x, pca_fit) {
  xc <- sweep(x, 2, pca_fit$center, "-")
  xc %*% pca_fit$rotation
}

assign_soft_pca <- function(x, pca_fit, centers_pc,
                            temp, features, ambiguous_cut = 0.6) {
  x <- x[, features, drop = FALSE]
  x_pc <- project_pca(x, pca_fit)
  d <- sapply(seq_len(nrow(centers_pc)), function(k) {
    sqrt(rowSums((x_pc - matrix(centers_pc[k, ], nrow(x_pc),
                                ncol(x_pc), byrow = TRUE))^2))
  })
  colnames(d) <- rownames(centers_pc)
  w <- exp(-d^2 / (2 * temp^2))
  prob <- w / rowSums(w)
  colnames(prob) <- rownames(centers_pc)
  max_prob <- apply(prob, 1, max)
  out <- data.frame(
    Tag       = rownames(x),
    subtype   = colnames(prob)[max.col(prob, ties.method = "first")],
    max_prob  = max_prob,
    ambiguous = max_prob < ambiguous_cut,
    stringsAsFactors = FALSE
  )
  cbind(out, as.data.frame(prob))
}

# ---- Load the TARGET expression matrix (batch-corrected, tumour samples)
load(file.path(dir_input, "TARGET.Rdata"))   # Object name TARGET_tumor_combat
rt <- TARGET_tumor_combat
cat("TARGET 维度：", dim(rt), "\n")

# ========== 1. Extract the patient ID
get_patient_id <- function(x) {
  x2 <- sub("^TARGET-[0-9]+-", "", x)   # Strip the TARGET-XX- prefix
  parts <- strsplit(x2, "-")
  sapply(parts, function(p) p[1])       # Take the first part as the patient ID
}

sample_ids  <- colnames(rt)
patient_ids <- get_patient_id(sample_ids)

cat("原始样本数：", length(sample_ids), "\n")
cat("唯一病人数：", length(unique(patient_ids)), "\n")

# ========== 2. Average per patient
patient_factor <- factor(patient_ids)

expr_merged <- t(apply(rt, 1, function(x) {
  tapply(x, patient_factor, mean, na.rm = TRUE)
}))
expr_merged <- as.matrix(expr_merged)

rownames(expr_merged) <- rownames(rt)
colnames(expr_merged) <- levels(patient_factor)

expr_merged[is.nan(expr_merged)] <- NA

# ---- 1. Wnt scoring
WNT_Score <- Calculate_bioactivity_scores(
  file_paths         = dir_input,
  expression_profile = as.data.frame(expr_merged),
  foundation         = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt       = geneSets,
  min.sz             = 1,
  max.sz             = 10000,
  geneset_direction  = c(rep("pos", sum(grepl("WPAGS", names(geneSets)))),
                         rep("neg", sum(grepl("WPIGS", names(geneSets)))))
)
sm <- as.matrix(WNT_Score$single_score_matrix)
cat("得分矩阵维度：", dim(sm), "\n")

# ---- 2. Align features
missing <- setdiff(panWSS_features, colnames(sm))
cat("缺失特征数：", length(missing), " / ", length(panWSS_features), "\n")
if (length(missing) > 0) {
  print(missing)
  for (m in missing) sm <- cbind(sm, setNames(0, m))
}
sm <- sm[, panWSS_features, drop = FALSE]

# ---- 3. Standardise with the training-set parameters
sm_z <- sweep(sm, 2, train_mean_pan, "-")
sm_z <- sweep(sm_z, 2, train_sd_pan, "/")

# Note: the PCA projection uses the training pca_fit (slim), which already centres the data,
# The standardised sm_z is therefore centred once more by pca_fit$center here.
# If prcomp used the raw (unstandardised) genesets during training, pass the raw sm here.
# See "Notes on standardisation" below.

# ---- 4. Assign the panWSS subtype (PCA + RBF soft assignment)
pred_df <- assign_soft_pca(
  x         = sm_z,
  pca_fit   = pca_fit,
  centers_pc = centers_pc,
  temp      = temp_pan,
  features  = panWSS_features
)
print(table(pred_df$subtype))
cat("边界样本比例：",
    round(mean(pred_df$ambiguous) * 100, 2), "%\n")

# ---- 5. Save
TARGET_panWSS <- data.frame(
  Tag       = pred_df$Tag,
  subtype   = pred_df$subtype,
  max_prob  = pred_df$max_prob,
  ambiguous = pred_df$ambiguous,
  pred_df[, c("panWSS1", "panWSS2", "panWSS3", "panWSS4")],
  dataset   = "TARGET",
  stringsAsFactors = FALSE
)
save(TARGET_panWSS,
     file = file.path(dir_input, "TARGET_panWSS_subtype.Rdata"))
cat("已保存：TARGET_panWSS_subtype.Rdata\n")
cat("样本数：", nrow(TARGET_panWSS), "\n")

# ---- 6. Also save sm_z; it is needed later for PCA and biological scoring
save(sm_z, file = file.path(dir_input, "TARGET_sm_z.Rdata"))
cat("已保存：TARGET_sm_z.Rdata（标准化得分矩阵）\n")
