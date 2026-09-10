# ==============================================================================
# WntQuant-Paper7 -- Single-cell and spatial analyses
# ------------------------------------------------------------------------------
# Purpose:    Resolve Wnt pathway activity at single-cell and spatial resolution using Seurat, harmony, and SCP-based workflows.
# Inputs:     Single-cell RNA-seq object and cell metadata.
# Outputs:    Cell-level Wnt activity annotations and visualizations.
# ==============================================================================

library(Seurat)
library(harmony)
library(mixOmics)
library(GSEABase)
library(ggplot2)
library(SCP)
library(patchwork)
library(scales)

set.seed(888)

setwd("data/")

scRNA <- readRDS("seurat_obj.rds")
metadata <- read.csv("Epithelial_metadata.csv")
metadata$iCMS_msi <- paste(metadata$iCMS, metadata$msi, sep = "_")

rownames(metadata) <- metadata$cell.ID
metadata <- metadata[colnames(scRNA), ]

scRNA <- AddMetaData(
  object = scRNA,
  metadata = metadata
)

scRNA$iCMS_msi <- as.character(scRNA$iCMS_msi)

scRNA$iCMS_msi[
  grepl("Normal", scRNA$iCMS_msi, ignore.case = TRUE)
] <- "Normal"

group_levels <- c(
  "Normal",
  "iCMS2_MSS",
  "iCMS3_MSI-H",
  "iCMS3_MSS"
)

scRNA$iCMS_msi <- factor(
  scRNA$iCMS_msi,
  levels = group_levels
)

Idents(scRNA) <- "iCMS_msi"
scRNA@assays$RNA <- as(
  object = scRNA@assays$RNA,
  Class = "Assay"
)

scRNA <- NormalizeData(
  scRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

scRNA <- FindVariableFeatures(
  scRNA,
  selection.method = "vst",
  nfeatures = 5000
)

scRNA <- ScaleData(
  scRNA,
  features = VariableFeatures(scRNA)
)
scRNA <- RunPCA(
  scRNA,
  features = VariableFeatures(scRNA),
  npcs = 50
)

scRNA <- RunHarmony(
  scRNA,
  group.by.vars = "dataset",
  reduction.use = "pca",
  dims.use = 1:20,
  nclust = 50,
  max_iter = 10,
  early_stop = TRUE,
  plot_convergence = FALSE
)
plsda_matrix <- t(
  as.matrix(scRNA@assays$RNA$scale.data)
)

group_factor <- factor(
  scRNA$iCMS_msi,
  levels = group_levels
)

keep_cells <- !is.na(group_factor)

plsda_matrix <- plsda_matrix[keep_cells, , drop = FALSE]
group_factor <- droplevels(group_factor[keep_cells])

plsda_result <- plsda(
  X = plsda_matrix,
  Y = group_factor,
  ncomp = 4
)

plsda_embeddings <- plsda_result$variates$X

scRNA_plsda <- subset(
  scRNA,
  cells = rownames(plsda_embeddings)
)

scRNA_plsda[["plsda"]] <- CreateDimReducObject(
  embeddings = plsda_embeddings,
  key = "PLSDA_",
  assay = DefaultAssay(scRNA_plsda)
)
scRNA_plsda <- RunTSNE(
  scRNA_plsda,
  reduction = "plsda",
  dims = 1:4,
  perplexity = 30,
  seed.use = 888
)

scRNA <- scRNA_plsda
work_path <- "/home/user/Fanglab1/ymh/dk/"

source(
  file.path(
    work_path,
    "Calculate_bioactivity_scores.R"
  )
)

geneSets <- getGmt(
  file.path(
    work_path,
    "wpags_wpigs.gmt"
  )
)

exp_matrix <- GetAssayData(
  scRNA,
  assay = "RNA",
  layer = "data"
)

WNT_Score <- Calculate_bioactivity_scores(
  file_paths = work_path,
  expression_profile = exp_matrix,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(
    rep(
      "pos",
      sum(grepl("WPAGS", names(geneSets)))
    ),
    rep(
      "neg",
      sum(grepl("WPIGS", names(geneSets)))
    )
  )
)
scRNA$WNT_Activity <- WNT_Score$final_activity_score$activity_score
scRNA$WNT_Pos <- WNT_Score$final_activity_score$pos_score
scRNA$WNT_Neg <- WNT_Score$final_activity_score$neg_score
target_genes <- c(
  "AXIN2",
  "RNF43",
  "ZNRF3",
  "LGR5",
  "TCF7"
)

for (gene in target_genes) {
  
  p <- FeaturePlot(
    scRNA,
    features = gene,
    reduction = "tsne",
    pt.size = 0.1,
    order = TRUE,
    max.cutoff = "q99"
  ) +
    scale_color_gradientn(
      colors = c(
        "grey90",
        "#afb4db",
        "#9b95c9",
        "#6950a1"
      )
    ) +
    theme_classic(base_size = 14) +
    ggtitle(gene)
  
  ggsave(
    filename = paste0(
      "FeaturePlot_tSNE_",
      gene,
      ".pdf"
    ),
    plot = p,
    width = 6,
    height = 5
  )
}
wnt_scores <- scRNA$WNT_Activity
wnt_scores <- wnt_scores[!is.na(wnt_scores)]

min_val <- min(wnt_scores)
q90 <- quantile(wnt_scores, 0.90)
q93 <- quantile(wnt_scores, 0.93)
q95 <- quantile(wnt_scores, 0.95)
q99 <- quantile(wnt_scores, 0.99)

q_values <- scales::rescale(
  c(
    min_val,
    q90,
    q93,
    q95,
    q99
  ),
  to = c(0, 1)
)

p_wnt <- FeaturePlot(
  scRNA,
  features = "WNT_Activity",
  reduction = "tsne",
  pt.size = 0.1,
  order = TRUE,
  max.cutoff = "q99"
) +
  scale_color_gradientn(
    colors = c(
      "grey90",
      "grey90",
      "#f5929d",
      "#f47983",
      "#f05b72"
    ),
    values = q_values,
    limits = c(min_val, q99),
    oob = scales::squish
  ) +
  labs(
    title = "WntQuant"
  ) +
  theme_classic(
    base_size = 14
  )

ggsave(
  "scRNA_WntQuant_tSNE.pdf",
  p_wnt,
  width = 6,
  height = 5
)
features_to_plot <- c(
  "WNT_Pos",
  "WNT_Neg"
)

color_list <- list(
  WNT_Pos = c(
    "grey90",
    "grey90",
    "#ffadad",
    "#ff7979",
    "#eb4d4b"
  ),
  WNT_Neg = c(
    "grey90",
    "grey90",
    "#a2c2e6",
    "#73a5d8",
    "#4589c9"
  )
)

for (feature in features_to_plot) {
  
  scores <- scRNA[[feature]][, 1]
  scores <- scores[!is.na(scores)]
  
  min_val <- min(scores)
  q90 <- quantile(scores, 0.90)
  q93 <- quantile(scores, 0.93)
  q95 <- quantile(scores, 0.95)
  q99 <- quantile(scores, 0.99)
  
  q_values <- scales::rescale(
    c(
      min_val,
      q90,
      q93,
      q95,
      q99
    ),
    to = c(0, 1)
  )
  
  p <- FeaturePlot(
    scRNA,
    features = feature,
    reduction = "tsne",
    pt.size = 0.1,
    order = TRUE,
    max.cutoff = "q99"
  ) +
    scale_color_gradientn(
      colors = color_list[[feature]],
      values = q_values,
      limits = c(min_val, q99),
      oob = scales::squish
    ) +
    theme_classic(
      base_size = 14
    )
  
  ggsave(
    paste0(
      "scRNA_",
      feature,
      "_tSNE.pdf"
    ),
    p,
    width = 6,
    height = 5
  )
}
p_gene <- DotPlot(
  scRNA,
  features = target_genes,
  group.by = "iCMS_msi"
) +
  scale_color_gradientn(
    colors = c(
      "grey90",
      "#f7acbc",
      "#f47983",
      "#f05b72"
    )
  ) +
  RotatedAxis() +
  theme_classic(
    base_size = 14
  )

ggsave(
  "DotPlot_Wnt_Target_Genes.pdf",
  p_gene,
  width = 8,
  height = 5
)

wnt_features <- c(
  "WNT_Activity",
  "WNT_Pos",
  "WNT_Neg"
)

p_score <- DotPlot(
  scRNA,
  features = wnt_features,
  group.by = "iCMS_msi"
) +
  scale_color_gradientn(
    colors = c(
      "grey90",
      "#DDA0DD",
      "#9370DB",
      "#4B0082"
    )
  ) +
  RotatedAxis() +
  theme_classic(
    base_size = 14
  )

ggsave(
  "DotPlot_WntQuant_Scores.pdf",
  p_score,
  width = 7,
  height = 5
)
saveRDS(
  scRNA,
  file = "scRNA_WntQuant_final.rds"
)

saveRDS(
  WNT_Score,
  file = "WntQuant_scores.rds"
)

#spatial
# =======================================================
# Spatial transcriptomics analysis for LUAD1
# Fang AST + WntQuant + conventional Wnt signatures
# =======================================================

library(Seurat)
library(GSEABase)
library(GSVA)
library(UCell)
library(Matrix)

# =======================================================
# 1. Paths
# =======================================================

sample_name <- "LUAD1"

spatial_path <- "/home/user/Fanglab1/ymh/LUAD_LUSC/one_lung_adenocarcinoma/AXB-2431-A1-20241204/square_008um"
work_path <- "/home/user/Fanglab1/ymh/dk/"

marker_file <- file.path(work_path, "marker.Rdata")
wntquant_script <- file.path(work_path, "Calculate_bioactivity_scores.R")
wntquant_gmt <- file.path(work_path, "wpags_wpigs.gmt")
usual_wnt_gmt <- file.path(work_path, "wnt_usual_using.gmt")

setwd(work_path)

# =======================================================
# 2. Load Visium HD sample
# =======================================================

load_sample <- function(sample_name, base_dir) {
  matrix_dir <- file.path(base_dir, "filtered_feature_bc_matrix")
  spatial_dir <- file.path(base_dir, "spatial")
  
  counts <- Read10X(data.dir = matrix_dir)
  image_obj <- Read10X_Image(
    image.dir = spatial_dir,
    filter.matrix = TRUE
  )
  
  obj <- CreateSeuratObject(
    counts = counts,
    assay = "Spatial"
  )
  
  image_obj <- image_obj[Cells(obj)]
  DefaultAssay(image_obj) <- "Spatial"
  obj[[sample_name]] <- image_obj
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = "^MT-"
  )
  
  obj$orig.ident <- sample_name
  return(obj)
}

luad1 <- load_sample(
  sample_name = sample_name,
  base_dir = spatial_path
)

DefaultAssay(luad1) <- "Spatial"

# =======================================================
# 3. Normalize expression
# =======================================================

luad1 <- NormalizeData(
  luad1,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# =======================================================
# 4. Load Fang ADC/SCC signatures from marker.Rdata
# =======================================================

extract_fang_signatures <- function(marker_file) {
  loaded_objects <- load(marker_file)
  
  if ("markers" %in% loaded_objects && is.list(get("markers"))) {
    marker_list <- get("markers")
  } else {
    list_objects <- loaded_objects[
      vapply(
        loaded_objects,
        function(x) is.list(get(x)),
        logical(1)
      )
    ]
    
    if (length(list_objects) == 1) {
      marker_list <- get(list_objects[1])
    } else {
      marker_list <- NULL
    }
  }
  
  if (!is.null(marker_list)) {
    marker_names <- names(marker_list)
    
    adc_match <- grep(
      "ADC|LUAD|ADENO",
      marker_names,
      ignore.case = TRUE,
      value = TRUE
    )
    
    scc_match <- grep(
      "SCC|LUSC|SQUAM",
      marker_names,
      ignore.case = TRUE,
      value = TRUE
    )
    
    if (length(adc_match) == 0 || length(scc_match) == 0) {
      stop("Cannot identify ADC/LUAD and SCC/LUSC signatures in marker.Rdata.")
    }
    
    adc_genes <- unique(as.character(marker_list[[adc_match[1]]]))
    scc_genes <- unique(as.character(marker_list[[scc_match[1]]]))
  } else {
    adc_objects <- grep(
      "ADC|LUAD|ADENO",
      loaded_objects,
      ignore.case = TRUE,
      value = TRUE
    )
    
    scc_objects <- grep(
      "SCC|LUSC|SQUAM",
      loaded_objects,
      ignore.case = TRUE,
      value = TRUE
    )
    
    if (length(adc_objects) == 0 || length(scc_objects) == 0) {
      stop("Cannot identify ADC/LUAD and SCC/LUSC signature objects in marker.Rdata.")
    }
    
    adc_genes <- unique(as.character(get(adc_objects[1])))
    scc_genes <- unique(as.character(get(scc_objects[1])))
  }
  
  list(
    ADC = adc_genes,
    SCC = scc_genes
  )
}

fang_signatures <- extract_fang_signatures(marker_file)

fang_signatures$ADC <- intersect(
  fang_signatures$ADC,
  rownames(luad1)
)

fang_signatures$SCC <- intersect(
  fang_signatures$SCC,
  rownames(luad1)
)

if (length(fang_signatures$ADC) < 2) {
  stop("Too few Fang ADC genes are present in the spatial expression matrix.")
}

if (length(fang_signatures$SCC) < 2) {
  stop("Too few Fang SCC genes are present in the spatial expression matrix.")
}

# =======================================================
# 5. Calculate UCell lineage scores
# =======================================================

ucell_signatures <- list(
  ADC_Fang = fang_signatures$ADC,
  SCC_Fang = fang_signatures$SCC
)

luad1 <- AddModuleScore_UCell(
  luad1,
  features = ucell_signatures
)

# =======================================================
# 6. Calculate Fang-style AST score
# =======================================================

mean_standardized_signature_score <- function(expression_matrix, genes) {
  genes <- intersect(genes, rownames(expression_matrix))
  mat <- expression_matrix[genes, , drop = FALSE]
  
  n <- ncol(mat)
  gene_sum <- Matrix::rowSums(mat)
  gene_sq_sum <- Matrix::rowSums(mat ^ 2)
  gene_mean <- gene_sum / n
  
  gene_var <- (
    gene_sq_sum - n * gene_mean ^ 2
  ) / max(n - 1, 1)
  
  gene_sd <- sqrt(
    pmax(as.numeric(gene_var), 0)
  )
  
  keep <- is.finite(gene_sd) & gene_sd > 0
  
  mat <- mat[keep, , drop = FALSE]
  gene_mean <- gene_mean[keep]
  gene_sd <- gene_sd[keep]
  
  if (nrow(mat) == 0) {
    stop("No variable genes remain for signature scoring.")
  }
  
  scaled_sparse <- Diagonal(x = 1 / gene_sd) %*% mat
  
  score <- Matrix::colMeans(scaled_sparse) - mean(gene_mean / gene_sd)
  names(score) <- colnames(mat)
  
  return(score)
}

exp_all <- GetAssayData(
  luad1,
  assay = "Spatial",
  layer = "data"
)

adc_score <- mean_standardized_signature_score(
  expression_matrix = exp_all,
  genes = fang_signatures$ADC
)

scc_score <- mean_standardized_signature_score(
  expression_matrix = exp_all,
  genes = fang_signatures$SCC
)

ast_score <- scc_score - adc_score

luad1$ADC_Scaled <- adc_score[colnames(luad1)]
luad1$SCC_Scaled <- scc_score[colnames(luad1)]
luad1$AST_Score <- ast_score[colnames(luad1)]

fang_ast_df <- data.frame(
  Barcode = colnames(luad1),
  ADC_Scaled = luad1$ADC_Scaled,
  SCC_Scaled = luad1$SCC_Scaled,
  AST_Score = luad1$AST_Score
)

write.csv(
  fang_ast_df,
  file = file.path(work_path, "LUAD1_Fang_AST_Scores.csv"),
  row.names = FALSE
)

# =======================================================
# 7. Retain bins with ADC or SCC lineage signal
# =======================================================

luad1_analysis <- subset(
  luad1,
  subset = ADC_Fang_UCell > 0 | SCC_Fang_UCell > 0
)

exp_matrix <- GetAssayData(
  luad1_analysis,
  assay = "Spatial",
  layer = "data"
)

# =======================================================
# 8. Calculate WntQuant activity
# =======================================================

source(wntquant_script)
geneSets <- getGmt(wntquant_gmt)

WNT_Score <- Calculate_bioactivity_scores(
  file_paths = work_path,
  expression_profile = exp_matrix,
  foundation = "relative_ssGSEA",
  activation_geneset = NA,
  inhibition_geneset = NA,
  geneSets_gmt = geneSets,
  min.sz = 1,
  max.sz = 10000,
  geneset_direction = c(
    rep("pos", sum(grepl("WPAGS", names(geneSets)))),
    rep("neg", sum(grepl("WPIGS", names(geneSets))))
  )
)

# =======================================================
# 9. Calculate conventional Wnt signature scores
# =======================================================

geneSets2 <- getGmt(usual_wnt_gmt)

param <- ssgseaParam(
  exp_matrix,
  geneSets2,
  minSize = 2,
  maxSize = 10000
)

raw_ssgsea <- gsva(
  param,
  verbose = FALSE
)

scaled_ssgsea_df <- as.data.frame(
  scale(
    as.data.frame(
      t(raw_ssgsea)
    )
  )
)

# =======================================================
# 10. Export all spatial scores
# =======================================================

barcodes <- colnames(luad1_analysis)
coords <- GetTissueCoordinates(luad1_analysis)
coords <- coords[barcodes, , drop = FALSE]

if (all(c("x", "y") %in% colnames(coords))) {
  x_coord <- coords$x
  y_coord <- coords$y
} else {
  x_coord <- coords[, 1]
  y_coord <- coords[, 2]
}

result_df <- data.frame(
  Barcode = barcodes,
  x = x_coord,
  y = y_coord,
  ADC_Fang_UCell = luad1_analysis$ADC_Fang_UCell,
  SCC_Fang_UCell = luad1_analysis$SCC_Fang_UCell,
  ADC_Scaled = luad1_analysis$ADC_Scaled,
  SCC_Scaled = luad1_analysis$SCC_Scaled,
  AST_Score = luad1_analysis$AST_Score,
  WPAGS_WPIGS_Activity = WNT_Score$final_activity_score$activity_score,
  WPAGS_WPIGS_Pos = WNT_Score$final_activity_score$pos_score,
  WPAGS_WPIGS_Neg = WNT_Score$final_activity_score$neg_score,
  check.names = FALSE
)

scaled_ssgsea_df <- scaled_ssgsea_df[barcodes, , drop = FALSE]
result_df <- cbind(result_df, scaled_ssgsea_df)

write.csv(
  result_df,
  file = file.path(work_path, "LUAD1_All_Scores.csv"),
  row.names = FALSE
)

saveRDS(
  luad1,
  file = file.path(work_path, "LUAD1_VisiumHD_Fang_AST.rds")
)

saveRDS(
  luad1_analysis,
  file = file.path(work_path, "LUAD1_VisiumHD_Analysis_Filtered.rds")
)

rm(
  exp_all,
  exp_matrix,
  raw_ssgsea,
  scaled_ssgsea_df,
  result_df,
  coords
)

gc()

message(
  "Finished. Outputs:\n",
  file.path(work_path, "LUAD1_Fang_AST_Scores.csv"), "\n",
  file.path(work_path, "LUAD1_All_Scores.csv")
)