# =======================================================
# LUAD1 Visium HD spatial analysis
# Fang AST + Scissor + WntQuant visualization/statistics
# =======================================================

from pathlib import Path
from PIL import Image
import json

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.stats import ranksums

# =======================================================
# 1. Global settings and paths
# =======================================================

sc.settings.verbosity = 0

plt.rcParams.update({
  "figure.dpi": 300,
  "font.size": 14,
  "font.family": "sans-serif",
  "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
  "text.color": "black",
  "axes.labelcolor": "black",
  "xtick.color": "black",
  "ytick.color": "black",
  "font.weight": "normal",
  "axes.titleweight": "normal"
})

sample = "LUAD1"

base_path = Path(
  "/home/user/Fanglab1/ymh/LUAD_LUSC/"
  "one_lung_adenocarcinoma/"
  "AXB-2431-A1-20241204/"
  "square_008um"
)

work_path = Path("/home/user/Fanglab1/ymh/dk/")

fang_ast_csv = work_path / "LUAD1_Fang_AST_Scores.csv"
all_scores_csv = work_path / "LUAD1_All_Scores.csv"
scissor_csv = work_path / "luad1_Scissor_AllSpots_for_Python.csv"

# =======================================================
# 2. Load one Visium HD sample
# =======================================================

def load_visium_hd_sample(base_path, library_id):
  adata = sc.read_10x_mtx(
    base_path / "filtered_feature_bc_matrix",
    var_names="gene_symbols",
    cache=True
  )

adata.var_names_make_unique()
adata.obs_names = library_id + "_" + adata.obs_names.astype(str)

positions_path = base_path / "spatial" / "tissue_positions.parquet"
positions = pd.read_parquet(positions_path, engine="pyarrow")

if "barcode" in positions.columns:
  positions = positions.set_index("barcode")
else:
  positions = positions.set_index(positions.columns[0])

positions.index = library_id + "_" + positions.index.astype(str)

common = adata.obs_names.intersection(positions.index)
adata = adata[common].copy()
positions = positions.loc[common].copy()

json_path = base_path / "spatial" / "scalefactors_json.json"
with open(json_path, "r") as f:
  scalefactors = json.load(f)

scale_factor = scalefactors["tissue_hires_scalef"]

positions["imagecol"] = positions["pxl_col_in_fullres"] * scale_factor
positions["imagerow"] = positions["pxl_row_in_fullres"] * scale_factor
positions["sample"] = library_id

adata.obs = positions.copy()
adata.obsm["spatial"] = positions[
  ["pxl_col_in_fullres", "pxl_row_in_fullres"]
].to_numpy()

image_path = base_path / "spatial" / "tissue_hires_image.png"
image = np.asarray(Image.open(image_path)).astype(np.float32) / 255.0

adata.uns["spatial"] = {
  library_id: {
    "images": {"hires": image},
    "scalefactors": scalefactors
  }
}

return adata


adata = load_visium_hd_sample(
  base_path=base_path,
  library_id=sample
)

# =======================================================
# 3. Inject Fang AST scores
# =======================================================

fang_df = pd.read_csv(fang_ast_csv)
fang_df["Barcode_full"] = sample + "_" + fang_df["Barcode"].astype(str)
fang_df = fang_df.set_index("Barcode_full")

fang_cols = ["ADC_Scaled", "SCC_Scaled", "AST_Score"]

adata.obs = adata.obs.drop(
  columns=[c for c in fang_cols if c in adata.obs.columns],
  errors="ignore"
)

adata.obs = adata.obs.join(
  fang_df[fang_cols],
  how="left"
)

# =======================================================
# 4. Inject WntQuant and conventional Wnt scores
# =======================================================

scores_df = pd.read_csv(all_scores_csv)
scores_df["Barcode_full"] = sample + "_" + scores_df["Barcode"].astype(str)
scores_df = scores_df.set_index("Barcode_full")

excluded_columns = {"Barcode", "x", "y"}
score_columns = [c for c in scores_df.columns if c not in excluded_columns]

duplicate_columns = [c for c in score_columns if c in adata.obs.columns]
if duplicate_columns:
  adata.obs = adata.obs.drop(columns=duplicate_columns)

adata.obs = adata.obs.join(
  scores_df[score_columns],
  how="left"
)

# =======================================================
# 5. Inject Scissor labels
# =======================================================

scissor_df = pd.read_csv(scissor_csv, index_col=0)

if not str(scissor_df.index[0]).startswith(f"{sample}_"):
  scissor_df.index = sample + "_" + scissor_df.index.astype(str)

scissor_columns = [
  c for c in ["Scissor_Label", "Scissor_Coef"]
  if c in scissor_df.columns
]

adata.obs = adata.obs.drop(
  columns=[c for c in scissor_columns if c in adata.obs.columns],
  errors="ignore"
)

adata.obs = adata.obs.join(
  scissor_df[scissor_columns],
  how="left"
)

adata.obs["Scissor_Label"] = adata.obs["Scissor_Label"].fillna("Background")

# =======================================================
# 6. Prepare spatial coordinates and grayscale H&E
# =======================================================

library_info = adata.uns["spatial"][sample]
scale_factor = library_info["scalefactors"]["tissue_hires_scalef"]

x_coords = adata.obsm["spatial"][:, 0] * scale_factor
y_coords = adata.obsm["spatial"][:, 1] * scale_factor

img = library_info["images"]["hires"].copy()
gray = np.dot(img[..., :3], [0.299, 0.587, 0.114])
gray_3ch = np.stack([gray, gray, gray], axis=-1)

adata.uns["spatial"][sample]["images"]["hires_gray"] = gray_3ch.astype(np.float32)
adata.uns["spatial"][sample]["scalefactors"]["tissue_hires_gray_scalef"] = scale_factor

# =======================================================
# 7. Color maps
# =======================================================

color_red = "#d71345"
color_blue = "#009ad6"
color_white = "#ffffff"

cmap_blue = mcolors.LinearSegmentedColormap.from_list(
  "ADC_blue", [color_white, color_blue]
)

cmap_red = mcolors.LinearSegmentedColormap.from_list(
  "SCC_red", [color_white, color_red]
)

cmap_ast = mcolors.LinearSegmentedColormap.from_list(
  "AST_blue_white_red", [color_blue, color_white, color_red]
)

cmap_activity = mcolors.LinearSegmentedColormap.from_list(
  "Wnt_blue_white_red", [color_blue, color_white, color_red]
)

scissor_palette = {
  "LUSC_like": color_red,
  "LUAD_like": color_blue,
  "Scissor+": color_red,
  "Scissor-": color_blue,
  "Background": "#e0e0e0"
}

# =======================================================
# 8. Plot Scissor classification
# =======================================================

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

sc.pl.spatial(
  adata,
  color="Scissor_Label",
  library_id=sample,
  palette=scissor_palette,
  alpha_img=0.5,
  size=1.2,
  show=False,
  title="Scissor phenotype mapping",
  ax=axes[0]
)

if "Scissor_Coef" in adata.obs.columns:
  sc.pl.spatial(
    adata,
    color="Scissor_Coef",
    library_id=sample,
    cmap=cmap_ast,
    vcenter=0,
    alpha_img=0.5,
    size=1.2,
    show=False,
    title="Scissor coefficient",
    ax=axes[1]
  )

for ax in axes:
  ax.axis("off")
sns.kdeplot(
  x=x_coords,
  y=y_coords,
  ax=ax,
  levels=[0.01],
  color="#555555",
  linestyles="--",
  linewidths=1.5
)

plt.tight_layout()
plt.savefig(
  work_path / "LUAD1_Scissor_Phenotype_Mapping.pdf",
  dpi=300,
  bbox_inches="tight"
)
plt.close(fig)

# =======================================================
# 9. Plot Fang ADC, SCC and AST scores
# =======================================================

fig, axes = plt.subplots(1, 3, figsize=(24, 7))

sc.pl.spatial(
  adata,
  color="ADC_Scaled",
  library_id=sample,
  cmap=cmap_blue,
  alpha_img=0.6,
  size=1.2,
  show=False,
  title="Fang ADC score",
  ax=axes[0]
)

sc.pl.spatial(
  adata,
  color="SCC_Scaled",
  library_id=sample,
  cmap=cmap_red,
  alpha_img=0.6,
  size=1.2,
  show=False,
  title="Fang SCC score",
  ax=axes[1]
)

sc.pl.spatial(
  adata,
  color="AST_Score",
  library_id=sample,
  cmap=cmap_ast,
  vcenter=0,
  alpha_img=0.6,
  size=1.2,
  show=False,
  title="Fang AST score (SCC - ADC)",
  ax=axes[2]
)

for ax in axes:
  ax.axis("off")
sns.kdeplot(
  x=x_coords,
  y=y_coords,
  ax=ax,
  levels=[0.01],
  color="#333333",
  linestyles="--",
  linewidths=1.5
)

plt.tight_layout()
plt.savefig(
  work_path / "LUAD1_Fang_AST_Spatial.pdf",
  dpi=300,
  bbox_inches="tight"
)
plt.close(fig)

# =======================================================
# 10. Plot WntQuant activity
# =======================================================

activity_column = "WPAGS_WPIGS_Activity"

if activity_column in adata.obs.columns:
  valid_values = adata.obs[activity_column].dropna()

if len(valid_values) > 0:
  vmax = np.nanquantile(np.abs(valid_values), 0.99)

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

common_args = {
  "library_id": sample,
  "cmap": cmap_activity,
  "vmin": -vmax,
  "vmax": vmax,
  "vcenter": 0,
  "size": 1.3,
  "show": False
}

sc.pl.spatial(
  adata,
  color=activity_column,
  img_key="hires_gray",
  alpha_img=0.6,
  alpha=0.8,
  title="WntQuant activity",
  ax=axes[0],
  **common_args
)

sc.pl.spatial(
  adata,
  color=activity_column,
  alpha_img=0.0,
  alpha=0.9,
  title="WntQuant activity",
  ax=axes[1],
  **common_args
)

for ax in axes:
  ax.axis("off")
sns.kdeplot(
  x=x_coords,
  y=y_coords,
  ax=ax,
  levels=[0.01],
  color="#333333",
  linestyles="--",
  linewidths=1.5
)

plt.tight_layout()
plt.savefig(
  work_path / "LUAD1_WntQuant_Spatial.pdf",
  dpi=300,
  bbox_inches="tight"
)
plt.close(fig)

# =======================================================
# 11. Compare LUAD-like and LUSC-like regions
# =======================================================

comparison_df = adata.obs[
  adata.obs["Scissor_Label"].isin(["LUAD_like", "LUSC_like"])
].copy()

analysis_columns = [
  c for c in score_columns
  if c in comparison_df.columns
  and pd.api.types.is_numeric_dtype(comparison_df[c])
]

results = []

for pathway in analysis_columns:
  data_luad = comparison_df.loc[
    comparison_df["Scissor_Label"] == "LUAD_like",
    pathway
  ].dropna()

data_lusc = comparison_df.loc[
  comparison_df["Scissor_Label"] == "LUSC_like",
  pathway
].dropna()

if len(data_luad) > 5 and len(data_lusc) > 5:
  stat, p_value = ranksums(data_lusc, data_luad)

mean_luad = data_luad.mean()
mean_lusc = data_lusc.mean()
median_luad = data_luad.median()
median_lusc = data_lusc.median()

results.append({
  "Pathway": pathway,
  "Mean_LUAD": mean_luad,
  "Mean_LUSC": mean_lusc,
  "Mean_Diff(LUSC-LUAD)": mean_lusc - mean_luad,
  "Median_Diff(LUSC-LUAD)": median_lusc - median_luad,
  "P_Value": p_value,
  "Enrichment": "LUSC_like" if mean_lusc > mean_luad else "LUAD_like"
})

stat_df = pd.DataFrame(results).sort_values("P_Value")

stat_df.to_csv(
  work_path / "LUAD1_Pathway_Wilcoxon_Stats.csv",
  index=False
)

# =======================================================
# 12. Plot directional significance
# =======================================================

if len(stat_df) > 0:
  df_plot = stat_df.copy()
df_plot = df_plot[
  ~df_plot["Pathway"].str.contains("UCell", case=False, na=False)
]

p_min_cap = 1e-200
df_plot["log10P"] = -np.log10(
  df_plot["P_Value"].clip(lower=p_min_cap)
)

df_plot["Directional_LogP"] = np.where(
  df_plot["Mean_Diff(LUSC-LUAD)"] > 0,
  df_plot["log10P"],
  -df_plot["log10P"]
)

df_plot = df_plot.sort_values(
  "Directional_LogP",
  ascending=False
)

effect_values = df_plot["Median_Diff(LUSC-LUAD)"]
max_abs_effect = max(np.nanmax(np.abs(effect_values)), 1e-6)

norm = mcolors.TwoSlopeNorm(
  vmin=-max_abs_effect,
  vcenter=0,
  vmax=max_abs_effect
)

diff_cmap = mcolors.LinearSegmentedColormap.from_list(
  "diff_cmap",
  [color_blue, "#f0f0f0", color_red]
)

row_colors = [
  diff_cmap(norm(value))
  for value in effect_values
]

fig_height = max(6, 0.38 * len(df_plot))
fig, ax = plt.subplots(figsize=(8, fig_height))

ax.barh(
  df_plot["Pathway"],
  df_plot["Directional_LogP"],
  color=row_colors,
  edgecolor="black",
  linewidth=0.6
)

ax.axvline(0, color="black", linewidth=1.2)
ax.invert_yaxis()
ax.set_xlabel("LUAD-like  <---  -log10(P)  --->  LUSC-like")
ax.set_ylabel("")

for spine in ax.spines.values():
  spine.set_visible(True)

sm = cm.ScalarMappable(cmap=diff_cmap, norm=norm)
sm.set_array([])

cbar = fig.colorbar(
  sm,
  ax=ax,
  shrink=0.5,
  aspect=20,
  pad=0.05
)

cbar.set_label(
  "Median difference\n(LUSC-like - LUAD-like)",
  rotation=270,
  labelpad=28
)

cbar.set_ticks([-max_abs_effect, 0, max_abs_effect])

plt.tight_layout()
plt.savefig(
  work_path / "LUAD1_Pathway_Directional_Barplot.pdf",
  dpi=300,
  bbox_inches="tight"
)
plt.close(fig)

# =======================================================
# 13. Save final AnnData object
# =======================================================

adata.write_h5ad(
  work_path / "LUAD1_VisiumHD_Fang_AST.h5ad"
)

print(
  "Finished:\n"
  f"{work_path / 'LUAD1_Fang_AST_Spatial.pdf'}\n"
  f"{work_path / 'LUAD1_WntQuant_Spatial.pdf'}\n"
  f"{work_path / 'LUAD1_Pathway_Wilcoxon_Stats.csv'}"
)