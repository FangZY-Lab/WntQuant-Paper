# ==============================================================================
# WntQuant revision -- LUAD1 Visium HD spatial workflow
# ------------------------------------------------------------------------------
# Purpose:    Read one Visium HD sample, run spot-level QC, project the Scissor,
#             AST and WntQuant scores onto the tissue, and export one PDF per
#             feature plus the LUAD-like vs LUSC-like pathway comparison.
# Inputs:     Visium HD output (filtered_feature_bc_matrix/, spatial/) and the
#             Scissor / AST / WntQuant score tables under work_path.
# Outputs:    Numbered PDF figures LUAD1_01 ... LUAD1_29 and the 18-pathway
#             Wilcoxon statistics table.
# ==============================================================================

from pathlib import Path
import json

from PIL import Image
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scanpy as sc
from scipy.stats import ranksums

# =======================================================
# 1. Global settings and paths
# =======================================================
sc.settings.verbosity = 0

# Figure style: Arial 14 pt, plain black, never bold.
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
    "axes.titleweight": "normal",
    "figure.titlesize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 14,
    "legend.title_fontsize": 14
})

sample = "LUAD1"

# Adjust these two paths to your local layout before running.
base_path = Path("/home/user/Fanglab1/ymh/aaa_lung_cancer/data/AXB-2431-A1-20241204/square_008um/")
work_path = Path("/home/user/Fanglab1/ymh/dk/")


# =======================================================
# 2. Load one Visium HD sample
# =======================================================
def load_visium_hd_sample(base_path, library_id):
    adata = sc.read_10x_mtx(
        base_path / "filtered_feature_bc_matrix",
        var_names="gene_symbols",
        cache=True,
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
    if image_path.exists():
        image = np.asarray(Image.open(image_path)).astype(np.float32) / 255.0
    else:
        image = np.zeros((100, 100, 3), dtype=np.float32)

    adata.uns["spatial"] = {
        library_id: {"images": {"hires": image}, "scalefactors": scalefactors}
    }
    return adata


adata = load_visium_hd_sample(base_path=base_path, library_id=sample)

# =======================================================
# 3. Technical QC metrics and pre-filter QC figures
# =======================================================
print("\nComputing per-spot QC metrics...")

if sp.issparse(adata.X):
    adata.obs["nUMI"] = adata.X.sum(axis=1).A1
    adata.obs["nGene"] = (adata.X > 0).sum(axis=1).A1
else:
    adata.obs["nUMI"] = adata.X.sum(axis=1)
    adata.obs["nGene"] = (adata.X > 0).sum(axis=1)

nUMI_thresh = 100
nGene_thresh = 50

fail_mask = (adata.obs["nUMI"] <= nUMI_thresh) | (adata.obs["nGene"] <= nGene_thresh)

# Spots below the threshold are drawn in deep blue; passing spots keep the
# Reds scale, so one figure shows both the gradient and the discarded spots.
color_fail = "#0b3d91"

adata.obs["nUMI_plot"] = adata.obs["nUMI"].astype(float)
adata.obs.loc[fail_mask, "nUMI_plot"] = np.nan

adata.obs["nGene_plot"] = adata.obs["nGene"].astype(float)
adata.obs.loc[fail_mask, "nGene_plot"] = np.nan

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.spatial(
    adata,
    color="nUMI_plot",
    library_id=sample,
    cmap="Reds",
    na_color=color_fail,
    size=1.2,
    alpha_img=0.1,
    title=f"nUMI (fail <= {nUMI_thresh})",
    show=False,
    ax=ax,
)
ax.axis("off")
plt.tight_layout()
plt.savefig(work_path / "LUAD1_01_QC_nUMI.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.spatial(
    adata,
    color="nGene_plot",
    library_id=sample,
    cmap="Reds",
    na_color=color_fail,
    size=1.2,
    alpha_img=0.1,
    title=f"nGene (fail <= {nGene_thresh})",
    show=False,
    ax=ax,
)
ax.axis("off")
plt.tight_layout()
plt.savefig(work_path / "LUAD1_02_QC_nGene.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

adata.obs = adata.obs.drop(columns=["nUMI_plot", "nGene_plot"])

# =======================================================
# 4. Technical hard filtering
# =======================================================
print(f"Filtering before: {adata.n_obs} spots (whole slide)")

mask_tech = (adata.obs["nUMI"] > nUMI_thresh) & (adata.obs["nGene"] > nGene_thresh)
adata = adata[mask_tech].copy()
sc.pp.filter_genes(adata, min_cells=1)

print(f"Filtering after: {adata.n_obs} spots (valid tissue retained)")

# =======================================================
# 5. Inject Scissor, AST and WntQuant scores
# =======================================================
scissor_csv = work_path / "luad1_Scissor_AllSpots_for_Python.csv"
ast_csv_path = work_path / "luad1_AST_ssGSEA_Scores.csv"
all_scores_csv = work_path / "LUAD1_All_Scores.csv"

# Scissor phenotype labels and coefficients.
if scissor_csv.exists():
    scissor_df = pd.read_csv(scissor_csv, index_col=0)
    if not str(scissor_df.index[0]).startswith(f"{sample}_"):
        scissor_df.index = sample + "_" + scissor_df.index.astype(str)
    sc_cols = [c for c in ["Scissor_Label", "Scissor_Coef"] if c in scissor_df.columns]
    adata.obs = adata.obs.drop(
        columns=[c for c in sc_cols if c in adata.obs.columns], errors="ignore"
    )
    adata.obs = adata.obs.join(scissor_df[sc_cols], how="left")
    adata.obs["Scissor_Label"] = adata.obs["Scissor_Label"].fillna("Background")

# AST (adeno-to-squamous transition) scores.
if ast_csv_path.exists():
    df_ast = pd.read_csv(ast_csv_path, index_col=0)
    if not str(df_ast.index[0]).startswith(f"{sample}_"):
        df_ast.index = sample + "_" + df_ast.index.astype(str)
    ast_cols = ["ADC_Scaled", "SCC_Scaled", "AST_Score"]
    adata.obs = adata.obs.drop(
        columns=[c for c in ast_cols if c in adata.obs.columns], errors="ignore"
    )
    adata.obs = adata.obs.join(df_ast[ast_cols], how="left")

# WntQuant activity and its component scores.
if all_scores_csv.exists():
    df_wnt = pd.read_csv(all_scores_csv)
    df_wnt["Barcode_full"] = sample + "_" + df_wnt["Barcode"].astype(str)
    df_wnt = df_wnt.set_index("Barcode_full")
    wnt_cols = [c for c in df_wnt.columns if c not in {"Barcode", "x", "y"}]
    adata.obs = adata.obs.drop(
        columns=[c for c in wnt_cols if c in adata.obs.columns], errors="ignore"
    )
    adata.obs = adata.obs.join(df_wnt[wnt_cols], how="left")

# =======================================================
# 6. Post-filter spatial maps
# =======================================================
print("\nRendering post-QC spatial maps (one PDF per feature)...")

library_id = list(adata.uns["spatial"].keys())[0]
scale_factor = adata.uns["spatial"][library_id]["scalefactors"]["tissue_hires_scalef"]

# Tissue outline drawn from the retained spots.
x_coords = adata.obsm["spatial"][:, 0] * scale_factor
y_coords = adata.obsm["spatial"][:, 1] * scale_factor

color_red, color_blue, color_white, color_bg = "#d71345", "#009ad6", "#ffffff", "#e0e0e0"
scissor_palette = {"LUSC_like": color_red, "LUAD_like": color_blue, "Background": color_bg}
cmap_diverge = mcolors.LinearSegmentedColormap.from_list(
    "BlueWhiteRed", [color_blue, color_white, color_red]
)


def outline(ax):
    """Draw the dashed tissue outline used by every spatial panel."""
    sns.kdeplot(
        x=x_coords,
        y=y_coords,
        ax=ax,
        levels=[0.01],
        color="#333333",
        linestyles="--",
        linewidths=1.5,
    )


# -- Scissor phenotype --
fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.spatial(
    adata,
    color="Scissor_Label",
    library_id=sample,
    palette=scissor_palette,
    size=1.2,
    alpha_img=0.1,
    title="Scissor Phenotypes",
    show=False,
    ax=ax,
)
ax.axis("off")
outline(ax)
plt.tight_layout()
plt.savefig(work_path / "LUAD1_03_PostQC_Scissor.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

# -- AST score --
fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.spatial(
    adata,
    color="AST_Score",
    library_id=sample,
    cmap=cmap_diverge,
    vcenter=0,
    size=1.2,
    alpha_img=0.1,
    title="AST Score (SCC - ADC)",
    show=False,
    ax=ax,
)
ax.axis("off")
outline(ax)
plt.tight_layout()
plt.savefig(work_path / "LUAD1_04_PostQC_AST.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

# -- WntQuant activity --
activity_col = "WPAGS_WPIGS_Activity"
if activity_col not in adata.obs.columns:
    adata.obs[activity_col] = 0.0

fig, ax = plt.subplots(figsize=(8, 8))
sc.pl.spatial(
    adata,
    color=activity_col,
    library_id=sample,
    cmap=cmap_diverge,
    vmin=-1.0,
    vmax=1.0,
    vcenter=0,
    size=1.2,
    alpha_img=0.1,
    title="WntQuant Activity",
    show=False,
    ax=ax,
)
ax.axis("off")
outline(ax)
plt.tight_layout()
plt.savefig(work_path / "LUAD1_05_PostQC_WntQuant.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"[*] 5 post-filter PDFs written to {work_path}")

# =======================================================
# 7. Normalisation and Wnt target-gene expression maps
# =======================================================
print("\nNormalising (target_sum = 1e4) and applying log1p...")

# Keep the raw counts for later use.
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

target_genes = ["AXIN2", "ZNRF3", "RNF43", "LGR5", "TCF7"]

print("\nRendering Wnt target-gene maps...")

for i, gene in enumerate(target_genes):
    if gene not in adata.var_names:
        print(f"[!] {gene} is not in the matrix, skipped")
        continue

    # Fraction of spots with non-zero expression, shown in the title.
    if sp.issparse(adata.X):
        expr = adata[:, gene].X.toarray().flatten()
    else:
        expr = adata[:, gene].X.flatten()

    non_zero_pct = (expr > 0).sum() / len(expr) * 100

    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.spatial(
        adata,
        color=gene,
        library_id=sample,
        cmap="Reds",
        size=1.2,
        alpha_img=0.1,
        title=f"Non-zero: {non_zero_pct:.2f}%",
        show=False,
        ax=ax,
    )
    ax.axis("off")
    ax.text(
        0.5,
        -0.05,
        f"Marker gene: {gene}",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=14,
        color="black",
    )
    outline(ax)

    plt.tight_layout()
    out_pdf = work_path / f"LUAD1_{i + 6:02d}_PostQC_Gene_{gene}.pdf"
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[*] {gene} -> {out_pdf.name} (non-zero: {non_zero_pct:.2f}%)")

# =======================================================
# 8. Wnt pathway score maps
# =======================================================
print("\nRendering pathway score maps...")

# Display name -> adata.obs column.
features_to_plot = {
    # WntQuant activation and inhibition
    "WntQuant_Activation": "WPAGS_WPIGS_Pos",
    "WntQuant_Inhibition": "WPAGS_WPIGS_Neg",
    # Adenocarcinoma vs squamous phenotype (UCell, falling back to Scaled)
    "LUAD_like": "LUAD_UCell" if "LUAD_UCell" in adata.obs.columns else "ADC_Scaled",
    "LUSC_like": "LUSC_UCell" if "LUSC_UCell" in adata.obs.columns else "SCC_Scaled",
    # 15 Wnt-related pathway scores
    "BIOCARTA_Wnt": "BIOCARTA_Wnt",
    "FEVR_Wnt_Down": "FEVR_Wnt_Down",
    "FEVR_Wnt_Up": "FEVR_Wnt_Up",
    "FLIER_Wnt": "FLIER_Wnt",
    "GOBP_Wnt": "GOBP_Wnt",
    "HALLMARK_Wnt": "HALLMARK_Wnt",
    "KEGG_Wnt": "KEGG_Wnt",
    "PID_Wnt": "PID_Wnt",
    "PROGENy_Wnt": "PROGENy_Wnt",
    "REACTOME_Wnt": "REACTOME_Wnt",
    "WP_Wnt": "WP_Wnt",
    "Geneformer_Wnt_Top100": "Geneformer_Wnt_Top100",
    "Geneformer_Wnt_Top200": "Geneformer_Wnt_Top200",
    "GenePT_Wnt_Top100": "GenePT_Wnt_Top100",
    "GenePT_Wnt_Top200": "GenePT_Wnt_Top200",
}

plot_idx = 11

for display_name, col_name in features_to_plot.items():
    # A missing column is filled with zeros so one gap does not stop the run.
    if col_name not in adata.obs.columns:
        print(f"[!] column {col_name} is missing, filled with 0.0")
        adata.obs[col_name] = 0.0

    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.spatial(
        adata,
        color=col_name,
        library_id=sample,
        cmap=cmap_diverge,
        vmin=-1.0,
        vmax=1.0,
        vcenter=0,
        size=1.2,
        alpha_img=0.1,
        title=display_name,
        show=False,
        ax=ax,
    )
    ax.axis("off")
    outline(ax)
    plt.tight_layout()

    out_pdf = work_path / f"LUAD1_{plot_idx:02d}_PostQC_{display_name}.pdf"
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[*] {display_name} -> {out_pdf.name}")
    plot_idx += 1

# =======================================================
# 9. Pathway comparison, LUAD-like vs LUSC-like spots
# =======================================================
print(f"\nTesting 18 pathways between Scissor phenotypes in {sample}...")

adata_sub = adata[adata.obs["sample"] == sample].copy()
comparison_df = adata_sub.obs[
    adata_sub.obs["Scissor_Label"].isin(["LUAD_like", "LUSC_like"])
].copy()

# Display name -> adata.obs column.
pathway_mapping = {
    "WntQuant": "WPAGS_WPIGS_Activity",
    "GenePT_Wnt_Top200": "GenePT_Wnt_Top200",
    "WP_Wnt": "WP_Wnt",
    "GenePT_Wnt_Top100": "GenePT_Wnt_Top100",
    "HALLMARK_Wnt": "HALLMARK_Wnt",
    "PID_Wnt": "PID_Wnt",
    "REACTOME_Wnt": "REACTOME_Wnt",
    "BIOCARTA_Wnt": "BIOCARTA_Wnt",
    "KEGG_Wnt": "KEGG_Wnt",
    "WntQuant_Activation": "WPAGS_WPIGS_Pos",
    "GOBP_Wnt": "GOBP_Wnt",
    "Geneformer_Wnt_Top100": "Geneformer_Wnt_Top100",
    "Geneformer_Wnt_Top200": "Geneformer_Wnt_Top200",
    "FLIER_Wnt": "FLIER_Wnt",
    "FEVR_Wnt_Up": "FEVR_Wnt_Up",
    "PROGENy_Wnt": "PROGENy_Wnt",
    "WntQuant_Inhibition": "WPAGS_WPIGS_Neg",
    "FEVR_Wnt_Down": "FEVR_Wnt_Down",
}

results = []

for display_name, col_name in pathway_mapping.items():
    if col_name not in comparison_df.columns:
        print(f"[!] column {col_name} is missing, skipped")
        continue

    data_luad = comparison_df[comparison_df["Scissor_Label"] == "LUAD_like"][col_name].dropna()
    data_lusc = comparison_df[comparison_df["Scissor_Label"] == "LUSC_like"][col_name].dropna()

    if len(data_luad) > 5 and len(data_lusc) > 5:
        stat, p_val = ranksums(data_lusc, data_luad)
        p_val = max(p_val, 1e-300)

        med_luad, med_lusc = np.median(data_luad), np.median(data_lusc)
        med_diff = med_lusc - med_luad

        # Signed significance: -log10(p) carrying the direction of the median shift.
        log_p = -np.log10(p_val)
        composite_score = log_p * np.sign(med_diff)

        results.append({
            "Pathway": display_name,
            "Median_Diff(LUSC-LUAD)": med_diff,
            "P_Value": p_val,
            "Composite_Score": composite_score,
            "Enrichment": "LUSC" if med_diff > 0 else "LUAD",
        })

stat_df = pd.DataFrame(results).sort_values("Composite_Score", ascending=False)

print("Rendering the signed bar chart...")

plt.rcParams.update({
    "figure.dpi": 300,
    "font.size": 16,
    "font.family": "sans-serif",
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "font.weight": "normal",
    "axes.titleweight": "normal",
})

fig, ax = plt.subplots(figsize=(15, 8))

x = np.arange(len(stat_df))
heights = stat_df["Composite_Score"].values
med_diffs = stat_df["Median_Diff(LUSC-LUAD)"].values

norm = mcolors.Normalize(vmin=-0.6, vmax=0.6)
cmap = mcolors.LinearSegmentedColormap.from_list(
    "BlueWhiteRed", ["#009ad6", "#ffffff", "#d71345"]
)
bar_colors = [cmap(norm(val)) for val in med_diffs]

ax.bar(x, heights, color=bar_colors, edgecolor="#333333", linewidth=0.6)
ax.axhline(0, color="black", linewidth=1.2)

ax.set_xticks(x)
ax.set_xticklabels(stat_df["Pathway"].values, rotation=90, ha="center", fontsize=14)
ax.set_ylabel("-log10(p) * sign(Median diff)", fontsize=16)

ax.text(
    -1.5, stat_df["Composite_Score"].max() * 0.6, "LUSC\n(UP)",
    ha="center", va="center", color="#d71345", fontsize=18, fontweight="bold",
)
ax.text(
    -1.5, stat_df["Composite_Score"].min() * 0.6, "LUAD\n(DOWN)",
    ha="center", va="center", color="#009ad6", fontsize=18, fontweight="bold",
)

# Colour bar for the median difference encoded by the bar colours.
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=0.02, aspect=30)
cbar.set_label("Median Difference (LUSC - LUAD)", fontsize=14, rotation=270, labelpad=20)
cbar.ax.tick_params(labelsize=12)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

out_pdf = work_path / f"{sample}_18Pathways_Wilcoxon_Barplot.pdf"
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)

out_csv = work_path / f"{sample}_18Pathways_Wilcoxon_Stats.csv"
stat_df.to_csv(out_csv, index=False)

print(f"[*] bar chart -> {out_pdf}")
print(f"[*] statistics -> {out_csv}")
