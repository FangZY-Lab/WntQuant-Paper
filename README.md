# WntQuant-Paper

![WntQuant-Paper overview](https://github.com/user-attachments/assets/216f94f6-23ce-4a21-a294-3a39412a9267)

## Overview

WntQuant-Paper is an end-to-end, multi-omics analytical framework for quantifying
Wnt/β-catenin pathway activity and translating it into an interpretable,
clinically actionable biomarker. It implements a **context-insensitive and
directionally resolved activity score** that separates pathway activation from
inhibition, enabling robust, cross-cohort measurement of Wnt signaling at
multiple biological resolutions.

## Highlights

- **Directional activity quantification.** WntQuant-Paper models activation
  (WPAGS, Wnt Pathway Activation Gene Sets) and inhibition (WPIGS, Wnt Pathway
  Inhibition Gene Sets) separately and combines them into a single signed
  activity score, preserving the biological direction of Wnt signaling.
- **Context-insensitive design.** The score is derived and validated across
  dozens of independent cohorts and cancer types, minimizing batch- and
  context-specific bias for a generalizable measurement.
- **Multi-cohort, multi-omic integration.** Bulk transcriptomes, pan-cancer
  atlases (TCGA and TARGET), single-cell RNA-seq, spatial transcriptomics, and
  metabolic flux are unified into one coherent framework.
- **Systematic downstream characterization.** The framework links Wnt activity
  to molecular subtypes, immune and stemness phenotypes, and metabolic
  reprogramming, yielding a multidimensional functional landscape.
- **Translational potential.** Pan-cancer subtyping and single-cell/spatial
  resolution position the Wnt activity score as a candidate predictive
  biomarker for Wnt-directed therapies and precision oncology.

## Core method

1. **Directional gene sets.** Two curated, non-redundant gene sets are used:
   `WPAGS` (activation / positive direction) and `WPIGS` (inhibition / negative
   direction).
2. **Activity scoring.** Single-sample enrichment (relative ssGSEA) is applied
   per direction, then the signed directions are combined into the final Wnt
   activity score.
3. **Quality control.** PCA-based inspection confirms group separation in each
   training dataset before score computation.
4. **Validation.** The score is benchmarked in large colorectal-cancer cohorts
   and pan-cancer atlases, then dissected by subtype, immunity, stemness, and
   metabolism.
5. **Resolution scaling.** The same score is projected onto single-cell and
   spatial data for cellular and spatial interpretation.

## Repository structure

The analyses are organized as a numbered, ordered pipeline:

| Script | Purpose |
| --- | --- |
| `WntQuant-Paper1_integrated_Wnt_signature_construction.R` | PCA-based training-data quality control and construction of the integrated directional Wnt signature. |
| `WntQuant-Paper2_analyses_in_colorectal_cancer.R` | Large-scale validation in colorectal cancer, including CMS molecular-subtype stratification. |
| `WntQuant-Paper3_initial_explorations_in_pan_cancer.R` | Pan-cancer (TCGA and TARGET) Wnt activity scoring and initial cross-cancer comparisons. |
| `WntQuant-Paper4_deep_explorations_in_pan_cancer.R` | Pan-cancer subtyping (k-means clustering) and deep multi-omic characterization. |
| `WntQuant-Paper5_immunological_and_stemness_analyses.R` | Immunological subtypes, MSI status, and stemness analyses. |
| `WntQuant-Paper6_metabolic_flux_analyses.R` | Metabolic flux analysis linking Wnt activity to metabolic reprogramming. |
| `WntQuant-Paper7_single_cell_spatial_analyses.R` | Single-cell RNA-seq and spatial analyses at cellular resolution. |
| `WntQuant-Paper8_spatial_transcriptomics_visualization.py` | Visium HD spatial transcriptomics visualization and statistics. |

### Revision materials

`revision/` holds the code and supplementary tables produced during revision.
It is kept separate from the numbered pipeline because these scripts extend
and re-run the analyses above rather than adding new pipeline stages. The
file names follow the same convention as the root scripts.

| File | Purpose |
| --- | --- |
| `revision/WntQuant-Revision1_wnt_subtype_comparison_and_validation.R` | Benchmarks the Wnt signature against published signatures, runs fGSEA, validates the subtypes in large CMS and iCMS colorectal cohorts, assembles the TCGA pan-cancer clinical table, and renders the multivariable OS/DSS/DFI/PFI Cox forest plots. |
| `revision/WntQuant-Revision2_per_cohort_immunotherapy_analyses.R` | Runs the classification and survival workflow cohort by cohort across the public immunotherapy datasets. |
| `revision/WntQuant-Revision3_immunotherapy_response_and_prognosis.R` | Relates Wnt activity to response (R vs NR) and to OS/PFS across immunotherapy cohorts, including ESTIMATE-adjusted analyses. |
| `revision/WntQuant-Revision4_metabolic_analysis_and_robustness_tests.R` | Links Wnt activity to TCGA metabolic programmes and quantifies robustness to gene and cohort dropout and to added noise. |
| `revision/WntQuant-Revision5_pan_cancer_wnt_score_comparison.R` | Compares Wnt activity across the TCGA and TARGET pan-cancer cohorts and correlates it with immune, stemness and Hallmark programmes. |
| `revision/WntQuant-Revision6_revised_immunotherapy_prognostic_analysis.R` | Re-runs the immunotherapy response analyses with the revised Hot/Cold assignment and per-cohort KNN + surv_cutpoint evaluation. |
| `revision/WntQuant-Revision7_wss_and_panwss_subtype_definition.R` | Defines the transferable WSS and panWSS subtypes, tests their internal validity, and projects them onto TCGA, TARGET and additional cohorts. |
| `revision/WntQuant-Revision8_spatial_visium_hd_luad1_workflow.py` | Reads one Visium HD sample, runs spot-level QC, projects the Scissor, AST and WntQuant scores onto the tissue, and exports the per-feature PDFs plus the LUAD-like vs LUSC-like pathway comparison. |
| `revision/supplementary material/` | Five workbooks: Wnt pathway gene sets, biological signatures and hallmarks, immune infiltration and stemness signatures, immunotherapy signatures, and the selected Wnt activity / iCMS signatures. |

The scripts were written against the authors' local layout and point at an
absolute `input` directory, so adjust those paths (and the working directory)
before running them elsewhere.

## Data requirements

The scripts read expression matrices, survival/group annotations, and the
directional GMT file (`wpags_wpigs.gmt`) from a local `input` directory. Raw
datasets are obtained from public repositories (GEO, TCGA, TARGET, etc.) and
are not redistributed here. Paths in the scripts are written relative to the
project root; adjust the working directory to your local environment.

## How to run

Run the scripts in numerical order, adjusting input paths to your local setup:

```r
source("WntQuant-Paper1_integrated_Wnt_signature_construction.R")
source("WntQuant-Paper2_analyses_in_colorectal_cancer.R")
# ... continue through WntQuant-Paper7
```

Spatial visualization is run with Python:

```bash
python WntQuant-Paper8_spatial_transcriptomics_visualization.py
```

## Notes

Part of the code involves the core interests of the laboratory and is not open
source at this time; the corresponding sections are clearly marked at the top
of each script.


## Author

Dingkang Zhao (赵定康)

Email: <dingkang.25@intl.zju.edu.cn>
