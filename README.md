# tnbc-kynurenine-pipeline

Reproducible R pipeline for:
**"Kynurenine-Centered Metabolic-Immune Trap Characterizes a Cold Tumor Microenvironment in Triple-Negative Breast Cancer"**

Every quantitative result in the manuscript is regenerated end-to-end by `run_all.R`
(19 scripts) from the raw inputs listed below. No manual intervention, no hand-edited
intermediate tables. Archived on Zenodo (DOI: 10.5281/zenodo.22747180).

## Pipeline structure (strictly ordered, one-way data flow)

    R/00_config.R            paths | frozen parameters (seed = 42; k = 3) | 16 pathway
                             gene sets | utility functions | integrity checks
    R/01_data.R              TCGA expression + survival -> patient-level cohort (n = 123)
    R/02_scores.R            16 pathway scores (ssGSEA for 4 metabolic pathways; mean-z
                             for the 12 immune/state pathways) + Kynurenine cutoff
    R/03_cluster.R           ConsensusClusterPlus (k-means, euclidean, 100 reps,
                             pItem = 0.8) + PAC/silhouette/gap; subtypes pre-specified
                             C1/C2/C3 by decreasing Kynurenine mean (locked before any
                             association analysis); Figure 1
    R/04_mechanism.R         Figure 2 (A-D: score-expression decoupling) + coupling panels
    R/04b_fig2_withE.R       Figure 2 five-panel version (adds E: swap-in clustering test)
    R/05_metabolic.R         Figure 3 (PGC-1a / OXPHOS axis)
    R/06_validation.R        METABRIC cohort assembly + survival objects (06_metabric.rds)
    R/07_figures.R           Figures 4 / 5 / 7 (immune desert, coupling spectrum,
                             single-cell patient-level association)
    R/08_fig9.R              Figure 8 TIDE (2x2 layout, panels A-D)
    R/09_tcga_survival.R     TCGA OS (19 events, median follow-up 35 months; log-rank p = 0.17)
    R/10_tide_supp.R         TIDE component scores (Figure 8D data)
    R/11_extra_analyses.R    Fig S1 volcano (KP genes highlighted), Fig S2 (SLC7A5
                             disconnection), Table S3
    R/12_k2_gap.R            k=2 profiles + gap statistic (Table S2)
    R/13_metabric_cluster.R  METABRIC independent subtyping (Figure S3) + OS
    R/14_sensitivity.R       Scoring-system sensitivity analysis (Table S4)
    R/15_cluster_abs_sensitivity.R  Swap-in clustering: absolute-expression z replaces
                             Kyn ssGSEA in the clustering input (original C1: 0% in the
                             absolute-high Kyn cluster)
    R/16_impress_validation.R       IMPRES consistency check for the TIDE direction (Figure S5)
    R/17_scrna_kp_heatmap.R         Per-cell KP gene expression atlas from the raw
                             GSE176078 count matrix (Figure S4; 13/14 genes, KYAT1 absent
                             from this annotation build)
    R/18_reviewer_round9.R   Robustness bundle: TIDE effect sizes (epsilon^2, Cohen's d),
                             global-background partial correlations, single-cell
                             cell-count robustness, all-ssGSEA unified-metric clustering
                             (original C1 recovered at 100%)
    R/19_final_robustness.R  14-gene C1-vs-C3 statistics with BH q-values
                             (results/19_kp_gene_stats.csv) + immune-infiltration-adjusted
                             IDO1 Cox models (METABRIC)

## Data acquisition (not redistributed; download required)

    TCGA expression    TNBC subset matrix (log2 FPKM+1) from UCSC Xena TCGA-BRCA.htseq_fpkm
    TCGA survival      Xena TCGA-BRCA.survival.tsv (OS, OS.time)
    METABRIC           cBioPortal brca_metabric study archive
    scRNA (modules)    Fig7_metadata_with_scores.csv (processed GSE176078 metadata with
                       published module scores)
    scRNA (counts)     GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz from GEO
                       (required by script 17; script searches ~/Downloads)
    TIDE               TIDEpy v1.3.9 run locally (Python >= 3.9); see below

Update paths in R/00_config.R to point to your local copies.

## Setup & run

    R >= 4.4
    install.packages(c("tidyverse", "survival", "survminer", "patchwork", "gridExtra", "Matrix"))
    BiocManager::install("ConsensusClusterPlus")
    source("run_all.R")

Seeded (params$seed = 42) and deterministic. sessionInfo.txt is written after each full run.

## TIDEpy (Figure 8)

    python3 -m venv tidepy-env && tidepy-env/bin/pip install tidepy pandas
    tidepy-env/bin/tidepy <expression.tsv> -o <output_file> -c Other

The expression file is exported from the cohort matrix (genes x samples, log2 FPKM+1);
the per-sample table (Patient, TIDE, Dysfunction, Exclusion, MDSC, CAF, TAM M2, Responder)
is the input at paths$response.

## Supplementary materials (as in the manuscript)

    Figure S1   C1 vs C3 volcano, KP genes highlighted (three downstream nodes with a
                consistent up-direction in C1; upstream/mid-pathway enzymes significantly
                lower, BH q <= 4e-4)
    Figure S2   SLC7A5 disconnection
    Figure S3   Independent subtyping in METABRIC (M1, n = 97)
    Figure S4   Per-cell KP gene expression atlas (GSE176078, n = 100,064 cells)
    Figure S5   IMPRES consistency check
    Figure S6   Dual-arm concept framework
    Table S2    Stability (PAC / silhouette / gap, k = 2-6, plus k=2 profile)
    Table S4    Sensitivity analysis (absolute-expression scores, partial correlations,
                swap-in and all-ssGSEA clustering, global-background correction)
    Table S5    Positioning vs prior kynurenine-pathway studies

## Notes

- TIDE/IMPRES are used as exploratory, cross-cancer prediction frameworks; all therapy
  stratification statements in the manuscript are directional.
- The three-subtype split is a descriptive framework; only the C1-like group (M1 in
  METABRIC) is treated as a robust entity.

## License

MIT.
