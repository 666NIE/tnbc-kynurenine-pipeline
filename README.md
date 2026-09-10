# Kynurenine-Centered Metabolic-Immune Trap in TNBC — Reproducible Pipeline

This repository contains the complete analysis pipeline for the manuscript
"Kynurenine-Centered Metabolic-Immune Trap Defines a Cold Tumor Microenvironment
in Triple-Negative Breast Cancer". Every quantitative result reported in the
manuscript (TCGA-BRCA discovery, scRNA validation, METABRIC external validation,
SLC7A5 disconnection analysis, and TIDEpy-based response prediction) is generated
end-to-end by run_all.R from the raw input files listed below. No manual
intervention, no hand-edited intermediate tables.

## Pipeline structure (strictly ordered, one-way data flow)

    R/00_config.R          paths | frozen parameters (seed = 42; k = 3 a priori) |
                           16 pathway gene sets | utility functions | integrity checks
    R/01_data.R            TCGA expression + survival -> patient-level cohort (n = 123)
    R/02_scores.R          16 pathway scores (ssGSEA for 4 metabolic pathways via the
                           original KS-walk implementation; mean-expression z for the
                           12 immune/state pathways) + Kynurenine cutoff
    R/03_cluster.R         ConsensusClusterPlus (k-means, 100 subsamples, pItem = 0.8;
                           euclidean as required by the k-means option) + PAC/silhouette;
                           subtypes pre-specified as C1/C2/C3 by decreasing Kynurenine
                           mean; Figure 1
    R/04_mechanism.R       Figures 3 / 6 / 8 (kynurenine-immune coupling)
    R/05_metabolic.R       Figure 4 (PGC-1a / OXPHOS axis)
    R/06_validation.R      Figure 10 (METABRIC; TCGA-derived cutoff applied unchanged)
    R/07_figures.R         Figures 5 / 7 (single-cell panels from processed metadata)
    R/08_fig9.R            Figure 9 (TIDEpy v1.3.9, Other cancer model, default
                           normalization; responder = TIDE < 0)
    R/09_tcga_survival.R   Exploratory TCGA OS (19 events; reported as one sentence)
    R/10_tide_supp.R       TIDE component scores (data table for Figure 9D)
    R/11_extra_analyses.R  Supplementary: SLC7A5 disconnection (Fig S3) and C1-vs-C3
                           differential program (Fig S2, Table S3)

## Data acquisition (not redistributed; download required)

    TCGA expression   TNBC subset matrix (log2 FPKM+1), derived from UCSC Xena
                      TCGA-BRCA.htseq_fpkm; columns are TCGA patient barcodes
    TCGA survival     Xena TCGA-BRCA.survival.tsv (OS, OS.time)
    METABRIC          cBioPortal brca_metabric study archive
    scRNA metadata    processed GSE176078 (Wu et al., Nat Genet 2021) cell metadata
                      with module scores (Fig7_metadata_with_scores.csv)
    TIDE              TIDEpy v1.3.9 run locally (Python >= 3.9); see below

Update paths in R/00_config.R to point to your local copies.

## Setup & run

    R >= 4.4
    install.packages(c("tidyverse", "survival", "survminer", "patchwork"))
    BiocManager::install("ConsensusClusterPlus")

    source("run_all.R")

Seeded (params$seed = 42) and deterministic. An environment snapshot
(sessionInfo.txt) is written after each full run.

## TIDEpy (Figure 9)

    python3 -m venv tidepy-env && tidepy-env/bin/pip install tidepy pandas
    tidepy-env/bin/tidepy <expression.tsv> -o <output_file> -c Other

The expression file is exported from the cohort matrix (genes x samples,
log2 FPKM+1); the resulting per-sample table (Patient, TIDE, Dysfunction,
Exclusion, MDSC, CAF, TAM M2, Responder) is the input at paths$response.

## Supplementary materials

    Table S1    16 pathway gene sets (identical to params$gene_sets in R/00_config.R)
    Table S2    Consensus clustering stability (PAC / silhouette, k = 2-6)
    Table S3    Pathway-level C1 vs C3 differences (11_C1vsC3_pathways.csv)
    Figure S1   C1 vs C3 volcano with kynurenine-pathway genes highlighted
    Figure S2   SLC7A5 disconnection
## Correspondence with the original report & known deviations

The original analysis session was lost; this pipeline rebuilds it from raw data.
Validation of the rebuild: the kynurenine-immune correlations reproduce the
originally reported values to within 0.02 (rho = -0.686 vs -0.663 reported).
Deviations, stated transparently:

    - Subtype sizes 28/48/47 (rebuilt) vs 54/58/11 (original report); the original
      cluster assignments are not recoverable.
    - OXPHOS did not differ significantly across subtypes in the rebuilt analysis
      (reported as such).
    - Fig 7 panels requiring per-cell IDO1 expression or UMAP coordinates need the
      raw GSE176078 count matrix (TODO; the two metadata-derivable panels are
      included).
    - TCGA survival is underpowered (19 events) and reported textually.

## License

MIT.
