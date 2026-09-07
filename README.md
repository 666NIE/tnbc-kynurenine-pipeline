# Kynurenine-Centered Metabolic-Immune Trap in TNBC — Reproducible Pipeline

End-to-end R pipeline for the TCGA-BRCA discovery, scRNA validation, and METABRIC external cohort analyses.

## Structure
    R/00_config.R     paths, frozen parameters, 16 pathway gene sets, simple_ssgsea
    R/01_data.R       data loading, ID alignment, cohort definition (n = 123)
    R/02_scores.R     16-pathway ssGSEA scores + Kynurenine cutoff
    R/03_cluster.R    consensus clustering (k=3) + stability + Figure 1
    R/04_mechanism.R  Figures 3 / 6 / 8
    R/05_metabolic.R  Figure 4 (PGC-1a / OXPHOS)
    R/06_validation.R Figure 10 (METABRIC); Fig 7 (scRNA) and Fig 9 (TIDE) are TODO
    run_all.R         sources all scripts in order

## Data (download required; not redistributed)
    TCGA expr     Xena TCGA-BRCA.htseq_fpkm-derived TNBC matrix (log2 FPKM+1), rds
    TCGA survival Xena TCGA-BRCA.survival.tsv
    METABRIC      cBioPortal brca_metabric
    scRNA         GSE176078 (Wu et al. 2021), processed metadata with module scores
    TIDE          tide.dfci.harvard.edu per-sample output (Fig 9)

## Setup
    R >= 4.4; install.packages(c("tidyverse","survival","survminer","patchwork"))
    BiocManager::install("ConsensusClusterPlus")

## Run
    source("run_all.R")

Seeded (params$seed = 42); deterministic. NOTE: ConsensusClusterPlus kmeans requires euclidean distance; Methods text "Pearson + K-means" implemented as kmeans + euclidean.

## Known limitations
    Immune-pathway gene sets are reconstructions per Methods 2.3; original definitions were lost with an intermediate object. Swap-in point: params$gene_sets in R/00_config.R.
    Fig 7 / Fig 9 require external inputs (scRNA metadata file; TIDE web output) and are not yet automated.

