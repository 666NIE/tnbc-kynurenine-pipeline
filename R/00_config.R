# 00_config.R — paths, parameters, utilities
paths <- list(
  tcga_expr  = "~/Downloads/01_data/TCGA_TNBC_expr.rds",
  tcga_surv  = "~/Downloads/TCGA-BRCA.survival.tsv",
  metabric   = "~/Downloads/brca_metabric",
  scrna_meta = "~/Downloads/fig7_output/Fig7_metadata_with_scores.csv"
)
dirs <- list(results = "results", figures = "figures")
invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

params <- list(
  seed      = 42,
  maxK      = 6,
  k_final   = 3,
  reps      = 100,
  pItem     = 0.8,
  down_genes  = c("KYNU","KMO","CCBL1","CCBL2","AFMID","AHR"),
  cyto_genes  = c("GZMA","PRF1","IFNG"),
  cd8_genes   = c("CD8A","CD8B"),
  etc_genes   = c("MT-CO1","MT-CYB","COX5B","COX4I1","COX6C","NDUFB8","SDHB","UQCRC2",
                  "ATP5F1A","ATP5F1B","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND5",
                  "MT-CO2","MT-CO3","MT-ATP6","MT-ATP8")
)

params$gene_sets <- list(
  Kynurenine = c("IDO1","IDO2","TDO2","AFMID","KYNU","KMO","KYAT1","AADAT","CCBL2","GPT2","HAAO","QPRT","NADSYN1","IL4I1"),
  Glycolysis = c("HK2","PFKP","ALDOA","GAPDH","PGK1","PGAM1","ENO1","PKM","LDHA","SLC2A1","SLC2A3","P4HA1","PFKFB3","PGAM2","TPI1","GPI","FBP1","G6PD"),
  Glutamine  = c("SLC1A5","SLC38A1","SLC38A2","SLC7A5","SLC3A2","SLC7A11","GLS","GLS2","GLUD1","GOT1","GOT2","GPT2","ASNS","MYC"),
  Adenosine  = c("ENTPD1","NT5E","ADORA2A","ADORA2B","ADORA1","CD38","ENPP1","AK3","AK4"),
  Tryptophan_metabolism = c("TAT","IDO1","IDO2","TDO2","KYNU","KMO","CCBL1","CCBL2","AFMID","AHR","IL4I1","DDC"),
  Hypoxia    = c("HIF1A","EPAS1","VEGFA","SLC2A1","PGK1","CA9","BNIP3","EGLN1","LDHA","ENO1"),
  EMT        = c("VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","FN1","MMP2","MMP9"),
  Stemness   = c("ALDH1A1","CD44","PROM1","SOX2","NANOG","MYC","ABCG2","EPCAM"),
  T_cell     = c("CD3D","CD3E","CD3G","CD2","LCK","ZAP70","TRAC"),
  NK_cell    = c("NKG7","GNLY","PRF1","GZMB","GZMA","KLRD1","NCR1"),
  Treg       = c("FOXP3","IL2RA","CTLA4","TNFRSF18","IKZF2","LAG3"),
  Macrophage = c("CD68","CD163","MSR1","MRC1","CSF1R","ITGAM","C1QA","C1QB"),
  Cytotoxicity = c("GZMA","GZMB","GZMK","PRF1","GNLY","IFNG","NKG7"),
  Exhaustion = c("PDCD1","HAVCR2","LAG3","TIGIT","TOX","LAYN","CXCL13"),
  Checkpoint = c("CD274","PDCD1","PDCD1LG2","CTLA4","HAVCR2","LAG3","TIGIT","ENTPD1"),
  IFNG_response = c("IRF1","STAT1","JAK1","GBP1","CXCL9","CXCL10","IDO1")
)

patient_id <- function(x) substr(x, 1, 12)
fmt_p <- function(p) vapply(p, function(pi) {
  if (is.null(pi) || length(pi) == 0 || is.na(pi)) return("p = NA")
  if (pi < 2.2e-16) return("p < 2.2e-16")
  sprintf("p = %.3g", pi)
}, character(1))
kw_p <- function(formula, dat) kruskal.test(formula, dat)$p.value
calc_z <- function(genes, mat) {
  genes <- intersect(genes, rownames(mat))
  as.numeric(scale(colMeans(mat[genes, , drop = FALSE], na.rm = TRUE)))
}
theme_paper <- function() theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 13, face = "bold"))

simple_ssgsea <- function(expr_matrix, gene_list) {
  gene_list <- intersect(gene_list, rownames(expr_matrix))
  if (length(gene_list) < 5) {
    warning("Only ", length(gene_list), " genes found")
    return(setNames(rep(NA, ncol(expr_matrix)), colnames(expr_matrix)))
  }
  vapply(seq_len(ncol(expr_matrix)), function(s) {
    r <- rank(expr_matrix[, s], ties.method = "average")
    n <- length(r)
    ri <- r[names(r) %in% gene_list]; ro <- r[!(names(r) %in% gene_list)]
    x <- sort(c(ri, ro))
    # NOTE: indexing mirrors the original implementation; recycling quirk kept for exact reproduction
    suppressWarnings(sum(diff(c(0, x)) * (ecdf(ri)(x) - ecdf(ro)(x))[-length(x)]))
  }, numeric(1))
}

need <- c("seed","maxK","k_final","reps","pItem","gene_sets","etc_genes")
miss <- setdiff(need, names(params))
if (length(miss) > 0) stop("config incomplete: ", paste(miss, collapse = ", "))
message("config loaded.")
paths$response <- "~/Downloads/TNBC_output/Table_S2_immunotherapy_scores.csv"
paths$response <- "完整路径"
paths$response <- "~/Downloads/TIDE_output.csv"
