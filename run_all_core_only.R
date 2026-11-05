# run_all.R — core-only driver for the TNBC kynurenine-pipeline repository
# Scope:
#   Build the core TCGA/METABRIC score/cluster objects and regenerate final Figure 2/Figure 3.
#   It intentionally does NOT run manual-download, spatial, CPTAC, DepMap, I-SPY2, or legacy modules.
# Run from repository root:  Rscript run_all.R

source("R/00_config.R")
stopifnot(exists("dirs"), exists("paths"))

dir.create(dirs$results, recursive = TRUE, showWarnings = FALSE)
dir.create(dirs$figures, recursive = TRUE, showWarnings = FALSE)

core_inputs <- c(
  tcga_expr = paths$tcga_expr,
  tcga_surv = paths$tcga_surv,
  metabric  = paths$metabric
)
missing <- core_inputs[!vapply(core_inputs, function(p) file.exists(path.expand(p)), logical(1))]
if (length(missing)) {
  stop(
    "Missing core input(s):\n  ", paste(missing, collapse = "\n  "),
    "\nPlace them under data-raw/manual_downloads/ as described in data-raw/manual_downloads/README.md"
  )
}

core_scripts <- c(
  "01_data.R",
  "02_scores.R",
  "03_cluster.R",
  "make_figure2_final_style.R",
  "make_figure3_final_style.R"
)

for (s in core_scripts) {
  message("\n===== run_all core: ", s, " =====")
  source(file.path("R", s))
}

message("\nrun_all core-only finished. Core outputs: results/01_*.rds, results/02_*.rds, results/03_*.rds, figures/Figure2.png, figures/Figure3.png")
message("Optional modules (not run here): 04-46 legacy/manual modules, 17/21/37 spatial/scRNA, 40/41/30/31 ICB, 42/43/44/45 pan-cancer/CPTAC.")
