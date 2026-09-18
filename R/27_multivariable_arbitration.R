# 27_multivariable_arbitration.R -- Kyn ~ subtype + CD8 + library size + global background
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
coh <- readRDS(file.path(dirs$results, "01_cohort.rds"))
scr <- readRDS(file.path(dirs$results, "02_scores.rds"))
cl  <- readRDS(file.path(dirs$results, "03_cluster.rds"))
scores <- scr$scores; cluster <- cl$cluster
cs <- read.delim(path.expand("~/Downloads/CIBERSORTx_Results.txt"), check.names = FALSE)
lib <- colSums(coh$expr); bg <- apply(coh$expr, 2, median)
df <- scores %>% mutate(cluster = cluster, CD8 = cs[["T cells CD8"]],
                        lib = as.numeric(scale(log10(lib))), bg = as.numeric(scale(bg)))
m <- lm(Kynurenine ~ cluster + CD8 + lib + bg, data = df)
print(anova(m))
cat(sprintf(">>> subtype omnibus p = %.4g\n", anova(m)$"Pr(>F)"[1]))
pcor <- function(x, y, covs) {
  rx <- resid(lm(x ~ ., data = as.data.frame(covs)))
  ry <- resid(lm(y ~ ., data = as.data.frame(covs)))
  suppressWarnings(cor.test(rx, ry, method = "spearman"))
}
r <- pcor(df$Kynurenine, df$Cytotoxicity, data.frame(CD8 = df$CD8, lib = df$lib, bg = df$bg))
cat(sprintf("partial corr Kyn~Cyto | CD8+lib+bg: rho = %.3f, p = %s\n",
            unname(r$estimate), format.pval(r$p.value)))
for (v in c("CD8", "lib", "bg")) {
  cc <- suppressWarnings(cor.test(df$Kynurenine, df[[v]], method = "spearman"))
  cat(sprintf("Kyn vs %s: rho = %.3f\n", v, unname(cc$estimate)))
}
message("27 done")
