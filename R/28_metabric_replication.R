# 28_metabric_replication.R -- METABRIC replication of matched null + multivariable arbitration
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
x <- readRDS(file.path(dirs$results, "06_metabric.rds"))
mdf <- x$metabric
mexpr <- readRDS(file.path(dirs$results, "06_metabric_expr.rds"))
asg <- read.csv(file.path(dirs$results, "13_metabric_assignments.csv"), check.names = FALSE)
mcl <- asg$cluster[match(mdf$PATIENT_ID, asg$patient)]
kp <- intersect(params$gene_sets$Kynurenine, rownames(mexpr))
gmean <- rowMeans(mexpr); qs <- quantile(gmean, seq(0, 1, 0.1))
kp_bin <- cut(rowMeans(mexpr[kp, , drop = FALSE]), qs, include.lowest = TRUE)
bin_pool <- split(rownames(mexpr), cut(gmean, qs, include.lowest = TRUE))
set.seed(99); z_c1 <- numeric(1000)
for (i in 1:1000) {
  gs <- vapply(seq_along(kp), function(j) sample(bin_pool[[as.character(kp_bin[j])]], 1), character(1))
  rss <- as.numeric(scale(as.numeric(simple_ssgsea(mexpr, gs))))
  z_c1[i] <- mean(rss[mcl == "M1"])
  if (i %% 200 == 0) cat(i, "\n")
}
real <- as.numeric(scale(as.numeric(simple_ssgsea(mexpr, kp))))
cat(sprintf("METABRIC M1 z = %.3f | null %.3f (SD %.3f) | empirical p = %.4f\n",
            mean(real[mcl == "M1"]), mean(z_c1), sd(z_c1),
            mean(abs(z_c1) >= abs(mean(real[mcl == "M1"])))))
cs <- read.delim(path.expand("~/Downloads/CIBERSORTx_METABRIC_Results.txt"), check.names = FALSE)
lib <- colSums(mexpr); bg <- apply(mexpr, 2, median)
df <- mdf %>% mutate(cluster = mcl, CD8 = cs[["T cells CD8"]],
                     lib = as.numeric(scale(log10(lib))), bg = as.numeric(scale(bg)))
m <- lm(Kyn_score ~ cluster + CD8 + lib + bg, data = df)
print(anova(m))
cat(sprintf(">>> METABRIC subtype omnibus p = %.4g\n", anova(m)$"Pr(>F)"[1]))
message("28 done")
