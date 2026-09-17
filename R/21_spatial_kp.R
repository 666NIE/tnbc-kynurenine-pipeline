# 21_spatial_kp.R -- 空间转录组定位：KP 模块在 TNBC 空间切片中的分布（图 S7）
# 数据：Zenodo 10.5281/zenodo.4739739（Wu et al. 2021；6 张 Visium，含 4 张 TNBC）
source("R/00_config.R")
suppressPackageStartupMessages({library(tidyverse); library(Matrix); library(gridExtra)})

kp <- c("IDO1","IDO2","TDO2","AFMID","KMO","KYNU","KYAT1","AADAT","CCBL2","GPT2","HAAO","QPRT","NADSYN1","IL4I1")
tcell_markers <- c("CD8A","CD8B","CD3D","CD3E","CD4","IL7R")

sp_dir <- path.expand("~/GSE176078_spatial")
dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)
for (f in c("filtered_count_matrices.tar.gz", "spatial.tar.gz", "metadata.tar.gz")) {
  src <- path.expand(file.path("~/Downloads", f))
  tag <- file.path(sp_dir, paste0(".done_", gsub("[.]", "_", f)))
  if (!file.exists(src)) stop("缺文件: ", src)
  if (!file.exists(tag)) { message("解压: ", f); untar(src, exdir = sp_dir); file.create(tag) }
}

mtx_files <- list.files(sp_dir, pattern = "matrix[.]mtx", recursive = TRUE, full.names = TRUE)
pos_files <- list.files(sp_dir, pattern = "tissue_positions", recursive = TRUE, full.names = TRUE)
meta_files <- list.files(sp_dir, pattern = "_metadata[.]csv$", recursive = TRUE, full.names = TRUE)
cat("矩阵:", length(mtx_files), "| 坐标:", length(pos_files), "| metadata:", length(meta_files), "\n")
stopifnot(length(mtx_files) > 0, length(pos_files) > 0, length(meta_files) > 0)

spot_all <- list()
for (i in seq_along(mtx_files)) {
  d <- dirname(mtx_files[i])
  sample_id <- regmatches(mtx_files[i], regexpr("1142243F|1160920F|CID4290|CID4465|CID44971|CID4535", mtx_files[i]))
  if (sample_id %in% c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")) sample_id <- basename(d)
  bc <- readLines(gzfile(file.path(d, list.files(d, pattern = "barcodes[.]tsv")[1])))
  ft <- read.delim(gzfile(file.path(d, list.files(d, pattern = "features[.]tsv")[1])), header = FALSE, stringsAsFactors = FALSE)
  sym <- if (ncol(ft) >= 2) ft[[2]] else ft[[1]]
  m <- readMM(gzfile(mtx_files[i]))
  if (nrow(m) == length(sym)) { colnames(m) <- bc; rownames(m) <- sym } else { m <- t(m); colnames(m) <- bc; rownames(m) <- sym }
  pf <- pos_files[grepl(sample_id, pos_files)]; pf <- if (length(pf) > 0) pf[1] else NA
  if (is.na(pf)) { message("跳过（无坐标）: ", sample_id); next }
  pos <- read.csv(pf, header = FALSE)
  coord <- data.frame(barcode = pos[[1]], x = as.numeric(pos[[5]]), y = as.numeric(pos[[6]]))
  mf <- meta_files[grepl(sample_id, meta_files)][1]
  meta <- read.csv(mf, check.names = FALSE)
  if (!"barcode" %in% colnames(meta)) meta$barcode <- meta[[1]]
  cells <- intersect(colnames(m), coord$barcode)
  message(sprintf("[%s] %d spots x %d genes | 坐标匹配 %d | metadata %d 列",
                  sample_id, ncol(m), nrow(m), length(cells), ncol(meta)))
  if (length(cells) < 500) next
  m <- m[, cells]; lib <- colSums(m)
  sub <- kp[kp %in% rownames(m)]
  cnt <- t(m[sub, , drop = FALSE])
  ids <- rownames(cnt)  # Diagonal %*% 会吞行名，先保存
  cnt <- Diagonal(x = 10000 / lib[rownames(cnt)]) %*% cnt
  cnt <- as(cnt, "dgCMatrix"); cnt@x <- log1p(cnt@x)
  rownames(cnt) <- ids
  spot <- as.data.frame(as.matrix(cnt))
  z <- function(x) as.numeric(scale(x))
  spot$KP <- rowMeans(sapply(sub, function(g) z(spot[[g]])))
  gz <- function(genes) { gg <- intersect(genes, rownames(m))
    rowMeans(sapply(gg, function(g) z(as.numeric(log1p(m[g, cells] * 10000 / lib[cells]))))) }
  spot$Tcell <- gz(tcell_markers)
  spot$barcode <- rownames(spot); spot$sample <- sample_id
  spot <- left_join(spot, coord, by = "barcode") %>% left_join(meta, by = "barcode")
  spot_all[[sample_id]] <- spot
}
spots <- bind_rows(spot_all)
cat("\nsample 分布:\n"); print(table(spots$sample))
saveRDS(spots, file.path(dirs$results, "21_spots.rds"))

# ================= 第二阶段：TNBC 空间分析（图 S7） =================
tnbc <- c("1142243F", "1160920F", "CID4465", "CID44971")
sp <- spots %>% filter(sample %in% tnbc, !is.na(Classification),
                       Classification != "", Classification != "Artefact")
cat("\nTNBC spots by sample x class:\n"); print(table(sp$sample, sp$Classification))

grad <- function(hi) scale_color_gradient2(low = "grey85", mid = "white", high = hi, midpoint = 0)
smap <- function(smp, score, hi, ttl) ggplot(sp %>% filter(sample == smp), aes(x, y, color = .data[[score]])) +
  geom_point(size = 1.0) + grad(hi) + scale_y_reverse() + coord_fixed() +
  theme_void(base_size = 11) + labs(title = ttl, color = NULL)

ps <- c("1142243F", "CID44971")
pA1 <- smap(ps[1], "KP",   "#B22222", paste0(ps[1], " — KP score"))
pA2 <- smap(ps[1], "Tcell","#1F4E8C", paste0(ps[1], " — T-cell score"))
pA3 <- smap(ps[2], "KP",   "#B22222", paste0(ps[2], " — KP score"))
pA4 <- smap(ps[2], "Tcell","#1F4E8C", paste0(ps[2], " — T-cell score"))

pB <- ggplot(sp, aes(reorder(Classification, KP, median), KP, fill = sample)) +
  geom_boxplot(outlier.shape = NA, alpha = .7) + coord_flip() + theme_bw(base_size = 11) +
  labs(x = NULL, y = "KP spot score", title = "KP score across pathological compartments (TNBC sections)")

cor_res <- sp %>% group_by(sample) %>%
  summarise(rho = round(unname(suppressWarnings(cor.test(KP, Tcell, method = "spearman")$estimate)), 3),
            p = format.pval(suppressWarnings(cor.test(KP, Tcell, method = "spearman")$p.value)),
            n = n(), .groups = "drop")
print(cor_res)
ra <- suppressWarnings(cor.test(sp$KP, sp$Tcell, method = "spearman"))
cat(sprintf("Pooled TNBC spots: rho = %.3f, p = %s, n = %d\n",
            unname(ra$estimate), format.pval(ra$p.value), nrow(sp)))
cat(sprintf("KP across pathological classes: KW p = %s\n",
            format.pval(kruskal.test(KP ~ Classification, data = sp)$p.value)))
pC <- ggplot(sp, aes(KP, Tcell, color = sample)) + geom_point(size = .7, alpha = .5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") + theme_bw(base_size = 11) +
  labs(title = "Spot-level KP vs T-cell score (TNBC sections)", x = "KP score", y = "T-cell score")

ggsave(file.path(dirs$figures, "FigureS7_spatial.pdf"),
       arrangeGrob(pA1, pA2, pA3, pA4, pB, pC, ncol = 2,
                   layout_matrix = rbind(c(1,2), c(3,4), c(5,5), c(6,6))),
       width = 13, height = 17)
write.csv(sp %>% group_by(sample, Classification) %>%
            summarise(KP = mean(KP), Tcell = mean(Tcell), n = n(), .groups = "drop"),
          file.path(dirs$results, "21_spatial_by_class.csv"), row.names = FALSE)
write.csv(cor_res, file.path(dirs$results, "21_spatial_correlations.csv"), row.names = FALSE)
message("21 done -> figures/FigureS7_spatial.pdf")
