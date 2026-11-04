#!/usr/bin/env Rscript
## ============================================================================
## 39_supplementary_composites.R — S1-S4 自动拼图；S5 由三张 panel 于投稿准备阶段手工拼合
# 注：S5a=分层耦合(28号产物) / S5b=TNBC ranking图(34_swap_noAFMID目录) / S5c=Lehmann组成 —— 期刊版式拼图
## 标准：174mm 图宽；panel 按内容比例分栏；小写粗体 a/b/c 标号（与正文一致）；
##       图内不放标题（图注在稿件 Supplementary Materials 清单中）
## 运行前准备（同 43 号脚本）：
##   S1c 源图命名 figures/FigureS3_METABRIC_composite.pdf
##   tif 转 pdf/png：I-SPY2 -> Figure_ICB_validation.pdf；蛋白 -> *_v2.png
## ============================================================================
suppressPackageStartupMessages({
  library(pdftools); library(grid); library(png); library(patchwork); library(ggplot2); library(tidyverse)
})
dir.create("figures", showWarnings = FALSE)
fig_dir  <- "figures"
proj_dir <- path.expand("~/afmid_project")

files <- list(
  S1a = file.path(fig_dir, "FigureS1_KynGenes_C1vsC3.pdf"),
  S1b = file.path(fig_dir, "FigureS2_SLC7A5.pdf"),
  S1c = file.path(fig_dir, "FigureS3_METABRIC_composite.pdf"),
  S2a = file.path(fig_dir, "FigureS4_KP_scRNA.pdf"),
  S2b = file.path(fig_dir, "FigureS5_IMPRES_C1vsnonC1.pdf"),
  S3a = file.path(fig_dir, "FigureS6_dual_model_v5_printsize.pdf"),
  S3b = file.path(fig_dir, "FigureS7_spatial_v2.pdf"),
  S3c = file.path(fig_dir, "FigureS8_cohort_flow.pdf"),
  S4a = file.path(fig_dir, "Figure_ICB_validation.pdf"),
  S4b = file.path(fig_dir, "Figure_S10_protein_level_evidence_v2.png"),
  S4c = file.path(proj_dir, "DepMap_AFMID_dependency.pdf"),
  S5a = file.path(fig_dir, "FigS14_scrna_stratified.pdf"),
  S5b = file.path(fig_dir, "Additional_file_21_TNBC_ranking.pdf"),
  S5c = file.path(fig_dir, "Additional_file_23_lehmann_crossclass.pdf"),
  S6a = file.path(fig_dir, "42c_pancancer.pdf"),
  S6b = file.path(fig_dir, "42e_LUAD_KM.pdf")
)

## ---------- 工具：渲染 + 裁白边 + 包成 grob ----------
trim_white <- function(img, thr = 0.98) {          # img: h*w*(3|4) 数组
  g <- img[, , 1] * 0.299 + img[, , 2] * 0.587 + img[, , 3] * 0.114
  rows <- which(rowMeans(g > thr) < 1)
  cols <- which(colMeans(g > thr) < 1)
  if (length(rows) < 10 || length(cols) < 10) return(img)
  img[min(rows):max(rows), min(cols):max(cols), , drop = FALSE]
}

to_grob <- function(path, dpi = 220) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") {
    pg  <- pdf_render_page(path, dpi = dpi, numeric = TRUE)
    tmp <- tempfile(fileext = ".png"); png::writePNG(pg / 255, tmp)
    img <- png::readPNG(tmp)
  } else if (ext == "png") {
    img <- png::readPNG(path)
  } else stop("请先转成 pdf/png：", path)
  img <- trim_white(img)
  wrap_elements(grid::rasterGrob(grDevices::as.raster(img), interpolate = TRUE))
}

## ---------- 工具：按期刊版式保存 ----------
## width/height 单位 mm；tag 标号样式与主文一致（小写粗体）
save_composite <- function(p, outfile, w_mm = 174, h_mm) {
  ggsave(file.path(fig_dir, outfile), p,
         width = w_mm, height = h_mm, units = "mm", dpi = 300)
  message("done: ", outfile, "  (", w_mm, " x ", h_mm, " mm)")
}

tag_theme <- theme(plot.tag = element_text(face = "bold", size = 14),
                   plot.margin = margin(2, 2, 2, 2))

## ---------- S1：a 通栏 + b/c 并排（b、c 为横向图，高度较小）----------
S1 <- (to_grob(files$S1a)) /
      (to_grob(files$S1b) | to_grob(files$S1c)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_composite(S1, "Additional file 1.pdf", h_mm = 200)

## ---------- S2：两张横图上下排 ----------
S2 <- to_grob(files$S2a) / to_grob(files$S2b) +
  plot_annotation(tag_levels = "a") & tag_theme
save_composite(S2, "Additional file 2.pdf", h_mm = 165)

## ---------- S3：a 概念框架（方形）左 + b/c 上下排右 ----------
S3 <- to_grob(files$S3a) | (to_grob(files$S3b) / to_grob(files$S3c)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_composite(S3, "Additional file 3.pdf", h_mm = 190)

## ---------- S4：a 横图通栏 + b（方形）/c（横图）并排 ----------
S4 <- to_grob(files$S4a) /
      (to_grob(files$S4b) | to_grob(files$S4c)) +
  plot_annotation(tag_levels = "a") & tag_theme
save_composite(S4, "Additional file 4.pdf", h_mm = 210)

## ---------- S5：a/b 并排（近方形）+ c 通栏 ----------
S5 <- (to_grob(files$S5a) | to_grob(files$S5b)) /
      to_grob(files$S5c) +
  plot_annotation(tag_levels = "a") & tag_theme
# save_composite(S5, "Supplementary_Figure_S5.pdf", h_mm = 210)

## ---------- S6：a（竖长条形图）左 + b（KM 竖图）右 ----------
S6 <- to_grob(files$S6a) | to_grob(files$S6b) +
  plot_annotation(tag_levels = "a") & tag_theme
save_composite(S6, "Additional file 6.pdf", h_mm = 185)

cat("\n全部完成 -> figures/Supplementary_Figure_S1.pdf ... S6.pdf\n")
