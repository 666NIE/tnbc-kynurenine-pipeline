# 36_gse21653_decomposition.R —— GSE21653 trap profile 成分分解分析（R 版，与 Python 验证口径一致）
# 产出: Figure_S11 三面板 KM + 36_gse21653_decomposition_summary.csv
# 依赖: survival, survminer, ggplot2, patchwork
source("R/00_config.R")
suppressPackageStartupMessages({ library(survival); library(survminer)
  library(ggplot2); library(patchwork) })

find_file <- function(pattern){
  for (dd in c(".", "data-raw/manual_downloads", "data", "results")) {
    if (dir.exists(dd)) {
      f <- list.files(dd, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
      if (length(f)) { message("found: ", f[1]); return(f[1]) }
    }
  }
  stop("找不到 ", pattern)
}
f_expr <- find_file("GSE21653_series_matrix")
f_ann  <- find_file("GPL570")

# ---- 1. 读矩阵 ----
con <- file(f_expr, open = "r"); lines <- readLines(con); close(con)
lines <- gsub("\r$", "", lines)
tb <- which(lines == "!series_matrix_table_begin"); te <- which(lines == "!series_matrix_table_end")
tbl <- lines[(tb+1):(te-1)]
first <- strsplit(tbl[1], "\t")[[1]]
if (first[1] %in% c("ID_REF","ID") || grepl("^GSM", first[2])) tbl <- tbl[-1]
con2 <- textConnection(tbl)
m <- read.table(con2, sep = "\t", quote = "\"", comment.char = "", fill = TRUE,
                na.strings = c("","NA"), stringsAsFactors = FALSE)
close(con2)
probe_ids <- as.character(m[[1]]); mat <- as.matrix(m[, -1]); storage.mode(mat) <- "numeric"
rownames(mat) <- probe_ids
colnames(mat) <- sub('^"|"$', '', strsplit(grep("^!Sample_geo_accession", lines, value = TRUE), "\t")[[1]][-1])

# ---- 2. 临床特征（逐样本字典法）----
ch_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
ch_list <- lapply(ch_lines, function(l)
  sub('^"|"$', '', strsplit(sub('^!Sample_characteristics_ch1[ \t]+"|"[ \t]*$', '', l), "\t")[[1]]))
n <- ncol(mat); clin <- vector("list", n)
for (j in seq_len(n)) {
  kv <- list()
  for (v in sapply(ch_list, `[`, j)) {
    if (is.na(v) || v == "") next
    s <- strsplit(v, ": ", fixed = TRUE)[[1]]
    if (length(s) >= 2) kv[[s[1]]] <- paste(s[-1], collapse = ": ") else kv[["misc"]] <- v
  }
  clin[[j]] <- kv
}
names(clin) <- colnames(mat)
getv <- function(k) vapply(clin, function(x) if (!is.null(x[[k]])) x[[k]] else NA_character_, character(1))
clean <- function(x) gsub('^[ \t]+|[ \t]+$|\"', '', x)
er  <- clean(getv("er ihc")); pr <- clean(getv("pr ihc")); her <- clean(getv("erbb2"))
negv <- function(x){ x <- tolower(x)
  ifelse(x %in% c("0","neg","negative","n","-","no"), "0",
  ifelse(x %in% c("1","pos","positive","p","+","yes"), "1", x)) }
er <- negv(er); pr <- negv(pr); her <- negv(her)
dtime <- suppressWarnings(as.numeric(clean(getv("dfs time (months)"))))
devt  <- suppressWarnings(as.numeric(clean(getv("dfs evt"))))
message("er 取值: "); print(table(er, useNA = "ifany"))
message("pr 取值: "); print(table(pr, useNA = "ifany"))
message("her 取值: "); print(table(her, useNA = "ifany"))

# ---- 3. 探针 -> 基因 ----
ann <- read.delim(f_ann, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
gcol <- grep("Gene Symbol|GENE_SYMBOL", colnames(ann), ignore.case = TRUE, value = TRUE)[1]
pcol <- grep("^(ID|ID_REF)$", colnames(ann), value = TRUE)[1]
ann <- ann[ann[[pcol]] %in% rownames(mat), c(pcol, gcol)]
colnames(ann) <- c("probe","symbol"); ann <- ann[ann$symbol != "", ]
ann$symbol <- sapply(strsplit(ann$symbol, " /// "), `[`, 1)
em <- rowMeans(mat, na.rm = TRUE)
best <- tapply(seq_len(nrow(ann)), ann$symbol, function(ii) ann$probe[ii][which.max(em[ann$probe[ii]])])
e <- mat[unlist(best), , drop = FALSE]; rownames(e) <- names(best)
message("基因矩阵: ", nrow(e), " x ", ncol(e))

# ---- 4. TNBC + z 标准化（全队列口径，复现 trap=24/85）----
mids <- colnames(e)[er == "0" & pr == "0" & her == "0"]
stopifnot(length(mids) > 0)
ez <- t(scale(t(e)))
tm <- intersect(c("CD8A","CD8B","CD3D","CD3E","CD4","IL7R"), rownames(ez))
tscore <- colMeans(ez[tm, mids, drop = FALSE]); ido1z <- ez["IDO1", mids]
d <- data.frame(gsm = mids,
                t = dtime[match(mids, names(clin))] / 12,
                e = devt[match(mids, names(clin))],
                tscore = tscore, ido1z = ido1z, stringsAsFactors = FALSE)
d <- d[!is.na(d$t) & !is.na(d$e), ]
d$cd8low <- d$tscore < 0; d$ido1low <- d$ido1z < 0; d$trap <- d$cd8low & d$ido1low
message(sprintf("n=%d, events=%d | trap=%d, CD8low=%d, IDO1low=%d",
                nrow(d), sum(d$e), sum(d$trap), sum(d$cd8low), sum(d$ido1low)))
stopifnot(sum(d$trap) == 24)   # 复现校验
write.csv(d, "results/36_gse21653_perpatient.csv", row.names = FALSE)  # ★ 导出逐患者数据

# ---- 5. log-rank + Cox + 增量 LRT ----
lr <- list(); cx <- list()
for (g in c("cd8low","ido1low","trap")) {
  sd <- survdiff(Surv(t, e) ~ d[[g]], data = d)
  lr[[g]] <- 1 - pchisq(sd$chisq, 1)
  cx[[g]] <- coxph(Surv(t, e) ~ d[[g]], data = d)
}
c2 <- coxph(Surv(t, e) ~ cd8low + ido1low, data = d)
inc <- anova(c2, cx$cd8low)
cat(sprintf("CD8low HR=%.2f p=%.2e | IDO1low HR=%.2f p=%.2e | trap HR=%.2f p=%.2e | IDO1增量LRT p=%.3f\n",
            exp(coef(cx$cd8low))[1], lr$cd8low, exp(coef(cx$ido1low))[1], lr$ido1low,
            exp(coef(cx$trap))[1], lr$trap, inc$`Pr(>|Chi|)`[2]))

# ---- 6. 三面板 KM（Figure S11）----
mk <- function(g, title){
  fit <- survfit(Surv(t, e) ~ d[[g]], data = d)
  ggsurvplot(fit, data = d, pval = TRUE, risk.table = FALSE, legend.title = "",
             palette = c("#377EB8","#E41A1C"),
             legend.labs = c("Others", gsub("low","-low",title)),
             title = title, xlab = "Time (years)", ylab = "Disease-free survival")$plot
}
g <- (mk("cd8low","CD8-low alone") | mk("ido1low","IDO1-low alone") | mk("trap","Joint profile")) +
     plot_annotation(title = sprintf("GSE21653 TNBC (n=%d, %d events): prognostic decomposition", nrow(d), sum(d$e)))
ggsave("results/36_gse21653_decomposition.png", g, width = 14, height = 5, dpi = 300)

write.csv(data.frame(test = names(lr), logrank_p = unlist(lr),
                     HR = sapply(cx, function(x) exp(coef(x))[1])),
          "results/36_gse21653_decomposition_summary.csv", row.names = FALSE)
message("done")
