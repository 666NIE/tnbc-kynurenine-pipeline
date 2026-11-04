# 40_gse25066_validation.R —— GSE25066 SuperSeries（508 例）trap 表型验证 v2
# 终点：pCR（Fisher）+ DRFS（KM/Cox）。口径与 I-SPY2 / GSE21653 一致。
source("R/00_config.R")
suppressPackageStartupMessages({ library(survival); library(survminer) })

# ---- 0. 定位文件 ----
find_file <- function(fname){
  for (d in c(".", "data-raw/manual_downloads", "results", "data")) {
    if (dir.exists(d)) {
      f <- list.files(d, pattern = fname, full.names = TRUE, ignore.case = TRUE)
      if (length(f)) { message("found: ", f[1]); return(f[1]) }
    }
  }
  message("未找到 ", fname, "，请弹窗选择："); file.choose()
}
f_expr <- find_file("GSE25066_series_matrix.txt")
f_ann  <- find_file("GPL96-57554.txt")

# ---- 1. 读 series matrix（read.table 稳健版）----
con <- file(f_expr, open = "r"); lines <- readLines(con); close(con)
lines <- gsub("\r$", "", lines)   # 关键修复：剥掉 Windows 换行符（行尾隐藏的 \r 会摧毁引号解析）
tb <- which(lines == "!series_matrix_table_begin")
te <- which(lines == "!series_matrix_table_end")
stopifnot(length(tb) == 1, length(te) == 1)
tbl <- lines[(tb+1):(te-1)]
con2 <- textConnection(tbl)
m <- read.table(con2, sep = "\t", quote = "\"", comment.char = "", fill = TRUE,
                na.strings = c("", "NA"), stringsAsFactors = FALSE)
close(con2)
probe_ids <- as.character(m[[1]])
mat <- as.matrix(m[, -1]); storage.mode(mat) <- "numeric"
if (probe_ids[1] %in% c("ID_REF", "ID")) {   # 万一首行是表头则丢弃
  mat <- mat[-1, ]; probe_ids <- probe_ids[-1]
}
rownames(mat) <- probe_ids
gsm_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
colnames(mat) <- sub('^"|"$', '', strsplit(gsm_line, "\t")[[1]][-1])
message("表达矩阵: ", nrow(mat), " 探针 x ", ncol(mat), " 样本")

# ---- 2. 临床特征 ----
ch_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
ch_split <- lapply(strsplit(sub('^!Sample_characteristics_ch1[ \t]+"|"[ \t]*$', '', ch_lines), "\t"),
                   function(x) sub('^"|"$', '', x))
parse_kv <- function(x){ s <- strsplit(x, ": ", fixed = TRUE)[[1]]
  if (length(s) >= 2) list(name = s[1], value = paste(s[-1], collapse = ": ")) else list(name = NA, value = x) }
ch_names <- sapply(ch_split, function(v) parse_kv(v[1])$name)
ch_df <- as.data.frame(setNames(lapply(ch_split, function(v) vapply(v, function(x) parse_kv(x)$value, character(1))),
                                ch_names), stringsAsFactors = FALSE)
rownames(ch_df) <- colnames(mat)
ch_df[] <- lapply(ch_df, function(col) gsub('"', "", as.character(col)))  # 全局清残留引号
message("===== 临床字段（全部）====="); print(ch_names)
pick <- function(patterns){
  hit <- ch_names[Reduce(`|`, lapply(patterns, function(p) grepl(p, ch_names, ignore.case = TRUE)))]
  if (length(hit)) hit[1] else NA
}
f_er  <- pick(c("^er", "estrogen")); f_pr <- pick(c("^pr", "progesterone"))
f_her <- pick(c("her2", "erb"))
f_pcr <- pick(c("pcr", "path"))
f_ev  <- pick(c("drfs.*event", "event.*drfs", "^event", "recurrence", "dmfs.*event"))
f_tm  <- pick(c("drfs.*time", "time.*drfs", "^time", "month", "dmfs.*time"))
message("自动匹配: er=", f_er, " pr=", f_pr, " her2=", f_her,
        " pcr=", f_pcr, " event=", f_ev, " time=", f_tm)
for (ff in c(f_er, f_pr, f_her, f_pcr, f_ev)) if (!is.na(ff)) print(table(ch_df[[ff]], useNA = "ifany"))

# ---- 3. 探针 → 基因 ----
ann <- read.delim(f_ann, comment.char = "#", check.names = FALSE)
gcol <- grep("Gene Symbol|GENE_SYMBOL", colnames(ann), ignore.case = TRUE, value = TRUE)[1]
pcol <- grep("^(ID|ID_REF)$", colnames(ann), value = TRUE)[1]
ann <- ann[ann[[pcol]] %in% rownames(mat), c(pcol, gcol)]
colnames(ann) <- c("probe", "symbol"); ann <- ann[ann$symbol != "", ]
ann$symbol <- sapply(strsplit(as.character(ann$symbol), " /// "), `[`, 1)
best <- tapply(seq_len(nrow(ann)), ann$symbol,
               function(ii) ii[which.max(rowMeans(mat[ann$probe[ii], , drop = FALSE]))])
e <- mat[ann$probe[unlist(best)], , drop = FALSE]; rownames(e) <- names(unlist(best))
message("基因矩阵: ", nrow(e), " x ", ncol(e), " | IDO1 在矩阵: ", "IDO1" %in% rownames(e))

# ---- 4. 圈 TNBC ----
auto_neg <- function(x){
  x <- tolower(as.character(x))
  x <- gsub('"$', '', x)                      # 清残留引号
  ifelse(x %in% c("n", "neg", "negative", "0", "-", "-1") | grepl("neg", x), TRUE,
  ifelse(x %in% c("p", "pos", "positive", "1", "+") | grepl("pos", x), FALSE, NA))
}
er_n  <- auto_neg(ch_df[[f_er]]); pr_n <- auto_neg(ch_df[[f_pr]])
her_n <- if (!is.na(f_her)) auto_neg(ch_df[[f_her]]) else rep(TRUE, length(er_n))
message("er 取值: "); print(table(ch_df[[f_er]], useNA = "ifany"))
message("pr 取值: "); print(table(ch_df[[f_pr]], useNA = "ifany"))
if (!is.na(f_her)) { message("her2 取值: "); print(table(ch_df[[f_her]], useNA = "ifany")) }
sel  <- which(er_n == TRUE & pr_n == TRUE & (is.na(her_n) | her_n == TRUE))  # which 自动剔除 NA
mids <- colnames(e)[sel]
message("er 判定 NA 数: ", sum(is.na(er_n)), " | pr NA 数: ", sum(is.na(pr_n)))
message("TNBC 样本数: ", length(mids), " / ", ncol(e))
if (length(mids) == 0) stop("TNBC 过滤后样本数为 0——把 er/pr/her2 的 table() 输出发我")
e <- e[, mids]; ch_df <- ch_df[mids, ]

# ---- 5. 16 通路评分 ----
ss_methods <- c("Kynurenine", "Glycolysis", "Glutamine", "Adenosine")
sc <- sapply(names(params$gene_sets), function(nm){
  v <- if (nm %in% ss_methods) simple_ssgsea(e, params$gene_sets[[nm]])
       else calc_z(params$gene_sets[[nm]], e)
  as.numeric(scale(v))
}); rownames(sc) <- colnames(e)
tcol <- grep("t.?cell|tcell|cd8|lymph|t_cell", colnames(sc), ignore.case = TRUE, value = TRUE)[1]
message("评分矩阵列名: "); print(colnames(sc))
message("T 细胞列名: ", tcol)
if (is.na(tcol)) stop("未匹配到 T 细胞列——请把上面打印的列名发给我")

# ---- 6. trap 表型（同口径：组内 z 双低）----
trap <- sc[, tcol] < 0 & as.numeric(scale(as.numeric(e["IDO1", ]))) < 0
message("trap-like: ", sum(trap), " / ", length(trap), " (", round(100*mean(trap),1), "%)")

# ---- 7. pCR ----
pcr <- tolower(ch_df[[f_pcr]])
pcr_bin <- ifelse(pcr %in% c("pcr", "pos", "positive", "1", "yes", "y"), "pCR",
           ifelse(pcr %in% c("rd", "neg", "negative", "0", "no", "n", "residual"), "RD", NA))
keep_p <- !is.na(pcr_bin)
message("pCR 可用样本: ", sum(keep_p), " | 取值分布: "); print(table(pcr, useNA = "ifany"))
tb <- table(factor(trap[keep_p], levels = c(FALSE, TRUE)),
            factor(pcr_bin[keep_p], levels = c("RD", "pCR")))
print(addmargins(tb))
fp <- fisher.test(tb)
cat(sprintf("pCR Fisher: p = %.4f, OR = %.2f (95%%CI %.2f-%.2f)\n",
            fp$p.value, unname(fp$estimate), fp$conf.int[1], fp$conf.int[2]))

# ---- 8. DRFS ----
message("drfs event 取值: "); print(table(ch_df[[f_ev]], useNA = "ifany"))
ev <- tolower(ch_df[[f_ev]])
ev_bin <- as.integer(ev %in% c("1","yes","y","recurrence","event","true","relapse","distant"))
tm <- suppressWarnings(as.numeric(ch_df[[f_tm]]))
keep <- !is.na(ev_bin) & !is.na(tm) & tm > 0
message("DRFS 可用样本: ", sum(keep))
fit <- survdiff(Surv(tm[keep], ev_bin[keep]) ~ trap[keep])
lp <- 1 - pchisq(fit$chisq, 1)
cox <- coxph(Surv(tm[keep], ev_bin[keep]) ~ trap[keep])
hr <- exp(coef(cox)); ci <- exp(confint(cox))
cat(sprintf("DRFS log-rank p = %.5f | HR = %.2f (95%%CI %.2f-%.2f)\n", lp, hr, ci[1], ci[2]))

# ---- 9. 图 + 输出 ----
dd <- data.frame(tm = tm[keep], ev = ev_bin[keep],
                 trap = factor(trap[keep], levels = c(TRUE, FALSE),
                               labels = c("CD8-low/IDO1-low", "Others")))
km <- survfit(Surv(tm, ev) ~ trap, data = dd)
g <- ggsurvplot(km, data = dd, pval = TRUE, risk.table = TRUE,
                title = "GSE25066 TNBC: trap-like profile vs others",
                legend.title = "", palette = c("#E41A1C", "#377EB8"))
ggsave("results/40_gse25066_KM.png", print(g), width = 7, height = 6)
write.csv(data.frame(sample = colnames(e), trap = trap, pcr = pcr_bin,
                     drfs_event = ev_bin, drfs_time = tm),
          "results/40_gse25066_trap_assignment.csv", row.names = FALSE)
write.csv(data.frame(test = c("trap_n","trap_pct","pCR_fisher","DRFS_logrank","DRFS_HR","DRFS_HR_lo","DRFS_HR_hi"),
                     value = c(sum(trap), round(100*mean(trap),1), fp$p.value, lp, hr, ci[1], ci[2])),
          "results/40_gse25066_summary.csv", row.names = FALSE)
message("done")
