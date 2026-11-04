# 33_verify_metabric_tnbc.R —— 独立版（无需 source config）
# 验证：13_metabric_assignments.csv 的 320 例是否为 ER-/PR-/HER2- 阴性（TNBC 亚队列）

# 0. 先关掉可能开着的 sink
sink()

# 1. 定位 METABRIC 临床文件
mdir <- "data-raw/manual_downloads/brca_metabric"   # 已确认的实际路径
if (!dir.exists(mdir)) mdir <- dirname(file.choose())
ms <- list.files(mdir, pattern = "data_clinical_sample",  full.names = TRUE)[1]
mp <- list.files(mdir, pattern = "data_clinical_patient", full.names = TRUE)[1]
message("clinical sample file: ", ms)
message("clinical patient file: ", mp)

# 2. 合并临床表
s <- read.delim(ms, comment.char = "#", check.names = FALSE)
p <- read.delim(mp, comment.char = "#", check.names = FALSE)
mclin <- merge(s, p, by = "PATIENT_ID", all.x = TRUE)

# 3. 自动找 ER/PR/HER2 状态列
er <- grep("ER_STATUS|ER_IHC|IHC_ER",       colnames(mclin), value = TRUE, ignore.case = TRUE)[1]
pr <- grep("PR_STATUS|PR_IHC|IHC_PR",       colnames(mclin), value = TRUE, ignore.case = TRUE)[1]
he <- grep("HER2_STATUS|HER2_IHC|IHC_HER2", colnames(mclin), value = TRUE, ignore.case = TRUE)[1]
message("ER 列: ", er, " | PR 列: ", pr, " | HER2 列: ", he)

# 4. 三阴性过滤 + 与分型文件对照
mtnbc   <- mclin[mclin[[er]] == "Negative" & mclin[[pr]] == "Negative" & mclin[[he]] == "Negative", ]
asg     <- read.csv("results/13_metabric_assignments.csv", check.names = FALSE)
in_tnbc <- asg$patient %in% mtnbc$PATIENT_ID
message("METABRIC 全队列中三阴性例数: ", nrow(mtnbc))
message("分型文件 ", nrow(asg), " 例中，属于三阴性的: ", sum(in_tnbc))
message("不属于三阴性的: ", sum(!in_tnbc), "  <- 必须为 0")

# 5. 320 例受体状态逐表核对
chk <- mclin[match(asg$patient, mclin$PATIENT_ID), c("PATIENT_ID", er, pr, he)]
print(table(chk[[er]], useNA = "ifany"))
print(table(chk[[pr]], useNA = "ifany"))
print(table(chk[[he]], useNA = "ifany"))
print(head(chk, 5))
