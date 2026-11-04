# 35_spatial_slice_level.R — 终版 v3：与原文 21_spatial_kp.R 完全同口径
# TNBC 4 切片 + Classification 过滤（去 NA/空/Artefact）+ KP 列为 KP
source("R/00_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(lmerTest) })

spots <- readRDS("results/21_spots.rds") %>%
  filter(sample %in% c("1142243F", "1160920F", "CID4465", "CID44971"),
         !is.na(Classification), Classification != "", Classification != "Artefact") %>%
  transmute(.slice = sample, .comp = Classification, .kp = KP, .t = Tcell)

cat("口径核对（应与论文一致: n = 11,877）: n =", nrow(spots), "\n")
print(table(spots$.slice))

# ① 原 spots 级 KW（复现论文 p）
kw <- kruskal.test(.kp ~ .comp, data = spots)
cat(sprintf("spots 级 KW: chi2 = %.1f, p = %.3g\n", unname(kw$statistic), kw$p.value))

# ② 混合模型：slice 随机截距
m <- lmer(.kp ~ .comp + (1 | .slice), data = spots)
print(anova(m))

# ③ 切片内置换检验
slice_comp <- spots %>% group_by(.slice, .comp) %>%
  summarise(kp_mean = mean(.kp), .groups = "drop")
obs_stat <- slice_comp %>% group_by(.slice) %>% summarise(mm = mean(kp_mean), .groups = "drop")
obs <- sum((obs_stat$mm - mean(obs_stat$mm))^2)
set.seed(42)
perm <- replicate(2000, {
  sp <- spots
  sp$.pc <- unlist(tapply(seq_len(nrow(sp)), sp$.slice,
                          function(ii) sample(sp$.comp[ii])))
  pc <- sp %>% group_by(.slice, comp = .pc) %>%
    summarise(kp_mean = mean(.kp), .groups = "drop")
  st <- pc %>% group_by(.slice) %>% summarise(mm = mean(kp_mean), .groups = "drop")
  sum((st$mm - mean(st$mm))^2)
})
cat(sprintf("切片内置换检验 p = %.4f\n", (sum(perm >= obs) + 1) / (length(perm) + 1)))

# ④ 每切片 KP–T 细胞 Spearman + 合并值（复现论文 rho = 0.054）
spots %>% group_by(.slice) %>%
  summarise(rho = cor(.kp, .t, method = "spearman"),
            n = n(), .groups = "drop") %>% print()
ra <- cor.test(spots$.kp, spots$.t, method = "spearman")
cat(sprintf("合并 spots: rho = %.3f, p = %.3g, n = %d\n",
            unname(ra$estimate), ra$p.value, nrow(spots)))

# ⑤ 切片 × 区室均值表（写正文/回应信用）
spots %>% group_by(.slice, .comp) %>%
  summarise(kp_mean = mean(.kp), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = .comp, values_from = c(kp_mean, n)) %>% print()
