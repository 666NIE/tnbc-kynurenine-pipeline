# 10_tide_supp.R — Supplementary: Kynurenine vs TIDE component scores
source("R/00_config.R")
suppressPackageStartupMessages(library(tidyverse))
f9  <- read.csv(file.path(dirs$results, "08_tide.csv"))
rdf <- read.csv(path.expand(paths$response), check.names = FALSE)
colnames(rdf)[1] <- "Patient"
comp <- bind_rows(lapply(c("Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2"), function(cc) {
  if (!cc %in% colnames(rdf)) return(NULL)
  v <- as.numeric(rdf[[cc]])[match(f9$patient, patient_id(as.character(rdf$Patient)))]
  ct <- suppressWarnings(cor.test(f9$Kynurenine, v, method = "spearman"))
  data.frame(component = cc, rho = unname(ct$estimate), p = ct$p.value)
})) %>% mutate(label = sprintf("rho=%.2f
%s", rho, fmt_p(p)))
p <- ggplot(comp, aes(reorder(component, rho), rho, fill = rho > 0)) +
  geom_col(width = .7, color = "black", linewidth = .3) +
  geom_text(aes(label = label), size = 3, hjust = ifelse(comp$rho > 0, -0.05, 1.08)) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5")) +
  coord_flip() + theme_paper() + theme(legend.position = "none") +
  ylim(-0.6, 0.75) +
  labs(title = "Kynurenine associates with exclusion and suppressive myeloid infiltration",
       subtitle = "Spearman correlation with TIDE component scores",
       x = NULL, y = "Spearman rho")
ggsave(file.path(dirs$figures, "FigureS_TIDE_components.pdf"), p, width = 7, height = 4.5)
write.csv(comp, file.path(dirs$results, "10_tide_components.csv"), row.names = FALSE)
message("10 done")
