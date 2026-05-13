#############################################################
# Correlation between 263-bp marker and 331-bp marker
#############################################################

library(vegan)
library(ggplot2)
library(dplyr)

# ---------- Paths and load data ----------
project_dir <- normalizePath(".", mustWork = TRUE)
dada_dir <- file.path(project_dir, "DADA2_outputs")
outdir <- file.path(project_dir, "ASV_correlation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

seqtab263 <- readRDS(file.path(dada_dir, "seqtab_read_red263.rds"))
seqtab331 <- readRDS(file.path(dada_dir, "seqtab_read_red331.rds"))

# ---------- Convert 331 marker to core263 ----------
seqs331 <- colnames(seqtab331)

core263 <- substr(
  seqs331,
  nchar(seqs331) - 263 + 1,
  nchar(seqs331)
)

collapse_counts <- table(core263)
unique_core <- unique(core263)

seqtab_core263 <- sapply(unique_core, function(core) {
  cols <- which(core263 == core)
  rowSums(seqtab331[, cols, drop = FALSE])
}, simplify = "matrix")

seqtab_core263 <- as.matrix(seqtab_core263)
colnames(seqtab_core263) <- unique_core
rownames(seqtab_core263) <- rownames(seqtab331)

# ---------- Shared ASV/core and shared samples ----------
asv_shared <- intersect(colnames(seqtab263), colnames(seqtab_core263))
samples_shared <- intersect(rownames(seqtab263), rownames(seqtab_core263))

# ---------- Relative abundance ----------
rel263 <- seqtab263 / rowSums(seqtab263)
rel331 <- seqtab_core263 / rowSums(seqtab_core263)

abd_shared_263 <- rel263[samples_shared, asv_shared, drop = FALSE]
abd_shared_331 <- rel331[samples_shared, asv_shared, drop = FALSE]

abd_shared_331 <- abd_shared_331[
  match(rownames(abd_shared_263), rownames(abd_shared_331)),
  match(colnames(abd_shared_263), colnames(abd_shared_331)),
  drop = FALSE
]

stopifnot(all(rownames(abd_shared_263) == rownames(abd_shared_331)))
stopifnot(all(colnames(abd_shared_263) == colnames(abd_shared_331)))

# ---------- Per-sample correlation ----------
sample_corr <- sapply(seq_len(nrow(abd_shared_263)), function(i) {
  cor(abd_shared_263[i, ], abd_shared_331[i, ], method = "spearman")
})

sample_corr_df <- data.frame(
  sample = rownames(abd_shared_263),
  spearman_rho = as.numeric(sample_corr),
  stringsAsFactors = FALSE
)

print(sample_corr_df)
print(median(sample_corr, na.rm = TRUE))

# ---------- Mantel test for overall ASV-level community consistency ----------
dist263 <- vegdist(abd_shared_263, method = "bray")
dist331 <- vegdist(abd_shared_331, method = "bray")

mantel_res <- mantel(dist263, dist331, method = "spearman")
print(mantel_res)

# ---------- Scatter plot data ----------
abd_scatter <- data.frame(
  ASV263 = rep(asv_shared, each = length(samples_shared)),
  sample = rep(samples_shared, times = length(asv_shared)),
  abd_263 = as.numeric(abd_shared_263),
  abd_331 = as.numeric(abd_shared_331),
  collapse_n = rep(as.numeric(collapse_counts[asv_shared]), each = length(samples_shared)),
  stringsAsFactors = FALSE
)

abd_scatter$collapse_n[is.na(abd_scatter$collapse_n)] <- 1

abd_scatter <- abd_scatter %>%
  mutate(
    abd_263 = abd_263 + 1e-6,
    abd_331 = abd_331 + 1e-6
  )

p_abd <- ggplot(abd_scatter, aes(x = abd_263, y = abd_331)) +
  geom_point(
    alpha = 0.35,
    shape = 21,
    stroke = 0.3,
    size = 1.5,
    fill = "grey40",
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_log10(
    breaks = c(1e-6, 1e-4, 1e-2, 1),
    labels = c("0", "1e-4", "1e-2", "1")
  ) +
  scale_y_log10(
    breaks = c(1e-6, 1e-4, 1e-2, 1),
    labels = c("0", "1e-4", "1e-2", "1")
  ) +
  labs(
    x = "marker 263-bp",
    y = "marker 331-bp"
  ) +
  theme_classic(base_size = 14)

print(p_abd)

# ---------- Outputs ----------
write.csv(sample_corr_df, file.path(outdir, "Fig1_sample_spearman_correlations.csv"), row.names = FALSE)
write.csv(abd_scatter, file.path(outdir, "Fig1_shared_ASV_abundance_scatter_data.csv"), row.names = FALSE)
saveRDS(mantel_res, file.path(outdir, "Fig1_mantel_shared_ASV_bray.rds"))
saveRDS(seqtab_core263,file.path(outdir, "seqtab_core263_from_331.rds"))

ggsave(file.path(outdir, "Fig1_correl_persample_abd_marker.png"), p_abd, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(outdir, "Fig1_correl_persample_abd_marker.pdf"), p_abd, width = 5.5, height = 5)
