############################################################
# Rarefaction and procrustes
############################################################

library(vegan)

# ---------- Paths ----------
project_dir <- normalizePath(".", mustWork = TRUE)
dada_dir <- file.path(project_dir, "DADA2_outputs")
outdir <- file.path(project_dir, "rarefaction_procrustes")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Load data ----------
seqtab263 <- readRDS(file.path(dada_dir, "seqtab_read_red263.rds"))
seqtab331 <- readRDS(file.path(dada_dir, "seqtab_read_red331.rds"))

# ---------- Rarefaction to minimum depth ----------
min_depth_263 <- min(rowSums(seqtab263))
min_depth_331 <- min(rowSums(seqtab331))
min_depth <- min(c(min_depth_263, min_depth_331))

set.seed(123)
seqtab263_rare <- rrarefy(seqtab263, sample = min_depth)
set.seed(123)
seqtab331_rare <- rrarefy(seqtab331, sample = min_depth)

# ---------- Relative abundance matrices ----------
rel263      <- decostand(seqtab263, method = "total")
rel263_rare <- decostand(seqtab263_rare, method = "total")

rel331      <- decostand(seqtab331, method = "total")
rel331_rare <- decostand(seqtab331_rare, method = "total")

# ---------- Bray-Curtis distances between original and rarefied datasets ----------
bray_263_pair <- sapply(seq_len(nrow(rel263)), function(i) {
  vegdist(rbind(rel263[i, ], rel263_rare[i, ]), method = "bray")[1]
})

bray_331_pair <- sapply(seq_len(nrow(rel331)), function(i) {
  vegdist(rbind(rel331[i, ], rel331_rare[i, ]), method = "bray")[1]
})

# ---------- NMDS ordinations ----------
set.seed(123)
nmds263 <- metaMDS(rel263, distance = "bray", trace = FALSE)

set.seed(123)
nmds263_rare <- metaMDS(rel263_rare, distance = "bray", trace = FALSE)

set.seed(123)
nmds331 <- metaMDS(rel331, distance = "bray", trace = FALSE)

set.seed(123)
nmds331_rare <- metaMDS(rel331_rare, distance = "bray", trace = FALSE)

# ---------- Procrustes comparison ----------
proc_263 <- procrustes(nmds263, nmds263_rare, symmetric = TRUE)
proc_331 <- procrustes(nmds331, nmds331_rare, symmetric = TRUE)

prot_263 <- protest(nmds263, nmds263_rare, permutations = 999)
prot_331 <- protest(nmds331, nmds331_rare, permutations = 999)

# ---------- Summary table ----------
# Mantel test was removed because it was a redundant exploratory analysis;
# the manuscript reports the Procrustes/protest comparison for the rarefaction effect.
rare_summary <- data.frame(
  amplicon = c("263", "331"),
  min_depth = c(min_depth, min_depth),
  median_bray_pair = c(median(bray_263_pair), median(bray_331_pair)),
  max_bray_pair = c(max(bray_263_pair), max(bray_331_pair)),
  procrustes_corr = c(1 - proc_263$ss, 1 - proc_331$ss),
  protest_corr = c(prot_263$t0, prot_331$t0),
  protest_p = c(prot_263$signif, prot_331$signif)
)

print(rare_summary)

write.csv(
  rare_summary,
  file.path(outdir, "rarefaction_comparison_summary.csv"),
  row.names = FALSE
)

# ---------- Procrustes plots ----------
png(file.path(outdir, "FigS2_Procrustes_263_rarefied_vs_original.png"), width = 1600, height = 1600, res = 300)
plot(proc_263, main = "Procrustes: 263-bp original vs rarefied")
dev.off()

png(file.path(outdir, "FigS2_Procrustes_331_rarefied_vs_original.png"), width = 1600, height = 1600, res = 300)
plot(proc_331, main = "Procrustes: 331-bp original vs rarefied")
dev.off()
