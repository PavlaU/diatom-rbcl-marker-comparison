#############################################################
# Primer mismatch analysis for 331 ASV
# shared vs. only331
# best primer / best mismatch / mismatch positions

# Required input:
# - group_info.csv
#   ASV group information manually curated based on BLAST inspection; included in GitHub
#  repository folder "input".
############################################################

library(dplyr)
library(ggplot2)
library(Biostrings)
library(tidyr)

# ---------- Paths and load data ----------
project_dir <- normalizePath(".", mustWork = TRUE)
dada_dir <- file.path(project_dir, "DADA2_outputs")
outdir <- file.path(project_dir, "primer_mismatch")
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

# ---------- Detail: which 331 ASV form each unmatched core ----------

core_map_df <- data.frame(
  seq331 = seqs331,
  core263 = core263,
  reads331_original = colSums(seqtab331),
  stringsAsFactors = FALSE)


# ---------- Load group info (Bacil/Eustig/Chryso/...) # ---------- 

group_map <- read.csv(
  file.path(project_dir, "input", "group_info.csv"),
  stringsAsFactors = FALSE
)

# ---------- Merge core_map with group info and remove non-target taxa ----------

core_map_df <- core_map_df %>%
  left_join(group_map, by = c("core263" = "sequence_263")) %>%
  filter(tax_group == "Bacil")


# ---------- primer-binding region in all 331 ASVs ----------

core_map_df$pfor <- substr(core_map_df$seq331, 42, 68)

primers <- c(
  p1 = "AGGTGAAGTAAAAGGTTCWTACTTAAA",
  p2 = "AGGTGAAGTTAAAGGTTCWTAYTTAAA",
  p3 = "AGGTGAAACTAAAGGTTCWTACTTAAA"
)

# ---------- helper: mismatch count + positions ----------

get_degenerate_mismatch_info <- function(region_seq, primer_seq) {
  r <- strsplit(region_seq, "", fixed = TRUE)[[1]]
  p <- strsplit(primer_seq, "", fixed = TRUE)[[1]]
  
  mismatch_idx <- which(!mapply(function(base, code) {
    base %in% strsplit(IUPAC_CODE_MAP[[code]], "", fixed = TRUE)[[1]]
  }, r, p))
  
  list(
    n = length(mismatch_idx),
    pos_region = if (length(mismatch_idx) == 0) "" else paste(mismatch_idx, collapse = "-")
  )
}

# ---------- mismatch counts for all 331 ASVs ----------

mm_mat <- do.call(cbind, lapply(primers, function(pr) {
  vapply(core_map_df$pfor, function(x) get_degenerate_mismatch_info(x, pr)$n, integer(1))
}))
colnames(mm_mat) <- names(primers)

# ---------- mismatch positions for all 331 ASVs ----------

pos_region_list <- lapply(primers, function(pr) {
  vapply(core_map_df$pfor, function(x) get_degenerate_mismatch_info(x, pr)$pos_region, character(1))
})
names(pos_region_list) <- names(primers)

# ---------- best primer ----------

best_idx <- max.col(-mm_mat, ties.method = "first")

core_map_df$best_primer <- colnames(mm_mat)[best_idx]
core_map_df$best_mismatches <- mm_mat[cbind(seq_len(nrow(core_map_df)), best_idx)]
core_map_df$best_mismatch_pos <- cbind(
  pos_region_list[["p1"]],
  pos_region_list[["p2"]],
  pos_region_list[["p3"]]
)[cbind(seq_len(nrow(core_map_df)), best_idx)]

core_map_df$primer_tie <- apply(mm_mat, 1, function(x) sum(x == min(x)) > 1)

# ---------- classify 331 ASV as shared vs only331 ----------
asv_shared <- intersect(colnames(seqtab_core263), colnames(seqtab263))
core_only_331 <- setdiff(colnames(seqtab_core263), colnames(seqtab263))

core_map_df$group <- ifelse(core_map_df$core263 %in% asv_shared, "shared", 
                            ifelse(core_map_df$core263 %in% core_only_331, "only331", NA))

core_map_mismatch <- core_map_df %>%
  filter(!is.na(group)) %>%
  mutate(
    group = factor(group, levels = c("shared", "only331")),
    collapse_n = as.numeric(collapse_counts[core263]),
    collapse_n = ifelse(is.na(collapse_n), 1, collapse_n)
  )


# ---------- Mismatch counts: shared vs only331 ----------

# Mismatch summary
mismatch_summary <- core_map_mismatch %>%
  group_by(group) %>%
  summarise(
    n_ASV = n(),
    median_mm = median(best_mismatches),
    mean_mm = mean(best_mismatches),
    max_mm = max(best_mismatches),
    n_zero = sum(best_mismatches == 0),
    .groups = "drop"
  )

print(mismatch_summary)

# Wilcoxon test
mm_test <- wilcox.test(best_mismatches ~ group, data = core_map_mismatch)
print(mm_test)


# ---------- Location of mismatches ----------

# ASV number in each group
mm_pos_summary_complete <- core_map_mismatch %>%
  filter(best_mismatch_pos != "") %>%
  separate_rows(best_mismatch_pos, sep = "-") %>%
  mutate(pos = as.integer(best_mismatch_pos)) %>%
  count(group, pos, name = "n_ASV") %>%
  complete(
    group,
    pos = 1:27,
    fill = list(n_ASV = 0)
  )

n_group <- core_map_mismatch %>%
  group_by(group) %>%
  summarise(total_ASV = n(), .groups = "drop")

# relative proportion
mm_pos_prop <- mm_pos_summary_complete %>%
  left_join(n_group, by = "group") %>%
  mutate(prop_ASV = n_ASV / total_ASV)

# graph
p_miss <- ggplot(mm_pos_prop, aes(x = pos, y = prop_ASV, color = group)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_continuous(breaks = seq(1, 27, by = 2)) +
  theme_bw(base_size = 14) +
  labs(
    x = "Position in primer-binding region",
    y = "Proportion of ASVs with mismatch",
    color = "Group"
  )

p_miss


# ---------- Outputs ----------

write.csv(core_map_mismatch, file.path(outdir, "primer_mismatch_ASV_table.csv"), row.names = FALSE)
write.csv(mismatch_summary, file.path(outdir, "primer_mismatch_summary.csv"), row.names = FALSE)
saveRDS(mm_test, file.path(outdir, "primer_mismatch_wilcox.rds"))

ggsave(file.path(outdir, "FigS4_primer_mismatch_location.png"), p_miss, width = 9, height = 3, dpi = 300)