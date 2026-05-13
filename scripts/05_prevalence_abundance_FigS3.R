#############################################################
# Prevalence vs maximum relative abundance plot
#############################################################

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# ---------- 0. Paths and load data ----------
project_dir <- normalizePath(".", mustWork = TRUE)
outdir <- file.path(project_dir, "prevalence_abundance")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


seqtab263 <- readRDS(file.path(project_dir, "DADA2_outputs", "seqtab_read_red263.rds"))
seqtab_core263 <- readRDS(file.path(project_dir, "ASV_correlation", "seqtab_core263_from_331.rds"))

# ---------- 1. Convert to long format ----------


make_long <- function(seqtab_obj, marker_name) {
  as.data.frame(seqtab_obj) %>%
    rownames_to_column("sample") %>%
    pivot_longer(
      cols = -sample,
      names_to = "core263",
      values_to = "reads"
    ) %>%
    filter(reads > 0) %>%
    mutate(marker = marker_name)
}

long331 <- make_long(seqtab_core263, "331")
long263 <- make_long(seqtab263, "263")


#  ---------- 2. Relative abundance within sample ----------


long331 <- long331 %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = reads / sum(reads)) %>%
  ungroup()

long263 <- long263 %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = reads / sum(reads)) %>%
  ungroup()


# ---------- 3. Core-level summary per marker ----------

summ331 <- long331 %>%
  group_by(core263) %>%
  summarise(
    marker = "331",
    total_reads = sum(reads),
    prevalence = n_distinct(sample),
    max_rel_abund = max(rel_abund_sample),
    mean_rel_abund = mean(rel_abund_sample),
    .groups = "drop"
  )

summ263 <- long263 %>%
  group_by(core263) %>%
  summarise(
    marker = "263",
    total_reads = sum(reads),
    prevalence = n_distinct(sample),
    max_rel_abund = max(rel_abund_sample),
    mean_rel_abund = mean(rel_abund_sample),
    .groups = "drop"
  )

core_summary <- bind_rows(summ331, summ263)


# ---------- 4. Shared vs marker-specific ----------


shared_core <- intersect(summ331$core263, summ263$core263)

core_summary <- core_summary %>%
  mutate(
    shared_status = ifelse(core263 %in% shared_core, "shared", "marker_specific"),
    shared_status = factor(shared_status, levels = c("marker_specific", "shared"))
  )



# ---------- 5. Plot: black-and-white version ----------


p1 <- ggplot(core_summary, aes(x = max_rel_abund, y = prevalence)) +
  geom_point(
    colour = "black",
    alpha = 0.35,
    size = 1.8
  ) +
  scale_x_log10() +
  facet_wrap(~ shared_status, nrow = 1) +
  labs(
    x = "Maximum relative abundance",
    y = "Prevalence (number of samples)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_line(linewidth = 0.2),
    panel.spacing.x = unit(1.2, "lines")
  )

print(p1)

# ---------- 6. Save outputs ----------
write.csv(core_summary, file.path(outdir, "FigS3_prevalence_abundance_core_summary.csv"), row.names = FALSE)

ggsave(file.path(outdir, "FigS3_shared_specific_prev_abd_bw.png"), p1, width = 8, height = 6, dpi = 300)
