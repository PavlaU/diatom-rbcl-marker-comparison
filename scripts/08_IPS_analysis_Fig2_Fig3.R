############################################################
# IPS analysis and IPS decomposition
#
# Required inputs:   
# - 2025-09-04-Diat.barcode_release-version 15.2.xlsx
#   Diat.Barcode v15.2 release file, available from:
#   https://doi.org/10.15454/TOMBYZ
############################################################

library(readxl)
library(stringr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)

# ---------- 0. Paths and load inputs ----------
project_dir <- normalizePath(".", mustWork = TRUE)
outdir <- file.path(project_dir, "ips-level_analyses")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Local copy of Diat.Barcode v15.2.
diatbarcode_xlsx <- file.path(project_dir, "input", "2025-09-04-Diat.barcode_release-version 15.2.xlsx")

seqtab263 <- readRDS(file.path(project_dir, "DADA2_outputs", "seqtab_read_red263.rds"))
seqtab331 <- readRDS(file.path(project_dir, "DADA2_outputs", "seqtab_read_red331.rds"))

seqtab_core263 <- readRDS(file.path(project_dir, "ASV_correlation", "seqtab_core263_from_331.rds"))

tax331_df <- readRDS(file.path(project_dir, "species-level_analyses", "tax331_df.rds"))
tax263_common_df <- readRDS(file.path(project_dir, "species-level_analyses", "tax263_common_df.rds"))

# ---------- 0.1 Helper functions ----------

normalize_species_name <- function(x) {
  x %>%
    str_replace_all(",", ".") %>%
    str_replace_all("\\baff\\.?\\b", "cf.") %>%
    str_replace_all("group([0-9]+)", "group \\1") %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\bsp\\b\\.?", "sp.") %>%
    str_squish()
}

make_long_from_seqtab <- function(seqtab_obj) {
  seqtab_obj %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(
      cols = -sample,
      names_to = "sequence",
      values_to = "reads"
    ) %>%
    filter(reads > 0)
}

make_species_table <- function(seq_long, tax_df, traits_df = NULL) {
  out <- seq_long %>%
    left_join(tax_df, by = "sequence")

  if (!is.null(traits_df)) {
    out <- out %>%
      left_join(traits_df, by = "Species")
  }

  out
}

make_taxon_reads_traits <- function(dat) {
  dat %>%
    filter(!is.na(Species), !is.na(IPSS), !is.na(IPSVV), !is.na(CFv2)) %>%
    group_by(sample, Species, IPSS, IPSVV, CFv2) %>%
    summarise(reads = sum(reads), .groups = "drop")
}

cf_correct <- function(df) {
  df %>%
    mutate(Reads_CF = reads / CFv2) %>%
    group_by(sample) %>%
    mutate(RelAbund_CF = Reads_CF / sum(Reads_CF, na.rm = TRUE)) %>%
    ungroup()
}

calc_ips <- function(df, ips_col = "IPS_20") {
  df %>%
    group_by(sample) %>%
    summarise(
      IPS = sum(RelAbund_CF * IPSS * IPSVV, na.rm = TRUE) /
        sum(RelAbund_CF * IPSVV, na.rm = TRUE),
      n_taxa = n_distinct(Species),
      .groups = "drop"
    ) %>%
    mutate(
      "{ips_col}" := 4.75 * IPS - 3.75
    )
}

assign_ips_class <- function(x) {
  case_when(
    x >= 17 & x <= 20 ~ "High",
    x >= 13 & x < 17 ~ "Good",
    x >= 9  & x < 13 ~ "Moderate",
    x >= 5  & x < 9  ~ "Poor",
    x >= 1  & x < 5  ~ "Bad",
    TRUE ~ NA_character_
  )
}

# ---------- 1. Import traits ----------

traits <- read_excel(
  diatbarcode_xlsx,
  sheet = "traits",
  col_types = "text"
) %>%
  transmute(
    Species = normalize_species_name(Species),
    IPSS  = parse_number(as.character(IPSS)),
    IPSVV = parse_number(as.character(IPSVV)),
    CFv2  = parse_number(as.character(`CF v2`))
  ) %>%
  filter(!is.na(IPSS), !is.na(IPSVV), !is.na(CFv2))

# ---------- 2. Long ASV tables ----------

seq331_long <- make_long_from_seqtab(seqtab331)
seq263_long <- make_long_from_seqtab(seqtab263)
seq_core263_long <- make_long_from_seqtab(seqtab_core263)

# ---------- 3. Join species and traits ----------

dat331 <- make_species_table(seq331_long, tax331_df, traits)
dat263 <- make_species_table(seq263_long, tax263_common_df, traits)
dat_core263 <- make_species_table(seq_core263_long, tax263_common_df, traits)

# ---------- 4. Collapse ASVs to species ----------

dat331_taxon <- make_taxon_reads_traits(dat331)
dat263_taxon <- make_taxon_reads_traits(dat263)
dat_core263_taxon <- make_taxon_reads_traits(dat_core263)

# ---------- 5. CF correction ----------

dat331_taxon <- cf_correct(dat331_taxon)
dat263_taxon <- cf_correct(dat263_taxon)
dat_core263_taxon <- cf_correct(dat_core263_taxon)

# ---------- 6. Full IPS ----------

ips331_all <- calc_ips(dat331_taxon, "IPS_20_331_all")
ips263_all <- calc_ips(dat263_taxon, "IPS_20_263_all")
ips263_all_core <- calc_ips(dat_core263_taxon, "IPS_20_263core_all")

ips_compare <- ips331_all %>%
  select(sample, IPS_20_331_all) %>%
  inner_join(
    ips263_all %>% select(sample, IPS_20_263_all),
    by = "sample"
  ) %>%
  mutate(
    class_331 = assign_ips_class(IPS_20_331_all),
    class_263 = assign_ips_class(IPS_20_263_all),
    class_changed = class_331 != class_263
  )

ips_core_compare <- ips263_all_core %>%
  select(sample, IPS_20_263core_all) %>%
  inner_join(
    ips263_all %>% select(sample, IPS_20_263_all),
    by = "sample"
  )

# ---------- 7. Statistics ----------

ips_wilcox <- wilcox.test(
  ips_compare$IPS_20_331_all,
  ips_compare$IPS_20_263_all,
  paired = TRUE,
  exact = FALSE
)
print(ips_wilcox)

cor_test <- cor.test(
  ips_compare$IPS_20_331_all,
  ips_compare$IPS_20_263_all
)
print(cor_test)

lm_fit <- lm(IPS_20_263_all ~ IPS_20_331_all, data = ips_compare)
print(summary(lm_fit))

ips_wilcox_core <- wilcox.test(
  ips_core_compare$IPS_20_263core_all,
  ips_core_compare$IPS_20_263_all,
  paired = TRUE,
  exact = FALSE
)
print(ips_wilcox_core)

cor_test_core <- cor.test(
  ips_core_compare$IPS_20_263core_all,
  ips_core_compare$IPS_20_263_all
)
print(cor_test_core)

lm_fit_core <- lm(IPS_20_263_all ~ IPS_20_263core_all, data = ips_core_compare)
print(summary(lm_fit_core))

# ---------- 8. Plots: full IPS comparisons ----------

label_txt <- paste0(
  "R = ", round(cor_test$estimate, 2),
  ", p = ", format.pval(cor_test$p.value, digits = 2), "\n",
  "y = ", round(coef(lm_fit)[1], 3), " + ",
  round(coef(lm_fit)[2], 2), "x"
)

p_ips <- ggplot(ips_compare, aes(x = IPS_20_331_all, y = IPS_20_263_all)) +
  annotate("rect", xmin = 1,  xmax = 5,  ymin = 1,  ymax = 5,  fill = "#e41a1c", alpha = 0.35) +
  annotate("rect", xmin = 5,  xmax = 9,  ymin = 5,  ymax = 9,  fill = "#ff7f00", alpha = 0.35) +
  annotate("rect", xmin = 9,  xmax = 13, ymin = 9,  ymax = 13, fill = "#ffd92f", alpha = 0.35) +
  annotate("rect", xmin = 13, xmax = 17, ymin = 13, ymax = 17, fill = "#4daf4a", alpha = 0.35) +
  annotate("rect", xmin = 17, xmax = 20, ymin = 17, ymax = 20, fill = "#377eb8", alpha = 0.35) +
  geom_point(size = 1.8, alpha = 0.85, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  annotate("text", x = 2, y = 19, label = label_txt, hjust = 0, size = 4.5) +
  scale_x_continuous(limits = c(1, 20), breaks = seq(1, 20, 2)) +
  scale_y_continuous(limits = c(1, 20), breaks = seq(1, 20, 2)) +
  labs(
    x = "IPS score 331-bp (Kelly)",
    y = "IPS score 263-bp (Vasselon)"
  ) +
  theme_classic(base_size = 14)

print(p_ips)

label_txt_core <- paste0(
  "R = ", round(cor_test_core$estimate, 2),
  ", p = ", format.pval(cor_test_core$p.value, digits = 2), "\n",
  "y = ", round(coef(lm_fit_core)[1], 3), " + ",
  round(coef(lm_fit_core)[2], 2), "x"
)

p_ips_core <- ggplot(ips_core_compare, aes(x = IPS_20_263core_all, y = IPS_20_263_all)) +
  annotate("rect", xmin = 1,  xmax = 5,  ymin = 1,  ymax = 5,  fill = "#e41a1c", alpha = 0.35) +
  annotate("rect", xmin = 5,  xmax = 9,  ymin = 5,  ymax = 9,  fill = "#ff7f00", alpha = 0.35) +
  annotate("rect", xmin = 9,  xmax = 13, ymin = 9,  ymax = 13, fill = "#ffd92f", alpha = 0.35) +
  annotate("rect", xmin = 13, xmax = 17, ymin = 13, ymax = 17, fill = "#4daf4a", alpha = 0.35) +
  annotate("rect", xmin = 17, xmax = 20, ymin = 17, ymax = 20, fill = "#377eb8", alpha = 0.35) +
  geom_point(size = 1.8, alpha = 0.85, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  annotate("text", x = 2, y = 19, label = label_txt_core, hjust = 0, size = 4.5) +
  scale_x_continuous(limits = c(1, 20), breaks = seq(1, 20, 2)) +
  scale_y_continuous(limits = c(1, 20), breaks = seq(1, 20, 2)) +
  labs(
    x = "IPS score 263core (Kelly)",
    y = "IPS score 263-bp (Vasselon)"
  ) +
  theme_classic(base_size = 14)

print(p_ips_core)

# ---------- 9. IPS decomposition ----------

species_shared_core <- intersect(
  unique(dat263_taxon$Species[!is.na(dat263_taxon$Species)]),
  unique(dat_core263_taxon$Species[!is.na(dat_core263_taxon$Species)])
)

# 9.1 Taxa shared globally by 263 and 263core
ips263_shared <- dat263_taxon %>%
  filter(Species %in% species_shared_core)

ips263_core_shared <- dat_core263_taxon %>%
  filter(Species %in% species_shared_core)

ips263_shared <- calc_ips(ips263_shared, "IPS_20_263_shared")
ips263_core_shared <- calc_ips(ips263_core_shared, "IPS_20_263core_shared")

ips_core_compare_share <- ips263_core_shared %>%
  select(sample, IPS_20_263core_shared) %>%
  inner_join(
    ips263_shared %>% select(sample, IPS_20_263_shared),
    by = "sample"
  )

# 9.2 Species shared within each sample
shared_species_by_sample <- inner_join(
  dat_core263_taxon %>%
    distinct(sample, Species) %>%
    mutate(in_263core = TRUE),
  dat263_taxon %>%
    distinct(sample, Species) %>%
    mutate(in_263 = TRUE),
  by = c("sample", "Species")
) %>%
  select(sample, Species)

ips263core_sample_shared <- dat_core263_taxon %>%
  semi_join(shared_species_by_sample, by = c("sample", "Species"))

ips263_sample_shared <- dat263_taxon %>%
  semi_join(shared_species_by_sample, by = c("sample", "Species"))

ips263core_sample_shared <- calc_ips(ips263core_sample_shared, "IPS_20_263core_sample_shared")
ips263_sample_shared <- calc_ips(ips263_sample_shared, "IPS_20_263_sample_shared")

ips_core_sample_compare <- ips263core_sample_shared %>%
  select(sample, IPS_20_263core_sample_shared) %>%
  inner_join(
    ips263_sample_shared %>% select(sample, IPS_20_263_sample_shared),
    by = "sample"
  )

# 9.3 Combine IPS levels
compare_all <- ips_compare %>%
  select(
    sample,
    IPS_20_331_all,
    IPS_20_263_all
  ) %>%
  left_join(
    ips_core_compare %>%
      select(
        sample,
        IPS_20_263core_all
      ),
    by = "sample"
  ) %>%
  left_join(
    ips_core_compare_share %>%
      select(
        sample,
        IPS_20_263_shared,
        IPS_20_263core_shared
      ),
    by = "sample"
  ) %>%
  left_join(
    ips_core_sample_compare %>%
      select(
        sample,
        IPS_20_263core_sample_shared,
        IPS_20_263_sample_shared
      ),
    by = "sample"
  )

compare_all <- compare_all %>%
  mutate(
    # 1) full difference: 263 vs 331
    delta_full_263_vs_331 = IPS_20_263_all - IPS_20_331_all,

    # 2) length/annotation effect: 263 vs 263core
    delta_263_vs_263core = IPS_20_263_all - IPS_20_263core_all,

    # 3) unique taxa effect: 263 vs 263core species found by both marker
    delta_263_unique_effect = IPS_20_263_shared - IPS_20_263core_shared,

    # 4) read abundance effect: species shared within sample
    delta_shared_abundance_effect = IPS_20_263_sample_shared - IPS_20_263core_sample_shared
  )

# 9.4 IPS decomposition boxplot
plot_delta_df <- compare_all %>%
  select(
    sample,
    delta_full_263_vs_331,
    delta_263_vs_263core,
    delta_263_unique_effect,
    delta_shared_abundance_effect
  ) %>%
  pivot_longer(
    cols = -sample,
    names_to = "comparison",
    values_to = "delta"
  ) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "delta_full_263_vs_331",
        "delta_263_vs_263core",
        "delta_263_unique_effect",
        "delta_shared_abundance_effect"
      ),
      labels = c("A", "B", "C", "D")
    )
  )

ips_box <- ggplot(plot_delta_df, aes(x = comparison, y = delta)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.7) +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Difference in IPS(20)"
  )

print(ips_box)

# ---------- 10. Save outputs ----------

write.csv(dat331_taxon, file.path(outdir, "ips_taxa_331.csv"), row.names = FALSE)
write.csv(dat263_taxon, file.path(outdir, "ips_taxa_263.csv"), row.names = FALSE)
write.csv(dat_core263_taxon, file.path(outdir, "ips_taxa_263core.csv"), row.names = FALSE)

write.csv(ips_compare, file.path(outdir, "ips_compare_331_vs_263.csv"), row.names = FALSE)
write.csv(ips_core_compare, file.path(outdir, "ips_compare_263core_vs_263.csv"), row.names = FALSE)
write.csv(compare_all, file.path(outdir, "ips_decomposition_values.csv"), row.names = FALSE)
write.csv(plot_delta_df, file.path(outdir, "ips_decomposition_long.csv"), row.names = FALSE)

saveRDS(ips_wilcox, file.path(outdir, "ips_wilcox_331_vs_263.rds"))
saveRDS(cor_test, file.path(outdir, "ips_cor_331_vs_263.rds"))
saveRDS(lm_fit, file.path(outdir, "ips_lm_331_vs_263.rds"))
saveRDS(ips_wilcox_core, file.path(outdir, "ips_wilcox_263core_vs_263.rds"))
saveRDS(cor_test_core, file.path(outdir, "ips_cor_263core_vs_263.rds"))
saveRDS(lm_fit_core, file.path(outdir, "ips_lm_263core_vs_263.rds"))

ggsave(file.path(outdir, "Fig2_IPS_331_vs_263.png"), p_ips, width = 7, height = 6, dpi = 300)
ggsave(file.path(outdir, "IPS_263core_vs_263.png"), p_ips_core, width = 7, height = 6, dpi = 300)
ggsave(file.path(outdir, "Fig3_IPS_decomposition_boxplot.png"), ips_box, width = 7, height = 5, dpi = 300)