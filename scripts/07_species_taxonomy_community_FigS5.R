############################################################
# species-level taxonomy and community analyses:
# 1) reference preparation
# 2) taxonomy + traits
# 3) helper functions
# 4) species composition per marker comparison 
# 5) outputs used by downstream IPS script

# Required inputs:
# - fromAliView.fasta
#   Manually curated rbcL reference file prepared by trimming the Diat.Barcode
#   v15.2 alignment to the 331-bp target region in AliView; included in the
#   GitHub repository folder "input".
#
# - 2025-09-04-Diat.barcode_release-version 15.2.xlsx
#   Diat.Barcode v15.2 release file, available from:
#   https://doi.org/10.15454/TOMBYZ
#
# - diat_barcode_v15_2_tax_assign_dada2.fa
#   DADA2-formatted taxonomic reference file derived from Diat.Barcode;
#   available from:
#   https://doi.org/10.57745/SE6GJH
#############################################################

library(dada2)
library(readxl)
library(Biostrings)
library(DECIPHER)
library(stringr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)
library(vegan)
library(ggrepel)


# 0. PATHS AND OUTPUTS

project_dir <- normalizePath(".", mustWork = TRUE)
dada_dir <- file.path(project_dir, "DADA2_outputs")
outdir <- file.path(project_dir, "species-level_analyses")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

diatbarcode_xlsx <- file.path(project_dir, "input", "2025-09-04-Diat.barcode_release-version 15.2.xlsx")
aliview_curated_fasta <- file.path(project_dir, "input", "fromAliView.fasta")
dada2_tax_reference_fasta <- file.path( project_dir, "input", "diat_barcode_v15_2_tax_assign_dada2.fa")

out_aliview <- file.path(outdir, "forAliView.fasta")
out_ref <- file.path(outdir, "reference.csv")

out_assTax331 <- file.path(outdir, "assTax_331.fasta")
out_assTax263 <- file.path(outdir, "assTax_263.fasta")

out_core_map <- file.path(outdir, "331_to_core263_mapping.csv")


# 1. REFERENCE PREPARATION
# ---------- 1.1 Import barcode table ----------
dat <- read_excel(
  diatbarcode_xlsx,
  sheet = "sequences_info",
  col_types = "text"
)

# aligned rbcL sequences form column 10
seqs <- dat[[10]]

keep <- !is.na(seqs) & nzchar(seqs)
dat2 <- dat[keep, ]

dna <- DNAStringSet(toupper(seqs[keep]))

# sequence names: curated species name + accession
ids <- paste0(dat2[[31]], "_", dat2[[9]])
names(dna) <- ids

# export for manual curation in AliView
writeXStringSet(dna, out_aliview)

# ---------- Manual curation in AliView  ----------

# Sequence ends close to 331R binding site were manually inspected and problematic terminal regions edited, particularly in sequences originating from Medlin et al. (1988), where sequencing errors near sequence termini were common. 5' end, including the 331F primer region, and the 3' end, including the 331R primer region, were removed. All gaps were removed and sequences shorter than 331 bp were discarded.

# ---------- 1.2 Read manually curated 331 bp alignment  ----------
aln <- readDNAStringSet(aliview_curated_fasta)

seq_names   <- names(aln)
seq_strings <- as.character(aln)


# ---------- 1.3 Parse sequence names ----------
parts_list <- strsplit(seq_names, "_", fixed = TRUE)

accession <- sapply(parts_list, function(x) x[length(x)])
genus     <- sapply(parts_list, function(x) x[1])

# everything between genus and accession = species part
species_raw <- sapply(parts_list, function(x) {
  paste(x[2:(length(x) - 1)], collapse = "_")
})

species <- gsub("_", " ", species_raw)
final_taxonomy <- paste(genus, species)

meta <- data.frame(
  original_name  = seq_names,
  sequence       = seq_strings,
  accession      = accession,
  genus          = genus,
  species        = species,
  final_taxonomy = final_taxonomy,
  stringsAsFactors = FALSE
)


# ---------- 1.4 Import genus taxonomy from dada2 reference ----------
 
ref <- readDNAStringSet(dada2_tax_reference_fasta)
ref_names <- names(ref)

ref_split <- strsplit(ref_names, ";", fixed = TRUE)
ref_split <- lapply(ref_split, function(x) x[x != ""])

ref_genus <- sapply(ref_split, function(x) x[length(x) - 1])
ref_tax_no_species <- sapply(ref_split, function(x) paste(x[-length(x)], collapse = ";"))

map_gentax <- data.frame(
  genus = ref_genus,
  taxonomy = ref_tax_no_species,
  stringsAsFactors = FALSE
) %>%
  distinct()


# ---------- 1.5 Merge taxonomy with curated reference ----------

meta2 <- merge(
  meta,
  map_gentax,
  by = "genus",
  all.x = TRUE,
  sort = FALSE
) %>%
  mutate(
    final_taxonomy = final_taxonomy %>%
      str_replace_all("\\baff\\.?\\b", "cf.") %>%
      str_replace_all("group([0-9]+)", "group \\1") %>%
      str_squish(),
    taxname = paste0(taxonomy, ";", final_taxonomy, ";"),
    sequence_263 = substr(sequence, 331 - 263 + 1, 331)
  )

# ---------- 1.6 Export FASTA references ----------

seq331_tax <- DNAStringSet(meta2$sequence)
names(seq331_tax) <- meta2$taxname
writeXStringSet(seq331_tax, out_assTax331)

seq263_tax <- DNAStringSet(meta2$sequence_263)
names(seq263_tax) <- meta2$taxname
writeXStringSet(seq263_tax, out_assTax263)


# 1.7 Reference exact-match lookup tables
# NOTE:
# These are NOT the main taxonomic assignments used downstream.
# They are helper lookup tables used only for diagnostics:
#   - exact sequence match checks
#   - ambiguous reference matches (same sequence linked to >1 species)
#   - distinguishing "no exact match" vs "multiple exact species"
#
# Species_type here distinguishes only:
#   - species
#   - explicit_sp (already present in reference as "Genus sp.")
############################
normalize_species_name <- function(x) {
  x %>%
    str_replace_all(",", ".") %>%          
    str_replace_all("\\baff\\.?\\b", "cf.") %>%
    str_replace_all("group([0-9]+)", "group \\1") %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\bsp\\b\\.?", "sp.") %>%  
    str_squish()
}

ref331_lookup <- meta2 %>%
  transmute(
    sequence = sequence,
    Species = normalize_species_name(final_taxonomy),
    Species_type = ifelse(
      grepl(" sp\\.$", str_squish(final_taxonomy)),
      "explicit_sp",
      "species"
    )
  ) %>%
  distinct()

ref263_lookup <- meta2 %>%
  transmute(
    sequence = sequence_263,
    Species = normalize_species_name(final_taxonomy),
    Species_type = ifelse(
      grepl(" sp\\.$", str_squish(final_taxonomy)),
      "explicit_sp",
      "species"
    )
  ) %>%
  distinct()

# exact-match summaries for diagnostics
ref331_exact_summary <- ref331_lookup %>%
  group_by(sequence) %>%
  summarise(
    ref_n_species_331 = n_distinct(Species),
    ref_species_list_331 = paste(sort(unique(Species)), collapse = " | "),
    ref_species_type_list_331 = paste(sort(unique(Species_type)), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    ref_match_status_331 = case_when(
      ref_n_species_331 == 1 ~ "unique_exact_match",
      ref_n_species_331 > 1 ~ "multiple_exact_species",
      TRUE ~ "no_exact_match"
    )
  )

ref263_exact_summary <- ref263_lookup %>%
  group_by(sequence) %>%
  summarise(
    ref_n_species_263 = n_distinct(Species),
    ref_species_list_263 = paste(sort(unique(Species)), collapse = " | "),
    ref_species_type_list_263 = paste(sort(unique(Species_type)), collapse = " | "),
    .groups = "drop"
  ) %>%
  mutate(
    ref_match_status_263 = case_when(
      ref_n_species_263 == 1 ~ "unique_exact_match",
      ref_n_species_263 > 1 ~ "multiple_exact_species",
      TRUE ~ "no_exact_match"
    )
  )



# 2. TAXONOMY QC + TRAITS

# ----------  2.1 Read seqtabs ---------- 

seqtab331 <- readRDS(file.path(dada_dir, "seqtab_read_red331.rds"))
seqtab263 <- readRDS(file.path(dada_dir, "seqtab_read_red263.rds"))

# derive 263-bp cores from full 331-bp ASV
seqs331 <- colnames(seqtab331)
core263 <- substr(seqs331, nchar(seqs331) - 263 + 1, nchar(seqs331))
unique_core <- unique(core263)

seqtab_core263 <- sapply(unique_core, function(core) {
  cols <- which(core263 == core)
  rowSums(seqtab331[, cols, drop = FALSE])
}, simplify = "matrix")

seqtab_core263 <- as.matrix(seqtab_core263)
colnames(seqtab_core263) <- unique_core
rownames(seqtab_core263) <- rownames(seqtab331)


#  ---------- 2.2 assignTaxonomy  ----------
# 331 is classified separately against 331-bp reference
# 263 and core263 are classified together in one shared run against 263-bp reference (identical 263-bp ASVs get identical names)

normalize_species_name <- function(x) {
  x %>%
    str_replace_all(",", ".") %>%          
    str_replace_all("\\baff\\.?\\b", "cf.") %>%
    str_replace_all("group([0-9]+)", "group \\1") %>%
    str_replace_all("\\s+", " ") %>%
    str_replace_all("\\bsp\\b\\.?", "sp.") %>%  
    str_squish()
}

make_tax_qc <- function(seqs, ref_fasta) {
  tax_out <- assignTaxonomy(
    seqs,
    ref_fasta,
    minBoot = 85,
    outputBootstraps = TRUE,
    taxLevels = c("Empire", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
  
  tax_df <- as.data.frame(tax_out$tax) %>%
    rownames_to_column("sequence")
  
  boot_df <- as.data.frame(tax_out$boot) %>%
    rownames_to_column("sequence") %>%
    rename_with(~ paste0(.x, "_boot"), -sequence)
  
  tax_df %>%
    left_join(boot_df, by = "sequence") %>%
    mutate(
      Species_original = Species,
      Species_original_norm = normalize_species_name(Species_original),
      
      Species = case_when(
        is.na(Species_original) & !is.na(Genus) ~ paste0(Genus, " sp."),
        TRUE ~ Species_original
      ),
      Species = normalize_species_name(Species)
      )
}

# 331 assignment
set.seed(1234)
tax331 <- make_tax_qc(colnames(seqtab331), out_assTax331)

# 263 + core263 shared assignment
seqs_263_common <- unique(c(colnames(seqtab263), colnames(seqtab_core263)))

set.seed(1234)
tax263_common <- make_tax_qc(seqs_263_common, out_assTax263)


saveRDS(tax331,       file.path(outdir, "tax331.rds"))
saveRDS(tax263_common,file.path(outdir, "tax263_common.rds"))


tax331_df <- tax331 %>%
  transmute(
    sequence,
    Species,
    Genus,
    Genus_boot,
    Species_boot
  ) %>%
  distinct()


tax263_common_df <- tax263_common %>%
  transmute(
    sequence,
    Species,
    Genus,
    Genus_boot,
    Species_boot
  ) %>%
  distinct()

# outputs used by downstream IPS script
saveRDS(tax331_df, file.path(outdir, "tax331_df.rds"))
saveRDS(tax263_common_df, file.path(outdir, "tax263_common_df.rds"))

# 331 -> core263 mapping table
map331_core263 <- tibble(
  sequence_331 = colnames(seqtab331),
  sequence_core263 = substr(colnames(seqtab331), nchar(colnames(seqtab331)) - 263 + 1, nchar(colnames(seqtab331)))
) %>%
  left_join(
    tax331_df %>%
      rename(
        sequence_331 = sequence,
        Species_331 = Species,
        Genus_331 = Genus,
        Genus_boot_331 = Genus_boot,
        Species_boot_331 = Species_boot
      ),
    by = "sequence_331"
  ) %>%
  left_join(
    tax263_common_df %>%
      rename(
        sequence_core263 = sequence,
        Species_core263_common = Species,
        Genus_core263_common = Genus,
        Genus_boot_core263_common = Genus_boot,
        Species_boot_core263_common = Species_boot
      ),
    by = "sequence_core263"
  ) %>%
  left_join(
    ref331_exact_summary %>%
      rename(sequence_331 = sequence),
    by = "sequence_331"
  ) %>%
  left_join(
    ref263_exact_summary %>%
      rename(sequence_core263 = sequence),
    by = "sequence_core263"
  ) %>%
  mutate(
    ref_match_status_331 = ifelse(is.na(ref_match_status_331), "no_exact_match", ref_match_status_331),
    ref_match_status_263 = ifelse(is.na(ref_match_status_263), "no_exact_match", ref_match_status_263),
    ref_n_species_331 = ifelse(is.na(ref_n_species_331), 0L, ref_n_species_331),
    ref_n_species_263 = ifelse(is.na(ref_n_species_263), 0L, ref_n_species_263),
    ref_species_list_331 = ifelse(is.na(ref_species_list_331), "", ref_species_list_331),
    ref_species_list_263 = ifelse(is.na(ref_species_list_263), "", ref_species_list_263)
  ) %>%
  distinct()


#  ---------- 2.3 Import traits ---------- 

traits <- read_excel(
  diatbarcode_xlsx,
  sheet = "traits",
  col_types = "text"
)%>%
  transmute(
    Species = normalize_species_name(Species),
    IPSS  = parse_number(as.character(IPSS)),
    IPSVV = parse_number(as.character(IPSVV)),
    CFv2  = parse_number(as.character(`CF v2`))
  ) %>%
  filter(!is.na(IPSS), !is.na(IPSVV), !is.na(CFv2))


# 3. HELPER FUNCTIONS

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

make_taxon_reads <- function(dat) {
  dat %>%
    filter(!is.na(Species)) %>%
    group_by(sample, Species) %>%
    summarise(reads = sum(reads), .groups = "drop")
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


# 4. SPECIES COMPOSITION COMPARISON

# ----------  4.1 Long tables ----------

seq331_long      <- make_long_from_seqtab(seqtab331)
seq263_long      <- make_long_from_seqtab(seqtab263)
seq_core263_long <- make_long_from_seqtab(seqtab_core263)


# ----------  4.2 Species-level tables without traits filter ---------- 

dat331_species <- make_species_table(seq331_long, tax331_df, traits_df = NULL)
dat263_species <- make_species_table(seq263_long, tax263_common_df, traits_df = NULL)

dat331_species_taxon <- make_taxon_reads(dat331_species) %>%
  mutate(marker = "331")
dat263_species_taxon <- make_taxon_reads(dat263_species) %>%
  mutate(marker = "263")

comp_long <- bind_rows(dat331_species_taxon, dat263_species_taxon) %>%
  mutate(sample_marker = paste(sample, marker, sep = "_"))

meta_comp <- comp_long %>%
  distinct(sample_marker, sample, marker)


# ----------  4.3 Presence / absence - 331 vs. 263 ----------

pa_mat <- comp_long %>%
  mutate(pa = 1L) %>%
  distinct(sample_marker, Species, .keep_all = TRUE) %>%
  select(sample_marker, Species, pa) %>%
  pivot_wider(
    names_from = Species,
    values_from = pa,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(pa_mat) <- pa_mat$sample_marker
pa_mat$sample_marker <- NULL

meta_pa <- meta_comp[match(rownames(pa_mat), meta_comp$sample_marker), ]

dist_pa <- vegdist(pa_mat, method = "jaccard", binary = TRUE)
adonis_pa <- adonis2(dist_pa ~ marker, data = meta_pa, strata = meta_pa$sample)
print(adonis_pa)

dist_pa_mat <- as.matrix(dist_pa)

pair_pa <- meta_pa %>%
  select(sample, marker, sample_marker) %>%
  distinct() %>%
  group_by(sample) %>%
  summarise(
    sample_331 = sample_marker[marker == "331"][1],
    sample_263 = sample_marker[marker == "263"][1],
    jaccard_pair = dist_pa_mat[sample_331, sample_263],
    .groups = "drop"
  ) %>%
  drop_na()

print(pair_pa)
summary(pair_pa$jaccard_pair)

all_dist_pa <- as.vector(dist_pa)
wilcox.test(pair_pa$jaccard_pair, all_dist_pa, alternative = "less", exact = FALSE)

set.seed(123)
nmds_pa <- metaMDS(pa_mat, distance = "jaccard", binary = TRUE, trace = FALSE)

scores_pa <- as.data.frame(scores(nmds_pa, display = "sites"))
scores_pa$sample_marker <- rownames(scores_pa)

plot_pa <- scores_pa %>%
  left_join(meta_pa, by = "sample_marker")

p_pa <- ggplot(plot_pa, aes(NMDS1, NMDS2, color = marker)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample), size = 3) +
  theme_bw(base_size = 14) +
  labs(
    title = "NMDS - presence/absence",
    subtitle = "Jaccard distance"
  )

print(p_pa)


# ---------- 4.4 Relative abundance - 331 vs. 263 ----------

rel_wide <- comp_long %>%
  select(sample_marker, Species, reads) %>%
  pivot_wider(
    names_from = Species,
    values_from = reads,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(rel_wide) <- rel_wide$sample_marker
rel_wide$sample_marker <- NULL

rel_mat <- decostand(rel_wide, method = "total")
meta_rel <- meta_comp[match(rownames(rel_mat), meta_comp$sample_marker), ]

dist_rel <- vegdist(rel_mat, method = "bray")
adonis_rel <- adonis2(dist_rel ~ marker, data = meta_rel, strata = meta_rel$sample)
print(adonis_rel)

dist_rel_mat <- as.matrix(dist_rel)

pair_rel <- meta_rel %>%
  select(sample, marker, sample_marker) %>%
  distinct() %>%
  group_by(sample) %>%
  summarise(
    sample_331 = sample_marker[marker == "331"][1],
    sample_263 = sample_marker[marker == "263"][1],
    bray_pair = dist_rel_mat[sample_331, sample_263],
    .groups = "drop"
  ) %>%
  drop_na()

print(pair_rel)
summary(pair_rel$bray_pair)

all_dist_rel <- as.vector(dist_rel)
wilcox.test(pair_rel$bray_pair, all_dist_rel, alternative = "less", exact = FALSE)

set.seed(123)
nmds_rel <- metaMDS(rel_mat, distance = "bray", trace = FALSE)

scores_rel <- as.data.frame(scores(nmds_rel, display = "sites"))
scores_rel$sample_marker <- rownames(scores_rel)

plot_rel <- scores_rel %>%
  left_join(meta_rel, by = "sample_marker")

p_rel <- ggplot(plot_rel, aes(NMDS1, NMDS2, color = marker)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sample), size = 3) +
  theme_bw(base_size = 14) +
  labs(
    title = "NMDS - relative abundance",
    subtitle = "Bray-Curtis distance"
  )

print(p_rel)



# 5. SAVE OUTPUTS
# The IPS calculation and decomposition are performed in 08_IPS_analysis_Fig2_Fig3.R.

write.csv(meta2, out_ref, row.names = FALSE)
write.csv(map331_core263, out_core_map, row.names = FALSE)

write.csv(comp_long, file.path(outdir, "species_composition_long.csv"), row.names = FALSE)
write.csv(pair_pa, file.path(outdir, "pairwise_jaccard_same_sample.csv"), row.names = FALSE)
write.csv(pair_rel, file.path(outdir, "pairwise_bray_same_sample.csv"), row.names = FALSE)

saveRDS(pa_mat, file.path(outdir, "pa_matrix.rds"))
saveRDS(rel_mat, file.path(outdir, "rel_matrix.rds"))
saveRDS(dist_pa, file.path(outdir, "dist_jaccard_pa.rds"))
saveRDS(dist_rel, file.path(outdir, "dist_bray_rel.rds"))
saveRDS(nmds_pa, file.path(outdir, "nmds_pa.rds"))
saveRDS(nmds_rel, file.path(outdir, "nmds_rel.rds"))

ggsave(file.path(outdir, "FigS5_NMDS_presence_absence.png"), p_pa, width = 7, height = 6, dpi = 300)
ggsave(file.path(outdir, "FigS5_NMDS_relative_abundance.png"), p_rel, width = 7, height = 6, dpi = 300)