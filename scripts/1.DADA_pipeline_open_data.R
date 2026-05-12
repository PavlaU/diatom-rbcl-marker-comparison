#####################################################################
# 1. DADA2 ASV generation
# Chimera-free ASV tables with merged technical extractions
#
# Run twice:
#   marker_263 <- TRUE   for the 263-bp marker
#   marker_263 <- FALSE  for the 331-bp marker
#
# Outputs are written to ~/DADA2_outputs and are used by scripts 2-8.
#####################################################################

library(Biostrings)
library(dada2)
library(stringr)

set.seed(1234)

# ---------- Marker selection ----------
# Set TRUE for 263-bp marker, FALSE for 331-bp marker.
marker_263 <- FALSE

if (marker_263) {
  fwd_primer <- c(
    "AGGTGAAGTAAAAGGTTCWTACTTAAA",
    "AGGTGAAGTTAAAGGTTCWTAYTTAAA",
    "AGGTGAAACTAAAGGTTCWTACTTAAA"
  )
  rev_primer <- c(
    "CCTTCTAATTTACCWACWACTG",
    "CCTTCTAATTTACCWACAACAG"
  )
  trunc_len <- c(260, 240)
  exp_len <- 263:266
  red_len <- 263
} else {
  fwd_primer <- "ATGCGTTGGAGAGARCGTTTC"
  rev_primer <- "GATCACCTTCTAATTTACCWACAACTG"
  trunc_len <- c(280, 250)
  exp_len <- 331:334
  red_len <- 331
}

marker_label <- if (marker_263) "263" else "331"
message("Running marker: ", marker_label)

dir.create(path.expand("~/DADA2_outputs"), showWarnings = FALSE, recursive = TRUE)

# ---------- Cutadapt ----------
cutadapt_bin <- "/path/to/cutadapt-env/bin/cutadapt"

# ---------- Input FASTQ ----------
# Raw FASTQ files must be placed in ~/NGS_raw/

fwd_files <- sort(list.files(path.expand("~/NGS_raw"), pattern = "_1.fastq.gz$", full.names = TRUE))
rev_files <- sort(list.files(path.expand("~/NGS_raw"), pattern = "_2.fastq.gz$", full.names = TRUE))

samples <- str_extract(basename(fwd_files), "^[^_]+_[^_]+_[^_]+")
names(fwd_files) <- samples
names(rev_files) <- samples

if (marker_263) {
  # 263-bp marker: samples without "J" in sample ID
  keep_idx <- !grepl("J", samples)
} else {
  # 331-bp marker: samples with "J" in sample ID
  keep_idx <- grepl("J", samples)
}

fwd_files <- fwd_files[keep_idx]
rev_files <- rev_files[keep_idx]
samples <- samples[keep_idx]

# ---------- Primer trimming ----------
fwd_primer_rc <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rc <- as.character(reverseComplement(DNAStringSet(rev_primer)))

trim_dir <- file.path(path.expand("~/DADA2_outputs"), paste0("01_trimmed_", marker_label))
dir.create(trim_dir, showWarnings = FALSE, recursive = TRUE)

fwd_trimmed <- file.path(trim_dir, basename(fwd_files))
rev_trimmed <- file.path(trim_dir, basename(rev_files))

make_args <- function(flag, seqs) {
  unlist(lapply(seqs, function(s) c(flag, s)))
}

cutadapt_args <- c(
  make_args("-g", paste0("^", fwd_primer)),
  make_args("-a", rev_primer_rc),
  make_args("-G", paste0("^", rev_primer)),
  make_args("-A", fwd_primer_rc),
  "-n", "2",
  "--discard-untrimmed"
)

for (i in seq_along(samples)) {
  system2(
    cutadapt_bin,
    args = c(
      cutadapt_args,
      "-o", fwd_trimmed[i],
      "-p", rev_trimmed[i],
      fwd_files[i],
      rev_files[i]
    ),
    stdout = TRUE,
    stderr = TRUE
  )
}

# ---------- Split reads by sequencing run ----------
split_dir <- file.path(path.expand("~/DADA2_outputs"), paste0("02_split_by_run_", marker_label))
dir.create(split_dir, showWarnings = FALSE, recursive = TRUE)

get_run_name <- function(header) {
  header <- sub("^@", "", header)
  header <- strsplit(header, " ", fixed = TRUE)[[1]][1]
  parts <- strsplit(header, ":", fixed = TRUE)[[1]]
  paste(parts[1], parts[2], parts[3], sep = "_")
}

closeAllConnections()

for (file in list.files(trim_dir, pattern = "\\.fastq.gz$", full.names = TRUE)) {
  con_in <- gzfile(file, "rt")
  out_connections <- list()

  repeat {
    lines <- readLines(con_in, n = 4)
    if (length(lines) == 0) break
    if (length(lines) < 4) break

    run_name <- get_run_name(lines[1])

    if (!run_name %in% names(out_connections)) {
      run_folder <- file.path(split_dir, run_name)
      dir.create(run_folder, showWarnings = FALSE, recursive = TRUE)
      out_connections[[run_name]] <- gzfile(file.path(run_folder, basename(file)), "wt")
    }

    writeLines(lines, out_connections[[run_name]])
  }

  close(con_in)
  for (x in out_connections) close(x)
}

# ---------- Filtering and DADA2 inference by run ----------
get_sample_name <- function(x) {
  str_extract(basename(x), "^[^_]+_[^_]+_[^_]+")
}

seqtabs_by_run <- list()
filter_stats_by_run <- list()

filt_base_dir <- file.path(path.expand("~/DADA2_outputs"), paste0("03_filtered_by_run_", marker_label))
dir.create(filt_base_dir, showWarnings = FALSE, recursive = TRUE)

for (run_dir in list.dirs(split_dir, full.names = TRUE, recursive = FALSE)) {
  run_name <- basename(run_dir)
  message("Processing run: ", run_name)

  run_fwd <- sort(list.files(run_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
  run_rev <- sort(list.files(run_dir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))

  run_samples <- get_sample_name(run_fwd)
  names(run_fwd) <- run_samples
  names(run_rev) <- run_samples

  run_filt_dir <- file.path(filt_base_dir, run_name)
  dir.create(run_filt_dir, showWarnings = FALSE, recursive = TRUE)

  filtFs <- file.path(run_filt_dir, basename(run_fwd))
  filtRs <- file.path(run_filt_dir, basename(run_rev))

  filt_stats <- filterAndTrim(
    fwd = run_fwd,
    filt = filtFs,
    rev = run_rev,
    filt.rev = filtRs,
    truncLen = trunc_len,
    maxEE = c(2, 2),
    truncQ = 2,
    maxN = 0,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = TRUE
  )

  errF <- learnErrors(filtFs, multithread = TRUE)
  errR <- learnErrors(filtRs, multithread = TRUE)

  derepF <- derepFastq(filtFs)
  derepR <- derepFastq(filtRs)
  names(derepF) <- run_samples
  names(derepR) <- run_samples

  dadaF <- dada(derepF, err = errF, multithread = TRUE)
  dadaR <- dada(derepR, err = errR, multithread = TRUE)

  mergers <- mergePairs(dadaF, derepF, dadaR, derepR, maxMismatch = 1)
  seqtab_run <- makeSequenceTable(mergers)

  seqtabs_by_run[[run_name]] <- seqtab_run
  filter_stats_by_run[[run_name]] <- filt_stats

  saveRDS(filt_stats, file.path(run_filt_dir, paste0("filter_stats_", run_name, ".rds")))
  saveRDS(errF, file.path(run_filt_dir, paste0("errF_", run_name, ".rds")))
  saveRDS(errR, file.path(run_filt_dir, paste0("errR_", run_name, ".rds")))
  saveRDS(seqtab_run, file.path(run_filt_dir, paste0("seqtab_", run_name, ".rds")))
}

# ---------- Merge runs ----------
seqtab <- mergeSequenceTables(tables = seqtabs_by_run, repeats = "sum")
seqtab <- seqtab[, nchar(colnames(seqtab)) %in% exp_len, drop = FALSE]

# ---------- Chimera removal ----------
seqtab_nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE
)

# ---------- Merge technical replicates ----------
base_id <- sub("_[12]$", "", rownames(seqtab_nochim))
unique_ids <- unique(base_id)

seqtab_final <- matrix(
  0,
  nrow = length(unique_ids),
  ncol = ncol(seqtab_nochim),
  dimnames = list(unique_ids, colnames(seqtab_nochim))
)

for (id in unique_ids) {
  idx <- which(base_id == id)
  x1 <- seqtab_nochim[idx[1], ]
  x2 <- seqtab_nochim[idx[2], ]

  merged <- x1 + x2
  merged[merged < 2] <- 0
  seqtab_final[id, ] <- merged
}

# ---------- Keep expected marker length and remove extremely rare ASVs ----------
read_red <- seqtab_final[, nchar(colnames(seqtab_final)) == red_len, drop = FALSE]
abd_red <- read_red / rowSums(seqtab_final)

keep_asv <- apply(abd_red, 2, function(x) any(x >= 0.0001))
read_red_filt <- read_red[, keep_asv, drop = FALSE]

# ---------- Remove internal DNA standard sequences ----------
if (marker_263) {
  ids <- as.character(readDNAStringSet("~/DADA2_outputs/ids263.fasta"))
} else {
  ids <- as.character(readDNAStringSet("~/DADA2_outputs/ids331.fasta"))
}

is_ids_match <- function(asv, ids_ref, max_mismatch = 5) {
  hits <- ids_ref[nchar(ids_ref) == nchar(asv)]
  if (length(hits) == 0) return(FALSE)

  asv_chars <- strsplit(asv, "", fixed = TRUE)[[1]]

  any(vapply(hits, function(ref) {
    ref_chars <- strsplit(ref, "", fixed = TRUE)[[1]]
    sum(asv_chars != ref_chars) <= max_mismatch
  }, logical(1)))
}

ids_cols <- vapply(
  colnames(read_red_filt),
  is_ids_match,
  logical(1),
  ids_ref = ids,
  max_mismatch = 5
)

read_ids <- read_red_filt[, ids_cols, drop = FALSE]
read_red <- read_red_filt[, !ids_cols, drop = FALSE]

# ---------- ASV metadata ----------
asv_meta <- data.frame(
  ASV_sequence = colnames(read_red),
  length = nchar(colnames(read_red)),
  total_reads = colSums(read_red, na.rm = TRUE),
  prevalence_nSamples = colSums(read_red > 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)

rownames(asv_meta) <- paste0("ASV", marker_label, "_", seq_len(nrow(asv_meta)))

# ---------- Outputs ----------
saveRDS(asv_meta, file.path(path.expand("~/DADA2_outputs"), paste0("asv_meta", marker_label, ".rds")))
saveRDS(seqtab_nochim, file.path(path.expand("~/DADA2_outputs"), paste0("seqtab_nochim", marker_label, ".rds")))
saveRDS(read_ids, file.path(path.expand("~/DADA2_outputs"), paste0("seqtab_read_ids", marker_label, ".rds")))
saveRDS(read_red, file.path(path.expand("~/DADA2_outputs"), paste0("seqtab_read_red", marker_label, ".rds")))

write.csv(asv_meta, file.path(path.expand("~/DADA2_outputs"), paste0("asv_meta", marker_label, ".csv")), row.names = FALSE)
