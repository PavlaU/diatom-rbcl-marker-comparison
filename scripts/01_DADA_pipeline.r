############################################################
# DADA_pipeline
# Chimera-free ASV generation with merged technical extractions
#
# Required inputs:
# - NGS_raw/
#   Folder with FASTQ files; raw reads are deposited in ENA under PRJEB106074.
#
# - input/ids263.fasta and input/ids331.fasta
#   Sequences of internal DNA spike-ins; included in this GitHub repository.
############################################################

library(Biostrings)
library(dada2)
library(stringr)

set.seed(1234)

# ---------- Paths ----------
project_dir <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(project_dir, "NGS_raw")
input_dir <- file.path(project_dir, "input")

# ---------- Output directory ----------
out_dir <- file.path(project_dir, "DADA2_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- Cutadapt ----------
cutadapt_bin <- "/path/to/cutadapt-env/bin/cutadapt"
if (!file.exists(cutadapt_bin)) stop("cutadapt not found: ", cutadapt_bin)

# ---------- Marker selection ----------
# Set TRUE for 263-bp, FALSE for 331-bp
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

fwd_primer_rc <- as.character(reverseComplement(DNAStringSet(fwd_primer)))
rev_primer_rc <- as.character(reverseComplement(DNAStringSet(rev_primer)))

marker_label <- if (marker_263) "263" else "331"
message("Running marker: ", marker_label)

# ---------- Input FASTQ ----------
fwd_files <- sort(list.files(raw_dir, pattern = "_1.fastq.gz$", full.names = TRUE))
rev_files <- sort(list.files(raw_dir, pattern = "_2.fastq.gz$", full.names = TRUE))

sample_names_raw <- str_extract(basename(fwd_files), "^[^_]+_[^_]+_[^_]+")

if (marker_263) {
  # 263-bp marker: samples without "J" in raw sample ID
  keep_idx <- !grepl("J", sample_names_raw)
} else {
  # 331-bp marker: samples with "J" in raw sample ID
  keep_idx <- grepl("J", sample_names_raw)
}

fwd_files <- fwd_files[keep_idx]
rev_files <- rev_files[keep_idx]
sample_names_raw <- sample_names_raw[keep_idx]

sample_names <- sample_names_raw
sample_names <- sub("^.*_([0-9]+J?_[12])$", "\\1", sample_names)
sample_names <- sub("J(?=_[12]$)", "", sample_names, perl = TRUE)

names(fwd_files) <- sample_names
names(rev_files) <- sample_names

# ---------- Primer trimming ----------
trim_dir <- file.path(out_dir, paste0("01_trimmed_", marker_label))
dir.create(trim_dir, showWarnings = FALSE, recursive = TRUE)

# Save trimmed files under standardised sample names.
fwd_trimmed <- file.path(trim_dir, paste0(sample_names, "_1.fastq.gz"))
rev_trimmed <- file.path(trim_dir, paste0(sample_names, "_2.fastq.gz"))

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

for (i in seq_along(sample_names)) {
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

# ---------- Split by run ----------
split_dir <- file.path(out_dir, paste0("02_split_by_run_", marker_label))
dir.create(split_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(trim_dir, pattern = "\\.fastq.gz$", full.names = TRUE)

get_run_name <- function(header) {
  header <- sub("^@", "", header)
  header <- strsplit(header, " ", fixed = TRUE)[[1]][1]
  parts <- strsplit(header, ":", fixed = TRUE)[[1]]
  paste(parts[1], parts[2], parts[3], sep = "_")
}

closeAllConnections()

for (file in files) {
  cat("Processing:", basename(file), "\n")
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
      
      out_file <- file.path(run_folder, basename(file))
      out_connections[[run_name]] <- gzfile(out_file, "wt")
    }
    
    writeLines(lines, out_connections[[run_name]])
  }
  
  close(con_in)
  
  for (x in out_connections) {
    close(x)
  }
  
  cat("Finished:", basename(file), "\n")
}

# ---------- Filtering + DADA2 by run ----------
getN <- function(x) sum(getUniques(x))

run_dirs <- list.dirs(split_dir, full.names = TRUE, recursive = FALSE)

seqtabs_by_run <- list()
filter_stats_by_run <- list()

filt_base_dir <- file.path(out_dir, paste0("03_filtered_by_run_", marker_label))
dir.create(filt_base_dir, showWarnings = FALSE, recursive = TRUE)

for (run_dir in run_dirs) {
  run_name <- basename(run_dir)
  
  run_fwd <- sort(list.files(run_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
  run_rev <- sort(list.files(run_dir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))

  run_samples <- sub("_1\\.fastq\\.gz$", "", basename(run_fwd))
  

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

# ---------- Keep red-lineage target-length ASVs with >=0.01% relative abundance ----------
read_red <- seqtab_final[, nchar(colnames(seqtab_final)) == red_len, drop = FALSE]
abd_red <- read_red / rowSums(seqtab_final)

keep_asv <- apply(abd_red, 2, function(x) any(x >= 0.0001))

abd_red_filt <- abd_red[, keep_asv, drop = FALSE]
read_red_filt <- read_red[, keep_asv, drop = FALSE]

# ---------- IDS / biological split ----------
if (marker_263) {
  ids <- as.character(readDNAStringSet(file.path(input_dir, "ids263.fasta")))
} else {
  ids <- as.character(readDNAStringSet(file.path(input_dir, "ids331.fasta")))
}

max_ids_mismatch <- 5

is_ids_match <- function(asv, ids_ref, max_mismatch = 2) {
  hits <- ids_ref[nchar(ids_ref) == nchar(asv)]
  
  if (length(hits) == 0) {
    return(FALSE)
  }
  
  asv_chars <- strsplit(asv, "", fixed = TRUE)[[1]]
  
  any(vapply(hits, function(ref) {
    ref_chars <- strsplit(ref, "", fixed = TRUE)[[1]]
    sum(asv_chars != ref_chars) <= max_mismatch
  }, logical(1)))
}

ids_cols <- vapply(
  colnames(abd_red_filt),
  is_ids_match,
  logical(1),
  ids_ref = ids,
  max_mismatch = max_ids_mismatch
)

read_ids <- read_red_filt[, ids_cols, drop = FALSE]
read_red <- read_red_filt[, !ids_cols, drop = FALSE]

# ---------- ASV metadata ----------
asv_sequences <- colnames(read_red)
asv_total_reads <- colSums(read_red, na.rm = TRUE)
asv_length <- nchar(asv_sequences)
asv_prevalence <- colSums(read_red > 0, na.rm = TRUE)

asv_meta <- data.frame(
  ASV_sequence = asv_sequences,
  length = asv_length,
  total_reads = asv_total_reads,
  prevalence_nSamples = asv_prevalence,
  stringsAsFactors = FALSE
)

rownames(asv_meta) <- paste0("ASV", marker_label, "_", seq_len(nrow(asv_meta)))

# ---------- Outputs ----------
saveRDS(asv_meta, file.path(out_dir, paste0("asv_meta", marker_label, ".rds")))
saveRDS(seqtab_nochim, file.path(out_dir, paste0("seqtab_nochim", marker_label, ".rds")))
saveRDS(read_ids, file.path(out_dir, paste0("seqtab_read_ids", marker_label, ".rds")))
saveRDS(read_red, file.path(out_dir, paste0("seqtab_read_red", marker_label, ".rds")))

write.csv(
  asv_meta,
  file.path(out_dir, paste0("asv_meta", marker_label, ".csv")),
  row.names = FALSE
)