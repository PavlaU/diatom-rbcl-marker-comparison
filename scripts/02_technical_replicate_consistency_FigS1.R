########################################
# TECHNICAL REPLICATE CONSISTENCY
# repeated rarefaction + binomial model
########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

# ---------- Paths and input ----------
project_dir <- normalizePath(".", mustWork = TRUE)
dada_dir <- file.path(project_dir, "DADA2_outputs")
out_dir <- file.path(project_dir, "tech_rep_consistency")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

nch263 <- readRDS(file.path(dada_dir, "seqtab_nochim263.rds"))
nch331 <- readRDS(file.path(dada_dir, "seqtab_nochim331.rds"))

target_depth <- min(min(rowSums(nch263)),min(rowSums(nch331)))
n_iter <- 100

# ---------- FUNCTION: get ASV-level reproducibility from one rarefied matrix ----------
get_repro_df <- function(mat, marker, iter_id) {
  base_id <- sub("_[12]$", "", rownames(mat))
  unique_ids <- unique(base_id)
  
  out_list <- vector("list", length(unique_ids))
  
  for (i in seq_along(unique_ids)) {
    id <- unique_ids[i]
    idx <- which(base_id == id)
    if (length(idx) != 2) next
    
    x1 <- as.numeric(mat[idx[1], ])
    x2 <- as.numeric(mat[idx[2], ])
    
    total <- x1 + x2
    present_in_both <- as.integer(x1 > 0 & x2 > 0)
    present_in_any  <- as.integer(x1 > 0 | x2 > 0)
    
    df <- data.frame(
      sample_id = id,
      ASV = colnames(mat),
      total_reads = total,
      present_in_both = present_in_both,
      present_in_any = present_in_any,
      marker = marker,
      iter = iter_id,
      stringsAsFactors = FALSE
    )
    
    df <- df[df$present_in_any == 1, , drop = FALSE]
    out_list[[i]] <- df
  }
  
  do.call(rbind, out_list)
}


# ---------- 100× rarefaction ----------
all_res <- vector("list", n_iter * 2)
k <- 1

for (i in seq_len(n_iter)) {
  set.seed(100 + i)
  
  rare_m1 <- rrarefy(nch263, sample = target_depth)
  rare_m2 <- rrarefy(nch331, sample = target_depth)
  
  rare_m1 <- rare_m1[rowSums(rare_m1) > 0, , drop = FALSE]
  rare_m2 <- rare_m2[rowSums(rare_m2) > 0, , drop = FALSE]
  
  all_res[[k]] <- get_repro_df(rare_m1, marker = "263", iter_id = i); k <- k + 1
  all_res[[k]] <- get_repro_df(rare_m2, marker = "331", iter_id = i); k <- k + 1
}

repro_df <- bind_rows(all_res)


# ---------- Binomial model ----------
glm_df <- repro_df %>%
  group_by(marker, sample_id, ASV) %>%
  summarise(
    mean_total_reads = mean(total_reads, na.rm = TRUE),
    prop_present_in_both = mean(present_in_both, na.rm = TRUE),
    n_iter = n(),
    .groups = "drop"
  ) %>%
  filter(mean_total_reads > 0)

# weighted binomial GLM
glm_fit <- glm(
  prop_present_in_both ~ log10(mean_total_reads) * marker,
  data = glm_df,
  weights = n_iter,
  family = binomial
)

print(summary(glm_fit))

# ---------- Prediction curves ----------
pred_df <- expand.grid(
  mean_total_reads = exp(seq(log(1), log(max(glm_df$mean_total_reads)), length.out = 200)),
  marker = c("263", "331")
)

pred_df$pred <- predict(glm_fit, newdata = pred_df, type = "response")


# ---------- Plots  ----------
# 1. observed by bins
repro_df <- repro_df %>%
  mutate(
    bin = cut(
      total_reads,
      breaks = c(0, 1, 2, 3, 5, 10, 20, 50, 100, Inf),
      labels = c("1", "2", "3", "4-5", "6-10", "11-20", "21-50", "51-100", "100+"),
      right = TRUE
    )
  )

bin_summary <- repro_df %>%
  group_by(marker, iter, bin) %>%
  summarise(
    prop_present_in_both = mean(present_in_both),
    n_ASV = n(),
    .groups = "drop"
  ) %>%
  group_by(marker, bin) %>%
  summarise(
    mean_prop_present_in_both = mean(prop_present_in_both, na.rm = TRUE),
    sd_prop_present_in_both   = sd(prop_present_in_both, na.rm = TRUE),
    mean_n_ASV               = mean(n_ASV, na.rm = TRUE),
    sd_n_ASV                 = sd(n_ASV, na.rm = TRUE),
    .groups = "drop"
  )

p_obs <- ggplot(bin_summary, aes(x = bin, y = mean_prop_present_in_both, color = marker, group = marker)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  geom_point(size = 2, alpha = 0.9) +
  theme_bw(base_size = 14) +
  labs(
    x = "Total reads across technical replicates",
    y = "Proportion of ASVs detected in both replicates",
    color = "Marker"
  )

print(p_obs)

# 2. model-based prediction (Fig S1)
p_pred <- ggplot() +
  geom_jitter(
    data = glm_df,
    aes(x = mean_total_reads, y = prop_present_in_both, color = marker),
    width = 0,
    height = 0.03,
    alpha = 0.08,
    size = 0.8
  ) +
  geom_line(
    data = pred_df,
    aes(x = mean_total_reads, y = pred, color = marker),
    linewidth = 1.2
  ) +
  scale_x_log10() +
  theme_bw(base_size = 14) +
  labs(
    x = "Mean total reads across technical replicates (log scale)",
    y = "Predicted probability of detection in both replicates",
    color = "Marker"
  )

print(p_pred)


# ---------- Outputs ---------- 
write.csv(repro_df, file.path(out_dir, "reproducibility_all_iterations.csv"), row.names = FALSE)
write.csv(bin_summary, file.path(out_dir, "reproducibility_bin_summary.csv"), row.names = FALSE)
write.csv(glm_df, file.path(out_dir, "reproducibility_glm_input.csv"), row.names = FALSE)

saveRDS(glm_fit, file.path(out_dir,"reproducibility_glm_fit.rds"))

ggsave(file.path(out_dir, "reproducibility_observed_bins.png"),p_obs,width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "FigS1_reproducibility_glm_prediction.png"), p_pred, width = 7, height = 5, dpi = 300)
