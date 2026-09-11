# HNSCC PDX transcription-factor activity and master-regulator analysis
# Adam McClatchey
#
# The full TF activity pipeline is shown for the representative average
# 3-6-week R vs NR analysis. Cross-classification sections use the saved
# outputs from the other five analyses, for which the same pipeline was
# previously repeated.

rm(list = ls(all.names = TRUE))
gc()

############################################################
# Packages
############################################################

library(DESeq2)
library(tidyverse)
library(dorothea)
library(viper)
library(limma)
library(ggrepel)
library(ComplexUpset)

############################################################
# Output folders
############################################################

root_output <- "Results_TF_Activity_Clean"
representative_output_dir <- file.path(
  root_output,
  "Average_3-6wk_RvsNR"
)
comparison_output_dir <- file.path(
  root_output,
  "Cross_Classification"
)

dir.create(representative_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_output_dir, recursive = TRUE, showWarnings = FALSE)

output <- paste0(representative_output_dir, "/")
prefix <- "DoRothEA_3-6wk_RvsNR"

message(
  "All outputs will be saved under: ",
  normalizePath(root_output, mustWork = FALSE)
)

############################################################
# REPRESENTATIVE AVERAGE 3-6-WEEK R VS NR TF PIPELINE
############################################################

############################################################
# Load raw counts and metadata
############################################################

counts <- read.csv(
  "Basal_RNA/PDX/All_RemovedSamples_Unfiltered_RawCounts_GraftHNSCC_270424.csv",
  check.names = FALSE
)

meta <- read.csv(
  "Basal_RNA/PDX/All_RemovedSamples_Unfiltered_GraftHNSCC_colData_240724.csv",
  check.names = FALSE
)

colnames(meta)[1] <- "SampleName"

rownames(counts) <- counts[, 1]
counts <- counts[, -1]

############################################################
# Select average 3-6-week R vs NR samples
############################################################

meta_avg <- meta[!is.na(meta$Definitive_Response_3_6wk), ]
counts_avg <- counts[, meta_avg$Sample_ID]

stopifnot(all(colnames(counts_avg) == meta_avg$Sample_ID))

# Remove the same low-depth outlier as in the original analysis.
outlier_sample <- "HNC0157PRX0A"

meta_avg <- meta_avg %>%
  filter(Sample_ID != outlier_sample)

counts_avg <- counts_avg[
  ,
  colnames(counts_avg) != outlier_sample
]

stopifnot(all(colnames(counts_avg) == meta_avg$Sample_ID))

# Explicit sample identifiers ensure correct count/metadata alignment.
rownames(meta_avg) <- meta_avg$Sample_ID

############################################################
# DESeq2 object and variance-stabilising transformation
############################################################

meta_avg$Definitive_Response_3_6wk <- factor(
  meta_avg$Definitive_Response_3_6wk,
  levels = c("NonResponder", "Responder")
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_avg,
  colData = meta_avg,
  design = ~ Definitive_Response_3_6wk
)

dds <- dds[rowSums(counts(dds)) >= 10, ]

vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

############################################################
# DoRothEA A-C regulons and VIPER TF activity
############################################################

data(dorothea_hs, package = "dorothea")

regulons <- dorothea_hs %>%
  filter(confidence %in% c("A", "B", "C"))

regulon <- df2regulon(regulons)

tf_activity <- viper(
  expr_mat,
  regulon,
  verbose = FALSE
)

write.csv(
  tf_activity,
  paste0(output, prefix, "_TF_activity_matrix.csv")
)

write.csv(
  regulons,
  paste0(output, prefix, "_DoRothEA_ABC_regulons.csv"),
  row.names = FALSE
)

############################################################
# Differential TF activity using limma
############################################################

# Positive activity_difference indicates greater inferred activity in
# responders; a negative value indicates greater activity in non-responders.
group <- factor(
  meta_avg$Definitive_Response_3_6wk,
  levels = c("NonResponder", "Responder")
)

design <- model.matrix(~ group)

fit <- lmFit(tf_activity, design)
fit <- eBayes(fit)

tf_results <- topTable(
  fit,
  coef = "groupResponder",
  number = Inf,
  sort.by = "P"
)

tf_results$TF <- rownames(tf_results)
tf_results <- tf_results %>% relocate(TF)

colnames(tf_results)[
  colnames(tf_results) == "logFC"
] <- "activity_difference"

tf_results$Significant <- tf_results$adj.P.Val < 0.05

write.csv(
  tf_results,
  paste0(output, prefix, "_TF_activity_results.csv"),
  row.names = FALSE
)

tf_results_sig_ranked <- tf_results %>%
  arrange(adj.P.Val)

tf_results_up <- tf_results %>%
  filter(Significant) %>%
  arrange(desc(activity_difference))

tf_results_down <- tf_results %>%
  filter(Significant) %>%
  arrange(activity_difference)

write.csv(
  tf_results_sig_ranked,
  paste0(output, prefix, "_TF_ranked_by_significance.csv"),
  row.names = FALSE
)

write.csv(
  tf_results_up,
  paste0(output, prefix, "_TF_highest_in_responders.csv"),
  row.names = FALSE
)

write.csv(
  tf_results_down,
  paste0(output, prefix, "_TF_highest_in_nonresponders.csv"),
  row.names = FALSE
)

############################################################
# Map significant TF regulons to significant DEGs
############################################################

sig_tfs <- tf_results %>%
  filter(Significant)

degs <- read.csv(
  "Results_Avg_RvsNR/3-6wkAvg_RvsNR_DEGs.csv"
)

mapped_interactions <- regulons %>%
  filter(
    tf %in% sig_tfs$TF,
    target %in% degs$Gene
  )

write.csv(
  mapped_interactions,
  paste0(output, prefix, "_TF_DEG_mapping.csv"),
  row.names = FALSE
)

tf_target_counts <- mapped_interactions %>%
  group_by(tf) %>%
  summarise(
    DEG_targets = n_distinct(target),
    .groups = "drop"
  ) %>%
  arrange(desc(DEG_targets))

write.csv(
  tf_target_counts,
  paste0(output, prefix, "_TF_target_counts.csv"),
  row.names = FALSE
)

gene_tf_counts <- mapped_interactions %>%
  group_by(target) %>%
  summarise(
    TF_count = n_distinct(tf),
    .groups = "drop"
  ) %>%
  arrange(desc(TF_count))

write.csv(
  gene_tf_counts,
  paste0(output, prefix, "_Shared_Target_Genes.csv"),
  row.names = FALSE
)

tf_summary <- tf_results %>%
  dplyr::select(
    TF,
    activity_difference,
    adj.P.Val
  ) %>%
  left_join(
    tf_target_counts,
    by = c("TF" = "tf")
  ) %>%
  arrange(desc(DEG_targets))

write.csv(
  tf_summary,
  paste0(output, prefix, "_TF_summary_table.csv"),
  row.names = FALSE
)

############################################################
# Figure 15 - candidate master regulators
############################################################

# This is the original Figure 15 selection and plotting code.
top_tfs_targets <- tf_target_counts %>%
  slice_max(
    order_by = DEG_targets,
    n = 15
  )

write.csv(
  top_tfs_targets,
  file.path(representative_output_dir, "Figure15_Master_Regulator_Data.csv"),
  row.names = FALSE
)

top_tfs_targets$tf <- factor(
  top_tfs_targets$tf,
  levels = rev(top_tfs_targets$tf)
)

master_regulator_plot <- ggplot(
  top_tfs_targets,
  aes(x = tf, y = DEG_targets)
) +
  geom_col(
    fill = "firebrick",
    width = 0.8
  ) +
  coord_flip() +
  labs(
    title = "Candidate Master Regulators",
    subtitle = "Number of mapped DEG targets per transcription factor",
    x = NULL,
    y = "Number of DEG Targets"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(
    representative_output_dir,
    "Figure15_Candidate_Master_Regulators.png"
  ),
  master_regulator_plot,
  width = 10,
  height = 8,
  dpi = 300
)

############################################################
# CROSS-CLASSIFICATION TF RESULTS
############################################################

analysis_order <- c(
  "3wk R vs NR",
  "6wk R vs NR",
  "Avg R vs NR",
  "3wk PR vs PD",
  "6wk PR vs PD",
  "Avg PR vs PD"
)

# The other five analyses are read-only inputs. The newly generated average
# R-vs-NR tables above are used for the representative comparison.
tf_3wk_rnr <- read.csv(
  "Results_DoRothEA_3wk_RvsNR/DoRothEA_3wk_RvsNR_TF_activity_results.csv"
)
tf_6wk_rnr <- read.csv(
  "Results_DoRothEA_6wk_RvsNR/DoRothEA_6wk_RvsNR_TF_activity_results.csv"
)
tf_avg_rnr <- tf_results
tf_3wk_prpd <- read.csv(
  "Results_DoRothEA_3wk_PRvsPD/DoRothEA_3wk_PRvsPD_TF_activity_results.csv"
)
tf_6wk_prpd <- read.csv(
  "Results_DoRothEA_6wk_PRvsPD/DoRothEA_6wk_PRvsPD_TF_activity_results.csv"
)
tf_avg_prpd <- read.csv(
  "Results_DoRothEA_3-6wkAvg_PRvsPD/DoRothEA_3-6wkAvg_PRvsPD_TF_activity_results.csv"
)

tf_result_tables <- list(
  "3wk R vs NR" = tf_3wk_rnr,
  "6wk R vs NR" = tf_6wk_rnr,
  "Avg R vs NR" = tf_avg_rnr,
  "3wk PR vs PD" = tf_3wk_prpd,
  "6wk PR vs PD" = tf_6wk_prpd,
  "Avg PR vs PD" = tf_avg_prpd
)

############################################################
# Significant TF counts, recurrence and supporting tables
############################################################

tf_lists <- lapply(
  tf_result_tables,
  function(x) {
    x %>%
      filter(Significant == TRUE) %>%
      pull(TF) %>%
      unique()
  }
)

significant_counts <- tibble(
  Analysis = names(tf_lists),
  Significant_TFs = lengths(tf_lists)
)

all_unique_tfs <- unique(unlist(tf_lists))

tf_frequency <- tibble(TF = all_unique_tfs) %>%
  rowwise() %>%
  mutate(
    Analyses_Present = sum(
      vapply(tf_lists, function(x) TF %in% x, logical(1))
    )
  ) %>%
  ungroup() %>%
  arrange(desc(Analyses_Present), TF)

frequency_summary <- tf_frequency %>%
  dplyr::count(
    Analyses_Present,
    name = "Number_of_TFs"
  ) %>%
  tidyr::complete(
    Analyses_Present = 1:6,
    fill = list(Number_of_TFs = 0)
  ) %>%
  arrange(desc(Analyses_Present))

results_3_4_summary <- tibble(
  Metric = c(
    "Unique significant TFs",
    "Significant in at least 2 analyses",
    "Significant in at least 3 analyses",
    "Significant in at least 4 analyses",
    "Significant in 5 analyses"
  ),
  Value = c(
    nrow(tf_frequency),
    sum(tf_frequency$Analyses_Present >= 2),
    sum(tf_frequency$Analyses_Present >= 3),
    sum(tf_frequency$Analyses_Present >= 4),
    sum(tf_frequency$Analyses_Present == 5)
  )
)

recurrent_tfs <- tf_frequency %>%
  filter(Analyses_Present >= 3)

highly_recurrent_tfs <- tf_frequency %>%
  filter(Analyses_Present >= 4)

unique_by_analysis <- map_dfr(
  names(tf_lists),
  function(analysis) {
    other_tfs <- unique(
      unlist(tf_lists[names(tf_lists) != analysis])
    )

    tibble(
      Analysis = analysis,
      TF = setdiff(tf_lists[[analysis]], other_tfs)
    )
  }
)

unique_counts <- unique_by_analysis %>%
  dplyr::count(
    Analysis,
    name = "Unique_TFs"
  ) %>%
  right_join(
    tibble(Analysis = analysis_order),
    by = "Analysis"
  ) %>%
  mutate(Unique_TFs = replace_na(Unique_TFs, 0L)) %>%
  arrange(match(Analysis, analysis_order))

tf_all <- bind_rows(
  lapply(names(tf_result_tables), function(x) {
    tf_result_tables[[x]] %>% mutate(Analysis = x)
  })
)

tf_significant_all <- tf_all %>%
  filter(Significant == TRUE)

tf_recurrence_table <- tf_significant_all %>%
  group_by(TF) %>%
  summarise(
    Analyses_Present = n_distinct(Analysis),
    Mean_Activity_Difference = mean(activity_difference, na.rm = TRUE),
    Mean_Absolute_Activity_Difference = mean(
      abs(activity_difference),
      na.rm = TRUE
    ),
    Max_Absolute_Activity_Difference = max(
      abs(activity_difference),
      na.rm = TRUE
    ),
    Analyses = paste(sort(unique(Analysis)), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(
    desc(Analyses_Present),
    desc(Mean_Absolute_Activity_Difference)
  )

tf_recurrence_top <- tf_recurrence_table %>%
  filter(Analyses_Present >= 4)

write.csv(
  significant_counts,
  file.path(comparison_output_dir, "Significant_TF_Counts.csv"),
  row.names = FALSE
)
write.csv(
  tf_frequency,
  file.path(comparison_output_dir, "TF_Frequency_Across_Analyses.csv"),
  row.names = FALSE
)
write.csv(
  frequency_summary,
  file.path(comparison_output_dir, "TF_Recurrence_Counts.csv"),
  row.names = FALSE
)
write.csv(
  results_3_4_summary,
  file.path(comparison_output_dir, "Results3.4_TF_Summary.csv"),
  row.names = FALSE
)
write.csv(
  recurrent_tfs,
  file.path(comparison_output_dir, "TFs_Recurrent_3plus.csv"),
  row.names = FALSE
)
write.csv(
  highly_recurrent_tfs,
  file.path(comparison_output_dir, "TFs_Recurrent_4plus.csv"),
  row.names = FALSE
)
write.csv(
  unique_by_analysis,
  file.path(comparison_output_dir, "TFs_Unique_To_Analyses.csv"),
  row.names = FALSE
)
write.csv(
  unique_counts,
  file.path(comparison_output_dir, "Unique_TF_Counts.csv"),
  row.names = FALSE
)
write.csv(
  tf_significant_all,
  file.path(comparison_output_dir, "Significant_TF_Results_All_Analyses.csv"),
  row.names = FALSE
)
write.csv(
  tf_recurrence_table,
  file.path(comparison_output_dir, "TF_Recurrence_Table_All.csv"),
  row.names = FALSE
)
write.csv(
  tf_recurrence_top,
  file.path(comparison_output_dir, "TF_Recurrence_Table_4plus.csv"),
  row.names = FALSE
)

############################################################
# Load TF-DEG summary and mapping tables for Tables 5-6
############################################################

tf_3wk_rnr_summary <- read.csv(
  "Results_DoRothEA_3wk_RvsNR/DoRothEA_3wk_RvsNR_TF_summary_table.csv"
)
tf_6wk_rnr_summary <- read.csv(
  "Results_DoRothEA_6wk_RvsNR/DoRothEA_6wk_RvsNR_TF_summary_table.csv"
)
tf_avg_rnr_summary <- tf_summary
tf_3wk_prpd_summary <- read.csv(
  "Results_DoRothEA_3wk_PRvsPD/DoRothEA_3wk_PRvsPD_TF_summary_table.csv"
)
tf_6wk_prpd_summary <- read.csv(
  "Results_DoRothEA_6wk_PRvsPD/DoRothEA_6wk_PRvsPD_TF_summary_table.csv"
)
tf_avg_prpd_summary <- read.csv(
  "Results_DoRothEA_3-6wkAvg_PRvsPD/DoRothEA_3-6wkAvg_PRvsPD_TF_summary_table.csv"
)

tf_summary_tables <- list(
  "3wk R vs NR" = tf_3wk_rnr_summary,
  "6wk R vs NR" = tf_6wk_rnr_summary,
  "Avg R vs NR" = tf_avg_rnr_summary,
  "3wk PR vs PD" = tf_3wk_prpd_summary,
  "6wk PR vs PD" = tf_6wk_prpd_summary,
  "Avg PR vs PD" = tf_avg_prpd_summary
)

mapping_3wk_rnr <- read.csv(
  "Results_DoRothEA_3wk_RvsNR/DoRothEA_3wk_RvsNR_TF_DEG_mapping.csv"
)
mapping_6wk_rnr <- read.csv(
  "Results_DoRothEA_6wk_RvsNR/DoRothEA_6wk_RvsNR_TF_DEG_mapping.csv"
)
mapping_avg_rnr <- mapped_interactions
mapping_3wk_prpd <- read.csv(
  "Results_DoRothEA_3wk_PRvsPD/DoRothEA_3wk_PRvsPD_TF_DEG_mapping.csv"
)
mapping_6wk_prpd <- read.csv(
  "Results_DoRothEA_6wk_PRvsPD/DoRothEA_6wk_PRvsPD_TF_DEG_mapping.csv"
)
mapping_avg_prpd <- read.csv(
  "Results_DoRothEA_3-6wkAvg_PRvsPD/DoRothEA_3-6wkAvg_PRvsPD_TF_DEG_mapping.csv"
)

mapping_tables <- list(
  "3wk R vs NR" = mapping_3wk_rnr,
  "6wk R vs NR" = mapping_6wk_rnr,
  "Avg R vs NR" = mapping_avg_rnr,
  "3wk PR vs PD" = mapping_3wk_prpd,
  "6wk PR vs PD" = mapping_6wk_prpd,
  "Avg PR vs PD" = mapping_avg_prpd
)

############################################################
# Table 5 - significant TF and TF-DEG mapping summary
############################################################

table5 <- tibble(
  Comparison = analysis_order,
  Significant_TFs = lengths(tf_lists[analysis_order]),
  Mapped_DEGs = vapply(
    mapping_tables[analysis_order],
    function(x) length(unique(x$target)),
    integer(1)
  ),
  TF_DEG_Interactions = vapply(
    mapping_tables[analysis_order],
    nrow,
    integer(1)
  )
)

write.csv(
  table5,
  file.path(comparison_output_dir, "Table5_TF_and_DEG_Mapping_Summary.csv"),
  row.names = FALSE
)

############################################################
# Figure 13 - TF overlap UpSet plot
############################################################

tf_summary_all <- bind_rows(
  lapply(names(tf_summary_tables), function(x) {
    tf_summary_tables[[x]] %>% mutate(Comparison = x)
  })
)

tf_summary_sig <- tf_summary_all %>%
  filter(adj.P.Val < 0.05)

tf_upset <- tf_summary_sig %>%
  dplyr::select(TF, Comparison) %>%
  distinct() %>%
  mutate(Present = TRUE) %>%
  pivot_wider(
    names_from = Comparison,
    values_from = Present,
    values_fill = FALSE
  )

tf_intersection_summary <- tf_upset %>%
  group_by(across(-TF)) %>%
  summarise(Intersection_Size = n(), .groups = "drop") %>%
  arrange(desc(Intersection_Size))

write.csv(
  tf_intersection_summary,
  file.path(comparison_output_dir, "Figure13_TF_Intersection_Summary.csv"),
  row.names = FALSE
)

tf_upset_plot <- upset(
  tf_upset,
  intersect = colnames(tf_upset)[-1],
  name = "Comparison",
  width_ratio = 0.2,
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  set_sizes = upset_set_size()
) +
  ggtitle("Overlap of Significant Transcription Factors Across Comparisons")

ggsave(
  file.path(comparison_output_dir, "Figure13_TF_UpSet_plot.png"),
  tf_upset_plot,
  width = 10,
  height = 6,
  dpi = 300
)

############################################################
# Figure 14 - top-20 differential TF activity heatmap
############################################################

analysis_labels <- c(
  "3wk R vs NR" = "3wk",
  "6wk R vs NR" = "6wk",
  "Avg R vs NR" = "Avg",
  "3wk PR vs PD" = "3wk",
  "6wk PR vs PD" = "6wk",
  "Avg PR vs PD" = "Avg"
)

selected_tfs <- tf_recurrence_table %>%
  arrange(desc(Mean_Absolute_Activity_Difference)) %>%
  slice_head(n = 20)

write.csv(
  selected_tfs,
  file.path(comparison_output_dir, "Figure14_Top20_TF_Data.csv"),
  row.names = FALSE
)

tf_order <- selected_tfs %>%
  arrange(
    desc(Analyses_Present),
    desc(Mean_Absolute_Activity_Difference)
  ) %>%
  pull(TF)

tf_heatmap_values <- tf_all %>%
  filter(
    TF %in% selected_tfs$TF,
    Significant == TRUE
  ) %>%
  dplyr::select(
    TF,
    Analysis,
    activity_difference
  )

tf_heatmap_df <- tf_heatmap_values %>%
  right_join(
    expand_grid(
      TF = selected_tfs$TF,
      Analysis = analysis_order
    ),
    by = c("TF", "Analysis")
  ) %>%
  mutate(
    Analysis = factor(Analysis, levels = analysis_order),
    Analysis_Label = recode(as.character(Analysis), !!!analysis_labels),
    Analysis_Label = factor(Analysis_Label, levels = c("3wk", "6wk", "Avg")),
    Comparison = case_when(
      Analysis %in% c(
        "3wk R vs NR",
        "6wk R vs NR",
        "Avg R vs NR"
      ) ~ "R vs NR",
      TRUE ~ "PR vs PD"
    ),
    Comparison = factor(Comparison, levels = c("R vs NR", "PR vs PD")),
    TF = factor(TF, levels = rev(tf_order))
  )

tf_heatmap_final <- ggplot(
  tf_heatmap_df,
  aes(
    x = Analysis_Label,
    y = TF,
    fill = activity_difference
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.55,
    width = 0.96,
    height = 0.92
  ) +
  facet_grid(
    cols = vars(Comparison),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_distiller(
    name = "Activity\ndifference",
    palette = "RdBu",
    direction = -1,
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    breaks = c(-1.5, -0.75, 0, 0.75, 1.5),
    na.value = "grey95"
  ) +
  labs(
    title = "Differential TF Activity Across Response Classifications",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "grey94", colour = "grey70"),
    panel.spacing.x = grid::unit(0.7, "cm"),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.4),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(4.5, "cm"),
      barwidth = grid::unit(0.5, "cm")
    )
  )

ggsave(
  file.path(comparison_output_dir, "Figure14_TF_Activity_Heatmap_Top20.png"),
  tf_heatmap_final,
  width = 8.8,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(comparison_output_dir, "Figure14_TF_Activity_Heatmap_Top20.pdf"),
  tf_heatmap_final,
  width = 8.8,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

############################################################
# Table 6 - leading candidate master regulators
############################################################

table6_detail <- tf_summary_sig %>%
  filter(!is.na(DEG_targets)) %>%
  group_by(Comparison) %>%
  slice_max(
    order_by = DEG_targets,
    n = 1,
    with_ties = TRUE
  ) %>%
  arrange(TF, .by_group = TRUE) %>%
  ungroup()

write.csv(
  table6_detail,
  file.path(comparison_output_dir, "Table6_Leading_Master_Regulators_Detail.csv"),
  row.names = FALSE
)

table6 <- table6_detail %>%
  group_by(Comparison) %>%
  summarise(
    Leading_TFs = paste(TF, collapse = " / "),
    DEG_Targets = if_else(
      n() > 1,
      paste0(dplyr::first(DEG_targets), " each"),
      as.character(dplyr::first(DEG_targets))
    ),
    Activity_Difference = paste(
      sprintf("%+.2f", activity_difference),
      collapse = " / "
    ),
    .groups = "drop"
  ) %>%
  right_join(
    tibble(Comparison = analysis_order),
    by = "Comparison"
  ) %>%
  mutate(
    Leading_TFs = replace_na(Leading_TFs, "-"),
    DEG_Targets = replace_na(DEG_Targets, "-"),
    Activity_Difference = replace_na(Activity_Difference, "-")
  ) %>%
  arrange(match(Comparison, analysis_order))

write.csv(
  table6,
  file.path(comparison_output_dir, "Table6_Leading_Master_Regulators.csv"),
  row.names = FALSE
)

############################################################
# Additional master-regulator comparison tables
############################################################

tf_repeat <- tf_summary_sig %>%
  group_by(TF) %>%
  summarise(
    Analyses_present = n(),
    Comparisons = paste(Comparison, collapse = ", "),
    Max_DEG_targets = if (
      all(is.na(DEG_targets))
    ) {
      NA_integer_
    } else {
      max(DEG_targets, na.rm = TRUE)
    },
    Max_activity_difference = max(abs(activity_difference)),
    Best_padj = min(adj.P.Val),
    .groups = "drop"
  ) %>%
  arrange(
    desc(Analyses_present),
    desc(Max_DEG_targets),
    Best_padj
  )

top_recurrent_tfs <- tf_repeat %>%
  slice_head(n = 10)

top_deg_target_tfs <- tf_summary_sig %>%
  filter(!is.na(DEG_targets)) %>%
  arrange(desc(DEG_targets)) %>%
  dplyr::select(
    Comparison,
    TF,
    activity_difference,
    adj.P.Val,
    DEG_targets
  ) %>%
  slice_head(n = 10)

top_activity_tfs <- tf_summary_sig %>%
  arrange(desc(abs(activity_difference))) %>%
  dplyr::select(
    Comparison,
    TF,
    activity_difference,
    adj.P.Val,
    DEG_targets
  ) %>%
  slice_head(n = 10)

write.csv(
  tf_repeat,
  file.path(comparison_output_dir, "Master_Regulator_Recurrence_All.csv"),
  row.names = FALSE
)
write.csv(
  top_recurrent_tfs,
  file.path(comparison_output_dir, "Top_Recurrent_TFs.csv"),
  row.names = FALSE
)
write.csv(
  top_deg_target_tfs,
  file.path(comparison_output_dir, "Top_DEG_Target_TFs.csv"),
  row.names = FALSE
)
write.csv(
  top_activity_tfs,
  file.path(comparison_output_dir, "Top_Activity_Difference_TFs.csv"),
  row.names = FALSE
)

############################################################
# Reproducibility and output manifest
############################################################

writeLines(
  capture.output(sessionInfo()),
  file.path(root_output, "Session_Info.txt")
)

generated_files <- list.files(
  root_output,
  recursive = TRUE,
  full.names = TRUE
)

output_manifest <- tibble(
  File = sub(paste0("^", root_output, "/"), "", generated_files),
  Size_bytes = file.info(generated_files)$size
)

write.csv(
  output_manifest,
  file.path(root_output, "Output_Manifest.csv"),
  row.names = FALSE
)

message("Finished. All newly written results are in: ", root_output)
