# HNSCC PDX cohort, response trajectory and baseline QC
# Adam McClatchey

rm(list = ls(all.names = TRUE))
gc()

############################################################
# Packages
############################################################

library(DESeq2)
library(tidyverse)
library(ggalluvial)
library(ggrepel)
library(scales)
library(gtsummary)
library(gt)
library(flextable)
library(officer)

############################################################
# Paths
############################################################

output <- "./Results_Cohort_Response_Clean/"
prefix <- "3-6wkAvg_RvsNR"

dir.create(output, showWarnings = FALSE, recursive = TRUE)

message(
  "All outputs will be saved to: ",
  normalizePath(output, mustWork = FALSE)
)

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

stopifnot(all(meta$Sample_ID %in% colnames(counts)))

outlier_sample <- "HNC0157PRX0A"

############################################################
# TABLE 1 - SAMPLE COMPOSITION
############################################################

# The outlier was excluded from the six baseline analyses.
meta_analysis <- meta %>%
  filter(Sample_ID != outlier_sample)

sample_composition <- tibble(
  Comparison = c(
    "3wk R vs NR",
    "6wk R vs NR",
    "Avg 3-6 R vs NR",
    "3wk PR vs PD",
    "6wk PR vs PD",
    "Avg 3-6 PR vs PD"
  ),
  Group_1 = c("R", "R", "R", "PR", "PR", "PR"),
  Group_1_n = c(
    sum(meta_analysis$Definitive_Response_3wk == "Responder", na.rm = TRUE),
    sum(meta_analysis$Definitive_Response_6wk == "Responder", na.rm = TRUE),
    sum(meta_analysis$Definitive_Response_3_6wk == "Responder", na.rm = TRUE),
    sum(meta_analysis$Response_3wk == "PR", na.rm = TRUE),
    sum(meta_analysis$Response_6wk == "PR", na.rm = TRUE),
    sum(meta_analysis$Response_3_6wk == "PR", na.rm = TRUE)
  ),
  Group_2 = c("NR", "NR", "NR", "PD", "PD", "PD"),
  Group_2_n = c(
    sum(meta_analysis$Definitive_Response_3wk == "NonResponder", na.rm = TRUE),
    sum(meta_analysis$Definitive_Response_6wk == "NonResponder", na.rm = TRUE),
    sum(meta_analysis$Definitive_Response_3_6wk == "NonResponder", na.rm = TRUE),
    sum(meta_analysis$Response_3wk == "PD", na.rm = TRUE),
    sum(meta_analysis$Response_6wk == "PD", na.rm = TRUE),
    sum(meta_analysis$Response_3_6wk == "PD", na.rm = TRUE)
  )
) %>%
  mutate(Total_Samples_Included = Group_1_n + Group_2_n)

# Exact values reported in thesis Table 1.
stopifnot(
  identical(sample_composition$Group_1_n, c(39L, 34L, 40L, 4L, 7L, 4L)),
  identical(sample_composition$Group_2_n, c(15L, 9L, 14L, 15L, 9L, 14L)),
  identical(
    sample_composition$Total_Samples_Included,
    c(54L, 43L, 54L, 19L, 16L, 18L)
  )
)

write.csv(
  sample_composition,
  paste0(output, "Table1_Sample_Composition.csv"),
  row.names = FALSE
)

############################################################
# FIGURE 1 - RESPONSE CLASSIFICATION AND TRAJECTORY
############################################################

# The outlier is intentionally retained in this section.

missing_6wk <- meta %>%
  filter(!is.na(Response_3wk)) %>%
  mutate(
    Missing_6wk_response = is.na(Response_6wk)
  ) %>%
  dplyr::select(
    Sample_ID,
    Response_3wk,
    Response_6wk,
    Definitive_Response_3wk,
    Definitive_Response_6wk,
    percentage_vol_change_3wk,
    percentage_vol_change_6wk,
    Missing_6wk_response
  )

write.csv(
  missing_6wk,
  paste0(output, "Missing_6wk_Response_Table.csv"),
  row.names = FALSE
)

missing_summary <- missing_6wk %>%
  dplyr::count(Response_3wk, Missing_6wk_response)

write.csv(
  missing_summary,
  paste0(output, "Missing_6wk_By_3wk_Response.csv"),
  row.names = FALSE
)

# This object is retained from the original response-dynamics script.
response_dynamic <- meta %>%
  mutate(
    GrowShrink_3wk = case_when(
      is.na(percentage_vol_change_3wk) ~ NA_character_,
      percentage_vol_change_3wk < 0 ~ "Shrink",
      percentage_vol_change_3wk >= 0 ~ "Grow"
    ),
    GrowShrink_6wk = case_when(
      is.na(percentage_vol_change_6wk) ~ NA_character_,
      percentage_vol_change_6wk < 0 ~ "Shrink",
      percentage_vol_change_6wk >= 0 ~ "Grow"
    )
  )

write.csv(
  response_dynamic,
  paste0(output, "Response_Dynamics_Metadata.csv"),
  row.names = FALSE
)

############################################################
# Figure 1A - response transition Sankey plot
############################################################

sankey_response <- response_dynamic %>%
  filter(!is.na(Response_3wk)) %>%
  mutate(
    Response_6wk_plot = ifelse(
      is.na(Response_6wk),
      "Missing",
      Response_6wk
    )
  ) %>%
  dplyr::count(Response_3wk, Response_6wk_plot)

write.csv(
  sankey_response,
  paste0(output, "Response_Transition_Counts.csv"),
  row.names = FALSE
)

p_sankey_response <- ggplot(
  sankey_response,
  aes(
    axis1 = Response_3wk,
    axis2 = Response_6wk_plot,
    y = n
  )
) +
  geom_alluvium(
    aes(fill = Response_6wk_plot),
    width = 0.15,
    alpha = 0.8
  ) +
  geom_stratum(
    width = 0.15,
    fill = "white",
    colour = "black"
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 4
  ) +
  scale_x_discrete(
    limits = c("3-week response", "6-week response"),
    expand = c(0.1, 0.1)
  ) +
  labs(
    title = "Transition of response categories from 3 to 6 weeks",
    y = "Number of PDX samples",
    x = NULL,
    fill = "6-week response"
  ) +
  theme_classic(base_size = 12)

ggsave(
  paste0(output, "Sankey_Response_3wk_to_6wk.png"),
  p_sankey_response,
  width = 8,
  height = 6,
  dpi = 500
)

############################################################
# Figure 1B - tumour-volume change at 3 and 6 weeks
############################################################

scatter_data <- response_dynamic %>%
  filter(
    !is.na(percentage_vol_change_3wk),
    !is.na(percentage_vol_change_6wk)
  )

cor_test <- cor.test(
  scatter_data$percentage_vol_change_3wk,
  scatter_data$percentage_vol_change_6wk,
  method = "spearman"
)

rho <- round(cor_test$estimate, 2)
pval <- signif(cor_test$p.value, 3)

write.csv(
  data.frame(
    N = nrow(scatter_data),
    Spearman_rho = unname(cor_test$estimate),
    P_value = cor_test$p.value
  ),
  paste0(output, "Spearman_Correlation_3wk_vs_6wk.csv"),
  row.names = FALSE
)

p_scatter <- ggplot(
  scatter_data,
  aes(
    x = percentage_vol_change_3wk,
    y = percentage_vol_change_6wk,
    colour = Response_3wk
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = min(scatter_data$percentage_vol_change_3wk),
    y = max(scatter_data$percentage_vol_change_6wk),
    label = paste0(
      "Spearman rho = ", rho,
      "\np = ", pval
    ),
    hjust = 0,
    size = 5
  ) +
  scale_colour_manual(
    values = c(
      "PD" = "#D55E00",
      "SD" = "#0072B2",
      "PR" = "#009E73"
    )
  ) +
  labs(
    title = "Relationship between 3-week and 6-week tumour volume change",
    x = "Tumour volume change at 3 weeks (%)",
    y = "Tumour volume change at 6 weeks (%)",
    colour = "3-week response"
  ) +
  theme_classic(base_size = 12)

ggsave(
  paste0(output, "Scatter_3wk_vs_6wk_Tumour_Change.png"),
  p_scatter,
  width = 7,
  height = 6,
  dpi = 500
)

############################################################
# Response trajectory table
############################################################

trajectory_table <- meta %>%
  dplyr::select(
    Sample_ID,
    Response_3wk,
    Response_6wk
  ) %>%
  filter(!is.na(Response_3wk)) %>%
  mutate(
    Response_6wk = ifelse(
      is.na(Response_6wk),
      "Missing",
      Response_6wk
    ),
    Trajectory_Group = case_when(
      Response_3wk == "PD" & Response_6wk == "PD" ~ "Stable_PD",
      Response_3wk == "SD" & Response_6wk == "SD" ~ "Stable_SD",
      Response_3wk == "PR" & Response_6wk == "PR" ~ "Stable_PR",
      Response_3wk == "PD" & Response_6wk == "Missing" ~ "PD_Missing",
      Response_3wk == "SD" & Response_6wk == "Missing" ~ "SD_Missing",
      Response_3wk == "PR" & Response_6wk == "Missing" ~ "PR_Missing",
      Response_3wk == "PD" & Response_6wk %in% c("SD", "PR") ~ "Improved",
      Response_3wk == "SD" & Response_6wk == "PR" ~ "Improved",
      Response_3wk == "PR" & Response_6wk %in% c("SD", "PD") ~ "Worsened",
      Response_3wk == "SD" & Response_6wk == "PD" ~ "Worsened",
      TRUE ~ "Other"
    )
  )

write.csv(
  trajectory_table,
  paste0(output, "Response_Trajectory_Table.csv"),
  row.names = FALSE
)

paired_trajectory <- trajectory_table %>%
  filter(Response_6wk != "Missing")

trajectory_summary <- data.frame(
  Measure = c(
    "Evaluable at 3 weeks",
    "Evaluable at both timepoints",
    "Missing 6-week classification",
    "Same category at both timepoints",
    "Changed category between timepoints"
  ),
  N = c(
    nrow(trajectory_table),
    nrow(paired_trajectory),
    sum(trajectory_table$Response_6wk == "Missing"),
    sum(paired_trajectory$Response_3wk == paired_trajectory$Response_6wk),
    sum(paired_trajectory$Response_3wk != paired_trajectory$Response_6wk)
  )
)

stopifnot(
  identical(trajectory_summary$N, c(55L, 43L, 12L, 34L, 9L))
)

write.csv(
  trajectory_summary,
  paste0(output, "Response_Trajectory_Summary.csv"),
  row.names = FALSE
)

############################################################
# FIGURE 2 - BASELINE RNA-SEQ QUALITY ASSESSMENT
############################################################

# The outlier is intentionally retained in both Figure 2 panels.

############################################################
# Figure 2A - PCA before outlier removal
############################################################

meta_temp <- meta[!is.na(meta$Definitive_Response_3_6wk), ]
counts_temp <- counts[, meta_temp$Sample_ID]

meta_temp$Definitive_Response_3_6wk <- factor(
  meta_temp$Definitive_Response_3_6wk
)

dds_temp <- DESeqDataSetFromMatrix(
  countData = counts_temp,
  colData = meta_temp,
  design = ~ Definitive_Response_3_6wk
)

dds_temp <- dds_temp[rowSums(counts(dds_temp)) >= 10, ]
vsd_temp <- vst(dds_temp, blind = TRUE)

pca_temp <- plotPCA(
  vsd_temp,
  intgroup = "Definitive_Response_3_6wk",
  returnData = TRUE
)

percentVar_temp <- round(100 * attr(pca_temp, "percentVar"))
pca_temp$Sample_ID <- rownames(pca_temp)

write.csv(
  pca_temp,
  paste0(output, prefix, "_PCA_BeforeOutlierRemoval_Coordinates.csv"),
  row.names = FALSE
)

p_pca_before <- ggplot(
  pca_temp,
  aes(
    x = PC1,
    y = PC2,
    colour = Definitive_Response_3_6wk
  )
) +
  geom_point(size = 4, alpha = 0.85) +
  geom_point(
    data = subset(pca_temp, Sample_ID == outlier_sample),
    aes(x = PC1, y = PC2),
    shape = 21,
    size = 6,
    stroke = 1.2,
    fill = NA,
    colour = "black",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Responder" = "#4DBBD5",
      "NonResponder" = "#F8766D"
    )
  ) +
  labs(
    title = "PCA of Variance-Stabilised Gene Expression",
    subtitle = "3-6 week Average Responders vs NonResponders",
    x = paste0("PC1: ", percentVar_temp[1], "% variance"),
    y = paste0("PC2: ", percentVar_temp[2], "% variance"),
    colour = "Response"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  paste0(output, prefix, "_PCA_BeforeOutlierRemoval.png"),
  p_pca_before,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
# Figure 2B - sequencing depth for all baseline samples
############################################################

library_sizes <- data.frame(
  Sample_ID = colnames(counts),
  Total_Reads = colSums(counts)
) %>%
  arrange(Total_Reads) %>%
  mutate(
    Rank = row_number(),
    Outlier = Sample_ID == outlier_sample
  )

write.csv(
  library_sizes,
  paste0(output, prefix, "_SequencingDepth_Table.csv"),
  row.names = FALSE
)

p_depth <- ggplot(library_sizes, aes(Rank, Total_Reads)) +
  geom_point(
    data = subset(library_sizes, !Outlier),
    colour = "grey65",
    size = 2.6,
    alpha = 0.9
  ) +
  geom_point(
    data = subset(library_sizes, Outlier),
    colour = "#D55E00",
    size = 4
  ) +
  geom_text_repel(
    data = subset(library_sizes, Outlier),
    aes(label = Sample_ID),
    colour = "#D55E00",
    fontface = "bold",
    size = 3.6,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.colour = "#D55E00",
    max.overlaps = Inf
  ) +
  scale_y_continuous(
    labels = label_number(scale = 1e-6, suffix = " M"),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "Library Size Across PDX Samples",
    subtitle = "Samples ranked by total sequencing reads",
    x = "Samples (ranked by sequencing depth)",
    y = "Total sequencing reads"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave(
  paste0(output, prefix, "_SequencingDepth_Ranked.png"),
  p_depth,
  width = 7,
  height = 5,
  dpi = 600
)

quality_summary <- data.frame(
  Measure = c(
    "Number of baseline RNA-seq samples",
    "Cohort median library size",
    "HNC0157PRX0A library size"
  ),
  Value = c(
    ncol(counts),
    median(library_sizes$Total_Reads),
    library_sizes$Total_Reads[library_sizes$Sample_ID == outlier_sample]
  )
)

write.csv(
  quality_summary,
  paste0(output, "Baseline_QC_Summary.csv"),
  row.names = FALSE
)

############################################################
# TABLE 2 - CLINICOPATHOLOGICAL CHARACTERISTICS
############################################################

# Average 3-6-week R vs NR samples only; outlier excluded.

clinic_meta <- meta %>%
  mutate(
    across(
      where(is.character),
      ~ na_if(trimws(.), "")
    )
  ) %>%
  filter(
    !is.na(Definitive_Response_3_6wk),
    Sample_ID != outlier_sample
  )

clinic_meta$Definitive_Response_3_6wk <- factor(
  clinic_meta$Definitive_Response_3_6wk,
  levels = c("Responder", "NonResponder")
)

stopifnot(
  sum(clinic_meta$Definitive_Response_3_6wk == "Responder") == 40L,
  sum(clinic_meta$Definitive_Response_3_6wk == "NonResponder") == 14L
)

clinic_meta <- clinic_meta %>%
  mutate(
    T_stage_group = case_when(
      T_Stage %in% c("T1", "T2") ~ "Early",
      T_Stage %in% c("T3", "T4", "T4a", "T4b") ~ "Advanced",
      TRUE ~ NA_character_
    ),
    N_stage_group = case_when(
      N_Stage == "N0" ~ "N0",
      !is.na(N_Stage) ~ "N+",
      TRUE ~ NA_character_
    ),
    Grade_group = case_when(
      HISTOLOGICAL.GRADE %in% c("g1-g2", "g2") ~ "Low/Intermediate",
      HISTOLOGICAL.GRADE %in% c("g2-g3", "g3") ~ "High",
      TRUE ~ NA_character_
    )
  )

table2_data <- clinic_meta %>%
  dplyr::select(
    Definitive_Response_3_6wk,
    Age,
    Sex,
    Site.of.Primary,
    T_stage_group,
    N_stage_group,
    Grade_group,
    ULCERATION,
    VASCULAR.INVASION,
    PERINEURAL.INVASION,
    LYMPH.NODE.MTS,
    EXTRA.CAPSULAR.EXTENSION,
    RESECTION.MARGINS
  )

table2 <- table2_data %>%
  tbl_summary(
    by = Definitive_Response_3_6wk,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    )
  ) %>%
  bold_labels()

table2

table2_word <- table2 %>%
  as_flex_table()

save_as_docx(
  "Baseline clinicopathological characteristics" = table2_word,
  path = paste0(output, "Baseline_Clinicopathological_Characteristics.docx")
)

writeLines(
  capture.output(sessionInfo()),
  paste0(output, "SessionInfo.txt")
)

message(
  "Finished. All outputs are in: ",
  normalizePath(output, mustWork = FALSE)
)
