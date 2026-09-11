# Longitudinal Analysis
# Adam McClatchey

############################################################
# Longitudinal RNA-seq Analysis
############################################################

rm(list = ls())

library(tidyverse)
library(variancePartition)
library(edgeR)
library(limma)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(grid)

# This path is relative to the working directory shown by getwd().
output_dir <- "Long_RNA/Longitudinal_Results_Clean"
hallmark_output <- file.path(output_dir, "Hallmark_GSEA")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hallmark_output, recursive = TRUE, showWarnings = FALSE)

message(
  "Clean longitudinal outputs will be saved to: ",
  normalizePath(output_dir, mustWork = FALSE)
)

############################################################
# Load longitudinal RNA-seq matrices
############################################################

raw_counts <- read.csv(
  "Long_RNA/FilteredLowCounts97_HPVRem_HNSCC_10_graft_raw_counts_matrix.csv",
  check.names = FALSE
)


############################################################
# Load longitudinal metadata
############################################################

long_meta <- read.csv(
  "Long_RNA/FilteredLowCounts97_HPVRem_GraftHNSCC_colData.csv",
  check.names = FALSE
)

############################################################
# Prepare raw count matrix
############################################################

# First column contains gene names
gene_names <- raw_counts[[1]]

# Remove gene-name column
counts <- raw_counts[, -1]

# Assign gene names as row names
rownames(counts) <- gene_names

# Convert to matrix
counts <- as.matrix(counts)

############################################################
# Align metadata and count matrix
############################################################

stopifnot(setequal(colnames(counts), long_meta$Sample_ID))

long_meta <- long_meta[
  match(colnames(counts), long_meta$Sample_ID),
]

rownames(long_meta) <- long_meta$Sample_ID

stopifnot(identical(rownames(long_meta), colnames(counts)))

############################################################
# Confirm raw counts are integer-valued
############################################################

storage.mode(counts) <- "integer"

############################################################
# Longitudinal Differential Expression Analysis
############################################################

############################################################
# Set longitudinal analysis variables
############################################################

# Sampling stage
long_meta$DefinitiveTreatment <- factor(
  long_meta$DefinitiveTreatment,
  levels = c(
    "Untreated",
    "ShortTreatment",
    "EndPointTreatment"
  )
)

# Average 3-6 week cetuximab response
long_meta$Definitive_Response_3_6wk <- factor(
  long_meta$Definitive_Response_3_6wk,
  levels = c(
    "NonResponder",
    "Responder"
  )
)

# PDX model identifier
long_meta$Case_ID <- factor(long_meta$Case_ID)

############################################################
# Create edgeR object from raw counts
############################################################

dge <- DGEList(
  counts = counts,
  samples = long_meta
)

############################################################
# Expression filtering for longitudinal analysis
############################################################

# Create sampling-stage × response groups
filter_group <- interaction(
  long_meta$DefinitiveTreatment,
  long_meta$Definitive_Response_3_6wk,
  drop = TRUE
)

############################################################
# edgeR expression filtering
############################################################

keep <- filterByExpr(
  dge,
  group = filter_group
)

message(sum(keep), " of ", length(keep), " genes retained by filterByExpr().")

############################################################
# Apply expression filtering
############################################################

dge <- dge[
  keep,
  ,
  keep.lib.sizes = FALSE
]

############################################################
# TMM normalisation
############################################################

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

############################################################
# Create stage × response groups
############################################################

long_meta$StageResponse <- interaction(
  long_meta$DefinitiveTreatment,
  long_meta$Definitive_Response_3_6wk,
  sep = "_"
)

############################################################
# Longitudinal model using explicit groups
############################################################

form_group <- ~ 0 + StageResponse + (1 | Case_ID)

############################################################
# Longitudinal contrasts
############################################################

L <- makeContrastsDream(
  form_group,
  long_meta,
  contrasts = c(
    
    NR_Short_vs_Untreated =
      "StageResponseShortTreatment_NonResponder -
       StageResponseUntreated_NonResponder",
    
    NR_Endpoint_vs_Untreated =
      "StageResponseEndPointTreatment_NonResponder -
       StageResponseUntreated_NonResponder",
    
    NR_Endpoint_vs_Short =
      "StageResponseEndPointTreatment_NonResponder -
       StageResponseShortTreatment_NonResponder",
    
    R_Short_vs_Untreated =
      "StageResponseShortTreatment_Responder -
       StageResponseUntreated_Responder",
    
    R_Endpoint_vs_Untreated =
      "StageResponseEndPointTreatment_Responder -
       StageResponseUntreated_Responder",
    
    R_Endpoint_vs_Short =
      "StageResponseEndPointTreatment_Responder -
       StageResponseShortTreatment_Responder"
  )
)

############################################################
# Voom weights for final longitudinal model
############################################################

vobj_group <- voomWithDreamWeights(
  dge,
  form_group,
  long_meta
)

############################################################
# Fit final longitudinal mixed model with all contrasts
############################################################

fit_long <- dream(
  vobj_group,
  form_group,
  long_meta,
  L = L
)

fit_long <- eBayes(fit_long)

############################################################
# Save final model
############################################################

saveRDS(
  fit_long,
  file.path(output_dir, "final_longitudinal_dream_fit.rds")
)

saveRDS(
  vobj_group,
  file.path(output_dir, "final_longitudinal_voom_object.rds")
)

############################################################
# Extract all six longitudinal contrasts
############################################################

contrast_names <- c(
  "NR_Short_vs_Untreated",
  "NR_Endpoint_vs_Untreated",
  "NR_Endpoint_vs_Short",
  "R_Short_vs_Untreated",
  "R_Endpoint_vs_Untreated",
  "R_Endpoint_vs_Short"
)

long_results <- lapply(contrast_names, function(x) {
  
  res <- topTable(
    fit_long,
    coef = x,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  res$Gene <- rownames(res)
  
  res <- res %>%
    dplyr::select(Gene, everything())
  
  return(res)
})

names(long_results) <- contrast_names

############################################################
# Summarise significant longitudinal DEGs
############################################################

long_summary <- data.frame(
  Comparison = contrast_names,
  
  Significant = sapply(
    long_results,
    function(x) sum(x$adj.P.Val < 0.05)
  ),
  
  Up = sapply(
    long_results,
    function(x) sum(x$adj.P.Val < 0.05 & x$logFC > 0)
  ),
  
  Down = sapply(
    long_results,
    function(x) sum(x$adj.P.Val < 0.05 & x$logFC < 0)
  )
)

long_summary

############################################################
# Significant early longitudinal DEGs in NonResponders
############################################################

nr_early_degs <- long_results$NR_Short_vs_Untreated %>%
  filter(adj.P.Val < 0.05)

############################################################
# Significant responder gene sets
############################################################

R_early <- long_results$R_Short_vs_Untreated %>%
  filter(adj.P.Val < 0.05) %>%
  pull(Gene)

R_overall <- long_results$R_Endpoint_vs_Untreated %>%
  filter(adj.P.Val < 0.05) %>%
  pull(Gene)

R_late <- long_results$R_Endpoint_vs_Short %>%
  filter(adj.P.Val < 0.05) %>%
  pull(Gene)

responder_endpoint_overlap <- length(intersect(R_overall, R_late))
message(responder_endpoint_overlap, " responder DEGs shared between endpoint contrasts.")

############################################################
# SAVE FINAL LONGITUDINAL RESULTS
############################################################


############################################################
# 1. SAVE COMPLETE R OBJECT
# Contains the FULL results tables for ALL SIX comparisons.
# Each table contains all 14,902 tested genes with:
# Gene, logFC, average expression, test statistic,
# raw P-value and adjusted P-value.
#
# This is the main object to reload into R later so the
# results do not need to be regenerated.
############################################################

saveRDS(
  long_results,
  file.path(output_dir, "longitudinal_all_results.rds")
)


############################################################
# 2. CREATE FOLDER FOR CSV RESULTS
############################################################

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
# 3. SAVE ALL GENES FOR EACH OF THE SIX COMPARISONS
# Creates SIX CSV files.
#
# Each file contains ALL 14,902 genes, not just significant
# genes. These are useful for:
# - GSEA
# - volcano plots
# - checking individual genes
# - using different significance thresholds later
############################################################

for (nm in names(long_results)) {
  
  write.csv(
    long_results[[nm]],
    file.path(output_dir, paste0(nm, "_ALL_genes.csv")),
    row.names = FALSE
  )
}


############################################################
# 4. SAVE SIGNIFICANT DEGs FOR EACH COMPARISON
# Creates SIX additional CSV files.
#
# Each contains ONLY genes with BH-adjusted P < 0.05.
# Therefore these are the final significant DEG lists.
#
# Some NR files will be empty because no genes reached
# FDR < 0.05 in those comparisons.
############################################################

for (nm in names(long_results)) {
  
  sig <- long_results[[nm]] %>%
    dplyr::filter(adj.P.Val < 0.05)
  
  write.csv(
    sig,
    file.path(output_dir, paste0(nm, "_FDR05_DEGs.csv")),
    row.names = FALSE
  )
}


############################################################
# 5. SAVE DEG COUNT SUMMARY
# One small CSV summarising the number of significant genes
# in each comparison, including numbers increased and
# decreased over time.
#
# e.g. R Endpoint vs Untreated:
# 1366 total = 893 increased + 473 decreased.
############################################################

write.csv(
  long_summary,
  file.path(output_dir, "DEG_summary.csv"),
  row.names = FALSE
)


############################################################


############################################################
# HALLMARK GSEA - LONGITUDINAL ANALYSIS
############################################################

############################################################
# Retrieve Hallmark gene sets
############################################################

h_t2g <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
) %>%
  dplyr::select(gs_name, gene_symbol)

############################################################
# Create Hallmark GSEA output folder
############################################################

dir.create(hallmark_output, recursive = TRUE, showWarnings = FALSE)

############################################################
# Run Hallmark GSEA for all six longitudinal contrasts
############################################################

hallmark_results <- list()

for (comparison in names(long_results)) {
  
  # Get complete dream results for this longitudinal contrast
  res_gsea <- long_results[[comparison]]
  
  # Remove genes without a valid t-statistic
  res_gsea <- res_gsea[
    !is.na(res_gsea$t) &
      !is.na(res_gsea$Gene),
  ]
  
  # Build ranked gene list using dream t-statistic
  geneList <- res_gsea$t
  names(geneList) <- res_gsea$Gene
  
  # Rank genes from most positive to most negative change
  geneList <- sort(geneList, decreasing = TRUE)
  
  # Run Hallmark GSEA
  gsea_result <- GSEA(
    geneList,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 10000,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = h_t2g,
    verbose = FALSE,
    seed = TRUE
  )
  
  # Convert result to dataframe
  gsea_df <- as.data.frame(gsea_result)
  
  # Store in R
  hallmark_results[[comparison]] <- gsea_df
  
  # Save complete Hallmark results for this comparison
  write.csv(
    gsea_df,
    file.path(
      hallmark_output,
      paste0(comparison, "_Hallmark_GSEA_ALL.csv")
    ),
    row.names = FALSE
  )
  
  # Save significant Hallmark pathways only
  gsea_sig <- gsea_df %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(desc(NES))
  
  write.csv(
    gsea_sig,
    file.path(
      hallmark_output,
      paste0(comparison, "_Hallmark_GSEA_FDR05.csv")
    ),
    row.names = FALSE
  )
}

############################################################
# Save all six Hallmark GSEA results as one R object
############################################################

saveRDS(
  hallmark_results,
  file.path(hallmark_output, "Hallmark_GSEA_all_six_comparisons.rds")
)

############################################################
# Number of significant Hallmark pathways in each comparison
############################################################

hallmark_summary <- data.frame(
  Comparison = names(hallmark_results),
  
  Significant_Pathways = sapply(
    hallmark_results,
    function(x) sum(x$p.adjust < 0.05)
  ),
  
  Positive_NES = sapply(
    hallmark_results,
    function(x) sum(x$p.adjust < 0.05 & x$NES > 0)
  ),
  
  Negative_NES = sapply(
    hallmark_results,
    function(x) sum(x$p.adjust < 0.05 & x$NES < 0)
  )
)

hallmark_summary

write.csv(
  hallmark_summary,
  file.path(hallmark_output, "Hallmark_GSEA_summary.csv"),
  row.names = FALSE
)

# OVERALL LONGITUDINAL ANALYSIS - ALL PDX MODELS
#
# Biological question:
# What transcriptional changes occur broadly across the
# treatment course, irrespective of eventual response?
#
# Reuses the already filtered and TMM-normalised dge object.
############################################################


############################################################
# 1. DEFINE OVERALL LONGITUDINAL MODEL
#
# DefinitiveTreatment = sampling stage:
# Untreated / ShortTreatment / EndPointTreatment
#
# (1 | Case_ID) = random intercept for each PDX model,
# accounting for repeated measurements from the same model.
############################################################

form_all <- ~ 0 + DefinitiveTreatment + (1 | Case_ID)


############################################################
# 2. DEFINE THE THREE LONGITUDINAL CONTRASTS
#
# Short vs Untreated:
# early treatment-associated change
#
# Endpoint vs Untreated:
# overall change from baseline to 6-week endpoint
#
# Endpoint vs Short:
# later change occurring after the early sampling stage
############################################################

L_all <- makeContrastsDream(
  form_all,
  long_meta,
  contrasts = c(
    
    All_Short_vs_Untreated =
      "DefinitiveTreatmentShortTreatment -
       DefinitiveTreatmentUntreated",
    
    All_Endpoint_vs_Untreated =
      "DefinitiveTreatmentEndPointTreatment -
       DefinitiveTreatmentUntreated",
    
    All_Endpoint_vs_Short =
      "DefinitiveTreatmentEndPointTreatment -
       DefinitiveTreatmentShortTreatment"
  )
)

# Check contrast matrix
L_all


############################################################
# 3. VOOM TRANSFORMATION + DREAM PRECISION WEIGHTS
#
# Uses the filtered, TMM-normalised count data and estimates
# the RNA-seq mean-variance relationship and precision weights
# for the overall repeated-measures model.
############################################################

vobj_all <- voomWithDreamWeights(
  dge,
  form_all,
  long_meta
)


############################################################
# 4. FIT REPEATED-MEASURES MIXED MODEL
#
# Fits all three contrasts simultaneously while accounting
# for repeated measurements within Case_ID.
############################################################

fit_all <- dream(
  vobj_all,
  form_all,
  long_meta,
  L = L_all
)

# Empirical Bayes moderation appropriate for dream models
fit_all <- variancePartition::eBayes(fit_all)


############################################################
# 5. SAVE EXPENSIVE MODEL OBJECTS
#
# Allows model/results to be reloaded without rerunning dream.
############################################################

saveRDS(
  fit_all,
  file.path(output_dir, "overall_longitudinal_dream_fit.rds")
)

saveRDS(
  vobj_all,
  file.path(output_dir, "overall_longitudinal_voom_object.rds")
)


############################################################
# 6. EXTRACT COMPLETE RESULTS FOR ALL THREE COMPARISONS
#
# Each table contains all 14,902 tested genes with effect
# size, test statistic, raw P-value and adjusted P-value.
############################################################

overall_contrast_names <- c(
  "All_Short_vs_Untreated",
  "All_Endpoint_vs_Untreated",
  "All_Endpoint_vs_Short"
)

overall_results <- lapply(
  overall_contrast_names,
  function(x) {
    
    res <- variancePartition::topTable(
      fit_all,
      coef = x,
      number = Inf,
      sort.by = "P"
    )
    
    res$Gene <- rownames(res)
    
    res %>%
      dplyr::select(Gene, everything())
  }
)

names(overall_results) <- overall_contrast_names


############################################################
# 7. SUMMARISE SIGNIFICANT OVERALL LONGITUDINAL DEGs
#
# FDR < 0.05
# Positive logFC = expression increased at the later stage
# Negative logFC = expression decreased at the later stage
############################################################

overall_summary <- data.frame(
  
  Comparison = overall_contrast_names,
  
  Significant = sapply(
    overall_results,
    function(x) sum(x$adj.P.Val < 0.05)
  ),
  
  Up = sapply(
    overall_results,
    function(x)
      sum(x$adj.P.Val < 0.05 & x$logFC > 0)
  ),
  
  Down = sapply(
    overall_results,
    function(x)
      sum(x$adj.P.Val < 0.05 & x$logFC < 0)
  )
)

overall_summary


############################################################
# 8. SAVE COMPLETE GENE-LEVEL RESULTS
############################################################

saveRDS(
  overall_results,
  file.path(output_dir, "overall_longitudinal_all_results.rds")
)

write.csv(
  overall_summary,
  file.path(output_dir, "overall_DEG_summary.csv"),
  row.names = FALSE
)

for (nm in names(overall_results)) {
  
  # All tested genes
  write.csv(
    overall_results[[nm]],
    file.path(output_dir, paste0(nm, "_ALL_genes.csv")),
    row.names = FALSE
  )
  
  # FDR-significant genes only
  write.csv(
    overall_results[[nm]] %>%
      dplyr::filter(adj.P.Val < 0.05),
    file.path(output_dir, paste0(nm, "_FDR05_DEGs.csv")),
    row.names = FALSE
  )
}

############################################################
# OVERALL LONGITUDINAL HALLMARK GSEA
#
# Biological question:
# Which biological programmes change across treatment
# irrespective of Responder / NonResponder phenotype?
#
# Uses ALL genes ranked by the dream t-statistic.
############################################################

overall_hallmark <- list()

for (comparison in names(overall_results)) {
  
  res_gsea <- overall_results[[comparison]]
  
  # Keep genes with valid gene names and t-statistics
  res_gsea <- res_gsea[
    !is.na(res_gsea$t) &
      !is.na(res_gsea$Gene),
  ]
  
  # Rank all genes by dream t-statistic
  geneList <- res_gsea$t
  names(geneList) <- res_gsea$Gene
  
  geneList <- sort(
    geneList,
    decreasing = TRUE
  )
  
  # Run Hallmark GSEA
  gsea_result <- GSEA(
    geneList,
    exponent = 1,
    minGSSize = 10,
    maxGSSize = 10000,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    TERM2GENE = h_t2g,
    verbose = FALSE,
    seed = TRUE
  )
  
  overall_hallmark[[comparison]] <-
    as.data.frame(gsea_result)
}


############################################################
# SUMMARISE SIGNIFICANT HALLMARK PATHWAYS
############################################################

overall_hallmark_summary <- data.frame(
  
  Comparison = names(overall_hallmark),
  
  Significant_Pathways = sapply(
    overall_hallmark,
    function(x) sum(x$p.adjust < 0.05)
  ),
  
  Positive_NES = sapply(
    overall_hallmark,
    function(x)
      sum(x$p.adjust < 0.05 & x$NES > 0)
  ),
  
  Negative_NES = sapply(
    overall_hallmark,
    function(x)
      sum(x$p.adjust < 0.05 & x$NES < 0)
  )
)

overall_hallmark_summary

############################################################
# SAVE OVERALL LONGITUDINAL HALLMARK GSEA RESULTS
#
# Saves complete results, significant pathways only,
# summary counts, and the complete R object.
############################################################

# Complete GSEA results for each comparison
for (nm in names(overall_hallmark)) {
  
  write.csv(
    overall_hallmark[[nm]],
    file.path(
      hallmark_output,
      paste0(nm, "_Hallmark_GSEA_ALL.csv")
    ),
    row.names = FALSE
  )
}


# Significant pathways only (FDR < 0.05)
for (nm in names(overall_hallmark)) {
  
  write.csv(
    overall_hallmark[[nm]] %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::arrange(NES),
    file.path(
      hallmark_output,
      paste0(nm, "_Hallmark_GSEA_SIGNIFICANT.csv")
    ),
    row.names = FALSE
  )
}


# Summary of significant pathway counts and directions
write.csv(
  overall_hallmark_summary,
  file.path(
    hallmark_output,
    "Overall_Longitudinal_Hallmark_GSEA_SUMMARY.csv"
  ),
  row.names = FALSE
)


# Complete GSEA objects for easy reloading later
saveRDS(
  overall_hallmark,
  file.path(
    hallmark_output,
    "Overall_Longitudinal_Hallmark_GSEA_all_results.rds"
  )
)





############################################################
# FIGURE 17 - COHORT-WIDE LONGITUDINAL HALLMARK HEATMAP
#
# Includes every Hallmark pathway significant (FDR < 0.05)
# in at least one of the three overall longitudinal contrasts.
#
# Rows are ordered by:
#   1. Number of comparisons in which pathway is significant
#   2. Mean absolute NES across the three comparisons
#
# * indicates FDR < 0.05
############################################################



############################################################
# 1. Combine the three complete GSEA result tables
############################################################

hallmark_long_df <- bind_rows(
  
  overall_hallmark$All_Short_vs_Untreated %>%
    mutate(Comparison = "Untreated \u2192 Short"),
  
  overall_hallmark$All_Endpoint_vs_Untreated %>%
    mutate(Comparison = "Untreated \u2192 Endpoint"),
  
  overall_hallmark$All_Endpoint_vs_Short %>%
    mutate(Comparison = "Short \u2192 Endpoint")
) %>%
  
  select(
    Description,
    Comparison,
    NES,
    p.adjust
  )


############################################################
# 2. Identify pathways significant in >=1 comparison
############################################################

pathway_summary <- hallmark_long_df %>%
  
  group_by(Description) %>%
  
  summarise(
    
    Analyses_Present = sum(
      p.adjust < 0.05,
      na.rm = TRUE
    ),
    
    Mean_Abs_NES = mean(
      abs(NES),
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  filter(Analyses_Present >= 1) %>%
  
  arrange(
    desc(Analyses_Present),
    desc(Mean_Abs_NES)
  )


############################################################
# 3. Keep those pathways and prepare labels
############################################################

heatmap_df <- hallmark_long_df %>%
  
  filter(
    Description %in% pathway_summary$Description
  ) %>%
  
  mutate(
    
    # Clean Hallmark pathway names
    Pathway = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    
    # Mark statistically significant cells
    Significance = ifelse(
      p.adjust < 0.05,
      "*",
      ""
    ),
    
    # Fix chronological column order
    Comparison = factor(
      Comparison,
      levels = c(
        "Untreated \u2192 Short",
        "Untreated \u2192 Endpoint",
        "Short \u2192 Endpoint"
      )
    )
  )


############################################################
# 4. Set pathway row order
############################################################

pathway_order <- pathway_summary %>%
  
  mutate(
    Pathway = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title()
  ) %>%
  
  pull(Pathway)


heatmap_df$Pathway <- factor(
  heatmap_df$Pathway,
  levels = rev(pathway_order)
)


############################################################
# 5. Determine symmetrical NES scale
############################################################

nes_limit <- max(
  abs(heatmap_df$NES),
  na.rm = TRUE
)


############################################################
# 6. Create thesis-ready heatmap
############################################################

p_hallmark_long <- ggplot(
  heatmap_df,
  aes(
    x = Comparison,
    y = Pathway,
    fill = NES
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +
  
  # Asterisk = significant FDR < 0.05
  geom_text(
    aes(label = Significance),
    size = 5,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-nes_limit, nes_limit),
    name = "NES"
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    title = "Longitudinal Hallmark Pathway Enrichment Across Treatment"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 15)
    ),
    
    axis.text.x = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    
    axis.ticks = element_blank(),
    
    axis.line = element_blank(),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


p_hallmark_long

############################################################
# SAVE FIGURE
############################################################

ggsave(
  file.path(output_dir, "Overall_Longitudinal_Hallmark_Heatmap.png"),
  p_hallmark_long,
  width = 10,
  height = 11,
  dpi = 500
)

############################################################
# FIGURE 18 - SAMPLE-LEVEL HALLMARK ssGSEA TRAJECTORIES
#
# Purpose:
# Generate one Hallmark enrichment score for every pathway
# in every RNA sample.
#
# This creates genuine sample-level pathway measurements
# that can be followed across:
# Untreated -> ShortTreatment -> EndPointTreatment
############################################################



############################################################
# 1. CREATE HALLMARK GENE-SET LIST
#
# h_t2g was already created for the longitudinal GSEA.
############################################################

hallmark_list <- split(
  h_t2g$gene_symbol,
  h_t2g$gs_name
)

length(hallmark_list)

############################################################
# 2. PREPARE EXPRESSION MATRIX
#
# vobj_all$E:
# genes x 152 samples
# voom log2-expression values from the overall model.
############################################################

ssgsea_expr <- vobj_all$E

dim(ssgsea_expr)

# Confirm sample alignment
all(colnames(ssgsea_expr) == rownames(long_meta))

############################################################
# 3. CALCULATE ssGSEA SCORES
#
# Output:
# Hallmark pathways x 152 individual RNA samples
############################################################

ssgsea_par <- ssgseaParam(
  exprData = ssgsea_expr,
  geneSets = hallmark_list,
  minSize = 10,
  maxSize = 10000,
  normalize = TRUE
)

ssgsea_scores <- gsva(
  ssgsea_par,
  verbose = TRUE
)

dim(ssgsea_scores)

############################################################
# 4. SAVE COMPLETE SAMPLE-LEVEL ssGSEA SCORES
############################################################

write.csv(
  ssgsea_scores,
  file.path(output_dir, "Hallmark_ssGSEA_SampleScores.csv")
)

saveRDS(
  ssgsea_scores,
  file.path(output_dir, "Hallmark_ssGSEA_SampleScores.rds")
)

############################################################
# 5. SELECT FOUR REPRESENTATIVE HALLMARK PATHWAYS
############################################################

selected_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
)

selected_names <- c(
  "Interferon Alpha Response",
  "E2F Targets",
  "TGF-beta Signalling",
  "PI3K/AKT/mTOR Signalling"
)


############################################################
# 6. CONVERT ssGSEA SCORES TO LONG FORMAT
############################################################

trajectory_df <- as.data.frame(
  t(ssgsea_scores[selected_pathways, ])
)

trajectory_df$Sample_ID <- rownames(trajectory_df)

trajectory_df <- trajectory_df %>%
  
  pivot_longer(
    cols = all_of(selected_pathways),
    names_to = "Pathway",
    values_to = "ssGSEA_Score"
  ) %>%
  
  left_join(
    long_meta %>%
      dplyr::select(
        Sample_ID,
        Case_ID,
        DefinitiveTreatment
      ),
    by = "Sample_ID"
  ) %>%
  
  mutate(
    
    DefinitiveTreatment = factor(
      DefinitiveTreatment,
      levels = c(
        "Untreated",
        "ShortTreatment",
        "EndPointTreatment"
      ),
      labels = c(
        "Untreated",
        "Short",
        "Endpoint"
      )
    ),
    
    Pathway = factor(
      Pathway,
      levels = selected_pathways,
      labels = selected_names
    )
  )

############################################################
# 7. CALCULATE MEAN ssGSEA TRAJECTORY
############################################################

trajectory_mean <- trajectory_df %>%
  
  group_by(
    Pathway,
    DefinitiveTreatment
  ) %>%
  
  summarise(
    Mean = mean(ssGSEA_Score),
    SE = sd(ssGSEA_Score) / sqrt(n()),
    .groups = "drop"
  )

############################################################
# 8. PLOT SAMPLE-LEVEL PATHWAY TRAJECTORIES
#
# Thin grey lines = individual PDX models
# Thick black line = cohort mean
# Error bars = mean +/- standard error
############################################################

p_ssgsea_trajectory <- ggplot(
  trajectory_df,
  aes(
    x = DefinitiveTreatment,
    y = ssGSEA_Score
  )
) +
  
  ##########################################################
# Individual PDX trajectories
##########################################################

geom_line(
  aes(group = Case_ID),
  colour = "grey65",
  alpha = 0.25,
  linewidth = 0.45
) +
  
  geom_point(
    colour = "grey55",
    alpha = 0.25,
    size = 1
  ) +
  
  ##########################################################
# Overall cohort mean trajectory
##########################################################

geom_line(
  data = trajectory_mean,
  aes(
    x = DefinitiveTreatment,
    y = Mean,
    group = 1
  ),
  inherit.aes = FALSE,
  linewidth = 1.25
) +
  
  geom_point(
    data = trajectory_mean,
    aes(
      x = DefinitiveTreatment,
      y = Mean
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  
  ##########################################################
# Standard error around cohort mean
##########################################################

geom_errorbar(
  data = trajectory_mean,
  aes(
    x = DefinitiveTreatment,
    ymin = Mean - SE,
    ymax = Mean + SE
  ),
  inherit.aes = FALSE,
  width = 0.08,
  linewidth = 0.7
) +
  
  ##########################################################
# Four pathway panels
##########################################################

facet_wrap(
  ~ Pathway,
  ncol = 2,
  scales = "free_y"
) +
  
  ##########################################################
# Labels
##########################################################

labs(
  x = "Sampling stage",
  y = "ssGSEA enrichment score",
  title = "Longitudinal Hallmark Pathway Trajectories"
) +
  
  ##########################################################
# Formatting
##########################################################

theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.title = element_text(
      face = "bold"
    )
  )


############################################################
# DISPLAY FIGURE
############################################################

p_ssgsea_trajectory

ggsave(
  file.path(output_dir, "Overall_Hallmark_ssGSEA_Trajectories.png"),
  p_ssgsea_trajectory,
  width = 10,
  height = 8,
  dpi = 500
)



############################################################
# FIGURE 19
# Longitudinal expression of top DEGs in Responders
#
# Genes selected as the four largest absolute log2FC values
# among FDR-significant genes in:
# Responder Endpoint vs Untreated
#
# PIGR    log2FC = +2.98
# MFSD4A  log2FC = +2.81
# PSCA    log2FC = +2.68
# IGFN1   log2FC = -2.54
############################################################


############################################################
# 1. SELECT GENES
############################################################

selected_genes <- c(
  "PIGR",
  "MFSD4A",
  "PSCA",
  "IGFN1"
)

############################################################
# 2. EXTRACT TMM-NORMALISED LOG2-CPM EXPRESSION
#    FROM EXISTING VOOM OBJECT
############################################################

gene_expression <- as.data.frame(
  t(vobj_group$E[selected_genes, , drop = FALSE])
)

gene_expression$Sample_ID <- rownames(gene_expression)

############################################################
# 3. CONVERT TO LONG FORMAT AND ADD METADATA
############################################################

gene_plot_df <- gene_expression %>%
  
  pivot_longer(
    cols = all_of(selected_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  
  left_join(
    long_meta %>%
      dplyr::select(
        Sample_ID,
        Case_ID,
        DefinitiveTreatment,
        Definitive_Response_3_6wk
      ),
    by = "Sample_ID"
  )

############################################################
# 4. KEEP RESPONDERS ONLY
############################################################

gene_plot_df_R <- gene_plot_df %>%
  
  filter(
    Definitive_Response_3_6wk == "Responder"
  ) %>%
  
  mutate(
    
    # Treatment stage
    Treatment = factor(
      DefinitiveTreatment,
      levels = c(
        "Untreated",
        "ShortTreatment",
        "EndPointTreatment"
      ),
      labels = c(
        "Untreated",
        "Short treatment",
        "Endpoint"
      )
    ),
    
    # Gene order
    Gene = factor(
      Gene,
      levels = selected_genes
    )
  ) %>%
  
  filter(
    !is.na(Treatment),
    !is.na(Expression)
  )

############################################################
# 5. CHECK SAMPLE NUMBERS
############################################################

gene_plot_df_R %>%
  dplyr::count(Gene, Treatment)

############################################################
# 6. CREATE FIGURE
############################################################

p_gene_trajectory <- ggplot(
  gene_plot_df_R,
  aes(
    x = Treatment,
    y = Expression,
    fill = Treatment
  )
) +
  
  ##########################################################
# Transparent coloured boxplots
##########################################################

geom_boxplot(
  width = 0.58,
  alpha = 0.40,
  outlier.shape = NA,
  linewidth = 0.7,
  colour = "black"
) +
  
  ##########################################################
# Individual samples
##########################################################

geom_jitter(
  width = 0.10,
  height = 0,
  size = 1.25,
  alpha = 0.65,
  colour = "black"
) +
  
  ##########################################################
# Four gene panels
##########################################################

facet_wrap(
  ~ Gene,
  ncol = 2,
  scales = "free_y"
) +
  
  ##########################################################
# Treatment-stage colours
##########################################################

scale_fill_manual(
  values = c(
    "Untreated" = "#5B8FF9",
    "Short treatment" = "#61B88A",
    "Endpoint" = "#B37AC5"
  )
) +
  
  ##########################################################
# Labels
##########################################################

labs(
  title = "Longitudinal Expression of Top Differentially Expressed Genes in Responders",
  x = "Treatment stage",
  y = "Expression (TMM-normalised log2-CPM)",
  fill = "Treatment stage"
) +
  
  ##########################################################
# Thesis formatting
##########################################################

theme_classic(
  base_size = 14
) +
  
  theme(
    
    # Main title
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 14)
    ),
    
    # Gene labels
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.7
    ),
    
    strip.text = element_text(
      face = "bold.italic",
      size = 13,
      margin = margin(
        t = 5,
        b = 5
      )
    ),
    
    # Axis titles
    axis.title.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 13,
      margin = margin(r = 10)
    ),
    
    # Axis text
    axis.text.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    # Axis lines/ticks
    axis.line = element_line(
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    # Spacing between panels
    panel.spacing = unit(
      1.2,
      "lines"
    ),
    
    # Legend
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    )
  )

############################################################
# 7. DISPLAY
############################################################

p_gene_trajectory

############################################################
# 8. SAVE
############################################################

ggsave(
  file.path(output_dir, "Responder_Top_DEGs_Longitudinal_Expression.png"),
  plot = p_gene_trajectory,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white"
)

############################################################
# FIGURE 20
# Longitudinal Hallmark pathway activity in
# Responders and NonResponders
#
# Uses existing sample-level Hallmark ssGSEA scores
#
# Selected pathways:
# 1. Interferon-alpha response - early immune response
# 2. E2F targets              - proliferation/cell cycle
# 3. mTORC1 signalling        - growth/metabolic signalling
# 4. TGF-beta signalling      - growth/phenotypic signalling
############################################################


############################################################
# 1. SELECT HALLMARK PATHWAYS
############################################################

selected_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_TGF_BETA_SIGNALING"
)

############################################################
# 2. CHECK THAT ALL FOUR EXIST IN ssGSEA MATRIX
############################################################

selected_pathways %in% rownames(ssgsea_scores)

# Should return:
# TRUE TRUE TRUE TRUE

############################################################
# 3. EXTRACT SAMPLE-LEVEL ssGSEA SCORES
############################################################

pathway_scores <- as.data.frame(
  t(
    ssgsea_scores[
      selected_pathways,
      ,
      drop = FALSE
    ]
  )
)

pathway_scores$Sample_ID <- rownames(pathway_scores)

############################################################
# 4. CONVERT TO LONG FORMAT AND ADD METADATA
############################################################

pathway_plot_df <- pathway_scores %>%
  
  pivot_longer(
    cols = all_of(selected_pathways),
    names_to = "Pathway",
    values_to = "ssGSEA_score"
  ) %>%
  
  left_join(
    long_meta %>%
      dplyr::select(
        Sample_ID,
        DefinitiveTreatment,
        Definitive_Response_3_6wk
      ),
    by = "Sample_ID"
  ) %>%
  
  mutate(
    
    ########################################################
    # Treatment order
    ########################################################
    
    Treatment = factor(
      DefinitiveTreatment,
      levels = c(
        "Untreated",
        "ShortTreatment",
        "EndPointTreatment"
      ),
      labels = c(
        "Untreated",
        "Short treatment",
        "Endpoint"
      )
    ),
    
    ########################################################
    # Response order
    ########################################################
    
    Response = factor(
      Definitive_Response_3_6wk,
      levels = c(
        "Responder",
        "NonResponder"
      ),
      labels = c(
        "Responder",
        "Non-responder"
      )
    ),
    
    ########################################################
    # Clean pathway names
    ########################################################
    
    Pathway = recode(
      Pathway,
      "HALLMARK_INTERFERON_ALPHA_RESPONSE" =
        "Interferon-α Response",
      "HALLMARK_E2F_TARGETS" =
        "E2F Targets",
      "HALLMARK_MTORC1_SIGNALING" =
        "mTORC1 Signalling",
      "HALLMARK_TGF_BETA_SIGNALING" =
        "TGF-β Signalling"
    ),
    
    Pathway = factor(
      Pathway,
      levels = c(
        "Interferon-α Response",
        "E2F Targets",
        "mTORC1 Signalling",
        "TGF-β Signalling"
      )
    )
  ) %>%
  
  filter(
    !is.na(Treatment),
    !is.na(Response),
    !is.na(ssGSEA_score)
  )

############################################################
# 5. CALCULATE GROUP MEANS AND 95% CONFIDENCE INTERVALS
############################################################

pathway_summary <- pathway_plot_df %>%
  
  group_by(
    Pathway,
    Response,
    Treatment
  ) %>%
  
  summarise(
    
    n = n(),
    
    mean_ssGSEA = mean(
      ssGSEA_score,
      na.rm = TRUE
    ),
    
    sd_ssGSEA = sd(
      ssGSEA_score,
      na.rm = TRUE
    ),
    
    se_ssGSEA = sd_ssGSEA / sqrt(n),
    
    CI_lower = mean_ssGSEA -
      qt(0.975, df = n - 1) * se_ssGSEA,
    
    CI_upper = mean_ssGSEA +
      qt(0.975, df = n - 1) * se_ssGSEA,
    
    .groups = "drop"
  )

############################################################
# 6. CHECK SUMMARY VALUES
############################################################

print(pathway_summary, n = Inf)

############################################################
# 7. CREATE LONGITUDINAL TRAJECTORY FIGURE
############################################################

p_pathway_trajectory <- ggplot(
  pathway_summary,
  aes(
    x = Treatment,
    y = mean_ssGSEA,
    colour = Response,
    group = Response
  )
) +
  
  ##########################################################
# 95% confidence intervals
##########################################################

geom_errorbar(
  aes(
    ymin = CI_lower,
    ymax = CI_upper
  ),
  width = 0.06,
  linewidth = 0.5,
  alpha = 0.55,
  position = position_dodge(width = 0.10)
) +
  
  ##########################################################
# Mean trajectory
##########################################################

geom_line(
  linewidth = 1.1,
  position = position_dodge(width = 0.10)
) +
  
  ##########################################################
# Mean at each treatment stage
##########################################################

geom_point(
  size = 3.0,
  position = position_dodge(width = 0.10)
) +
  
  ##########################################################
# Four pathway panels
##########################################################

facet_wrap(
  ~ Pathway,
  ncol = 2,
  scales = "free_y"
) +
  
  ##########################################################
# Response colours
##########################################################

scale_colour_manual(
  values = c(
    "Responder" = "#D55E00",
    "Non-responder" = "#0072B2"
  )
) +
  
  ##########################################################
# Labels
##########################################################

labs(
  title = "Longitudinal Hallmark Pathway Activity by Treatment Response",
  x = "Treatment stage",
  y = "Hallmark ssGSEA enrichment score",
  colour = "Response"
) +
  
  ##########################################################
# Thesis formatting
##########################################################

theme_classic(
  base_size = 14
) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16,
      margin = margin(b = 14)
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.7
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12.5,
      margin = margin(
        t = 5,
        b = 5
      )
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 10)
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 13,
      margin = margin(r = 10)
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 10.5
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.line = element_line(
      linewidth = 0.6,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    panel.spacing = unit(
      1.2,
      "lines"
    ),
    
    legend.position = "top",
    
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    
    legend.text = element_text(
      size = 10.5
    )
  )

############################################################
# 8. DISPLAY
############################################################

p_pathway_trajectory

############################################################
# 9. SAVE
############################################################

ggsave(
  file.path(output_dir, "Response_Specific_Hallmark_Trajectories.png"),
  plot = p_pathway_trajectory,
  width = 10,
  height = 8,
  dpi = 500,
  bg = "white"
)

# TOP 2 UP + TOP 2 DOWN
# Overall Endpoint vs Untreated
############################################################

overall_endpoint <- overall_results[["All_Endpoint_vs_Untreated"]] %>%
  dplyr::filter(adj.P.Val < 0.05)

# Top 2 upregulated by logFC
top2_up <- overall_endpoint %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 2)

# Top 2 downregulated by logFC
top2_down <- overall_endpoint %>%
  arrange(logFC) %>%
  slice_head(n = 2)

selected_overall_genes <- bind_rows(
  top2_up,
  top2_down
)

selected_overall_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val)

############################################################
# FIGURE 16 - TOP 2 UP + TOP 2 DOWN OVERALL LONGITUDINAL DEGs
############################################################


############################################################
# 1. Selected genes
############################################################

selected_genes_overall <- c(
  "PIGR",
  "FCGBP",
  "NRG1",
  "IGFN1"
)


############################################################
# 2. Extract normalised expression from voom object
############################################################

overall_gene_expression <- as.data.frame(
  t(
    vobj_all$E[
      rownames(vobj_all$E) %in% selected_genes_overall,
      ,
      drop = FALSE
    ]
  )
)

overall_gene_expression$Sample_ID <- rownames(overall_gene_expression)


############################################################
# 3. Convert to long format
############################################################

overall_gene_plot_df <- overall_gene_expression %>%
  
  pivot_longer(
    cols = all_of(selected_genes_overall),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  
  left_join(
    long_meta %>%
      dplyr::select(
        Sample_ID,
        DefinitiveTreatment
      ),
    by = "Sample_ID"
  ) %>%
  
  filter(
    !is.na(DefinitiveTreatment),
    !is.na(Expression)
  )


############################################################
# 4. Clean treatment labels and ordering
############################################################

overall_gene_plot_df <- overall_gene_plot_df %>%
  
  mutate(
    
    Treatment = factor(
      DefinitiveTreatment,
      levels = c(
        "Untreated",
        "ShortTreatment",
        "EndPointTreatment"
      ),
      labels = c(
        "Untreated",
        "Short treatment",
        "Endpoint"
      )
    ),
    
    # Keep genes in desired 2 x 2 order
    Gene = factor(
      Gene,
      levels = c(
        "PIGR",
        "FCGBP",
        "NRG1",
        "IGFN1"
      )
    )
  )


############################################################
# 5. Check sample numbers
############################################################

overall_gene_plot_df %>%
  dplyr::count(Gene, Treatment)


############################################################
# 6. Create figure
############################################################

p_gene_overall <- ggplot(
  overall_gene_plot_df,
  aes(
    x = Treatment,
    y = Expression,
    fill = Treatment
  )
) +
  
  geom_boxplot(
    width = 0.55,
    alpha = 0.55,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  
  geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.65,
    colour = "black"
  ) +
  
  facet_wrap(
    ~ Gene,
    ncol = 2,
    scales = "free_y",
    labeller = label_parsed
  ) +
  
  scale_fill_manual(
    values = c(
      "Untreated" = "#9DBAF2",
      "Short treatment" = "#A9D7BF",
      "Endpoint" = "#D1AED8"
    )
  ) +
  
  labs(
    title = "Longitudinal Expression of Top Differentially Expressed Genes",
    x = "Treatment stage",
    y = "Expression (TMM-normalised log2-CPM)",
    fill = "Treatment stage"
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    plot.title = element_text(
      size = 19,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    
    axis.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 12,
      colour = "black"
    ),
    
    strip.background = element_rect(
      fill = "grey95",
      colour = "black",
      linewidth = 0.7
    ),
    
    strip.text = element_text(
      size = 15,
      face = "bold.italic"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    panel.spacing = unit(1.2, "lines")
  )


############################################################
# 7. Display
############################################################

p_gene_overall


############################################################
# 8. Save PNG
############################################################

ggsave(
  file.path(output_dir, "Overall_Longitudinal_Top4_DEGs.png"),
  p_gene_overall,
  width = 11,
  height = 9,
  dpi = 300
)
