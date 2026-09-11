# HNSCC PDX baseline differential expression and pathway analysis
# Adam McClatchey
#
# The full DESeq2 and GSEA pipeline is the representative
# average 3-6-week R vs NR analysis.

rm(list = ls(all.names = TRUE))
gc()

############################################################
# Packages
############################################################

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(apeglm)
library(ggVennDiagram)
library(ComplexUpset)
library(grid)

############################################################
# Output folders
############################################################

root_output <- "Results_Baseline_DE_Pathway_Clean"
representative_output_dir <- file.path(
  root_output,
  "Average_3-6wk_RvsNR"
)
deg_comparison_output_dir <- file.path(
  root_output,
  "Cross_Classification_DEG"
)
pathway_output_dir <- file.path(
  root_output,
  "Cross_Classification_Pathways"
)

dir.create(representative_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_comparison_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pathway_output_dir, recursive = TRUE, showWarnings = FALSE)

# Retain the original output/prefix pattern inside the core pipeline.
output <- paste0(representative_output_dir, "/")
prefix <- "3-6wkAvg_RvsNR"

message(
  "All outputs will be saved under: ",
  normalizePath(root_output, mustWork = FALSE)
)

############################################################
# REPRESENTATIVE AVERAGE 3-6-WEEK R VS NR ANALYSIS
############################################################

############################################################
# Load raw count matrix and metadata
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

# Remove the low-depth sample exactly as in the original script.
outlier_sample <- "HNC0157PRX0A"

meta_avg <- meta_avg %>%
  filter(Sample_ID != outlier_sample)

counts_avg <- counts_avg[
  ,
  colnames(counts_avg) != outlier_sample
]

stopifnot(all(colnames(counts_avg) == meta_avg$Sample_ID))

# Make the sample identifiers explicit for DESeq2. This does not alter the
# samples or model; it guarantees that count columns and metadata rows align.
rownames(meta_avg) <- meta_avg$Sample_ID

############################################################
# DESeq2 model and low-count filtering
############################################################

meta_avg$Definitive_Response_3_6wk <- factor(
  meta_avg$Definitive_Response_3_6wk
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_avg,
  colData = meta_avg,
  design = ~ Definitive_Response_3_6wk
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Positive log2FC indicates higher expression in responders.
res <- results(dds)
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df <- res_df[, c(
  "Gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

write.csv(
  res_df,
  paste0(output, prefix, "_DESeq2_unshrunken_results.csv"),
  row.names = FALSE
)

############################################################
# Variance-stabilising transformation used for Figure 5
############################################################

vsd <- vst(dds, blind = TRUE)

############################################################
# Log2 fold-change shrinkage
############################################################

res_shrunk <- lfcShrink(
  dds,
  coef = "Definitive_Response_3_6wk_Responder_vs_NonResponder",
  type = "apeglm"
)

res_shrunk_df <- as.data.frame(res_shrunk)
res_shrunk_df$Gene <- rownames(res_shrunk_df)
res_shrunk_df <- res_shrunk_df[, c(
  "Gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj"
)]

write.csv(
  res_shrunk_df,
  paste0(output, prefix, "_DESeq2_shrunk_results.csv"),
  row.names = FALSE
)

############################################################
# DEG tables
############################################################

padj_threshold <- 0.05
log2fc_threshold <- 0.585

all_degs <- res_shrunk_df %>%
  filter(
    !is.na(padj),
    padj < padj_threshold,
    abs(log2FoldChange) > log2fc_threshold
  ) %>%
  mutate(
    Direction = case_when(
      log2FoldChange > 0 ~ "Up in Responder",
      log2FoldChange < 0 ~ "Up in NonResponder"
    )
  ) %>%
  arrange(desc(abs(log2FoldChange)))

degs_up <- all_degs %>%
  filter(log2FoldChange > 0)

degs_down <- all_degs %>%
  filter(log2FoldChange < 0)


write.csv(
  all_degs,
  paste0(output, prefix, "_DEGs.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  degs_up,
  paste0(output, prefix, "_DEGs_UpInResponder.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  degs_down,
  paste0(output, prefix, "_DEGs_UpInNonResponder.csv"),
  row.names = FALSE,
  quote = FALSE
)

top20_up_padj <- degs_up %>%
  arrange(padj) %>%
  slice_head(n = 20)

top20_down_padj <- degs_down %>%
  arrange(padj) %>%
  slice_head(n = 20)

top20_up_fc <- degs_up %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 20)

top20_down_fc <- degs_down %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 20)

write.csv(
  top20_up_padj,
  paste0(output, prefix, "_Top20_UpInResponder_padj.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  top20_down_padj,
  paste0(output, prefix, "_Top20_UpInNonResponder_padj.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  top20_up_fc,
  paste0(output, prefix, "_Top20_UpInResponder_FC.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  top20_down_fc,
  paste0(output, prefix, "_Top20_UpInNonResponder_FC.csv"),
  row.names = FALSE,
  quote = FALSE
)

top_degs_fc <- all_degs %>%
  slice_head(n = 20)

top_degs_padj <- all_degs %>%
  arrange(padj) %>%
  slice_head(n = 20)

write.csv(
  top_degs_fc,
  paste0(output, prefix, "_Top20_DEGs_FC.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  top_degs_padj,
  paste0(output, prefix, "_Top20_DEGs_padj.csv"),
  row.names = FALSE,
  quote = FALSE
)

############################################################
# Figure 4 - MA plot
############################################################

ma_padj_thresh <- 0.05
ma_fc_thresh <- 0.585

ma_data <- res_shrunk_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange), baseMean > 0) %>%
  mutate(
    status = case_when(
      padj < ma_padj_thresh & log2FoldChange > ma_fc_thresh ~ "UP",
      padj < ma_padj_thresh & log2FoldChange < -ma_fc_thresh ~ "DOWN",
      TRUE ~ "NS"
    )
  )

label_genes <- bind_rows(
  ma_data %>%
    filter(status == "UP") %>%
    arrange(padj) %>%
    slice_head(n = 10),
  ma_data %>%
    filter(status == "DOWN") %>%
    arrange(padj) %>%
    slice_head(n = 10)
)

n_up <- sum(ma_data$status == "UP")
n_down <- sum(ma_data$status == "DOWN")
n_ns <- sum(ma_data$status == "NS")


p_ma <- ggplot(
  ma_data,
  aes(
    x = log10(baseMean),
    y = log2FoldChange,
    color = status
  )
) +
  geom_point(
    data = filter(ma_data, status == "NS"),
    alpha = 0.3,
    size = 0.8
  ) +
  geom_point(
    data = filter(ma_data, status != "NS"),
    alpha = 0.8,
    size = 1.2
  ) +
  scale_color_manual(
    values = c(
      "NS" = "grey60",
      "UP" = "#D55E00",
      "DOWN" = "#0072B2"
    ),
    labels = c(
      "NS" = paste0("Not significant (n=", n_ns, ")"),
      "UP" = paste0("Up in Responder (n=", n_up, ")"),
      "DOWN" = paste0("Up in NonResponder (n=", n_down, ")")
    ),
    name = NULL
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(
    yintercept = ma_fc_thresh,
    color = "#D55E00",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -ma_fc_thresh,
    color = "#0072B2",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 2.8,
    max.overlaps = 20,
    segment.color = "grey50",
    segment.size = 0.3,
    box.padding = 0.3,
    show.legend = FALSE
  ) +
  labs(
    title = "MA Plot — 3-6 week Avg Responder vs NonResponder",
    subtitle = paste0(
      "Shrunken log2FC | padj < ", ma_padj_thresh,
      " and |log2FC| > ", ma_fc_thresh,
      " | Top 10 up/down labelled"
    ),
    x = "Mean expression [log10(baseMean)]",
    y = "log2 Fold Change (Responder / NonResponder)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
  )

ggsave(
  paste0(output, prefix, "_MAPlot_Shrink.png"),
  p_ma,
  width = 10,
  height = 7,
  dpi = 500
)

############################################################
# Figure 3 - volcano plot
############################################################

padj_thresh <- 0.05
fc_thresh <- 0.585

volcano_data <- as.data.frame(res_shrunk) %>%
  filter(!is.na(padj)) %>%
  mutate(
    diffexpressed = case_when(
      padj < padj_thresh & log2FoldChange > fc_thresh ~ "UP",
      padj < padj_thresh & log2FoldChange < -fc_thresh ~ "DOWN",
      TRUE ~ "NO"
    ),
    Gene = rownames(.),
    delabel = ifelse(diffexpressed != "NO", rownames(.), NA)
  )

total_up <- sum(volcano_data$diffexpressed == "UP")
total_down <- sum(volcano_data$diffexpressed == "DOWN")

top_up <- volcano_data %>%
  filter(diffexpressed == "UP") %>%
  arrange(padj, desc(log2FoldChange)) %>%
  slice_head(n = 20)

top_down <- volcano_data %>%
  filter(diffexpressed == "DOWN") %>%
  arrange(padj, log2FoldChange) %>%
  slice_head(n = 20)

top_genes <- rbind(top_up, top_down)
ymax <- max(-log10(volcano_data$padj), na.rm = TRUE) + 1

p_volcano <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = diffexpressed
  )
) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(
    data = top_genes,
    aes(label = delabel),
    max.overlaps = 12,
    size = 4
  ) +
  scale_color_manual(
    values = c(
      "DOWN" = "blue",
      "NO" = "black",
      "UP" = "red"
    )
  ) +
  geom_hline(
    yintercept = -log10(padj_thresh),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(
    xintercept = c(-fc_thresh, fc_thresh),
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    "text",
    x = Inf,
    y = ymax,
    label = paste("UP:", total_up),
    color = "red",
    hjust = 1.1
  ) +
  annotate(
    "text",
    x = -Inf,
    y = ymax,
    label = paste("DOWN:", total_down),
    color = "blue",
    hjust = -0.1
  ) +
  labs(
    title = "Volcano Plot — 3-6 week Avg Responder vs NonResponder",
    subtitle = "Shrunken log2FC | padj < 0.05 | |log2FC| > 0.585",
    x = "log2 Fold Change (Responder / NonResponder)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  paste0(output, prefix, "_VolcanoPlot_Shrink.png"),
  p_volcano,
  width = 11,
  height = 11,
  dpi = 500
)

############################################################
# Figure 5 - heatmap of all significant DEGs
############################################################

norm_counts <- counts(dds, normalized = TRUE)

write.csv(
  norm_counts,
  paste0(output, prefix, "_Normalized_counts.csv"),
  quote = FALSE
)

vsd_mat <- assay(vsd)

sample_annotation <- as.data.frame(colData(dds)) %>%
  dplyr::select(Definitive_Response_3_6wk)

colnames(sample_annotation) <- "Response"
rownames(sample_annotation) <- colnames(vsd_mat)

ann_colors <- list(
  Response = c(
    "Responder" = "#009E73",
    "NonResponder" = "#D55E00"
  )
)

heatmap_colours <- colorRampPalette(
  c("navy", "white", "firebrick")
)(100)

heatmap_breaks <- seq(-3, 3, length.out = 101)

all_deg_genes <- all_degs$Gene
all_deg_mat <- vsd_mat[rownames(vsd_mat) %in% all_deg_genes, ]

all_deg_z <- t(scale(t(all_deg_mat)))
all_deg_z[is.na(all_deg_z)] <- 0

sample_order <- order(sample_annotation$Response)
all_deg_z <- all_deg_z[, sample_order]
sample_annotation_ordered <- sample_annotation[
  sample_order,
  ,
  drop = FALSE
]

pdf(
  paste0(output, prefix, "_Heatmap_All_DEGs.pdf"),
  width = 10,
  height = 10
)

pheatmap(
  all_deg_z,
  scale = "none",
  color = heatmap_colours,
  breaks = heatmap_breaks,
  annotation_col = sample_annotation_ordered,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  main = "All DEGs — 3-6 week Avg Responder vs NonResponder"
)

dev.off()

############################################################
# GSEA ranked list
############################################################

DESeq_results <- as.data.frame(res)
DESeq_results <- DESeq_results[!is.na(DESeq_results$padj), ]

geneList <- DESeq_results$stat
names(geneList) <- rownames(DESeq_results)
geneList <- sort(geneList, decreasing = TRUE)

############################################################
# HALLMARK GSEA
############################################################

h_t2g <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
) %>%
  dplyr::select(gs_name, gene_symbol)

h_gseaResult <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = h_t2g,
  verbose = TRUE,
  seed = TRUE
)

h_df <- as.data.frame(h_gseaResult)

write.csv(
  h_df,
  paste0(output, prefix, "_Hallmarks_GSEA_RESULTS.csv"),
  row.names = FALSE
)

h_significant <- h_df[h_df$p.adjust < 0.05, ]
h_sorted <- h_significant[order(-h_significant$NES), ]


write.csv(
  h_sorted,
  paste0(output, prefix, "_Hallmarks_GSEA_Pathways_Sorted.csv"),
  row.names = FALSE
)

############################################################
# KEGG GSEA AND FIGURE 10
############################################################

kegg_t2g <- msigdbr(
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:KEGG_LEGACY"
) %>%
  dplyr::select(gs_name, gene_symbol)

kegg_gseaResult <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = kegg_t2g,
  verbose = TRUE,
  seed = TRUE
)

kegg_df <- as.data.frame(kegg_gseaResult)

write.csv(
  kegg_df,
  paste0(output, prefix, "_KEGG_GSEA_RESULTS.csv"),
  row.names = FALSE
)

kegg_significant <- kegg_df[kegg_df$p.adjust < 0.05, ]
kegg_sorted <- kegg_significant[order(-kegg_significant$NES), ]


write.csv(
  kegg_sorted,
  paste0(output, prefix, "_KEGG_GSEA_Pathways_Sorted.csv"),
  row.names = FALSE
)

kegg_significant <- kegg_significant %>%
  mutate(
    Pathway = Description %>%
      stringr::str_remove("^KEGG_") %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_to_title() %>%
      stringr::str_replace_all("\\bJak\\b", "JAK") %>%
      stringr::str_replace_all("\\bStat\\b", "STAT") %>%
      stringr::str_replace_all("\\bMapk\\b", "MAPK") %>%
      stringr::str_replace_all("\\bVegf\\b", "VEGF") %>%
      stringr::str_replace_all("\\bErbb\\b", "ERBB") %>%
      stringr::str_replace_all("\\bFc\\b", "Fc") %>%
      stringr::str_replace_all("\\bNk\\b", "NK") %>%
      stringr::str_replace_all("\\bDna\\b", "DNA") %>%
      stringr::str_replace_all("\\bRna\\b", "RNA") %>%
      stringr::str_replace_all("\\bEcm\\b", "ECM") %>%
      stringr::str_replace_all("\\bP450\\b", "P450") %>%
      stringr::str_replace_all("\\bH Pylori\\b", "H. pylori") %>%
      stringr::str_replace_all("\\bE Coli\\b", "E. coli")
  )

p_kegg_sig <- ggplot(
  kegg_significant,
  aes(reorder(Pathway, NES), NES)
) +
  geom_col(fill = "#29BF4E", width = 0.7) +
  coord_flip() +
  labs(
    x = "Pathway",
    y = "Normalized Enrichment Score",
    title = paste0(
      "Representative KEGG pathway enrichment\n",
      "(Average responder vs non-responder comparison)"
    )
  ) +
  theme(
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, colour = "black", face = "bold"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 14, colour = "black", face = "bold"),
    axis.text.x = element_text(size = 12, colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  paste0(output, prefix, "_KEGG_GSEA_barplot_Significant.png"),
  p_kegg_sig,
  width = 11,
  height = 8.5,
  dpi = 500
)

############################################################
# GO GSEA AND FIGURE 11
############################################################

c5_t2g <- msigdbr(
  species = "Homo sapiens",
  collection = "C5"
) %>%
  dplyr::select(gs_name, gene_symbol)

c5_gseaResult <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  TERM2GENE = c5_t2g,
  verbose = TRUE,
  seed = TRUE
)

c5_df <- as.data.frame(c5_gseaResult)

write.csv(
  c5_df,
  paste0(output, prefix, "_GO_GSEA_RESULTS.csv"),
  row.names = FALSE
)

c5_df <- c5_df[order(-abs(c5_df$NES)), ]


write.csv(
  c5_df,
  paste0(output, prefix, "_GO_GSEA_RESULTS_Sorted.csv"),
  row.names = FALSE
)

top_go_positive <- c5_df %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

top_go_negative <- c5_df %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  slice_head(n = 10)

top_go_plot <- bind_rows(
  top_go_positive,
  top_go_negative
) %>%
  mutate(
    Description = Description %>%
      str_remove("^GOBP_") %>%
      str_remove("^GOCC_") %>%
      str_remove("^GOMF_") %>%
      str_remove("^HP_") %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_replace_all("\\bDna\\b", "DNA") %>%
      str_replace_all("\\bRna\\b", "RNA") %>%
      str_replace_all("\\bEr\\b", "ER") %>%
      str_replace_all("\\bGolgi\\b", "Golgi")
  )

write.csv(
  top_go_plot,
  paste0(output, prefix, "_GO_Top10_Positive_and_Negative.csv"),
  row.names = FALSE
)

p_go <- ggplot(
  top_go_plot,
  aes(reorder(Description, NES), NES)
) +
  geom_col(fill = "#29BF4E", width = 0.7) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score",
    title = paste0(
      "Representative Functional Term Enrichment\n",
      "(Average responder vs non-responder comparison)"
    )
  ) +
  theme(
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 15, face = "bold", colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank()
  )

ggsave(
  filename = paste0(output, prefix, "_GO_GSEA_barplot.png"),
  plot = p_go,
  width = 12,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = paste0(output, prefix, "_GO_GSEA_barplot.pdf"),
  plot = p_go,
  width = 12,
  height = 8.5,
  device = cairo_pdf,
  bg = "white"
)

############################################################
# C6 ONCOGENIC SIGNATURE GSEA
############################################################

c6_t2g <- msigdbr(
  species = "Homo sapiens",
  collection = "C6"
) %>%
  dplyr::select(gs_name, gene_symbol)

c6_gseaResult <- GSEA(
  geneList,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = c6_t2g,
  verbose = TRUE,
  seed = TRUE
)

c6_df <- as.data.frame(c6_gseaResult)

write.csv(
  c6_df,
  paste0(output, prefix, "_C6_GSEA_RESULTS.csv"),
  row.names = FALSE
)

c6_significant <- c6_df[c6_df$p.adjust < 0.05, ]
c6_sorted <- c6_significant[order(-c6_significant$NES), ]


write.csv(
  c6_sorted,
  paste0(output, prefix, "_C6_GSEA_Pathways_Sorted.csv"),
  row.names = FALSE
)

############################################################
# CROSS-CLASSIFICATION DEG ANALYSIS
############################################################

# The newly generated average R vs NR DEG table is used here.
# The other five tables are read-only inputs from the repeated
# applications of the same DE pipeline.

rnr_3wk <- read.csv(
  "Results_3wk_RvsNR/3wk_RvsNR_DEGs.csv"
)

rnr_6wk <- read.csv(
  "Results_6wk_RvsNR/6wk_RvsNR_DEGs.csv"
)

rnr_avg <- read.csv(
  paste0(output, prefix, "_DEGs.csv")
)

prpd_3wk <- read.csv(
  "Results_3wk_PRvsPD/3wk_PRvsPD_DEGs.csv"
)

prpd_6wk <- read.csv(
  "Results_6wk_PRvsPD/6wk_PRvsPD_DEGs.csv"
)

prpd_avg <- read.csv(
  "Results_Avg_PRvsPD/3-6wkAvg_PRvsPD_DEGs.csv"
)

deg_tables <- list(
  RNR_3wk = rnr_3wk,
  RNR_6wk = rnr_6wk,
  RNR_Avg = rnr_avg,
  PRPD_3wk = prpd_3wk,
  PRPD_6wk = prpd_6wk,
  PRPD_Avg = prpd_avg
)

############################################################
# Table 3 - DEGs across the six response classifications
############################################################

table3 <- tibble(
  Comparison = c(
    "3wk R vs NR",
    "6wk R vs NR",
    "Avg 3-6 wk R vs NR",
    "3wk PR vs PD",
    "6wk PR vs PD",
    "Avg 3-6 wk PR vs PD"
  ),
  Group_1 = c("39 R", "34 R", "40 R", "4 PR", "7 PR", "4 PR"),
  Group_2 = c("15 NR", "9 NR", "14 NR", "15 PD", "9 PD", "14 PD"),
  Total_Samples = c(54L, 43L, 54L, 19L, 16L, 18L),
  Significant_DEGs = vapply(deg_tables, nrow, integer(1)),
  Higher_in_Group_1 = vapply(
    deg_tables,
    function(x) sum(x$log2FoldChange > 0),
    integer(1)
  ),
  Higher_in_Group_2 = vapply(
    deg_tables,
    function(x) sum(x$log2FoldChange < 0),
    integer(1)
  )
)


write.csv(
  table3,
  file.path(deg_comparison_output_dir, "Table3_DEG_Counts.csv"),
  row.names = FALSE
)

############################################################
# Gene sets used in Figures 6-7 and Table 4
############################################################

genes_rnr_3wk <- unique(rnr_3wk$Gene)
genes_rnr_6wk <- unique(rnr_6wk$Gene)
genes_rnr_avg <- unique(rnr_avg$Gene)

genes_prpd_3wk <- unique(prpd_3wk$Gene)
genes_prpd_6wk <- unique(prpd_6wk$Gene)
genes_prpd_avg <- unique(prpd_avg$Gene)

common_rnr <- Reduce(
  intersect,
  list(
    genes_rnr_3wk,
    genes_rnr_6wk,
    genes_rnr_avg
  )
)

common_prpd <- Reduce(
  intersect,
  list(
    genes_prpd_3wk,
    genes_prpd_6wk,
    genes_prpd_avg
  )
)

common_all <- Reduce(
  intersect,
  list(
    genes_rnr_3wk,
    genes_rnr_6wk,
    genes_rnr_avg,
    genes_prpd_3wk,
    genes_prpd_6wk,
    genes_prpd_avg
  )
)


write.csv(
  data.frame(Gene = common_rnr),
  file.path(deg_comparison_output_dir, "Common_Genes_RvsNR.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Gene = common_prpd),
  file.path(deg_comparison_output_dir, "Common_Genes_PRvsPD.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Gene = common_all),
  file.path(deg_comparison_output_dir, "Common_Genes_All6.csv"),
  row.names = FALSE
)

all_gene_lists_internal <- list(
  RNR_3wk = genes_rnr_3wk,
  RNR_6wk = genes_rnr_6wk,
  RNR_Avg = genes_rnr_avg,
  PRPD_3wk = genes_prpd_3wk,
  PRPD_6wk = genes_prpd_6wk,
  PRPD_Avg = genes_prpd_avg
)

all_genes <- unique(unlist(all_gene_lists_internal))

gene_frequency <- data.frame(
  Gene = all_genes,
  Analyses_Present = vapply(
    all_genes,
    function(g) {
      sum(vapply(
        all_gene_lists_internal,
        function(x) g %in% x,
        logical(1)
      ))
    },
    integer(1)
  )
) %>%
  arrange(desc(Analyses_Present))

gene_recurrence_counts <- gene_frequency %>%
  dplyr::count(Analyses_Present) %>%
  tidyr::complete(
    Analyses_Present = 1:6,
    fill = list(n = 0)
  ) %>%
  arrange(Analyses_Present)

write.csv(
  gene_frequency,
  file.path(deg_comparison_output_dir, "Consensus_Gene_Frequency.csv"),
  row.names = FALSE
)

write.csv(
  gene_recurrence_counts,
  file.path(deg_comparison_output_dir, "DEG_Recurrence_Counts.csv"),
  row.names = FALSE
)

############################################################
# Figure 6 - R vs NR and PR vs PD Venn diagrams
############################################################

venn_rnr <- list(
  `3 Week` = genes_rnr_3wk,
  `6 Week` = genes_rnr_6wk,
  `3-6 Week Avg` = genes_rnr_avg
)

p_rnr <- ggVennDiagram(
  venn_rnr,
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#D6ECFF",
    high = "#4DA6FF"
  ) +
  ggtitle("DEG Overlap: Responder vs NonResponder")

ggsave(
  file.path(deg_comparison_output_dir, "Venn_RvsNR.png"),
  p_rnr,
  width = 8,
  height = 7,
  dpi = 500
)

venn_prpd <- list(
  `3 Week` = genes_prpd_3wk,
  `6 Week` = genes_prpd_6wk,
  `3-6 Week Avg` = genes_prpd_avg
)

p_prpd <- ggVennDiagram(
  venn_prpd,
  label_alpha = 0
) +
  scale_fill_gradient(
    low = "#D6ECFF",
    high = "#4DA6FF"
  ) +
  ggtitle("DEG Overlap: PR vs PD")

ggsave(
  file.path(deg_comparison_output_dir, "Venn_PRvsPD.png"),
  p_prpd,
  width = 8,
  height = 7,
  dpi = 500
)

############################################################
# Figure 7 - DEG UpSet plot across all six analyses
############################################################

all_gene_lists <- list(
  `RvsNR 3 Week` = genes_rnr_3wk,
  `RvsNR 6 Week` = genes_rnr_6wk,
  `RvsNR Avg` = genes_rnr_avg,
  `PRvsPD 3 Week` = genes_prpd_3wk,
  `PRvsPD 6 Week` = genes_prpd_6wk,
  `PRvsPD Avg` = genes_prpd_avg
)

all_genes <- unique(unlist(all_gene_lists))
upset_df <- data.frame(Gene = all_genes)

for (set_name in names(all_gene_lists)) {
  upset_df[[set_name]] <- as.integer(
    all_genes %in% all_gene_lists[[set_name]]
  )
}

intersection_summary <- upset_df %>%
  group_by(across(-Gene)) %>%
  summarise(Intersection_Size = n(), .groups = "drop") %>%
  arrange(desc(Intersection_Size))

write.csv(
  intersection_summary,
  file.path(deg_comparison_output_dir, "DEG_Intersection_Summary_All6.csv"),
  row.names = FALSE
)

p_upset <- upset(
  upset_df,
  intersect = names(all_gene_lists),
  name = "DEGs",
  min_size = 20,
  width_ratio = 0.25
) +
  ggplot2::labs(
    title = "Overlap of Differentially Expressed Genes Across Six Analyses",
    subtitle = paste0(
      "Genes significant in 3-week, 6-week and average ",
      "response classifications"
    )
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)
  )

ggsave(
  file.path(deg_comparison_output_dir, "DEG_UpSet_All6.png"),
  p_upset,
  width = 14,
  height = 8,
  dpi = 600
)

############################################################
# Consensus DEG summary and Table 4
############################################################

consensus_genes <- gene_frequency %>%
  filter(Analyses_Present >= 4)

consensus_summary <- lapply(
  consensus_genes$Gene,
  function(gene) {
    fc_values <- c()

    for (tbl in deg_tables) {
      hit <- tbl[tbl$Gene == gene, ]

      if (nrow(hit) > 0) {
        fc_values <- c(fc_values, hit$log2FoldChange[1])
      }
    }

    data.frame(
      Gene = gene,
      Analyses_Present = length(fc_values),
      Mean_log2FC = mean(fc_values),
      Mean_abs_log2FC = mean(abs(fc_values)),
      Positive_Count = sum(fc_values > 0),
      Negative_Count = sum(fc_values < 0)
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    Direction = case_when(
      Positive_Count == Analyses_Present ~ "Up in favourable response",
      Negative_Count == Analyses_Present ~ "Up in poor response",
      TRUE ~ "Mixed direction"
    )
  ) %>%
  arrange(
    desc(Analyses_Present),
    desc(Mean_abs_log2FC)
  )

write.csv(
  consensus_summary,
  file.path(
    deg_comparison_output_dir,
    "Consensus_Gene_Summary_with_Direction.csv"
  ),
  row.names = FALSE
)

table4 <- consensus_summary %>%
  filter(Analyses_Present >= 5) %>%
  arrange(desc(Mean_abs_log2FC)) %>%
  slice_head(n = 10) %>%
  dplyr::select(
    Gene,
    Analyses_Present,
    Mean_log2FC,
    Direction
  )


write.csv(
  table4,
  file.path(deg_comparison_output_dir, "Table4_Top10_Consensus_DEGs.csv"),
  row.names = FALSE
)

############################################################
# CROSS-CLASSIFICATION PATHWAY ANALYSIS
############################################################

analysis_folders <- c(
  "3wk_RvsNR"  = "Results_3wk_RvsNR",
  "6wk_RvsNR"  = "Results_6wk_RvsNR",
  "Avg_RvsNR"  = representative_output_dir,
  "3wk_PRvsPD" = "Results_3wk_PRvsPD",
  "6wk_PRvsPD" = "Results_6wk_PRvsPD",
  "Avg_PRvsPD" = "Results_Avg_PRvsPD"
)

analysis_prefixes <- c(
  "3wk_RvsNR",
  "6wk_RvsNR",
  "3-6wkAvg_RvsNR",
  "3wk_PRvsPD",
  "6wk_PRvsPD",
  "3-6wkAvg_PRvsPD"
)

hallmark <- map2(
  analysis_folders,
  analysis_prefixes,
  ~ read_csv(
    file.path(.x, paste0(.y, "_Hallmarks_GSEA_Pathways_Sorted.csv")),
    show_col_types = FALSE
  )
)

kegg <- map2(
  analysis_folders,
  analysis_prefixes,
  ~ read_csv(
    file.path(.x, paste0(.y, "_KEGG_GSEA_Pathways_Sorted.csv")),
    show_col_types = FALSE
  )
)

go <- map2(
  analysis_folders,
  analysis_prefixes,
  ~ read_csv(
    file.path(.x, paste0(.y, "_GO_GSEA_RESULTS_Sorted.csv")),
    show_col_types = FALSE
  )
)

oncogenic <- map2(
  analysis_folders,
  analysis_prefixes,
  ~ read_csv(
    file.path(.x, paste0(.y, "_C6_GSEA_Pathways_Sorted.csv")),
    show_col_types = FALSE
  )
)

names(hallmark) <- names(analysis_folders)
names(kegg) <- names(analysis_folders)
names(go) <- names(analysis_folders)
names(oncogenic) <- names(analysis_folders)

############################################################
# Per-analysis pathway counts used in Results 3.3
############################################################

summarise_gsea <- function(gsea_list) {
  bind_rows(
    lapply(names(gsea_list), function(x) {
      df <- gsea_list[[x]]

      tibble(
        Analysis = x,
        Positive = sum(df$NES > 0),
        Negative = sum(df$NES < 0),
        Total = nrow(df),
        Max_Positive_NES = max(df$NES),
        Max_Negative_NES = min(df$NES)
      )
    })
  )
}

hallmark_summary <- summarise_gsea(hallmark)
kegg_summary <- summarise_gsea(kegg)
go_summary <- summarise_gsea(go)
onc_summary <- summarise_gsea(oncogenic)

write.csv(
  hallmark_summary,
  file.path(pathway_output_dir, "Hallmark_Per_Analysis_Counts.csv"),
  row.names = FALSE
)
write.csv(
  kegg_summary,
  file.path(pathway_output_dir, "KEGG_Per_Analysis_Counts.csv"),
  row.names = FALSE
)
write.csv(
  go_summary,
  file.path(pathway_output_dir, "GO_Per_Analysis_Counts.csv"),
  row.names = FALSE
)
write.csv(
  onc_summary,
  file.path(pathway_output_dir, "C6_Per_Analysis_Counts.csv"),
  row.names = FALSE
)


summarise_recurrence <- function(gsea_list) {
  bind_rows(
    lapply(names(gsea_list), function(x) {
      gsea_list[[x]] %>% mutate(Analysis = x)
    })
  ) %>%
    group_by(Description) %>%
    summarise(
      Analyses_Present = n(),
      Analyses = paste(sort(Analysis), collapse = ", "),
      Positive = sum(NES > 0),
      Negative = sum(NES < 0),
      Mean_NES = mean(NES),
      Mean_Abs_NES = mean(abs(NES)),
      Max_NES = max(NES),
      Min_NES = min(NES),
      Best_padj = min(p.adjust),
      .groups = "drop"
    ) %>%
    arrange(desc(Analyses_Present), desc(Mean_Abs_NES))
}

recurrence_counts <- function(recurrence_table) {
  recurrence_table %>%
    dplyr::count(Analyses_Present) %>%
    tidyr::complete(Analyses_Present = 1:6, fill = list(n = 0)) %>%
    arrange(desc(Analyses_Present))
}

############################################################
# Hallmark recurrence and Figure 8
############################################################

hallmark_all <- bind_rows(
  lapply(names(hallmark), function(x) hallmark[[x]] %>% mutate(Analysis = x))
)
hallmark_recurrence <- summarise_recurrence(hallmark)
hallmark_counts <- recurrence_counts(hallmark_recurrence)

write.csv(
  hallmark_recurrence,
  file.path(pathway_output_dir, "Hallmark_Pathway_Recurrence.csv"),
  row.names = FALSE
)
write.csv(
  hallmark_counts,
  file.path(pathway_output_dir, "Hallmark_Recurrence_Counts.csv"),
  row.names = FALSE
)


analysis_order <- c(
  "3wk_RvsNR",
  "6wk_RvsNR",
  "Avg_RvsNR",
  "3wk_PRvsPD",
  "6wk_PRvsPD",
  "Avg_PRvsPD"
)

analysis_labels <- c(
  "3wk_RvsNR" = "3wk",
  "6wk_RvsNR" = "6wk",
  "Avg_RvsNR" = "Avg",
  "3wk_PRvsPD" = "3wk",
  "6wk_PRvsPD" = "6wk",
  "Avg_PRvsPD" = "Avg"
)

selected_hallmarks <- hallmark_recurrence %>%
  filter(Analyses_Present >= 3) %>%
  arrange(desc(Analyses_Present), desc(Mean_Abs_NES)) %>%
  mutate(
    Pathway = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title() %>%
      str_replace_all("\\bTgf Beta\\b", "TGF-β") %>%
      str_replace_all("\\bTnfa\\b", "TNFα") %>%
      str_replace_all("\\bNfkb\\b", "NF-κB") %>%
      str_replace_all("\\bIl6\\b", "IL6") %>%
      str_replace_all("\\bIl2\\b", "IL2") %>%
      str_replace_all("\\bJak\\b", "JAK") %>%
      str_replace_all("\\bStat3\\b", "STAT3") %>%
      str_replace_all("\\bStat5\\b", "STAT5") %>%
      str_replace_all("\\bPi3k\\b", "PI3K") %>%
      str_replace_all("\\bAkt\\b", "AKT") %>%
      str_replace_all("\\bMtor\\b", "mTOR") %>%
      str_replace_all("\\bKras\\b", "KRAS") %>%
      str_replace_all("\\bP53\\b", "p53") %>%
      str_replace_all("\\bMyc\\b", "MYC") %>%
      str_replace_all("\\bE2f\\b", "E2F") %>%
      str_replace_all("\\bG2m\\b", "G2M") %>%
      str_replace_all("\\bUv\\b", "UV") %>%
      str_replace_all("\\bV1\\b", "V1") %>%
      str_replace_all("\\bV2\\b", "V2") %>%
      str_replace_all("\\bDn\\b", "DN") %>%
      str_replace_all("\\bUp\\b", "UP")
  )


pathway_order <- rev(selected_hallmarks$Pathway)

hallmark_heatmap_df <- hallmark_all %>%
  filter(Description %in% selected_hallmarks$Description) %>%
  dplyr::select(Description, Analysis, NES) %>%
  right_join(
    expand_grid(
      Description = selected_hallmarks$Description,
      Analysis = analysis_order
    ),
    by = c("Description", "Analysis")
  ) %>%
  left_join(
    selected_hallmarks %>% dplyr::select(Description, Pathway),
    by = "Description"
  ) %>%
  mutate(
    Analysis = factor(Analysis, levels = analysis_order),
    Analysis_Label = recode(as.character(Analysis), !!!analysis_labels),
    Analysis_Label = factor(Analysis_Label, levels = c("3wk", "6wk", "Avg")),
    Comparison = case_when(
      Analysis %in% c("3wk_RvsNR", "6wk_RvsNR", "Avg_RvsNR") ~ "R vs NR",
      TRUE ~ "PR vs PD"
    ),
    Comparison = factor(Comparison, levels = c("R vs NR", "PR vs PD")),
    Pathway = factor(Pathway, levels = pathway_order)
  )

hallmark_heatmap_final <- ggplot(
  hallmark_heatmap_df,
  aes(x = Analysis_Label, y = Pathway, fill = NES)
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
    name = "NES",
    palette = "RdBu",
    direction = -1,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    breaks = c(-2, 0, 2),
    na.value = "grey95"
  ) +
  labs(
    title = "Recurrent Hallmark Pathway Enrichment Across Response Classifications",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
    axis.text.y = element_text(size = 9.5, face = "bold", colour = "black"),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey94", colour = "grey70"),
    panel.spacing.x = unit(0.7, "cm"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.6),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  ) +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(4.5, "cm"),
      barwidth = unit(0.5, "cm")
    )
  )

ggsave(
  file.path(pathway_output_dir, "Figure8_Hallmark_Recurrence_Heatmap.png"),
  hallmark_heatmap_final,
  width = 8.8,
  height = 11,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(pathway_output_dir, "Figure8_Hallmark_Recurrence_Heatmap.pdf"),
  hallmark_heatmap_final,
  width = 8.8,
  height = 11,
  device = cairo_pdf,
  bg = "white"
)

############################################################
# KEGG recurrence and Figure 9
############################################################

kegg_recurrence <- summarise_recurrence(kegg)
kegg_counts <- recurrence_counts(kegg_recurrence) %>%
  mutate(
    Recurrence = factor(
      paste0(Analyses_Present, "/6 analyses"),
      levels = paste0(1:6, "/6 analyses")
    )
  )

write.csv(
  kegg_recurrence,
  file.path(pathway_output_dir, "KEGG_Pathway_Recurrence.csv"),
  row.names = FALSE
)
write.csv(
  kegg_counts,
  file.path(pathway_output_dir, "KEGG_Recurrence_Counts.csv"),
  row.names = FALSE
)


kegg_recurrence_plot <- ggplot(
  kegg_counts,
  aes(y = Recurrence, x = n, fill = Analyses_Present)
) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.25, size = 4.2, fontface = "bold") +
  scale_fill_gradient(low = "grey80", high = "#1F4E79", guide = "none") +
  expand_limits(x = max(kegg_counts$n) + 3) +
  labs(
    x = "Number of pathways",
    y = NULL,
    title = "Recurrence of enriched KEGG pathways across response classifications"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.9, size = 12),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(colour = "black", size = 11),
    axis.ticks.y = element_blank()
  )

ggsave(
  file.path(pathway_output_dir, "Figure9_KEGG_Pathway_Recurrence.png"),
  kegg_recurrence_plot,
  width = 6.2,
  height = 4.5,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(pathway_output_dir, "Figure9_KEGG_Pathway_Recurrence.pdf"),
  kegg_recurrence_plot,
  width = 6.2,
  height = 4.5,
  device = cairo_pdf,
  bg = "white"
)

############################################################
# GO recurrence supporting tables
############################################################

go_recurrence <- summarise_recurrence(go)
go_counts <- recurrence_counts(go_recurrence)

write.csv(
  go_recurrence,
  file.path(pathway_output_dir, "GO_Recurrent_Pathways.csv"),
  row.names = FALSE
)
write.csv(
  go_counts,
  file.path(pathway_output_dir, "GO_Recurrence_Counts.csv"),
  row.names = FALSE
)
write.csv(
  go_recurrence %>%
    filter(Analyses_Present == 6, Mean_NES < 0) %>%
    dplyr::select(Description, Mean_NES) %>%
    arrange(Mean_NES),
  file.path(pathway_output_dir, "GO_Negative_in_All6_Analyses.csv"),
  row.names = FALSE
)



############################################################
# Oncogenic-signature recurrence and Figure 12
############################################################

onc_recurrence <- summarise_recurrence(oncogenic)
onc_counts <- recurrence_counts(onc_recurrence) %>%
  mutate(
    Recurrence = factor(
      paste0(Analyses_Present, "/6 analyses"),
      levels = paste0(6:1, "/6 analyses")
    )
  )

write.csv(
  onc_recurrence,
  file.path(pathway_output_dir, "Oncogenic_Recurrent_Signatures.csv"),
  row.names = FALSE
)
write.csv(
  onc_counts,
  file.path(pathway_output_dir, "C6_Recurrence_Counts.csv"),
  row.names = FALSE
)



onc_recurrence_plot <- ggplot(
  onc_counts,
  aes(y = Recurrence, x = n, fill = Analyses_Present)
) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.25, size = 4.2, fontface = "bold") +
  scale_fill_gradient(low = "grey80", high = "#1F4E79", guide = "none") +
  expand_limits(x = max(onc_counts$n) + 4) +
  labs(
    x = "Number of signatures",
    y = NULL,
    title = "Oncogenic signature recurrence across response classifications"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.47, size = 12),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(colour = "black", size = 11),
    axis.ticks.y = element_blank()
  )

ggsave(
  file.path(pathway_output_dir, "Figure12_Oncogenic_Signature_Recurrence.png"),
  onc_recurrence_plot,
  width = 6.2,
  height = 4.5,
  dpi = 600,
  bg = "white"
)
ggsave(
  file.path(pathway_output_dir, "Figure12_Oncogenic_Signature_Recurrence.pdf"),
  onc_recurrence_plot,
  width = 6.2,
  height = 4.5,
  device = cairo_pdf,
  bg = "white"
)

############################################################
# Reproducibility record
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
