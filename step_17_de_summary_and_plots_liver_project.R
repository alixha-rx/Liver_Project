# =============================================================================
# STEP 17: DE summary tables, diagnostics and figure panels  -->  FIGURE 2c, 2d
# =============================================================================
#
# By Nancy R. Zhang
#
# Turns the six DE tables from step 16 into the four summary folders the paper
# figures are assembled from, runs the p-value diagnostics that say whether the
# analysis behaved, and enriches each gene set.
#
# ---------------------------------------------------------------------------
# THE FOUR CATEGORIES
# ---------------------------------------------------------------------------
# Every gene-zone result passing FDR < 0.2 falls into one of:
#     female_up_in_old     female_up_in_young
#     male_up_in_old       male_up_in_young
#
# A gene can appear in more than one zone within a category, so the number of
# ROWS exceeds the number of unique GENES. Both counts are reported, and they
# answer different questions: unique genes is "how many genes respond to age",
# rows is "how many zone-specific effects were detected".
#
# ---------------------------------------------------------------------------
# READ THE P-VALUE HISTOGRAMS FIRST
# ---------------------------------------------------------------------------
# Block 2 writes a raw p-value histogram per sex-zone. This is the quickest
# check that the model is sound, and it should be looked at before any gene
# list is taken seriously:
#
#   flat with a spike near zero   what you want: a uniform null plus real signal
#   flat with no spike            no detectable effect in that zone
#   rising towards 1              the model is misspecified or over-corrected
#   spike at both ends            usually a variance estimation problem
#
# ---------------------------------------------------------------------------
# INPUT   RESULTS_DIR/pipeline_E/DE_*.csv                       (step 16)
#         FIGURES_DIR/gene_heatmaps/Zone_heatmap_*.png          (step 10)
#         CACHE_DIR/expr_pseudobulk.rds                         (step 15)
# OUTPUT  RESULTS_DIR/pipeline_E/summary/<category>/
#             <category>.csv     the significant gene-zone rows
#             <category>.pdf     the same, formatted as a table (Fig 2d)
#             <GENE>_heatmap.png one per gene, copied from step 10 (Fig 2e-h)
#         RESULTS_DIR/pipeline_E/summary/enrichment_plots/       (Fig 2c)
#         FIGURES_DIR/limma_diagnostics/pvalue_histograms.png
#
# RUN TIME: a few minutes; the enrichment queries the KEGG API.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

source("config.R")

library(data.table)
library(ggplot2)
library(gridExtra)

de_dir       <- file.path(RESULTS_DIR, "pipeline_E")
summary_base <- file.path(de_dir, "summary")
heatmaps_dir <- file.path(FIGURES_DIR, "gene_heatmaps")
diag_dir     <- file.path(FIGURES_DIR, "limma_diagnostics")
dir.create(summary_base, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

zones <- ZONES
genders <- c("female", "male")


# --- BLOCK 1: LOAD THE SIX DE TABLES -----------------------------------------
#
# Zone and sex are recovered from the filename, which is why step 16 writes
# them in the fixed form DE_<zone>_<sex>WT_old_vs_young.csv.

de_files <- list.files(de_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
if (length(de_files) == 0) stop("No DE files in ", de_dir, " -- run step 16 first.")

all_de <- list()
for (f in de_files) {
  tt <- fread(f)
  bn <- basename(f)
  parts <- sub("^DE_", "", sub("WT_old_vs_young\\.csv$", "", bn))
  tt$zone   <- sub("_[^_]+$", "", parts)
  tt$gender <- sub("^.*_", "", parts)
  all_de[[bn]] <- tt
}
de <- rbindlist(all_de)
cat("[step 17] loaded", length(de_files), "DE tables,", nrow(de), "gene-zone rows\n")


# --- BLOCK 2: P-VALUE DIAGNOSTICS --------------------------------------------

plots <- list()
for (sex in genders) {
  for (z in zones) {
    sub <- de[de$gender == sex & de$zone == z, ]
    if (nrow(sub) == 0) next
    n_sig <- sum(sub$adj.P.Val < FDR_THRESHOLD)
    plots[[paste(sex, z)]] <- ggplot(sub, aes(x = P.Value)) +
      geom_histogram(breaks = seq(0, 1, 0.05), fill = "steelblue", colour = "white") +
      geom_hline(yintercept = nrow(sub) / 20, linetype = "dashed", colour = "red") +
      labs(title = paste(sex, z), subtitle = paste(n_sig, "at FDR <", FDR_THRESHOLD),
           x = "raw p-value", y = "genes") +
      theme_bw(base_size = 9)
  }
}
png(file.path(diag_dir, "pvalue_histograms.png"), width = 1200, height = 700)
do.call(grid.arrange, c(plots, ncol = 3))
dev.off()
cat("[step 17] p-value histograms ->", diag_dir, "\n")
cat("[step 17] the dashed line is the uniform-null expectation; a spike in the\n")
cat("          leftmost bin above it is real signal.\n")


# --- BLOCK 3: SUMMARY FOLDERS ------------------------------------------------

categories <- list(
  list(gender = "female", direction = "OLD",   label = "female_up_in_old"),
  list(gender = "female", direction = "YOUNG", label = "female_up_in_young"),
  list(gender = "male",   direction = "OLD",   label = "male_up_in_old"),
  list(gender = "male",   direction = "YOUNG", label = "male_up_in_young")
)

category_genes <- list()

for (cat_info in categories) {
  sex <- cat_info$gender
  direction <- cat_info$direction
  label <- cat_info$label

  sub_de <- de[de$gender == sex & de$direction == direction &
               de$adj.P.Val < FDR_THRESHOLD, ]

  if (nrow(sub_de) == 0) {
    cat("[step 17]", label, ": no genes at FDR <", FDR_THRESHOLD, "\n")
    next
  }

  sig_genes <- unique(sub_de$Gene)
  category_genes[[label]] <- sig_genes

  out_dir <- file.path(summary_base, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  all_rows <- sub_de[order(sub_de$Gene, sub_de$zone), ]
  fwrite(all_rows, file.path(out_dir, paste0(label, ".csv")))

  # Formatted table, the basis of the Figure 2d panel.
  pdf(file.path(out_dir, paste0(label, ".pdf")),
      width = 14, height = max(4, 0.3 * nrow(all_rows) + 2))
  display <- all_rows[, .(Gene, zone, logFC = round(logFC, 4),
                          AveExpr = round(AveExpr, 4), t = round(t, 3),
                          P.Value = signif(P.Value, 3),
                          adj.P.Val = signif(adj.P.Val, 3),
                          B = round(B, 2), direction)]
  gridExtra::grid.table(display, rows = NULL,
                        theme = gridExtra::ttheme_minimal(base_size = 9))
  dev.off()

  # Pull in the matching per-gene zonation heatmap from step 10, so each
  # category folder is self-contained for figure assembly.
  missing <- character(0)
  for (gene in sig_genes) {
    src <- file.path(heatmaps_dir, paste0("Zone_heatmap_across_age_gender_", gene, ".png"))
    if (file.exists(src)) {
      file.copy(src, file.path(out_dir, paste0(gene, "_heatmap.png")), overwrite = TRUE)
    } else {
      missing <- c(missing, gene)
    }
  }
  if (length(missing) > 0) {
    cat("[step 17]  ", label, ": no heatmap for", length(missing),
        "genes (run step 10):", paste(head(missing, 5), collapse = ", "), "\n")
  }

  cat("[step 17]", label, ":", length(sig_genes), "genes,", nrow(all_rows),
      "gene-zone rows ->", out_dir, "\n")
}


# --- BLOCK 4: ENRICHMENT PER CATEGORY  -->  FIGURE 2c ------------------------
#
# GO Biological Process and KEGG on each of the four gene sets. Figure 2c
# places the up-in-old terms and the up-in-young terms on opposite sides of a
# shared axis, which is why the two directions are enriched separately rather
# than as one list.
#
# The same caveat as step 13 applies: the background is the full mouse
# annotation, not the 1,000-gene panel, so the p-values rank terms rather than
# calibrate them. With sets of 25-40 genes, few terms survive multiple-testing
# correction, and the plots below show terms at p < 0.05 uncorrected.

if (requireNamespace("clusterProfiler", quietly = TRUE)) {

  library(clusterProfiler)
  library(org.Mm.eg.db)

  enrich_dir <- file.path(summary_base, "enrichment_plots")
  dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)

  for (label in names(category_genes)) {
    sig_genes <- category_genes[[label]]
    gene_df <- suppressWarnings(
      bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    )
    if (nrow(gene_df) < 3) {
      cat("[step 17]", label, ": too few mapped genes for enrichment\n")
      next
    }

    ego <- enrichGO(gene = gene_df$ENTREZID, OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                    pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    ekegg <- tryCatch(
      enrichKEGG(gene = gene_df$ENTREZID, organism = "mmu", pvalueCutoff = 0.05),
      error = function(e) { cat("[step 17] KEGG unavailable:", conditionMessage(e), "\n"); NULL }
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego),
                file.path(enrich_dir, paste0(label, "_GO_BP_enrichment.csv")),
                row.names = FALSE)
      ggsave(file.path(enrich_dir, paste0(label, "_GO_dotplot.png")),
             dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP:", label)),
             width = 7, height = 6, bg = "white")
      ggsave(file.path(enrich_dir, paste0(label, "_GO_barplot.png")),
             barplot(ego, showCategory = 15) + ggtitle(paste("GO BP:", label)),
             width = 7, height = 6, bg = "white")
    }
    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      write.csv(as.data.frame(ekegg),
                file.path(enrich_dir, paste0(label, "_KEGG_enrichment.csv")),
                row.names = FALSE)
      ggsave(file.path(enrich_dir, paste0(label, "_KEGG_dotplot.png")),
             dotplot(ekegg, showCategory = 15) + ggtitle(paste("KEGG:", label)),
             width = 7, height = 6, bg = "white")
    }
    cat("[step 17]", label, ": enrichment written\n")
  }
} else {
  cat("[step 17] clusterProfiler not installed; skipping enrichment.\n")
}

cat("\n[step 17] summary ->", summary_base, "\n")
