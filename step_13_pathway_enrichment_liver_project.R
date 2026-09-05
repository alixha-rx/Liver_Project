# =============================================================================
# STEP 13: GO and KEGG enrichment of senescence-associated genes  -->  FIG 4b, 4d
# =============================================================================
#
# By Nancy R. Zhang
#
# Takes the genes whose expression correlates POSITIVELY and reproducibly with
# the senescence score (step 12) and asks which biological processes and KEGG
# pathways they belong to.
#
# WHAT IS TESTED, AND AGAINST WHAT
#   The input set is the positively-associated k-means cluster from step 12,
#   with the score's own constituent genes removed (SASP, cell-cycle inhibitor
#   and anti-apoptosis panels). Removing them matters: they would otherwise
#   guarantee enrichment for inflammation and cell-cycle terms by construction,
#   since they are the genes the score is built from. What survives is the set
#   of genes that TRACK senescence without having been used to define it.
#
#   enrichGO and enrichKEGG both test against the full annotated background for
#   Mus musculus, not against the 1,000-gene Visium panel. The panel is itself
#   enriched for immune and senescence genes, so the reported p-values are
#   optimistic in absolute terms. They are used here to RANK pathways, which is
#   how the figure presents them, and not as calibrated significance.
#
# SELECTING THE POSITIVE CLUSTER
#   The original code took cluster 1. For males the k-means centres are
#   supplied in step 12, so cluster 1 is by construction the positively
#   associated one. For females the centres are found from the data, and the
#   label the positive cluster receives depends on initialization -- cluster 1
#   is not guaranteed to be the positive one. This step therefore selects the
#   cluster with the HIGHEST MEAN CORRELATION rather than trusting the label,
#   and prints which cluster that turned out to be.
#
# KEGG GeneRatio FILTER
#   The KEGG dot plots are filtered to pathways where more than a threshold
#   fraction of the input genes map to the pathway (5% male, 8% female) and
#   capped at 30 terms. Without it the plot fills with large generic pathways
#   that a handful of genes hit. The two thresholds differ because the female
#   input set is larger.
#
# INPUT   male_cors, female_cors (step 12)
# OUTPUT  RESULTS_DIR/enrichment/*_{GO,KEGG}_enrichment.csv
#         FIGURES_DIR/pathway_plots/*.png   (dot, bar, emap, tree, upset)
#
# RUN TIME: a few minutes. enrichKEGG queries the KEGG REST API, so this step
# needs a network connection and can fail if that service is unavailable.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("male_cors")) source("step_12_senescence_association_liver_project.R")

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

enrich_res_dir <- file.path(RESULTS_DIR, "enrichment")
enrich_fig_dir <- file.path(FIGURES_DIR, "pathway_plots")
dir.create(enrich_res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(enrich_fig_dir, showWarnings = FALSE, recursive = TRUE)


# --- BLOCK 1: PICK THE POSITIVELY ASSOCIATED CLUSTER -------------------------

positive_cluster_genes <- function(cors_df, label) {
  # Correlation columns only: drop `cluster` and `known_senescence`.
  cor_cols <- setdiff(colnames(cors_df), c("cluster", "known_senescence"))
  cluster_means <- tapply(rowMeans(cors_df[, cor_cols, drop = FALSE]),
                          cors_df$cluster, mean)
  pos <- names(which.max(cluster_means))
  cat("[step 13]", label, "cluster mean correlations:\n")
  print(round(cluster_means, 4))
  cat("[step 13]", label, "-> using cluster", pos, "as the positive cluster\n")
  if (pos != "1") {
    cat("[step 13] NOTE:", label, "positive cluster is not cluster 1;",
        "the label differs from the original run but the gene set is the same.\n")
  }
  rownames(cors_df)[cors_df$cluster == pos]
}

cor_sen_genes_pos_male   <- positive_cluster_genes(male_cors, "male")
cor_sen_genes_pos_female <- positive_cluster_genes(female_cors, "female")


# --- BLOCK 2: ENRICHMENT ------------------------------------------------------

run_enrichment <- function(gene_symbols, exclude_genes, sex, ratio_thresh, top_n = 30) {

  # Drop the score's own constituent genes, then map symbols to Entrez IDs.
  # bitr warns about symbols it cannot map; those are dropped silently by
  # clusterProfiler, so the input count is printed for comparison.
  input <- setdiff(gene_symbols, exclude_genes)
  gene_df <- bitr(input, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)
  entrez_genes <- gene_df$ENTREZID
  cat("[step 13]", sex, ":", length(input), "genes in,",
      length(entrez_genes), "mapped to Entrez\n")

  ego <- enrichGO(gene = entrez_genes,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)

  ekegg <- enrichKEGG(gene = entrez_genes,
                      organism = "mmu",
                      pvalueCutoff = 0.05)

  write.csv(as.data.frame(ego),
            file.path(enrich_res_dir,
                      paste0("senescence_DE_genes_pos_", sex, "_allzones_GO_enrichment.csv")),
            row.names = FALSE)
  write.csv(as.data.frame(ekegg),
            file.path(enrich_res_dir,
                      paste0("senescence_DE_genes_pos_", sex, "_allzones_KEGG_enrichment.csv")),
            row.names = FALSE)

  # Filter KEGG to the informative pathways (see header).
  res2 <- ekegg@result %>%
    mutate(GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    })) %>%
    filter(GeneRatio_num > ratio_thresh) %>%
    arrange(pvalue) %>%
    slice_head(n = top_n)
  res2$GeneRatio_num <- NULL

  enr_filt <- ekegg
  enr_filt@result <- res2

  list(ego = ego, ekegg = ekegg, enr_filt = enr_filt, n_kegg = nrow(res2))
}


# --- BLOCK 3: PLOTS -----------------------------------------------------------

save_enrichment_plots <- function(e, sex) {
  tag <- paste0("(positive cluster, ", sex, ")")
  fig <- function(name) file.path(enrich_fig_dir, paste0(sex, "_", name, ".png"))

  dotplot(e$ego, showCategory = 15) + ggtitle(paste("GO top pathways", tag))
  ggsave(fig("go_dotplot"), height = 6, width = 6)

  # Fig 4b / 4d: KEGG pathways ranked by enrichment p-value, dot size = gene
  # count, colour = adjusted p-value.
  dotplot(e$enr_filt, showCategory = e$n_kegg) + ggtitle(paste("KEGG top pathways", tag))
  ggsave(fig("kegg_dotplot"), height = 5.5, width = 6)

  barplot(e$ego, showCategory = 15, title = paste("GO top pathways", tag))
  ggsave(fig("go_barplot"), height = 6, width = 6)

  barplot(e$enr_filt, showCategory = e$n_kegg) + ggtitle(paste("KEGG top pathways", tag))
  ggsave(fig("kegg_barplot"), height = 5.5, width = 6)

  # Enrichment map and tree: terms linked by shared genes, which shows when a
  # long list of significant terms is really one redundant cluster.
  edo <- pairwise_termsim(e$ego)
  emapplot(edo, showCategory = 15) + ggtitle(paste("GO enrichment map", tag))
  ggsave(fig("go_graphplot"), height = 6, width = 6)
  treeplot(edo, showCategory = 15) + ggtitle(paste("GO enrichment tree", tag))
  ggsave(fig("go_treeplot"), height = 6, width = 6)

  edo <- pairwise_termsim(e$ekegg)
  emapplot(edo, showCategory = 15) + ggtitle(paste("KEGG enrichment map", tag))
  ggsave(fig("kegg_graphplot"), height = 6, width = 6)
  treeplot(edo, showCategory = 15) + ggtitle(paste("KEGG enrichment tree", tag))
  ggsave(fig("kegg_treeplot"), height = 6, width = 6)

  upsetplot(e$ego, n = 15) + ggtitle(paste("GO upset plot", tag))
  ggsave(fig("go_upsetplot"), height = 6, width = 6)
}


# --- BLOCK 4: RUN -------------------------------------------------------------

enrich_male <- run_enrichment(
  cor_sen_genes_pos_male,
  exclude_genes = c(genes_cellcycle, genes_antiapoptosis, genes_sasp_male),
  sex = "male", ratio_thresh = 0.05
)
save_enrichment_plots(enrich_male, "male")

enrich_female <- run_enrichment(
  cor_sen_genes_pos_female,
  exclude_genes = c(genes_cellcycle, genes_antiapoptosis, genes_sasp_female),
  sex = "female", ratio_thresh = 0.08
)
save_enrichment_plots(enrich_female, "female")

# Kept under the original names for anything downstream that expects them.
ego_male <- enrich_male$ego; ekegg_male <- enrich_male$ekegg
ego_female <- enrich_female$ego; ekegg_female <- enrich_female$ekegg

cat("[step 13] enrichment tables ->", enrich_res_dir, "\n")
cat("[step 13] enrichment plots  ->", enrich_fig_dir, "\n")
