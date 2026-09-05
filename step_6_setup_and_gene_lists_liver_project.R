# =============================================================================
# STEP 6: Setup -- libraries, gene lists, and gene categories
# =============================================================================
#
# By Nancy R. Zhang
#
# This is the first step of the downstream analysis. Steps 1-5 (Alisha Aristel)
# take raw Visium output through cropping, clipping, masking, normalization and
# zonation scoring. Steps 6-17 start from that output and produce the Visium
# panels of Figures 1-4.
#
# This step defines no data and computes no results. It loads the packages and
# builds the gene lists that every later step refers to. It is sourced, not
# run on its own.
#
# INPUT
#   gene_lists/mfuzz.csv                     zonation clusters from GeoMx
#                                            (gene name, cluster id 1-5)
#   gene_lists/zonation_genes.txt            21 canonical zonation markers
#   gene_lists/senescence_genes_20250307.txt senescence panel, curated by the
#                                            Niedernhofer lab
#   gene_lists/fibroblast_genes.txt          fibroblast markers
#
# OUTPUT (variables placed in the global environment)
#   genes_pericentral, genes_midlobular, genes_periportal
#   genes_senescence_all, genes_sasp, genes_cellcycle, genes_antiapoptosis
#   genes_housekeeping, genes_immune, genes_endothelial, genes_fibroblast
#   genes_all                 union of everything above
#   gene_category_colors      colour key reused by every heatmap
#
# RUN TIME: a few seconds.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

source("config.R")


# --- BLOCK 1: LIBRARIES ------------------------------------------------------
#
# Plotting:   ggplot2, pheatmap, ComplexHeatmap, gplots, gridExtra, colorspace
# Data:       data.table, dplyr, magrittr, matrixStats
# Statistics: limma
# Spatial:    fields, MBA (surface interpolation for the tissue maps)
# Progress:   pbapply, tictoc

library_paths <- c(
  "pheatmap", "fields", "ggplot2", "gplots", "gridExtra",
  "data.table", "tictoc", "limma"
)
invisible(lapply(library_paths, library, character.only = TRUE))

library(MBA)
library(matrixStats)
library(magrittr)
library(dplyr)
library(colorspace)
library(pbapply)

# Diverging palette for the expression heatmaps, lightened so that the extreme
# reds and blues do not dominate at small panel sizes.
RdBucols <- lighten(RColorBrewer::brewer.pal(11, "RdBu"), amount = 0.3)
immune_color_scheme <- lighten(RColorBrewer::brewer.pal(9, "Oranges"), amount = 0.1)


# --- BLOCK 2: ZONATION CLUSTERS FROM GeoMx -----------------------------------
#
# mfuzz.csv assigns each gene to one of five soft clusters describing its
# expression profile across the three zones, from the GeoMx analysis in Fig 1b.
# Rows with a missing cluster id are genes that were not assigned; drop them.

tab <- read.table(file.path(GENE_LISTS_DIR, "mfuzz.csv"), sep = ",", header = TRUE)
sel <- which(is.na(tab$x))
tab <- tab[-sel, ]
colnames(tab) <- c("genename", "clusterid")
genes_zonation_clusters <- tab


# --- BLOCK 3: GENE LISTS -----------------------------------------------------
#
# Three lists come from files; the rest are short enough to state inline, which
# keeps the exact panel used for the paper visible in the code.

genes_senescence_all <- readLines(file.path(GENE_LISTS_DIR, "senescence_genes_20250307.txt"))
genes_zonation       <- readLines(file.path(GENE_LISTS_DIR, "zonation_genes.txt"))
genes_fibroblast     <- readLines(file.path(GENE_LISTS_DIR, "fibroblast_genes.txt"))

# Zonation markers, by zone. Zone 1 = periportal, zone 2 = midlobular,
# zone 3 = pericentral. These define the expected direction of the zonation
# score and are the genes plotted in Figure 1c.
genes_midlobular  <- c("Hamp", "Hamp2", "Ccnd1")
genes_pericentral <- c("Cyp2e1", "Gsta3", "Cyp27a1", "Mup17", "Nt5e", "Axin2",
                       "Cyp1a2", "Oat", "Gstm3", "Axin2", "Lgr5")
genes_periportal  <- c("Alb", "Cyp2f2", "Asl", "Gls2", "Cdh1", "Cps1", "Pck1", "Sdhd")

# Housekeeping genes: the denominator of the per-pixel normalization applied in
# preprocessing, and the basis of the PC1 artifact correction.
genes_housekeeping <- c("Gapdh", "Actb", "Rplp0", "Vcl", "L3mbtl2", "Rbck1",
                        "Vamp7", "Wdr55", "Hprt")

genes_immune      <- c("Cd45", "Cd68", "Cd64", "Cd3", "Cd4", "Cd8")
genes_endothelial <- c("Cd31", "Cd144", "Vwf")

# The senescence panel splits into three modules that enter the senescence
# score (step 11) with different exponents.
genes_cellcycle <- c("Cdkn1a", "Cdkn2a", "Cdkn2b", "Cdkn2d", "Cdkn1b", "Cdkn2c",
                     "Cdkn1c", "Ccna2", "Ccnb1", "Ccnd1", "Ccnd2", "Ccng1")
genes_antiapoptosis <- c("Bcl2l1", "Bcl2l2", "Mcl1")

# SASP = everything in the senescence panel that is not a cell-cycle inhibitor,
# not anti-apoptotic, and not one of the two markers that go DOWN in senescence.
genes_sasp <- setdiff(genes_senescence_all,
                      c(genes_cellcycle, genes_antiapoptosis, "Mki67", "Lmnb1"))

genes_downsenescent <- c("Mki67", "Lmnb1")
genes_proliferation <- c("Mki67")

genes_all <- c(genes_zonation, genes_senescence_all, genes_fibroblast,
               genes_housekeeping, genes_immune, genes_endothelial,
               genes_midlobular, genes_periportal, genes_pericentral,
               genes_zonation_clusters[, 1])
genes_all <- unique(genes_all)


# --- BLOCK 4: COLOUR KEY -----------------------------------------------------
#
# One colour per gene category, used consistently across every annotated
# heatmap so that a reader who learns the key on one figure can carry it to
# the next.

gene_category_colors <- c(
  "Pericentral"       = "darkgreen",
  "Midlobular"        = "forestgreen",
  "Periportal"        = "lightgreen",
  "Housekeeping"      = "gray",
  "Fibroblast"        = "purple",
  "Cell_Cycle"        = "orange",
  "SASP"              = "red",
  "Anti_Apoptosis"    = "brown",
  "Endothelial"       = "cornflowerblue",
  "Immune"            = "hotpink",
  "Down_in_Senescent" = "cyan"
)

cat("[step 6] setup complete:", length(genes_all), "genes across",
    length(gene_category_colors), "categories\n")
