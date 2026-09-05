# =============================================================================
# STEP 11: Per-pixel senescence scoring  -->  FIGURES 3a, 3b, 3c, 3g
# =============================================================================
#
# By Nancy R. Zhang
#
# Assigns every pixel a senescence score, then produces the per-sample
# diagnostic and figure panels built on it.
#
# ---------------------------------------------------------------------------
# THE SCORE
# ---------------------------------------------------------------------------
# Three gene modules are summed within each pixel:
#
#     totalSASP            secreted / inflammatory SASP genes
#     totalCellCycle       cell cycle inhibitors (Cdkn*)
#     totalAntiApoptosis   Bcl2l1, Bcl2l2, Mcl1
#
# Each sum is then converted to a PERCENTILE RANK within a pooled reference
# distribution built from every wildtype pixel of that sex, across all ages
# (blocks 2 and 3). The score multiplies the three percentiles, with exponents:
#
#     senscore = (pct_SASP ^ 1) * (pct_CellCycle ^ 2) * (pct_AntiApoptosis ^ 1)
#
# Two properties follow from that form, and both matter for interpretation:
#
#   * Because it is a PRODUCT, a pixel scores high only if all three modules
#     are elevated together. A pixel with strong SASP but no cell-cycle arrest
#     scores near zero. This is deliberate: it encodes the requirement that
#     senescence show several hallmarks at once, rather than any one of them.
#
#   * The exponent of 2 on the cell-cycle term makes it the dominant factor.
#     Cdkn2a in particular drives most of the between-sample variation in the
#     final score.
#
# Ranks are computed against a POOLED reference across ages, not within each
# sample. That is what makes the score comparable between samples: a value of
# 0.9 means the same thing in a young and an old animal. It also means adding
# or removing a sample changes the reference and shifts every score slightly.
#
# ---------------------------------------------------------------------------
# SEX-SPECIFIC GENE PANELS
# ---------------------------------------------------------------------------
# The SASP and cell-cycle panels DIFFER between female and male (block 1).
# They were derived separately from the female and male data. As a result the
# female and male scores are computed from different genes against different
# reference distributions, and the two are NOT on a common scale. Figure 3a
# and 3b are read within a sex, never across, and the senescence-high
# thresholds are set per sex for the same reason.
#
# ---------------------------------------------------------------------------
# Mki67 GATING
# ---------------------------------------------------------------------------
# After the product is formed, pixels in the top 5% of Mki67 are zeroed. Mki67
# marks proliferating cells, which by definition are not senescent; without the
# gate, dense proliferative regions pick up high scores through the SASP term.
#
# ---------------------------------------------------------------------------
# INPUT   samples, metatab, genes (step 7)
# OUTPUT  senscore4_list   list, one numeric vector of per-pixel scores per
#                          sample (masked pixels already removed)
#         genecors_list    list, per-sample correlation of every non-panel gene
#                          with the senescence score; the input to step 12
#         CACHE_DIR/wildtype_senescence_results.Rdata   both of the above
#         FIGURES_DIR/per_sample_senescence/*.png       per-sample panels
#
# RUN TIME: 30-60 minutes over 13 samples.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("samples")) source("step_7_load_samples_liver_project.R")

do_mask <- DO_MASK

# The original analysis did not set a seed. It is set here so that repeated
# runs agree; see README, "Known caveats", for what this means for exactly
# reproducing the published k-means cluster labels in step 12.
set.seed(1)

out_dir <- file.path(FIGURES_DIR, "per_sample_senescence")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# The exploratory version of this script also built a permutation null: 100
# random gene sets of the same sizes as the three modules, summed per pixel.
# It was never used -- the block that would have turned it into permuted
# scores stayed commented out -- and it cost 100 passes over every pixel of
# every sample, so it is not reproduced here. Nothing in the published results
# depends on it.
#
# Zone-specific senescence correlation is run separately, in step 14, against
# the full 5k-gene tables rather than the harmonized panel.
do_zone_specific_correlation <- FALSE


# --- BLOCK 1: SEX-SPECIFIC GENE PANELS ---------------------------------------

genes_sasp_female <- c("Apobec3", "Ccl4", "Ccl9", "Ccl24", "Cxcl1", "Cxcl14", "Ddit4",
                       "Eda2r", "Esm1", "Gadd45a", "Hpse", "Igfbp3", "Igfbp5", "Il1b",
                       "Il7", "Iqgap2", "Irf7", "Mmp9", "Mmp12", "Mmp12", "Mif", "Nrg1",
                       "Pecam1", "Plk2", "Ptges", "Serpine1", "Tnfaip2")
genes_sasp_male   <- c("Ccl2", "Ccl3", "Ccl4", "Ccl9", "Ccl25", "Cd9", "Csf1", "Ctnnb1",
                       "Ctsd", "Cxcl1", "Ddit4", "Eda2r", "Egfr", "Hpse", "Ifi44", "Ifih1",
                       "Lcp1", "Lgals3bp", "Mif", "Mmp12", "Ptbp1", "Ptges", "Serpine2",
                       "Serpine1", "Stat1", "Susd6", "Usp18", "Wnt16")

genes_cellcycle_female <- c("Cdkn2a", "Cdkn2b")
genes_cellcycle_male   <- c("Cdkn1b", "Cdkn2a")

# Genes NOT in the score. Only these are scanned in block 5 for correlation
# with the score, so that a gene cannot appear correlated simply because it is
# one of the terms being correlated.
non_senscore_genes_male <- genes[!(genes %in% genes_sasp_male) &
                                 !(genes %in% genes_cellcycle) &
                                 !(genes %in% genes_antiapoptosis)]
non_senscore_genes_female <- genes[!(genes %in% genes_sasp_female) &
                                   !(genes %in% genes_cellcycle) &
                                   !(genes %in% genes_antiapoptosis)]


# --- BLOCK 2: BUILD THE POOLED REFERENCE DISTRIBUTIONS -----------------------
#
# Concatenate the three module sums over every wildtype pixel of one sex, then
# rank once. Ranking the pooled vector a single time here, rather than once per
# sample later, is what puts all samples of a sex on one common scale.
#
# The "_control" vectors hold the young (<12mo) subset only. They are retained
# for diagnostics and are not used by the published score.

build_reference <- function(sex, genes_sasp_sel, genes_cellcycle_sel) {

  ref <- list(
    totalSASP = rep(0, 0), totalCellCycle = rep(0, 0),
    totalAntiApoptosis = rep(0, 0), totalProliferation = rep(0, 0),
    totalSASP_control = rep(0, 0), totalCellCycle_control = rep(0, 0),
    totalAntiApoptosis_control = rep(0, 0), totalProliferation_control = rep(0, 0),
    pixtot = rep(0, 0), source_handle = rep(0, 0)
  )

  for (handle in which(metatab$gender == sex & metatab$mutant == "wildtype")) {
    cat("[step 11] reference", sex, "- sample", handle, "\n")
    dat <- samples[[handle]]$datfm_inbox[, ]
    maskpix <- samples[[handle]]$maskpix
    if (do_mask && length(maskpix) == nrow(dat)) {
      dat <- dat[!maskpix, , drop = FALSE]
    }

    totalSASP <- rowSums(dat[, colnames(dat) %in% genes_sasp_sel])
    totalCellCycle <- rowSums(dat[, colnames(dat) %in% genes_cellcycle_sel])
    totalAntiApoptosis <- rowSums(dat[, colnames(dat) %in% genes_antiapoptosis, drop = FALSE])
    totalProliferation <- rowSums(dat[, colnames(dat) %in% genes_proliferation, drop = FALSE])
    pixtot <- rowSums(dat)

    if (metatab$age_in_months[handle] < 12) {
      ref$totalSASP_control <- c(ref$totalSASP_control, totalSASP)
      ref$totalCellCycle_control <- c(ref$totalCellCycle_control, totalCellCycle)
      ref$totalAntiApoptosis_control <- c(ref$totalAntiApoptosis_control, totalAntiApoptosis)
      ref$totalProliferation_control <- c(ref$totalProliferation_control, totalProliferation)
    }

    ref$totalSASP <- c(ref$totalSASP, totalSASP)
    ref$totalCellCycle <- c(ref$totalCellCycle, totalCellCycle)
    ref$totalAntiApoptosis <- c(ref$totalAntiApoptosis, totalAntiApoptosis)
    ref$totalProliferation <- c(ref$totalProliferation, totalProliferation)
    ref$pixtot <- c(ref$pixtot, pixtot)
    ref$source_handle <- c(ref$source_handle, rep(handle, length(totalSASP)))
  }

  # Rank once, here. Every per-sample score is a lookup into these ranks.
  ref$totalSASP_ranks <- rank(ref$totalSASP)
  ref$totalCellCycle_ranks <- rank(ref$totalCellCycle)
  ref$totalAntiApoptosis_ranks <- rank(ref$totalAntiApoptosis)
  ref
}

ref_f <- build_reference("female", genes_sasp_female, genes_cellcycle_female)
ref_m <- build_reference("male",   genes_sasp_male,   genes_cellcycle_male)

cat("[step 11] female reference:", length(ref_f$totalSASP), "pixels\n")
cat("[step 11] male reference:  ", length(ref_m$totalSASP), "pixels\n")


# --- BLOCK 3: PER-SAMPLE SCORING AND PLOTS -----------------------------------

senscore4_thresh_female <- SENSCORE_THRESH_FEMALE
senscore4_thresh_male   <- SENSCORE_THRESH_MALE

genecors_list  <- vector("list", length(samples))
senscore4_list <- vector("list", length(samples))

for (handle in which(metatab$mutant == "wildtype")) {
  cat("[step 11] scoring sample", handle, ":", metatab$handle[handle], "\n")

  dat <- samples[[handle]]$datfm_inbox[, ]
  zone <- samples[[handle]]$zone31rqt
  maskpix <- samples[[handle]]$maskpix
  pixtot <- samples[[handle]]$pixtot
  coords <- samples[[handle]]$loc

  if (do_mask && length(maskpix) == nrow(dat)) {
    dat <- dat[!maskpix, , drop = FALSE]
    zone <- zone[!maskpix]
    pixtot <- pixtot[!maskpix]
    coords <- coords[!maskpix, , drop = FALSE]
  }

  # Pick the sex-matched panels, reference ranks and threshold.
  if (metatab$gender[handle] == "female") {
    genes_sasp_sel <- genes_sasp_female
    genes_cellcycle_sel <- genes_cellcycle_female
    non_senscore_genes <- non_senscore_genes_female
    ref <- ref_f
    senscore4_thresh <- senscore4_thresh_female
  } else {
    genes_sasp_sel <- genes_sasp_male
    genes_cellcycle_sel <- genes_cellcycle_male
    non_senscore_genes <- non_senscore_genes_male
    ref <- ref_m
    senscore4_thresh <- senscore4_thresh_male
  }

  # Module sums for this sample.
  totalSASP <- rowSums(dat[, colnames(dat) %in% genes_sasp])          # full panel
  totalSASP_sel <- rowSums(dat[, colnames(dat) %in% genes_sasp_sel])  # sex panel
  totalCellCycle <- rowSums(dat[, colnames(dat) %in% genes_cellcycle_sel])
  totalAntiApoptosis <- rowSums(dat[, colnames(dat) %in% genes_antiapoptosis, drop = FALSE])
  totalImmune <- rowSums(dat[, colnames(dat) %in% genes_immune, drop = FALSE])
  totalFibroblast <- rowSums(dat[, colnames(dat) %in% genes_fibroblast, drop = FALSE])

  # Diagnostic: how the signatures relate to each other and to zonation. A
  # strong pixtot correlation here would mean the score is tracking total
  # capture rather than biology.
  geneSignatures <- cbind(pixtot, totalSASP_sel, totalCellCycle, totalAntiApoptosis,
                          totalImmune, totalFibroblast, zone)
  colnames(geneSignatures)[ncol(geneSignatures)] <- "zonescore"
  png(file.path(out_dir, paste0(metatab$handle[handle], "_diagnostic_cor_geneSignatures.png")),
      height = 500, width = 500)
  pheatmap(cor(geneSignatures), cluster_rows = FALSE, cluster_cols = FALSE)
  dev.off()

  # --- the score ---
  cat("  computing senescence score\n")
  alpha_sasp <- ALPHA_SASP
  alpha_cell_cycle <- ALPHA_CELL_CYCLE
  alpha_antiapoptosis <- ALPHA_ANTIAPOPTOSIS

  # Look up each pixel value in the pooled reference to get its percentile.
  mm <- match(totalSASP_sel, ref$totalSASP)
  totalSASP_rank_in_all <- ref$totalSASP_ranks[mm]
  mm <- match(totalCellCycle, ref$totalCellCycle)
  totalCellCycle_rank_in_all <- ref$totalCellCycle_ranks[mm]
  mm <- match(totalAntiApoptosis, ref$totalAntiApoptosis)
  totalAntiApoptosis_rank_in_all <- ref$totalAntiApoptosis_ranks[mm]
  npix_all <- length(ref$totalSASP_ranks)

  senscore4 <- ((totalSASP_rank_in_all / npix_all)^alpha_sasp) *
               ((totalCellCycle_rank_in_all / npix_all)^alpha_cell_cycle) *
               ((totalAntiApoptosis_rank_in_all / npix_all)^alpha_antiapoptosis)

  # Zero out proliferating pixels: a dividing cell is not senescent.
  mki67 <- dat[, "Mki67"]
  cdkn2a <- dat[, "Cdkn2a"]
  softblue <- rgb(red = 0, green = 0, blue = 1, alpha = 0.5)
  senscore4 <- senscore4 * (mki67 < quantile(mki67, 0.95))

  senscore4_list[[handle]] <- senscore4

  # --- scatter of the three modules, coloured by the resulting score ---
  df <- data.frame(totalSASP_sel = totalSASP_sel, totalCellCycle = totalCellCycle,
                   totalAntiApoptosis = totalAntiApoptosis, senscore4)
  p <- ggplot(df, aes(x = totalSASP_sel, y = totalCellCycle,
                      color = senscore4, size = totalAntiApoptosis)) +
    geom_point(alpha = 0.5) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(0.1, 4)) +
    theme_minimal() +
    labs(x = "SASP", y = "Cell cycle inhibitors", color = "Senescence score",
         size = "Anti-apoptosis", title = metatab$description[handle])
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_senescence_scatter_senscore4.png")),
         plot = p, width = 6, height = 6, dpi = 300, bg = "white")


  # --- BLOCK 4: SPATIAL PANELS ---------------------------------------------
  #
  # zonescore is the within-sample percentile of the zonation score, so the
  # tissue map is on a 0-1 scale regardless of section size. Senescence-high
  # pixels are overlaid as black squares on the zonation map: this is the
  # panel used for Figures 3d and 3h.

  zonescore <- (rank(zone) / length(zone))
  senscore4_high_pixels <- coords[senscore4 > senscore4_thresh, ]

  psens <- spatial_map_value(senscore4, coords, label = "Senescence", color_theme = "plasma")
  pzonation <- spatial_map_value(zonescore, coords, label = "Zone",
                                 points_to_highlight = senscore4_high_pixels,
                                 color_theme = RdBucols, col_highlight = "black")
  pzonation_nosens <- spatial_map_value(zonescore, coords, label = "Zone", color_theme = RdBucols)

  pimmune <- spatial_map_value(totalImmune, coords, label = "Immune",
                               points_to_highlight = senscore4_high_pixels,
                               color_theme = immune_color_scheme, col_highlight = "black")
  pimmune2 <- spatial_map_value(totalImmune, coords, label = "Immune",
                                color_theme = immune_color_scheme)
  pfibroblast <- spatial_map_value(totalFibroblast, label = "Fibroblast", coords,
                                   points_to_highlight = senscore4_high_pixels,
                                   color_theme = immune_color_scheme, col_highlight = "black")
  pfibroblast2 <- spatial_map_value(totalFibroblast, label = "Fibroblast", coords,
                                    color_theme = immune_color_scheme)

  pgene1 <- spatial_map_gene("Cyp2f2", dat, coords)     # periportal marker, Fig 3e
  pgene2 <- spatial_map_gene("Cyp2e1", dat, coords)     # pericentral marker, Fig 3f
  pgene3 <- spatial_map_gene("Serpine1", dat, coords)   # SASP
  pgene4 <- spatial_map_gene("Spp1", dat, coords)       # SASP

  gg <- grid.arrange(pzonation, pgene1, pgene2, pimmune, psens, pgene3, pgene4,
                     pfibroblast, ncol = 4)
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_senescence_vs_zone_grid.png")),
         gg, width = 18, height = 8)

  gg2 <- grid.arrange(pzonation_nosens, pgene1, pgene2, pimmune2, psens, pgene3, pgene4,
                      pfibroblast2, ncol = 4)
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_senescence_vs_zone_grid2.png")),
         gg2, width = 18, height = 8)

  # Senescence score across zonation bins: Figures 3c and 3g.
  psens_by_zone <- plot_senescence_by_zone_violin(zonescore, senscore4,
                                                  title = metatab$description[handle])
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_senscore4_by_zone_violin.png")),
         psens_by_zone, width = 6, height = 4)

  # SASP against proliferation, restricted to midlobular pixels. If SASP were
  # simply reporting proliferation, this would show a positive trend.
  sel <- which(zonescore > 0.3 & zonescore < 0.7)
  corr <- cor(totalSASP_sel[sel], mki67[sel])
  sasp_vs_mki67 <- ggplot(data.frame(x = totalSASP_sel[sel], y = mki67[sel]), aes(x, y)) +
    geom_hex(bins = 60) +
    theme_minimal() +
    ggtitle(paste(metatab$handle[handle], ", corr=", format(corr, digits = 2))) +
    labs(x = "SASP gene expression", y = "Mki67")
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_sasp_vs_mki67_midlobular.png")),
         sasp_vs_mki67, width = 4, height = 4, bg = "white")

  # Wnt16 / Wnt2 maps. Motivates the WNT signalling argument in the discussion.
  wnt2 <- dat[, "Wnt2"]
  wnt16 <- dat[, "Wnt16"]
  p1 <- spatial_map_value(wnt16, coords, label = "Wnt16")
  gg <- grid.arrange(psens, p1, ncol = 2)
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_wnt16.png")),
         gg, width = 8, height = 4, bg = "white")

  sel <- which(wnt16 > quantile(wnt16, 0.99))
  p1 <- spatial_map_value(pmin(wnt2, 0.1), coords, label = "Wnt2",
                          points_to_highlight = coords[sel, ], col_highlight = "limegreen")
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_wnt2_vs_wnt16.png")),
         p1, width = 4, height = 4, bg = "white")
  p1 <- spatial_map_value(senscore4, coords, label = "Senescence",
                          points_to_highlight = coords[sel, ], col_highlight = "limegreen")
  ggsave(file.path(out_dir, paste0(metatab$handle[handle], "_senescore4_vs_wnt16.png")),
         p1, width = 4, height = 4, bg = "white")


  # --- BLOCK 5: WHICH OTHER GENES TRACK THE SCORE? -------------------------
  #
  # Correlate every gene that is NOT part of the score against the score,
  # within this sample. Step 12 asks which of these correlations reproduce
  # across biological replicates; that reproducibility is Figure 4a and 4c,
  # and the reproducible set is what goes into pathway enrichment.

  corsens <- function(z) cor(z, senscore4)
  genecors <- apply(dat, 2, corsens)
  genecors <- genecors[non_senscore_genes]

  png(file.path(out_dir, paste0(metatab$handle[handle], "_senescence_topgenes.png")),
      width = 500, height = 500)
  plot(sort(genecors), col = "blue", pch = 18, xlab = "Rank", ylab = "Correlation",
       main = "Correlation with senescence score")
  grid(); abline(h = 0, col = "black")
  text(50, max(genecors, na.rm = TRUE) / 2,
       paste(c("Top genes:", names(genecors)[order(genecors, decreasing = TRUE)[1:20]]),
             collapse = "\n"), col = "blue", cex = 0.8)
  dev.off()

  genecors_list[[handle]] <- genecors
}


# --- BLOCK 6: SAVE ------------------------------------------------------------

save(file = file.path(CACHE_DIR, "wildtype_senescence_results.Rdata"),
     genecors_list, senscore4_list)

cat("[step 11] wrote", file.path(CACHE_DIR, "wildtype_senescence_results.Rdata"), "\n")
