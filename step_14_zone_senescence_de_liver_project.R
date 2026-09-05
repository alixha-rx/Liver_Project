# =============================================================================
# STEP 14: Zone-stratified senescence DE  -->  FIGURES 4f, 4g, 4h
# =============================================================================
#
# By Nancy R. Zhang
#
# Asks, WITHIN each hepatic zone of each old animal, which genes differ between
# senescence-high and senescence-low pixels. Figures 4f-h show the answer for
# the zonation markers: whether zone identity genes are lost specifically in
# the senescent pixels of their own zone.
#
# ---------------------------------------------------------------------------
# THE COMPARISON
# ---------------------------------------------------------------------------
# For one animal and one zone, pixels are split into senescence-high and
# senescence-low, and every gene is tested between the two groups with a
# Wilcoxon rank-sum test, BH-adjusted within the animal-zone.
#
# This is a WITHIN-ANIMAL, WITHIN-ZONE comparison, so it is not confounded by
# age, batch, or section quality: the two groups being compared are pixels from
# the same tissue section. What it cannot do is establish that the difference
# generalizes, which is why the results are only reported where they replicate
# across every old animal of a sex (block 3).
#
# The test is on pixels, not animals, so the p-values are pixel-level. With
# tens of thousands of pixels per zone almost anything reaches significance;
# the asterisks on the figures are therefore read as effect direction and
# consistency across replicates, not as evidence of a population effect. The
# log2 fold change is the quantity that carries the biology.
#
# ---------------------------------------------------------------------------
# TWO PARTS
# ---------------------------------------------------------------------------
# PART A (block 1, off by default) computes the per-animal per-zone tables.
#        It runs against the FULL ~5,000-gene superpixel tables, not the
#        1,000-gene harmonized panel used elsewhere, because the point is to
#        discover genes rather than to test a chosen panel. That data is a
#        separate deposit; set RUN_PART_A and DATA_5K_DIR to regenerate.
#
# PART B (blocks 2-4) reads those tables and builds the figures. This is the
#        part that reproduces from the deposited intermediate files.
#
# ---------------------------------------------------------------------------
# INPUT   DE_TABLES_DIR/DE_perm_OLD_<handle>_<zone>.tsv
#           columns: Gene, log2FC, pwilcox_adj_BH
#         metatab, gene lists (step 7)
# OUTPUT  FIGURES_DIR/zone_senescence_de/*.png
#         RESULTS_DIR/zone_senescence_de/*.csv
#
# RUN TIME: Part B, under a minute. Part A, hours.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("metatab")) source("step_7_load_samples_liver_project.R")

library(pheatmap)

# DE_TABLES_DIR is set in config.R.
RUN_PART_A <- FALSE          # see Part A note above
DATA_5K_DIR <- NULL          # full 5k-gene superpixel tables, if running Part A

out_fig <- file.path(FIGURES_DIR, "zone_senescence_de")
out_res <- file.path(RESULTS_DIR, "zone_senescence_de")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)
dir.create(out_res, showWarnings = FALSE, recursive = TRUE)

# NOTE the zone order here. It is periportal -> midlobular -> pericentral,
# the REVERSE of ZONES in config.R, and it is the order the DE table filenames
# use. Kept as-is so the filenames match; do not reorder without renaming.
zonestr <- c("periportal", "midlobular", "pericentral")

oldhandles_female <- metatab$handle[metatab$gender == "female" & metatab$age_in_months > 12]
oldhandles_male   <- metatab$handle[metatab$gender == "male"   & metatab$age_in_months > 12]


# --- BLOCK 1 (PART A): GENERATE THE PER-ANIMAL TABLES ------------------------
#
# compute_corr_senescence_per_zone() (R/loading.R) does the work: it takes the
# full pixel table plus the senescence score from step 11, splits by zone, and
# runs the Wilcoxon test gene by gene.

if (RUN_PART_A) {
  if (is.null(DATA_5K_DIR)) stop("Set DATA_5K_DIR to the 5k-gene superpixel tables.")
  if (!exists("senscore4_list")) source("step_11_senescence_scoring_liver_project.R")

  dir.create(DE_TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
  infiles <- list.files(DATA_5K_DIR, pattern = "_withpca\\.RData$",
                        full.names = TRUE, recursive = TRUE)

  for (handle in which(metatab$mutant == "wildtype" & metatab$age_in_months > 12)) {
    fid <- grep(paste(metatab$handle[handle], "_", sep = ""), infiles)
    if (length(fid) == 0) { cat("  no 5k file for", metatab$handle[handle], "\n"); next }
    cat("[step 14] Part A:", metatab$handle[handle], "\n")

    load(infiles[fid[1]])                       # provides sample_list
    stopifnot(exists("sample_list"))
    s <- sample_list

    dfbig <- as.data.table(s$datfm_inbox)
    z <- s$zone31rqt
    stopifnot(nrow(dfbig) == length(z))
    gene_cols <- colnames(dfbig)[-c(1, 2)]
    dfbig[, zone31rqt := z]
    if (!is.null(s$maskpix) & DO_MASK) dfbig <- dfbig[!s$maskpix, ]

    # Align the 5k table to the pixels the score was computed on, by matching
    # on the (x, y) coordinate pair.
    coords <- samples[[handle]]$loc[!samples[[handle]]$maskpix, ]
    coords2 <- dfbig[, c(1, 2)]
    basenum <- 100000
    key1 <- coords[, 1] * basenum + coords[, 2]
    key2 <- coords2[[1]] * basenum + coords2[[2]]
    idx <- match(key2, key1)
    dfbig[, senscore := senscore4_list[[handle]][idx]]

    compute_corr_senescence_per_zone(dfbig, gene_cols, DE_TABLES_DIR,
                                     metatab$handle[handle])
  }
}


# --- BLOCK 2 (PART B): READ AND COMBINE THE TABLES ---------------------------
#
# For each zone, assemble a genes x animals matrix of adjusted p-values and one
# of log2 fold changes, restricted to genes present in every animal.

read_zone_tables <- function(handles) {
  out <- vector("list", 3)
  for (zid in 1:3) {
    combinedpval <- NULL; combinedeff <- NULL
    for (i in seq_along(handles)) {
      f <- file.path(DE_TABLES_DIR,
                     paste("DE_perm_OLD_", handles[i], "_", zonestr[zid], ".tsv", sep = ""))
      if (!file.exists(f)) stop("missing DE table: ", f)
      temp <- read.table(f, sep = "\t", header = TRUE)

      pvals <- temp$pwilcox_adj_BH; names(pvals) <- temp$Gene
      log2FC <- temp$log2FC;        names(log2FC) <- temp$Gene

      if (i == 1) {
        combinedpval <- data.frame(pvals); rownames(combinedpval) <- temp$Gene
        combinedeff  <- data.frame(log2FC); rownames(combinedeff) <- temp$Gene
      } else {
        all_genes <- intersect(rownames(combinedpval), temp$Gene)
        combinedpval <- cbind(combinedpval[match(all_genes, rownames(combinedpval)), , drop = FALSE],
                              pvals[match(all_genes, names(pvals))])
        combinedeff  <- cbind(combinedeff[match(all_genes, rownames(combinedeff)), , drop = FALSE],
                              log2FC[match(all_genes, names(log2FC))])
        rownames(combinedpval) <- all_genes
        rownames(combinedeff)  <- all_genes
      }
    }
    colnames(combinedpval) <- handles
    colnames(combinedeff)  <- handles
    out[[zid]] <- list(p_adj_wilcoxon = combinedpval, log2FC = combinedeff)
  }
  out
}

sende_male   <- read_zone_tables(oldhandles_male)
sende_female <- read_zone_tables(oldhandles_female)


# --- BLOCK 3: DO THE EFFECTS AGREE ACROSS ZONES? -----------------------------
#
# Average each gene's log2FC over animals within a zone, then correlate the
# three zone profiles. High correlation means the senescence-associated
# expression change is largely shared across the lobule; low correlation means
# it is zone-specific. This is the summary behind the "senotype" comparison.

zone_effect_matrix <- function(sende, n_animals) {
  m <- cbind(rowSums(sende[[1]]$log2FC) / n_animals,
             rowSums(sende[[2]]$log2FC) / n_animals,
             rowSums(sende[[3]]$log2FC) / n_animals)
  colnames(m) <- zonestr
  m
}

effects_across_zones_male   <- zone_effect_matrix(sende_male,   length(oldhandles_male))
effects_across_zones_female <- zone_effect_matrix(sende_female, length(oldhandles_female))

for (sx in c("male", "female")) {
  m <- if (sx == "male") effects_across_zones_male else effects_across_zones_female
  png(file.path(out_fig, paste0("Zone_senotype_heatmap_", sx, ".png")),
      height = 300, width = 300)
  pheatmap(cor(m, use = "complete.obs"), cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("white", "orange", "red"))(21))
  dev.off()
  write.csv(m, file.path(out_res, paste0("zone_effects_", sx, ".csv")))
}


# --- BLOCK 4: ZONATION MARKERS BY ZONE  -->  FIGURES 4f, 4g, 4h --------------
#
# Rows are zonation markers, grouped periportal / midlobular / pericentral.
# Columns are animal x zone, grouped by zone. Colour is log2FC between
# senescence-high and senescence-low pixels; asterisks mark the Wilcoxon
# adjusted p-value (*** < 0.001, ** < 0.01, * < 0.05).
#
# The question the figure answers: in the pixels that look senescent, do the
# markers OF THAT ZONE go down? A negative (blue) block on the diagonal --
# zone 1 markers down within zone 1, and so on -- is loss of zone identity in
# senescent tissue.
#
# COLUMN ORDER: the female matrix is reordered to c(5,6,3,4,1,2). That is not
# arbitrary -- it puts the six old female animals in ascending age order, which
# the raw handle order does not. It is asserted below rather than assumed.

sig_stars <- function(pval_mat) {
  ifelse(pval_mat < 0.001, "***",
         ifelse(pval_mat < 0.01, "**",
                ifelse(pval_mat < 0.05, "*", "")))
}

make_marker_figure <- function(sende, handles, col_order, sex) {

  has_genes <- rownames(sende[[2]]$p_adj_wilcoxon)
  marker_sets <- list(intersect(genes_periportal,  has_genes),
                      intersect(genes_midlobular,  has_genes),
                      intersect(genes_pericentral, has_genes))
  breaks <- seq(-1.2, 1.2, length.out = 101)

  log2fc_big <- NULL; sig_big <- NULL

  for (i in 1:3) {
    sel_genes <- marker_sets[[i]]
    log2fc_slice <- matrix(nrow = length(sel_genes), ncol = 0)
    sig_slice <- matrix(nrow = length(sel_genes), ncol = 0)

    for (zid in 1:3) {
      pval_mat <- sende[[zid]]$p_adj_wilcoxon[sel_genes, col_order, drop = FALSE]
      log2fc_mat <- sende[[zid]]$log2FC[sel_genes, col_order, drop = FALSE]
      sig_mat <- sig_stars(pval_mat)

      log2fc_slice <- cbind(log2fc_slice, log2fc_mat)
      sig_slice <- cbind(sig_slice, sig_mat)

      # Individual panel for the diagonal cases (markers of zone i shown within
      # zone i), plus midlobular markers in every zone since zone 2 expansion
      # is the effect of interest.
      if (zid == i | i == 2) {
        row_ord <- order(rowSums(log2fc_mat), decreasing = TRUE)
        png(file.path(out_fig,
                      paste("sen_log2fc_zone_", i, "_genes_in_", zid, "_", sex, ".png", sep = "")),
            height = 30 + 20 * length(sel_genes), width = 100 + 30 * ncol(pval_mat))
        pheatmap(log2fc_mat[row_ord, , drop = FALSE], cluster_rows = FALSE, cluster_cols = FALSE,
                 display_numbers = sig_mat[row_ord, , drop = FALSE], number_color = "black",
                 fontsize_number = 12, breaks = breaks)
        dev.off()
      }
    }
    log2fc_big <- rbind(log2fc_big, log2fc_slice)
    sig_big <- rbind(sig_big, sig_slice)
  }

  n <- length(col_order)
  gaps_row <- c(length(marker_sets[[1]]),
                length(marker_sets[[1]]) + length(marker_sets[[2]]))
  gaps_col <- c(n, 2 * n)

  png(file.path(out_fig, paste0(sex, "_senescence_zonation.png")),
      height = 600, width = 300)
  pheatmap(log2fc_big, cluster_rows = FALSE, cluster_cols = FALSE,
           gaps_row = gaps_row, gaps_col = gaps_col,
           display_numbers = sig_big, number_color = "black", fontsize_number = 12)
  dev.off()

  write.csv(log2fc_big, file.path(out_res, paste0(sex, "_senescence_zonation_log2FC.csv")))
  cat("[step 14] wrote", sex, "marker figure:", nrow(log2fc_big), "genes x",
      ncol(log2fc_big), "animal-zones\n")
}

# Female: reorder columns to ascending age, and check that is what we get.
female_col_order <- c(5, 6, 3, 4, 1, 2)
female_ages <- metatab$age_in_months[match(oldhandles_female, metatab$handle)]
if (!identical(order(female_ages[female_col_order]), seq_along(female_col_order))) {
  warning("female_col_order no longer sorts the old female animals by age; ",
          "check the handle order before reading Figure 4f-h column groups.")
}
make_marker_figure(sende_female, oldhandles_female, female_col_order, "female")

# Male: three animals, all ~30 months, so no reordering is needed.
make_marker_figure(sende_male, oldhandles_male, seq_along(oldhandles_male), "male")


# --- BLOCK 5: GENES THAT REPLICATE ACROSS EVERY ANIMAL -----------------------
#
# Summarize each gene by its WORST adjusted p-value across animals within a
# zone, and its mean log2FC. Taking the maximum p-value is a deliberately
# conservative "replicates in all animals" filter: a gene passes only if it was
# significant in every one of them.

summarize_replication <- function(sende, sex) {
  n_genes <- nrow(sende[[1]]$p_adj_wilcoxon)
  max_pval <- matrix(nrow = n_genes, ncol = 3)
  avg_log2FC <- matrix(nrow = n_genes, ncol = 3)
  for (zid in 1:3) {
    max_pval[, zid] <- apply(sende[[zid]]$p_adj_wilcoxon, 1, max, na.rm = TRUE)
    avg_log2FC[, zid] <- apply(sende[[zid]]$log2FC, 1, mean, na.rm = TRUE)
  }
  colnames(max_pval) <- zonestr;   rownames(max_pval) <- rownames(sende[[1]]$log2FC)
  colnames(avg_log2FC) <- zonestr; rownames(avg_log2FC) <- rownames(sende[[1]]$log2FC)

  write.csv(max_pval,   file.path(out_res, paste0("senescence_max_wilcoxon_p_", sex, "_old.csv")))
  write.csv(avg_log2FC, file.path(out_res, paste0("senescence_avg_log2FC_", sex, "_old.csv")))
  list(max_pval = as.data.frame(max_pval), avg_log2FC = as.data.frame(avg_log2FC))
}

repl_male   <- summarize_replication(sende_male, "male")
repl_female <- summarize_replication(sende_female, "female")

cat("[step 14] done. Figures ->", out_fig, "\n")
