# =============================================================================
# STEP 8: Zonation binning -- gene x bin x sample expression arrays
# =============================================================================
#
# By Nancy R. Zhang
#
# This is the core of the zone-resolved analysis. Every zonation figure in the
# paper (Fig 1c, 2e-h, 3c, 4e) is a slice of one of the arrays built here.
#
# THE IDEA
#   Each pixel carries a continuous zonation score. Rank the pixels of a sample
#   on that score and cut the ranks into bins of EQUAL PIXEL COUNT, then average
#   expression within each bin. The result is a profile of each gene along the
#   porto-central axis, on a common x-axis across samples of very different
#   size and shape.
#
#   Equal COUNT rather than equal WIDTH is the important choice. It makes bin i
#   mean the same thing in every sample -- "the i-th slice of the tissue as you
#   walk from central vein to portal triad" -- regardless of how the zonation
#   score happens to be distributed in that particular section.
#
#   Bin 1 = pericentral (zone 3).  Bin nbins = periportal (zone 1).
#
# THREE NORMALIZATIONS, AND WHY THERE ARE THREE
#   The same bin averages are produced under three normalizations, because
#   different figures need different things:
#
#   1. zone_average                          RAW bin means, no normalization.
#      Use when absolute expression matters, or as the input to a later
#      normalization. This is what feeds the cross-sample version below.
#
#   2. zone_average_persample_genenormed     each gene divided by its own max
#      WITHIN that sample. Every gene runs 0-1 in every sample, so the SHAPE of
#      the zonation profile is comparable across samples but the LEVEL is not.
#      Use for "where along the axis is this gene expressed", not "how much".
#
#   3. zone_average_crosssample_genenormed   each gene divided by the mean
#      young-wildtype maximum FOR THAT SEX. Level is preserved and comparable
#      across ages, because every sample is measured against the same young
#      reference. This is the one used for the age-comparison heatmaps, since
#      it is the only one in which "the old sample is higher than the young
#      sample" is a meaningful statement.
#
#      The reference is computed separately for female and male because
#      baseline expression differs substantially between sexes; using a pooled
#      reference would push a sex difference into every age comparison.
#
# TWO RESOLUTIONS
#   NBINS_COARSE (20)  per-gene heatmaps and violins -- enough pixels per bin
#                      for a stable mean.
#   NBINS_FINE  (100)  zonation marker heatmaps, where the shape of the
#                      boundary is the point of the figure.
#
# INPUT   samples, metatab, genes  (from step 7)
# OUTPUT  six arrays of dim (gene x bin x sample), named as above with a
#         _fine suffix for the 100-bin versions, plus a diagnostic scatter
#         written to FIGURES_DIR.
#
# RUN TIME: a few minutes.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("samples")) source("step_7_load_samples_liver_project.R")

do_mask <- DO_MASK
nbins <- NBINS_COARSE
nbins_fine <- NBINS_FINE
axis_transform <- sqrt

ngenes <- length(genes)


# --- BLOCK 1: COARSE BINS (20) -----------------------------------------------

zone_average <- array(dim = c(ngenes, nbins, nsamples))
zone_average_persample_genenormed <- array(dim = c(ngenes, nbins, nsamples))

dimnames(zone_average) <- list(
  gene = genes,
  bin = paste("bin:", c(1:nbins), sep = ""),
  sample = metatab$handle
)
dimnames(zone_average_persample_genenormed) <- dimnames(zone_average)

for (i in 1:nsamples) {
  cat("[step 8] binning sample", i, "of", nsamples, "\n")
  dat <- samples[[i]]$datfm_inbox[, genes]   # pixel-by-gene table
  zone <- samples[[i]]$zone31rqt             # continuous zonation score
  maskpix <- samples[[i]]$maskpix            # TRUE = pixel should be dropped
  colnames(dat) <- genes

  if (do_mask) {
    dat <- dat[!maskpix, ]
    zone <- zone[!maskpix]
  }

  # Raw bin means. normalize_per_gene = FALSE because normalization happens
  # below, against a reference shared across samples.
  res <- compute_bin_average(dat, zone = zone, nbins = nbins,
                             normalize_per_gene = FALSE)
  zone_average[, , i] <- t(res)              # note the transpose: the array is
                                             # gene x bin, res is bin x gene

  # Per-sample per-gene version: divide each gene by its own max across bins.
  gene_max <- apply(res, 2, max)
  res <- sweep(res, 2, gene_max, "/")
  zone_average_persample_genenormed[, , i] <- t(res)
}


# --- BLOCK 2: CROSS-SAMPLE NORMALIZATION (COARSE) ----------------------------
#
# For each gene, take the mean profile of the young wildtype samples of a given
# sex, find its maximum across bins, and divide every sample of that sex by
# that single number. After this, 1.0 means "as high as the young reference
# peak", and values above 1 mean genuinely elevated relative to young.

femaleyoung <- which(metatab$gender == "female" & metatab$age_in_months < 6 &
                     metatab$mutant == "wildtype")
maleyoung <- which(metatab$gender == "male" & metatab$age_in_months < 6 &
                   metatab$mutant == "wildtype")

gene_average_femaleyoung <- apply(zone_average[, , femaleyoung], c(1, 2), mean)
gene_average_maleyoung <- apply(zone_average[, , maleyoung], c(1, 2), mean)
gene_max_female_young_avg <- apply(gene_average_femaleyoung, 1, max)
gene_max_male_young_avg <- apply(gene_average_maleyoung, 1, max)

# Diagnostic: female vs male reference levels. Genes far off the diagonal have
# strongly sex-dependent baseline expression -- worth knowing before reading
# any sex comparison.
png(file.path(FIGURES_DIR, "step8_reference_female_vs_male.png"),
    width = 600, height = 600)
plot(log(gene_max_female_young_avg), log(gene_max_male_young_avg),
     xlab = "log max, young female reference",
     ylab = "log max, young male reference",
     main = "Young-wildtype reference level by sex")
abline(0, 1, col = "red")
grid()
dev.off()

zone_average_crosssample_genenormed <- zone_average
for (i in 1:nsamples) {
  if (metatab$gender[i] == "female") {
    gene_norm_factor <- gene_max_female_young_avg
  } else {
    gene_norm_factor <- gene_max_male_young_avg
  }
  zone_average_crosssample_genenormed[, , i] <-
    sweep(zone_average[, , i], 1, gene_norm_factor, "/")
}


# --- BLOCK 3: FINE BINS (100) ------------------------------------------------
#
# Identical to blocks 1-2 at higher resolution. Used for Figure 1c, where the
# question is exactly where the zone boundaries sit and how they move with age.

zone_average_fine <- array(dim = c(ngenes, nbins_fine, nsamples))
zone_average_persample_genenormed_fine <- array(dim = c(ngenes, nbins_fine, nsamples))

dimnames(zone_average_fine) <- list(
  gene = genes,
  bin = paste("bin:", c(1:nbins_fine), sep = ""),
  sample = metatab$handle
)
dimnames(zone_average_persample_genenormed_fine) <- dimnames(zone_average_fine)

for (i in 1:nsamples) {
  cat("[step 8] fine binning sample", i, "of", nsamples, "\n")
  dat <- samples[[i]]$datfm_inbox[, genes]
  zone <- samples[[i]]$zone31rqt
  maskpix <- samples[[i]]$maskpix
  if (do_mask) {
    dat <- dat[!maskpix, ]
    zone <- zone[!maskpix]
  }
  res <- compute_bin_average(dat, zone = zone, nbins = nbins_fine,
                             normalize_per_gene = FALSE)
  zone_average_fine[, , i] <- t(res)
  gene_max <- apply(res, 2, max)
  res <- sweep(res, 2, gene_max, "/")
  zone_average_persample_genenormed_fine[, , i] <- t(res)
}

gene_average_femaleyoung_fine <- apply(zone_average_fine[, , femaleyoung], c(1, 2), mean)
gene_average_maleyoung_fine <- apply(zone_average_fine[, , maleyoung], c(1, 2), mean)
gene_max_female_young_avg_fine <- apply(gene_average_femaleyoung_fine, 1, max)
gene_max_male_young_avg_fine <- apply(gene_average_maleyoung_fine, 1, max)

zone_average_crosssample_genenormed_fine <- zone_average_fine
for (i in 1:nsamples) {
  if (metatab$gender[i] == "female") {
    gene_norm_factor <- gene_max_female_young_avg_fine
  } else {
    gene_norm_factor <- gene_max_male_young_avg_fine
  }
  zone_average_crosssample_genenormed_fine[, , i] <-
    sweep(zone_average_fine[, , i], 1, gene_norm_factor, "/")
}


# --- BLOCK 4: CACHE ----------------------------------------------------------
#
# These arrays take a few minutes to build and are used by steps 9, 10 and 17.
# Saving them lets a later step be re-run on its own.

saveRDS(list(coarse = zone_average,
             coarse_persample = zone_average_persample_genenormed,
             coarse_crosssample = zone_average_crosssample_genenormed,
             fine = zone_average_fine,
             fine_persample = zone_average_persample_genenormed_fine,
             fine_crosssample = zone_average_crosssample_genenormed_fine,
             nbins = nbins, nbins_fine = nbins_fine),
        file.path(CACHE_DIR, "zone_averages.rds"))

cat("[step 8] wrote", file.path(CACHE_DIR, "zone_averages.rds"), "\n")
