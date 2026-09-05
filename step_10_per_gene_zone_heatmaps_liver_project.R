# =============================================================================
# STEP 10: Per-gene zonation heatmaps and line plots  -->  FIGURES 2e-h, 4e
# =============================================================================
#
# By Nancy R. Zhang
#
# Produces one small heatmap per gene, showing that gene alone across the
# porto-central axis, with one row per age group: four female rows on top, two
# male rows below. Figures 2f-h and 4e are assembled from these panels.
#
# HOW TO READ ONE PANEL
#   Six rows, six age-by-sex groups. The x-axis is the 20 coarse zonation bins
#   from step 8. Because the array is cross-sample normalized against the young
#   wildtype of the SAME sex, colour is comparable down a column within a sex:
#   a row that is warmer than the 4mo row above it means genuinely higher
#   expression in the older animals, not just a different dynamic range.
#
#   Female and male blocks are normalized against different references, so
#   comparing a female row to a male row is not meaningful. The 2-row gap
#   between the blocks is there to discourage exactly that.
#
# WHY EVERY GENE, NOT JUST THE SIGNIFICANT ONES
#   selgenes is the full harmonized gene set. Generating a panel for every gene
#   is cheap, and it means that when a gene later turns up as significant in
#   the limma analysis (steps 15-16) its panel already exists and can be
#   dropped straight into a figure. Step 17 copies the relevant subset into the
#   per-category summary folders.
#
# NON-FINITE VALUES
#   A gene whose young-wildtype reference maximum is zero divides to Inf or
#   NaN. Those cells are set to 1, i.e. "equal to the reference", which is the
#   neutral colour. This affects only genes with no detectable young signal;
#   such a gene should not be interpreted from this panel.
#
# INPUT   zone_average_crosssample_genenormed, metatab (steps 7-8)
# OUTPUT  FIGURES_DIR/gene_heatmaps/Zone_heatmap_across_age_gender_<GENE>.png
#         FIGURES_DIR/gene_lineplots/<GENE>_lineplot.png  (optional, block 4)
#
# RUN TIME: roughly one second per gene, so tens of minutes for the full set.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("zone_average_crosssample_genenormed")) {
  source("step_8_zonation_binning_liver_project.R")
}

library(ComplexHeatmap)
library(grid)

heatmap_dir <- file.path(FIGURES_DIR, "gene_heatmaps")
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

# Set to TRUE to also produce the per-gene line plots of block 4. They show the
# variation between biological replicates that a heatmap of group means hides.
MAKE_LINEPLOTS <- FALSE


# --- BLOCK 1: AGE-GROUP MEAN PROFILES ----------------------------------------
#
# Four female groups and two male groups, each averaged over its replicates.
# fwy2 is excluded from the young female group: it was flagged in QC.
#
# The matching standard deviations are computed alongside. They are not drawn
# on the heatmaps, but they are what block 4 uses and are worth inspecting
# before believing a small difference between two age rows.

group_mean <- function(idx) apply(zone_average_crosssample_genenormed[, , idx, drop = FALSE], c(1, 2), mean)
group_sd   <- function(idx) apply(zone_average_crosssample_genenormed[, , idx, drop = FALSE], c(1, 2), sd)

sel_f1 <- which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                metatab$age_in_months < 12 & metatab$handle != "fwy2")
sel_f2 <- which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                metatab$age_in_months == 18)
sel_f3 <- which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                metatab$age_in_months == 24)
sel_f4 <- which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                metatab$age_in_months > 24)
sel_m1 <- which(metatab$gender == "male" & metatab$mutant == "wildtype" &
                metatab$age_in_months < 12)
sel_m2 <- which(metatab$gender == "male" & metatab$mutant == "wildtype" &
                metatab$age_in_months > 12)

zone_average_f1 <- group_mean(sel_f1); zone_sd_f1 <- group_sd(sel_f1)
zone_average_f2 <- group_mean(sel_f2); zone_sd_f2 <- group_sd(sel_f2)
zone_average_f3 <- group_mean(sel_f3); zone_sd_f3 <- group_sd(sel_f3)
zone_average_f4 <- group_mean(sel_f4); zone_sd_f4 <- group_sd(sel_f4)
zone_average_m1 <- group_mean(sel_m1); zone_sd_m1 <- group_sd(sel_m1)
zone_average_m2 <- group_mean(sel_m2); zone_sd_m2 <- group_sd(sel_m2)

age_group_profiles <- list(zone_average_f1, zone_average_f2, zone_average_f3,
                           zone_average_f4, zone_average_m1, zone_average_m2)
age_group_labels <- c("4mo", "18mo", "24mo", ">28mo", "4mo", ">28mo")
age_group_sex <- c(1, 1, 1, 1, 2, 2)   # row_split: female block, male block


# --- BLOCK 2: WHICH GENES ----------------------------------------------------

selgenes <- genes
cat("[step 10] generating panels for", length(selgenes), "genes\n")


# --- BLOCK 3: DRAW ------------------------------------------------------------
#
# Deliberately small (400 x 110 px): these are figure components, sized to be
# stacked into the multi-gene panels of Figures 2 and 4 without rescaling.

for (i in 1:length(selgenes)) {
  plotmat <- matrix(nrow = 0, ncol = nbins)
  for (prof in age_group_profiles) {
    plotmat <- rbind(plotmat, prof[match(selgenes[i], rownames(prof)), ])
  }
  plotmat[which(!is.finite(plotmat), arr.ind = TRUE)] <- 1   # see header note

  rownames(plotmat) <- age_group_labels

  ht <- Heatmap(plotmat, cluster_rows = FALSE, cluster_columns = FALSE,
                name = " ",
                row_split = age_group_sex, gap = unit(2, "mm"),
                row_title = c("female", "male"),
                show_column_names = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 8),
                row_title_gp = gpar(fontsize = 8))

  png(file.path(heatmap_dir,
                paste0("Zone_heatmap_across_age_gender_", selgenes[i], ".png")),
      height = 110, width = 400)
  draw(ht, padding = unit(c(2, 2, 2, 1), "mm"))
  dev.off()
}

cat("[step 10] wrote", length(selgenes), "heatmaps to", heatmap_dir, "\n")


# --- BLOCK 4: LINE PLOTS (OPTIONAL) ------------------------------------------
#
# One line per SAMPLE rather than per age group, so that between-replicate
# spread is visible. A difference between age groups that is smaller than the
# spread between replicates of the same age should not be trusted, and these
# plots are the quickest way to see that.

if (MAKE_LINEPLOTS) {
  lineplot_dir <- file.path(FIGURES_DIR, "gene_lineplots")
  dir.create(lineplot_dir, showWarnings = FALSE, recursive = TRUE)

  for (sex in c("female", "male")) {
    sel <- which(metatab$gender == sex & metatab$mutant == "wildtype")
    for (g in selgenes) {
      gi <- match(g, genes)
      png(file.path(lineplot_dir, paste0(g, "_", sex, "_lineplot.png")),
          width = 500, height = 400)
      matplot(t(zone_average_crosssample_genenormed[gi, , sel]),
              type = "l", lty = metatab$plotlty[sel], col = metatab$plotcolor[sel],
              xlab = "zonation bin (pericentral -> periportal)",
              ylab = "expression / young reference",
              main = paste(g, "-", sex))
      legend("topleft", legend = metatab$handle[sel], col = metatab$plotcolor[sel],
             lty = metatab$plotlty[sel], cex = 0.7, bty = "n")
      grid()
      dev.off()
    }
  }
  cat("[step 10] wrote line plots to", lineplot_dir, "\n")
}
