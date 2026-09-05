# =============================================================================
# STEP 9: Zonation marker heatmaps across age  -->  FIGURE 1c
# =============================================================================
#
# By Nancy R. Zhang
#
# Produces the figure that shows zone contraction and expansion: for each
# canonical zonation marker, one row per age group, colour showing expression
# along the porto-central axis.
#
# HOW TO READ THE OUTPUT
#   Each gene occupies a block of rows, one row per age group, oldest at the
#   bottom of its block. The x-axis runs periportal (left) to pericentral
#   (right) across the 100 fine bins from step 8.
#
#   A periportal marker whose warm band widens rightward with age, or a
#   pericentral marker whose band narrows, is the loss of zonation the paper
#   reports. Midlobular markers (Hamp, Hamp2, Ccnd1) spreading outward is the
#   zone 2 expansion described in the abstract.
#
# TWO DISPLAY CHOICES WORTH KNOWING
#   1. The matrix is square-rooted before plotting. This is display only -- it
#      compresses the dynamic range so that a gene with one very bright bin
#      does not flatten the rest of its own row.
#   2. The columns are REVERSED (plotmat[, ncol:1]) so that the axis reads
#      periportal -> pericentral, matching the published figure. Bin 1 out of
#      step 8 is pericentral, so without this reversal the axis would run the
#      other way. The grid.text labels underneath state the direction.
#
# WHY FEMALES GET FOUR AGE GROUPS AND MALES TWO
#   The female series was sampled at 4, 18, 24 and >28 months; the male series
#   only at 4 and >28. The panels are therefore not symmetric, and the male
#   heatmap should not be read as a coarser version of the female one.
#
# INPUT   zone_average_crosssample_genenormed_fine, metatab (steps 7-8)
# OUTPUT  FIGURES_DIR/zonation_heatmaps/Zonation_markers_heatmap_{female,male}_WT.png
#
# RUN TIME: seconds.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("zone_average_crosssample_genenormed_fine")) {
  source("step_8_zonation_binning_liver_project.R")
}

library(ComplexHeatmap)
library(grid)

out_dir <- file.path(FIGURES_DIR, "zonation_heatmaps")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# --- BLOCK 1: THE MARKER PANEL -----------------------------------------------
#
# A reduced marker set: the genes with the cleanest, most canonical zonation in
# these data. Ordered periportal -> midlobular -> pericentral so that the row
# order matches the left-to-right reading of the heatmap.

genes_periportal_small  <- c("Alb", "Cyp2f2", "Asl", "Gls2", "Cdh1", "Pck1")
genes_midlobular_small  <- c("Hamp", "Hamp2", "Ccnd1")
genes_pericentral_small <- c("Cyp2e1", "Cyp1a2", "Oat", "Gstm3")

selgenes <- c(genes_periportal_small, genes_midlobular_small, genes_pericentral_small)
gene_group <- c(rep("periportal",  length(genes_periportal_small)),
                rep("midlobular",  length(genes_midlobular_small)),
                rep("pericentral", length(genes_pericentral_small)))

# Wider gap where the zone group changes, so the three blocks read as blocks.
row_gap_vec <- unit(ifelse(gene_group[-1] != gene_group[-length(gene_group)], 6, 2), "mm")


# --- BLOCK 2: HELPER ---------------------------------------------------------
#
# Averages the fine-binned, cross-sample-normalized array over the samples in
# one age group, giving a gene x bin matrix for that group.

mean_profile <- function(sample_idx) {
  apply(zone_average_crosssample_genenormed_fine[, , sample_idx, drop = FALSE],
        c(1, 2), mean)
}

# Assembles the stacked matrix and draws it. `age_profiles` is a list of
# gene x bin matrices, one per age group; `age_labels` names them.
draw_zonation_heatmap <- function(age_profiles, age_labels, outfile) {

  plotmat <- matrix(nrow = 0, ncol = nbins_fine)
  row_groups <- rep(0, 0)
  row_title <- rep("", 0)

  for (i in 1:length(selgenes)) {
    for (a in seq_along(age_profiles)) {
      prof <- age_profiles[[a]]
      plotmat <- rbind(plotmat, prof[match(selgenes[i], rownames(prof)), ])
    }
    row_groups <- c(row_groups, rep(i, length(age_profiles)))
    row_title <- c(row_title, selgenes[i])
  }
  rownames(plotmat) <- rep(age_labels, length(selgenes))

  plotmat <- sqrt(plotmat)                      # display transform only
  plotmat <- plotmat[, ncol(plotmat):1]         # periportal on the left

  ht <- Heatmap(plotmat, cluster_rows = FALSE, cluster_columns = FALSE,
                name = "SQRT Expr.",
                row_split = row_groups, row_gap = row_gap_vec, row_title = row_title,
                show_column_names = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 10))

  png(outfile, height = 800, width = 600)
  draw(ht, padding = unit(c(10, 5, 10, 2), "mm"))
  grid.text("Periportal ----", x = unit(0.05, "npc"), y = unit(5, "mm"), just = "left")
  grid.text("--- Midlobular ---", x = unit(0.5, "npc"), y = unit(5, "mm"), just = "center")
  grid.text("--- Pericentral", x = unit(0.9, "npc"), y = unit(5, "mm"), just = "right")
  dev.off()
  cat("[step 9] wrote", outfile, "\n")
}


# --- BLOCK 3: FEMALE, FOUR AGE GROUPS ----------------------------------------

female_profiles <- list(
  mean_profile(which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                     metatab$age_in_months < 12)),
  mean_profile(which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                     metatab$age_in_months == 18)),
  mean_profile(which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                     metatab$age_in_months == 24)),
  mean_profile(which(metatab$gender == "female" & metatab$mutant == "wildtype" &
                     metatab$age_in_months > 24))
)

draw_zonation_heatmap(female_profiles,
                      c("4mo", "18mo", "24mo", ">28mo"),
                      file.path(out_dir, "Zonation_markers_heatmap_female_WT.png"))


# --- BLOCK 4: MALE, TWO AGE GROUPS -------------------------------------------

male_profiles <- list(
  mean_profile(which(metatab$gender == "male" & metatab$mutant == "wildtype" &
                     metatab$age_in_months < 12)),
  mean_profile(which(metatab$gender == "male" & metatab$mutant == "wildtype" &
                     metatab$age_in_months > 12))
)

draw_zonation_heatmap(male_profiles,
                      c("4mo", ">28mo"),
                      file.path(out_dir, "Zonation_markers_heatmap_male_WT.png"))
