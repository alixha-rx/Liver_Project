# =============================================================================
# STEP 12: Cross-sample senescence association  -->  FIGURES 3a-b, 3d-f, 4a, 4c
# =============================================================================
#
# By Nancy R. Zhang
#
# Step 11 produced, for each sample independently, the correlation of every
# non-panel gene with the senescence score. This step asks the question that
# makes those correlations worth anything: DO THEY REPRODUCE ACROSS BIOLOGICAL
# REPLICATES?
#
# THE LOGIC OF FIGURES 4a AND 4c
#   Each old animal gives one correlation per gene. Plot animal against animal
#   and you get a scatter with one point per gene. If the senescence score were
#   picking up sample-specific noise, the points would form a cloud with no
#   structure. They do not -- the correlations agree between animals, which is
#   what the reported Pearson r and its p-value quantify. Only after that check
#   is it reasonable to take the genes that correlate consistently and ask what
#   pathways they belong to (step 13).
#
#   k-means with k = 3 splits the genes into positively associated, unassociated
#   and negatively associated groups, and colours the scatter by that grouping.
#   For males the cluster centres are SUPPLIED (0.2 / 0.0 / -0.1) rather than
#   found, which fixes the meaning of the three clusters; for females they are
#   found from the data.
#
# WHAT ELSE IS HERE
#   Block 4 writes the per-sample senescence score summaries that become
#   Figures 3a and 3b: box plots of the pixel-level score distribution for each
#   animal, with the 90th / 95th / 99th percentiles marked.
#   Block 3 writes the spatial expression maps of the genes Linshan selected.
#
# THRESHOLDS
#   Block 5 recomputes the senescence-high thresholds from the data: the 98th
#   percentile of the score in the oldest group of each sex. These are the
#   values hard-coded as SENSCORE_THRESH_* in config.R and used for the black
#   overlay dots in Figures 3d and 3h. They are recomputed here so that the
#   provenance of those two numbers is visible, and so that a change in the
#   sample set is caught rather than silently inherited.
#
# INPUT   senscore4_list, genecors_list (step 11); samples, metatab (step 7)
#         gene_lists/senescence_genes_{female,male}.txt  genes selected by
#         Linshan Laux for the spatial panels
# OUTPUT  RESULTS_DIR/senescence_association/{male,female}_senescence_association.csv
#         FIGURES_DIR/senescence_association/*_gene_pairs.png     (Fig 4a, 4c)
#         FIGURES_DIR/senescence_summary/*                        (Fig 3a, 3b)
#         FIGURES_DIR/spatial_maps_senescence/{female,male}/*.png (Fig 3e-f)
#         male_cors, female_cors  data frames, the input to step 13
#
# RUN TIME: 10-20 minutes, most of it the spatial maps.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

if (!exists("senscore4_list")) source("step_11_senescence_scoring_liver_project.R")

library(GGally)

set.seed(1)   # k-means below; see README, "Known caveats"

assoc_fig_dir <- file.path(FIGURES_DIR, "senescence_association")
assoc_res_dir <- file.path(RESULTS_DIR, "senescence_association")
summary_dir   <- file.path(FIGURES_DIR, "senescence_summary")
for (d in c(assoc_fig_dir, assoc_res_dir, summary_dir,
            file.path(FIGURES_DIR, "spatial_maps_senescence", "female"),
            file.path(FIGURES_DIR, "spatial_maps_senescence", "male"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ggally_cor groups by colour by default, which would report a correlation per
# k-means cluster. Strip the aesthetic so the printed r is the full-sample
# correlation -- the number quoted in the figure legend.
my_cor <- function(data, mapping, ...) {
  mapping$colour <- NULL
  mapping$color  <- NULL
  mapping$fill   <- NULL
  GGally::ggally_cor(data = data, mapping = mapping, ...)
}

genecors <- genecors_list


# --- BLOCK 1: MALE REPLICATES  -->  FIGURE 4c --------------------------------

sel <- which((metatab$gender == "male") & (metatab$age_in_months > 12))
male_cors <- do.call("cbind", genecors[sel])
colnames(male_cors) <- metatab$handle[sel]
male_cors[which(is.na(male_cors), arr.ind = TRUE)] <- 0

# Fixed centres: positively associated, unassociated, negatively associated.
centers <- matrix(data = c(0.2, 0.2, 0.2, 0, 0, 0, -0.1, -0.1, -0.1),
                  nrow = 3, byrow = TRUE)
kk <- kmeans(male_cors, center = centers)
male_cors <- as.data.frame(male_cors)
male_cors$cluster <- as.factor(kk$cluster)

p1 <- ggpairs(male_cors,
              columns = colnames(male_cors)[1:(ncol(male_cors) - 1)],
              mapping = aes(color = cluster),
              lower = list(continuous = wrap("points", alpha = 0.7, size = 1.5)),
              upper = list(continuous = my_cor),
              legend = 1)
p1 <- p1 + theme(legend.position = "right")
ggsave(file.path(assoc_fig_dir, "male_senescence_association_gene_pairs.png"),
       p1, width = 7, height = 5.5, bg = "white")

male_cors$known_senescence <- rownames(male_cors) %in% genes_senescence_all
write.csv(male_cors, file = file.path(assoc_res_dir, "male_senescence_association.csv"))

# How many genes clear r > 0.1 in ALL male replicates, and how many of those
# were already known senescence genes? The novel ones are the subject of Fig 4e.
temp <- rowSums(male_cors[, c(1:3)] > 0.1)
cat("[step 12] male: genes above r>0.1 in n replicates:\n"); print(table(temp))
cat("[step 12] male: consistent genes already in the senescence panel:",
    sum(rownames(male_cors[which(temp == 3), ]) %in% genes_senescence_all), "\n")


# --- BLOCK 2: FEMALE REPLICATES  -->  FIGURE 4a ------------------------------
#
# Ordered by age so the pair panels read youngest to oldest.

sel <- which((metatab$gender == "female") & (metatab$age_in_months > 12))
o <- order(metatab$age_in_months[sel])
sel <- sel[o]
female_cors <- do.call("cbind", genecors[sel])
colnames(female_cors) <- metatab$handle[sel]
female_cors[which(is.na(female_cors), arr.ind = TRUE)] <- 0

kk <- kmeans(female_cors, center = 3)
female_cors <- as.data.frame(female_cors)
female_cors$cluster <- as.factor(kk$cluster)

p1 <- ggpairs(female_cors,
              columns = colnames(female_cors)[1:(ncol(female_cors) - 1)],
              mapping = aes(color = cluster),
              legend = 1,
              lower = list(continuous = wrap("points", alpha = 0.7, size = 0.8)),
              upper = list(continuous = my_cor))
p1 <- p1 + theme(legend.position = "right")
ggsave(file.path(assoc_fig_dir, "female_senescence_association_gene_pairs.png"),
       p1, width = 7, height = 5.5, bg = "white")

female_cors$known_senescence <- rownames(female_cors) %in% genes_senescence_all
write.csv(female_cors, file = file.path(assoc_res_dir, "female_senescence_association.csv"))

temp <- rowSums(female_cors[, 1:6] > 0.2)
cat("[step 12] female: genes above r>0.2 in n replicates:\n"); print(table(temp))


# --- BLOCK 3: SPATIAL MAPS OF SELECTED GENES  -->  FIGURE 3e-f, 3i-j ---------
#
# The gene lists are the ones Linshan Laux selected for the figure panels.
# One map per gene per old animal, so that a spatial pattern can be checked
# against every replicate rather than shown for a single favourable section.

plot_spatial_panel <- function(gene_file, sex, out_folder) {
  genes_to_plot <- readLines(gene_file)
  old_samps <- which(metatab$mutant == "wildtype" & metatab$gender == sex &
                     metatab$age_in_months > 12)
  for (i in 1:length(genes_to_plot)) {
    genename <- genes_to_plot[i]
    for (handle in old_samps) {
      dat <- samples[[handle]]$datfm_inbox[, ]
      coords <- samples[[handle]]$loc
      maskpix <- samples[[handle]]$maskpix
      zone <- samples[[handle]]$zone31rqt

      if (DO_MASK && length(maskpix) == nrow(dat)) {
        dat <- dat[!maskpix, , drop = FALSE]
        zone <- zone[!maskpix]
        coords <- coords[!maskpix, , drop = FALSE]
      }
      if (!(genename %in% colnames(dat))) next   # gene absent from the panel

      pgene <- spatial_map_gene(genename, dat, coords)
      ggsave(file.path(out_folder,
                       paste(genename, "_", metatab$handle[handle], "_spatial.png", sep = "")),
             pgene, width = 3, height = 3, bg = "white")
    }
  }
}

plot_spatial_panel(file.path(GENE_LISTS_DIR, "senescence_genes_female.txt"), "female",
                   file.path(FIGURES_DIR, "spatial_maps_senescence", "female"))
plot_spatial_panel(file.path(GENE_LISTS_DIR, "senescence_genes_male.txt"), "male",
                   file.path(FIGURES_DIR, "spatial_maps_senescence", "male"))


# --- BLOCK 4: PER-SAMPLE SCORE DISTRIBUTIONS  -->  FIGURE 3a, 3b -------------
#
# One box per animal, ordered and coloured by age group, with the 90th, 95th
# and 99th percentiles of the pixel-level score marked. The paper reads the
# UPPER TAIL of these distributions rather than the median: senescent pixels
# are rare, so a shift in the 95th and 99th percentiles is the signal, and the
# median barely moves.

p1 <- plot_faceted_violins(senscore4_list, metatab, gender = "female")
ggsave(file.path(summary_dir, "senscore4_faceted_violins_female.png"),
       plot = p1, height = 3, width = 6, bg = "white")

p1 <- plot_faceted_percentile_boxes(senscore4_list, metatab, gender = "female",
                                    transform = NULL)
ggsave(file.path(summary_dir, "senscore4_faceted_percentiles_female.png"),
       plot = p1, height = 3, width = 6, bg = "white")

p1 <- plot_faceted_violins(senscore4_list, metatab, gender = "male")
ggsave(file.path(summary_dir, "senscore4_faceted_violins_male.png"),
       plot = p1, height = 3, width = 3, bg = "white")

p1 <- plot_faceted_percentile_boxes(senscore4_list, metatab, gender = "male",
                                    transform = NULL)
ggsave(file.path(summary_dir, "senscore4_faceted_percentiles_male.png"),
       plot = p1, height = 3, width = 6, bg = "white")


# --- BLOCK 5: SENESCENCE-HIGH THRESHOLDS -------------------------------------
#
# Defined as the 98th percentile of the score in the OLDEST group of each sex,
# then applied to every animal of that sex. The fraction of pixels in the
# younger animals that clears the same bar is printed: that fraction is the
# quantitative claim behind "senescence-high pixels accumulate with age".
#
# These are the numbers stored as SENSCORE_THRESH_FEMALE / _MALE in config.R.
# If they disagree with the values printed here, the sample set has changed.

sel <- which(metatab$gender == "male" & metatab$age_in_months > 24)
senscore4_thresh_male <- quantile(unlist(senscore4_list[sel]), 0.98)
cat("[step 12] male threshold (98th pct of >24mo):",
    round(senscore4_thresh_male, 4),
    " config.R has:", SENSCORE_THRESH_MALE, "\n")

sel <- which(metatab$gender == "female" & metatab$age_in_months == 24)
senscore4_thresh_female <- quantile(unlist(senscore4_list[sel]), 0.98)
cat("[step 12] female threshold (98th pct of 24mo):",
    round(senscore4_thresh_female, 4),
    " config.R has:", SENSCORE_THRESH_FEMALE, "\n")

for (grp in list(list("female", 18), list("female", 30), list("male", 4))) {
  sex <- grp[[1]]; age <- grp[[2]]
  thresh <- if (sex == "female") senscore4_thresh_female else senscore4_thresh_male
  sel <- if (age == 30) which(metatab$gender == sex & metatab$age_in_months > 24)
         else which(metatab$gender == sex & metatab$age_in_months == age)
  if (length(sel) == 0) next
  v <- unlist(senscore4_list[sel])
  cat(sprintf("[step 12] %s %smo: %.3f%% of pixels above threshold\n",
              sex, age, 100 * sum(v > thresh) / length(v)))
}
