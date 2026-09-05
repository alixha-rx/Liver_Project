# =============================================================================
# STEP 15: Pseudobulk expression by animal and zone
# =============================================================================
#
# By Nancy R. Zhang
#
# Collapses the pixel-level data into one expression value per gene, per
# animal, per zone. This matrix is the input to the limma analysis in step 16
# and is the point at which the analysis stops being spatial.
#
# ---------------------------------------------------------------------------
# WHY PSEUDOBULK
# ---------------------------------------------------------------------------
# A Visium section contains tens of thousands of pixels, but they are not
# independent observations: they come from one animal. Testing at the pixel
# level would treat a single mouse as tens of thousands of replicates and
# return p-values that reflect pixel count rather than biological effect.
#
# Averaging to one value per animal-zone puts the analysis at the level of the
# experimental unit -- the mouse -- so that n is the number of animals, which
# is what the limma model in step 16 assumes. It is a large loss of data and a
# necessary one.
#
# ---------------------------------------------------------------------------
# THE THREE-ZONE SPLIT
# ---------------------------------------------------------------------------
# Masked pixels are removed FIRST, then the remaining pixels are ranked by
# zonation score and cut into equal thirds. Masking before ranking matters:
# masked pixels are not randomly distributed (vessels sit pericentrally), so
# ranking first and masking second would leave the three zones with unequal
# pixel counts and shifted boundaries.
#
#   bottom third = pericentral    middle = midlobular    top = periportal
#
# ---------------------------------------------------------------------------
# SAMPLE SELECTION
# ---------------------------------------------------------------------------
# Wildtype only -- the Ercc1 mutants are a separate comparison. fwy3 is
# excluded for quality (EXCLUDE_HANDLES in config.R).
#
# Age groups: young below 18 months, old at or above it for females and above
# it for males. The asymmetry is an artifact of the sampling -- there is an
# 18-month female group but no 18-month male group -- and has no effect on
# which animals land in which group.
#
# ---------------------------------------------------------------------------
# INPUT   DATA_RDS_DIR/<handle>.rds     per-animal pixel tables
#         DATA_RDS_DIR/metatab.csv      sample sheet
# OUTPUT  CACHE_DIR/expr_pseudobulk.rds   genes x (animal x zone) matrix
#         CACHE_DIR/pheno_pseudobulk.csv  one row per column of the above
#         FIGURES_DIR/limma_diagnostics/pseudobulk_diagnostic_heatmap.png
#
# RUN TIME: 10-20 minutes on the first run; instant afterwards, since the
# result is cached and reused.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

source("config.R")

library(limma)
library(data.table)
library(pheatmap)

zones <- ZONES
pseudobulk_file <- file.path(CACHE_DIR, "expr_pseudobulk.rds")
pheno_file      <- file.path(CACHE_DIR, "pheno_pseudobulk.csv")
diag_dir        <- file.path(FIGURES_DIR, "limma_diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)


# --- BLOCK 1: BUILD OR LOAD THE PSEUDOBULK MATRIX ----------------------------

if (file.exists(pseudobulk_file) && file.exists(pheno_file)) {

  cat("[step 15] loading cached pseudobulk from", CACHE_DIR, "\n")
  expr <- readRDS(pseudobulk_file)
  pheno <- as.data.frame(fread(pheno_file))
  rownames(pheno) <- pheno$col
  stopifnot(identical(colnames(expr), pheno$col))

} else {

  # --- 1a. sample sheet and selection ---
  metatab_rds <- read.csv(file.path(DATA_RDS_DIR, "metatab.csv"), stringsAsFactors = FALSE)
  cat("[step 15] metatab:", nrow(metatab_rds), "samples\n")

  wt <- metatab_rds[metatab_rds$mutant == "wildtype" &
                    !metatab_rds$handle %in% EXCLUDE_HANDLES, ]
  cat("[step 15] wildtype samples used:\n")
  print(wt[, c("handle", "gender", "age_months")])

  wt$age_grp <- NA
  fem <- wt$gender == "female"
  wt$age_grp[fem & wt$age_months <  AGE_YOUNG_MAX] <- "young"
  wt$age_grp[fem & wt$age_months >= AGE_YOUNG_MAX] <- "old"
  mal <- wt$gender == "male"
  wt$age_grp[mal & wt$age_months <  AGE_YOUNG_MAX] <- "young"
  wt$age_grp[mal & wt$age_months >  AGE_YOUNG_MAX] <- "old"

  cat("[step 15] sample counts by sex and age group:\n")
  print(table(wt$gender, wt$age_grp))

  # --- 1b. first pass: the gene set common to all animals ---
  # Done as a separate pass so that only one animal is in memory at a time;
  # the pixel tables are several GB each.
  cat("[step 15] finding common genes...\n")
  all_gene_sets <- list()
  for (i in seq_len(nrow(wt))) {
    h <- wt$handle[i]
    s <- readRDS(file.path(DATA_RDS_DIR, paste0(h, ".rds")))
    all_gene_sets[[h]] <- colnames(s$datfm_inbox)[-(1:2)]   # drop x, y
    rm(s); gc(verbose = FALSE)
  }
  common_genes <- Reduce(intersect, all_gene_sets)
  cat("[step 15]  ", length(common_genes), "genes common to all", nrow(wt), "animals\n")

  # --- 1c. second pass: zone means ---
  cat("[step 15] computing pseudobulk zone averages...\n")
  pseudobulk_list <- list()
  for (i in seq_len(nrow(wt))) {
    h <- wt$handle[i]
    cat("  ", h, "...")
    s <- readRDS(file.path(DATA_RDS_DIR, paste0(h, ".rds")))

    dat_genes <- s$datfm_inbox[, common_genes]
    zone_raw  <- s$zone31rqt
    maskpix   <- s$maskpix

    keep <- !maskpix                      # mask BEFORE ranking -- see header
    dat_genes <- dat_genes[keep, ]
    zone_raw  <- zone_raw[keep]

    zone_label <- assign_three_zones(zone_raw)   # R/zonation.R

    for (z in zones) {
      sel <- zone_label == z
      pseudobulk_list[[paste0(h, ".", z)]] <- colMeans(dat_genes[sel, , drop = FALSE])
    }
    cat(" done (", sum(keep), "unmasked pixels )\n")
    rm(s, dat_genes); gc(verbose = FALSE)
  }

  expr <- do.call(cbind, pseudobulk_list)
  cat("[step 15] pseudobulk:", nrow(expr), "genes x", ncol(expr), "animal-zones\n")

  # --- 1d. phenotype table, in the same column order as expr ---
  pheno <- data.frame(col = colnames(expr), stringsAsFactors = FALSE)
  pheno$handle <- sub("\\.[^.]+$", "", pheno$col)
  pheno$zone   <- sub("^.*\\.", "", pheno$col)
  pheno$zone   <- factor(pheno$zone, levels = zones)
  col_order <- pheno$col
  pheno <- merge(pheno, wt[, c("handle", "gender", "age_months", "age_grp")],
                 by = "handle", all.x = TRUE)
  pheno <- pheno[match(col_order, pheno$col), ]   # merge reorders; restore it
  rownames(pheno) <- pheno$col

  saveRDS(expr, pseudobulk_file)
  fwrite(pheno, pheno_file)
  cat("[step 15] cached to", CACHE_DIR, "\n")
}


# --- BLOCK 2: CHECK THE ZONE ASSIGNMENT --------------------------------------
#
# The single most consequential thing that can silently go wrong upstream is an
# inverted zonation score, which would swap the pericentral and periportal
# labels and reverse the interpretation of every zone-specific result. Verify
# it directly against marker genes rather than trusting the label.

pc_markers <- intersect(c("Cyp2e1", "Oat", "Gstm3", "Cyp1a2"), rownames(expr))
pp_markers <- intersect(c("Cyp2f2", "Cps1", "Pck1", "Asl", "Gls2"), rownames(expr))

cat("\n[step 15] zone assignment check (marker mean expression):\n")
zone_check <- sapply(zones, function(z) {
  s <- pheno$col[pheno$zone == z]
  c(pericentral_markers = mean(expr[pc_markers, s]),
    periportal_markers  = mean(expr[pp_markers, s]))
})
print(round(zone_check, 4))

if (!(zone_check["pericentral_markers", "pericentral"] >
      zone_check["pericentral_markers", "periportal"] &&
      zone_check["periportal_markers", "periportal"] >
      zone_check["periportal_markers", "pericentral"])) {
  stop("Zone labels look inverted: pericentral markers are not highest in the ",
       "pericentral zone. Check the sign of zone31rqt before going further.")
}
cat("[step 15] zone assignment OK\n\n")


# --- BLOCK 3: DIAGNOSTIC HEATMAP ---------------------------------------------
#
# A handful of genes with known behaviour, row-scaled, with animals ordered by
# sex then zone then age. Read it before trusting step 16: the zone blocks
# should be visibly different from each other, and replicates of the same
# sex-zone-age should look alike. A column that does not resemble its
# neighbours is an animal worth investigating.

diag_genes <- intersect(c("Aldh1b1", "Nxpe2", "Ang", "Cd36", "Serpina7"), rownames(expr))

ann <- pheno[, c("gender", "zone", "age_months", "age_grp")]
ann$zone <- factor(ann$zone, levels = zones)
ord <- order(ann$gender, match(ann$zone, zones), ann$age_months)
expr_diag <- expr[diag_genes, ord, drop = FALSE]
ann_diag  <- ann[ord, , drop = FALSE]

gender_rle <- rle(as.character(ann_diag$gender))
gaps_col <- cumsum(gender_rle$lengths)
gaps_col <- gaps_col[-length(gaps_col)]

col_labels <- paste0(pheno$handle[match(colnames(expr_diag), pheno$col)], " (",
                     substr(as.character(ann_diag$zone), 1, 4), ")")

ann_colors <- list(
  gender  = c(female = "pink", male = "lightblue"),
  zone    = c(pericentral = "darkgreen", midlobular = "gold", periportal = "salmon"),
  age_grp = c(young = "chartreuse3", old = "indianred3")
)

png(file.path(diag_dir, "pseudobulk_diagnostic_heatmap.png"), width = 1400, height = 500)
pheatmap::pheatmap(expr_diag,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   annotation_col = ann_diag[, c("gender", "zone", "age_grp"), drop = FALSE],
                   annotation_colors = ann_colors,
                   labels_col = col_labels, fontsize_col = 8, fontsize_row = 12,
                   gaps_col = gaps_col,
                   main = "Pseudobulk diagnostic: selected genes by zone and age",
                   scale = "row")
dev.off()

cat("[step 15] diagnostic heatmap ->", diag_dir, "\n")
cat("[step 15] pseudobulk ready:", nrow(expr), "genes x", ncol(expr), "animal-zones\n")
