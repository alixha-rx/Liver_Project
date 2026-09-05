# =============================================================================
# STEP 16: Zone-specific differential expression with limma  -->  FIGURE 2b, 2d
# =============================================================================
#
# By Nancy R. Zhang
#
# The published differential expression analysis. Compares old against young
# wildtype animals WITHIN each hepatic zone, separately for each sex.
#
# ---------------------------------------------------------------------------
# THE MODEL
# ---------------------------------------------------------------------------
#     ~ zone + age_grp + zone:age_grp
#
# The interaction term is the whole point. Without it the model would fit a
# single aging effect shared by all three zones; with it, the aging effect is
# allowed to differ by zone, which is the biological question -- does the liver
# age differently at the portal triad than at the central vein?
#
# Three contrasts are then extracted per sex, one per zone:
#     pericentral   age_grpold                                (the baseline zone)
#     midlobular    age_grpold + zonemidlobular:age_grpold
#     periportal    age_grpold + zoneperiportal:age_grpold
#
# Each estimates the old-vs-young log fold change within that zone. Because
# pericentral is the reference level of the zone factor, its contrast is the
# main effect alone; the other two add their interaction term back in.
#
# ---------------------------------------------------------------------------
# duplicateCorrelation: THE PART THAT MATTERS MOST
# ---------------------------------------------------------------------------
# Each animal contributes THREE rows to this matrix -- one per zone. They are
# not independent: a mouse with globally high expression is high in all three
# zones. Treating them as 33 independent observations would badly overstate
# the sample size.
#
# duplicateCorrelation(block = handle) estimates the within-animal correlation
# and lmFit uses it to down-weight accordingly. The consensus correlation comes
# out around 0.87-0.89, which is high, and confirms that the three zones of one
# animal are close to being one observation for many genes. This is what makes
# the resulting p-values defensible.
#
# ---------------------------------------------------------------------------
# WHY 30-MONTH FEMALES ARE EXCLUDED  (this is "pipeline E")
# ---------------------------------------------------------------------------
# Five variants of this analysis were run during development, differing in
# whether the sexes were modelled jointly, whether quantile normalization was
# applied, and whether the two 30-month females were included. This script is
# variant E: separate sexes, quantile normalized, 30+ month females dropped.
#
# Dropping them raised the female count from roughly 7 to roughly 50 genes per
# zone at FDR < 0.2, and left the male results unchanged. The reason is reduced
# heterogeneity in the old female group: with only four old females, two of
# them six months older than the other two, the within-group variance was large
# enough to absorb most of the aging signal.
#
# The female old group in this pipeline is therefore 18- and 24-month animals.
# Any statement about 30-month females rests on Figures 3 and 4, not on this
# differential expression analysis.
#
# ---------------------------------------------------------------------------
# OTHER CHOICES
# ---------------------------------------------------------------------------
# Gene filtering    the lowest-expressed 10% of genes are dropped before
#                   fitting. Low-expression genes have unstable variance and
#                   inflate the multiple-testing burden without contributing.
# Quantile norm     normalizeBetweenArrays puts every animal-zone on a common
#                   distribution, removing per-section differences in overall
#                   capture that survive the per-pixel normalization.
# eBayes trend      trend = TRUE models the mean-variance relationship;
#                   robust = TRUE stops a few extreme genes from dominating the
#                   variance prior. Both matter at this sample size.
# FDR < 0.2         a permissive threshold, chosen because the experiment has
#                   2-4 animals per group. It is the threshold used throughout
#                   the paper and should be read as "worth following up", not
#                   as strong evidence for any single gene.
#
# ---------------------------------------------------------------------------
# INPUT   CACHE_DIR/expr_pseudobulk.rds, CACHE_DIR/pheno_pseudobulk.csv (step 15)
# OUTPUT  RESULTS_DIR/pipeline_E/DE_<zone>_<sex>WT_old_vs_young.csv
#           one row per gene: Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B,
#           direction ("OLD" if logFC > 0, else "YOUNG")
#
# RUN TIME: under a minute.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

source("config.R")

library(limma)
library(data.table)

zones <- ZONES
out_dir <- file.path(RESULTS_DIR, "pipeline_E")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# --- BLOCK 1: INPUT ----------------------------------------------------------

expr_raw <- readRDS(file.path(CACHE_DIR, "expr_pseudobulk.rds"))
pheno <- as.data.frame(fread(file.path(CACHE_DIR, "pheno_pseudobulk.csv")))
rownames(pheno) <- pheno$col
stopifnot(identical(colnames(expr_raw), pheno$col))

pheno$gender  <- factor(pheno$gender, levels = c("female", "male"))
pheno$zone    <- factor(pheno$zone, levels = zones)      # pericentral = reference
pheno$age_grp <- factor(pheno$age_grp, levels = c("young", "old"))

cat("[step 16] pseudobulk:", nrow(expr_raw), "genes x", ncol(expr_raw), "animal-zones\n")


# --- BLOCK 2: GENE FILTERING -------------------------------------------------

gene_means <- rowMeans(expr_raw)
min_expr <- quantile(gene_means, GENE_FILTER_Q)
keep_genes <- gene_means > min_expr
cat("[step 16] keeping", sum(keep_genes), "of", length(keep_genes),
    "genes (min mean expression =", round(min_expr, 4), ")\n")


# --- BLOCK 3: SAMPLE FILTERING (pipeline E) ----------------------------------

keep_E <- !(pheno$gender == "female" & pheno$age_months >= FEMALE_AGE_EXCLUDE)
pheno_E <- droplevels(pheno[keep_E, ])
pheno_E$gender  <- factor(pheno_E$gender, levels = c("female", "male"))
pheno_E$zone    <- factor(pheno_E$zone, levels = zones)
pheno_E$age_grp <- factor(pheno_E$age_grp, levels = c("young", "old"))

cat("[step 16] excluding females >=", FEMALE_AGE_EXCLUDE, "months\n")
cat("[step 16] remaining animal-zones by sex and age group:\n")
print(table(pheno_E$gender, pheno_E$age_grp))


# --- BLOCK 4: NORMALIZATION --------------------------------------------------

expr_E <- expr_raw[keep_genes, keep_E]
expr_E_qn <- normalizeBetweenArrays(expr_E, method = "quantile")


# --- BLOCK 5: FIT, PER SEX ---------------------------------------------------

for (sex in c("female", "male")) {

  keep_sex <- pheno_E$gender == sex
  expr_sex <- expr_E_qn[, keep_sex, drop = FALSE]
  P_sex <- droplevels(pheno_E[keep_sex, , drop = FALSE])
  P_sex$zone <- factor(P_sex$zone, levels = zones)
  P_sex$age_grp <- factor(P_sex$age_grp, levels = c("young", "old"))

  design <- model.matrix(~ zone + age_grp + zone:age_grp, data = P_sex)
  cat("\n[step 16]", sex, ": design", ncol(design), "coefficients,",
      nrow(design), "observations, rank", qr(design)$rank, "\n")
  if (qr(design)$rank < ncol(design)) {
    warning(sex, ": design matrix is rank deficient; contrasts are not estimable.")
  }

  # Within-animal correlation across the three zones. See header.
  corfit <- duplicateCorrelation(expr_sex, design, block = P_sex$handle)
  cat("[step 16]", sex, "consensus within-animal correlation:",
      round(corfit$consensus.correlation, 4), "\n")

  fit <- lmFit(expr_sex, design, block = P_sex$handle,
               correlation = corfit$consensus.correlation)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  coef_names <- colnames(design)

  # Build one contrast per zone. Pericentral is the reference level, so its
  # aging effect is the main effect alone; the other zones add their
  # interaction coefficient.
  contrast_list <- list()

  contrast_list[["pericentral"]] <- setNames(rep(0, length(coef_names)), coef_names)
  contrast_list[["pericentral"]]["age_grpold"] <- 1

  contrast_list[["midlobular"]] <- setNames(rep(0, length(coef_names)), coef_names)
  contrast_list[["midlobular"]]["age_grpold"] <- 1
  contrast_list[["midlobular"]]["zonemidlobular:age_grpold"] <- 1

  contrast_list[["periportal"]] <- setNames(rep(0, length(coef_names)), coef_names)
  contrast_list[["periportal"]]["age_grpold"] <- 1
  contrast_list[["periportal"]]["zoneperiportal:age_grpold"] <- 1

  contrast_mat <- do.call(cbind, contrast_list)
  fit2 <- contrasts.fit(fit, contrast_mat)
  fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

  # topTable adjusts p-values within each contrast independently, so the FDR
  # is per zone, not across all three.
  for (i in seq_along(contrast_list)) {
    z <- names(contrast_list)[i]
    tt <- topTable(fit2, coef = i, number = Inf, sort.by = "P")
    tt$direction <- ifelse(tt$logFC > 0, "OLD", "YOUNG")
    outfile <- file.path(out_dir, paste0("DE_", z, "_", sex, "WT_old_vs_young.csv"))
    fwrite(as.data.table(tt, keep.rownames = "Gene"), outfile)
    cat("[step 16]   ", z, ":", sum(tt$adj.P.Val < FDR_THRESHOLD),
        "genes at FDR <", FDR_THRESHOLD, "\n")
  }
}

cat("\n[step 16] results ->", out_dir, "\n")
