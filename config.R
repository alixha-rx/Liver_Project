# =============================================================================
# config.R -- The ONE file you need to edit before running steps 6-16.
# =============================================================================
#
# Every analysis script begins with `source("config.R")`. Nothing else in the
# repository contains a machine-specific path, so pointing the four variables
# below at your own copy of the data is the whole of the setup.
#
# WHAT THE ANALYSIS EXPECTS AS INPUT
# ----------------------------------
# Steps 6-16 start from the output of the preprocessing notebooks (steps 1-5),
# not from raw Visium output. Two forms of that output are used:
#
#   (A) One .RData file per sample, each holding a single list `s` with fields:
#         s$datfm_inbox  data.frame, pixels x genes. Columns 1-2 are x, y;
#                        the remaining columns are normalized gene expression.
#         s$zone31rqt    numeric, one value per pixel. The continuous zonation
#                        score, quantile-normalized to a common 0-1 scale.
#                        LOW = pericentral (zone 3), HIGH = periportal (zone 1).
#         s$maskpix      logical, one per pixel. TRUE = drop this pixel
#                        (intensity outlier or large blood vessel).
#         s$genes        character, the gene names.
#         s$meta         list(handle, gender, mutant, age)
#       Used by steps 6-13 (zonation binning, senescence scoring, spatial maps).
#
#   (B) One .rds file per sample with the same fields, plus a metatab.csv
#       sample sheet with columns handle, gender, mutant, age_months.
#       Used by steps 14-16 (pseudobulk + limma).
#
#   (A) and (B) hold the same measurements in two different containers; the
#   split is historical and is preserved here so that published results
#   reproduce exactly. If you are starting fresh, write both from step 5.
#
# The pixel tables are far too large for GitHub. They are deposited with the
# rest of the study data in the SenNet consortium data repository:
# https://data.sennetconsortium.org/  (sample IDs in Supplementary Table 7).
# =============================================================================


# --- 1. EDIT THESE FOUR PATHS ------------------------------------------------

# Directory holding the per-sample .RData files from preprocessing (form A).
# May be a character vector if the samples are split over several folders --
# every folder in the vector is searched, and each handle is taken from the
# first folder in which it is found.
DATA_RDATA_DIRS <- c(
  "~/liver_project_data/rdata_by_handle_batch1",
  "~/liver_project_data/rdata_by_handle_batch2"
)

# Directory holding the per-sample .rds files and metatab.csv (form B).
DATA_RDS_DIR <- "~/liver_project_data/rds_by_handle"

# Directory holding the curated gene lists that ship with this repo
# (mfuzz.csv, zonation_genes.txt, senescence_genes_20250307.txt, ...).
GENE_LISTS_DIR <- file.path(getwd(), "gene_lists")

# Directory holding the per-animal, per-zone senescence DE tables consumed by
# step 14 (files named DE_perm_OLD_<handle>_<zone>.tsv). These are an
# intermediate product; step 14 can also regenerate them from the full
# 5,000-gene superpixel tables. Only step 14 uses this.
DE_TABLES_DIR <- "~/liver_project_data/de_perm_tables"

# Where results and figures are written. Created if it does not exist.
OUT_DIR <- file.path(getwd(), "output")


# --- 2. DERIVED PATHS (no need to edit) --------------------------------------

RESULTS_DIR <- file.path(OUT_DIR, "results")   # tables: DE CSVs, correlations
FIGURES_DIR <- file.path(OUT_DIR, "figures")   # PNG/PDF panels
CACHE_DIR   <- file.path(OUT_DIR, "cache")     # intermediates (pseudobulk, ...)

for (d in c(OUT_DIR, RESULTS_DIR, FIGURES_DIR, CACHE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# --- 3. ANALYSIS PARAMETERS --------------------------------------------------
# These are the values used for the published analysis. They are collected
# here, rather than buried in the scripts, so that the choices are auditable.

# Masking. TRUE = exclude pixels flagged in preprocessing (intensity outliers
# and large blood vessels) from every downstream computation.
DO_MASK <- TRUE

# Zonation binning. The continuous zonation score is turned into ordered bins
# of equal pixel count. Two resolutions are used:
NBINS_COARSE <- 20    # violin plots and per-gene heatmaps (Fig 2e-h, 3c, 4e)
NBINS_FINE   <- 100   # zonation marker heatmaps (Fig 1c)

# Three-zone split for pseudobulk. Pixels are ranked by zonation score and cut
# into equal thirds. Because LOW score = pericentral, the bottom third is
# pericentral and the top third is periportal.
ZONES <- c("pericentral", "midlobular", "periportal")

# Senescence score exponents. See step 10 for the formula.
ALPHA_SASP           <- 1
ALPHA_CELL_CYCLE     <- 2
ALPHA_ANTIAPOPTOSIS  <- 1

# Senescence-high pixel thresholds, defined as the 98th percentile of the
# senescence score in the oldest age group of each sex (24mo female, 30mo male)
# and then applied to every sample of that sex. Used for Fig 3d and 3h.
SENSCORE_THRESH_FEMALE <- 0.86
SENSCORE_THRESH_MALE   <- 0.41

# Differential expression.
FDR_THRESHOLD  <- 0.2    # significance cutoff used throughout the paper
GENE_FILTER_Q  <- 0.1    # drop the lowest-expressed 10% of genes before limma

# Age grouping for the limma comparison (months).
AGE_YOUNG_MAX  <- 18     # young = strictly below this
# Female samples aged >= FEMALE_AGE_EXCLUDE are dropped (see step 15).
FEMALE_AGE_EXCLUDE <- 25

# Samples excluded from the limma analysis for quality reasons.
EXCLUDE_HANDLES <- c("fwy3")


# --- 4. LOAD THE SHARED FUNCTION LIBRARY -------------------------------------

for (f in list.files(file.path(getwd(), "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(f)
}

cat("[config] output ->", OUT_DIR, "\n")
