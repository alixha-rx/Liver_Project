# =============================================================================
# STEP 7: Load the preprocessed samples and harmonize genes across them
# =============================================================================
#
# By Nancy R. Zhang
#
# Reads one .RData file per sample, assembles the sample sheet (metatab), and
# reduces every sample to the set of genes present in ALL of them.
#
# WHY THE GENE HARMONIZATION MATTERS
#   The samples were processed in two batches, and the 1,000-gene panels built
#   in preprocessing are not identical between them. Any cross-sample
#   comparison therefore has to run on the intersection. This step computes
#   that intersection and subsets every sample to it, so that no later step has
#   to think about ragged gene sets. Expect the intersection to be noticeably
#   smaller than any individual panel -- the count is printed.
#
# INPUT
#   DATA_RDATA_DIRS/<handle>.RData   one list per sample (see config.R)
#
# OUTPUT (variables placed in the global environment)
#   samples   list of per-sample lists, each subset to the common gene set,
#             with an added $loc holding the x,y pixel coordinates
#   metatab   one row per sample: handle, gender, mutant, age, age_in_months,
#             npix, medpixtot, percmasked, description, plotcolor, plotlty
#   genes     character vector, the common gene set
#   nsamples  number of samples loaded
#
# A NOTE ON THE PLOT COLOURS
#   Colour encodes age and genotype, line type encodes sex. The scheme was set
#   in preprocessing and is carried through unchanged, so that a given sample
#   looks the same in every figure of the paper:
#     green    = young wildtype       purple = Ercc1 mutant
#     light -> dark red = 6-18mo, 18-26mo, >26mo wildtype
#     solid    = female               dashed = male
#
# RUN TIME: several minutes; the pixel tables are large.
# =============================================================================


# --- BLOCK 0: CONFIG ---------------------------------------------------------

source("step_6_setup_and_gene_lists_liver_project.R")


# --- BLOCK 1: FIND THE SAMPLE FILES ------------------------------------------
#
# Every handle is looked up in each directory of DATA_RDATA_DIRS in turn.
# Handles carry a suffix where preprocessing wrote a filtered version; both
# spellings are listed so that either is picked up.

handles_to_use <- c(
  "fwi3_filtered_withmaskpix", "fwi4_filtered_withmaskpix",
  "fwi_filtered_withmaskpix", "fwi2", "fwo", "fwo2", "fwy", "fwy3",
  "mwo_filtered_withmaskpix", "mwo3_filtered_withmaskpix",
  "mwy_filtered_withmaskpix", "mwo2", "mwy2"
)

rdata_files <- unlist(lapply(DATA_RDATA_DIRS, function(d) {
  fp <- file.path(d, paste0(handles_to_use, ".RData"))
  fp[file.exists(fp)]
}))
stopifnot(length(rdata_files) > 0)
cat("[step 7] found", length(rdata_files), ".RData files\n")


# --- BLOCK 2: LOAD -----------------------------------------------------------
#
# load_one_sample() (R/loading.R) validates each file and merges the intensity
# mask with the blood-vessel mask into a single logical maskpix.

samples <- lapply(rdata_files, load_one_sample)
nsamples <- length(samples)


# --- BLOCK 3: BUILD THE SAMPLE SHEET -----------------------------------------
#
# One row per sample, carrying the biological metadata plus three QC numbers
# (pixel count, median pixel intensity, fraction masked) that are worth a look
# before trusting any cross-sample comparison.

metatab <- data.frame(handle = rep("", nsamples),
                      gender = rep("", nsamples),
                      mutant = rep("", nsamples),
                      age = rep("", nsamples),
                      npix = rep(NA, nsamples),
                      medpixtot = rep(NA, nsamples))

for (i in 1:length(samples)) {
  metatab[i, 1:4] <- samples[[i]]$meta
  metatab$npix[i] <- nrow(samples[[i]]$datfm_inbox)
  metatab$medpixtot[i] <- median(samples[[i]]$pixtot)
  metatab$percmasked[i] <- sum(samples[[i]]$maskpix) / metatab$npix[i]
}

# Age arrives as a string ("4.5 months", "2.5 years"); convert to months.
temp <- strsplit(metatab$age, split = " ")
agenum <- unlist(lapply(temp, function(obj) { as.numeric(obj[1]) }))
agestr <- unlist(lapply(temp, function(obj) { obj[2] }))
sel <- which(agestr == "years")
agenum[sel] <- agenum[sel] * 12
metatab$age_in_months <- agenum

metatab$description <- paste(metatab$gender, ", ", metatab$mutant, ", ",
                             metatab$age_in_months, "mo", sep = "")


# --- BLOCK 4: PLOTTING KEY ---------------------------------------------------

plotcolor <- rep("chartreuse4", nsamples)                    # young wildtype
sel <- which(metatab$mutant == "mutant")
plotcolor[sel] <- "darkorchid2"                              # Ercc1 mutant
sel <- which(metatab$mutant == "wildtype" & metatab$age_in_months > 20)
plotcolor[sel] <- "indianred4"
sel <- which(metatab$mutant == "wildtype" & metatab$age_in_months > 6 &
             metatab$age_in_month <= 18)
plotcolor[sel] <- "indianred1"
sel <- which(metatab$mutant == "wildtype" & metatab$age_in_months > 18 &
             metatab$age_in_month <= 26)
plotcolor[sel] <- "indianred3"

plotlty <- ifelse(metatab$gender == "female", 1, 2)
plotlty[metatab$handle == "fwy2"] <- 4   # flagged as possibly problematic

metatab$plotcolor <- plotcolor
metatab$plotlty <- plotlty


# --- BLOCK 5: HARMONIZE GENES ACROSS SAMPLES ---------------------------------
#
# Take the intersection of the gene sets, then subset every sample to it. The
# hasgene check is belt-and-braces: after Reduce(intersect, ...) every gene is
# present everywhere by construction, but the counts are printed so that a
# surprising drop would be visible rather than silent.

ll <- vector("list", length(samples))
for (i in 1:length(samples)) {
  ll[[i]] <- samples[[i]]$genes
}
genes <- Reduce(intersect, ll)
cat("[step 7]", length(genes), "genes common to all", nsamples, "samples\n")

temp <- matrix(nrow = length(genes), ncol = length(samples))
rownames(temp) <- genes
colnames(temp) <- metatab$handle
for (i in 1:length(samples)) {
  temp[, i] <- genes %in% samples[[i]]$genes
}
hasgene <- rowSums(temp)
genes <- genes[hasgene == length(samples)]


# --- BLOCK 6: ASSIGN EACH GENE A CATEGORY ------------------------------------
#
# Order matters: a gene matching several lists keeps the LAST category
# assigned. Ccnd1, for instance, is both a midlobular marker and a cell-cycle
# gene, and ends up categorized as Cell_Cycle.

gene_category <- rep(NA, length(genes))
gene_category[genes %in% genes_pericentral]    <- "Pericentral"
gene_category[genes %in% genes_midlobular]     <- "Midlobular"
gene_category[genes %in% genes_periportal]     <- "Periportal"
gene_category[genes %in% genes_housekeeping]   <- "Housekeeping"
gene_category[genes %in% genes_fibroblast]     <- "Fibroblast"
gene_category[genes %in% genes_cellcycle]      <- "Cell_Cycle"
gene_category[genes %in% genes_sasp]           <- "SASP"
gene_category[genes %in% genes_antiapoptosis]  <- "Anti_Apoptosis"
gene_category[genes %in% genes_endothelial]    <- "Endothelial"
gene_category[genes %in% genes_immune]         <- "Immune"
gene_category[genes %in% genes_downsenescent]  <- "Down_in_Senescent"


# --- BLOCK 7: REDUCE EACH SAMPLE TO THE COMMON GENE SET ----------------------
#
# Pull the x,y coordinates into $loc before subsetting, since the first two
# columns of datfm_inbox are coordinates rather than genes. This cuts memory by
# roughly a quarter.

for (i in 1:length(samples)) {
  samples[[i]]$loc <- samples[[i]]$datfm_inbox[, c(1, 2)]
  samples[[i]]$datfm_inbox <- samples[[i]]$datfm_inbox[, genes]
  samples[[i]]$genes <- genes
}
gc()

cat("[step 7] loaded", nsamples, "samples x", length(genes), "genes\n")
print(metatab[, c("handle", "gender", "mutant", "age_in_months", "npix", "percmasked")])
