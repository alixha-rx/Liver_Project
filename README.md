# Liver Project Spatial Data Pipeline
## By Alisha Aristel and Nancy R Zhang, PhD
PrePrint Paper: doi: https://doi.org/10.64898/2026.08.08.743614

Analysis code for the Visium spatial transcriptomics component of *Cellular
senescence is associated with age-related loss of liver zonation and hepatocyte
function*.

The pipeline runs in two parts. **Steps 1-5 (Python, Alisha Aristel)** take raw
Visium output through cropping, clipping, masking, normalization and zonation
scoring. **Steps 6-17 (R, Nancy R. Zhang)** run the downstream analysis and
produce the Visium panels of Figures 1-4.

---

### STEPS 1-5: preprocessing before analysis

`step_1_preprocessing_liver_project.ipynb` — loading, per-gene clipping at the
98th percentile, intensity-based pixel masking, square-root variance
stabilization and housekeeping normalization.

---

### STEPS 6-17: downstream analysis

All R. Run them in order with:

```
Rscript run_all_analysis_liver_project.R
```

or run any single step on its own — each sources what it needs.

| Step | Script | What it does | Figure |
|------|--------|--------------|--------|
| 6 | `step_6_setup_and_gene_lists` | Libraries, gene lists, gene categories | — |
| 7 | `step_7_load_samples` | Load samples, build sample sheet, harmonize genes | — |
| 8 | `step_8_zonation_binning` | Equal-count zonation bins; three normalizations | — |
| 9 | `step_9_zonation_marker_heatmaps` | Zone markers across age | **1c** |
| 10 | `step_10_per_gene_zone_heatmaps` | One zonation heatmap per gene | **2e-h, 4e** |
| 11 | `step_11_senescence_scoring` | Per-pixel senescence score | **3a-c, 3g** |
| 12 | `step_12_senescence_association` | Replicate correlation; spatial maps | **3d-f, 4a, 4c** |
| 13 | `step_13_pathway_enrichment` | GO / KEGG on senescence-associated genes | **4b, 4d** |
| 14 | `step_14_zone_senescence_de` | Senescence-high vs -low within each zone | **4f-h** |
| 15 | `step_15_pseudobulk_by_zone` | Collapse pixels to animal x zone | — |
| 16 | `step_16_limma_pipeline_E` | Zone-specific old-vs-young limma | **2b, 2d** |
| 17 | `step_17_de_summary_and_plots` | Summary tables, diagnostics, enrichment | **2c** |

Steps 15-17 are independent of 6-14. To reproduce only the Figure 2
differential expression results, run those three.

---

### Repository layout

```
config.R                          the only file you edit: paths and parameters
run_all_analysis_liver_project.R  runs steps 6-17 in order
step_6..17_*.R                    the analysis pipeline
R/                                shared function library
  zonation.R                        binning, three-zone split, normalization
  loading.R                         reading one sample; per-zone senescence DE
  plotting.R                        spatial maps, violins, percentile boxes
gene_lists/                       curated gene panels used by the analysis
```

### Setup

1. **R 4.4.1** or later.

2. Install dependencies:

   ```r
   install.packages(c("ggplot2", "pheatmap", "gplots", "gridExtra", "data.table",
                      "dplyr", "magrittr", "matrixStats", "colorspace", "fields",
                      "MBA", "pbapply", "tictoc", "GGally", "RColorBrewer"))

   if (!require("BiocManager")) install.packages("BiocManager")
   BiocManager::install(c("limma", "ComplexHeatmap", "clusterProfiler",
                          "org.Mm.eg.db", "enrichplot"))
   ```

3. Edit the four paths at the top of `config.R`.

4. Run.

Steps 13 and 17 call `enrichKEGG`, which queries the KEGG REST API and needs a
network connection.

### Data

The pixel tables are far too large for GitHub. All datasets generated in this
study are deposited in the SenNet consortium data repository,
https://data.sennetconsortium.org/ — sample IDs are listed in Supplementary
Table 7 of the paper. `config.R` documents the exact structure each step
expects.

### Runtime

Two to three hours end to end on a workstation, dominated by the senescence
scoring (step 11), the per-gene heatmaps (step 10) and the first pseudobulk
pass (step 15). Peak memory around 32 GB. The pseudobulk matrix and the zone
average arrays are cached after the first run, so steps 16 and 17 rerun in
under a minute.

---

### Zone orientation

Throughout, the continuous zonation score `zone31rqt` runs:

```
LOW  score  =  pericentral  =  zone 3  (central vein)
HIGH score  =  periportal   =  zone 1  (portal triad)
```

Step 15 asserts this against marker genes and stops if it does not hold. If you
adapt this code to new data, check that assertion first — an inverted zonation
score reverses the interpretation of every zone-specific result in the paper
while producing output that looks entirely normal.

### Known caveats

These are limitations of the analysis, recorded here so that anyone reusing the
code knows where to be careful.

**Female and male senescence scores are not on a common scale.** They are
computed from different gene panels against different reference distributions
(step 11). Compare within a sex, not across.

**Pixel-level p-values are not animal-level evidence.** The Wilcoxon tests in
step 14 treat pixels as independent observations. With tens of thousands of
pixels almost any difference is significant, so those asterisks indicate
direction and consistency, not population-level effect. Step 16 is the analysis
that tests at the level of the animal.

**Enrichment uses the whole-genome background,** not the 1,000-gene Visium
panel, which is itself enriched for immune and senescence genes. The p-values
in steps 13 and 17 rank pathways; they are not calibrated significance.

**FDR < 0.2** is permissive, chosen for an experiment with 2-4 animals per
group. Individual genes at that threshold are candidates for follow-up.

**Random seeds.** The original analysis did not set one. `set.seed(1)` is now
set in steps 11 and 12, so repeated runs agree with each other. The k-means
cluster *labels* in step 12 may differ from the original run; step 13 therefore
selects the positively-associated cluster by its mean correlation rather than
by its label.
