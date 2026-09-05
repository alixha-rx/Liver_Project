# =============================================================================
# run_all_analysis_liver_project.R -- Run steps 6 through 17 in order.
# =============================================================================
#
# Usage, from the repository root:
#     Rscript run_all_analysis_liver_project.R
#
# Edit config.R first to point at your copy of the data.
#
# Each step can also be run on its own; every one begins by sourcing the steps
# it depends on if their results are not already in the environment. Running
# them all in one session is faster, because the pixel tables are then loaded
# once rather than repeatedly.
#
# Steps 15-17 (pseudobulk and limma) do NOT depend on steps 6-14. If you only
# want the Figure 2 differential expression results, run those three.
#
# TOTAL RUN TIME: two to three hours, dominated by step 11 (senescence
# scoring), step 10 (one heatmap per gene) and step 15 (first pseudobulk pass).
# Peak memory is around 32 GB; the pixel tables are the reason.
# =============================================================================

steps <- c(
  "step_6_setup_and_gene_lists_liver_project.R",
  "step_7_load_samples_liver_project.R",
  "step_8_zonation_binning_liver_project.R",
  "step_9_zonation_marker_heatmaps_liver_project.R",
  "step_10_per_gene_zone_heatmaps_liver_project.R",
  "step_11_senescence_scoring_liver_project.R",
  "step_12_senescence_association_liver_project.R",
  "step_13_pathway_enrichment_liver_project.R",
  "step_14_zone_senescence_de_liver_project.R",
  "step_15_pseudobulk_by_zone_liver_project.R",
  "step_16_limma_pipeline_E_liver_project.R",
  "step_17_de_summary_and_plots_liver_project.R"
)

t_start <- Sys.time()
for (s in steps) {
  cat("\n", strrep("=", 70), "\n>>> ", s, "\n", strrep("=", 70), "\n", sep = "")
  t0 <- Sys.time()
  source(s)
  cat(">>> ", s, " finished in ",
      format(round(difftime(Sys.time(), t0), 1)), "\n", sep = "")
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Pipeline complete in ", format(round(difftime(Sys.time(), t_start), 1)), "\n", sep = "")
cat("Results: ", RESULTS_DIR, "\n", sep = "")
cat("Figures: ", FIGURES_DIR, "\n", sep = "")

cat("\nSession info:\n")
print(sessionInfo())
