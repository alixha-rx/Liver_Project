# =============================================================================
# R/zonation.R -- Turning a continuous zonation score into ordered zone bins.
# =============================================================================
#
# The zonation score (`zone31rqt`, computed in preprocessing) is a continuous
# pseudotime along the porto-central axis of the liver lobule:
#
#     LOW  score  =  pericentral  =  zone 3  (around the central vein)
#     HIGH score  =  periportal   =  zone 1  (around the portal triad)
#
# Every zone-resolved analysis in this repository works by ranking pixels on
# this score and cutting the ranks into groups of EQUAL PIXEL COUNT. Equal
# count, rather than equal score width, is what makes bins comparable between
# samples: each bin then represents the same fraction of the tissue regardless
# of how the score happens to be distributed in that particular section.
# =============================================================================


#' Assign ranks to equal-count bins
#'
#' @param ranks  integer vector of ranks, as returned by `rank()`. Must run
#'               from 1 to length(ranks).
#' @param nbins  number of bins to produce.
#' @return integer vector, same length as `ranks`, giving the bin index 1..nbins.
#'
#' Each bin receives floor(n / nbins) elements; the final bin absorbs the
#' remainder, so it can be up to nbins-1 elements larger than the others.
#' Bin 1 = lowest ranks = pericentral. Bin nbins = highest = periportal.
bin <- function(ranks, nbins = 20) {
  n <- length(ranks)
  num_per_bin <- floor(n / nbins)
  bin_id <- rep(NA, n)
  for (i in 1:nbins) {
    bin_start <- (i - 1) * (num_per_bin) + 1
    bin_end <- i * num_per_bin
    if (i == nbins) bin_end <- n
    which_in_bin <- which((ranks >= bin_start) & (ranks <= bin_end))
    bin_id[which_in_bin] <- i
  }
  bin_id
}


#' Average expression within each zonation bin
#'
#' @param dat    pixels x genes matrix or data.frame of normalized expression.
#' @param zone   numeric vector of zonation scores, one per row of `dat`.
#' @param nbins  number of equal-count bins.
#' @param normalize_per_gene  if TRUE, divide each gene by its own maximum
#'               across bins, putting every gene on a 0-1 scale. Use FALSE when
#'               bins will later be normalized against a common reference
#'               (which is what steps 8 and 9 do).
#' @return nbins x genes matrix. Row 1 = pericentral, row nbins = periportal.
#'
#' Pixels must already be masked before calling: this function has no notion of
#' `maskpix` and will happily average over pixels you meant to drop.
compute_bin_average <- function(dat, zone, nbins, normalize_per_gene = FALSE) {
  zone_bin <- bin(rank(zone), nbins = nbins)
  bin_avg <- matrix(nrow = nbins, ncol = ncol(dat))
  colnames(bin_avg) <- colnames(dat)
  for (i in 1:nbins) {
    ids <- which(zone_bin == i)
    bin_avg[i, ] <- colMeans(dat[ids, , drop = FALSE])
  }
  if (normalize_per_gene) {
    gene_max <- apply(bin_avg, 2, max)
    bin_avg <- sweep(bin_avg, 2, gene_max, "/")
  }
  bin_avg
}


#' Cut pixels into three equal-count hepatic zones
#'
#' @param zone_score  numeric vector of zonation scores (already masked).
#' @return character vector of zone labels, same length as `zone_score`.
#'
#' The bottom third of the ranked score is pericentral, the top third
#' periportal, the remainder midlobular. This is the split used for the
#' pseudobulk matrix that feeds limma (steps 14-15).
#'
#' The orientation here is the one thing in this file worth checking if you
#' adapt it to new data: confirm that pixels labelled "pericentral" really do
#' have higher Cyp2e1 / Oat / Gstm3 / Cyp1a2 than those labelled "periportal",
#' and the reverse for Cyp2f2 / Cps1 / Pck1 / Asl / Gls2.
assign_three_zones <- function(zone_score) {
  zone_rank <- rank(zone_score) / length(zone_score)
  zone_label <- rep("midlobular", length(zone_rank))
  zone_label[zone_rank <= 1/3] <- "pericentral"
  zone_label[zone_rank >= 2/3] <- "periportal"
  zone_label
}


#' Square-root variance stabilization + housekeeping normalization
#'
#' @param mat         pixels x genes matrix of counts.
#' @param transformf  variance-stabilizing transform. sqrt is used throughout;
#'                    it is gentler than log for Poisson-like counts, where
#'                    variance tracks the mean.
#' @param scalegenes  if TRUE, additionally divide each gene by its median.
#' @param ctrgenes    housekeeping gene names; each pixel is divided by the sum
#'                    of these genes in that pixel, which removes per-pixel
#'                    differences in total capture efficiency.
#' @param ctrval      alternatively, supply the per-pixel denominator directly.
#'
#' Exactly one of `ctrgenes` or `ctrval` must be given.
normalizeData <- function(mat, transformf = sqrt, scalegenes = TRUE,
                          ctrgenes = NULL, ctrval = NULL) {
  mat <- as.matrix(mat)
  mat <- transformf(mat)

  if (!is.null(ctrgenes) & is.null(ctrval)) {
    ctrval <- rowSums(mat[, ctrgenes], na.rm = TRUE)
    mat <- sweep(mat, 1, ctrval, FUN = "/")
  } else if (!is.null(ctrval)) {
    ctrval <- transformf(ctrval)
    mat <- sweep(mat, 1, ctrval, FUN = "/")
  } else {
    stop("Exactly one of ctrgenes or ctrval has to be non-null.")
  }
  if (scalegenes) {
    mat <- scale(mat, center = FALSE, scale = matrixStats::colMedians(mat, na.rm = TRUE))
  }
  mat
}
