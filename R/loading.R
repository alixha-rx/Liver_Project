# =============================================================================
# R/loading.R -- Reading one preprocessed sample, and the per-zone senescence
#                correlation test.
# =============================================================================
#
# CONTENTS
#   load_one_sample()                Read one .RData file produced by the
#                                    preprocessing steps and validate it.
#   compute_corr_senescence_per_zone() For one sample, split pixels into three
#                                    zones and test, within each zone, which
#                                    genes are differentially expressed between
#                                    senescence-high and senescence-low pixels
#                                    (Fig 4f-h).
#
# WHAT load_one_sample() GUARANTEES
#   It fails loudly rather than silently returning a half-populated object.
#   After it returns, the sample is known to have datfm_inbox, zone31rqt,
#   genes and maskpix, and meta with handle/gender/mutant/age.
#
#   It also merges the two mask fields. Preprocessing may write `maskpix`
#   (intensity outliers) and `maskpix_bloodvessels` (large vessels) separately;
#   this function ORs them into a single `maskpix`, so downstream code has one
#   mask to think about and a pixel flagged by either rule is dropped.
# =============================================================================


load_one_sample <- function(f) {
  e <- new.env()
  obj_names <- load(f, envir = e)

  # detect  s or sample_list
  candidate <- c("s", "sample_list", "sample")
  name <- candidate[candidate %in% obj_names][1]
  if (is.na(name)) stop("Cannot detect main object in ", f)

  s <- get(name, envir = e)

  # check meta structure
  if (!"meta" %in% names(s))
    stop("Missing 'meta' field in ", f)
  if (!all(c("handle", "gender", "mutant", "age") %in% names(s$meta)))
    stop("meta must have handle, gender, mutant, age in ", f)

  # handle the new maskpix_bloodvessels field if present
  if ("maskpix_bloodvessels" %in% names(s)) {
    # combine masks so either TRUE = masked pixel
    if ("maskpix" %in% names(s)) {
      s$maskpix <- s$maskpix | s$maskpix_bloodvessels
    } else {
      s$maskpix <- s$maskpix_bloodvessels
    }
  } else if (!"maskpix" %in% names(s)) {
    # if both missing, initialize empty logical
    s$maskpix <- logical(0)
  }

  # sanity check this is to make sure required things are present
  needed <- c("datfm_inbox", "zone31rqt", "genes", "maskpix")
  missing <- setdiff(needed, names(s))
  if (length(missing) > 0)
    stop("Missing fields in ", f, ": ", paste(missing, collapse = ", "))

  s
}


compute_corr_senescence_per_zone<-function(dfbig, gene_cols, zone_specific_corr_outdir, handlestr){

  # define 3 zones by tertiles within this sample
  q1 <- quantile(dfbig$zone31rqt, 1/3, na.rm=TRUE)
  q2 <- quantile(dfbig$zone31rqt, 2/3, na.rm=TRUE)
  dfbig[, zone3 := fifelse(zone31rqt < q1, "periportal",
                           fifelse(zone31rqt < q2, "midlobular", "pericentral"))]
  dfbig[, zone3 := factor(zone3, levels=c("periportal","midlobular","pericentral"))]

  # bucket per zone: top 5% vs bottom 50%
  dfbig[, `:=`(
    t5 = quantile(senscore, 0.95, na.rm=TRUE),
    b5 = quantile(senscore, 0.50, na.rm=TRUE)
  ), by = zone3]
  dfbig[, group := fifelse(senscore >= t5, "top5",
                           fifelse(senscore <= b5, "bottom50", NA_character_)),
        by = zone3]

  dfbig <- dfbig[!is.na(group)]
  if (!nrow(dfbig)) {
    cat(sprintf("[SKIP] %s: no pixels in top5/bottom50 after bucketing\n", handlestr))
    next
  }

  # target genes = SASP ∪ AntiApop present in this sample
  gs_target <- intersect(unique(c(genes_sasp, genes_antiapoptosis)), gene_cols)

  zones <- levels(dfbig$zone3)
  res_list <- vector("list", length(zones))

  # do each zone and bind for this sample
  for (z in zones) {
    cat("--- Processing zone ", z, "\n")
    sub <- dfbig[zone3 == z & group %in% c("top5","bottom50")]
    if (!nrow(sub)) next

    # this is the gs for only the scene genes
    #gs <- intersect(gs_target, colnames(sub))
    # this is everything thats in the columns of sub the genes
    gs  <- intersect(gene_cols, colnames(sub))

    if (!length(gs)) next

    # genes x pixels matrix
    # X <- t(as.matrix(sub[, ..gs])) this was the orginal
    # X <- t(as.matrix(sub[, ..gs, with = FALSE]))
    # G <- nrow(X); P <- ncol(X)

    cat("Getting X matrix consisting of gene cols ...\n")
    X <- sub[, ..gs, with = FALSE]  #### NOTE that this is different from scene_score_nancy_2.R, which uses the above line.
    G <- ncol(X); P <- nrow(X)


    # normalize each pixel to library size.
    cat("Scaling X by pixtot ...\n")
    X_scaled = sweep(X, 1, sub$pixtot, "/")

    idx_top <- which(sub$group == "top5")
    idx_bot <- which(sub$group == "bottom50")

    #plot(sub$x[idx_top], sub$y[idx_top])  # there is "clumping"
    # NRZ: subsample pixels, this alleviates "clumping"
    #idx_top_sel = sample(idx_top, 100, replace=FALSE)
    #plot(sub$x[idx_top_sel], sub$y[idx_top_sel])  # is there still "clumping"
    #idx_top=idx_top_sel

    n_top <- length(idx_top); n_bot <- length(idx_bot)
    if (n_top < 2 || n_bot < 2) next

    # observed effect

    # NRZ: match the top and bottom by zonation
    zone_top = sub$zone31rqt[idx_top]
    zone_bot = sub$zone31rqt[idx_bot]
    dist_mat <- abs(outer(zone_top, zone_bot, "-"))
    nearest_idx <- max.col(-dist_mat)          # index of min distance
    nearest_dist <- dist_mat[cbind(seq_along(zone_top), nearest_idx)]

    ## ---- enrich with DE info ----
    #    cat("Computing per gene means and sds for top and bottom senescence brackets ...\n")
    #    mu_top <- colMeans(X[idx_top, ,drop=FALSE], na.rm=TRUE)
    #    mu_bot <- colMeans(X[idx_bot, ,drop=FALSE], na.rm=TRUE)
    #    sd_top <- apply(X[idx_top, ,drop=FALSE], 2, sd, na.rm=TRUE)
    #    sd_bot <- apply(X[idx_bot, ,drop=FALSE], 2, sd, na.rm=TRUE)

    #    diff_mean <- mu_top - mu_bot
    #    log2FC    <- log2((mu_top + 1e-8) / (mu_bot + 1e-8))
    #    sd_pool <- sqrt( ((n_top - 1) * sd_top^2 + (n_bot - 1) * sd_bot^2) /
    #                      max(n_top + n_bot - 2, 1) )
    #    d_cohen <- diff_mean / (sd_pool + 1e-12)


    # ---- NRZ: try t test and wilcoxon on the zone-matched ---------
    # mu1, mu2, s1, s2, s are essentially what Alisha computed on top,
    # but with zone-matched nearest_idx, instead of idx_bot.
    cat("Computing per gene means and sds for top and bottom senescence brackets ...\n")
    mu1 = apply(X_scaled[idx_top, , drop=FALSE], 2, mean, na.rm=TRUE)
    mu2 = apply(X_scaled[nearest_idx, , drop=FALSE], 2, mean, na.rm=TRUE)
    s1 = apply(X_scaled[idx_top, , drop=FALSE], 2, sd, na.rm=TRUE)
    s2 = apply(X_scaled[nearest_idx, , drop=FALSE], 2, sd, na.rm=TRUE)
    s_pool = sqrt(((s1^2)*(n_top-1)+(s2^2)*(n_top-1))/(n_top+n_top-2))
    tstat = (mu1-mu2)/(s_pool*sqrt(1/n_top + 1/n_top))
    dof=n_top+n_top-2  # since we matched top and bottom, there are only n_top points in bottom.
    p_t <- 2 * (1 - pt(abs(tstat), dof))
    log2FC_scaled    <- log2((mu1 + 1e-8) / (mu2 + 1e-8))
    d_cohen_scaled <- (mu1-mu2) / (s_pool + 1e-12)


    # NRZ: try wilcoxon rank-sum
    cat("Computing per gene wilcoxon ...\n")
    p_wilcoxon <- unlist(pblapply(seq_len(ncol(X_scaled)), function(i)
      wilcox.test(X_scaled[idx_top,i ], X_scaled[nearest_idx,i], paired = TRUE)$p.value
    ))

    # NRZ: diagnostic plots:
    cat("Outputing plots ...\n")
    png(sprintf("DE_perm_OLD_%s_%s.png", handlestr, z), height=800, width=800)
    par(mfrow=c(2,2))
    #hist(tstat)
    idx_top_sel = sample(idx_top, 100, replace=FALSE)
    plot(sub$x[idx_top_sel], sub$y[idx_top_sel])  # is there still "clumping"
    hist(p_wilcoxon, 100)
    # t and wilcoxon should roughly agree
    xmax=max(sqrt(-log(p_t)), na.rm=TRUE)
    if(!is.finite(xmax)) xmax=20
    ymax=max(sqrt(-log(p_wilcoxon)), na.rm=TRUE)
    if(!is.finite(ymax)) ymax=20

    plot(sqrt(-log(p_t)), sqrt(-log(p_wilcoxon)), xlim=c(0, xmax), ylim=c(0,ymax))
    # log fold change versus wilcoxon's p
    plot(log((mu1+ 1e-8)/(mu2+ 1e-8)), sqrt(-log(p_wilcoxon)))
    #sel_genes = which(abs(log(mu1/mu2))>0.15 & unlist(p_wilcoxon)<0.05/nrow(X))
    dev.off()


    # get_auc <- function(v_top, v_bot) {
    #   wt <- suppressWarnings(wilcox.test(v_top, v_bot, exact = FALSE))
    #   U  <- as.numeric(wt$statistic)
    #   U / (length(v_top) * length(v_bot))
    # }
    # AUC <- pblapply(seq_len(nrow(X)),
    #               function(i) get_auc(X[i, idx_top], X[i, idx_bot]))

    cat("Writing results ...\n")
    tab <- data.frame(
      sample    = handlestr,
      gender    = gender,
      age_raw   = age,
      zone      = z,
      Gene      = colnames(X),
      #      mean_top  = mu_top,
      mean_top = mu1,
      #      mean_bot  = mu_bot,
      mean_bot  = mu2,
      #      diff_mean = diff_mean,
      #      log2FC    = log2FC,
      #      sd_top    = sd_top,
      #      sd_bot    = sd_bot,
      log2FC = log2FC_scaled,
      sd_top = s1,
      sd_bot = s2,
      #      d_cohen   = d_cohen,
      d_cohen   = d_cohen_scaled,
      #      AUC       = AUC,
      #      p_perm    = p_perm,
      mean_bot_matched = mu2,
      log2FC_matched = log((mu1+ 1e-8)/(mu2+ 1e-8)),
      p_wilcoxon = p_wilcoxon,
      p_t = p_t,
      stringsAsFactors = FALSE
    )
    # NRZ substitute wilcoxon for perm
    #    tab$padj_BH <- p.adjust(tab$p_perm, method = "BH")
    tab$pwilcox_adj_BH <- p.adjust(tab$p_wilcoxon, method = "BH")
    tab$direction <- ifelse(tab$log2FC_matched > 0, "up_in_top5", "down_in_top5")
    tab$type <- ifelse(tab$Gene %in%  genes_sasp, "SASP", "Others")

    ## ---- write one TSV per sample × zone ----
    outfile_zone <- file.path(zone_specific_corr_outdir, sprintf("DE_perm_OLD_%s_%s.tsv", handlestr, z))
    fwrite(tab, outfile_zone, sep = "\t")
    cat("[OK] Wrote:", outfile_zone, "\n")
    wrote_any <- TRUE
  }  # end of for(z in zones)

  return(res_list)
} # end of compute_corr_senescence_per_zone<- function(...)
