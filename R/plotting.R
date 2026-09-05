# =============================================================================
# R/plotting.R -- Figure functions shared by steps 9-13.
# =============================================================================
#
# Every function here returns a ggplot object; nothing is written to disk from
# inside these functions. The calling step decides the filename and size, so
# that all output paths live in the step scripts and none are buried here.
#
# CONTENTS
#   spatial_map_value()             Paint any per-pixel vector onto tissue
#                                   coordinates. Used for senescence score maps
#                                   (Fig 3d, 3h) and for zonation maps.
#   spatial_map_gene()              Same, for one gene's expression
#                                   (Fig 3e-f, 3i-j: Cyp2f2 and Cyp2e1).
#   plot_senescence_by_zone_violin() Senescence score distribution across
#                                   zonation bins, one violin per bin, with
#                                   median / 95th / 99th percentile traces
#                                   (Fig 3c, 3g).
#   plot_faceted_percentile_boxes() Per-sample senescence score box plots with
#                                   90/95/99th percentile markers (Fig 3a, 3b).
#   plot_faceted_violins()          Per-sample senescence violins grouped by
#                                   age category.
#
# ORIENTATION NOTE
#   plot_senescence_by_zone_violin() REVERSES the bin numbering, so that its
#   x-axis reads periportal (left) -> pericentral (right). This is the opposite
#   of the raw `bin()` numbering in R/zonation.R, where bin 1 is pericentral.
#   The reversal is deliberate and matches the published figure orientation;
#   the axis label states the direction explicitly.
#
# The blood-vessel masking function that used to live in this file has been
# left with the preprocessing code (steps 1-5), where the mask is actually
# computed. Steps 6-16 consume `s$maskpix` and never recompute it.
# =============================================================================

spatial_map_value=function(vec, coords, label="value",
                           color_theme="plasma",points_to_highlight=NULL,  col_highlight="black"){
  temp=cbind(coords, vec)
  names(temp) = c("x","y", label)
  temp=as.data.frame(temp)
#  min_value=0
#  if(is.null(max_value)) max_value = max(vec, na.rm=TRUE)
  p1=ggplot(temp[,], aes_string(x="x", y="y", fill=label))+
    geom_raster()+
    theme_minimal() + theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank())
  if(!is.null(points_to_highlight)){
    p1=p1+geom_point(data=points_to_highlight, aes(x=x, y=y),
                     color=col_highlight,shape=22,  size=0.5, fill=col_highlight)

  }
  if(length(color_theme)>1){
    p1 = p1+scale_fill_gradientn(colors=color_theme)
  } else {
    p1 = p1+scale_fill_viridis_c(option=color_theme)
  }
  p1
}

spatial_map_gene=function(genename, dat, coords, color_theme="cividis",points_to_highlight=NULL,  col_highlight="black"){
  temp=cbind(coords, dat[,genename])
  names(temp) = c("x","y", genename)
  temp=as.data.frame(temp)
  pgene=ggplot(temp[,], aes_string(x="x", y="y", fill=genename))+
    geom_raster()+
    scale_fill_viridis_c(option="cividis")+
    theme_minimal() + theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank())
  if(!is.null(points_to_highlight)){
    pgene=pgene+geom_point(data=points_to_highlight, aes(x=x, y=y),
                     color=col_highlight,shape=22,  size=0.5, fill=col_highlight)

  }

  pgene
}

plot_senescence_by_zone_violin=function(zonescore, senscore, nbreaks=20, title=NULL){
  # `title` defaults to the global metatab/handle lookup that the original
  # analysis relied on, so behaviour is unchanged when called from step 10.
  # Pass it explicitly to use this function outside that loop.
  if (is.null(title)) title <- metatab$description[handle]
  df=data.frame(zonescore, senscore)
  # Reverse bin numbering so bin 1 = periportal, bin N = pericentral
  df <- df %>%
    mutate(x_bin = nbreaks + 1L - cut(zonescore, breaks = nbreaks, labels = FALSE))
  med_df <- df %>%
    group_by(x_bin) %>%
    summarise(
      x_mid = mean(range(zonescore)),    # midpoint of the bin
      y_med = median(senscore),
      .groups = "drop"
    )
  qt_df <- df %>%
    group_by(x_bin) %>%
    summarise(
      x_mid = mean(range(zonescore)),    # midpoint of the bin
      y_qt = quantile(senscore,0.95),
      y_qt99 = quantile(senscore, 0.99),
      .groups = "drop"
    )

  p1=ggplot(df, aes(x = factor(x_bin), y = senscore)) +
    geom_violin(scale = "width",fill = "lightblue") +
    geom_line(data = med_df, aes(x = x_bin, y = y_med),
              inherit.aes = FALSE, color = "blue", linewidth = 1) +
    geom_line(data = qt_df, aes(x = x_bin, y = y_qt),
              inherit.aes = FALSE, color = "orange", linewidth = 1) +
    geom_line(data = qt_df, aes(x = x_bin, y = y_qt99),
              inherit.aes = FALSE, color = "red", linewidth = 1) +
    labs(x = "Zonation periportal ----> pericentral", y = "Senescence score", title=title) +
    theme_bw() +  # <-- white background with grid lines
    theme(axis.text.x = element_blank())+  # hide cluttered x-axis labels
    scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)   # optional, removes extra padding
    )

    return(p1)
}




plot_faceted_percentile_boxes <- function(senscore_list, metatab, gender="female", do_transform=FALSE, transform=NULL){

  #  long-form data frame
  names(senscore_list)= metatab$handle
  sel= which(metatab$gender==gender & metatab$mutant=="wildtype" & metatab$handle != "fwy2")

  age_in_months = metatab$age_in_months[sel]
  mylist = senscore_list[sel]

  # Assign age categories
  categories = rep(0, length(sel))
  categories[age_in_months < 6] = "<6 months"
  categories[age_in_months == 18] = "18 months"
  categories[age_in_months == 24] = "24 months"
  categories[age_in_months > 28] = ">28 months"


  df <- data.frame(
    ind = rep(names(mylist), times = sapply(mylist, length)),
    values = unlist(mylist)
  )

  if(!is.null(transform)){
    df$values = transform(df$values)
  }

  # optional transform
  if (do_transform) {
    df$values <- sqrt(df$values)
  }


  df$category <- categories[match(df$ind, names(mylist))]

  df$category <- factor(
    df$category,
    levels = c("<6 months", "18 months", "24 months", ">28 months"))



  # Compute quantiles (median, 90%, 95%, 99%)
  qdf <- df %>%
    group_by(ind, category) %>%
    summarise(
      med = median(values),
      q90 = quantile(values, 0.90),
      q95 = quantile(values, 0.95),
      q99 = quantile(values, 0.99),
      .groups="drop"
    )

  #ymax <- max(df$values, na.rm = TRUE)
  #print(paste("YMAX for gender", gender, "=", ymax))
  ymax <- ifelse(gender == "female", 1, 1)






  # Plot
  p1 = ggplot(df, aes(x = ind, y = values)) +
    geom_boxplot(outlier.shape = NA, width=0.5, fill="lightgray") +
#    geom_point(data = qdf, aes(x = ind, y = med, color="50%"), shape="-", size=10) +
    geom_point(data = qdf, aes(x = ind, y = q90, color="90%"), shape="-",  size=10) +
    geom_point(data = qdf, aes(x = ind, y = q95, color="95%"),  shape="-", size=10) +
    geom_point(data = qdf, aes(x = ind, y = q99, color="99%"),  shape="-", size=10) +
    scale_color_manual(
      name = "Percentile",
      values = c("50%"="blue", "90%"="orange", "95%"="red", "99%"="purple")
    ) +  scale_y_continuous(limits = c(0, ymax), expand = c(0,0)) +
    facet_wrap(~category, scales="free_x", nrow=1) +
    theme_minimal() +
    labs(x = "", y = "Senescence Score") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
  #annotate("text", x=1, y=1.02, label="med=blue, 90%=orange, 95%=red, 99%=purple", size=3)

  return(p1)
}



plot_faceted_violins=function(senscore_list, metatab, gender="female", do_transform=FALSE){
  # Convert list to long-form data frame
  names(senscore_list)= metatab$handle
  sel= which(metatab$gender==gender & metatab$mutant=="wildtype" & metatab$handle != "fwy2")
  age_in_months = metatab$age_in_months[sel]
  mylist = senscore_list[sel]
  categories = rep(0, length(sel))
  categories[age_in_months<6] = "<6 months"
  categories[age_in_months==18] = "18 months"
  categories[age_in_months==24] = "24 months"
  categories[age_in_months>28] = ">28 months"

  df <- stack(mylist)
  if(do_transform){df$values = sqrt(df$values)}
  df$category <- categories[df$ind]  # add category info

  df$category <- factor(df$category, levels = c("<6 months", "18 months", "24 months", ">28 months"))

  # Plot grouped violins
  p1=ggplot(df, aes(x = ind, y = values, fill = category)) +
    geom_violin(scale = "width", trim = FALSE) +
    facet_wrap(~category, scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = "", y = "Senescence Score", title = "")

  p1
}
