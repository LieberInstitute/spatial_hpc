#' plotVisium
#' 
#' Plots for spatially resolved transcriptomics data from the 10x Genomics
#' Visium platform
#' 
#' Function to generate plots for spatially resolved transcriptomics datasets
#' from the 10x Genomics Visium spatially platform.
#' 
#' This function generates a plot for spot-based spatially resolved
#' transcriptomics data from the 10x Genomics Visium platform, with several
#' options available to adjust the plot type and style.
#' 
#' 
#' @param spe (SpatialExperiment) Input data object.
#' 
#' @param spots (logical) Whether to display spots (spatial barcodes) as points.
#'   Default = TRUE.
#' 
#' @param annotate (character) Column in \code{colData} to use to fill points by
#'   color. If \code{annotate} contains a numeric column (e.g. total UMI counts), a
#'   continuous color scale will be used. If \code{annotate} contains a factor (e.g.
#'   cluster labels), a discrete color scale will be used. Default = NULL.
#' 
#' @param highlight (character) Column in \code{colData} to use to highlight
#'   points by outlining them. For example, \code{in_tissue} will highlight
#'   spots overlapping with tissue. Default = NULL.
#' 
#' @param facets (character) Column in \code{colData} to use to facet plots,
#'   i.e. show multiple panels of plots. Default = "sample_id". Set to NULL to
#'   disable.
#' 
#' @param image (logical) Whether to show histology image as background. Default
#'   = TRUE.
#'   
#' @param zoom (logical) Whether to zoom to area of tissue containing spots.
#'   Default = FALSE
#'   
#' @param show_axes (logical) Whether to show axes and coordinates. Default =
#'   FALSE
#' 
#' @param assay (character) Name of assay data to use when \code{annotate} is in
#'   \code{rownames(spe)}. Should be one of \code{assayNames(spe)}.
#'   
#' @param trans Transformation to apply for continuous scales. Ignored unless
#'   \code{annotate} is numeric, e.g. feature expression. (See
#'   \code{\link{ggplot2}{continuous_scale}} for valid options.)
#' 
#' @param point_size (numeric) Point size. Default = 1.
#' 
#' @param legend_position Legend position for annotations. Options are "left",
#'   "right", "top", "bottom", and "none". Default = "right".
#' 
#' @param x_coord (character) Column in \code{spatialCoords} containing
#'   x-coordinates. Default = NULL, which selects the first column.
#' 
#' @param y_coord (character) Column in \code{spatialCoords} containing
#'   y-coordinates. Default = NULL, which selects the second column.
#' 
#' @param y_reverse (logical) Whether to reverse y coordinates, which is often 
#'   required for Visium data, depending on the orientation of the raw data.
#'   Default = TRUE.
#' 
#' @param pal (character) Color palette for points. Options for discrete
#'   labels are "libd_layer_colors", "Okabe-Ito", or a custom vector of hex
#'   color codes. Options for continuous values are "viridis", a single color
#'   name (e.g. "red", "navy", etc), or a vector of length two containing color
#'   names for each end of the scale. Default = "libd_layer_colors" for discrete
#'   data, and "viridis" for continuous data.
#' 
#' @param sample_ids (character) Samples to show, if multiple samples are
#'   available. Default = NULL (show all samples).
#' 
#' @param image_ids (character) Images to show, if multiple images are
#'   available. Default = NULL (show all images).
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords spatialCoordsNames imgData
#'   'imgData<-' imgRaster scaleFactors
#' @importFrom SummarizedExperiment colData assayNames
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual scale_fill_gradient
#'   scale_fill_gradientn scale_fill_viridis_c scale_color_identity
#'   scale_fill_identity facet_wrap guides guide_colorbar guide_legend
#'   theme_void element_text margin unit layer
#' @importFrom grid rasterGrob
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales hue_pal
#' @importFrom methods is as
#' @importFrom stats setNames
#' 
#' @export
#' 
#' @author Helena L. Crowell, with modifications by Lukas M. Weber and Yixing E.
#'   Dong
#' 
#' @examples
#' library(STexampleData)
#' 
#' spe <- Visium_mouseCoronal()
#' 
#' # color by x coordinate, highlight in-tissue spots
#' plotVisium(spe, annotate = "pxl_col_in_fullres", highlight = "in_tissue")
#' 
#' # subset in-tissue spots
#' sub <- spe[, as.logical(colData(spe)$in_tissue)]
#' 
#' # color by feature counts, don't include image
#' rownames(sub) <- make.names(rowData(sub)$gene_name)
#' plotVisium(sub, annotate = "Gad2", assay = "counts")
#' 
plotVisium <- function(spe, 
                       spots = TRUE, annotate = NULL, highlight = NULL, 
                       facets = "sample_id", image = TRUE, zoom = FALSE, show_axes = FALSE,
                       assay = "counts", trans = "identity", point_size = 1, legend_position = "right",
                       x_coord = NULL, y_coord = NULL, y_reverse = TRUE, 
                       sample_ids = NULL, image_ids = NULL, pal = NULL, stroke_size=0.25) {
  
  # check validity of input arguments
  stopifnot(
    is(spe, "SpatialExperiment"), 
    is.logical(spots), length(spots) == 1, 
    is.logical(image), length(image) == 1, 
    is.logical(y_reverse), length(y_reverse) == 1)
  
  stopifnot(legend_position %in% c("left", "right", "top", "bottom", "none"))
  
  if (!is.null(annotate)) {
    stopifnot(is.character(annotate))
  }
  
  if(is.null(x_coord)) x_coord <- spatialCoordsNames(spe)[1]
  if(is.null(y_coord)) y_coord <- spatialCoordsNames(spe)[2]
  
  # set up data for plotting
  df <- data.frame(colData(spe), spatialCoords(spe))
  if (!is.null(annotate)) {
    # check validity of 'annotate' argument
    stopifnot(is.character(annotate), length(annotate) == 1)
    if (!annotate %in% c(names(df), rownames(spe))) {
      stop("'annotate' should be in rownames(spe) or names(colData(spe))")
    }
    # (optionally) add feature assay data to 'df'
    if (annotate %in% rownames(spe)) {
      stopifnot(
        is.character(assay), 
        length(grep(assay, assayNames(spe))) == 1)
      df[[annotate]] <- assay(spe, assay)[annotate, ]
    }
    if (is.numeric(df[[annotate]]) & is.null(pal)) {
      # for continuous feature, ensure length(pal) == 1 (instead of 0 if NULL)
      pal <- "seuratlike"
    }
    # get color palette
    pal <- .get_pal(pal, df[[annotate]])
  } else {
    annotate <- "foo"
    df[[annotate]] <- "black"
  }
  
  if (is.null(sample_ids)) {
    # default to using all samples
    sample_ids <- unique(spe$sample_id)
  } else {
    # subset specified samples
    spe <- spe[, spe$sample_id %in% sample_ids]
  }
  
  # subset selected images
  img_df <- .sub_imgData(spe, sample_ids, image_ids)
  rownames(img_df) <- img_df$sample_id
  
  # construct image layers
  # note: images could also be plotted using 'annotation_custom()', 
  # however, this does not allow for faceting, so we instead 
  # construct a separate image layer for each sample
  if (image) {
    images <- lapply(sample_ids, function(s) {
      spi <- img_df[s, "data"]
      img <- imgRaster(spi[[1]])
      layer(
        data = data.frame(sample_id = s), 
        inherit.aes = FALSE, 
        stat = "identity", 
        position = "identity", 
        geom = ggplot2::GeomCustomAnn, 
        params = list(
          grob = rasterGrob(img), 
          xmin = 0, xmax = ncol(img), 
          ymin = 0, ymax = nrow(img))
      )
    })
    img <- img_df$data[[1]]
    xlim <- c(0, ncol(img))
    ylim <- c(0, nrow(img))
    
    if (zoom) {
      xlim <- ylim <- NULL
    }
  } else {
    img <- NULL
    images <- xlim <- ylim <- NULL
  }
  
  # scale spatial coordinates
  for (s in sample_ids) {
    ix <- df$sample_id == s
    xy <- c(x_coord, y_coord)
    sf <- img_df[s, "scaleFactor"]
    df[ix, xy] <- sf * df[ix, xy]
    # reverse y coordinates to match orientation of images 
    # (sometimes required for Visium data)
    if (y_reverse) df <- .y_reverse(df, ix, y_coord, img)
  }
  
  # construct points and highlights
  if (spots) {
    # check whether 'annotate' is continuous (numeric) or discrete (factor)
    guide <- ifelse(is.numeric(df[[annotate]]), guide_colorbar, guide_legend)
    points <- list(
      guides(fill = guide(
        title = annotate, order = 1, override.aes = list(col = NA, size = 3))), 
      geom_point(shape = 21, size = point_size, stroke = stroke_size, alpha = 0.8))
    if (!is.null(highlight)) {
      df$highlight <- as.factor(df[[highlight]])
      highlights <- list(
        scale_color_manual(highlight, values = c("#e0e0e0", "black")), 
        guides(col = guide_legend(override.aes = list(
          size = 2, stroke = 1, 
          col = c("#e0e0e0", "black")[seq_along(unique(df$highlight))]))))
    } else {
      df$highlight <- "transparent"
      highlights <- scale_color_identity()
    }
  } else {
    # this is required, else the image layer doesn't show
    points <- geom_point(col = "transparent")
    highlights <- NULL
  }
  
  # color scale
  scale <- if(annotate != "foo") {
    if (is.numeric(df[[annotate]])) {
      if (length(pal) == 1 && 
          pal %in% c("viridis", "magma", "inferno", "plasma", 
                     "cividis", "rocket", "mako", "turbo")) {
        scale_fill_viridis_c(trans = trans, option = pal)
      } else if (length(pal) == 1 && pal == "seuratlike") {
        scale_fill_gradientn(
          colors = colorRampPalette(
            colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100), 
          trans = trans, 
          limits = c(min(df[[annotate]]), max(df[[annotate]])))
      } else {
        scale_fill_gradient(low = pal[1], high = pal[2], trans = trans)
      }
    } else if (is.factor(df[[annotate]])) {
      # for categorical feature, automate palette
      if (is.null(pal)) {
        scale_fill_manual(
          name = annotate, 
          values = hue_pal()(length(unique(df[[annotate]]))))
      } else if (!is.null(pal)) {
        scale_fill_manual(values = pal)
      }
    }
  } else {
    scale_fill_identity()
  }
  
  # display plot
  p <- ggplot(df, 
              aes_string(x_coord, y_coord, fill = annotate, col = "highlight")) + 
    images + points + highlights + scale + 
    coord_fixed(xlim = xlim, ylim = ylim) 
  
  if (show_axes) {
    p <- p + 
      theme_bw() + 
      theme(strip.text = element_text(margin = margin(0, 0, 0.5, 0, "lines"), 
                                      size = 12), 
            legend.position = legend_position) +
      labs(x = paste0("pxl_col_in_", img_df[s, "image_id"]),
           y = paste0("pxl_col_in_", img_df[s, "image_id"])) + 
      if (!is.null(facets)) facet_wrap(facets)
  } else {
    p <- p + 
      theme_void() + 
      theme(strip.text = element_text(margin = margin(0, 0, 0.5, 0, "lines"), 
                                      size = 12), 
            legend.position = legend_position) + 
      if (!is.null(facets)) facet_wrap(facets)
  }
  
  p
}

.get_pal <- function(pal, val) {
  
  if (length(pal) == 1) {
    pal <- switch(pal, 
                  "libd_layer_colors" = c(
                    "#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
                    "#FFD700", "#FF7F00", "#1A1A1A", "#666666"), 
                  "Okabe-Ito" = c(
                    "#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 
                  # use 'scale_fill_viridis_c' for the following options
                  "viridis" = pal, 
                  "magma" = pal, 
                  "inferno" = pal, 
                  "plasma" = pal, 
                  "viridis" = pal, 
                  "cividis" = pal, 
                  "rocket" = pal, 
                  "mako" = pal, 
                  "turbo" = pal, 
                  "seuratlike" = pal, 
                  # for a single color name, combine with "gray95" for continuous color scale
                  c("gray95", pal)
    )
  }
  
  # if length(pal) == 0 (i.e. 'pal' is NULL), leave 'pal' unchanged and the
  # plotting functions will select default palettes instead
  
  # if length(pal) > 1, use 'pal' as provided (e.g. multiple colors for discrete
  # labels, or length 2 for continuous gradient)
  
  return(pal)
}


#' @importFrom SpatialExperiment imgData
.sub_imgData <- function(spe, sample_ids, image_ids) {
  .get_img_idx <- SpatialExperiment:::.get_img_idx
  if (is.null(image_ids)) {
    # default to first available image for each sample
    idx <- .get_img_idx(spe, TRUE, NULL)
  } else {
    if (length(image_ids) == 1) {
      idx <- .get_img_idx(spe, TRUE, image_ids)
    } else {
      stopifnot(length(image_ids) == length(sample_ids))
      idx <- mapply(s = sample_ids, i = image_ids,
                    function(s, i) .get_img_idx(spe, s, i))
    }
  }
  imgData(spe)[idx, ]
}


.y_reverse <- function(df, ix, y, img) {
  y_tmp <- df[ix, y]
  if (!is.null(img)) {
    y_tmp <- nrow(img) - y_tmp
  } else {
    y_tmp <- max(y_tmp) - y_tmp
  }
  df[ix, y] <- y_tmp
  return(df)
}
