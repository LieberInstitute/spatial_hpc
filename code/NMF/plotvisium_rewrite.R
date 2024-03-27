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
#' @param fill (character) Column in \code{colData} to use to fill points by
#'   color. If \code{fill} contains a numeric column (e.g. total UMI counts), a
#'   continuous color scale will be used. If \code{fill} contains a factor (e.g.
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
#' @param assay (character) Name of assay data to use when \code{fill} is in
#'   \code{rownames(spe)}. Should be one of \code{assayNames(spe)}.
#'
#' @param trans Transformation to apply for continuous scales. Ignored unless
#'   \code{fill} is numeric, e.g. feature expression. (See
#'   \code{\link{ggplot2}{continuous_scale}} for valid options.)
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
#' @param palette (character) Color palette for points. Options for discrete
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
#'   scale_fill_viridis_c scale_color_identity scale_fill_identity facet_wrap
#'   guides guide_colorbar guide_legend theme_void element_text margin unit
#'   layer
#' @importFrom grid rasterGrob
#' @importFrom methods is as
#' @importFrom stats setNames
#'
#' @export
#'
#' @author Helena L. Crowell with modifications by Lukas M. Weber
#'
#' @examples
#' library(STexampleData)
#'
#' spe <- Visium_mouseCoronal()
#'
#' # color by x coordinate, highlight in-tissue spots
#' plotVisium(spe, fill = "pxl_col_in_fullres", highlight = "in_tissue")
#'
#' # subset in-tissue spots
#' sub <- spe[, as.logical(colData(spe)$in_tissue)]
#'
#' # color by feature counts, don't include image
#' rownames(sub) <- make.names(rowData(sub)$gene_name)
#' plotVisium(sub, fill = "Gad2", assay = "counts")
#'
plotVisium <- function(spe,
                       spots = TRUE, fill = NULL, highlight = NULL,
                       facets = "sample_id", image = TRUE,
                       assay = "logcounts", trans = "identity",
                       x_coord = NULL, y_coord = NULL, y_reverse = TRUE,
                       sample_ids = NULL, image_ids = NULL, palette = NULL,
                       values=NULL) {

    # check validity of input arguments
    stopifnot(
        is(spe, "SpatialExperiment"),
        is.logical(spots), length(spots) == 1,
        is.logical(image), length(image) == 1,
        is.logical(y_reverse), length(y_reverse) == 1)

    if (is.null(x_coord)) x_coord <- spatialCoordsNames(spe)[1]
    if (is.null(y_coord)) y_coord <- spatialCoordsNames(spe)[2]

    # set up data for plotting
    plt_df <- data.frame(colData(spe), spatialCoords(spe))
    if (!is.null(fill)) {
        # check validity of 'fill' argument
        stopifnot(is.character(fill), length(fill) == 1)
        if (!fill %in% c(names(plt_df), rownames(spe))) {
            stop("'fill' should be in rownames(spe) or names(colData(spe))")
        }
        # (optionally) add feature assay data to 'plt_df'
        if (fill %in% rownames(spe)) {
            stopifnot(
                is.character(assay),
                length(grep(assay, assayNames(spe))) == 1)
            plt_df[[fill]] <- assay(spe, assay)[fill, ]
        }
        # get color palette
        palette <- .get_pal(palette, plt_df[[fill]])
    } else {
        fill <- "foo"
        plt_df[[fill]] <- "black"
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
    } else {
        img <- NULL
        images <- xlim <- ylim <- NULL
    }

    # scale spatial coordinates
    for (s in sample_ids) {
        ix <- plt_df$sample_id == s
        xy <- c(x_coord, y_coord)
        sf <- img_df[s, "scaleFactor"]
        plt_df[ix, xy] <- sf * plt_df[ix, xy]
        # reverse y coordinates to match orientation of images
        # (sometimes required for Visium data)
        if (y_reverse) plt_df <- .y_reverse(plt_df, ix, y_coord, img)
    }

    # construct points and highlights
    if (spots) {
        # check whether 'fill' is continuous (numeric) or discrete (factor)
        guide <- ifelse(is.numeric(plt_df[[fill]]), guide_colorbar, guide_legend)
        points <- list(
            guides(fill = guide(
                title = NULL, order = 1, override.aes = list(col = NA, size = 3))),
            geom_point(shape = 21, size = 1, stroke = 0.25, alpha = 1))
        if (!is.null(highlight)) {
            plt_df$highlight <- as.factor(plt_df[[highlight]])
            highlights <- list(
                scale_color_manual(highlight, values = values),
                guides(col = guide_legend(override.aes = list(
                    size = 2, stroke = 1,
                    col = values[seq_along(unique(plt_df$highlight))]))))
        } else {
            plt_df$highlight <- "transparent"
            highlights <- scale_color_identity()
        }
    } else {
        # this is required, else the image layer doesn't show
        points <- geom_point(col = "transparent")
        highlights <- NULL
    }

    # color scale
    scale <- if (fill != "foo") {
        if (is.numeric(plt_df[[fill]])) {
            if (length(palette) == 1 && palette == "viridis") {
                scale_fill_viridis_c(trans = trans)
            } else {
                scale_fill_gradient(low = palette[1], high = palette[2], trans = trans)
            }
        } else {
            scale_fill_manual(values = palette)
        }
    } else scale_fill_identity()

    # display plot
    ggplot(plt_df,
           aes_string(x_coord, y_coord, fill = fill, col = "highlight")) +
        images + points + highlights + scale +
        coord_fixed(xlim = xlim, ylim = ylim) +
        theme_void() +
        theme(legend.key.size = unit(0.5, "lines"),
              strip.text = element_text(margin = margin(0, 0, 0.5, 0, "lines"))) +
        if (!is.null(facets)) facet_wrap(facets,ncol = 2)
}

###################
# utility functions
###################

.get_pal <- function(pal, val) {
    # (if missing) default to 'viridis' for continuous
    # and 'libd_layer_colors' for discrete scales
    if (is.null(pal)) {
        pal <- ifelse(is.numeric(val), "viridis", "libd_layer_colors")
    }
    # if length(palette) > 1, use palette as provided
    # (either multiple colors for discrete labels,
    # or length 2 for continuous gradient)
    if (length(pal) == 1) {
        pal <- switch(pal,
                      "libd_layer_colors" = c(
                          "#F0027F", "#377EB8", "#4DAF4A", "#984EA3",
                          "#FFD700", "#FF7F00", "#1A1A1A", "#666666"),
                      "Okabe-Ito" = c(
                          "#000000", "#E69F00", "#56B4E9", "#009E73",
                          "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                      # use 'scale_fill_viridis_c' instead of hex codes or color names
                      "viridis" = pal,
                      # if providing a single color name (e.g. 'navy' or 'red'),
                      # combine with 'gray95' for continuous color scale
                      c("gray95", pal)
        )
    }
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
