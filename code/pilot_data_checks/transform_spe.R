#   Return a copy of the input SpatialExperiment where a set of geometric
#   transformations has been applied. These consist of an optional reflection
#   (across the vertical axis prior to any rotation), an optional translation,
#   and an optional rotation (about the center of the coordinates).
#
#   x: a SpatialExperiment object
#
#   sample_id: either a length-1 character vector, to subset 'x' to the
#           specified sample ID (colData(x)$sample_id must be a non-NULL
#           column!), or NULL to specify no subsetting.
#
#   refl: a length-1 logical vector indicating whether to reflect coordinates
#           about the vertical axis prior to any rotation (affecting
#           spatialCoords(x)[,1])
#
#   trans: a length-2 integer vector corresponding to values to add to the
#           existing coordinates (applying a 2D translation).
#
#   degrees: a length-1 float vector, corresponding to the clockwise
#           rotation, in degrees, to apply about the center of the coordinates.
#           This is currently expected to be a multiple of 90 degrees.
trans_geom <- function(x, sample_id = NULL, refl = FALSE, trans = c(0, 0), degrees = 0) {
    x_name <- "pxl_col_in_fullres"
    y_name <- "pxl_row_in_fullres"

    stopifnot(all(spatialCoordsNames(x) == c(x_name, y_name)))

    # 'degrees' should be a single numeric divisible by 90
    stopifnot(
        length(degrees) == 1,
        is.numeric(degrees),
        degrees %% 90 == 0
    )

    #   Subset to a particular sample if requested
    if (is.null(sample_id)) {
        warning(
            paste0(
                "Applying transformation to entire object. For multi-sample",
                " objects, this might not be desired."
            )
        )
    } else {
        #   Ensure 'sample_id' is a column name that exists
        if (is.null(colData(x)$sample_id)) {
            stop(
                paste0(
                    "'sample_id' must be a column name in the",
                    " SpatialExperiment in order to subset by sample."
                )
            )
        }

        #   Ensure the sample is present in the object
        if (!(sample_id %in% colData(x)$sample_id)) {
            stop(
                paste0(
                    "'", sample_id, "' is not a sample in",
                    "'colData(x)$sample_id': cannot subset."
                )
            )
        }

        x <- x[, colData(x)$sample_id == sample_id]
    }

    #   For consistency with 'rotateImg', we use degrees. Note that positive
    #   angles represent counter-clockwise rotations
    radians <- degrees * pi / 180

    if (refl) {
        refl_vec <- c(-1, 1)
    } else {
        refl_vec <- c(1, 1)
    }

    #   Determine the matrix by which left-multiplication represents rotation
    rotation_mat <- matrix(
        c(
            cos(radians), sin(radians), -1 * sin(radians), cos(radians)
        ),
        nrow = 2
    )

    #   Get the dimensions of the "rectangle" containing the set of
    #   spatialCoords within the object
    dim_max <- dim(imgRaster(x)) / scaleFactors(x)[1]

    new_coords <- refl_vec * t(spatialCoords(x)) # reflect across v. axis

    #   If reflecting across v. axis, return "rectangle" to its original
    #   location
    if (refl) {
        new_coords <- new_coords + c(dim_max[2], 0)
    }

    new_coords <- rotation_mat %*% new_coords # rotate about origin

    #   Since the rotation is about the origin, we'll need to return the
    #   "rectangle" such that its top left corner is at the spot where its
    #   previous top left corner was
    if (degrees %% 360 == 90) {
        new_coords <- new_coords + c(dim_max[1], 0)
    } else if (degrees %% 360 == 180) {
        new_coords <- new_coords + rev(dim_max)
    } else if (degrees %% 360 == 270) {
        new_coords <- new_coords + c(0, dim_max[2])
    }

    new_coords <- t(new_coords + trans) # transpose and translate
    #   Add names to spatialCoords of the new object
    colnames(new_coords) <- colnames(spatialCoords(x))

    #   Ensure points are at integer values
    new_coords <- round(new_coords)

    #   Return a copy of the SpatialExperiment with the new coordinates
    spatialCoords(x) <- new_coords
    return(x)
}


#   Transform both the 'spatialCoords' and 'imgData' of a SpatialExperiment
#   object, returning the transformed copy of the full object.
transform_spe <- function(x, refl = FALSE, trans = c(0, 0), degrees = 0) {
    x <- trans_geom(x, refl = FALSE, trans = c(0, 0), degrees = 0)

    if (refl) {
        x <- mirrorImg(x, axis = "v")
    }

    if (degrees != 0) {
        x <- rotateImg(x, degrees = degrees)
    }

    return(x)
}

#   Given a SpatialExperiment object 'x', plot the 'spatialCoords' in the same
#   orientation as 'plot(imgRaster(x))'. Note the names of the 'spatialCoords'
#   columns were swapped (a bug) until a recent update to SpatialExperiment
#   under Bioc 3.15, so this function works properly since Bioc 3.15.
plot_spatial_coords <- function(x, title) {
    ggplot(
        data.frame(spatialCoords(x)),
        aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)
    ) +
        geom_point() +
        coord_fixed() +
        scale_y_reverse() +
        labs(title = title)
}
