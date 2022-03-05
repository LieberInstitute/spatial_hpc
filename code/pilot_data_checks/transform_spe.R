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
trans_geom = function(x, sample_id = NULL, refl = FALSE, trans = c(0, 0), degrees = 0) {
  tol = 1e-8
  
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
    if (! (sample_id %in% colData(x)$sample_id)) {
      stop(
        paste0(
          "'", sample_id, "' is not a sample in",
          "'colData(x)$sample_id': cannot subset."
        )
      )
    }
    
    x = x[, colData(x)$sample_id == sample_id]
  }
  
  #   For consistency with 'rotateImg', we use degrees and consider the
  #   clockwise direction positive
  radians = -1 * degrees * pi / 180
  
  if (refl) {
    refl_vec = c(-1, 1)
  } else {
    refl_vec = c(1, 1)
  }
  
  #   Determine the matrix by which left-multiplication represents rotation
  rotation_mat = matrix(
    c(
      cos(radians), sin(radians), -1 * sin(radians), cos(radians)
    ),
    nrow = 2
  )
  
  #   When performing the rotation, we want the coordinates to be centered at
  #   the origin
  x_center = (max(spatialCoords(x)[,'pxl_col_in_fullres']) + min(spatialCoords(x)[,'pxl_col_in_fullres'])) / 2
  y_center = (max(spatialCoords(x)[,'pxl_row_in_fullres']) + min(spatialCoords(x)[,'pxl_row_in_fullres'])) / 2
  dist_to_origin = c(x_center, y_center)
  
  # Perform geometric transformation(s)
  new_coords = t(spatialCoords(x)) - dist_to_origin   # center at origin
  new_coords = refl_vec * new_coords                  # reflect across v. axis
  new_coords = rotation_mat %*% new_coords            # rotate about origin
  new_coords = t(new_coords + dist_to_origin + trans) # re-center + translate
  
  colnames(new_coords) = colnames(spatialCoords(x))
  
  # Verify the center of the points has moved exactly by 'trans'
  new_x_center = (max(new_coords[,1]) + min(new_coords[,1])) / 2
  new_y_center = (max(new_coords[,2]) + min(new_coords[,2])) / 2
  new_dist_to_origin = c(new_x_center, new_y_center)
  stopifnot(all(new_dist_to_origin - dist_to_origin - trans < tol))
  
  # Ensure points are at integer values
  new_coords = round(new_coords)
  
  #   Return a copy of the SpatialExperiment with the new coordinates
  spatialCoords(x) = new_coords
  return(x)
}

#   Transform both the 'spatialCoords' and 'imgData' of a SpatialExperiment
#   object, returning the transformed copy of the full object.
transform_spe = function(x, refl = FALSE, trans = c(0, 0), degrees = 0) {
  x = trans_geom(x, refl = FALSE, trans = c(0, 0), degrees = 0)
  
  if (refl) {
    x = mirrorImg(x, axis = "v")
  }
  
  if (degrees != 0) {
    x = rotateImg(x, degrees=degrees)
  }
  
  return(x)
}

xlab = "spatialCoords axis 2"
ylab = "spatialCoords axis 1"

#   Given a SpatialExperiment object 'x', plot the 'spatialCoords' in the same
#   orientation as 'plot(imgRaster(x))'. Note the names of the 'spatialCoords'
#   columns were swapped (a bug) until a recent update to SpatialExperiment
#   under Bioc 3.15, so this function works properly since Bioc 3.15.
plot_spatial_coords = function(x, title) {
  ggplot(
    data.frame(spatialCoords(x)),
    aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)
  ) + 
    geom_point() +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = title)
}

