# ####Function to rescale multiple continuous variables in one RGB space and plot
# ###Args:
# ##spe: a SpatialExperiment object
# ##vars: names of variables  to plot. Can be genes with names in rownames(spe) OR colData columns
# ##independent: a logical asking whether to rescale variables independently or not.
# #If FALSE, rescale all vars to 0 to 1 from min to max of all vars
# #if TRUE, rescale each var to 0 to 1 from min to max of each var independently
# #default TRUE
# ##assay: which assay in assays(spe) to plot
# #default is logcounts
# plotVisiumRGB <- function(spe, vars, independent=TRUE, assay='logcounts', ...) {
#     plt_df <- data.frame(colData(spe), spatialCoords(spe))
#     if (any(vars %in% rownames(spe))) {
#         df<-assay(spe,assay)[rownames(assay(spe,assay)) %in% vars,]
#         plt_df<-cbind(t(as.matrix(df)),plt_df)
#         rm(df)
#     }
#
#     if (any(!vars %in% names(plt_df))) {
#         stop("One or more variables not found in the data.")
#     }
#
#     if(independent==TRUE){
#     for (var in vars) {
#         plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
#     }
#     }
#     if(independent==FALSE){
#         # Find the minimum and maximum values across all variables in vars
#         min_value <- min(unlist(lapply(plt_df[vars], min)))
#         max_value <- max(unlist(lapply(plt_df[vars], max)))
#
#         # Scale each variable in vars dependently from 0 to 1
#         for (var in vars) {
#             plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1), from = c(min_value, max_value))
#         }
#     }
#     num_vars <- length(vars)
#
#     # Initialize RGB channels based on number of variables
#     if (num_vars >= 3) {
#         plt_df$R <- plt_df[[vars[1]]]
#         plt_df$G <- plt_df[[vars[2]]]
#         plt_df$B <- plt_df[[vars[3]]]
#     } else {
#         plt_df$R <- plt_df[[vars[1]]]
#         plt_df$B <- plt_df[[vars[2]]]
#         plt_df$G <- rep(0, nrow(plt_df))
#     }
#
#     if (num_vars >= 4) {
#         W <- plt_df[[vars[4]]]
#         plt_df$R <- plt_df$R + W * (1 - plt_df$R)
#         plt_df$G <- plt_df$G + W * (1 - plt_df$G)
#         plt_df$B <- plt_df$B + W * (1 - plt_df$B)
#     }
#
#     if (num_vars >= 5) {
#         C <- plt_df[[vars[5]]]
#         plt_df$G <- plt_df$G + C * (1 - plt_df$G)
#         plt_df$B <- plt_df$B + C * (1 - plt_df$B)
#     }
#
#     if (num_vars == 6) {
#         M <- plt_df[[vars[6]]]
#         plt_df$R <- plt_df$R + M * (1 - plt_df$R)
#         plt_df$B <- plt_df$B + M * (1 - plt_df$B)
#     }
#
#     spe$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
#     plotVisium(spe, fill = "RGB", ...)+scale_fill_identity()#+
#         #geom_point(aes(fill = RGB),size=size,stroke=stroke) #+
#       #  scale_color_identity()
# }

plotVisiumRGB<-function(spe, pink = NULL, yellow = NULL, blue = NULL, green = NULL, red = NULL, cyan = NULL, ...) {
    plt_df <- data.frame(colData(spe), spatialCoords(spe))

    plt_df <- data.frame(colData(spe), spatialCoords(spe))

    # Check for the presence of variables in colData and rownames and rescale them
    variables <- list(pink = pink, yellow = yellow, blue = blue, green = green, red = red, cyan = cyan)
    for (var_name in names(variables)) {
        var <- variables[[var_name]]
        if (!is.null(var)) {
            if (var %in% names(plt_df)) {
                # Variable found in colData
                plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
            } else if (var %in% rownames(spe)) {
                # Variable found in rownames, add corresponding logcounts row to plt_df
                logcounts_row <- logcounts(spe)[var, ]
                plt_df[[var]] <- scales::rescale(logcounts_row, to = c(0, 1))
            } else {
                stop(paste("Variable", var, "not found in the data."))
            }
        }
    }

    # Initialize RGB channels
    plt_df$R <- plt_df$G <- plt_df$B <- rep(0, nrow(plt_df))

    # Assign channels based on selected colors
    if (!is.null(pink)) {
        plt_df$R <- plt_df[[pink]]
        plt_df$B <- plt_df[[pink]]
    }

    if (!is.null(yellow)) {
        plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[yellow]]
        plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[yellow]]
    }

    if (!is.null(blue)) {
        plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[blue]]
    }

    if (!is.null(green)) {
        plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[green]]
    }

    if (!is.null(red)) {
        plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[red]]
    }

    if (!is.null(cyan)) {
        plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[cyan]]
        plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[cyan]]
    }

    # Create RGB color
    spe$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
    plotVisium(spe, fill = "RGB", ...)+scale_fill_identity()
}
