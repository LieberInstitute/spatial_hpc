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
    plotVisium(spe,fill='RGB',...)+scale_fill_identity()
}
