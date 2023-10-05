plotVisiumRGB <- function(spe, vars, ...) {
    plt_df <- data.frame(colData(spe), spatialCoords(spe))

    if (any(!vars %in% names(plt_df))) {
        stop("One or more variables not found in the data.")
    }

    for (var in vars) {
        plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
    }

    num_vars <- length(vars)

    # Initialize RGB channels based on number of variables
    if (num_vars >= 3) {
        plt_df$R <- plt_df[[vars[1]]]
        plt_df$G <- plt_df[[vars[2]]]
        plt_df$B <- plt_df[[vars[3]]]
    } else {
        plt_df$R <- plt_df[[vars[1]]]
        plt_df$B <- plt_df[[vars[2]]]
        plt_df$G <- rep(0, nrow(plt_df))
    }

    if (num_vars >= 4) {
        W <- plt_df[[vars[4]]]
        plt_df$R <- plt_df$R + W * (1 - plt_df$R)
        plt_df$G <- plt_df$G + W * (1 - plt_df$G)
        plt_df$B <- plt_df$B + W * (1 - plt_df$B)
    }

    if (num_vars >= 5) {
        C <- plt_df[[vars[5]]]
        plt_df$G <- plt_df$G + C * (1 - plt_df$G)
        plt_df$B <- plt_df$B + C * (1 - plt_df$B)
    }

    if (num_vars == 6) {
        M <- plt_df[[vars[6]]]
        plt_df$R <- plt_df$R + M * (1 - plt_df$R)
        plt_df$B <- plt_df$B + M * (1 - plt_df$B)
    }

    spe$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
    plotVisium(spe, fill = "RGB", ...)+scale_fill_identity()#+
        #geom_point(aes(fill = RGB),size=size,stroke=stroke) #+
      #  scale_color_identity()
}

plotVisiumRGB
function(spe, vars, ...) {
    plt_df <- data.frame(colData(spe), spatialCoords(spe))

    if (any(!vars %in% names(plt_df))) {
        stop("One or more variables not found in the data.")
    }

    if (length(vars) > 4) {
        stop("A maximum of 4 variables is allowed.")
    }

    for (var in vars) {
        plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))
    }

    num_vars <- length(vars)

    # Initialize channels based on the number of variables:
    # Magenta (R and B channels)
    # Yellow (R and G channels)
    # Green (G channel)
    # Blue (B channel)

    if (num_vars >= 1) {
        plt_df$R <- plt_df[[vars[1]]] # Part of Magenta and Yellow
        plt_df$B <- plt_df[[vars[1]]] # Part of Magenta
    }

    if (num_vars >= 2) {
        plt_df$R <- plt_df$R + (1 - plt_df$R) * plt_df[[vars[2]]] # Yellow component
        plt_df$G <- plt_df[[vars[2]]] # Green
    }

    if (num_vars >= 3) {
        plt_df$G <- plt_df$G + (1 - plt_df$G) * plt_df[[vars[3]]] # Green component
    }

    if (num_vars == 4) {
        plt_df$B <- plt_df$B + (1 - plt_df$B) * plt_df[[vars[4]]] # Blue component
    }

    spe$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)
    plotVisium(spe, fill = "RGB", ...)+scale_fill_identity()
}

