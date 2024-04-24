### Notes
# The following code was adapted from the patternMarkers() function from `CoGAPS` to accept as input a matrix instead of a `CoGAPSObject`
# Additionally, matrix is not normalized before pattern marker calculation because RcppML already normalizes patterns
# Original function: https://rdrr.io/bioc/CoGAPS/man/patternMarkers-methods.html
patternMarkers <- function(featureLoadingsMatrix, sampleFactorsMatrix, threshold, axis,n)
{
    ## check inputs to the function
    if (!(threshold %in% c("cut", "all",'flex')))
        stop("threshold must be either 'cut' or 'all' or 'flex")
    if (!(axis %in% 1:2))
        stop("axis must be either 1 or 2")
    # Validate the new argument 'n'
    if (!is.numeric(n) || n < 1)
        stop("n must be a positive integer")

    ## need to scale each row of the matrix of interest so that the maximum is 1
    resultMatrix <- if (axis == 1) featureLoadingsMatrix else stop("Invalid axis for this function.")
    library(scales)

    #normedMatrix <- t(apply(resultMatrix, 1, function(row) row / max(row)))
    normedMatrix <- resultMatrix

    ## default pattern marker calculation, each pattern has unit weight
    markerScores <- sapply(1:ncol(normedMatrix), function(patternIndex)
        apply(normedMatrix, 1, function(row)
        {
            lp <- rep(0, ncol(normedMatrix))
            lp[patternIndex] <- 1
            return(sqrt(sum((row-lp)^2)))
        })
    )

    markerRanks <- apply(markerScores, 2, rank)
    colnames(markerScores) <- colnames(markerRanks) <- colnames(normedMatrix)

    ## Define the simplicityGENES function
    simplicityGENES <- function(As, Ps) {
        # rescale p's to have max 1
        pscale <- apply(Ps,1,max)

        # rescale A in accordance with p's having max 1
        As <- sweep(As, 2, pscale, FUN="*")

        # find the A with the highest magnitude
        Arowmax <- t(apply(As, 1, function(x) x/max(x)))

        # determine which genes are most associated with each pattern
        ssl <- matrix(NA, nrow=nrow(As), ncol=ncol(As), dimnames=dimnames(As))
        for (i in 1:ncol(As)) {
            lp <- rep(0, ncol(As))
            lp[i] <- 1
            ssl.stat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp)%*%(x-lp)))
            ssl[order(ssl.stat),i] <- 1:length(ssl.stat)
        }

        return(ssl)
    }

    ## keep only a subset of markers for each pattern depending on the type of threshold
    if (threshold == "cut")
    {
        simGenes <- simplicityGENES(As = resultMatrix, Ps = sampleFactorsMatrix)
        patternMarkers <- list()
        nP <- ncol(simGenes)

        for (i in 1:nP) {
            sortSim <- names(sort(simGenes[,i], decreasing = FALSE))
            geneThresh <- min(which(simGenes[sortSim,i] > apply(simGenes[sortSim,], 1, min)))
            markerGenes <- sortSim[1:geneThresh]
            markerGenes <- unique(markerGenes)
            patternMarkers[[i]] <- markerGenes
        }

        markersByPattern <- patternMarkers
    }
    else if (threshold == "all") # only the markers with the lowest scores
    {
        min_indices <- apply(markerScores, 1, which.min)
        patternsByMarker <- colnames(markerScores)[sapply(min_indices, `[`, 1)]
        markersByPattern <- sapply(colnames(markerScores), USE.NAMES = TRUE, simplify = FALSE,
                                   function(pattern) rownames(markerScores)[which(patternsByMarker == pattern)])
     }

    else if (threshold == "flex") {
        # Assign genes to a pattern if they are among the three lowest values for that gene
        flexPatternsByGene <- apply(markerScores, 1, function(geneScores) {
            lowestThreePatterns <- order(geneScores)[1:2]
            return(colnames(markerScores)[lowestThreePatterns])
        })

        markersByPattern <- lapply(colnames(markerScores), function(pattern) {
            genesAssignedToPattern <- rownames(markerScores)[which(flexPatternsByGene == pattern, arr.ind = TRUE)]
            return(unique(genesAssignedToPattern))
        })
    }

    ## add TopRankedGenes
    topRankedGenes <- list()
    for (patternIndex in 1:ncol(markerScores)) {
        geneScores <- markerScores[, patternIndex]
        geneRanks <- markerRanks[, patternIndex]

        # Filter out genes with zero loading for this pattern
        nonZeroIndices <- which(featureLoadingsMatrix[, patternIndex] != 0)
        filteredGeneScores <- geneScores[nonZeroIndices]
        filteredGeneRanks <- geneRanks[nonZeroIndices]

        # Sort the filtered gene scores and ranks and take the top N
        sortedGeneIndices <- order(filteredGeneScores)
        topNIndices <- sortedGeneIndices[1:min(n, length(sortedGeneIndices))]

        # Extract the gene names for these top N indices
        topRankedGenes[[patternIndex]] <- rownames(markerScores)[nonZeroIndices[topNIndices]]
    }

    return(list(
        "PatternMarkers" = markersByPattern,
        "PatternMarkerRanks" = markerRanks,
        "PatternMarkerScores" = markerScores,
        "TopRankedGenes" = topRankedGenes
    ))
}

