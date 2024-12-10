# as of 2024/12/10 this CoGAPS patternMarkers rewrite to run with matrices is nearly identical to source code found here:
# https://github.com/FertigLab/CoGAPS/blob/5595e8c69120ff67991385398b5c8d0fd34780d1/R/methods-CogapsResult.R#L393
### lines 393 to 494
### I added the rescale argument which allows the user to skip the CoGAPS rescaling of the Wmatrix based on Hmatrix max
### I also added a dim check to see if the Hmatrix needs to be transposed
### tecnhically if 'rescale=FALSE' the Hmatrix isn't needed but I kept the function as similar to the source code as possible

patternMarkers_m <- function(Wmatrix, Hmatrix, threshold, rescale=FALSE, lp=NULL, axis){
            #look for features-markers of patterns with axis=1
            #or samples-markers of patterns with axis=2
            if(axis == 1){
              Amatrix <- Wmatrix
	      if(dim(Hmatrix)[1]==dim(Wmatrix)[2]) {Pmatrix = Hmatrix}
	      else {Pmatrix <- t(Hmatrix)}
            } else if(axis == 2){
              Amatrix <- Hmatrix
              Pmatrix <- t(Wmatrix)
            } else {
              stop("axis must be 1 or 2")
            }
            
            if(rescale==TRUE) {
            # determine norm for A if Ps were rescaled to have max 1
            pscale <- apply(Pmatrix,1,max)
            
            # rescale A in accordance with p's having max 1
            Amatrix <- sweep(Amatrix, 2, pscale, FUN="*")
            }
            # normalize each row of A to have max 1
            Arowmax <- t(apply(Amatrix, 1, function(x) x/max(x)))
            
            nP=dim(Amatrix)[2]
            
            if(!is.null(lp)){
              if(!(unique(lengths(lp)) > 1 && unique(lengths(lp)) == nP)){
                warning("lp length must equal the number of columns of the Amatrix")
              }
              # container for feature ranks by L2 distance from lp in case of provided lp
              ssranks<-matrix(NA, nrow=nrow(Amatrix),ncol=length(lp),
                              dimnames=list(rownames(Amatrix), names(lp)))
              if(max(sapply(lp,max))>1){
                stop("lp should be a list of vectors with max value of 1")
              }
            } else {
              lp <- list()
              for(i in 1:nP){
                lp_mock <- rep(0,nP)
                lp_mock[i] <- 1
                lp[[i]] <- lp_mock
              }
              names(lp) <- colnames(Amatrix)
              # container for feature ranks by L2 distance from lp in case of one-hot encoded lp
              ssranks<-matrix(NA, nrow=nrow(Amatrix), ncol=ncol(Amatrix),dimnames=dimnames(Amatrix))
            }
            
            #container for feature scores
            ssscores<-ssranks
            
            #for each lp, calculate the L2 distance from each row of A to lp[i], rank
            for (i in seq_along(lp)){
              sstat <- apply(Arowmax, 1, function(x) sqrt(t(x-lp[[i]])%*%(x-lp[[i]])))
              ssranks[,i] <- rank(sstat, ties.method="first")
              ssscores[,i] <- sstat
            }
            
            if(threshold=="all"){
              ssgenes.th <- .patternMarkers_all(ssranks)
            } else if(threshold=="cut"){
              ssgenes.th <- .patternMarkers_cut(ssranks)
            }
            
            return(list("PatternMarkers"=ssgenes.th,
                        "PatternRanks"=ssranks,
                        "PatternScores"=ssscores))
            
          }


.patternMarkers_all <- function(ssranks) {
  pIndx<-apply(ssranks,1,which.min)
  pNames<-setNames(seq_along(colnames(ssranks)), colnames(ssranks))
  ssgenes.th <- lapply(pNames,function(x) names(pIndx[pIndx==x]))
  
  #sort genes by rank for output
  for (i in seq_along(ssgenes.th)){
    order <- names(sort(ssranks[,i]))
    ssgenes.th[[i]] <- intersect(order, ssgenes.th[[i]])
  }
  return(ssgenes.th)
}


.patternMarkers_cut <- function(ssranks) {
  ssgenes.th <- list()
  for (i in seq_along(colnames(ssranks))){
    sortSim <- names(sort(ssranks[,i], decreasing = FALSE))
    #first intra-pattern rank that is worse than inter-pattern rank
    geneThresh <- min(which(ssranks[sortSim, i] > apply(ssranks[sortSim,], 1, min)))
    markerGenes <- sortSim[1:geneThresh-1]
    ssgenes.th[[i]] <- markerGenes
  }
  names(ssgenes.th) <- colnames(ssranks)
  
  return(ssgenes.th)
}
