temp = reducedDim(spe_pseudo, "PCA")
N = which(temp[,3]>50)
# [1]  43 103 150 160 177

colData(spe_pseudo[,N])

Br6522