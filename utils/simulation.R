
generateBasis <- function(genes,cellTypes,sd_=0.2) {
  
  data_1 <- c(rnorm(genes*2,mean=4,sd=0.75),rnorm(genes*3,mean=10,sd=1.5))
  basis <- matrix(0,nrow=genes,ncol=cellTypes)
  sds_ <- rnorm(cellTypes,mean=sd_,sd=0.02)
  
  for (i in 1:genes) {
    basis[i,1] <- data_1[sample(1:length(data_1),1)] * rnorm(1, mean=1, sd=sds_[1])
    for (j in 2:cellTypes) {
      basis[i,j] <- basis[i,1] * rnorm(1, mean=1, sd=sds_[j])  
    }
  }
  
  basis <- normalizeBetweenArrays(basis, method="quantile")
  basis <- 2^basis
  rownames(basis) <- paste("Gene",1:genes)
  colnames(basis) <- paste0("Cell type ", 1:cellTypes)
  
  basis
}



generateProportions <- function(samples,cellTypes) {
  proportions <- linseed:::sampleFromSimplexUniformly(samples, cellTypes, 100000)
  colnames(proportions) <- paste0("Sample ", 1:samples)
  rownames(proportions) <- paste0("Cell type ", 1:cellTypes)  
  
  proportions
}

## Main function to generate simulation dataset
# genes - number of genes
# samples - number of samples
# cell_types - number of cell types
# noiseDeviation - level of noise (sd)

generateData <- function(genes,samples,cellTypes,noiseDeviation=0) {
  basis <- generateBasis(genes,cellTypes)
  
  proportions <- generateProportions(samples,cellTypes)
  data <- basis %*% proportions
  
  data[data < 0] <- 0
  
  if (noiseDeviation>0) {
    noise <- matrix(rnorm(length(data), sd=noiseDeviation),
                    nrow = genes, ncol = samples)
    noised <- data + 2^noise
    noised[noised < 0] <- 0
    
    data <- noised
  }
  
  
  
  return(list(basis=basis,
              proportions=proportions,
              data=data))
}


