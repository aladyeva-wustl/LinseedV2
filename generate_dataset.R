
plotDistanceR <- function(obj) {
  ggplot(obj$distances_R,aes(x=idx,y=dist)) + geom_point()
}

plotDistanceS <- function(obj) {
  ggplot(obj$distances_S,aes(x=idx,y=dist)) + geom_point()
}

plotProjectionGenes <- function(obj) {
  ggplot(data.frame(obj$new_points),aes(x=X2,y=X3)) + geom_point()
}

plotProjectionSamples <- function(obj) {
  ggplot(data.frame(obj$new_samples_points),aes(x=X2,y=X3)) + geom_point()
}

plotSVD <- function(obj) {
  components <- nrow(obj$V__)
  vars <- obj$d_ ^ 2
  vars <- cumsum(vars / sum(vars))
  df <-
    data.frame(
      dimension = 1:length(vars),
      variance = vars,
      filtered = FALSE
    )
  colors <- 1
  
  ggplot(data = df, aes(x = dimension, y = variance)) +
    geom_point(aes(y = variance, color = filtered),
               size = 0.5,
               alpha = 1) +
    geom_line(aes(y = variance, color = filtered, group =
                    filtered),
              size = 0.5,
              alpha = 1) +
    scale_color_manual(values = c("#999999", "#E41A1C")[1:colors]) +
    theme_minimal(base_size = 8) +
    theme(
      axis.line.x = element_line(
        colour = 'black',
        size = 0.5,
        linetype = 'solid'
      ),
      axis.line.y = element_line(
        colour = 'black',
        size = 0.5,
        linetype = 'solid'
      ),
      legend.position = c(1, 0),
      legend.justification = c(1, 0),
      legend.background = element_rect(colour =
                                         "black", size = 0.2),
      legend.key.size = unit(0.1, "in")
    ) +
    scale_x_continuous(minor_breaks = 1:components,
                       limits = c(0, components))
}

generateExample <- function(N,genes,cell_types,sampleLogSd,sampleLogMean,
                            noiseDeviation,k, addPure = F) {
  generated_data <- generateMixedData(samples = N, genes = genes, cellTypes = cell_types, 
                                      sampleLogSd = sampleLogSd, sampleLogMean = sampleLogMean, 
                                      noiseDeviation = 0)
  gene_names <- rownames(generated_data$basis)
  
  if (addPure) {
    for (i in 1:cell_types) {
      tt <- matrix(0,ncol=1,nrow=cell_types)
      tt[i,] <- 1
      generated_data$proportions <- cbind(generated_data$proportions,tt)
      generated_data$basis <- rbind(t(tt),generated_data$basis)
    }
    rownames(generated_data$basis) <- c(paste0("Pure",c(1:cell_types)),gene_names)
    
    colnames(generated_data$proportions) <- paste("Sample",c(1:ncol(generated_data$proportions)))
    pure_data <- generated_data$basis %*% generated_data$proportions
    generated_data$data <- pure_data  
    genes <- genes+cell_types
    N <- N+cell_types
  } else {
    pure_data <- generated_data$data
  }
  
  if (noiseDeviation > 0) {
    noise <- matrix(rnorm(length(generated_data$data), sd=noiseDeviation),
                    nrow = genes, ncol = N)
    noised <- generated_data$data + 2^noise
    noised[noised < 0] <- 0
    generated_data$data <- noised      
  }
  
  V <- generated_data$data
  N <- ncol(V)
  M <- nrow(V)
  V_ <- V
  for (i in 1:3000) {
    V_ <- V_ / rowSums(V_)
    t1 <- sqrt(sum((matrix(1,nrow=1,ncol=N) - colSums(V_))^2))
    V_ <- t(t(V_) / rowSums(t(V_)))  #apply(V_1,2,function(x) {x / sum(x)})
    t2 <- sqrt(sum((matrix(1,nrow=1,ncol=M) - rowSums(V_))^2))
  }
  coef_ <- t1/t2
  
  data_ <- V_
  R_0 <- matrix(1/sqrt(N),ncol=N,nrow=1)
  data_t <- t(data_)-(t(R_0) %*% R_0 %*% t(data_))
  svd_ <- svd(data_t)
  svd_k <- matrix(0,ncol=N,nrow=k)
  svd_k[1,] <- R_0
  svd_k[2:k,] <- t(svd_$u[,1:(k-1)])
  R <- svd_k
  new_points <- data_ %*% t(R)
  
  distances_R <- data.frame(sort(sqrt(apply((t(V_) - t(R) %*% R %*% t(V_))^2,2,sum)),decreasing = T))
  colnames(distances_R) <- c("dist")
  distances_R$idx <- 1:nrow(distances_R)
  
  data_ <- t(V_)
  S_0 <- matrix(1/sqrt(M),ncol=M,nrow=1)
  data_t <- t(data_)-(t(S_0) %*% S_0 %*% t(data_))
  svd_ <- svd(data_t)
  svd_s_k <- matrix(0,ncol=M,nrow=k)
  svd_s_k[1,] <- S_0
  svd_s_k[2:k,] <- t(svd_$u[,1:(k-1)])
  
  S <- svd_s_k
  new_samples_points <- t(S %*% t(data_))
  
  distances_S <- data.frame(sort(sqrt(apply((V_ - t(S) %*% S %*% V_)^2,2,sum)),decreasing = T))
  colnames(distances_S) <- c("dist")
  distances_S$idx <- 1:nrow(distances_S)
  
  V__ <- S %*% V_ %*% t(R)
  
  d_ <- svd(V__)$d
  
  
  H_real <- generated_data$proportions / rowSums(generated_data$proportions)
  X_real <- H_real %*% t(R)
  if (addPure) {
    pure_data_real <- pure_data[,paste("Sample",(ncol(V_)-k+1):ncol(V_))]
  } else {
    pure_data_real <- generated_data$basis
  }
  Omega_real <- S %*% pure_data_real
  
  
  list(generated_data = generated_data,
       pure_data = pure_data,
       N = N, M = M, k = k, V_ = V_,
       coef_ = coef_, R = R, S = S,
       new_points = new_points,
       distances_R = distances_R,
       new_samples_points = new_samples_points,
       distances_S = distances_S,
       V__ = V__, d_ = d_,
       H_real = H_real,
       X_real = X_real,
       pure_data_real = pure_data_real,
       Omega_real = Omega_real)
}



