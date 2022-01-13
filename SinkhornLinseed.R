library(R6)
library(linseed)
library(NMF)
library(ggplot2)
library(combinat)
library(progress)
library(corpcor)
library(MASS)

SinkhornLinseed <- R6Class(
  "SinkhornLinseed",
  public = list(
    filtered_samples = NULL,
    filtered_dataset = NULL,
    raw_dataset = NULL,
    linseed_object = NULL,
    path_ = NULL,
    analysis_name = NULL,
    dataset = NULL,
    topGenes = NULL,
    cell_types = NULL,
    samples = NULL,
    coef_der_X = NULL,
    coef_der_Omega = NULL,
    coef_hinge_H = NULL,
    coef_hinge_W = NULL,
    coef_pos_D_w = NULL,
    coef_pos_D_h = NULL,
    global_iterations = NULL,
    new_points = NULL,
    new_samples_points = NULL,
    orig_full_proportions = NULL,
    full_proportions = NULL,
    orig_full_basis = NULL,
    full_basis = NULL,
    data = NULL,
    V_row = NULL,
    V_column = NULL,
    M = NULL,
    N = NULL,
    Sigma = NULL,
    R = NULL,
    S = NULL,
    A = NULL,
    B = NULL,
    X = NULL,
    Omega = NULL,
    W_ = NULL,
    H_ = NULL,

    D_h = NULL,
    D_w = NULL,
    D = NULL,
    D_v_row = NULL,
    D_v_column = NULL,

    init_D_w = NULL,
    init_D_h = NULL,
    init_D = NULL,
    init_X = NULL,
    init_H = NULL,
    init_W = NULL,
    init_Omega = NULL,

    unity = NULL,
    init_count_neg_props = NULL,
    init_count_neg_basis = NULL,
    count_neg_props = NULL,
    count_neg_basis = NULL,
    errors_statistics = NULL,
    genes_mean = NULL,
    genes_sd = NULL,
    genes_mad = NULL,
    top_genes = NULL,
    metric = NULL,
    
    getFoldChange = function(signatures) {
      cell_types_fc <-
        matrix(0,
               ncol = ncol(signatures),
               nrow = nrow(signatures))
      rownames(cell_types_fc) <-
        rownames(signatures)
      colnames(cell_types_fc) <-
        paste0("FC_", colnames(signatures))
      for (i in 1:ncol(cell_types_fc)) {
        cell_types_fc[, i] <- apply(signatures, 1, function(x) {
          x[i] / mean(x[-i])
        })
      }
      cbind(cell_types_fc, signatures)
    },
    
    fcnnls_coefs = function(H_0, filtered_dataset) {
      W <- (.fcnnls(H_0, filtered_dataset, pseudo = TRUE))$coef
      W
    },
    
    initialize = function(dataset,
                          path,
                          analysis_name,
                          cell_types,
                          filtered_samples = c(),
                          topGenes = 100000,
                          data = NULL,
                          metric="mad",
                          coef_der_X = 0.00001,
                          coef_der_Omega = 0.0000001,
                          coef_pos_D_w = 0.01,
                          coef_pos_D_h = 0.01,
                          coef_hinge_H = 100,
                          coef_hinge_W = 10,
                          global_iterations = 10000) {
      self$filtered_samples <- filtered_samples
      self$dataset <- dataset
      self$path_ <- path
      self$analysis_name <- analysis_name
      self$topGenes <- topGenes
      self$cell_types <- cell_types
      
      self$data <- data
      self$unity <- matrix(1, nrow = self$cell_types, ncol = 1)
      
      self$coef_der_X <- coef_der_X
      self$coef_der_Omega <- coef_der_Omega
      self$coef_hinge_H <- coef_hinge_H
      self$coef_hinge_W <- coef_hinge_W
      self$coef_pos_D_w <- coef_pos_D_w
      self$coef_pos_D_h <- coef_pos_D_h
      
      self$global_iterations <- global_iterations
      
      if (!is.null(data)) {
        input_data <- self$data
      } else {
        input_data <- self$dataset
      }
      
      if (length(self$filtered_samples) != 0) {
        self$linseed_object <- LinseedObject$new(
          input_data,
          topGenes = self$topGenes,
          samples = self$filtered_samples
        )
      } else {
        self$linseed_object <- LinseedObject$new(input_data,
                                                 topGenes = self$topGenes)
      }
      
      self$filtered_dataset <- self$linseed_object$exp$full$norm
      self$raw_dataset <- self$linseed_object$exp$full$raw
      
      self$genes_mean <- apply(self$raw_dataset,1,mean)
      self$genes_sd <- apply(self$raw_dataset,1,sd)
      self$genes_mad <- apply(self$raw_dataset,1,mad)
      
      self$samples <- ncol(self$filtered_dataset)
      self$metric <- metric

      self$N <- ncol(self$filtered_dataset)
      self$M <- nrow(self$filtered_dataset)
      
    },
    
    plotMetricHistogram = function(logScale=F,
                                   breaks=100) {
      if (self$metric == "mean") {
        toPlot <- self$genes_mean
      } else if (self$metric == "sd") {
        toPlot <- self$genes_sd
      } else if (self$metric == "mad") {
        toPlot <- self$genes_mad
      } else {
        stop("Metric not find. Available metrics: mean, sd, mad")
      }
      
      if (logScale) {
        toPlot <- log(toPlot,10)
      }
      binwidth <- round((max(toPlot)-min(toPlot))/breaks)
      toPlot <- data.frame(toPlot)
      colnames(toPlot) <- "metric"
      p <- ggplot(toPlot, aes(x=metric)) + 
        geom_histogram(binwidth = binwidth)
      p
    },
    
    selectTopGenes = function(genes_number=10000) {
      if (self$metric == "mean") {
        dataset_ <- self$genes_mean
      } else if (self$metric == "sd") {
        dataset_ <- self$genes_sd
      } else if (self$metric == "mad") {
        dataset_ <- self$genes_mad
      } else {
        stop("Metric not find. Available metrics: mean, sd, mad")
      }
      self$top_genes <- names(sort(dataset_,decreasing = T)[1:genes_number])
      self$filtered_dataset <- self$filtered_dataset[self$top_genes,]
      self$M <- nrow(self$filtered_dataset)
    },
    
    scaleDataset = function(iterations = 10000){
      V <- self$raw_dataset[rownames(self$filtered_dataset),]
      V_row <- V
      V_column <- V
      pb <- progress_bar$new(
        format = "Scaling dataset [:bar] :percent eta: :eta",
        total = iterations, clear = FALSE, width= 60)

      for (i in 1:iterations) {
        self$D_v_row <- diag(1/rowSums(V_column))
        V_row <- self$D_v_row %*% V_column
        self$D_v_column <- diag(1/rowSums(t(V_row)))
        V_column <- V_row %*% self$D_v_column
        pb$tick()
      }
      self$V_row <- V_row
      self$V_column <- V_column
    },
    
    getSvdProjections = function(k = self$cell_types) {
      V <- self$V_row
      
      #R
      R_0 <- matrix(1/sqrt(self$N),ncol=self$N,nrow=1)
      V_t <- t(V)-(t(R_0) %*% R_0 %*% t(V))
      svd_ <- svd(V_t)
      
      svd_k <- matrix(0,ncol=self$N,nrow=k)
      svd_k[1,] <- R_0
      svd_k[2:k,] <- t(svd_$u[,1:(k-1)])
      
      
      self$R <- svd_k
      self$A <- matrix(apply(self$R,1,sum),ncol=1,nrow=self$cell_types)
      self$new_points <- self$V_ %*% t(self$R)
      
      #S
      V_col <- t(self$V_column)
      S_0 <- matrix(1/sqrt(self$M),ncol=self$M,nrow=1)
      V_m <- t(V_col) - t(S_0) %*% S_0 %*% t(V_col)
      svd_ <- svd(V_m)
      
      svd_s_k <- matrix(0,ncol=self$M,nrow=k)
      svd_s_k[1,] <- S_0
      svd_s_k[2:k,] <- t(svd_$u[,1:(k-1)])
      
      self$S <- svd_s_k
      self$B <- matrix(apply(self$S,1,sum),ncol=1,nrow=self$cell_types)
      self$new_samples_points <- t(self$S %*% t(V_col))
    },

    getSvdProjectionsNew = function(k = self$cell_types){
      svd_ <- svd(self$V_row)
      self$S <- t(svd_$u[,1:k])
      self$R <- t(svd_$v[,1:k])
      self$Sigma <- diag(svd_$d[1:k])
      self$S[1,] <- -self$S[1,]
      self$R[1,] <- -self$R[1,]

      self$A <- matrix(apply(self$R,1,sum),ncol=1,nrow=self$cell_types)
      self$new_points <- self$V_row %*% t(self$R)

      self$B <- matrix(apply(self$S,1,sum),ncol=1,nrow=self$cell_types)
      self$new_samples_points <- t(self$S %*% self$V_column)
    },
    
    selectInitOmega = function() {
      restored <- t(self$S) %*% t(self$new_samples_points)
      p <- self$cell_types
      x <- t(self$new_samples_points)
      u <- rowMeans(x)
      y <- x / matrix(kronecker(colSums(x * u), rep(1, p)), nrow=p)
      
      indice <- rep(0, p)
      A <- matrix(0, nrow=p, ncol=p)
      A[p, 1] <- 1
      for (i in 1:p) {
        w <- matrix(runif(p), ncol=1)
        f <- w - A %*% pseudoinverse(A) %*% w
        f <- f / sqrt(sum(f^2))

        v <- t(f) %*% y
        indice[i] <- which.max(abs(v))
        A[, i] <- y[, indice[i]]
      }
      Ae <- restored[, indice]
      ## Omega
      self$init_Omega <- self$S %*% Ae

      ## D
      self$init_D_w <- ginv(self$init_Omega) %*% self$B
      self$init_D_h <- self$init_D_w * (self$N/self$M)
      
      ## X
      V__ <- self$S %*% self$V_row %*% t(self$R)
      self$init_X <- ginv(self$init_Omega %*% diag(self$init_D_w[,1])) %*% V__
    },

    selectInitX = function() {
      restored <- self$new_points %*% self$R
      p <- self$cell_types

      x <- t(self$new_points)
      u <- rowMeans(x)
      y <- x / matrix(kronecker(colSums(x * u), rep(1, p)), nrow=p)

      indice <- rep(0, p)
      A <- matrix(0, nrow=p, ncol=p)
      A[p, 1] <- 1
      for (i in 1:p) {
        w <- matrix(runif(p), ncol=1)
        f <- w - A %*% pseudoinverse(A) %*% w;
        f <- f / sqrt(sum(f^2))

        v <- t(f) %*% y
        indice[i] <- which.max(abs(v))
        A[, i] <- y[, indice[i]]
      }
      Ae <- restored[indice, ]
      ## X
      self$init_X <- Ae %*% t(self$R)
      ## D
      self$init_D_h <- ginv(t(self$init_X)) %*% self$A
      self$init_D_w <- self$init_D_h * (self$M/self$N)
      ## Omega
      V__ <- self$S %*% self$V_row %*% t(self$R)
      self$init_Omega <- V__ %*% ginv(diag(self$init_D_w[,1]) %*% self$init_X)

    },
    
    hinge = function(X) {
      N <- ncol(X)
      M <- nrow(X)
      h <- 0
      for (i in 1:N) {
        for (j in 1:M) {
          h <- h + max(-X[j,i],0)
        }  
      }
      return(h)
    },
    
    
    hinge_der_proportions = function(H,R){
      m <- nrow(H)
      n <- ncol(H)
      der_R <- list()
      for (c in 1:m) {
        der_loc_R <- matrix(0,nrow=m,ncol=n)
        for (i in 1:m) {
          for (j in 1:n) {
            if (H[i,j] < 0) {
              der_loc_R[i,j] <- -R[c,j]
            }
          }
        }
        der_R[[c]] <- der_loc_R
      }
      res <- matrix(0,nrow=m,ncol=m)
      for (c in 1:m) {
        mtx <- der_R[[c]]
        res[,c] <- apply(mtx,1,sum)
      }
      res[,1] <- 0
      res
    },
    
    hinge_der_basis = function(W,S){
      
      n <- ncol(W)
      res <- matrix(0,nrow=n,ncol=n)
      
      for (j in 1:n) {
        if (any(W[,j]>0)) {
          res[,j] <- -apply(t(S)[which(W[,j]<0),,drop=F],2,sum)
        }
      }
      res[1,] <- 0
      res
    },
    
    hinge_der = function(X,D) {
      N <- ncol(X)
      M <- nrow(X)
      h <- matrix(0,nrow=M,ncol=N)
      for (i in 1:N) {
        for (j in 1:M) {
          if (X[i,j] < 0) {
            h[i,j] <- -D[i,j]
          }
        }  
      }
      return(h)
    },

    plotPoints2D = function(points="init") {
      if (!points %in% c("init","current")) {
        stop("Allowed values for points are 'init', 'current'")
      }

      if (points == "init") {
        X <- self$init_X
        Omega <- self$init_Omega
        count_neg_props <- self$init_count_neg_props
        count_neg_basis <- self$init_count_neg_basis
      }
      if (points == "current") {
        X <- self$X
        Omega <- self$Omega
        count_neg_props <- self$count_neg_props
        count_neg_basis <- self$count_neg_basis
      }

      ## plot X
      toPlot <- as.data.frame(self$V_row %*% t(self$R))
      colnames(toPlot) <- c("X","Y","Z")
      colnames(X) <- c("X","Y","Z")
      pltX <- ggplot(toPlot, aes(x=Y, y=Z)) +
        geom_point() + 
        geom_polygon(data=as.data.frame(X), fill=NA, color = "green") 
        theme_minimal()
      if (!is.null(count_neg_props)) {
        pltX <- pltX + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_props / (self$cell_types*self$N),4)*100,"%"), vjust=1, hjust=1)
      }

      ## plot Omega  
      toPlot <- as.data.frame(t(self$S %*% self$V_column))
      colnames(toPlot) <- c("X","Y","Z")
      rownames(Omega) <- c("X","Y","Z")
      pltOmega <- ggplot(toPlot, aes(x=Y, y=Z)) +
        geom_point() + 
        geom_polygon(data=as.data.frame(t(Omega)), fill=NA, color = "green") +
        theme_minimal()
      
      if (!is.null(count_neg_basis)) {
        pltOmega <- pltOmega + annotate("text",  x=Inf, y = Inf, label = paste0(round(count_neg_basis / (self$cell_types*self$M),4)*100,"%"), vjust=1, hjust=1)
      }

      grid.arrange(pltX,pltOmega,nrow=1)
    },
    
    runOptimization = function(debug=FALSE, idx = NULL) {
      
      self$errors_statistics <- NULL
      V__ <- self$S %*% self$V_row %*% t(self$R)

      self$X <- self$init_X
      self$D_w <- self$init_D_w
      self$Omega <- self$init_Omega
      self$W_ <- t(self$S) %*% self$Omega %*% diag(self$D_w[,1])
      self$init_count_neg_basis <- sum(self$W_ < -1e-10)

      self$D_h <- self$init_D_h
      self$H_ <- self$X %*% self$R
      self$full_proportions <- diag(self$D_h[,1]) %*% self$H_
      self$init_count_neg_props <- sum(self$full_proportions < -1e-10)

      cnt <- 0
      error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
      lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
      beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega)
      D_h_error <- self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
      D_w_error <- self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
      prev_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
      self$errors_statistics <- rbind(self$errors_statistics, c(cnt,0,1,1,1,1,error_,
      lambda_error,D_h_error,beta_error,D_w_error,prev_error,self$init_count_neg_props,self$init_count_neg_basis))
      
      ## Stage I. Optimization of deconvolution and negative basis/proportions
      
      pb <- progress_bar$new(
        format = "Optimization. Stage I [:bar] :percent eta: :eta",
        total = self$global_iterations, clear = FALSE, width= 60)
for (t in seq(1,length.out=self$global_iterations)) {
  der_X <- -2*(t(diag(self$D_w[,1])) %*% t(self$Omega) %*% (V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X))
  der_X <- der_X + self$coef_hinge_H * hinge_der_proportions(self$X %*% self$R, self$R)
  der_X_f <- der_X
  der_X_f[,1] <- 0
  self$X <- self$X - (self$coef_der_X*der_X_f)

  self$H_ <- self$X %*% self$R
  
  self$full_proportions <- diag(self$D_h[,1]) %*% self$H_
  self$count_neg_props <- sum(self$full_proportions < -1e-10)
  
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,1,0,0,0,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  
  for (j in 1:100) {
    mult_D <-  ((t(self$Omega) %*% V__ %*% t(self$X)) / (t(self$Omega) %*% self$Omega %*% diag(self$D_w[,1]) %*% self$X %*% t(self$X)))
    self$D_w <- matrix(diag(diag(self$D_w[,1]) * mult_D),nrow=self$cell_types,ncol=1)
    self$D_w <- (self$D_w / sum(self$D_w)) * self$M
    self$D_h <- self$D_w * (self$N/self$M)
  }
  
  
  
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,0,1,0,0,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  
  der_Omega <- -2*(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X) %*% t(self$X) %*% t(diag(self$D_w[,1]))
  
  der_Omega <- der_Omega + self$coef_hinge_W * simulation_3_no_noise$hinge_der_basis(t(simulation_3_no_noise$S)%*%self$Omega, simulation_3_no_noise$S)
  der_Omega_f <- der_Omega
  der_Omega_f[1,] <- 0
  self$Omega <- self$Omega - (self$coef_der_Omega*der_Omega_f)
  
  self$W_ <- t(self$S) %*% self$Omega %*% diag(self$D_w[,1])
  self$count_neg_basis <- sum(self$W_ < -1e-10)
          
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,0,0,1,0,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  
  for (j in 1:100) {
    mult_D <-  ((t(self$Omega) %*% V__ %*% t(self$X)) / (t(self$Omega) %*% self$Omega %*% diag(self$D_w[,1]) %*% self$X %*% t(self$X)))
    self$D_w <- matrix(diag(diag(self$D_w[,1]) * mult_D),nrow=self$cell_types,ncol=1)
    self$D_w <- (self$D_w / sum(self$D_w)) * self$M
    self$D_h <- self$D_w * (self$N/self$M)
  }
  
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,0,0,0,1,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  if (debug) {
        cat(paste("\n",cnt,":",
                  "\nIteration:",t,
                  "\nTotal:",new_error,
                  "\nDeconv:",error_,
                  "\nLambda:",lambda_error,
                  "\nPos_D_h:",D_h_error,
                  "\nBeta:",beta_error,
                  "\nPos_D_w:",D_w_error,
                  "\nNegative proportions:",self$count_neg_props,
                  "\nNegative basis:",self$count_neg_basis,"\n"))
      }
  
  pb$tick()
}

      ## Stage II. Optimization of deconvolution and negative basis/proportions with sum-to-one constraints
      self$coef_der_X <- self$coef_der_X / 10

pb <- progress_bar$new(
        format = "Optimization. Stage II [:bar] :percent eta: :eta",
        total = self$global_iterations, clear = FALSE, width= 60)
for (t in seq(max(self$errors_statistics[,2])+1,length.out=self$global_iterations)) {
  der_X <- -2*(t(diag(self$D_w[,1])) %*% t(self$Omega) %*% (V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X))
  der_X <- der_X + self$coef_hinge_H * hinge_der_proportions(self$X %*% self$R, self$R) + self$coef_pos_D_h * self$D_h %*% t(t(self$X)%*%self$D_h-self$A)
  der_X_f <- der_X
  der_X_f[,1] <- 0
  self$X <- self$X - (self$coef_der_X*der_X_f)

  self$H_ <- self$X %*% self$R
  
  self$full_proportions <- diag(self$D_h[,1]) %*% self$H_
  self$count_neg_props <- sum(self$full_proportions < -1e-10)
  
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,1,0,0,0,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  
  for (j in 1:100) {
    mult_D <-  ((t(self$Omega) %*% V__ %*% t(self$X)) / (t(self$Omega) %*% self$Omega %*% diag(self$D_w[,1]) %*% self$X %*% t(self$X)))
    self$D_w <- matrix(diag(diag(self$D_w[,1]) * mult_D),nrow=self$cell_types,ncol=1)
    self$D_w <- (self$D_w / sum(self$D_w)) * self$M
    self$D_h <- self$D_w * (self$N/self$M)
  }
  
  
  
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,0,1,0,0,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  
  der_Omega <- -2*(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X) %*% t(self$X) %*% t(diag(self$D_w[,1]))
  
  der_Omega <- der_Omega + self$coef_hinge_W * simulation_3_no_noise$hinge_der_basis(t(simulation_3_no_noise$S)%*%self$Omega, simulation_3_no_noise$S) + self$coef_pos_D_w * (self$Omega%*%self$D_w-simulation_3_no_noise$B) %*% t(self$D_w)
  der_Omega_f <- der_Omega
  der_Omega_f[1,] <- 0
  self$Omega <- self$Omega - (self$coef_der_Omega*der_Omega_f)
  
  self$W_ <- t(self$S) %*% self$Omega %*% diag(self$D_w[,1])
  self$count_neg_basis <- sum(self$W_ < -1e-10)
          
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,0,0,1,0,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))
  
  for (j in 1:100) {
    mult_D <-  ((t(self$Omega) %*% V__ %*% t(self$X)) / (t(self$Omega) %*% self$Omega %*% diag(self$D_w[,1]) %*% self$X %*% t(self$X)))
    self$D_w <- matrix(diag(diag(self$D_w[,1]) * mult_D),nrow=self$cell_types,ncol=1)
    self$D_w <- (self$D_w / sum(self$D_w)) * self$M
    self$D_h <- self$D_w * (self$N/self$M)
  }
  
  
  
  error_ <- norm(V__ - self$Omega %*% diag(self$D_w[,1]) %*% self$X,"F")^2
  lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
  beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
  D_h_error <-
  self$coef_pos_D_h * norm(t(self$X)%*%self$D_h-self$A,"F")^2
  D_w_error <-
  self$coef_pos_D_w * norm(self$Omega%*%self$D_w-self$B,"F")^2
          
  new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
  cnt <- cnt + 1
  self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,0,0,0,1,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   self$count_neg_props,
                                                                   self$count_neg_basis))

  if (debug) {
        cat(paste("\n",cnt,":",
                  "\nIteration:",t,
                  "\nTotal:",new_error,
                  "\nDeconv:",error_,
                  "\nLambda:",lambda_error,
                  "\nPos_D_h:",D_h_error,
                  "\nBeta:",beta_error,
                  "\nPos_D_w:",D_w_error,
                  "\nNegative proportions:",self$count_neg_props,
                  "\nNegative basis:",self$count_neg_basis,"\n"))
  }
  
  pb$tick()
}

      colnames(self$errors_statistics) <- c("idx","iteration","is_X","is_D_X","is_Omega","is_D_Omega",
                                                    "deconv_error","lamdba_error","D_h_error",
                                                    "beta_error","D_w_error","total_error",
                                                    "neg_proportions","neg_basis")
      self$H_ <- self$X %*% self$R
      self$full_proportions <- diag(self$D_h[,1]) %*% self$H_
      self$orig_full_proportions <- self$full_proportions
      self$full_proportions[self$full_proportions < -1e-10] <- 0
      self$full_proportions <- t(t(self$full_proportions) / rowSums(t(self$full_proportions)))

      self$W_ <- t(self$S) %*% self$Omega
      self$full_basis <- self$W_ %*% diag(self$D_w[,1])
      self$orig_full_basis <- self$full_basis
      self$full_basis[self$full_basis < -1e-10] <- 0
      self$full_basis <- self$full_basis / rowSums(self$full_basis)
      
    },
    
    plotErrors = function(filter_var="is_Omega",
    variables = c("deconv_error","lamdba_error","beta_error","D_h_error","D_w_error","total_error")) {
      if (is.null(colnames(self$errors_statistics))) {
        colnames(self$errors_statistics) <- c("idx","iteration","is_X","is_D_X","is_Omega","is_D_Omega",
                                                    "deconv_error","lamdba_error","D_h_error",
                                                    "beta_error","D_w_error","total_error",
                                                    "neg_proportions","neg_basis")
      }
      toPlot <- melt(data.frame(self$errors_statistics[,c("iteration",filter_var,variables)]) %>% filter(.data[[filter_var]]==1),id.vars="iteration",measure.vars = variables)
      plt <- ggplot(toPlot,aes(x=iteration,y=value,color=variable)) +
      geom_point(size=0.2) +
      geom_line() + theme_minimal()
      
      plt
    },
    saveResults = function() {
      write.table(self$full_proportions,
                  file=paste0(self$path_,"/","markers_",self$analysis_name,"_proportions.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)
      colnames(self$full_basis) <- paste0("Cell_type_",1:self$cell_types)
      toSave <- self$full_basis
      toSave <- self$getFoldChange(toSave)
      toSave <- rbind(c(rep(NA,self$cell_types),round(apply(self$full_proportions,1,mean),4)),toSave)
      rownames(toSave) <- c("avg_proportions",rownames(self$filtered_dataset))
      write.table(toSave,file=paste0(self$path_,"/","markers_",self$analysis_name,"_basis_fc.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)
    }
  )
)