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
    inner_iterations_X = NULL,
    inner_iterations_Omega = NULL,
    global_iterations = NULL,
    new_points = NULL,
    new_samples_points = NULL,
    orig_full_proportions = NULL,
    full_proportions = NULL,
    orig_basis = NULL,
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
    init_Omega = NULL,
    init_proportions_rows = NULL,
    init_proportions_ = NULL,
    unity = NULL,
    inits_statistics = NULL,
    errors_statistics = NULL,
    optim_init_proportions_rows = NULL,
    optim_init_proportions_ = NULL,
    genes_mean = NULL,
    genes_sd = NULL,
    genes_mad = NULL,
    top_genes = NULL,
    metric = NULL,
    positive_dataset = NULL,
    
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
                          inner_iterations_X = 1000,
                          inner_iterations_Omega = 1000,
                          global_iterations = 10) {
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
      
      self$inner_iterations_X <- inner_iterations_X
      self$inner_iterations_Omega <- inner_iterations_Omega
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
      N <- ncol(V)
      M <- nrow(V)
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
      N <- ncol(V)
      R_0 <- matrix(1/sqrt(N),ncol=N,nrow=1)
      V_t <- t(V)-(t(R_0) %*% R_0 %*% t(V))
      svd_ <- svd(V_t)
      
      svd_k <- matrix(0,ncol=N,nrow=k)
      svd_k[1,] <- R_0
      svd_k[2:k,] <- t(svd_$u[,1:(k-1)])
      
      
      self$R <- svd_k
      self$A <- matrix(apply(self$R,1,sum),ncol=1,nrow=self$cell_types)
      self$new_points <- self$V_ %*% t(self$R)
      
      #S
      V_col <- t(self$V_column)
      M <- nrow(V)
      S_0 <- matrix(1/sqrt(M),ncol=M,nrow=1)
      V_m <- t(V_col) - t(S_0) %*% S_0 %*% t(V_col)
      svd_ <- svd(V_m)
      
      svd_s_k <- matrix(0,ncol=M,nrow=k)
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
    
    selectInit = function() {
      constraints_ = F
      cnt_ <- 0
      cnt_proportions <- 0
      cnt_coefficients <- 0
      select_k <- self$cell_types
      limit_ <- select_k*nrow(self$V_row)
      while (!constraints_) {
        cnt_ <- cnt_ + 1
        if (cnt_ > limit_) {
          self$init_X <- NULL
          stop(paste0("Couldn't find initial points"))
        }
        
        self$init_proportions_rows <- sample(nrow(self$V_row), select_k)
        self$init_proportions_ <- self$V_row[self$init_proportions_rows, ]
        self$init_X <- self$init_proportions_ %*% t(self$R)
        rownames(self$init_X) <- paste('Cell type', 1:self$cell_types)
        
        out <- tryCatch(solve(t(self$init_X),self$A)[,1], error = function(e) e)
        if (!any(class(out) == "error")) {
          if (all(out>=0)) {
            constraints_ <- T
            self$init_H <- self$init_X %*% self$R
            self$init_D_h <- diag(out)
            self$init_D <- self$init_D_h * (self$M/self$N)
            self$init_Omega <- self$Sigma%*%ginv(self$init_D%*%self$init_X)
          }
        }
      }
    },
    
    optimizeInitProportions = function(iterations_=5) {
      if (is.null(self$init_X)) {
        self$selectInit()
      }
      self$inits_statistics <- NULL

      new_init_X <- self$init_X
      new_init_H <- self$init_H
      new_D_h <- self$init_D_h
      new_D <- self$init_D
      new_init_Omega <- self$init_Omega

      V_row_ <- self$S %*% self$V_row %*% t(self$R)
      
      init_error <- norm(V_row_ - self$init_Omega %*% self$init_D %*% self$init_X,"F")
      lambda_error <- self$hinge(self$init_H) 
      beta_error <- self$hinge(t(self$S) %*% self$init_Omega)
      d_error <- self$hinge(self$init_D)
      total_init_error <- init_error + self$coef_hinge_H * lambda_error + self$coef_hinge_W * beta_error + d_error
      
      print(paste("Init error:",total_init_error))
      all_selections <- self$init_proportions_rows
      new_init_proportions_rows <- self$init_proportions_rows

      genes_ <- rownames(self$filtered_dataset[self$init_proportions_rows,])
      self$inits_statistics <- rbind(self$inits_statistics,c(genes_,init_error,lambda_error,beta_error,d_error,total_init_error))
      
      for (itr_ in 1:iterations_){
        for (change_point in 1:self$cell_types) {
          print(paste(itr_,change_point))
          left_points <- self$V_row[-all_selections, ]
          shuffle_set <- sample(nrow(left_points), nrow(left_points))
          for (elem in shuffle_set) {
            try_points <- new_init_proportions_rows
            try_points[change_point] <- elem
            
            X <- self$V_row[try_points,] %*% t(self$R)
            out <- tryCatch(solve(t(X),self$A)[,1], error = function(e) e)
              if (!any(class(out) == "error")) {
                
                  H <- X %*% self$R
                  D_h <- diag(as.vector(out))
                  D <- D_h * (self$M/self$N)
                  Omega <- self$Sigma %*% ginv(D %*% X)
                  
                  new_error <- norm(V_row_ - Omega %*% D %*% X,"F")
                  new_lambda_error <- self$hinge(H) 
                  new_beta_error <- self$hinge(t(self$S) %*% Omega) 
                  new_d_error <- self$hinge(D)
                  new_total_error <- new_error + self$coef_hinge_H * new_lambda_error + self$coef_hinge_W * new_beta_error + new_d_error

                  genes_ <- rownames(self$filtered_dataset[try_points,])
                  self$inits_statistics <- rbind(self$inits_statistics,c(genes_,new_error,new_lambda_error,new_beta_error,new_d_error,new_total_error))

                  if (all(out<0)) {
                    continue
                  }
                  
                  if (new_total_error < total_init_error) {
                    print(new_total_error)
                    
                    new_init_X <- X
                    new_init_H <- H
                    new_D_h <- D_h
                    new_D <- D
                    new_init_Omega <- Omega
                    
                    total_init_error <- new_total_error
                    all_selections <- c(all_selections,elem)
                    new_init_proportions_rows <- try_points

                    
                    break 
                  }
            }
          }
        }
      }
      self$optim_init_proportions_rows <- new_init_proportions_rows
      self$optim_init_proportions_ <- self$V_row[try_points,]
      
      self$init_X <- new_init_X
      self$init_H <- new_init_H
      
      self$init_D_h <- new_D_h
      self$init_D <- new_D
      self$init_Omega <- new_init_Omega
      
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
    
    runOptimization = function(debug=FALSE, idx = NULL) {
      
      B_2 <- matrix(apply(self$S,1,sum),ncol=1,nrow=self$cell_types)
      
      self$errors_statistics <- NULL
      unity_mtx_ <- self$unity %*% t(self$unity)
      V__ <- self$S %*% self$V_ %*% t(self$R)
      M <- nrow(self$V_)
      
      self$X <-  self$init_X
      self$D_w <-  self$init_D_w
      self$Omega <- self$init_Omega_
      
      self$D_h <- diag((ginv(t(self$X)) %*% self$A %*% t(self$unity))[,1]) #diag((ginv(self$X %*% t(self$X)) %*% self$unity)[,1])
      self$H_ <- self$X %*% self$R
      self$full_proportions <- self$D_h %*% self$H_
      count_neg_props <- sum(self$full_proportions<0)
      
      self$W_ <- t(self$S) %*% self$Omega
      count_neg_basis <- sum(self$W_<0)
      
      cnt <- 0
      
      error_ <- norm(V__ - self$Omega %*% self$D_w %*% self$X,"F")
      lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R) 
      beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
      D_h_error <- self$coef_pos_D_h * self$hinge(self$D_h)
      #D_w_error <- self$coef_pos_D_w * self$hinge(self$D_w)
      prev_error <- error_ + lambda_error + beta_error + D_h_error #+ D_w_error
      
      self$errors_statistics <- rbind(self$errors_statistics,c(cnt,0,0,1,1,0,
                                                               error_,
                                                               lambda_error,
                                                               D_h_error,
                                                               beta_error,
                                                               #D_w_error,
                                                               prev_error,
                                                               count_neg_props,
                                                               count_neg_basis))
      
      
      pb <- progress_bar$new(
        format = "Optimization [:bar] :percent eta: :eta",
        total = self$global_iterations, clear = FALSE, width= 60)
      
      changed_ct <- sample(seq(1:self$cell_types),1)
      
      
      if (debug) {
        cat(paste("\n",0,":",
                  "\nInitialization",
                  "\nSelected cell type:",0,
                  "\nTotal:",prev_error,
                  "\nDeconv:",error_,
                  "\nLambda:",lambda_error,
                  "\nPos_D_h:",D_h_error,
                  "\nBeta:",beta_error,
                  "\nPos_D_w:",D_w_error,
                  "\nNegative proportions:",count_neg_props,
                  "\nNegative basis:",count_neg_W__,"\n"))
      }
      
      
      
      for (t in (1:self$global_iterations)) {
        
        for (i in (1:self$inner_iterations)) {
          
          der_X <- -2*(t(self$D_w) %*% t(self$Omega) %*% (V__ - self$Omega %*% self$D_w %*% self$X))
          #der_D_h <- -self$hinge_der(matrix(rep(diag(self$D_h),self$cell_types),ncol=self$cell_types,nrow=self$cell_types),
          #                           t(ginv(self$X) %*% ginv(self$X) %*% A %*% B ))
          
          der_D_h <- matrix(0,nrow=self$cell_types,ncol=self$cell_types)
          
          for (idx__ in which(diag(self$D_h)<0)) {
            A <- matrix(0,nrow=1,ncol=self$cell_types)
            A[,idx__] <- 1
            der_D_h <- der_D_h + (ginv(self$X) %*% B %*% A %*% ginv(self$X))
          }
          
          der_X <- der_X + self$coef_hinge_H * self$hinge_der_proportions(self$X %*% self$R,self$R) + self$coef_pos_D_h * der_D_h
          
          der_X_f <- der_X / sqrt(sum(der_X^2))
          der_X_f[,1] <- 0
          
          prev_X <- self$X
          self$X[changed_ct,] <- self$X[changed_ct,] - (self$coef_der_X*der_X_f)[changed_ct,]
          self$H_ <- self$X %*% self$R
          self$D_h <- diag((ginv(t(self$X)) %*% B %*% matrix(1,nrow=1,ncol=self$cell_types))[,1])
          
          if (i == self$inner_iterations) {
            cur_try_D_h <- 0
            while (any(self$D_h<0) & (cur_try_D_h <= self$tries_D_h)) {
              der_X <- -2*(t(self$D_w) %*% t(self$Omega) %*% (V__ - self$Omega %*% self$D_w %*% self$X))
              #der_D_h <- -self$hinge_der(matrix(rep(diag(self$D_h),self$cell_types),ncol=self$cell_types,nrow=self$cell_types),
              #                           t(ginv(self$X) %*% ginv(self$X) %*% A %*% B ))
              
              der_D_h <- matrix(0,nrow=self$cell_types,ncol=self$cell_types)
              
              for (idx__ in which(diag(self$D_h)<0)) {
                A <- matrix(0,nrow=1,ncol=self$cell_types)
                A[,idx__] <- 1
                der_D_h <- der_D_h + (ginv(self$X) %*% B %*% A %*% ginv(self$X))
              }
              
              der_X <- der_X + self$coef_hinge_H * self$hinge_der_proportions(self$X %*% self$R,self$R) + self$coef_pos_D_h * der_D_h
              
              der_X_f <- der_X / sqrt(sum(der_X^2))
              der_X_f[,1] <- 0
              
              prev_X <- self$X
              self$X[changed_ct,] <- self$X[changed_ct,] - (self$coef_der_X*der_X_f)[changed_ct,]
              self$H_ <- self$X %*% self$R
              self$D_h <- diag((ginv(t(self$X)) %*% B %*% matrix(1,nrow=1,ncol=self$cell_types))[,1])
              cur_try_D_h <- cur_try_D_h + 1
            } 
            if (cur_try_D_h > self$tries_D_h) {
              stop("D_h matrix didn't converge (has negative elements)")
            }
          }
          
          prev_changed_ct <- changed_ct
          
          if (all(self$D_h>=0)) {
            changed_ct <- sample(seq(1:self$cell_types),1)
          }
          
          self$full_proportions <- self$D_h %*% self$H_
          count_neg_props <- sum(self$full_proportions<0)
          
          error_ <- norm(V__ - self$Omega %*% self$D_w %*% self$X,"F")
          lambda_error <- self$coef_hinge_H * self$hinge(self$X %*% self$R)
          D_h_error <- self$coef_pos_D_h * self$hinge(self$D_h)
          
          new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
          cnt <- cnt + 1
          self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,i,1,0,prev_changed_ct,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   count_neg_props,
                                                                   count_neg_W__))
          
          
          if (debug) {
            cat(paste("\n",cnt,":",
                      "\nUpdating X",
                      "\nSelected cell type:",prev_changed_ct,
                      "\nTotal:",new_error,
                      "\nDeconv:",error_,
                      "\nLambda:",lambda_error,
                      "\nPos_D_h:",D_h_error,
                      "\nBeta:",beta_error,
                      "\nPos_D_w:",D_w_error,
                      "\nNegative proportions:",count_neg_props,
                      "\nNegative basis:",count_neg_W__,"\n"))
          }
          
          if (!is.null(idx)){
            if (cnt==idx) {
              stop()
            }  
          }
          
        }
        
        
        for (j in (1:self$inner_iterations_Omega)) {
          
          der_Omega <- -2*(V__ - self$Omega %*% self$D_w %*% self$X) %*% t(self$X) %*% t(self$D_w)
          + 2*t(ginv(self$Omega))%*%diag(self$cell_types)*(t(self$Omega)%*%V__%*%t(self$X))%*%t(B_1)%*%t(ginv(self$Omega))
          + 2*t(ginv(self$Omega))%*%diag(self$cell_types)*(t(self$Omega)%*%self$Omega %*% self$D_w %*% self$X%*%t(self$X))%*%t(B_1)%*%t(ginv(self$Omega))
          
          der_Omega <- der_Omega + self$coef_hinge_W * self$hinge_der_basis(t(self$S)%*%self$Omega,self$S)
          
          #der_Omega <- der_Omega - self$coef_pos_D_w * self$hinge_der(matrix(rep(diag(self$D_w),self$cell_types),ncol=self$cell_types,nrow=self$cell_types),
          #                                                            ginv(self$Omega)%*%ginv(self$Omega)%*%B_1)
          
          der_D_w <- matrix(0,nrow=self$cell_types,ncol=self$cell_types)
          
          for (idx__ in which(diag(self$D_w)<0)) {
            A <- matrix(0,nrow=1,ncol=self$cell_types)
            A[,idx__] <- 1
            der_D_w <- der_D_w + t(ginv(self$Omega)) %*% t(A) %*% t(B_2) %*% t(ginv(self$Omega))
          }
          
          der_Omega <- der_Omega + self$coef_pos_D_w * der_D_w
          
          der_Omega_f <- der_Omega / sqrt(sum(der_Omega^2))
          der_Omega_f[1,] <- 0
          prev_Omega <- self$Omega
          self$Omega[,changed_ct] <- self$Omega[,changed_ct] - (self$coef_der_Omega*der_Omega_f)[,changed_ct]
          
          self$D_w <- diag(diag(ginv(self$Omega) %*% B_1)) #diag(diag(solve(self$Omega,B_1)))
          
          if (j == self$inner_iterations_Omega) {
            cur_try_D_w <- 0
            while (any(self$D_w<0) & (cur_try_D_w < self$tries_D_w)) {
              der_Omega <- -2*(V__ - self$Omega %*% self$D_w %*% self$X) %*% t(self$X) %*% t(self$D_w)
              + 2*t(ginv(self$Omega))%*%diag(self$cell_types)*(t(self$Omega)%*%V__%*%t(self$X))%*%t(B_1)%*%t(ginv(self$Omega))
              + 2*t(ginv(self$Omega))%*%diag(self$cell_types)*(t(self$Omega)%*%self$Omega %*% self$D_w %*% self$X%*%t(self$X))%*%t(B_1)%*%t(ginv(self$Omega))
              
              der_Omega <- der_Omega + self$coef_hinge_W * self$hinge_der_basis(t(self$S)%*%self$Omega,self$S)
              
              #der_Omega <- der_Omega - self$coef_pos_D_w * self$hinge_der(matrix(rep(diag(self$D_w),self$cell_types),ncol=self$cell_types,nrow=self$cell_types),
              #                                                            ginv(self$Omega)%*%ginv(self$Omega)%*%B_1)
              
              der_D_w <- matrix(0,nrow=self$cell_types,ncol=self$cell_types)
              
              for (idx__ in which(diag(self$D_w)<0)) {
                A <- matrix(0,nrow=1,ncol=self$cell_types)
                A[,idx__] <- 1
                der_D_w <- der_D_w + t(ginv(self$Omega)) %*% t(A) %*% t(B_2) %*% t(ginv(self$Omega))
              }
              
              der_Omega <- der_Omega + self$coef_pos_D_w * der_D_w
              
              der_Omega_f <- der_Omega / sqrt(sum(der_Omega^2))
              der_Omega_f[1,] <- 0
              prev_Omega <- self$Omega
              self$Omega[,changed_ct] <- self$Omega[,changed_ct] - (self$coef_der_Omega*der_Omega_f)[,changed_ct]
              
              self$D_w <- diag(diag(ginv(self$Omega) %*% B_1))
              cur_try_D_w <- cur_try_D_w + 1
            } 
            if (cur_try_D_w > self$tries_D_w) {
              stop("D_w matrix didn't converge (has negative elements)")
            }
          }
          
          prev_changed_ct <- changed_ct
          if (all(self$D_w>=0)) {
            changed_ct <- sample(seq(1:self$cell_types),1)
          }
          
          self$W__ <- t(self$S) %*% self$Omega
          count_neg_W__ <- sum(self$W__<0)
          
          error_ <- norm(V__ - self$Omega %*% self$D_w %*% self$X,"F")
          beta_error <- self$coef_hinge_W * self$hinge(t(self$S) %*% self$Omega) 
          D_w_error <- self$coef_pos_D_w * self$hinge(self$D_w)
          new_error <- error_ + lambda_error + beta_error + D_h_error + D_w_error
          
          cnt <- cnt + 1
          self$errors_statistics <- rbind(self$errors_statistics,c(cnt,t,j,0,1,prev_changed_ct,
                                                                   error_,
                                                                   lambda_error,
                                                                   D_h_error,
                                                                   beta_error,
                                                                   D_w_error,
                                                                   new_error,
                                                                   count_neg_props,
                                                                   count_neg_W__))
          
          if (debug) {
            cat(paste("\n",cnt,":",
                      "\nUpdating Omega",
                      "\nSelected cell type:",prev_changed_ct,
                      "\nTotal:",new_error,
                      "\nDeconv:",error_,
                      "\nLambda:",lambda_error,
                      "\nPos_D_h:",D_h_error,
                      "\nBeta:",beta_error,
                      "\nPos_D_w:",D_w_error,
                      "\nNegative proportions:",count_neg_props,
                      "\nNegative basis:",count_neg_W__,"\n"))
          }
          
          if (!is.null(idx)){
            if (cnt==idx) {
              stop()
            }  
          }
          
        }
        
        
        
        
        pb$tick()
      }
      
      colnames(self$errors_statistics) <- c("idx","global_idx",
                                            "inner_idx","is_X","is_Omega","cell_type",
                                            "deconv","lambda","d_h","beta","d_w","total",
                                            "negative_props","negative_basis")
      self$H_ <- self$X %*% self$R
      self$D_h <- diag(solve(self$X %*% t(self$X),self$unity)[,1])
      self$full_proportions <- self$D_h %*% self$H_
      self$orig_full_proportions <- self$full_proportions
      self$full_proportions[self$full_proportions<0] <- 0
      self$full_proportions <- t(t(self$full_proportions) / rowSums(t(self$full_proportions)))
      self$W__ <- t(self$S) %*% self$Omega
      self$orig_W__ <- self$W__
      self$W__[self$W__<0] <- 0
      self$W__ <- t(t(self$W__) / rowSums(t(self$W__)))
      self$W_ <- self$W__ %*% self$D_w
      
    },
    
    plotErrors = function(y_limit = NULL) {
      data_toPlot <- self$errors_statistics
      toPlot <- as.data.frame(data_toPlot[,c("idx","total")])
      colnames(toPlot) <- c("idx","error")
      toPlot$type <- "Total"
      toPlotE <- as.data.frame(data_toPlot[,c("idx","deconv")])
      colnames(toPlotE) <- c("idx","error")
      toPlotE$type <- "Deconv"
      toPlot <- rbind(toPlot,toPlotE)
      toPlotE <- as.data.frame(data_toPlot[,c("idx","lambda")])
      colnames(toPlotE) <- c("idx","error")
      toPlotE$type <- "Lambda"
      toPlot <- rbind(toPlot,toPlotE)
      toPlotE <- as.data.frame(data_toPlot[,c("idx","beta")])
      colnames(toPlotE) <- c("idx","error")
      toPlotE$type <- "Beta"
      toPlot <- rbind(toPlot,toPlotE)
      toPlotE <- as.data.frame(data_toPlot[,c("idx","d_h")])
      colnames(toPlotE) <- c("idx","error")
      toPlotE$type <- "D_h"
      toPlot <- rbind(toPlot,toPlotE)
      toPlotE <- as.data.frame(data_toPlot[,c("idx","d_w")])
      colnames(toPlotE) <- c("idx","error")
      toPlotE$type <- "D_w"
      toPlot <- rbind(toPlot,toPlotE)
      
      plt <- ggplot(data=toPlot, aes(x=idx, y=error,color=type)) +
        geom_line()  +
        #annotate("text",x=Inf,y=Inf,label=paste("Minimal error", min_error), hjust=1.5, vjust=10) + 
        theme_bw() + xlab("Iteration") + ylab("Error")
      
      if (!is.null(y_limit)) {
        plt <- plt + coord_cartesian(ylim = c(0,y_limit))
      }
      plt
    },
    saveResults = function() {
      write.table(self$full_proportions,
                  file=paste0(self$path_,"/","markers_",self$analysis_name,"_proportions.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)
      colnames(self$W_) <- paste0("Cell_type_",1:self$cell_types)
      toSave <- self$W_
      toSave <- self$getFoldChange(toSave)
      toSave <- rbind(c(rep(NA,self$cell_types),round(apply(self$full_proportions,1,mean),4)),toSave)
      rownames(toSave) <- c("avg_proportions",rownames(self$filtered_dataset))
      write.table(toSave,file=paste0(self$path_,"/","markers_",self$analysis_name,"_basis_fc.tsv"),
                  sep="\t",col.names = NA, row.names = T, quote = F)
    }
  )
)