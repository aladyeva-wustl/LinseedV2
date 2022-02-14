library(R6)
library(linseed)
library(NMF)
library(Rcpp)
library(RcppArmadillo)
library(corpcor)
library(progress)
library(ggplot2)
library(combinat)
library(collections)
library(plotly)
library(matlib)
library(matrixcalc)
library(reshape2)

MonteCarloLinseed <- R6Class(
  "MonteCarloLinseed",
  public = list(
    filtered_samples = NULL,
    linseed_object = NULL,
    radius_ = NULL,
    optimize_iterations_ = NULL,
    search_iterations_ = NULL,
    path_ = NULL,
    analysis_name = NULL,
    dataset = NULL,
    topGenes = NULL,
    markers_dict = NULL,
    init_markers_names = NULL,
    cell_types = NULL,
    markers_list = NULL,
    p = NULL,
    k = NULL,
    samples = NULL,
    filtered_dataset = NULL,
    new_points = NULL,
    init_points = NULL,
    full_proportions = NULL,
    full_basis = NULL,
    data = NULL,
    R = NULL,
    A = NULL,
    X = NULL,
    W = NULL,
    H = NULL,
    D = NULL,
    init_X = NULL,
    init_H = NULL,
    init_W = NULL,
    optim_init_X = NULL,
    optim_init_H = NULL,
    optim_init_W = NULL,
    init_proportions_rows = NULL,
    init_proportions_ = NULL,
    unity = NULL,
    errors_ = NULL,
    optim_init_proportions_rows = NULL,
    optim_init_proportions_ = NULL,
    step_R = NULL,
    random_R = NULL,
    step_metric = NULL,
    metric_changes = NULL,
    epsilon_limit = NULL,
    epsilon_genes = NULL,
    genes_mean = NULL,
    genes_sd = NULL,
    genes_mad = NULL,
    linear_genes = NULL,
    raw_dataset = NULL,
    metric = NULL,
    top_genes = NULL,
    genes_distances = NULL,
    positive_dataset = NULL,
    original_space_dataset = NULL,
    
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
    l2_norm = function(ds) {
      #ds_ <- t(hnsc_male_svd10_filter_ct4$filtered_dataset)-t(hnsc_male_svd10_filter_ct4$init_H)%*%t(hnsc_male_svd10_filter_ct4$init_W)
      sum(apply(ds,2,function(x){sqrt(sum(x^2))}))
    },
    
    initialize = function(dataset,
                          path,
                          analysis_name,
                          cell_types,
                          filtered_samples = c(),
                          topGenes = 100000,
                          radius = 0.01,
                          optimize_iterations = 5000,
                          search_iterations = 100,
                          epsilon = 0.00001,
                          data = NULL,
                          metric="mad") {
      self$filtered_samples <- filtered_samples
      self$dataset <- dataset
      self$path_ <- path
      self$radius_ <- radius
      self$analysis_name <- analysis_name
      self$topGenes <- topGenes
      self$optimize_iterations_ <- optimize_iterations
      self$search_iterations_ <- search_iterations
      self$init_markers_names <- c()
      self$cell_types <- cell_types
      self$p <- cell_types - 1
      
      self$data <- data
      self$unity <- matrix(1, nrow = self$cell_types, ncol = 1)
      
      
      
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
      self$epsilon_limit <- sqrt(self$samples*epsilon)
      self$markers_dict <- Dict()
      self$metric <- metric
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
    },
    
    calculateDistances = function(){
      if (self$metric == "mean") {
        dataset_norm <-t(apply(self$filtered_dataset,1,function(x) { x / mean(x)}))
      } else if (self$metric == "sd") {
        dataset_norm <-t(apply(self$filtered_dataset,1,function(x) { x / sd(x)}))
      } else if (self$metric == "mad") {
        dataset_norm <-t(apply(self$filtered_dataset,1,function(x) { x / mad(x)}))
      } else {
        stop("Metric not find. Available metrics: mean, sd, mad")
      }
      genes <- nrow(dataset_norm)
      genes_distances <- matrix(0,nrow=genes,ncol=genes)
      pb <- progress_bar$new(
        format = "Calculate distance [:bar] :percent eta: :eta",
        total = genes, clear = FALSE, width= 60)
      
      for (i in 1:genes) {
        pb$tick()
        x <- as.matrix(dataset_norm[i,])
        S <- t(x) %*% x
        R_x <- t(x) / sqrt(S[1,1])
        
        genes_distances[i,] <- sqrt(colSums((t(dataset_norm)-t(R_x)%*%R_x%*%t(dataset_norm))^2))
        genes_distances[i,i] <- NaN
      }
      rownames(genes_distances) <- rownames(dataset_norm)
      colnames(genes_distances) <- rownames(dataset_norm)
      self$genes_distances <- genes_distances
    },
    
    getSvdR = function(minus_one = F, add_bias = F, k = self$cell_types) {
      self$k <- k
      V <- self$filtered_dataset
      N <- ncol(V)
      V_t <- t(V)
      if (minus_one) {
        V <- V - matrix(1, ncol=ncol(V), nrow=nrow(V))
      }
      if (add_bias) {
        
        V_t_added <- rbind(V_t,matrix(1,nrow=1,ncol=ncol(V_t)))
        svd_ <- svd(V_t_added)
        svd_k <- t(svd_$u[,1:k])
        a <- apply(svd_k,1,sum)[-1]
        a_1 <- sqrt(N - sum(a^2))
        svd_k[1,] <- (matrix(1,nrow=1,ncol=N+1) - a%*%svd_k[-1,]) / a_1 
        svd_k <- svd_k[,1:N]
        
      } else {
        
        svd_ <- svd(V_t)
        svd_k <- t(svd_$u[,1:k])
        a <- apply(svd_k,1,sum)[-1]
        a_1 <- sqrt(N - sum(a^2))
        svd_k[1,] <- (matrix(1,nrow=1,ncol=N) - a%*%svd_k[-1,]) / a_1  
      }
      
      self$R <- svd_k
      self$A <- apply(self$R,1,sum)
      self$new_points <- self$filtered_dataset %*% t(self$R)
    },
    
    setReconstructedSpace = function(){
      if (is.null(self$R)) {
        stop("No initial R. Run getSvdR() first")
      }
      self$original_space_dataset <- self$filtered_dataset
      self$filtered_dataset <- t(self$R) %*% self$R %*% t(self$filtered_dataset)
    },
    
    filterGenes = function(remove_genes=2000) {
      noise_distance <- apply(t(self$filtered_dataset)-t(self$R)%*%self$R%*%t(self$filtered_dataset),2,function(x){norm(as.matrix(x),"F")})
      selected_genes <- names(sort(noise_distance)[1:(length(noise_distance)-remove_genes)])
      self$filtered_dataset <- self$filtered_dataset[selected_genes,]
    },
    
    getX = function(point, a) {
      X <- (1 - point%*%a[-1]) / a[1]
      cbind(X,point)
    },
    
    getProjectionsPlot = function(points, dims = 2:3) {
      toPlot <- self$new_points[, dims]
      toPlot <- as.data.frame(toPlot)
      colnames(toPlot) <- c("x", "y")
      toPlot$type <- "gene"
      
      toPlotE <- (points %*% t(self$R))[, dims]
      toPlotE <- as.data.frame(toPlotE)
      colnames(toPlotE) <- c("x", "y")
      toPlotE$corner <- 1:self$cell_types
      toPlotE$type <- "predicted"
      toPlot$corner <- NA
      toPlot <- rbind(toPlot, toPlotE)
      toPlot$corner <- as.factor(toPlot$corner)
      
      axisNames <- paste0("Projection ", dims)
      pl <- ggplot(data=toPlot, aes(x=x, y=y)) +
        theme_bw() +
        labs(x=axisNames[1], y=axisNames[2])
      pl <- pl + geom_point(aes(shape=type, size=type))
      
      pl <- pl + 
        scale_shape_manual(values=c(20, 18), labels=c("gene", "predicted")) +
        scale_size_manual(values=c(1, 3), labels=c("gene", "predicted")) +
        geom_polygon(data=dplyr::filter(toPlot, type=="predicted"),
                     mapping=aes(x, y), fill=NA, color="red", lty=2) +
        geom_label(data=dplyr::filter(toPlot, type=="predicted"),
                   mapping=aes(x, y,label=corner))
      pl
    },
    guessOrder = function(predicted, actual) {
      ctn <- nrow(predicted)
      allPerms <- permn(ctn)
      
      vals <- sapply(allPerms, function(perm) {
        sum(diag(cor(t(predicted[perm, ]), t(actual))))
      })
      perm <- allPerms[[which.max(vals)]]
      return(perm)
    },
    dotPlotPropotions = function(actual, guess=T, main=NULL, showR2=T) {
      predicted <- as.matrix(self$full_proportions)
      actual <- as.matrix(actual)
      
      if (guess) {
        predicted <- predicted[self$guessOrder(predicted, actual), ]
      }
      
      colnames(predicted) <- colnames(actual)
      rownames(predicted) <- rownames(actual)
      
      xmelt <- melt(predicted)
      ymelt <- melt(actual)
      
      colnames(ymelt) <- c("Cell Type", "Sample", "Actual")
      colnames(xmelt) <- c("Cell Type", "Sample", "Predicted")
      
      total <- cbind(ymelt, xmelt[, 3, drop=F])
      
      pred <- as.numeric(predicted)
      act <- as.numeric(actual)
      
      r2 <- summary(lm(pred ~ act))$adj.r.squared
      
      pl <- ggplot(data=total, aes(x=Actual, y=Predicted, color=`Cell Type`)) +
        geom_point() + theme_bw(base_size=8) +
        theme(aspect.ratio = 1) + geom_abline(slope=1, intercept = 0, lty=2) +
        xlim(c(0, 1)) + ylim(c(0, 1))
      if (!is.null(main)) {
        pl <- pl + labs(title=main)
      }
      if (showR2) {
        subs <- substitute(italic(R)^2~"="~r2, list(r2=r2))
        pl <- pl + annotate("text", label=as.character(as.expression(subs)), parse=T, x = 0.2, y=0.9)
      }
      pl
    },
    plotErrors = function() {
      min_error <- min(self$errors_[,2])
      toPlot <- as.data.frame(self$errors_)
      plt <- ggplot(data=toPlot, aes(x=V1, y=V2)) +
        geom_line()  +
        annotate("text",x=Inf,y=Inf,label=paste("Minimal error", min_error), hjust=1.5, vjust=10) + 
        theme_bw() + xlab("Iteration") + ylab("Error")
      plt
    },
    selectInit = function() {
      constraints_ = F
      cnt_ <- 0
      cnt_proportions <- 0
      cnt_coefficients <- 0
      only_positive_genes <- rownames(self$filtered_dataset)[apply(self$new_points %*% self$R,1,function(x){all(x>0)})]
      self$positive_dataset <- self$filtered_dataset[only_positive_genes,]
      while (!constraints_) {
        cnt_ <- cnt_ + 1
        if (cnt_ > (self$cell_types*nrow(self$positive_dataset))) {
          self$init_X <- NULL
          stop(paste0("Couldn't find initial points. Negative proportions: ",cnt_proportions,
                      " Negative coefficients: ",cnt_coefficients))
        }
        self$init_proportions_rows <- sample(nrow(self$positive_dataset), self$cell_types)
        self$init_proportions_ <- self$positive_dataset[self$init_proportions_rows, ]
        rownames(self$init_proportions_) <- paste('Cell type', 1:self$cell_types)
        self$init_X <- self$init_proportions_ %*% t(self$R)
        if (all(self$init_X %*% self$R > 0)) {
          out <- tryCatch(solve(t(self$init_X %*% t(self$init_X)),self$unity), error = function(e) e)
          if (!any(class(out) == "error")) {
            #print(t(out))
            if (all(out > 0)) {
              #if (all(self$init_X %*% t(self$init_X) %*% self$unity > 0)) {
              self$init_H <- self$init_X %*% self$R
              self$init_W <- t(self$fcnnls_coefs(t(self$init_H), t(self$filtered_dataset)))
              self$init_W <- self$init_W + 0.0000000001
              constraints_ <- T  
            } else {
              cnt_coefficients <- cnt_coefficients + 1
            }
          }
        } else {
          cnt_proportions <- cnt_proportions + 1
        }
      }
    },
    optimizeInitProportions = function(iterations_=5) {
      if (is.null(self$init_X)) {
        self$selectInit()
      }
      new_init_X <- self$init_X
      new_init_H <- self$init_H
      new_init_W <- self$init_W
      
      init_error <- norm(self$filtered_dataset-self$init_W%*%self$init_H,"F")
      #init_error <- self$l2_norm(t(self$filtered_dataset-self$init_W%*%self$init_H))
      all_selections <- self$init_proportions_rows
      new_init_proportions_rows <- self$init_proportions_rows
      for (itr_ in 1:iterations_){
        #change_point <- sample(self$cell_types,1)
        for (change_point in 1:self$cell_types) {
          print(paste(itr_,change_point))
          left_points <- self$positive_dataset[-all_selections, ]
          shuffle_set <- sample(nrow(left_points), nrow(left_points))
          for (elem in shuffle_set) {
            try_points <- new_init_proportions_rows
            try_points[change_point] <- elem
            
            init_X <- self$positive_dataset[try_points,] %*% t(self$R)
            if (all(init_X %*% self$R > 0)) {
              out <- tryCatch(solve(t(init_X %*% t(init_X)),self$unity), error = function(e) e)
              if (!any(class(out) == "error")) {
                if (all(out > 0)) {
                  #if (all(init_X %*% t(init_X) %*% self$unity > 0)) {
                  init_H <- init_X %*% self$R
                  init_W <- t(self$fcnnls_coefs(t(init_H), t(self$filtered_dataset)))
                  init_W <- init_W + 0.0000000001
                  new_error <- norm(t(self$filtered_dataset)-t(init_H)%*%t(init_W),"F")
                  #new_error <- self$l2_norm(t(self$filtered_dataset)-t(init_H)%*%t(init_W))
                  if (new_error < init_error) {
                    print(new_error)
                    new_init_X <- init_X
                    new_init_H <- init_H
                    new_init_W <- init_W
                    init_error <- new_error
                    all_selections <- c(all_selections,elem)
                    new_init_proportions_rows <- try_points
                    break 
                  }
                }
              }
            }
          }
        }
      }
      self$optim_init_proportions_rows <- new_init_proportions_rows
      self$optim_init_proportions_ <- self$filtered_dataset[try_points,]
      
      self$optim_init_X <- new_init_X
      self$optim_init_W <- new_init_W
      self$optim_init_H <- new_init_H
      
      self$init_X <- new_init_X
      self$init_W <- new_init_W
      self$init_H <- new_init_H
    },
    plotRDistance = function(){
      if (is.null(self$R)) {
        stop("No initial R. Run getSvdR() first")
      }
      
      noise_distance_before <- apply(t(self$filtered_dataset)-t(self$R)%*%self$R%*%t(self$filtered_dataset),2,function(x){norm(as.matrix(x),"F")})
      
      
      toPlot <- as.data.frame(sort(noise_distance_before, decreasing = T))
      colnames(toPlot) <- "distance"
      toPlot$position <- 1:nrow(toPlot)
      toPlot$type <- "SVD"
      
      if (!is.null(self$step_R)) {
        noise_distance_after <- apply(t(self$filtered_dataset)-t(self$step_R)%*%self$step_R%*%t(self$filtered_dataset),2,function(x){norm(as.matrix(x),"F")})
        toPlotE <-  as.data.frame(sort(noise_distance_after, decreasing = T))
        colnames(toPlotE) <- "distance"
        toPlotE$position <- 1:nrow(toPlotE)
        toPlotE$type <- "Optim SVD"
        #diff <- merge(toPlot,toPlotE,by=0)
        
        toPlot <- rbind(toPlot, toPlotE) 
      }
      
      ggplot(toPlot, aes(x=position, y=distance, color=type)) + geom_point() + 
        theme_bw() + geom_hline(aes(yintercept=self$epsilon_limit,
                                    linetype = "limit"), colour= 'blue') +
        scale_linetype_manual(name = "", values = c(2, 2), 
                              guide = guide_legend(override.aes = list(color = c("blue"))))
    },
    svdPlot = function(components = 50) {
      dataFull <- self$filtered_dataset - matrix(1, ncol=ncol(self$filtered_dataset),
                                                 nrow=nrow(self$filtered_dataset))
      
      if (is.null(dataFull)) {
        stop("Full dataset appears to be NULL, something is wrong")
      }
      
      components <-
        min(components, ncol(dataFull))
      
      
      vars <- svd(dataFull)$d ^ 2
      vars <- cumsum(vars / sum(vars))
      df <-
        data.frame(
          dimension = 1:length(vars),
          variance = vars,
          filtered = FALSE
        )
      colors <- 1
      
      if (!is.null(self$linear_genes)) {
        dataFiltered <- dataFull[self$linear_genes,]
        vars <- svd(dataFiltered)$d ^ 2
        vars <- cumsum(vars / sum(vars))
        dfFiltered <-
          data.frame(
            dimension = 1:length(vars),
            variance = vars,
            filtered = TRUE
          )
        df <- rbind(df, dfFiltered)
        colors <- 2
      }
      
      ggplot(data = df, aes(x = dimension, y = variance)) +
        geom_point(aes(y = variance, color = filtered),
                   size = 0.5,
                   alpha = 1) +
        geom_line(
          aes(y = variance, color = filtered, group = filtered),
          size = 0.5,
          alpha = 1
        ) +
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
    },
    runSearch = function() {
      if (is.null(self$init_X)) {
        stop("No initial points found. Run selectInit() first")
      }
      self$X <- self$init_X
      self$W <- self$init_W
      self$H <- self$init_H
      error_ <- norm(t(self$filtered_dataset) - t(self$H) %*% t(self$W),"F")
      #error_ <- self$l2_norm(t(self$filtered_dataset) - t(self$H) %*% t(self$W))
      print(error_)
      
      pb <- progress_bar$new(
        format = "Optimization [:bar] :percent eta: :eta",
        total = self$optimize_iterations_, clear = FALSE, width= 60)
      
      constraints_ <- T
      self$errors_ <- matrix(0,ncol=2,nrow=self$optimize_iterations_+1)
      self$errors_[1,] <- c(1,error_)
      
      for (i in 1:self$optimize_iterations_) {
        pb$tick()
        points_ <- self$X[,-1]
        constraints_ <- T
        
        for (iter_ in 1:self$search_iterations_) {
          ct <- sample(self$cell_types, 1)
          x <- rnorm((self$k-1), mean = 0, sd = 1)
          u <- runif(1,0,1)^(1/(self$k-1))
          tmp_x <- (1/sqrt(sum(x^2)))
          changes_ <- ((self$radius_ * u) / tmp_x) * x
          points_[ct,] <- points_[ct,] + changes_
          X_ <- self$getX(points_, self$A)
          if (any(X_ %*% self$R < 0)) {
            constraints_ <- F
            next
          }
          if (any(solve(t(X_ %*% t(X_)),self$unity) < 0)) {
            #if (any(X_ %*% t(X_) %*% self$unity < 0)) {
            constraints_ <- F
            next
          } else {
            constraints_ <- T
          }
          
          H_ <- X_ %*% self$R
          new_error <- norm(t(self$filtered_dataset) - t(H_) %*% t(self$W),"F")
          #new_error <- self$l2_norm(t(self$filtered_dataset) - t(H_) %*% t(self$W))
          
          if (new_error < error_) {
            self$X <- X_
            self$H <- self$X %*% self$R
            break
          }
        }
        if(constraints_) {
          gradW_ <- (self$filtered_dataset - self$W%*%self$H) %*% t(self$H)
          lrW <- self$W / (self$W %*% self$H %*% t(self$H))
          self$W <- self$W + lrW * gradW_
        }
        error_ <- norm(self$filtered_dataset - self$W%*%self$H,"F")
        #error_ <- self$l2_norm(t(self$filtered_dataset - self$W%*%self$H))
        self$errors_[i+1,] <- c(i+1,error_)
        #self$D <- solve(t(self$X),self$A)
        self$D <- t(solve(t(self$X %*% t(self$X)),self$unity))
        self$full_proportions <- diag(self$D[1,]) %*% self$H
      }
    },
    
    getEpsilonApproximation = function(min_=-5,max_=-2,values_=6) {
      if (is.null(self$R)) {
        stop("No initial R. Run getSvdR() first")
      }
      self$epsilon_genes <- matrix(0,nrow=values_,ncol=3)
      colnames(self$epsilon_genes) <- c("epsilon","limit","linear_genes")
      i <- 1
      for (epsilon_ in (10^linspace(min_, max_, values_))){
        genes_ <- self$getLinearMetric(self$R, epsilon_)
        epsilon_limit <- sqrt(self$samples*epsilon_)
        self$epsilon_genes[i,] <- c(epsilon_,epsilon_limit, genes_)
        i <- i+1
      }
      self$epsilon_genes
    },
    
    getLinearMetric = function(R, epsilon_ = 0.00001){
      N <- ncol(self$filtered_dataset)
      self$epsilon_limit <- sqrt(N*epsilon_)
      vector_distance <- t(self$filtered_dataset)-t(R)%*%R%*%t(self$filtered_dataset)
      genes_list_ <- apply(vector_distance,2,function(x){
        sqrt(sum(x^2))<self$epsilon_limit
      })
      metrics <- sum(genes_list_)
      metrics  
    },
    
    updateR = function(outer_loops = 10, tries_ = 10, epsilon_ = 0.00001, option_R = 'R',
                       alpha = 0.1) {
      
      init_R <- get(option_R,self)
      init_metric <- self$getLinearMetric(init_R, epsilon_)
      self$step_R <- init_R
      self$step_metric <- init_metric
      self$metric_changes <- c(self$step_metric)
      components <- 2:self$cell_types 
      N <- self$samples
      pb <- progress_bar$new(
        format = "Updating R [:bar] :percent eta: :eta",
        total = outer_loops, clear = FALSE, width= 60)
      
      for (l in 1:outer_loops) {
        pb$tick()
        for (i in components) {
          for (j in 1:N) {
            for (try_ in 1:tries_) {
              tmp_R <- t(self$step_R)
              alpha_ <- runif(1,min=-alpha,max=alpha)
              
              #update selected j element on i component 
              tmp_R[j,i] <- tmp_R[j,i] * (1+alpha_)
              sum_const_i_ <- tmp_R[,i] %*% tmp_R[,i]
              tmp_R[,i] <- tmp_R[,i]/c(sqrt(sum_const_i_))
              
              #update rest of the R components (Gram-Schmidt)
              prev_components <- c(i)
              for (k in components[-(i-1)]) {
                sum_projection <- matrix(0,nrow=1,ncol=N)
                for (comp_ in prev_components) {
                  sum_projection <- sum_projection + Proj(t(self$step_R)[,k], tmp_R[,comp_])
                }
                tmp_R[,k] <- t(self$step_R)[,k] - sum_projection
                prev_components <- c(prev_components,k)
                #tmp_R[k,] <- (self$step_R[i,]*self$step_R[k,]) / tmp_R[i,]
                #sum_const_k_ <- tmp_R[k,] %*% tmp_R[k,]
                #tmp_R[k,] <- tmp_R[k,]/c(sqrt(sum_const_k_))
              }
              tmp_R <- tmp_R %*% diag(1 / len(tmp_R))
              tmp_R  <- t(tmp_R)
              #update R1
              a <- apply(tmp_R,1,sum)[-1]
              a_1 <- sqrt(N - sum(a^2))
              tmp_R[1,] <- (matrix(1,nrow=1,ncol=N) - a%*%tmp_R[-1,]) / a_1
              tmp_R <- t(t(tmp_R) %*% diag(1 / len(t(tmp_R))))
              
              #check metrics value
              tmp_metrics <- self$getLinearMetric(tmp_R, epsilon_)
              if (tmp_metrics > self$step_metric) {
                self$step_metric <- tmp_metrics
                self$metric_changes <- c(self$metric_changes,self$step_metric)
                self$step_R <- tmp_R
                break
              }
            }
          }
        }
      }
    }
  )
)
