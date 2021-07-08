library(R6)

LinseedMetadata <- R6Class(
  "LinseedMetadata",
  public = list(
    filtered_dataset = NULL,
    init_X = NULL,
    init_Omega = NULL,
    init_H = NULL,
    init_W = NULL,
    init_W__ = NULL,
    init_D_h = NULL,
    init_D_w = NULL,
    final_X = NULL,
    final_Omega = NULL,
    final_H = NULL,
    final_W = NULL,
    final_W__ = NULL,
    D_h = NULL,
    D_w = NULL,
    deconv_errors = NULL,
    deconv_errors_X = NULL,
    deconv_errors_Omega = NULL,
    lambda_errors = NULL,
    beta_errors = NULL,
    total_errors = NULL,
    R = NULL,
    S = NULL,
    filtered_samples = NULL,
    dataset = NULL,
    path_ = NULL,
    radius_ = NULL,
    analysis_name = NULL,
    topGenes = NULL,
    optimize_iterations_ = NULL,
    search_iterations_ = NULL,
    cell_types = NULL,
    lambda = NULL,
    beta = NULL,
    full_proportions = NULL,
    metric = NULL,
    j = NULL,
    
    initialize = function(linseed_object) {
      self$filtered_samples <- linseed_object$filtered_samples
      self$dataset <- linseed_object$dataset
      self$path_ <- linseed_object$path_
      self$radius_ <- linseed_object$radius
      self$analysis_name <- linseed_object$analysis_name
      self$topGenes <- linseed_object$topGenes
      self$optimize_iterations_ <- linseed_object$optimize_iterations
      self$search_iterations_ <- linseed_object$search_iterations
      self$cell_types <- linseed_object$cell_types
      
      self$filtered_dataset <- linseed_object$filtered_dataset
      self$init_X <- linseed_object$init_X
      self$init_Omega <- linseed_object$init_Omega_
      self$init_H <- linseed_object$init_H
      self$init_W <- linseed_object$init_W
      self$init_W__ <- linseed_object$init_W__
      self$final_X <- linseed_object$X
      self$final_Omega <- linseed_object$Omega
      self$final_H <- linseed_object$H
      self$final_W <- linseed_object$W
      self$final_W__ <- linseed_object$W__
      self$deconv_errors <- linseed_object$deconv_errors
      self$deconv_errors_X <- linseed_object$deconv_errors_X
      self$deconv_errors_Omega <- linseed_object$deconv_errors_Omega
      self$lambda_errors <- linseed_object$lambda_errors
      self$beta_errors <- linseed_object$beta_errors
      self$total_errors <- linseed_object$errors_
      self$R <- linseed_object$R
      self$S <- linseed_object$S
      self$init_D_h <- linseed_object$init_D_h
      self$init_D_w <- linseed_object$init_D_w
      self$D_h <- linseed_object$D_h
      self$D_w <- linseed_object$D_w
      self$full_proportions <- linseed_object$full_proportions
      self$metric <- linseed_object$metric
      self$lambda <- linseed_object$lambda_
      self$beta <- linseed_object$beta_
      self$j <- linseed_object$j
    }
  )
)
