
#########################
# 1. Simulate Data
#########################



# Function to generate a sparse matrix with positive numbers
generateSparseMatrix <- function(nrows, ncols, density, minVal = 0.1, maxVal = 1) {
  set.seed(123)
  library(Matrix)
  # Ensure that density is between 0 and 1
  if (density < 0 || density > 1) {
    stop("Density must be between 0 and 1.")
  }
  
  # Use rsparsematrix to create the sparse matrix
  # The rand.x function generates positive numbers (from runif)
  sparse_mat <- rsparsematrix(nrow = nrows,
                              ncol = ncols,
                              density = density,
                              symmetric = TRUE,
                              rand.x = function(n) runif(n, min = minVal, max = maxVal))
  
  diag(sparse_mat) <- max(sparse_mat)*runif(nrows, 1.2, 1.5)
  return(sparse_mat)
}


simulate_claims <- function(n, years, spatial_type, additive, mixing, area = 5, model_type = "poisson", density = 0.4, exposure_lambda = 0, seed = 123){
  
  A1_true <- as.matrix(generateSparseMatrix(area, area, density))
  
  set.seed(seed)
  
  # store claim data
  while(TRUE){
    locs <- sample(1:area, n, replace = TRUE)
    if(length(unique(locs)) == area){
      break
    }
  }
  claims <- matrix(0, nrow = n, ncol = years+1)
  mu <- matrix(0, nrow = n, ncol = years)
  phi <- matrix(0, nrow = n, ncol = years)
  z <- matrix(0, nrow = n, ncol = years)
  years_vec <- rep(1:years, each = n)
  locs_vec <- c()
  
  exposure <- rpois(years*n, lambda = exposure_lambda) + 1
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent < 0.4] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),
                  rnorm(n*years, 0.1, sd = 0.1),
                  rnorm(n*years, -0.3, sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1, -1, 1)
  
  
  
  psi_true <- runif(area, 1, 3)
  
  prop_true <- 0.7
  p <- length(unique(locs))
  
  if(mixing == "none"){
    beta2_true <- 0
  } else if(mixing == "gamma"){
    beta2_true <- 1
  } else if(mixing == "ig"){
    beta2_true <- 0.5
  } else if(mixing == "ln"){
    beta2_true <- 0
  }
  
  # Generate data over time
  for(t in 1:years){
    if(spatial_type == "psi"){
      if(additive){
        mu[, t] <- (exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) +
                      psi_true[locs] * ( (A1_true %*% y_latent[, t]) )[locs]) *
          exposure[((t-1)*n + 1):(t*n)]
      } else {
        mu[, t] <- (1 + psi_true[locs] * ( (A1_true %*% y_latent[, t]) )[locs]) *
          exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
          exposure[((t-1)*n + 1):(t*n)]
      }
    } else if(spatial_type == "graph") {
      if(additive){
        mu[, t] <- (exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) +
                      ( (A1_true %*% y_latent[, t]))[locs]) *
          exposure[((t-1)*n + 1):(t*n)]
      } else {
        mu[, t] <- (1 + ( (A1_true %*% y_latent[, t]))[locs]) *
          exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
          exposure[((t-1)*n + 1):(t*n)]
      }
    } else if(spatial_type == "none") {
      mu[, t] <- exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
        exposure[((t-1)*n + 1):(t*n)]
    }
    
    if(mixing == "none"){
      z[, t] <- 1
    } else if(mixing == "gamma"){
      z[, t] <- rgamma(n, exp(beta2_true), exp(beta2_true))
    } else if(mixing == "ig"){
      z[, t] <- actuar::rinvgauss(n, mean = 1, shape = exp(beta2_true)^2)
    } else if(mixing == "ln"){
      z[, t] <- rlnorm(n, meanlog = -exp(beta2_true)^2/2, sdlog = exp(beta2_true))
    }
    
    u <- runif(n)
    if(model_type == "poisson"){
      claims[, t+1] <- rpois(n, mu[, t] * z[, t])
    }else if(model_type == "zip"){
      claims[, t+1] <- (u >= prop_true) * rpois(n, mu[, t] * z[, t])
    } else if(model_type == "hurdle"){
      for(i in 1:n){
        # Decide if observation is a hurdle (y=0) or positive outcome (y>0)
        if (runif(1) <= prop_true) {
          claims[i,t+1] <- 0
        } else {
          lambda_val <- mu[i,t]* z[i,t]
          # Rejection sampling for truncated Poisson (ensure y > 0)
          repeat {
            y_sample <- rpois(1, lambda = lambda_val)
            if (y_sample > 0) break
          }
          claims[i,t+1] <- y_sample
        }
      }
    }

    locs_vec <- c(locs_vec, locs)
  }
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  return(list(X = X_mat1, claims = claims, z = z, mu = mu, phi = phi, A = A1_true, psi = psi_true, 
              beta2 = beta2_true, beta1 = beta1_true, agg_claims = y_latent,
              years = years_vec, locs = locs_vec, exposure = exposure, 
              a = A1_true[upper.tri(A1_true, TRUE)] ))
}

simulate_hurdle_claims <- function(n, years, spatial_type, additive, mixing, area = 5){
  set.seed(123)
  
  # store claim data
  locs <- sample(1:area, n, replace = TRUE)
  claims <- matrix(0, nrow = n, ncol = years+1)
  mu <- matrix(0, nrow = n, ncol = years)
  phi <- matrix(0, nrow = n, ncol = years)
  z <- matrix(0, nrow = n, ncol = years)
  years_vec <- rep(1:years, each = n)
  locs_vec <- c()
  
  exposure <- rep(1, years*n)
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent < 0.4] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),
                  rnorm(n*years, 0.1, sd = 0.1),
                  rnorm(n*years, -0.3, sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1, -1, 1)
  
  A1_true <- outer(0.2*runif(area, 0, 2), 0.2*runif(area, 0, 2))
  mask_true <- matrix(rbinom(area*area, size = 1, 0.3), nrow = area, ncol = area)
  diag(mask_true) <- TRUE
  A1_true <- A1_true * mask_true
  A1_true <- A1_true + t(A1_true)
  diag(A1_true) <- 0.1
  
  psi_true <- runif(area, 1, 3)
  
  prop_true <- 0.7
  p <- length(unique(locs))
  beta2_true <- 0
  
  # Generate data over time
  for(t in 1:years){
    if(spatial_type == "psi"){
      if(additive){
        mu[, t] <- (exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) +
                      psi_true[locs] * ( (A1_true %*% y_latent[, t]) )[locs]) *
          exposure[((t-1)*n + 1):(t*n)]
      } else {
        mu[, t] <- (1 + psi_true[locs] * ( (A1_true %*% y_latent[, t]) )[locs]) *
          exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
          exposure[((t-1)*n + 1):(t*n)]
      }
    } else if(spatial_type == "graph") {
      if(additive){
        mu[, t] <- (exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) +
                      (15 * (A1_true %*% y_latent[, t]))[locs]) *
          exposure[((t-1)*n + 1):(t*n)]
      } else {
        mu[, t] <- (1 + (15 * (A1_true %*% y_latent[, t]))[locs]) *
          exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
          exposure[((t-1)*n + 1):(t*n)]
      }
    } else if(spatial_type == "none") {
      mu[, t] <- exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
        exposure[((t-1)*n + 1):(t*n)]
    }
    
    if(mixing == "none"){
      z[, t] <- 1
    } else if(mixing == "gamma"){
      z[, t] <- rgamma(n, exp(beta2_true), exp(beta2_true))
    } else if(mixing == "ig"){
      z[, t] <- actuar::rinvgauss(n, mean = 1, shape = exp(beta2_true)^2)
    } else if(mixing == "ln"){
      z[, t] <- rlnorm(n, meanlog = -exp(beta2_true)^2/2, sdlog = exp(beta2_true))
    }
    
    u <- runif(n)
    for(i in 1:n){
      # Decide if observation is a hurdle (y=0) or positive outcome (y>0)
      if (runif(1) <= prop_true) {
        claims[i,t+1] <- 0
      } else {
        lambda_val <- mu[i,t]* z[i,t]
        # Rejection sampling for truncated Poisson (ensure y > 0)
        repeat {
          y_sample <- rpois(1, lambda = lambda_val)
          if (y_sample > 0) break
        }
        claims[i,t+1] <- y_sample
      }
    }
    
    
    #claims[,t+1] <- ((u>=prop_true)*qpois(runif(n, dpois(0,mu[,t]* z[,t])+1e-5, 1), mu[,t]* z[,t]))
    
    
    
    
    
    locs_vec <- c(locs_vec, locs)
  }
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  return(list(X = X_mat1, claims = claims, z = z, mu = mu, phi = phi, A = A1_true, psi = psi_true, 
              beta2 = beta2_true, beta1 = beta1_true, agg_claims = y_latent,
              years = years_vec, locs = locs_vec, exposure = exposure, 
              a = A1_true[upper.tri(A1_true, TRUE)] ))
}



evaluate_matrix_recovery <- function(A, A_est, threshold = 1e-3) {
  # Check if matrices have identical dimensions
  if (!all(dim(A) == dim(A_est))) {
    stop("A and A_est must have the same dimensions.")
  }
  
  # --- Method 1: Relative Reconstruction Error ---
  # Frobenius norm error between A and its estimate A_est
  relative_error <- norm(A - A_est, type = "F") / norm(A, type = "F")
  
  # --- Method 2: Link Recovery Metrics ---
  # In A, consider an entry to represent a "link" if it is nonzero
  true_links <- (A != 0)
  
  # For A_est, count an entry as a recovered link if its value exceeds the threshold.
  # (If negative values might occur and you only care about magnitude, you can use abs(A_est).)
  estimated_links <- (A_est > threshold)
  
  # Compute confusion counts
  TP <- sum(true_links & estimated_links)       # True links correctly uncovered
  FP <- sum(!true_links & estimated_links)      # False links (non-links predicted as links)
  FN <- sum(true_links & !estimated_links)      # Missed links
  TN <- sum(!true_links & !estimated_links)     # Correctly absent links
  
  # Calculate precision, recall, and F1 score
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1_score  <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall)
  else NA
  
  # Create a confusion matrix.
  # Rows are actual classes, and columns are predicted classes.
  confusion_matrix <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE)
  rownames(confusion_matrix) <- c("Actual Positive", "Actual Negative")
  colnames(confusion_matrix) <- c("Predicted Positive", "Predicted Negative")
  
  # Return all computed metrics as a list
  list(
    relative_error = relative_error,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    confusion_matrix = confusion_matrix
  )
}


plot_graphs <- function(A_est, A_true){
  A_est[abs(A_est) <= 1e-3 ]<- 0
  # suppose A_true and A_est are your two N×N 0/1 matrices
  df <- expand.grid(i = 1:nrow(A_true), j = 1:ncol(A_true)) %>%
    mutate(
      true = as.vector(A_true),
      est  = as.vector(A_est),
      type = case_when(
        abs(true) > 1e-3 & abs(est) > 1e-3 ~ "TP",
        abs(true) > 1e-3 & abs(est) <= 1e-3 ~ "FN",
        abs(true) <= 1e-3  & abs(est) > 1e-3 ~ "FP",
        TRUE                ~ "TN"
      )
    )
  
  true_false_graph <- ggplot(df, aes(x = j, y = i, fill = type)) +
    geom_tile(color = "grey80") +
    # flip the y-axis so row 1 appears at the top
    scale_y_reverse(breaks = 1:nrow(A_true)) +
    coord_equal() +
    scale_fill_manual(
      values = c(TP = "#1b9e77", FN = "#d95f02", FP = "#DE3163", TN = "white")
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Column (j)",
      y = "Row (i)",
      fill = "Edge type",
      title = "True vs. Estimated Adjacency (matrix layout)"
    )
  
  # suppose A_true and A_est are your two N×N 0/1 matrices
  # turn them into a long data frame with a “which” indicator
  df_true <- expand.grid(i = 1:nrow(A_true), j = 1:ncol(A_true)) %>%
    mutate(type = "True", value = as.vector(A_true))
  
  df_est  <- expand.grid(i = 1:nrow(A_est), j = 1:ncol(A_est)) %>%
    mutate(type = "Estimated", value = as.vector(A_est))
  
  df_long <- bind_rows(df_true, df_est)
  
  # now plot with facet
  side_by_side <- ggplot(df_long, aes(x = j, y = i, fill = factor(value))) +
    geom_tile(color = "grey80") +
    scale_y_reverse(breaks = 1:nrow(A_true)) +
    coord_equal() +
    scale_fill_manual(
      values = c(`0` = "white", `1` = "#1b9e77"),
      labels = c("0" = "no edge", "1" = "edge")
    ) +
    facet_wrap(~type, ncol = 2) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Column (j)",
      y = "Row (i)",
      fill = "Presence",
      title = "True vs. Estimated Adjacency Matrices"
    )
  
  
  
  return(list(true_false_graph = true_false_graph, side_by_side = side_by_side))
  
}


