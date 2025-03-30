
#########################
# 1. Simulate Data
#########################

simulate_claims <- function(n, years, spatial_type, additive, mixing, area = 5, model_type = "poisson"){
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
  
  
  if(spatial_type== "graph"){
    A1_true <- 15 * A1_true
  }
  
  psi_true <- runif(area, 1, 3)
  
  prop_true <- 0.7
  p <- length(unique(locs))
  beta2_true <- 1
  
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
