# Example objective function for L-BFGS-B that uses integration
# to compute the marginal likelihood f(y) for each observation.
#
# Parameters:
#  - par: a vector of parameters, where
#       par = c(beta, psi, pi)
#    beta: regression coefficients (length d)
#    psi: spatial multipliers (length p)
#    pi: zero-inflation probability (scalar)
#
#  - claims: vector of observed claims (y_i)
#  - X: design matrix for covariates (n x d)
#  - locs: a vector of region indicators (each in {1, ..., p})
#  - exposure: exposure values (vector, length n)
#  - phi: mixing parameter for the Gamma density (assumed known)






# Set seed for reproducibility
set.seed(123)
source("utils.R")

#########################
# 1. Simulate Data
#########################

simulate_zip_claims <- function(n, years, spatial_type, additive, mixing, area = 5){
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
  beta2_true <- -0.5
  
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
    claims[, t+1] <- (u >= prop_true) * rpois(n, mu[, t] * z[, t])
    
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

data_sim <- simulate_zip_claims(100, 20, "psi", TRUE, mixing = "gamma")




zip_obj <- function(y, mu, mu_old, pi_val, pi_old, phi){
  epsilon <- 1e-8
  
  gamma_density <- function(z) dgamma(z, shape = phi, rate = phi)
  
  if (y == 0) {
    integrand_num <- function(z) {
      h_val <- log(pi_val + (1-pi_val)*exp(-z*mu) + epsilon)
      h_val <- ifelse(is.finite(h_val), h_val, 0)
      lik <- pi_old + (1 - pi_old) * exp(-z * mu_old)
      val <- h_val * lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
    integrand_denom <- function(z) {
      lik <- pi_old + (1 - pi_old) * exp(-z * mu_old)
      val <- lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
  } else {
    integrand_num <- function(z) {
      h_val <- log(1-pi_val) + y*log(mu) - z*mu
      lik <- (1 - pi_old) * ((z * mu_old)^y * exp(-z * mu_old)) / gamma(y + 1)
      val <- h_val * lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
    integrand_denom <- function(z) {
      lik <- (1 - pi_old) * ((z * mu_old)^y * exp(-z * mu_old)) / gamma(y + 1)
      val <- lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
  }
  
  num <- integrate(integrand_num, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  denom <- integrate(integrand_denom, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  if (denom == 0) return(0)
  return(num / denom)
}
zip_obj <- Vectorize(zip_obj, vectorize.args = c("y", "mu", "mu_old"))

objective_function <- function(par, claims, X, years, locs, exposure, agg_claims, phi, mu_old, pi_old) {
  # Determine dimensions
  d <- ncol(X)                        # number of covariates
  p <- length(unique(locs))           # number of regions
  
  # Extract parameters:
  beta <- par[1:d]
  #psi  <- par[(d+1):(d+p)]
  #pi_val <- par[length(par)]
  
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type = "learn_psi")
  nu <- as.numeric(exp(X %*% beta))  # baseline mu from covariates
  
  # Get full mu by combining covariate effects (nu), spatial effects, and exposure
  mu <- get_mu(nu, se$spatial_effect, exposure, additive = TRUE)
  
  lik <- zip_obj(claims, mu, mu_old, pi_val, pi_old, phi)

  
  
  print(sum(lik))
  
  return(-sum(lik))
  

}



expected_h_lfbgs <- function(y, beta, phi, pi, mu, mu_old, pi_old) {
  gamma_density <- function(z) dgamma(z, shape = phi, rate = phi)
  
  if (y == 0) {
    integrand_num <- function(z) {
      den <-  pi * exp(z * mu) + (1 - pi)
        
      h_val <- ifelse(is.finite(den), -((1 - pi) * z) / den, 0)
      lik <- pi_old + (1 - pi_old) * exp(-z * mu_old)
      val <- h_val * lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
    integrand_denom <- function(z) {
      lik <- pi_old + (1 - pi_old) * exp(-z * mu_old)
      val <- lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
  } else {
    integrand_num <- function(z) {
      h_val <- (y / mu) - z
      lik <- (1 - pi_old) * ((z * mu_old)^y * exp(-z * mu_old)) / gamma(y + 1)
      val <- h_val * lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
    integrand_denom <- function(z) {
      lik <- (1 - pi_old) * ((z * mu_old)^y * exp(-z * mu_old)) / gamma(y + 1)
      val <- lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
  }
  
  num <- integrate(integrand_num, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  denom <- integrate(integrand_denom, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  if (denom == 0) return(0)
  return(num / denom)
}
expected_h_lfbgs <- Vectorize(expected_h_lfbgs, vectorize.args = c("y", "mu", "mu_old"))


expected_pi_grad <- function(y, beta, phi, pi, mu, mu_old, pi_old) {
  gamma_density <- function(z) dgamma(z, shape = phi, rate = phi)
  
  if (y == 0) {
    integrand_num <- function(z) {
      h_val <- (1-exp(-mu*z))/(pi + (1+pi)*exp(-mu*z))
      h_val <- ifelse(is.finite(h_val), h_val, 0)
      lik <- pi_old + (1 - pi_old) * exp(-z * mu_old)
      val <- h_val * lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
    integrand_denom <- function(z) {
      lik <- pi_old + (1 - pi_old) * exp(-z * mu_old)
      val <- lik * gamma_density(z)
      val[!is.finite(val)] <- 0
      return(val)
    }
    num <- integrate(integrand_num, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
    denom <- integrate(integrand_denom, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
    if (denom == 0) return(0)
    return(num / denom)
  } else {
    return(-1 / (1 - pi))
  }
}
expected_pi_grad <- Vectorize(expected_pi_grad, vectorize.args = c("y", "mu", "mu_old"))


get_grad_beta <- function(grad_mu, nu, exposure, X){
  as.numeric(t(grad_mu * nu * exposure) %*% X)
}

get_grad_psi <- function(grad_mu, exposure, agg_effect, locs, nr_regions){
  grad_psi <- grad_mu * exposure * agg_effect
  sapply(1:nr_regions, function(r) sum(grad_psi[locs == r]))
}

get_grad_pi <- function(claims, beta, phi, pi, mu, mu_old, pi_old){
  sum(expected_pi_grad(y = claims, beta = beta, phi = phi, pi = pi, mu = mu,  mu_old = mu_old, pi_old = pi_old))
}



gradient_function <- function(par, claims, X, years, locs, exposure, agg_claims, phi, mu_old, pi_old){
  
  # Determine dimensions
  d <- ncol(X)                        # number of covariates
  nr_regions <- length(unique(locs))           # number of regions
  
  # Extract parameters:
  beta <- par[1:d]
  #psi  <- par[(d+1):(d+p)]
  #pi_val <- par[length(par)]
  
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type = "learn_psi")
  nu <- as.numeric(exp(X %*% beta))  # baseline mu from covariates
  
  # Get full mu by combining covariate effects (nu), spatial effects, and exposure
  mu <- get_mu(nu, se$spatial_effect, exposure, additive = TRUE)
  grad_mu <- expected_h_lfbgs(y = claims, beta = beta, phi = phi, 
                              pi = pi_val, mu = mu, mu_old = mu_old, pi_old = pi_old)
  
  
  grad_beta_total <- get_grad_beta(grad_mu, nu, exposure, X)
  #grad_psi_total <- get_grad_psi(grad_mu, exposure, se$agg_effect, locs, nr_regions)
  #grad_pi_total <- get_grad_pi(claims, beta, phi, pi_val, mu, mu_old, pi_old)

  
  
  return(- c(grad_beta_total))
  
  
  
}






# Extract variables from simulation
claims <- data_sim$claims
X <- data_sim$X
years <- data_sim$years
locs <- data_sim$locs
agg_claims <- data_sim$agg_claims
A <- data_sim$A
exposure <- data_sim$exposure
model_type <- "learn_psi"
additive <- TRUE


# Suppose:
# d = number of covariates, p = number of regions.
d <- ncol(data_sim$X)
p <- length(unique(data_sim$locs))
# The total number of parameters is d + p + 1.
par_init <- c(rep(0, d), rep(1, p), 0.5)  # initial guess

# Lower and upper bounds:
lower_bounds <- c(rep(-Inf, d), rep(1e-6, p), 1e-6)
upper_bounds <- c(rep(Inf, d), rep(Inf, p), 1 - 1e-6)





# Extract parameters:
beta <- par_init[1:d]
psi  <- data_sim$psi#par_init[(d+1):(d+p)]
pi_old <- 0.7#par_init[length(par_init)]
pi_val <- 0.7


par <- par_init
for(iter in 1:200){
  
  # Extract parameters:
  beta <- par[1:d]
  #psi  <- par[(d+1):(d+p)]
  #pi_old <- par[length(par)]
  
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type = "learn_psi")
  nu <- as.numeric(exp(X %*% beta))  # baseline mu from covariates
  # Get full mu by combining covariate effects (nu), spatial effects, and exposure
  mu_old <- get_mu(nu, se$spatial_effect, exposure, additive = TRUE)
  
  
  #par, claims, X, years, locs, exposure, agg_claims, phi, mu_old, pi_old
  # Now call optim with method "L-BFGS-B"
  res <- optim(par = beta, 
               fn = objective_function, 
               gr = gradient_function,
               method = "L-BFGS-B", 
               lower = lower_bounds, 
               upper = upper_bounds,
               claims = claims, 
               years = years,
               X = X, 
               locs = locs, 
               agg_claims = agg_claims,
               exposure = exposure,
               phi = exp(-0.5),
               mu_old = mu_old, 
               pi_old = pi_old,
               control = list(maxit = 2, trace = 3))
  
  par <- res$par
  print(par)
  
}

# Extract the optimized parameters:
opt_params <- res$par
beta_opt <- opt_params[1:d]
psi_opt  <- opt_params[(d+1):(d+p)]
pi_opt   <- opt_params[length(opt_params)]
cat("Optimized beta:", beta_opt, "\n")
cat("Optimized psi:", psi_opt, "\n")
cat("Optimized pi:", pi_opt, "\n")
