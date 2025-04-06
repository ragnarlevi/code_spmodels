# Set seed for reproducibility
set.seed(123)
source("utils.R")
source("zip.R")
source("simulate_data.R")

#########################
# 1. Define Expected Gradient Functions
#########################

# Function inside the E-step
h_mu_func <- function(z, y, pi, mu) {
  if (y == 0) {
    den <- pi * exp(z * mu) + (1 - pi)
    h_val <- ifelse(is.finite(den), -((1 - pi) * z) / den, 0)
    return(h_val)
  } else {
    return((y / mu) - z)
  }
}

h_pi_func <- function(z, y, pi, mu) {
  if (y == 0) {
    den <- pi * exp(z * mu) + (1 - pi)
    h_val <- ifelse(is.finite(den), -((1 - pi) * z) / den, 0)
    return(h_val)
  } else {
    return((y / mu) - z)
  }
}

# Function to compute the likelihood
zip_likelihood_func <- function(z, y, pi, mu) {
  if (y == 0) {
    return(pi + (1 - pi) * exp(-z * mu))
  } else {
    return((1 - pi) * ((z * mu)^y * exp(-z * mu)) / gamma(y + 1))
  }
}




expected_pi_grad <- function(y, beta, phi, pi, mu, z_density, do_integral = TRUE) {

  # THE DO_INTEGRAL IS JUST TO TEST DERIVATIVES
  if(do_integral){
    if (y == 0) {
      integrand_num <- function(z) {
        num <- (1 - exp(-z * mu))
        denom <- pi + (1 - pi) * exp(-z * mu)
        val <- (num / denom) * (pi + (1 - pi) * exp(-z * mu)) * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
      integrand_denom <- function(z) {
        lik <- pi + (1 - pi) * exp(-z * mu)
        val <- lik * z_density(z, phi)
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
  }else{
    if (y == 0) {
        num <- (1 - exp(-1 * mu))
        denom <- pi + (1 - pi) * exp(-1 * mu)
        val <- (num / denom) 
        val[!is.finite(val)] <- 0
      return(val)
    } else {
      return(-1 / (1 - pi))
    } 
    
  }


}
expected_pi_grad <- Vectorize(expected_pi_grad, vectorize.args = c("y", "mu"))



zzip_log_lik_function <- function(y,  phi, pi, mu, z_density) {
  
  # f(y|z) for y > 0 under the zip  model:
  f_y_given_z <- function(z) {
    zip_likelihood_func(z, y, pi, mu)
  }
  
  # Numerator: integrate h(z)*f(y|z)*f(z) over z
  numerator_integrand <- function(z) {
    f_y_given_z(z) * z_density(z, phi)
  }
  
  
  
  num_result <- integrate(numerator_integrand, lower = 1e-3, subdivisions = 1000,  upper = Inf, rel.tol = 1e-8)
  
  
  lik <- num_result$value 
  return(lik)
}
zzip_log_lik_function <- Vectorize(zzip_log_lik_function, vectorize.args = c("y", "mu"))






zip_mixed <- function(claims, X, years, locs, agg_claims, A, exposure, model_type, additive, mixing, Emethod = "integration",
                         n_iter = 100, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd", 
                         optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                         batch_size = 500, param_tol = 1e-5, Q_tol = 1e-5, verbose = 0, ...){
  
  
  control_list <- list(...)
  control_list$S <- control_list$S %||%  1000
  control_list$n_nodes <- control_list$n_nodes %||%  50
  
  
  if(mixing == "gamma"){
    z_density <- gamma_density  
    rzdensity <- rgamma_density
  }else if(mixing == "ln"){
    z_density <- ln_density  
    rzdensity <- rln_density
  }else if(mixing == "ig"){
    z_density <- ig_density 
    rzdensity <- rig_density
  }
  
  
  # Adam parameters for beta and pi
  d <- ncol(X)  # dimension of beta
  
  # Get number of regions
  nr_regions <- length(unique(locs))
  nr_edges <- nr_regions*(nr_regions+1)/2
  
  
  
  # Initialize parameters with non mixing Poisson
  beta_phi <- 0
  phi <- exp(beta_phi)
  
  out_zip <- zip(claims, X, locs, years,  agg_claims, 
                 A, additive, "learn_psi", lambda = 0, exposure = exposure, max_itr = 300)
  
  beta_est <- out_zip$beta1    # beta (d-dimensional)
  pi_est <- out_zip$prop
  psi_est <- out_zip$psi     # psi (vector of length nr_regions)
  a_est <- out_zip$a
  
  
  # Default controls
  beta1_mom <- 0.9
  beta2_mom <- 0.999
  epsilon <- 1e-8
  momentum <- 0.07
  
  
  # Set gradient descent controls 
  controls <- list()
  controls$beta <- get_controls(optimizer_beta,  
                                lr = control_list$beta$lr %||%  0.01, 
                                len = d, 
                                beta1 = control_list$beta$beta1 %||%  beta1_mom, 
                                beta2 = control_list$beta$beta2 %||%  beta2_mom , 
                                epsilon = epsilon, 
                                iter = 1, 
                                momentum = control_list$beta$momentum %||%  momentum )
  
  
  controls$psi <- get_controls(optimizer_psi,  
                               lr = control_list$psi$lr %||%  0.01, 
                               len = nr_regions, 
                               beta1 = control_list$psi$beta1 %||%  beta1_mom, 
                               beta2 = control_list$psi$beta2 %||%  beta2_mom , 
                               epsilon = epsilon, 
                               iter = 1, 
                               momentum = control_list$psi$momentum %||%  momentum )
  
  
  controls$beta_phi <- get_controls(optimizer_beta_phi,  
                                    lr = control_list$beta_phi$lr %||%  0.05, 
                                    len = 1, 
                                    beta1 = control_list$beta_phi$beta1 %||%  beta1_mom, 
                                    beta2 = control_list$beta_phi$beta2 %||%  beta2_mom , 
                                    epsilon = epsilon, 
                                    iter = 1, 
                                    momentum = control_list$beta_phi$momentum %||%  momentum )
  
  controls$a <- get_controls(optimizer_a,  
                             lr = control_list$a$lr %||%  0.01, 
                             len = nr_edges, 
                             beta1 = control_list$a$beta1 %||%  beta1_mom, 
                             beta2 = control_list$a$beta2 %||%  beta2_mom , 
                             epsilon = epsilon, 
                             iter = 1, 
                             momentum = control_list$a$momentum %||%  momentum )
  
  
  controls$pi <- get_controls(optimizer_pi, 
                              lr = control_list$pi$lr %||% 0.0001, 
                              len = 1, 
                              beta1 = control_list$pi$beta1 %||%  beta1_mom,
                              beta2 = control_list$pi$beta2 %||%  beta2_mom, 
                              epsilon = epsilon, 
                              iter = 1, 
                              momentum = control_list$pi$momentum %||%  momentum)
  
  # Start the loop
  beta_nrom <- c()
  a_norm <- c()
  beta_phi_norm <- c()
  log_lik <- -Inf
  
  N <- length(claims)
  # (Assuming get_spatial_aggregate() and get_mu() are defined in utils.R)
  cat("\nStarting updates for beta, pi (Adam) and psi (Adagrad):\n")
  
  for (iter in 1:n_iter) {
    
    # 0. Prepare updates
    # Store old parameters and update
    beta_old <- beta_est
    beta_phi_old <- beta_phi
    phi <- exp(beta_phi)
    
    if(model_type == "learn_graph"){
      a_old <- a_est
      A <- get_W_from_array(a_est, nr_regions)
      psi_est <- NA
      psi_old <- NA
    }else if(model_type == "learn_psi"){
      a_est <- 1
      a_old <- 1
      psi_old <- psi_est
    }

    
    # Compute spatial aggregate effects and baseline nu
    se <- get_spatial_aggregate(locs, A, psi_est, agg_claims, years, model_type)
    nu <- as.numeric(exp(X %*% beta_est))  # baseline mu from covariates
    
    # Get full mu by combining covariate effects (nu), spatial effects, and exposure
    mu <- get_mu(nu, se$spatial_effect, exposure, additive)
    
    
    # If using SGD, sample a mini-batch
    if (sgd) {
      idx <- c()
      for(i in 1:nr_regions){
        idx <- c(idx, sample(which(locs == i), batch_size, replace = FALSE))
      }
      claims_batch   <- claims[idx]
      X_batch        <- X[idx, , drop = FALSE]
      exposure_batch <- exposure[idx]
      locs_batch     <- locs[idx]
      years_batch <- years[idx]
      mu_batch <- mu[idx]
      nu_batch <- nu[idx]
      spatial_effect_batch <- se$spatial_effect[idx]
      agg_effect_batch <- se$agg_effect[idx]
    } else {
      claims_batch   <- claims
      X_batch        <- X
      exposure_batch <- exposure
      locs_batch     <- locs
      years_batch <- years
      mu_batch <- mu
      nu_batch <- nu
      spatial_effect_batch <- se$spatial_effect
      agg_effect_batch <- se$agg_effect
    }
    
    
    # 1. E-steps
    # Compute expected derivative with respect to mu
    grad_mu <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, 
                          h_func = h_mu_func, likelihood_func = zip_likelihood_func, 
                          method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    
    grad_pi <- expected_pi_grad(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, TRUE)
    
    
    if(mixing == "gamma"){
      w1 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, 
                       h_func = function(z, y, pi, mu) identity(z), likelihood_func = zip_likelihood_func, method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
      w2 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, 
                       h_func = function(z, y, pi, mu) log(z), likelihood_func = zip_likelihood_func, 
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    }else if(mixing == "ig"){
      w1 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, 
                       h_func = function(z, y, pi, mu) identity(z), likelihood_func = zip_likelihood_func, 
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
      w3 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, 
                       h_func = function(z, y, pi, mu) 1/z, likelihood_func = zip_likelihood_func, 
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    } else if(mixing == "ln"){
      w4 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density, 
                       h_func = function(z, y, pi, mu) log(z)^2, likelihood_func = zip_likelihood_func, 
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    }
    
    
    
    
    # Compute gradient for beta via chain rule: sum_i (grad_mu_i * nu_i * x_i)
    grad_beta_total <- get_grad_beta(additive, grad_mu, nu_batch, exposure_batch, X_batch, spatial_effect_batch, locs_batch)
    
    # Compute psi/a gradient
    if(model_type == "learn_graph"){
      grad_a_total <- get_grad_a(additive, grad_mu, agg_claims, years_batch, locs_batch, exposure_batch, nu_batch, lambda)
    }else if (model_type == "learn_psi"){
      grad_psi_total <- get_grad_psi(additive, grad_mu, agg_effect_batch, nu_batch, exposure_batch, locs_batch, nr_regions)
    }
    
    # Compute gradient for pi (scalar)
    grad_pi_total <- get_grad_pi(grad_pi)
    
    # phi grad
    grad_beta_phi_total <- get_grad_beta_phi(phi, w1, w2, w3, w4)
    
    # scale gradient for sgd
    if(sgd){
      grad_beta_total <- grad_beta_total* (N / batch_size)  
      if(model_type == "learn_graph"){
        grad_a_total <- grad_a_total* (N / batch_size)  
      }else if (model_type == "learn_psi"){
        grad_psi_total <- grad_psi_total* (N / batch_size)  
      }
      grad_pi_total <- grad_pi_total* (N / batch_size)  
      grad_beta_phi_total <- grad_beta_phi_total* (N / batch_size)  
    }
    
    # 2. M-steps - update according to gradient method
    # update beta
    beta_out <- update_gradient(optimizer_beta, beta_est, grad_beta_total, controls$beta)
    beta_est <- beta_out$param
    controls$beta <- beta_out$controls
    
    beta_converged <-  isTRUE(all.equal(beta_est, beta_old, tolerance = param_tol)) 
    
    # update psi/a
    if(model_type == "learn_graph"){
      a_out <- update_gradient(optimizer_a, a_est, grad_a_total, controls$a)
      a_est <- a_out$param
      controls$a <- a_out$controls
      a_est <- pmax(a_est, 1e-6)
      spatial_converged <- isTRUE(all.equal(a_est, a_old, tolerance = param_tol)) 
    }else if(model_type == "learn_psi"){
      psi_out <- update_gradient(optimizer_psi, psi_est, grad_psi_total, controls$psi)
      psi_est <- psi_out$param
      controls$psi <- psi_out$controls
      psi_est <- pmax(psi_est, 1e-6)
      spatial_converged <- isTRUE(all.equal(psi_est, psi_old, tolerance = param_tol))
    }
    
    
    # update pi
    pi_out <- update_gradient(optimizer_pi, pi_est, grad_pi_total, controls$pi)
    pi_est <- pi_out$param
    controls$pi <- pi_out$controls
    pi_est <- pmin(pmax(pi_est, 1e-6), 1-1e-6)
    
    
    # update phi
    beta_phi_out <- update_gradient(optimizer_beta_phi, beta_phi, grad_beta_phi_total, controls$beta_phi)
    beta_phi <- beta_phi_out$param
    controls$beta_phi <- beta_phi_out$controls
    
    
    
    # 3. Check Convergence
    if(beta_converged & spatial_converged  ){
      print("Breaking because parameters have stopped changing")
      break
    }
    
    
    mu_old <- mu
    # stop if expected likelihood has stopped decreasing
    se <- get_spatial_aggregate(locs, A, psi_est, agg_claims, years, model_type)
    nu <- as.numeric(exp(X %*% beta_est))  # baseline mu from covariates
    
    # Get full mu by combining covariate effects (nu), spatial effects, and exposure
    mu <- get_mu(nu, se$spatial_effect, exposure, additive)
    
    
    
    if(iter >1){
      log_lik_old <- log_lik
      log_lik <- sum(log(zzip_log_lik_function(claims,  exp(beta_phi), pi_est,  mu,  z_density)))
      
      if(iter > 2){
        if(isTRUE(all.equal(log_lik, log_lik_old, tolerance = Q_tol))){
          print("Breaking because log-likelihood has stopped changing")
          break
        } 
      }
      if(log_lik < log_lik_old){
        
        print("There is a decrease in the likelihood. Some numerical instabilities occuring. Maybe change SGD batch if using SGD")
      }
      
    }
    
    if(verbose >= 1){
      print(paste0("log-likelihood value  ",sum(log_lik)))
    }
    
    if(verbose >=2){
      if(model_type == "learn_graph"){
        cat(sprintf("Iteration %d: beta = [%s], a = [%s], beta_phi = %f, pi = %f \n,",
                    iter,
                    paste(round(beta_est, 4), collapse = ", "),
                    paste(round(a_est, 4), collapse = ", "),
                    beta_phi,
                    pi_est))
      }else if (model_type == "learn_psi"){
        cat(sprintf("Iteration %d: beta = [%s], psi = [%s], beta_phi = %f, pi = %f \n,",
                    iter,
                    paste(round(beta_est, 4), collapse = ", "),
                    paste(round(psi_est, 4), collapse = ", "),
                    beta_phi,
                    pi_est))
      } 
    }
    
    
    
    
    
  }
  
  
  
  
}


#########################
# 4. Test
#########################
data_sim <- simulate_claims(100, 20, "psi", TRUE, mixing = "none", model_type = "zip")

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
mixing <- "ig"


zip_mixed (claims, X, years, locs, agg_claims, A, exposure, model_type, additive, mixing, 
              n_iter = 100, lambda = 0, optimizer_beta = "adam", optimizer_psi = "adam", 
              optimizer_a = "adam", optimizer_pi = "adam", optimizer_beta_phi = "adam", sgd = FALSE,
              batch_size = 100, param_tol = 1e-9, verbose = 1)



# # Adam parameters for beta and pi
# learning_rate_beta <- 0.05
# learning_rate_pi <- 0.01
# n_iter <- 500
# d <- length(data_sim$beta1)  # dimension of beta
# beta1_mom <- 0.9
# beta2_mom <- 0.999
# epsilon <- 1e-8
# momentum <- 0.07
# 
# # Get number of regions
# nr_regions <- length(unique(data_sim$locs))
# nr_edges <- nr_regions*(nr_regions+1)/2
# beta_phi <- 0
# phi <- exp(beta_phi)
# 
# 
# lambda <- 0
# 
# 
# # Initialize parameters
# out_zip <- zip(data_sim$claims, data_sim$X, data_sim$locs, data_sim$years,  data_sim$agg_claims, 
#                data_sim$A, additive, "learn_psi", lambda = 0, exposure = data_sim$exposure, max_itr = 300)
# 
# 
# beta_est <- c(0,0,0) #out_zip$beta1    # beta (d-dimensional)
# pi_est <- 0.5#out_zip$prop           # pi (scalar)
# psi_est <- rep(0, nr_regions)# out_zip$psi     # psi (vector of length nr_regions)
# a_est <- rep(0, nr_regions*(nr_regions + 1)/2)# out_zip$a
# 
# # Initialize controls
# optimizer_beta <- "gd"
# optimizer_psi <- "gd"
# optimizer_a <- "gd"
# optimizer_pi <- "gd"
# optimizer_beta_phi <- "gd"
# sgd <- FALSE
# batch_size <- 100
# 
# 
# 
# controls <- list()
# controls$beta <- get_controls(optimizer_beta, 0.001, d, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
# controls$psi <- get_controls(optimizer_psi, 0.1, nr_regions, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
# controls$pi <- get_controls(optimizer_pi, 0.0001, 1, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
# controls$beta_phi <- get_controls(optimizer_beta_phi, 0.001, 1, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
# controls$a <- get_controls(optimizer_a, 0.001, nr_edges, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
# 
# 
# 
# 
# 
# 
# 
# t1 <- Sys.time()
# 
# 
# N <- length(claims)
# # (Assuming get_spatial_aggregate() and get_mu() are defined in utils.R)
# cat("\nStarting updates for beta, pi (Adam) and psi (Adagrad):\n")
# for (iter in 1:n_iter) {
#   t1 <- Sys.time()
#   if(model_type == "learn_graph"){
#     A <- get_W_from_array(a_est, nr_regions)
#     psi <- NA
#   }else if(model_type == "learn_psi"){
#     a_est <- 1
#   }
#   
#   
#   # If using SGD, sample a mini-batch
#   if (sgd) {
#     idx <- c()
#     for(i in 1:nr_regions){
#       idx <- c(idx, sample(which(locs == i), batch_size, replace = FALSE))
#     }
#     claims_batch   <- claims[idx]
#     X_batch        <- X[idx, , drop = FALSE]
#     exposure_batch <- exposure[idx]
#     locs_batch     <- locs[idx]
#     years_batch <- years[idx]
#   } else {
#     idx <- 1:N
#     claims_batch   <- claims
#     X_batch        <- X
#     exposure_batch <- exposure
#     locs_batch     <- locs
#     years_batch <- years
#   }
#   
#   
#   phi <- exp(beta_phi)
# 
#   # Compute spatial aggregate effects and baseline nu
#   se <- get_spatial_aggregate(locs_batch, A, psi_est, agg_claims, years_batch, model_type)
#   nu <- as.numeric(exp(X_batch %*% beta_est))  # baseline mu from covariates
#   
#   # Get full mu by combining covariate effects (nu), spatial effects, and exposure
#   mu <- get_mu(nu, se$spatial_effect, exposure_batch, additive)
#   
#   mu_batch <- mu[idx]
#   nu_batch <- nu[idx]
#   
#   # E-steps
#   # Compute expected derivative with respect to mu
#   grad_mu <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
#                         h_func = h_mu_func, likelihood_func = zip_likelihood_func, method  = "MCEM", S = 10, n_nodes = 10, rzdensity = rzdensity)
#   grad_pi <- expected_pi_grad(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, TRUE)
# 
#   
#   if(mixing == "gamma"){
#     w1 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
#                      h_func = function(z, y, pi, mu) identity(z), likelihood_func = zip_likelihood_func,  method  = "MCEM", S = 10, n_nodes = 10, rzdensity = rzdensity)
#     w2 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
#                      h_func = function(z, y, pi, mu) log(z), likelihood_func = zip_likelihood_func, method  = "MCEM", S = 10, n_nodes = 10, rzdensity = rzdensity)
#   }else if(mixing == "ig"){
#     w1 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
#                      h_func = function(z, y, pi, mu) identity(z), likelihood_func = zip_likelihood_func, method  = "MCEM", S = 10, n_nodes = 10, rzdensity = rzdensity)
#     w3 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
#                      h_func = function(z, y, pi, mu) 1/z, likelihood_func = zip_likelihood_func, method  = "MCEM", S = 10, n_nodes = 10, rzdensity = rzdensity)
#   } else if(mixing == "ln"){
#     w4 <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
#                      h_func = function(z, y, pi, mu) log(z)^2, likelihood_func = zip_likelihood_func, method  = "MCEM", S = 10, n_nodes = 10, rzdensity = rzdensity)
#   }
#   
# 
#   
#   
#   # Compute gradient for beta via chain rule: sum_i (grad_mu_i * nu_i * x_i)
#   grad_beta_total <- get_grad_beta(additive, grad_mu, nu_batch, exposure_batch, X_batch, se$spatial_effect[idx], locs_batch)
#   
#   # Compute psi/a gradient
#   if(model_type == "learn_graph"){
#     grad_a_total <- get_grad_a(additive, grad_mu, agg_claims, years_batch, locs_batch, exposure_batch, nu_batch, lambda)
#   }else if (model_type == "learn_psi"){
#     grad_psi_total <- get_grad_psi(additive,grad_mu, se$agg_effect[idx], nu_batch, exposure_batch, locs_batch)
#   }
#   
#   # Compute gradient for pi (scalar)
#   grad_pi_total <- get_grad_pi()
#   
#   # phi grad
#   grad_beta_phi_total <- get_grad_beta_phi(phi, w1, w2, w3, w4)
#   
#   # scale gradient for sgd
#   if(sgd){
#     grad_beta_total <- grad_beta_total* (N / batch_size)  
#     if(model_type == "learn_graph"){
#       grad_a_total <- grad_a_total* (N / batch_size)  
#     }else if (model_type == "learn_psi"){
#       grad_psi_total <- grad_psi_total* (N / batch_size)  
#     }
#     grad_pi_total <- grad_pi_total* (N / batch_size)  
#     grad_beta_phi_total <- grad_beta_phi_total* (N / batch_size)  
#   }
#   
#   
#   # update beta
#   beta_out <- update_gradient(optimizer_beta, beta_est, grad_beta_total, controls$beta)
#   beta_est <- beta_out$param
#   controls$beta <- beta_out$controls
#   
#   # update psi/a
#   if(model_type == "learn_graph"){
#     a_out <- update_gradient(optimizer_a, a_est, grad_a_total, controls$a)
#     a_est <- a_out$param
#     controls$a <- a_out$controls
#     a_est <- pmax(a_est, 1e-6)
#   }else if(model_type == "learn_psi"){
#     psi_out <- update_gradient(optimizer_psi, psi_est, grad_psi_total, controls$psi)
#     psi_est <- psi_out$param
#     controls$psi <- psi_out$controls
#     psi_est <- pmax(psi_est, 1e-6)
#   }
# 
#   
#   # update pi
#   pi_out <- update_gradient(optimizer_pi, pi_est, grad_pi_total, controls$pi)
#   pi_est <- pi_out$param
#   controls$pi <- pi_out$controls
#   pi_est <- pmin(pmax(pi_est, 1e-6), 1-1e-6)
#   
#   
#   # update phi
#   beta_phi_out <- update_gradient(optimizer_beta_phi, beta_phi, grad_beta_phi_total, controls$beta_phi)
#   beta_phi <- beta_phi_out$param
#   controls$beta_phi <- beta_phi_out$controls
#   
#   
#   
#   ##### Check convergence
#   # stop if expected likelihood has stopped decreasing
#   se <- get_spatial_aggregate(locs, A, psi_est, agg_claims, years, model_type)
#   nu <- as.numeric(exp(X %*% beta_est))  # baseline mu from covariates
#   
#   # Get full mu by combining covariate effects (nu), spatial effects, and exposure
#   mu <- get_mu(nu, se$spatial_effect, exposure, additive)
#   log_lik <- sum(log(log_lik_function(claims,  exp(beta_phi), pi_est,  mu,  z_density)))
#   print(paste0("log-likelihood value  ",sum(log_lik)))
#   
#   
# 
#   
#   if(model_type == "learn_graph"){
#     cat(sprintf("Iteration %d: beta = [%s], pi = %f, a = [%s], beta_phi = %f\n,",
#                 iter,
#                 paste(round(beta_est, 4), collapse = ", "),
#                 pi_est,
#                 paste(round(a_est, 4), collapse = ", "),
#                 beta_phi))
#   }else if (model_type == "learn_psi"){
#     cat(sprintf("Iteration %d: beta = [%s], pi = %f, psi = [%s], beta_phi = %f\n,",
#                 iter,
#                 paste(round(beta_est, 4), collapse = ", "),
#                 pi_est,
#                 paste(round(psi_est, 4), collapse = ", "),
#                 beta_phi))
#   }
#   
#   t2 <- Sys.time()
#   print(paste0("Time  ",t2-t1))
# 
# }
# 

