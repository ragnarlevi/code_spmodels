# Set seed for reproducibility
set.seed(123)
source("utils.R")
source("zip.R")
source("simulate_data.R")

#########################
# 1. Define Expected Gradient Functions
#########################

# Function inside the E-step
h_mu_func <- function(z, y, pi, mu, phi) {
  if (y==0) {
    num   <- -(1 - pi) * z * exp(-z * mu)
    den   <-  pi + (1 - pi) * exp(-z * mu)
    return(num/den)
  } else {
    return((y/mu) - z)
  }
}



#h_mu_func <- Vectorize(h_mu_func, vectorize.args = c("z", "y", "pi", "mu"))

h_pi_func <- function(z, y, pi, mu, phi) {
  if (y == 0) {
    num <- (1 - exp(-z * mu))
    denom <- pi + (1 - pi) * exp(-z * mu)
    h_val <- ifelse(is.finite(denom), num/denom, 0)
    return(h_val)
  } else {
    return(-1 / (1 - pi))
  }
}
#h_pi_func <- Vectorize(h_pi_func, vectorize.args = c("z", "y", "pi", "mu"))

# Function to compute the likelihood
zip_likelihood_func <- function(z, y, pi, mu, phi) {
  if (y == 0) {
    return(pi + (1 - pi) * exp(-z * mu))
  } else {
    return((1 - pi) * dpois(y, z*mu))
  }
}
#zip_likelihood_func <- Vectorize(zip_likelihood_func, vectorize.args = c("z", "y", "pi", "mu"))

log_zip_likelihood_func <- function(z, y, pi, mu, phi) {
  if (y == 0) {
    return(log(pi + (1 - pi) * exp(-z * mu)))
  } else {
    return(log(1 - pi) +  dpois(y, z*mu, log = TRUE))
  }
}

#log_zip_likelihood_func <- Vectorize(log_zip_likelihood_func, vectorize.args = c("z", "y", "pi", "mu"))




expected_pi_grad <- function(y, beta, phi, pi, mu, z_density, phi_old = phi, pi_old = pi) {


    if (y == 0) {
      integrand_num <- function(z) {
        num <- (1 - exp(-z * mu))
        denom <- pi + (1 - pi) * exp(-z * mu)
        val <- (num / denom) * (pi_old + (1 - pi_old) * exp(-z * mu)) * z_density(z, phi_old)
        val[!is.finite(val)] <- 0
        return(val)
      }
      integrand_denom <- function(z) {
        lik <- pi_old + (1 - pi_old) * exp(-z * mu)
        val <- lik * z_density(z, phi_old)
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
expected_pi_grad <- Vectorize(expected_pi_grad, vectorize.args = c("y", "mu"))



zzip_log_lik_function <- function(y, phi, pi, mu, z_density) {
  
  # f(y|z) for y > 0 under the ZIP model:
  f_y_given_z <- function(z) {
    zip_likelihood_func(z, y, pi, mu)
  }
  
  # lik <-     integrate(zip_likelihood_func(z, y, pi, mu)
  #   function(z) f_y_given_z(z) * z_density(z, phi),
  #   lower      = 1e-3,
  #   upper      = Inf,
  #   subdivisions = 1000,
  #   rel.tol    = 1e-8
  # )$value
  
  # safe integration: try integrate(), otherwise fallback to Legendre
  lik <- tryCatch({
    # numerator: ∫ f(y|z) f(z) dz
    integrate(
      function(z) f_y_given_z(z) * z_density(z, phi),
      lower      = 1e-3,
      upper      = Inf,
      subdivisions = 1000,
      rel.tol    = 1e-8
    )$value
  }, error = function(e) {
    
        warning("`integrate()` failed: ", paste0(y, " ", mu, " ", pi, " ", phi),
                "\n  falling back to Gauss–Legendre quadrature.")

    # — fallback Gauss–Legendre —
    quad     <- statmod::gauss.quad(50, kind = "legendre")
    t_nodes  <- (quad$nodes + 1) / 2
    t_weights<- quad$weights / 2

    # transform [0,1] → [0,∞): z = t/(1−t), dz = 1/(1−t)^2 dt
    z_vals   <- t_nodes / (1 - t_nodes)
    jacobian <- 1 / (1 - t_nodes)^2

    num_quad <- sum(
      t_weights *
        ( sapply(z_vals, function(z) f_y_given_z(z) * z_density(z, phi)) * jacobian ),
      na.rm = TRUE
    )

    num_quad
  })

  return(lik)
}

zzip_log_lik_function <- Vectorize(zzip_log_lik_function, vectorize.args = c("y", "mu"))






zip_mixed <- function(claims, X, years, locs, agg_claims, A, exposure, model_type, additive, mixing, Emethod = "integration",
                         n_iter = 100, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd", 
                         optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                         batch_size = 500, param_tol = 1e-5, Q_tol = 1e-5, verbose = 0, control_list = list(), 
                      do_optim = FALSE, a_known = FALSE, calc_se = TRUE, beta2_start = NULL,
                      max_itr_poisson = 300){
  
  
  
  if(a_known & do_optim){
    stop("Optim with a_known = TRUE, not supported")
  }
  
  # Start by defining two functions that go into the optim base package to
  # 1. Perform L-BFGS-B method
  # 2. Obtain Hessian for standard error calculations
  beta_optim <- function(param){
    
    beta_phi_new <- param[1]
    pi_new <- param[2]
    phi_new <- exp(beta_phi_new)
    
    
    if(model_type == "ordinary"){
      beta_new <- param[3:(ncol(X)+2)]
    }else if(model_type == "learn_psi"){
      beta_new <- param[3:(ncol(X) +2)]
      psi_new <- param[(ncol(X)+3):length(param)]
      a_new <- 1
    }else if(model_type == "learn_graph"){
      beta_new <- param[3:(ncol(X) +2)]
      a_new <- param[(ncol(X)+3):length(param)]
      A <- get_W_from_array(a_new, nr_regions)
      psi_new <- NA
    }
    
    
    # Compute spatial aggregate effects and baseline nu
    se_new <- get_spatial_aggregate(locs_batch, A, psi_est, agg_claims, years_batch, model_type)
    nu_new <- as.numeric(exp(X_batch %*% beta_new))  # baseline mu from covariates
    
    # Get full mu by combining covariate effects (nu), spatial effects, and exposure
    mu_new <- get_mu(nu_new, se_new$spatial_effect, exposure_batch, additive)
    
    
    
    lik <- expected_h(y = claims_batch, phi = phi_new, pi = pi_new, mu = mu_new, z_density = z_density,
                      h_func = log_zip_likelihood_func, likelihood_func = zip_likelihood_func,
                      method  = "integration", S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity,
                      mu_old = mu_batch, phi_old = exp(beta_phi_old), pi_old = pi_old)
    
    
    return(-sum(lik))
    
  }
  
  beta_optim_grad <- function(param){
    beta_phi_new <- param[1]
    pi_new <- param[2]
    phi_new <- exp(beta_phi_new)
    
    
    if(model_type == "ordinary"){
      beta_new <- param[3:(ncol(X)+2)]
    }else if(model_type == "learn_psi"){
      beta_new <- param[3:(ncol(X) +2)]
      psi_new <- param[(ncol(X)+3):length(param)]
      a_new <- 1
    }else if(model_type == "learn_graph" ){
      beta_new <- param[3:(ncol(X) +2)]
      a_new <- param[(ncol(X)+3):length(param)]
      A <- get_W_from_array(a_new, nr_regions)
      psi_new <- NA
    }
    
    
    
    # Compute spatial aggregate effects and baseline nu
    se_new <- get_spatial_aggregate(locs_batch, A, psi_new, agg_claims, years_batch, model_type)
    nu_new <- as.numeric(exp(X_batch %*% beta_new))  # baseline mu from covariates
    
    # Get full mu by combining covariate effects (nu), spatial effects, and exposure
    mu_new <- get_mu(nu_new, se_new$spatial_effect, exposure_batch, additive)
    
    grad_mu <- expected_h(y = claims_batch, phi = phi_new, pi = pi_new, mu = mu_new, z_density = z_density,
                          h_func = h_mu_func, likelihood_func = zip_likelihood_func,
                          method  = "integration", S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity,
                          mu_old = mu_batch, phi_old = exp(beta_phi_old), pi_old = pi_old)
    
    grad_pi <- expected_pi_grad(y = claims_batch, phi = phi_new, pi = pi_new, mu = mu_new, z_density = z_density, 
                                phi_old = exp(beta_phi_old), pi_old = pi_old)
    
    grad_beta_total <- get_grad_beta(additive, grad_mu, nu_new, exposure_batch, X_batch, se_new$spatial_effect, locs_batch)
    
    # Compute psi/a gradient
    if(model_type == "learn_graph"){
      grad_a_total <- get_grad_a(additive, grad_mu, agg_claims, years_batch, locs_batch, exposure_batch, nu_new, lambda, nr_regions)
    }else if (model_type == "learn_psi"){
      grad_psi_total <- get_grad_psi(additive, grad_mu, se_new$agg_effect, nu_new, exposure_batch, locs_batch, nr_regions)
    }
    
    # Compute gradient for pi (scalar)
    grad_pi_total <- get_grad_pi(grad_pi)
    
    # phi grad
    grad_beta_phi_total <- get_grad_beta_phi(phi_new, mixing, w1, w2, w3, w4)
    
    
    if(model_type == "learn_graph"){
      return(-c(grad_beta_phi_total, grad_pi_total, grad_beta_total, grad_a_total ))
    }else if (model_type == "learn_psi"){
      return(-c(grad_beta_phi_total, grad_pi_total, grad_beta_total, grad_psi_total))
    }
    
  }
  
  if(model_type != "learn_graph"){
    a_known <- FALSE
  }
  
  
  
  control_list$S <- control_list$S %||%  1000
  control_list$n_nodes <- control_list$n_nodes %||%  50
  
  

  

  # Get number of regions
  nr_regions <- length(unique(locs))
  nr_edges <- nr_regions*(nr_regions+1)/2
  
  
  
  # Initialize parameters with non mixing Poisson
  a1 <- Sys.time()
  out_zip <- zip(claims = claims, X1 = X, locs = locs, years = years,  agg_claims = agg_claims, 
                 A = A, additive = additive, model_type = model_type, exposure = exposure, 
                 lambda = lambda, max_itr = max_itr_poisson, a_known = a_known,  calc_hessian = FALSE)
  time_to_init <-  Sys.time() - a1
  
  # Initialize mixing parameters

  if(mixing == "gamma"){
    z_density <- gamma_density  
    rzdensity <- rgamma_density
    deriv_log_density <- deriv_log_gamma_density
  }else if(mixing == "ln"){
    z_density <- ln_density  
    rzdensity <- rln_density
    deriv_log_density <- deriv_log_ln_density
  }else if(mixing == "ig"){
    z_density <- ig_density 
    rzdensity <- rig_density
    deriv_log_density <- deriv_log_ig_density
  }

  
  # Initialize mixing parameters
  beta_phi_est <- beta2_start %||% 0.3
  phi <- exp(beta_phi_est)
  
  beta_est <- out_zip$beta1    # beta (d-dimensional)
  pi_est <- out_zip$prop
  psi_est <- out_zip$psi     # psi (vector of length nr_regions)
  a_est <- out_zip$a
  
  if(verbose>0){
    print("initial parameters:")
    print(beta_est)
    print(pi_est)    
  }
  


  # Adam parameters for beta and pi
  d <- ncol(X)  # dimension of beta  
  # Default controls for gradient descent methods
  beta1_mom <- 0.9
  beta2_mom <- 0.999
  epsilon <- 1e-8
  momentum <- 0.07

  nr_time <- length(unique(years))
  
  # Set gradient descent controls 
  controls <- list()
  controls$beta <- get_controls(optimizer_beta,  
                                lr = (control_list[["beta", exact = TRUE]]$lr )%||%  (1/nrow(X)), 
                                len = d, 
                                beta1 = control_list[["beta", exact = TRUE]]$beta1 %||%  beta1_mom, 
                                beta2 = control_list[["beta", exact = TRUE]]$beta2 %||%  beta2_mom , 
                                epsilon = epsilon, 
                                iter = 1, 
                                momentum = control_list[["beta", exact = TRUE]]$momentum %||%  momentum )

  
  
  controls$psi <- get_controls(optimizer_psi,  
                               lr = control_list$psi$lr %||%  (5/nr_time), 
                               len = nr_regions, 
                               beta1 = control_list$psi$beta1 %||%  beta1_mom, 
                               beta2 = control_list$psi$beta2 %||%  beta2_mom , 
                               epsilon = epsilon, 
                               iter = 1, 
                               momentum = control_list$psi$momentum %||%  momentum )
  
  

  controls$beta_phi <- get_controls(optimizer_beta_phi,  
                                    lr = control_list$beta_phi$lr %||%  (5/nrow(X)), 
                                    len = 1, 
                                    beta1 = control_list$beta_phi$beta1 %||%  beta1_mom, 
                                    beta2 = control_list$beta_phi$beta2 %||%  beta2_mom , 
                                    epsilon = epsilon, 
                                    iter = 1, 
                                    momentum = control_list$beta_phi$momentum %||%  momentum )
  
  controls$a <- get_controls(optimizer_a,  
                             lr = control_list$a$lr %||%  (0.5/nr_time), 
                             len = nr_edges, 
                             beta1 = control_list$a$beta1 %||%  beta1_mom, 
                             beta2 = control_list$a$beta2 %||%  beta2_mom , 
                             epsilon = epsilon, 
                             iter = 1, 
                             momentum = control_list$a$momentum %||%  momentum )
  
  
  controls$pi <- get_controls(optimizer_pi, 
                              lr = control_list$pi$lr %||% ( 0.01/nrow(X)), 
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
  cat("\nStarting EM updates:\n")
  
  iteration_times <- c()
  for (iter in 1:n_iter) {
    
    a1 <- Sys.time()
    
    # 0. Prepare updates
    # Store old parameters and update
    pi_old <- pi_est
    beta_old <- beta_est
    beta_phi_old <- beta_phi_est
    phi <- exp(beta_phi_est)
    
    if(model_type == "learn_graph" & !a_known){
      a_old <- a_est
      A <- get_W_from_array(a_est, nr_regions)
      psi_est <- NA
      psi_old <- NA
    }else if(model_type == "learn_psi"){
      a_est <- 1
      a_old <- 1
      psi_old <- psi_est
    }else if(a_known){
      a_est <- A[upper.tri(A, TRUE)]
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
    
    # E-steps that can be calculated seperately:
    
    if(mixing == "gamma"){
      w1 <- expected_h(y = claims_batch,  phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density,
                       h_func = function(z, y, pi, mu, phi) identity(z), likelihood_func = zip_likelihood_func, method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
      w2 <- expected_h(y = claims_batch, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density,
                       h_func = function(z, y, pi, mu, phi) log(z), likelihood_func = zip_likelihood_func,
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    }else if(mixing == "ig"){
      w1 <- expected_h(y = claims_batch, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density,
                       h_func = function(z, y, pi, mu, phi) identity(z), likelihood_func = zip_likelihood_func,
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
      w3 <- expected_h(y = claims_batch, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density,
                       h_func = function(z, y, pi, mu, phi) 1/z, likelihood_func = zip_likelihood_func,
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    } else if(mixing == "ln"){
      w4 <- expected_h(y = claims_batch, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density,
                       h_func = function(z, y, pi, mu, phi) log(z)^2, likelihood_func = zip_likelihood_func,
                       method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
    }
    
    
    # Next we either optimize via L-BFGS-B 
    # or some gradient descent method
    
    if(do_optim){

      if(model_type == "learn_graph"){
        print("doing optim")
        lower <- c(-Inf, 1e-3, rep(-Inf,ncol(X)), rep(1e-8, nr_regions*(nr_regions+1)/2))
        upper <- c(Inf, 1-1e-3, rep(Inf,ncol(X)), rep(Inf, nr_regions*(nr_regions+1)/2))
        
        out_optim <- optim(par = c(beta_phi_est, pi_est, beta_est, a_est), fn = beta_optim, gr = beta_optim_grad, 
                           control = list(maxit = 2, trace = 3), method = "L-BFGS-B",
                           lower = lower,
                           upper = upper
        )
        par <- out_optim$par
        
        beta_phi_est <- par[1]
        pi_est <- par[2]
        beta_est <- par[3:(ncol(X)+2)]
        a_est <- par[(ncol(X)+3):length(par)]
        print(out_optim$message)
        A <- get_W_from_array(a_est, nr_regions)
        
      }else if (model_type == "learn_psi"){
        lower <- c(-Inf, 1e-3, rep(-Inf,ncol(X)), rep(1e-8, nr_regions))
        upper <- c(Inf, 1-1e-3, rep(Inf,ncol(X)), rep(Inf, nr_regions))
        
        out_optim <- optim(par = c(beta_phi_est, pi_est, beta_est, psi_est), fn = beta_optim, gr = beta_optim_grad, 
                           control = list(maxit = 2, trace = 3), method = "L-BFGS-B",
                           lower = lower,
                           upper = upper
        )
        par <- out_optim$par
        
        beta_phi_est <- par[1]
        pi_est <- par[2]
        beta_est <- par[3:(ncol(X)+2)]
        psi_est <- par[(ncol(X)+3):length(par)]
        print(out_optim$message)
      }
      
      
      beta_converged <-  isTRUE(all.equal(beta_est, beta_old, tolerance = param_tol)) 
      if(model_type == "learn_graph"){
        spatial_converged <- isTRUE(all.equal(a_est, a_old, tolerance = param_tol)) 
      }else if(model_type == "learn_psi"){
        spatial_converged <- isTRUE(all.equal(psi_est, psi_old, tolerance = param_tol))
      }
      
      beta_phi_converged <- isTRUE(all.equal(beta_phi_est, beta_phi_old, tolerance = param_tol))
      
    }else{
      # 1. E-steps
      # Compute expected derivative with respect to mu
      grad_mu <- expected_h(y = claims_batch, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density,
                            h_func = h_mu_func, likelihood_func = zip_likelihood_func,
                            method  = Emethod, S = control_list$S, n_nodes = control_list$n_nodes, rzdensity = rzdensity)
      
      grad_pi <- expected_pi_grad(y = claims_batch, phi = phi, pi = pi_est, mu = mu_batch, z_density = z_density)
      
      
      # Compute gradient for beta via chain rule: sum_i (grad_mu_i * nu_i * x_i)
      grad_beta_total <- get_grad_beta(additive, grad_mu, nu_batch, exposure_batch, X_batch, spatial_effect_batch, locs_batch)
      
      # Compute psi/a gradient
      if(model_type == "learn_graph" & !a_known){
        grad_a_total <- get_grad_a(additive, grad_mu, agg_claims, years_batch, locs_batch, exposure_batch, nu_batch, lambda, nr_regions)
      }else if (model_type == "learn_psi"){
        grad_psi_total <- get_grad_psi(additive, grad_mu, agg_effect_batch, nu_batch, exposure_batch, locs_batch, nr_regions)
      }
      
      # Compute gradient for pi (scalar)
      grad_pi_total <- get_grad_pi(grad_pi)
      
      # phi grad
      grad_beta_phi_total <- get_grad_beta_phi(phi,mixing, w1, w2, w3, w4)
      
      # scale gradient for sgd
      if(sgd){
        grad_beta_total <- grad_beta_total* (N / batch_size)
        if(model_type == "learn_graph" & !a_known){
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
      
      # # update psi/a
      if(model_type == "learn_graph" & !a_known){
        a_out <- update_gradient(optimizer_a, a_est, grad_a_total, controls$a)
        a_est <- a_out$param
        controls$a <- a_out$controls
        a_est <- pmax(a_est, 1e-6)
        spatial_converged <- isTRUE(all.equal(a_est, a_old, tolerance = param_tol)) 
        A <- get_W_from_array(a_est, nr_regions)
      }else if(model_type == "learn_psi"){
        psi_out <- update_gradient(optimizer_psi, psi_est, grad_psi_total, controls$psi)
        psi_est <- psi_out$param
        controls$psi <- psi_out$controls
        psi_est <- pmax(psi_est, 1e-6)
        spatial_converged <- isTRUE(all.equal(psi_est, psi_old, tolerance = param_tol))
      }else if(a_known){
        spatial_converged <- TRUE
      }
      
      
      # update pi
      pi_out <- update_gradient(optimizer_pi, pi_est, grad_pi_total, controls$pi)
      pi_est <- pi_out$param
      controls$pi <- pi_out$controls
      pi_est <- pmin(pmax(pi_est, 1e-6), 1-1e-6)

      
      # update phi
      beta_phi_out <- update_gradient(optimizer_beta_phi, beta_phi_est, grad_beta_phi_total, controls$beta_phi)
      beta_phi_est <- beta_phi_out$param
      controls$beta_phi <- beta_phi_out$controls
      
      
      beta_phi_converged <- isTRUE(all.equal(beta_phi_est, beta_phi_old, tolerance = param_tol))
      
    }
    
    # 3. Check Convergence
    if(beta_converged & spatial_converged & beta_phi_converged  ){
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
      log_lik <- sum(log(zzip_log_lik_function(claims,  exp(beta_phi_est), pi_est,  mu,  z_density))) 
      if(model_type == "learn_graph"){
        log_lik <- log_lik - lambda*sum(a_est)
      }
      
      if(iter > 2){
        if(isTRUE(all.equal(log_lik, log_lik_old, tolerance = Q_tol))  ){
          
          print("Breaking because log-likelihood has stopped changing")
          break
        } 
      }
      if(log_lik < log_lik_old){
        
        print("There is a decrease in the likelihood. Some numerical instabilities occuring. Change step sizes or SGD batch if using SGD.")
      }
      
    }
    
    if(verbose >= 1){
      print(paste0("log-likelihood value  ", sum(log_lik)))
    }
    
    if(verbose >=2){
      if(model_type == "learn_graph"){
        cat(sprintf("Iteration %d: beta = [%s], a = [%s], beta_phi = %f, pi = %f \n,",
                    iter,
                    paste(round(beta_est, 4), collapse = ", "),
                    paste(round(a_est, 4), collapse = ", "),
                    beta_phi_est,
                    pi_est))
      }else if (model_type == "learn_psi"){
        cat(sprintf("Iteration %d: beta = [%s], psi = [%s], beta_phi = %f, pi = %f \n,",
                    iter,
                    paste(round(beta_est, 4), collapse = ", "),
                    paste(round(psi_est, 4), collapse = ", "),
                    beta_phi_est,
                    pi_est))
      } 
    }
    
    
    
    
    iteration_times[iter] <-  Sys.time() - a1 
  }
  

  
  
  
  # Finalize 
  # Hessian
  
  if(calc_se){
    y <- claims
    
    if(model_type == "learn_graph"){
      Hessian <- optimHess(c(beta_phi_est, pi_est, beta_est, a_est), fn = beta_optim, gr = beta_optim_grad)
    }else if (model_type == "learn_psi"){
      Hessian <- optimHess(c(beta_phi_est, pi_est, beta_est, psi_est), fn = beta_optim, gr = beta_optim_grad)
    }
    
    # Variance
    
    grad_mu_sq <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                             h_func = function(z, y, pi, mu, phi) h_mu_func(z, y, pi, mu, phi)^2, likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0,
                             mu_old = mu, phi_old = phi, method = "integration")
    
    grad_phi_sq <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                              h_func = function(z, y, pi, mu, phi) deriv_log_density(z, phi)^2, likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0,
                              mu_old = mu, phi_old = phi, method = "integration")
    
    grad_pi_sq <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                             h_func = function(z, y, pi, mu, phi) h_pi_func(z, y, pi, mu, phi)^2, likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0,
                             mu_old = mu, phi_old = phi, method = "integration")
    
    
    
    grad_mu_phi <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                              h_func = function(z, y, pi, mu, phi) h_mu_func(z, y, pi, mu, phi)* deriv_log_density(z, phi), likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0,
                              mu_old = mu, phi_old = phi, method = "integration")
    
    grad_mu_pi <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                             h_func = function(z, y, pi, mu, phi) h_mu_func(z, y, pi, mu, phi)* h_pi_func(z, y, pi, mu, phi), likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0,
                             mu_old = mu, phi_old = phi, method = "integration")
    
    grad_phi_pi <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                              h_func = function(z, y, pi, mu, phi) deriv_log_density(z, phi)* h_pi_func(z, y, pi, mu, phi), likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0,
                              mu_old = mu, phi_old = phi, method = "integration")
    
    
    
    
    
    a_known <- FALSE
    if(a_known){
      W11 <- matrix(0, nrow = ncol(X)+2, ncol = ncol(X)+2)
    }else if(model_type == "learn_graph"){
      W11 <- matrix(0, nrow = ncol(X)+2 + length(a_est), ncol = ncol(X)+2+ length(a_est))
    }else if(model_type == "learn_psi"){
      W11 <- matrix(0, nrow = ncol(X)+2 + length(psi_est), ncol = ncol(X)+2+ length(psi_est))
    }else{
      W11 <- matrix(0, nrow = ncol(X)+2, ncol = ncol(X)+2)
    }
    
    for(i in 1:nrow(X)){
      
      g1 <- get_grad_beta(additive, 1, nu[i], exposure[i], X[i, , drop = FALSE], se$spatial_effect[i], locs[i])
      
      if(a_known){
        g22 <- c()
      }else if(model_type == "learn_graph"){
        g22 <- get_grad_a_2(additive, 1, agg_claims, years[i],locs[i], exposure[i], nu[i], lambda, nr_regions)
      }else if(model_type == "learn_psi"){
        g22 <- get_grad_psi_2(additive, 1, se$agg_effect[i], nu[i], exposure[i], locs[i], nr_regions)
      }else{
        g22 <- c()
      }
      
      
      
      u1 <- grad_mu_sq[i] * outer(c(g1, g22), c(g1, g22))
      u2 <- grad_phi_sq[i]
      u3 <- grad_mu_phi[i] * c(g1, g22)
      u4 <- grad_pi_sq[i]
      u5 <- grad_mu_pi[i] * c(g1, g22)
      u6 <- grad_phi_pi[i]
      
      
      
      W11[1,1] <-  W11[1,1] + u2
      
      W11[1,2] <- W11[1,2] + u6
      W11[2,1] <- W11[2,1] + u6
      
      W11[1,3:ncol(W11)] <- W11[1,3:ncol(W11)] + u3
      W11[3:ncol(W11), 1] <-  W11[3:ncol(W11), 1] + u3
      
      W11[2,2] <- W11[2,2] + u4
      
      W11[2,3:ncol(W11)] <-  W11[2,3:ncol(W11)] + u5
      W11[3:ncol(W11), 2] <-  W11[3:ncol(W11), 2] + u5
      
      
      W11[3:ncol(W11), 3:ncol(W11)] <- W11[3:ncol(W11), 3:ncol(W11)] + u1
      
    }
    
    
    grad_mu <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                          h_func = function(z, y, pi, mu, phi) h_mu_func(z, y, pi, mu, phi), likelihood_func = zip_likelihood_func,
                          special_case = FALSE, special_fun = 0, method = "integration")
    
    grad_phi <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                           h_func = function(z, y, pi, mu, phi) deriv_log_density(z, phi),
                           likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0, method = "integration")
    
    grad_pi <- expected_h(y = y, phi = phi, pi = pi_est, mu = mu, z_density = z_density, 
                          h_func = function(z, y, pi, mu, phi) h_pi_func(z, y, pi, mu, phi),
                          likelihood_func = zip_likelihood_func, special_case = FALSE, special_fun = 0, method = "integration")
    
    
    
    if(a_known){
      W22 <- matrix(0, nrow = ncol(X)+2, ncol = ncol(X)+2)
    }else if(model_type == "learn_graph"){
      W22 <- matrix(0, nrow = ncol(X)+2 + length(a_est), ncol = ncol(X)+2+ length(a_est))
    }else if(model_type == "learn_psi"){
      W22 <- matrix(0, nrow = ncol(X)+2 + length(psi_est), ncol = ncol(X)+2+ length(psi_est))
    }else{
      W22 <- matrix(0, nrow = ncol(X)+2, ncol = ncol(X)+2)
    }
    
    for(i in 1:nrow(X)){
      
      g1 <- get_grad_beta(additive, 1, nu[i], exposure[i], X[i, , drop = FALSE], se$spatial_effect[i], locs[i])
      
      if(a_known){
        g22 <- c()
      }else if(model_type == "learn_graph"){
        g22 <- get_grad_a_2(additive, 1, agg_claims, years[i],locs[i], exposure[i], nu[i], 0, nr_regions)
      }else if(model_type == "learn_psi"){
        g22 <- get_grad_psi_2(additive, 1, se$agg_effect[i], nu[i], exposure[i], locs[i], nr_regions)
      }else{
        g22 <- c()
      }
      
      u1 <- grad_mu[i]^2 * outer(c(g1, g22), c(g1, g22))
      u2 <- (grad_phi[i])^2
      u3 <- grad_mu[i]*grad_phi[i] * c(g1, g22)
      
      u4 <- grad_pi[i]^2
      u5 <- grad_mu[i]*grad_pi[i] * c(g1, g22)
      u6 <- grad_pi[i] * grad_phi[i]
      
      W22[1,1] <-  W22[1,1] + u2
      
      W22[1,2] <-  W22[1,2] + u6
      W22[2,1] <-  W22[2,1] + u6
      
      W22[1,3:ncol(W22)] <- W22[1,3:ncol(W22)] + u3
      W22[3:ncol(W22), 1] <-  W22[3:ncol(W22), 1] + u3
      
      W22[2,2] <- W22[2,2] + u4
      
      W22[2,3:ncol(W22)] <-  W22[2,3:ncol(W22)] + u5
      W22[3:ncol(W22), 2] <-  W22[3:ncol(W22), 2] + u5
      
      W22[3:ncol(W22), 3:ncol(W22)] <- W22[3:ncol(W22), 3:ncol(W22)] + u1
      
    }
    
    
    
    var_loglik <- W11-W22
  }else{
    Hessian <- NA
    var_loglik <- NA
  }
  
  # Find number of params
  if(a_known){
    nr_param <- length(beta_est) + 1
  }else if(model_type == "learn_graph"){
    nr_param <- length(beta_est) + 1 + sum(abs(a_est) > 1e-3)
  }else if(model_type == "learn_psi"){
    nr_param <- length(beta_est) + 1 + length(psi_est)
  }
  
  
  
  return(list(beta1 = beta_est, pi = pi_est, a = a_est, psi = psi_est, beta2 = beta_phi_est, mu = mu, Hessian = Hessian, 
              log_lik = log_lik, var_loglik = var_loglik, nr_param = nr_param, time_to_init = time_to_init, iteration_times = iteration_times ))
  
}


#########################
# 4. Test
#########################
# n <- 500
# t <- 2
# data_sim <- simulate_claims(n, t, "psi", TRUE, mixing = "gamma", model_type = "zip",
#                             exposure_lambda = 0, area = 20, seed = 2, order_regions = TRUE, density = 0.2)
#
# # Extract variables from simulation
# claims <- data_sim$claims
# X <- data_sim$X
# years <- data_sim$years
# locs <- data_sim$locs
# agg_claims <- data_sim$agg_claims
# A <- data_sim$A
# exposure <- data_sim$exposure
# model_type <- "learn_psi"
# additive <- TRUE
# mixing <- "gamma"

#
# out_zip <- zip(claims, X, locs, years,  agg_claims,
#                A, additive, model_type, lambda = 0, exposure = exposure, max_itr = 1, calc_hessian = FALSE)
# #
# out <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, model_type, additive, mixing,  Emethod = "integration",
#               n_iter = 10, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
#               optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
#               batch_size = 100, verbose = 2, do_optim = FALSE, calc_se = FALSE, max_itr_poisson = 1)
# sum(abs(out_zip$beta1 - data_sim$beta1))
# sum(abs(out$beta1 - data_sim$beta1))


# diag(solve(out$Hessian - out$var_loglik))
# 
# diag(solve(out$Hessian))
