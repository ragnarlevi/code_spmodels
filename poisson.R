#### This file contains Possion models. That is - no mixing and no EM needed.

set.seed(123)
source("utils.R")
source("simulate_data.R")


##########
# Define a Poisson model with no mixing - base
######
# Objective of Poisson

Q_P_no_phi <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, eother, scaling_factor = 1, a_known = FALSE){
  
  etheta <- as.numeric(etheta)

    if(model_type == "ordinary"){
      beta1 <- param[1:(ncol(X_mat)+0)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[1:(ncol(X_mat)+0)]
      psi <- param[(ncol(X_mat)+1):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[1:(ncol(X_mat)+0)]
      a <- param[(ncol(X_mat)+1):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }
    
    if(a_known){
      beta1 <- param
      a <- 0
    }
    
    
    

  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  ll <- sum(claims*log(mu) -mu*etheta) - sum(a*lambda)*(model_type == "learn_graph")  # note minus in z
  

  return(-ll)
  
}


# deriviative of the Poisson
Q_P_deriv_no_phi <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, eother, scaling_factor = 1, a_known = FALSE){
  
  
  etheta <- as.numeric(etheta)
  
    if(model_type == "ordinary"){
      beta1 <- param[1:(ncol(X_mat)+0)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[1:(ncol(X_mat)+0)]
      psi <- param[(ncol(X_mat)+1):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[1:(ncol(X_mat)+0)]
      a <- param[(ncol(X_mat)+1):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }
    
    if(a_known){
      beta1 <- param
      a <- 0
    }
    

  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- exp(X_mat %*% beta1)
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  
  der1_logdp <- claims/mu - etheta  # same for all models
  
  
  
  
  # add gradients depending on model type
  if(model_type == "ordinary"){
    grad <- t( der1_logdp * nu * exp(log(exposure)) ) %*% X_mat
    
  }else if(model_type == "learn_psi"){
    
    if(additive){
      #beta
      g1 <- t( der1_logdp * nu * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- der1_logdp[,1]*exp(log(exposure)) * se$agg_effect
      g2 <- tapply(g2,locs,sum)
    }else{
      # beta
      g1 <- t( der1_logdp * nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- der1_logdp[,1]*exp(log(exposure)) * se$agg_effect * nu
      g2 <- tapply(g2,locs,sum)
    }
    
    grad <- c(g1, g2)
    
    
  }else if(model_type == "learn_graph"){
    if(additive){
      g1 <- t( der1_logdp * nu * exp(log(exposure)) ) %*% X_mat
      
      # a
      g2 <- der1_logdp[,1]*t(agg_claims[, years]) * exp(log(exposure))
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    }else{
      # beta
      g1 <- t( der1_logdp * nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      
      # A
      g2 <-  der1_logdp[,1]*t(agg_claims[, years]) * exp(log(exposure)) *as.numeric(nu)
      
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
      
    }
    
    if(a_known){
      g22 <- c()
    }
    
    grad <- c(g1, g22)*scaling_factor
  }
  
  grad <- -grad

  return(grad)
  
}



Q_P <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, eother, Q_beta_2, Q_beta_2_deriv, mixing, scaling_factor = 1, a_known = FALSE){
  
  etheta <- as.numeric(etheta)
  if(mixing){
    beta2 <- param[1]
    if(model_type == "ordinary"){
      beta1 <- param[2:(ncol(X_mat)+1)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[2:(ncol(X_mat)+1)]
      psi <- param[(ncol(X_mat)+2):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[2:(ncol(X_mat)+1)]
      a <- param[(ncol(X_mat)+2):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }    
    
    if(a_known){
      beta1 <- param[2:length(param)]
      a <- 0
    }
    
  }else{
    if(model_type == "ordinary"){
      beta1 <- param[1:(ncol(X_mat)+0)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[1:(ncol(X_mat)+0)]
      psi <- param[(ncol(X_mat)+1):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[1:(ncol(X_mat)+0)]
      a <- param[(ncol(X_mat)+1):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }
    
    if(a_known){
      beta1 <- param
      a <- 0
    }
    
    
    
  }

  

  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  ll <- sum(claims*log(mu) -mu*etheta) - sum(a*lambda)*(model_type == "learn_graph")  # note minus in z
  
  if(mixing){
    ll <- ll - Q_beta_2(beta2,  claims, exposure, etheta, eother)
  }
  
  return(-ll)
  
}


# deriviative of the Poisson
Q_P_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, eother, Q_beta_2, Q_beta_2_deriv, mixing, scaling_factor = 1, a_known = FALSE){
  

  etheta <- as.numeric(etheta)
  
  if(mixing){
    beta2 <- param[1]
    if(model_type == "ordinary"){
      beta1 <- param[2:(ncol(X_mat)+1)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[2:(ncol(X_mat)+1)]
      psi <- param[(ncol(X_mat)+2):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[2:(ncol(X_mat)+1)]
      a <- param[(ncol(X_mat)+2):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }    
    
    
    if(a_known){
      beta1 <- param[2:length(param)]
      a <- 0
    }
    
    
  }else{
    if(model_type == "ordinary"){
      beta1 <- param[1:(ncol(X_mat)+0)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[1:(ncol(X_mat)+0)]
      psi <- param[(ncol(X_mat)+1):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[1:(ncol(X_mat)+0)]
      a <- param[(ncol(X_mat)+1):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }
    
    if(a_known){
      beta1 <- param
      a <- 0
    }
    
  }
  

  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- exp(X_mat %*% beta1)
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  

  der1_logdp <- claims/mu - etheta  # same for all models
  
  
  
  
  # add gradients depending on model type
  if(model_type == "ordinary"){
    grad <- t( der1_logdp * nu * exp(log(exposure)) ) %*% X_mat
    
  }else if(model_type == "learn_psi"){
    
    if(additive){
      #beta
      g1 <- t( der1_logdp * nu * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- der1_logdp[,1]*exp(log(exposure)) * se$agg_effect
      g2 <- tapply(g2,locs,sum)
    }else{
      # beta
      g1 <- t( der1_logdp * nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- der1_logdp[,1]*exp(log(exposure)) * se$agg_effect * nu
      g2 <- tapply(g2,locs,sum)
    }
    
    grad <- c(g1, g2)
    
    
  }else if(model_type == "learn_graph"){
    if(additive){
      g1 <- t( der1_logdp * nu * exp(log(exposure)) ) %*% X_mat
      
      # a
      g2 <- der1_logdp[,1]*t(agg_claims[, years]) * exp(log(exposure))
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    }else{
      # beta
      g1 <- t( der1_logdp * nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      
      # A
      g2 <-  der1_logdp[,1]*t(agg_claims[, years]) * exp(log(exposure)) *as.numeric(nu)
      
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
      
    }
    
    if(a_known){
      g22 <- c()
    }
    
    grad <- c(g1, g22)*scaling_factor
  }
  
  grad <- -grad
  if(mixing){
    beta2_deriv <- Q_beta_2_deriv(beta2,  claims, exposure, etheta, eother)
    grad <- c(beta2_deriv, grad)
  }
  
  
  return(grad)
  
}

Q_P_deriv2 <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, scaling_factor = 1, a_known = FALSE){
  
  
  etheta <- as.numeric(etheta)
  
    beta2 <- param[1]
    if(model_type == "ordinary"){
      beta1 <- param[2:(ncol(X_mat)+1)]
    }else if(model_type == "learn_psi"){
      beta1 <- param[2:(ncol(X_mat)+1)]
      psi <- param[(ncol(X_mat)+2):length(param)]
      a <- 1
    }else if(model_type == "learn_graph" & !(a_known)){
      beta1 <- param[2:(ncol(X_mat)+1)]
      a <- param[(ncol(X_mat)+2):length(param)]
      A <- get_W_from_array(a, p)
      psi <- NA
    }    
  
  
  if(a_known){
    beta1 <- param
    a <- 0
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- exp(X_mat %*% beta1)
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  

  
  # add gradients depending on model type
  if(model_type == "ordinary"){
    grad <- t( nu * exp(log(exposure)) ) %*% X_mat
    
  }else if(model_type == "learn_psi"){
    
    if(additive){
      #beta
      g1 <- t(nu * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- exp(log(exposure)) * se$agg_effect
      g2 <- tapply(g2,locs,sum)
    }else{
      # beta
      g1 <- t(nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- exp(log(exposure)) * se$agg_effect * nu
      g2 <- tapply(g2,locs,sum)
    }
    
    grad <- c(g1, g2)
    
    
  }else if(model_type == "learn_graph"){
    if(additive){
      g1 <- t( nu * exp(log(exposure)) ) %*% X_mat
      
      # a
      g2 <- t(agg_claims[, years]) * exp(log(exposure))
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    }else{
      # beta
      g1 <- t(  nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      
      # A
      g2 <-  t(agg_claims[, years]) * exp(log(exposure)) *as.numeric(nu)
      
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) <- 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
      
    }
    
    if(a_known){
      g22 <- c()
    }
    
    grad <- c(g1, g22)*scaling_factor
  }
  
  
  
  return(-c(-1, grad))
  
}


# Poisson model
Poisson <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure,  lambda = 0,  max_itr = 1000, a_known = FALSE){
  
  
  if(model_type != "learn_graph"){
    a_known <- FALSE
  }
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  

  # add initial param depending on model type
  if(model_type == "ordinary"){
    theta0 <- c(beta1)
    lower <- c(rep(-20,ncol(X1)))
  }else if(model_type == "learn_psi"){
    psi <- rep(1,p)
    theta0 <- c(beta1, psi)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
  }else if(model_type == "learn_graph" & !a_known){
    a <- rep(1,p*(p+1)/2)*0.001
    A <- get_W_from_array(a, p)
    theta0 <- c(beta1, a)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2))
    psi <- NA
  }else if(a_known){
    theta0 <- c(beta1)
    lower <- c(rep(-20,ncol(X1)))
    psi <- NA
  }else{
    stop("model type not known")
  }
  
  if(a_known){
    theta0 <- c(beta1)
    a <- A[upper.tri(A, TRUE)]
  }
  

    
    
    se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
    nu <- as.numeric(exp(X1 %*% beta1))
    mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
    
    
    etheta <- 1
    
    ### poisson part
    out <- optim(par = theta0,
                 fn = Q_P,
                 gr = Q_P_deriv,
                 locs = locs,
                 claims = claims, 
                 exposure = exposure, 
                 X_mat = X1, 
                 agg_claims = agg_claims, 
                 years = years,
                 etheta = etheta,
                 A = A,
                 additive = additive,
                 model_type = model_type,
                 lambda = lambda,
                 p = p,
                 a_known = a_known,
                 eother = 0, 
                 Q_beta_2 = NA, 
                 Q_beta_2_deriv = NA,
                 mixing = FALSE,
                 method = 'L-BFGS-B',
                 control = list(maxit = max_itr),
                 lower = lower,
                 hessian = T)
    

  
  
  if(model_type == "ordinary"){
    beta1 <- out$par[1:ncol(X1)]
  }else if(model_type == "learn_psi"){
    beta1 <- out$par[1:ncol(X1)]
    psi <- out$par[(ncol(X1)+1):length(out$par)]
  }else if(model_type == "learn_graph" & !a_known){
    beta1 <- out$par[1:ncol(X1)]
    a <- out$par[(ncol(X1)+1):length(out$par)]
    A <- get_W_from_array(a, p)
  }
    
  if(a_known){
    beta1 <- out$par
  }
  
    
    # Find number of params
    if(a_known){
      nr_param <- length(beta1) + 1
    }else if(model_type == "learn_graph"){
      nr_param <- length(beta1) + 1 + sum(abs(A[upper.tri(A, diag = T)]) > 1e-3)
    }else if(model_type == "learn_psi"){
      nr_param <- length(beta1) + 1 + length(psi)
    }
    
  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], beta2 = 1 , H = out$hessian, H_beta2 = NA, mu = mu, 
              optim_obj = out, model_type = model_type, nr_param = nr_param))
  
  
}



##################################
# Define "ordinary Poisson" 
# with mixing
# We use L-BFGS-B to optimize this model
##################################



############## Poisson Gamma #############

# E-step for Poisson-Gamma
Etheta <- function(claims,lambda,phi) (phi + claims)/(phi + lambda)
Elogtheta <- function(claims,lambda,phi) digamma(phi + claims) - log(phi + lambda)

# Q function for the mixing Gamma
Q_PG_beta_2 <- function(param,  claims, exposure, etheta, eother){
  
  etheta = as.numeric(etheta)
  elogtheta = as.numeric(eother)
  
  nu2 <- exp(param)
  phi <- nu2 
  
  
  ll = sum(phi*log(phi) -lgamma(phi) + (phi-1)*elogtheta - phi*etheta )
  return(-ll)
  
}

Q_PG_beta_2_deriv <- function(param,  claims, exposure, etheta, eother){
  
  etheta = as.numeric(etheta)
  elogtheta = as.numeric(eother)
  
  nu2 <- exp(param)
  phi <- nu2 
  
  der1_logdga <- log(phi) + 1 - digamma(phi) + elogtheta - etheta
  
  
  g3 <- sum( der1_logdga*nu2 ) 
  
  return(-g3)
}



################ Poisson IG ############


# E-step for PIG
Etheta_gig <- function(claims, a, b){
  
  
  num <- (sqrt(b)/sqrt(a)) * besselK(sqrt(a*b), claims+1-0.5,expon.scaled = TRUE)
  denom <- besselK(sqrt(a*b), claims-0.5,expon.scaled = TRUE)
  
  return(num/denom)
}
# E-step for PIG
Einvtheta_gig <- function(claims, a, b){
  num <- (sqrt(a)/sqrt(b)) * besselK(sqrt(a*b), claims+1-0.5,expon.scaled = TRUE)
  denom <- besselK(sqrt(a*b), claims-0.5,expon.scaled = TRUE)
  
  
  return(num/denom- 2*(claims-0.5)/b)
  
}

# Q function for the mixing IG
Q_PIG_beta_2 <- function(param, claims, exposure, etheta, eother){
  
  etheta = as.numeric(etheta)
  einvtheta = as.numeric(eother)
  
  nu2 <- exp(param)
  phi <- nu2 
  
  
  ll <- sum(log(phi) + phi^2 - 0.5*phi^2*(einvtheta + etheta))
  
  return(-ll)
  
}

Q_PIG_beta_2_deriv <- function(param, claims, exposure, etheta, eother){
  
  etheta = as.numeric(etheta)
  einvtheta = as.numeric(eother)
  
  nu2 <- exp(param)
  phi <- nu2 
  
  der1_logdgauss <- 1/phi + 2*phi - phi*einvtheta - phi*etheta
  g3 <- sum(der1_logdgauss* nu2)
  
  return(-g3)
}







############# Poisson LN ##########




# Q function for the mixing log normal
Q_PLN_beta_2 <- function(param, claims, exposure, etheta = 1, eother){
  
  elogthetasq = as.numeric(eother)
  
  nu2 <- exp(param)
  phi <- nu2
  
  
  ll = sum(-log(phi) - elogthetasq/(2*phi^2)  - phi^2/8)
  return(-ll)
  
}

Q_PLN_beta_2_deriv <- function(param, claims, exposure, etheta = 1, eother){
  
  elogthetasq = as.numeric(eother)
  
  nu2 <- exp(param)
  phi <- nu2
  
  der1_logdig <- -1/phi + elogthetasq/(phi ^3) - 2*phi/8
  
  
  g3 <- sum( der1_logdig*nu2 )
  
  return(-g3)
}



 pois_likelihood_func <- function(z, y, pi, mu, phi){
   dpois(y, lambda = z*mu)
 }

 
 
 library(statmod)
 
 zpois_log_lik_function <- function(y, phi, pi, mu, z_density) {
   # f(y | z)
   f_y_given_z <- function(z) {
     pois_likelihood_func(z, y, pi, mu, phi)
   }
   
   # Unnormalized numerator integrand: f(y|z) * f(z)
   numerator_integrand <- function(z) {
     val <- f_y_given_z(z) * z_density(z, phi)
     val[!is.finite(val)] <- 0
     val
   }
   
   # 1) Try integrate()
   num_result <- tryCatch({
     integrate(
       numerator_integrand,
       lower        = 0,
       upper        = Inf,
       subdivisions = 1000,
       rel.tol      = 1e-8,
       stop.on.error= FALSE
     )$value
   }, error = function(e) {
     NA_real_
   }, warning = function(w) {
     NA_real_
   })
   
   # 2) If integrate() failed or gave NA/NaN/Inf, fall back to Gauss–Legendre
   if (!is.finite(num_result)) {
     quad      <- statmod::gauss.quad(50, kind = "legendre")
     t_nodes   <- (quad$nodes + 1) / 2
     t_weights <- quad$weights / 2
     
     # Transform [0,1] → [0,∞): z = t/(1-t), dz = 1/(1-t)^2 dt
     z_vals    <- t_nodes / (1 - t_nodes)
     jacobian  <- 1 / (1 - t_nodes)^2
     
     lik_vals   <- f_y_given_z(z_vals)
     prior_vals <- z_density(z_vals, phi)
     
     num_result <- sum(t_weights * lik_vals * prior_vals * jacobian)
   }
   
   num_result
 }
 
 # vectorize over y and mu as before
 zpois_log_lik_function <- Vectorize(
   zpois_log_lik_function,
   vectorize.args = c("y", "mu")
 )
 
 
 

# Poisson mixture model 
 Poisson_mixed <- function(claims, X, locs, years, agg_claims, A, additive, model_type, exposure, 
                          lambda = 0, param_tol = 1e-5, Q_tol = 1e-5, nr_em = 100, max_itr = 1000, 
                          mixing_var = "gamma", z= "", verbose = 0, sgd = FALSE, batch_size = 500,
                          a_known = FALSE, control_list = list(), beta2_start = NULL){
  

   # claims <- Y_claims_train
   # X <- X_mat_train
   # locs <- locs_train
   # years <- years_train
   # agg_claims <- agg_claims
   # A <- A
   # additive <- TRUE
   # model_type <- "learn_graph"
   # exposure <- offset_train
   # lambda = 0
   # param_tol = 1e-5
   # Q_tol = 1e-5
   # nr_em = 100
   # max_itr = 100
   # mixing_var = "ig"
   # z = ""
   # verbose = 2
   # sgd = FALSE
   # control_list = list()
   # a_known = FALSE
   # beta2_start = NULL


  
  if(model_type != "learn_graph"){
    a_known <- FALSE
  }
  
  N <- length(claims)
  p <- length(unique(locs))
  
  # Set initial parameters
  out_poisson <- Poisson(claims, X, locs, years, agg_claims, A, additive, model_type, lambda = lambda, exposure, 
                         max_itr = max_itr, a_known = a_known)
  beta1 <- out_poisson$beta1
  if(verbose>0){
    print("initial parameters:")
    print(beta1)
  }
  
  
  # Compute the moment‐ratios:
  res <- claims - out_poisson$mu
  M2  <- mean((res^2 - out_poisson$mu) / out_poisson$mu^2)   # for gamma & lognormal
  M3  <- mean((res^2 - out_poisson$mu) / out_poisson$mu^3)   # for inverse Gaussian
  

  if(model_type == "learn_graph"){
  a <- out_poisson$a
  A <- get_W_from_array(a, p)
  }else if(model_type == "learn_psi"){
    psi <- out_poisson$psi
  }

  # add initial param depending on model type
  if(model_type == "ordinary"){
    theta0 <- c(beta1)
    lower <- c(rep(-20,ncol(X)))
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi)
    lower <- c(rep(-20,ncol(X)), rep(1e-8, p))
  }else if(model_type == "learn_graph"){
    theta0 <- c(beta1, a)
    lower <- c(rep(-20,ncol(X)), rep(1e-8, p*(p+1)/2))
    psi <- NA
  }else if(a_known){
    theta0 <- c(beta1)
    lower <- c(rep(-20,ncol(X)))
    psi <- NA
  }else{
    stop("model type not known")
  }
  
  # set initial for the mixing
  if(mixing_var == "gamma"){
    beta20 <- beta2_start%||%  1
    z_density <- gamma_density
    rzdensity <- rgamma_density
    deriv_log_density <- deriv_log_gamma_density
    Q_beta_2 <- Q_PG_beta_2
    Q_beta_2_deriv <- Q_PG_beta_2_deriv
  }else if(mixing_var == "ig"){
    beta20 <- beta2_start %||%  0
    z_density <- ig_density
    rzdensity <- rig_density
    deriv_log_density <- deriv_log_ig_density
    Q_beta_2 <- Q_PIG_beta_2
    Q_beta_2_deriv <- Q_PIG_beta_2_deriv
  }else if(mixing_var == "ln"){
    beta20 <- beta2_start %||%  0
    phi <- exp(beta20)
    z_density <- ln_density
    rzdensity <- rln_density
    deriv_log_density <- deriv_log_ln_density
    Q_beta_2 <- Q_PLN_beta_2
    Q_beta_2_deriv <- Q_PLN_beta_2_deriv
  }else{
    stop("mixing type not known")
  }
  

  if(a_known){
    theta0 <- c(beta1)
    a <- A[upper.tri(A, diag = TRUE)]
  }
  
  ## Loop starts ##
  
  log_lik <- -Inf

  for(iter in 1:nr_em){
    # 0. Prepare updates

    se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
    nu <- as.numeric(exp(X %*% beta1))
    mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
    phi <- exp(beta20)
    
    # If using SGD, sample a mini-batch
    if (sgd) {
      idx <- c()
      for(i in 1:p){
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
      scaling_factor <- (N / batch_size)  
    }else{
      claims_batch   <- claims
      X_batch        <- X
      exposure_batch <- exposure
      locs_batch     <- locs
      years_batch <- years
      mu_batch <- mu
      nu_batch <- nu
      spatial_effect_batch <- se$spatial_effect
      agg_effect_batch <- se$agg_effect
      scaling_factor <- 1
    }

    
    
    
    
    if(mixing_var == "gamma"){
      
        etheta <- Etheta(claims_batch, mu_batch, phi)
        eother <- Elogtheta(claims_batch, mu_batch, phi)

    }else if(mixing_var == "ig"){
      etheta <- Etheta_gig(claims_batch, phi^2 + 2*mu_batch, phi^2)
      eother <- Einvtheta_gig(claims_batch, phi^2 + 2*mu_batch, phi^2)
      
    }else if(mixing_var == "ln"){
      
      eother <- expected_h(y = claims_batch, phi = phi, pi = NA, mu = mu_batch, z_density = z_density, 
                       h_func = function(z, y, pi, mu, phi) log(z)^2, likelihood_func = pois_likelihood_func, 
                       method  = "GaussLegendre")
      
      etheta <- expected_h(y = claims_batch, phi = phi, pi = NA, mu = mu_batch, z_density = z_density, 
                          h_func = function(z, y, pi, mu, phi) z, likelihood_func = pois_likelihood_func, 
                          method  = "GaussLegendre")

    }

    
    ### poisson part
    out <- optim(par = theta0,
                 fn = Q_P_no_phi,
                 gr = Q_P_deriv_no_phi,
                 locs = locs_batch,
                 claims = claims_batch, 
                 exposure = exposure_batch, 
                 X_mat = X_batch, 
                 agg_claims = agg_claims, 
                 years = years_batch,
                 etheta = etheta,
                 A = A,
                 additive = additive,
                 model_type = model_type,
                 lambda = lambda,
                 p = p,
                 eother = eother, 
                 scaling_factor = scaling_factor,
                 a_known = a_known,
                 method = 'L-BFGS-B',
                 control = list(maxit = 2),
                 lower = lower)

    ## mixing part

  
    out_beta2 = optim(beta20,
                       fn = Q_beta_2,
                       gr = Q_beta_2_deriv,
                       claims = claims_batch,
                       exposure =exposure_batch,
                       etheta = etheta,
                       eother = eother,
                       method = 'L-BFGS',
                       control = list(maxit = 2))

    beta2 <- out_beta2$par
    H_beta2 <- out_beta2$hessian
    

    theta <- out$par

    if(model_type == "ordinary"){
      beta1 <- out$par[1:(ncol(X)+0)]
    }else if(model_type == "learn_psi"){
      beta1 <- out$par[1:(ncol(X)+0)]
      psi <- out$par[(ncol(X)+1):length(out$par)]
    }else if(model_type == "learn_graph" & !(a_known) ){
      beta1 <- out$par[1:(ncol(X)+0)]
      a <- out$par[(ncol(X)+1):length(out$par)]
      A <- get_W_from_array(a, p)
    }
    if(a_known){
      beta1 <- out$par[1:length(out$par)]
      theta <- out$par
    }

   
    mu_old <- mu
    # stop if expected likelihood has stopped decreasing
    se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
    nu <- as.numeric(exp(X %*% beta1))  # baseline mu from covariates
    
    # Get full mu by combining covariate effects (nu), spatial effects, and exposure
    mu <- get_mu(nu, se$spatial_effect, exposure, additive)
    
    
    
    # check convergence
    if(isTRUE(all.equal(theta, theta0, tolerance = param_tol)) & isTRUE(all.equal(beta2, beta20, tolerance = param_tol)) ){
      print("Breaking because parameters have stopped changing")
      break
    }else{
      theta0 <- theta
      beta20 <- beta2
    }
    
    
    if(iter >1){
      log_lik_old <- log_lik
      log_lik <- sum(log(zpois_log_lik_function(claims,  exp(beta2), 0, mu,  z_density)))
      if(model_type == "learn_graph"){
        log_lik <- log_lik - lambda*sum(a)
      }
      
      if(iter > 2){
        if(isTRUE(all.equal(log_lik, log_lik_old, tolerance = Q_tol))){
          print(paste0("New: ",log_lik, " vs Old:", log_lik_old))
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
        cat(sprintf("Iteration %d: beta = [%s], a = [%s], beta_phi = %f \n,",
                    iter,
                    paste(round(beta1, 4), collapse = ", "),
                    paste(round(a, 4), collapse = ", "),
                    beta2))
      }else if (model_type == "learn_psi"){
        cat(sprintf("Iteration %d: beta = [%s], psi = [%s], beta_phi = %f \n,",
                    iter,
                    paste(round(beta1, 4), collapse = ", "),
                    paste(round(psi, 4), collapse = ", "),
                    beta2))
      } 
    }
  
    
    
  }
  
  # Finalize 
  # Find Hessian
  y <- claims
  nr_regions <- p

  # Find variance
  if(mixing_var == "gamma"){
    etheta <- Etheta(claims_batch, mu, phi)
    eother <- Elogtheta(claims_batch, mu, phi)
  }else if(mixing_var == "ig"){
    etheta <- Etheta_gig(claims_batch, phi^2 + 2*mu, phi^2)
    eother <- Einvtheta_gig(claims_batch, phi^2 + 2*mu, phi^2)
  }else if(mixing_var == "ln"){
    
    eother <- expected_h(y = claims_batch, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                         h_func = function(z, y, pi, mu, phi) log(z)^2, likelihood_func = pois_likelihood_func, 
                         method  = "GaussLegendre")
    
    etheta <- expected_h(y = claims_batch, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                         h_func = function(z, y, pi, mu, phi) z, likelihood_func = pois_likelihood_func, 
                         method  = "GaussLegendre")
  }
  
  
  
  # Hessian part
  Hessian <- optimHess(par = c(beta2, theta),
                       fn = Q_P,
                       gr = Q_P_deriv,
                       locs = locs,
                       claims = claims, 
                       exposure = exposure, 
                       X_mat = X_batch, 
                       agg_claims = agg_claims, 
                       years = years_batch,
                       etheta = etheta,
                       A = A,
                       additive = additive,
                       model_type = model_type,
                       lambda = 0,
                       p = p,
                       eother = eother, 
                       Q_beta_2 = Q_beta_2, 
                       Q_beta_2_deriv = Q_beta_2_deriv,
                       mixing = TRUE,
                       scaling_factor = scaling_factor,
                       a_known = a_known)

  
  # Variance

  grad_mu_sq <- expected_h(y = y, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                                       h_func = function(z, y, pi, mu, phi) (y/mu-z)^2, likelihood_func = pois_likelihood_func, special_case = FALSE, special_fun = 0,
                                       mu_old = mu, phi_old = phi, method = "integration")
  
  grad_phi_sq <- expected_h(y = y, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                                        h_func = function(z, y, pi, mu, phi) deriv_log_density(z, phi)^2, likelihood_func = pois_likelihood_func, special_case = FALSE, special_fun = 0,
                                        mu_old = mu, phi_old = phi, method = "integration")
  
  
  grad_mu_sq <- expected_h(y = y, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                                       h_func = function(z, y, pi, mu, phi) (y/mu-z)^2, likelihood_func = pois_likelihood_func, special_case = FALSE, special_fun = 0,
                                       mu_old = mu, phi_old = phi, method = "integration")
  
  grad_mu_phi <- expected_h(y = y, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                                        h_func = function(z, y, pi, mu, phi) (y/mu-z)* deriv_log_density(z, phi), likelihood_func = pois_likelihood_func, special_case = FALSE, special_fun = 0,
                                        mu_old = mu, phi_old = phi, method = "integration")
  
  
  
  
  
  if(a_known){
    W11 <- matrix(0, nrow = ncol(X)+1, ncol = ncol(X)+1)
  }else if(model_type == "learn_graph"){
    W11 <- matrix(0, nrow = ncol(X)+1 + length(a), ncol = ncol(X)+1+ length(a))
  }else if(model_type == "learn_psi"){
    W11 <- matrix(0, nrow = ncol(X)+1 + length(psi), ncol = ncol(X)+1+ length(psi))
  }else{
    W11 <- matrix(0, nrow = ncol(X)+1, ncol = ncol(X)+1)
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
    
    W11[1,1] <-  W11[1,1] + u2
    W11[1,2:ncol(W11)] <- W11[1,2:ncol(W11)] + u3
    W11[2:ncol(W11), 1] <-  W11[2:ncol(W11), 1] + u3
    W11[2:ncol(W11), 2:ncol(W11)] <- W11[2:ncol(W11), 2:ncol(W11)] + u1
    
  }
  
  
  grad_mu <- expected_h(y = y, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                                    h_func = function(z, y, pi, mu, phi) (y/mu-z), likelihood_func = pois_likelihood_func,
                                    special_case = FALSE, special_fun = 0, method = "integration")
  
  grad_phi <- expected_h(y = y, phi = phi, pi = NA, mu = mu, z_density = z_density, 
                                     h_func = function(z, y, pi, mu, phi) deriv_log_density(z, phi),
                                     likelihood_func = pois_likelihood_func, special_case = FALSE, special_fun = 0, method = "integration")
  
  
  if(a_known){
    W22 <- matrix(0, nrow = ncol(X)+1, ncol = ncol(X)+1)
  }else if(model_type == "learn_graph"){
    W22 <- matrix(0, nrow = ncol(X)+1 + length(a), ncol = ncol(X)+1+ length(a))
  }else if(model_type == "learn_psi"){
    W22 <- matrix(0, nrow = ncol(X)+1 + length(psi), ncol = ncol(X)+1+ length(psi))
  }else{
    W22 <- matrix(0, nrow = ncol(X)+1, ncol = ncol(X)+1)
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
    
    W22[1,1] <-  W22[1,1] + u2
    W22[1,2:ncol(W22)] <- W22[1,2:ncol(W22)] + u3
    W22[2:ncol(W22), 1] <-  W22[2:ncol(W22), 1] + u3
    W22[2:ncol(W22), 2:ncol(W22)] <- W22[2:ncol(W22), 2:ncol(W22)] + u1
    
  }
  

  
  var_loglik <- W11-W22
  
  
  # Find number of params
  if(a_known){
    nr_param <- length(beta1) + 1
  }else if(model_type == "learn_graph"){
    nr_param <- length(beta1) + 1 + sum(abs(A[upper.tri(A, diag = T)]) > 1e-3)
  }else if(model_type == "learn_psi"){
    nr_param <- length(beta1) + 1 + length(psi)
  }
  


  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], beta2 = beta2 , Hessian = Hessian, mu = mu, 
              optim_obj = out, model_type = model_type, log_lik = log_lik, var_loglik = var_loglik, nr_param = nr_param, optim_obj_beta2 = out_beta2))
  
  
}





##### Test some data set ###########

# sim <- simulate_claims(50, 200, "graph", FALSE, mixing = "gamma", model_type = "poisson", exposure_lambda = 0,area = 12)
# # 
# # 
# #  out_poisson <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, sim$A, FALSE, "learn_graph",
# #                         lambda = 0, sim$exposure,
# #                         nr_em = 100, max_itr = 1000, a_known = FALSE)
# # 
# # beta1 <- out_poisson$beta1
# # 
# # 
# # 
# out_mixed <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, sim$A, FALSE, "learn_graph", lambda = 100,
#               exposure = sim$exposure, max_itr = 0, mixing_var = "gamma", nr_em = 60,Q_tol = 0,
#               verbose = 2, sgd = FALSE, batch_size = 100, param_tol = 0, a_known = FALSE)
# # 
# 
# sum(abs(out_poisson$beta1 - sim$beta1))
# sum(abs(out_mixed$beta1 - sim$beta1))
# 
# 
# sum(abs(out_poisson$a - sim$a))
# sum(abs(out_mixed$a - sim$a))
# 
# 
# diag(solve(out_mixed$Hessian - out_mixed$var_loglik))
