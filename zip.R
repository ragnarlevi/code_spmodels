# Source code for zip models
source("utils.R")

# Objective of ZIP
log_dZIP <- function(x,mu,prop, z){
  lab <- (x==0)*1
  lab*log(prop + (1-prop)*exp(-mu*z)) + (1-lab)*(log((1-prop)) -mu*z  + x*log(mu*z) - lgamma(x+1))
}

log_dZIP_der <- function(x,mu,prop, z){
  lab <- (x == 0)*1
  
  const1 <- 1/(prop + (1-prop)*exp(-mu*z))
  
  
  fphi <- lab*const1*(1 - exp(-mu*z)) - (1-lab)/(1-prop)
  fmu <- -lab*z*const1*(1-prop)*exp(-mu*z) + (1-lab)*(x/mu - z)
  
  return(list(fmu = fmu,fphi = fphi))
}


 

Q_ZIP_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, a_known = FALSE){
  
  etheta <- as.numeric(etheta)
  
  prop <- param[length(param)]
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param)-1)]
    a <- 1
  }else if(model_type == "learn_graph" & !a_known){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param)-1)]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  if(a_known){
    beta1 <- param[1:ncol(X_mat)]
    psi <- NA
  }
  
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  
  derivs <-  log_dZIP_der(claims, mu, prop, 1)
  fmu <- as.matrix(derivs$fmu)
  fphi_deriv <- derivs$fphi
  
  # gradient for prop
  g3 <- sum(fphi_deriv)
  
  
  
  # add gradients depending on model type
  if(model_type == "ordinary"){
    grad <- t( fmu * nu * exp(log(exposure)) ) %*% X_mat
    
  }else if(model_type == "learn_psi"){
    
    if(additive){
      #beta
      g1 <- t( fmu * nu * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- fmu[,1]*exp(log(exposure)) * se$agg_effect
      g2 <- tapply(g2,locs,sum)
    }else{
      # beta
      g1 <- t( fmu * nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      # psi
      g2 <- fmu[,1]*exp(log(exposure)) * se$agg_effect * nu
      g2 <- tapply(g2,locs,sum)
    }
    
    grad <- c(g1, g2)
    
    
  }else if(model_type == "learn_graph"){
    if(additive){
      g1 <- t( fmu * nu * exp(log(exposure)) ) %*% X_mat
      
      # a
      g2 <- fmu[,1]*t(agg_claims[, years]) * exp(log(exposure))
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    }else{
      # beta
      g1 <- t( fmu * nu*(1+se$spatial_effect) * exp(log(exposure)) ) %*% X_mat
      
      # A
      g2 <-  fmu[,1]*t(agg_claims[, years]) * exp(log(exposure)) *as.numeric(nu)
      
      g2 <- by(g2,locs, FUN=colSums)
      G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
      g22 <- G2[upper.tri(G2, diag = T)]
      diag(G2) = 0  # make sure we do not double count the diagonal when we add
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
      
    }
    
    if(a_known){
      g22 <- c()
    }
    
    grad <- c(g1, g22)
  }
  
  
  return(-c(grad, g3))  # adding prop deriv
  
  
}



Q_ZIP <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, 
                  additive, model_type, lambda, p, a_known = FALSE){
  
  
  etheta <- as.numeric(etheta)
  
  prop <- param[length(param)]
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param)-1)]
    a <- 1
  }else if(model_type == "learn_graph" & !a_known){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param)-1)]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  if(a_known){
    beta1 <- param[1:ncol(X_mat)]
    psi <- NA
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  ll <- sum(log_dZIP(claims, mu, prop, etheta))
  
  return(-ll)
  
}



 

zip <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure,  lambda = 0,  max_itr = 1000, a_known = FALSE, calc_hessian = TRUE){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  
  
  if(model_type == "learn_graph"& !a_known){
    a <- rep(1,p*(p+1)/2)*0.001
    A <- get_W_from_array(a, p)
    lower_a <- rep(1e-4, p*(p+1)/2)
    upper_a <- rep(0.5, p*(p+1)/2)
  }else if(model_type == "learn_psi"){
    psi <- rep(1,p)
  }
  

  # add initial param depending on model type
  if(model_type == "ordinary"){
    theta0 <- c(beta1, 0.3)
    lower <- c(rep(-20,ncol(X1)), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)), 1-1e-3)
    psi <- NA
    A <- NA
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi, 0.3)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)+p), 1-1e-3)
  }else if(model_type == "learn_graph"& !a_known){
    theta0 <- c(beta1, a, 0.3)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)+p*(p+1)/2), 1-1e-3)
    psi <- NA
  }else if(a_known){
    theta0 <- c(beta1, 0.3)
    lower <- c(rep(-20,ncol(X1)), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)), 1-1e-3)
    psi <- NA
  }else{
    stop("model type not known")
  }
  
  
  
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X1 %*% beta1))
  mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
  
  
  etheta <- 1
  
  ### ZIP poisson part
  out <- optim(par = theta0,
               fn = Q_ZIP,
               gr = Q_ZIP_deriv,
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
               method = 'L-BFGS-B',
               control = list(maxit = max_itr),
               lower = lower,
               upper = upper,
               hessian = calc_hessian)
  

  if(model_type == "ordinary"){
    beta1 <- out$par[1:ncol(X1)]
  }else if(model_type == "learn_psi"){
    beta1 <- out$par[1:ncol(X1)]
    psi <- out$par[(ncol(X1)+1):(length(out$par)-1)]
  }else if(model_type == "learn_graph"& !a_known ){
    beta1 <- out$par[1:ncol(X1)]
    a <- out$par[(ncol(X1)+1):(length(out$par)-1)]
    A <- get_W_from_array(a, p)
  }else if(a_known){
    beta1 <- out$par[1:ncol(X1)]
    a <- A[upper.tri(A, TRUE)]
  }
  
  prop <- out$par[length(out$par)]
  
  
  # Find number of params
  if(a_known){
    nr_param <- length(beta1) + 1
  }else if(model_type == "learn_graph"){
    nr_param <- length(beta1) + 1 + sum(abs(A[upper.tri(A, diag = TRUE)]) > 1e-3)
  }else if(model_type == "learn_psi"){
    nr_param <- length(beta1) + 1 + length(psi)
  }
  
  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], 
              beta2 = 1 , H = out$hessian, H_beta2 = NA, mu = mu, optim_obj = out, 
              model_type = model_type, prop = prop, nr_param = nr_param))
  
  
}

########## Mixed Models ##############

# we need integral f(y|z) f(z)
fyzz <- function(y,mu,prop, phi){
  
  int_fun <- function(z, mu, y) exp(log_dZIP(y,mu,prop, z))*dgamma(z,shape = phi, rate = phi )
  
  y_max <- max(y)
  
  denominator <- mapply(function(y, mu){pracma::integral(int_fun, 1e-3, y_max*10, mu = mu, y = y )}, 
                        y = y, mu = mu)
  
  
  return(denominator)
  
}



E_fzy <- function(y,mu, mu_old, prop, phi){
  
  
  int_fun <- function(z, mu, mu_old, y) log_dZIP(y,mu,prop, z)*exp(log_dZIP(y,mu_old,prop, z))*dgamma(z,shape = phi, rate = phi )
  
  y_max <- max(y)
  
  num <- mapply(function(y, mu, mu_old){integrate(int_fun, 1e-3, y_max*10, mu = mu, mu_old = mu_old, y = y )$value}, 
                        y = y, mu = mu, mu_old = mu_old)
  
  
  return(num)
  
}


log_dZIP_der_mu <- function(y,mu, mu_old, prop){
  
  grad_fun <- function(z, y, mu, prop){
    lab <- (y == 0)*1
    
    const1 <- 1/(prop + (1-prop)*exp(-mu*z))
    
    
    fmu <- -lab*z*const1*(1-prop)*exp(-mu*z) + (1-lab)*(y/mu - z)
    
    return(fmu) 
  }
  
  # Integrate the Gradient
  int_fun <- function(z, mu, mu_old, y) grad_fun(z, y, mu, prop)*exp(log_dZIP(y,mu_old,prop, z))*dgamma(z,shape = phi, rate = phi )
  
  y_max <- max(y)
  
  out <- mapply(function(y, mu, mu_old){pracma::integral(int_fun, 1e-3, y_max*10, mu = mu, mu_old = mu_old, y = y )}, 
                y = y, mu = mu, mu_old = mu_old)
  
  
  denominator <- fyzz(y,mu,prop, phi)
  
  return(log_dZIP_der_mu/denominator)
  
  
}



Q_ZIP_mixed <- function(param, locs, claims, exposure, X_mat, agg_claims, years, A, additive, model_type, lambda, p, mu_old, phi, prop, a_known = FALSE){
  

  
  prop <- param[length(param)]
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param)-1)]
    a <- 1
  }else if(model_type == "learn_graph" & !(a_known)){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param)-1)]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  if(a_known){
    beta1 <- param[2:length(param)]
    a <- 0
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  Q_val <- E_fzy(y, mu,mu_old,prop, phi)
  
  ll <- sum(Q_val)
  
  return(-ll)
  
}


# 
# # Poisson mixture model 
# zip_mixed <- function(claims, X1, locs, years,  agg_claims, A, additive, model_type, exposure, lambda = 0, nr_em = 100, max_itr = 1000, mixing_var = "gamma", z= "", a_known = FALSE){
#   
#   
#   if(model_type != "learn_graph"){
#     a_known <- FALSE
#   }
#   
#   
#   p <- length(unique(locs))
#   
#   # Set initial parameters
#   out_zip <- zip(claims, X1, locs, years, agg_claims, A, additive, model_type, lambda = lambda, exposure, nr_em = 100, max_itr = 300)
#   beta1 <- out_zip$beta1
#   prop <- out_zip$prop
#   
#   
#   if(model_type == "learn_graph"){
#     a <- rep(1,p*(p+1)/2)*0.001
#     A <- get_W_from_array(a, p)
#     lower_a <- rep(1e-4, p*(p+1)/2)
#     upper_a <- rep(0.5, p*(p+1)/2)
#   }else if(model_type == "learn_psi" ){
#     psi <- rep(1,p)
#   }
#   
#   
#   
#   
#   # add gradients depending on model type
#   if(model_type == "ordinary"){
#     theta0 <- c(beta1, prop)
#     lower <- c(rep(-20,ncol(X1)))
#   }else if(model_type == "learn_psi"){
#     theta0 <- c(beta1, psi, prop)
#     lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
#   }else if(model_type == "learn_graph" & !a_known){
#     theta0 <- c(beta1, a, prop)
#     lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2))
#     psi <- NA
#   }else if(a_known){
#     theta0 <- c(beta1)
#     lower <- c(rep(-20,ncol(X1)))
#     psi <- NA
#   }else{
#     stop("model type not known")
#   }
#   
#   if(a_known){
#     theta0 <- c(beta1)
#     a <- A[upper.tri(A, TRUE)]
#   }
#   
#   # set initial for the mixing
#   if(mixing_var == "gamma"){
#     beta20 <- 1
#   }else if(mixing_var == "ig"){
#     beta20 <- 0
#   }else if(mixing_var == "ln"){
#     beta20 <- 0
#   }else{
#     stop("mixing type not known")
#   }
#   
#   
#   X2 <-  as.matrix(rep(1, nrow(X1)))
#   
#   
#   for(i in 1:nr_em){
#     
#     
#     se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
#     nu <- as.numeric(exp(X1 %*% beta1))
#     mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
#     phi <- exp(beta20)
#     
# 
#  
#     
#     ### zip part
#     out <- optim(par = theta0,
#                  fn = Q_ZIP_mixed,
#                  #gr = Q_P_deriv,
#                  locs = locs,
#                  claims = claims, 
#                  exposure = exposure, 
#                  X_mat = X1, 
#                  agg_claims = agg_claims, 
#                  years = years,
#                  A = A,
#                  additive = additive,
#                  model_type = model_type,
#                  lambda = lambda,
#                  p = p,
#                  mu_old = mu,
#                  phi = phi,
#                  prop = prop,
#                  a_known = a_known,
#                  method = 'L-BFGS-B',
#                  control = list(maxit = ifelse(mixing_var == "none", 1000, 1)),
#                  lower = lower,
#                  hessian = T)
#     
#     
#     ## mixing part
#     
#     out_beta2 = optim(beta20,
#                       Q_beta_2,
#                       gr = Q_beta_2_deriv,
#                       claims = claims, 
#                       exposure =exposure, 
#                       X2 = X2, 
#                       etheta = etheta,
#                       eother = eother,
#                       method = 'L-BFGS',
#                       control = list(maxit = 2),
#                       hessian = TRUE)
#     
#     beta2 <- out_beta2$par
#     H_beta2 <- out_beta2$hessian
#     
#     
#     
#     if(model_type == "ordinary"){
#       beta1 <- out$par[1:ncol(X1)]
#       theta <- beta1
#     }else if(model_type == "learn_psi"){
#       beta1 <- out$par[1:ncol(X1)]
#       psi <- out$par[(ncol(X1)+1):length(out$par)]
#       theta <- c(beta1, psi)
#     }else if(model_type == "learn_graph"  & !a_known ){
#       beta1 <- out$par[1:ncol(X1)]
#       a <- out$par[(ncol(X1)+1):length(out$par)]
#       A <- get_W_from_array(a, p)
#       theta <- c(beta1, a)
#     }
#     
#     if(a_known){
#       beta1 <- out$par
#       theta <- beta1
#     }
#     
#     # check convergence
#     
#     if(sum(abs(theta-theta0))/sum(abs(theta0)) < 1e-5 &  sum(abs(beta2-beta20))/sum(abs(beta20)) < 1e-5 ){
#       break
#     }else{
#       theta0 <- theta
#       beta20 <- beta2
#     }
#     
#     
#     
#   }
#   
#   
#   
#   if(model_type == "ordinary"){
#     beta1 <- out$par[1:ncol(X1)]
#   }else if(model_type == "learn_psi"){
#     beta1 <- out$par[1:ncol(X1)]
#     psi <- out$par[(ncol(X1)+1):length(out$par)]
#   }else if(model_type == "learn_graph"){
#     beta1 <- out$par[1:ncol(X1)]
#     a <- out$par[(ncol(X1)+1):length(out$par)]
#     A <- get_W_from_array(a, p)
#   }
#   
#   
#   return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], beta2 = beta2 , H = out$hessian, H_beta2 = H_beta2, mu = mu, optim_obj = out, model_type = model_type))
#   
#   
# }
# 
# 
# 
# 


