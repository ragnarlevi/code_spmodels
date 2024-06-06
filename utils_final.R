

# Utils for Spatial temporal models ----


min_max_scale <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}



log_dnb <- function(x,mu,phi){
  lgamma(phi + x) - lgamma(x+1) - lgamma(phi) + phi*log(phi) + x*log(mu) - (phi+x)*log(mu+phi)
}


# new spatial aggregate function, Now assumes that the first column of agg_claims is claims at time t=0
get_spatial_aggregate_new <- function(locs, w, psi, agg_claims, years){
  
  if(all(psi == 0)){
    return(data.frame(locs = locs,years = years, agg_effect = 0, spatial_effect = 0))
  }
  
  se <- data.frame(locs = locs,years = years)
  
  # index the psi through locations
  if(length(psi) == 1){
    se$psi <- 1
  }else{
    se$psi <- psi[locs] 
  }
  
  tmp <- as.matrix(w) %*% agg_claims
  
  # index the agg_effect through locations and years
  se$agg_effect <- tmp[cbind(locs,years)] 
  se$spatial_effect <-  se$agg_effect*se$psi
  return(se)
}


# Function to build matrix from vector
get_W_from_array = function(a, p){
  
  A = matrix(0, nrow =p, ncol = p)
  A[upper.tri(A, diag = T)] = a
  A = A +  t(A)
  diag(A) = diag(A)*0.5
  
  return(A)
  
}



# E-step for PIG
Etheta_gig <- function(claims, a, b){
  
  
  num <- (sqrt(b)/sqrt(a)) * besselK(sqrt(a*b), claims+1-0.5,expon.scaled = TRUE)
  denom <- besselK(sqrt(a*b), claims-0.5,expon.scaled = TRUE)
  
  return(num/denom)
}

Einvtheta_gig <- function(claims, a, b){
  num <- (sqrt(a)/sqrt(b)) * besselK(sqrt(a*b), claims+1-0.5,expon.scaled = TRUE)
  denom <- besselK(sqrt(a*b), claims-0.5,expon.scaled = TRUE)
  
  
  return(num/denom- 2*(claims-0.5)/b)
  
}

# Functions for Poisson Gamma ----


# E-step for Poisson-Gamma
Etheta <- function(claims,lambda,phi) (phi + claims)/(phi + lambda)
Elogtheta <- function(claims,lambda,phi) digamma(phi + claims) - log(phi + lambda)



# functions to perform integrations for pln
dlno<-function(z,phi) {
  
  (1/(sqrt(2*pi)*phi*z) ) *exp(  -((log(z)+((phi^2)/2))^2)  /(2*(phi^2)) )
  
}


jpdf<-function(y,mu,phi, max_int = Inf){
  fun=function(z) {
    dpois(y, mu*z)*dlno(z,phi)
  }
  
  out <- tryCatch({
    integrate(fun,0,max_int)$value 
  }, error = function(e) return(""))
  
  if(is.character(out)){
    x <- NA
  }else{
    x <- out
  }
  
  return(x)
}

jpdf <- Vectorize(jpdf, vectorize.args = c("y", "mu", "phi", "max_int"))

pm2<-function(y,mu,phi, numer, max_int = Inf){
  fun=function(z) {
    (log(z)^2)* dpois(y, mu*z)* dlno(z,phi)
  }
  
  
  out <- tryCatch({
    integrate(fun,0,max_int)$value/numer
  }, error = function(e) return(""))
  
  if(is.character(out)){
    x <- NA
  }else{
    x <- out
  }
  
  return(x)
}

pm2 <- Vectorize(pm2, vectorize.args = c("y", "mu", "phi","numer", "max_int"))

pm1<-function(y,mu,phi, numer, max_int = Inf){
  fun=function(z) {
    z * dpois(y, mu*z)* dlno(z,phi)
  }
  out <- tryCatch({
    integrate(fun,0,max_int)$value/numer
  }, error = function(e) return(""))
  
  if(is.character(out)){
    x <- NA
  }else{
    x <- out
  }
  
  return(x)
}
pm1 <- Vectorize(pm1, vectorize.args = c("y", "mu", "phi","numer", "max_int"))




log_dpln <- function(y,mu,phi, max_int = Inf){
  fun=function(z) {
    dpois(y, mu*z)* dlno(z,phi)
  }
  out <- tryCatch({
    integrate(fun,0,max_int)$value
  }, error = function(e) return(""))
  
  if(is.character(out)){
    x <- NA
  }else{
    x <- out
  }
  
  return(log(x))
  
}
log_dpln <- Vectorize(log_dpln, vectorize.args = c("y", "mu", "phi", "max_int"))


# Poisson models -----


########### Following 4 functions are Poisson models, no mixing #####################
# X1 is data matrix for mean
# locs - vector containing location of each observation from 1 to p, where p is number of regions
# years - vector specifying year/time of each observation goes from 1 to T, where T is max time
# exposure - exposure
# do_optim - TRUE then each M-step is optimization with L-BFGS, else beta_1 and beta_2 optimized with newton
# only_diag - only used if newton is performed, If true then only the diagonal of the Hessian is used to do newton update - much quicker
# max_itr - max EM steps
# step_size_beta_1 - only if newton, step size of newton for beta_1 (mean)
# step_size_beta_2 - only if newton, step size of newton for beta_2 (dispersion)



# Objective of the additive Poisson with known graph w.r.t to beta_1 treating psi_1 as known
Q_P_beta_1 <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, psi1, w1){
  
  beta1 <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w1, psi1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta1))
  mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  
  ll <- sum(claims*log(mu) -mu*etheta)
  
  
  
  return(-ll)
  
}




# deriviative of the objective of the additive Poisson with known graph w.r.t to beta_1 treating psi_1 as known
Q_P_beta_1_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, psi1, w1){
  
  beta1 <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w1, psi1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  
  der1_logdp <- claims/mu - etheta
  
  g1 <- t( der1_logdp * nu1 * exp(log(exposure)) ) %*% X_mat
  
  
  
  return(-g1)
  
}

# Poisson model with no spatial features
Poisson_no_spatial <- function(claims, X1, exposure = rep(1, length(claims)), do_optim =T, only_diag = T, max_itr = 1000, 
                               step_size_beta_1 = 0.1, step_size_beta_2 = 0.1){
  
  
  # Set initial parameters
  #beta1_pre <- glm(claims ~ -1 + X1 + offset(exposure), family = 'poisson')$coefficients
  beta1_pre <- c(log(sum(claims)/sum(exposure)),rep(0, ncol(X1)-1))
  
  nu1 <- exp(X1 %*% as.matrix(beta1_pre))
  mu <- nu1*exp(log(exposure))
  
  ll <-  sum(dpois(claims,mu,log = T))
  
  out_beta1_deriv_em = optim(beta1_pre,
                             Q_P_beta_1,
                             gr = Q_P_beta_1_deriv,
                             locs = NA,
                             claims = claims, 
                             exposure =exposure, 
                             X_mat = X1, 
                             agg_claims = 0, 
                             years = NA, 
                             etheta = 1, 
                             psi1 = 0,
                             w1 = 0,
                             method = 'L-BFGS',
                             control = list(maxit = 200), hessian = T)
  
  
  beta1 = out_beta1_deriv_em$par
  beta1_H <- out_beta1_deriv_em$hessian
  
  
  # Create new estimates
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- nu1*exp(log(exposure))
  
  
  
  
  
  
  
  return(list(beta1 = out_beta1_deriv_em$par, beta1_H = beta1_H,beta1_pre=beta1_pre ))
  
  
  
}

# Poisson model with additive psi
Poisson_additive_psi <- function(claims, X1, locs, years, agg_claims, w, exposure = rep(1, length(claims)), do_optim =T, only_diag = T, 
                                 max_itr = 1000, step_size_beta_1 = 0.1, step_size_beta_2 = 0.1, z = ""){
  
  # do_optim if use L-BFGS for beta1 and beta2 updates
  # only_diag if the only the diagonal of the Hessian should be used, only used if do_optim = FALSE
  
  p_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, w){
    
    beta1 <- param[1:ncol(X_mat)]
    psi1 <- param[(ncol(X_mat)+1):length(param)]
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
    nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta1)))
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    
    ll <- sum(claims*log(mu) -mu)
    
    
    
    return(-ll)
    
  }
  
  
  p_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, w){
    
    beta1 <- param[1:ncol(X_mat)]
    psi1 <- param[(ncol(X_mat)+1):length(param)]
    
    se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
    nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta1)))
    
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    
    der1_logdp <- claims/mu - 1
    
    g1 <- t( der1_logdp * nu1 * exp(log(exposure)) ) %*% X_mat
    
    # update psi1
    g2 <- (claims/mu - 1)*exp(log(exposure)) * se1$agg_effect
    g2 <- tapply(g2,locs,sum)
    
    
    
    return(-c(g1,g2))
    
  }
  
  
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  psi1 <- rep(1,p)*0.001
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  
  
  out_beta1_deriv_em = optim(c(beta1, psi1),
                             p_obj,
                             gr = p_obj_deriv,
                             locs = locs,
                             claims = claims, 
                             exposure = exposure, 
                             X_mat = X1, 
                             agg_claims = agg_claims, 
                             years = years,
                             w = w,
                             method = 'L-BFGS-B',
                             control = list(maxit = max_itr),
                             lower = c(rep(-20,ncol(X1)), rep(1e-8, p)),
                             hessian = T)
  
  
  
  
  print("finalising")
  
  beta1 <- out_beta1_deriv_em$par[1:ncol(X1)]
  psi1 <- out_beta1_deriv_em$par[(ncol(X1)+1):length(out_beta1_deriv_em$par)]
  
  
  
  
  return(list(beta1 = beta1, psi = psi1, H = out_beta1_deriv_em$hessian, mu = mu))
  
  
}

# Poisson model with multiplicative psi
Poisson_multiplicative_psi <- function(claims, X1, locs, years, agg_claims, w, exposure = rep(1, length(claims)), do_optim =T, only_diag = T, 
                                       max_itr = 1000, step_size_beta_1 = 0.1, step_size_beta_2 = 0.1, z = ""){
  
  # do_optim if use L-BFGS for beta1 and beta2 updates
  # only_diag if the only the diagonal of the Hessian should be used, only used if do_optim = FALSE
  
  p_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, w){
    
    beta1 <- param[1:ncol(X_mat)]
    psi1 <- param[(ncol(X_mat)+1):length(param)]
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
    nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta1)))
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    
    ll <- sum(claims*log(mu) -mu)
    
    
    
    return(-ll)
    
  }
  
  
  p_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, w){
    
    beta1 <- param[1:ncol(X_mat)]
    psi1 <- param[(ncol(X_mat)+1):length(param)]
    
    se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
    nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta1)))
    
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    der1_logdp <- claims/mu - 1
    
    
    g1 <- t( der1_logdp * nu1*(1+se1$spatial_effec) * exp(log(exposure)) ) %*% X_mat
    # update psi1
    g2 <- (claims/mu - 1)*exp(log(exposure)) * se1$agg_effect * nu1
    g2 <- tapply(g2,locs,sum)
    
    
    
    return(-c(g1,g2))
    
  }
  
  
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  # beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta1 <- c(log(sum(claims)/sum(exposure)),rep(0, ncol(X1)-1))
  psi1 <- rep(1,p)*1e-3
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  
  
  out_beta1_deriv_em = optim(c(beta1, psi1),
                             p_obj,
                             #gr = p_obj_deriv,
                             locs = locs,
                             claims = claims, 
                             exposure = exposure, 
                             X_mat = X1, 
                             agg_claims = agg_claims, 
                             years = years,
                             w = w,
                             method = 'L-BFGS-B',
                             control = list(maxit = max_itr),
                             lower = c(rep(-20,ncol(X1)), rep(1e-8, p)),
                             hessian = T)
  
  
  
  
  print("finalising")
  
  beta1 <- out_beta1_deriv_em$par[1:ncol(X1)]
  psi1 <- out_beta1_deriv_em$par[(ncol(X1)+1):length(out_beta1_deriv_em$par)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  
  
  
  return(list(beta1 = beta1, psi = psi1, H = out_beta1_deriv_em$hessian, mu = mu, out = out_beta1_deriv_em))
  
  
}



# Poisson model with additive spatial effect where graph is learnt
Poisson_additive_a <- function(claims, X1, locs, years, agg_claims, exposure = rep(1, length(claims)), do_optim =T, only_diag = T, lambda = 0,
                               max_itr = 1000, step_size_beta_1 = 0.1, step_size_beta_2 = 0.1, z = ""){
  
  # do_optim if use L-BFGS for beta1 and beta2 updates
  # only_diag if the only the diagonal of the Hessian should be used, only used if do_optim = FALSE
  
  p_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, lambda){
    
    beta1 <- param[1:ncol(X_mat)]
    a1 <- param[(ncol(X_mat)+1):length(param)]
    A1 = get_W_from_array(a1,p)
    
    
    spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
    nu1 <- exp(X_mat %*% as.matrix(beta1))
    mu <- (nu1 + spatial_effect_mu)*exp(log(exposure))
    
    
    
    
    
    ll <- sum(claims*log(mu) -mu) - sum(a1*lambda)
    
    
    
    return(-ll)
    
  }
  
  
  p_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, lambda){
    
    beta1 <- param[1:ncol(X_mat)]
    a1 <- param[(ncol(X_mat)+1):length(param)]
    A1 = get_W_from_array(a1,p)
    
    spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
    nu1 <- exp(X1 %*% as.matrix(beta1))
    mu <- (nu1 + spatial_effect_mu)*exp(log(exposure))
    
    
    
    
    
    der1_logdp <- claims/mu - 1
    
    g1 <- t( der1_logdp * nu1*spatial_effect_mu * exp(log(exposure)) ) %*% X_mat
    
    g2 <- der1_logdp[,1]*t(agg_claims[, years]) * exp(log(exposure))
    g2 <- by(g2,locs, FUN=colSums)
    G2 <- matrix(unlist(g2),nrow = nrow(A1), ncol = ncol(A1), byrow = T)
    g22 <- G2[upper.tri(G2, diag = T)]
    diag(G2) = 0  # make sure we do not double count the diagonal when we add
    g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    
    
    return(-c(g1,g22))
    
  }
  
  
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  #beta1 <- c(log(sum(claims)/sum(exposure)),rep(0, ncol(X1)-1))
  a1 = rep(0.1, p*(p+1)/2)
  A1 = get_W_from_array(a1, p)
  
  
  
  spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 + spatial_effect_mu)*exp(log(exposure))
  
  
  
  out_beta1_deriv_em = optim(c(beta1, a1),
                             p_obj,
                             gr = p_obj_deriv,
                             locs = locs,
                             claims = claims, 
                             exposure = exposure, 
                             X_mat = X1, 
                             agg_claims = agg_claims, 
                             years = years,
                             method = 'L-BFGS-B',
                             control = list(maxit = max_itr),
                             lambda =lambda,
                             lower = c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2)),
                             hessian = T)
  
  
  
  
  print("finalising")
  
  beta1 <- out_beta1_deriv_em$par[1:ncol(X1)]
  a1 <- out_beta1_deriv_em$par[(ncol(X1)+1):length(out_beta1_deriv_em$par)]
  A1 <- get_W_from_array(a1,p) 
  
  
  
  spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 + spatial_effect_mu)*exp(log(exposure))
  
  
  
  return(list(beta1 = beta1, A = A1, a =a1, H = out_beta1_deriv_em$hessian, mu = mu))
  
  
}




# Poisson model with multiplicative spatial effect where graph is learnt
Poisson_multiplicative_a <- function(claims, X1, locs, years, agg_claims, exposure = rep(1, length(claims)), do_optim =T, only_diag = T, lambda = 0,
                                     max_itr = 1000, step_size_beta_1 = 0.1, step_size_beta_2 = 0.1, z = ""){
  
  # do_optim if use L-BFGS for beta1 and beta2 updates
  # only_diag if the only the diagonal of the Hessian should be used, only used if do_optim = FALSE
  
  p_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, lambda){
    
    beta1 <- param[1:ncol(X_mat)]
    a1 <- param[(ncol(X_mat)+1):length(param)]
    A1 = get_W_from_array(a1,p)
    
    
    spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
    nu1 <- exp(X_mat %*% as.matrix(beta1))
    mu <- (nu1 * (1+spatial_effect_mu))*exp(log(exposure))
    
    ll <- sum(claims*log(mu) -mu) - sum(a1*lambda)
    
    return(-ll)
    
  }
  
  
  p_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, lambda){
    
    beta1 <- param[1:ncol(X_mat)]
    a1 <- param[(ncol(X_mat)+1):length(param)]
    A1 = get_W_from_array(a1,p)
    
    spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
    nu1 <- exp(X1 %*% as.matrix(beta1))
    mu <- (nu1 * (1+spatial_effect_mu))*exp(log(exposure))
    
    
    
    der1_logdp <- claims/mu - 1
    g1 <- t( der1_logdp * nu1 * exp(log(exposure))*(1+spatial_effect_mu) ) %*% X_mat
    
    
    # find derivative
    der1_logdp <- as.numeric(claims/mu - 1)
    
    # update w1. Note how I try to sum correclt for each parameter in the matrix
    # The goal is to have da_{ij} <-  da_{ij}+ da_{ji}, but for the diagonal we only have da_{ii} <- da_{ii}, thus the diag(G2)<-0
    g2 <- der1_logdp*t(agg_claims[, years]) * exp(log(exposure)) *as.numeric(nu1)
    g2 <- by(g2,locs, FUN=colSums)
    G2 = matrix(unlist(g2),nrow = nrow(A1), ncol = ncol(A1), byrow = T)
    g22 = G2[upper.tri(G2, diag = T)]
    diag(G2) = 0  # make sure we do not double count the diagonal when we add
    g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    
    
    return(-c(g1,g22))
    
  }
  
  
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  #beta1 <- c(log(sum(claims)/sum(exposure)),rep(0, ncol(X1)-1))
  a1 = rep(0.1, p*(p+1)/2)
  A1 = get_W_from_array(a1, p)
  
  
  
  spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 * (1+spatial_effect_mu))*exp(log(exposure))
  
  
  
  out_beta1_deriv_em = optim(c(beta1, a1),
                             p_obj,
                             gr = p_obj_deriv,
                             locs = locs,
                             claims = claims, 
                             exposure = exposure, 
                             X_mat = X1, 
                             agg_claims = agg_claims, 
                             years = years,
                             method = 'L-BFGS-B',
                             control = list(maxit = max_itr),
                             lambda =lambda,
                             lower = c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2)),
                             hessian = T)
  
  
  
  
  print("finalising")
  
  beta1 <- out_beta1_deriv_em$par[1:ncol(X1)]
  a1 <- out_beta1_deriv_em$par[(ncol(X1)+1):length(out_beta1_deriv_em$par)]
  A1 <- get_W_from_array(a1,p) 
  
  
  
  spatial_effect_mu = unname(rowSums(t(agg_claims[, years])*A1[locs,]))
  nu1 <- exp(X1 %*% as.matrix(beta1))
  mu <- (nu1 * spatial_effect_mu)*exp(log(exposure))
  
  
  
  return(list(beta1 = beta1, A1 = A1, a1 =a1, H = out_beta1_deriv_em$hessian, mu = mu, out = out_beta1_deriv_em))
  
  
}








# Ordindary Models -----



log_dP <- function(x,mu, z){
  #dpois(x, lambda = mu*z, log = TRUE)
  x*log(mu) -mu*z
}

log_dP_der <- function(x, mu, z){
  -z  + x/mu 
}



Q_P <- function(param, claims, exposure, X1, etheta){
  
  beta1 <- param
  etheta <- as.numeric(etheta)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  mu <- nu1*exp(log(exposure))
  ll_tmp <- log_dP(claims,mu,etheta)
  ll <- sum(log_dP(claims,mu,etheta))
  
  return(-ll)
  
}

Q_P_deriv <- function(param, claims, exposure, X1, etheta){
  
  
  beta1 <- param
  
  
  etheta <- as.numeric(etheta)
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  mu <- nu1*exp(log(exposure))
  
  
  # beta_deriv
  
  mu_deriv <- log_dP_der(claims, mu, etheta)
  
  beta_deriv <- as.numeric(t(as.numeric(mu_deriv)*nu1*exp(log(exposure))) %*% X1)
  
  
  
  
  
  return(-c(beta_deriv))
  
  
}



Q_P_psi <- function(param, locs, claims, exposure, X1, agg_claims, years, etheta, w, additive){
  
  beta1 <- param[1:ncol(X1)]
  psi <- param[(ncol(X1)+1):length(param)]
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll_tmp <- log_dP(claims,mu,etheta)
  
  
  ll <- sum(log_dP(claims,mu,etheta))
  
  return(-ll)
  
}

Q_P_a <- function(param, locs, claims, exposure, X1, agg_claims, years, etheta, additive, p, lambda){
  
  beta1 <- param[1:ncol(X1)]
  a <- param[(ncol(X1)+1):length(param)]
  
  etheta <- as.numeric(etheta)
  
  A <- get_W_from_array(a, p)  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  
  ll <- sum(log_dP(claims,mu,etheta)) - sum(lambda*a)
  
  return(-ll)
  
}


Q_P_psi_deriv <- function(param, locs, claims, exposure, X1, agg_claims, years, etheta, w, additive){
  
  
  beta1 <- param[1:ncol(X1)]
  psi <- param[(ncol(X1)+1):length(param)]
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  
  etheta <- as.numeric(etheta)
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  # beta_deriv
  
  mu_deriv <- log_dP_der(claims, mu, etheta)
  
  beta_deriv <- as.numeric(t(as.numeric(mu_deriv)*nu1*exp(log(exposure))) %*% X1)
  
  
  # psi deriv
  if(additive){
    psi_deriv <- mu_deriv*se1$agg_effect*exp(log(exposure))
  }else{
    psi_deriv <- mu_deriv*se1$agg_effect*exp(log(exposure))*nu1
  }
  
  psi_deriv <- tapply(psi_deriv,locs,sum)
  
  
  
  
  return(-c(beta_deriv, psi_deriv))
  
  
}


Q_P_a_deriv <- function(param, locs, claims, exposure, X1, agg_claims, years, etheta, additive, p, lambda){
  
  
  beta1 <- param[1:ncol(X1)]
  a <- param[(ncol(X1)+1):length(param)]
  
  
  A <- get_W_from_array(a, p)  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  
  
  
  etheta <- as.numeric(etheta)
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  # beta_deriv
  
  mu_deriv <- log_dP_der(claims, mu, etheta)
  
  beta_deriv <- as.numeric(t(as.numeric(mu_deriv)*nu1*exp(log(exposure))) %*% X1)
  
  
  
  
  # A
  if(additive){
    g2 <- as.numeric(mu_deriv)*t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <- as.numeric(mu_deriv)*t(agg_claims[, years]) * exp(log(exposure))*as.numeric(nu1)
  }
  
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  
  return(-c(beta_deriv, g22))
  
  
}



# PG -----
Q_PG_beta_2 <- function(param,  claims, exposure, X2,  etheta, elogtheta){
  
  etheta = as.numeric(etheta)
  elogtheta = as.numeric(elogtheta)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  
  ll = sum(phi*log(phi) -lgamma(phi) + (phi-1)*elogtheta - phi*etheta )
  return(-ll)
  
}

Q_PG_beta_2_deriv <- function(param,  claims, exposure, X2,  etheta, elogtheta){
  
  etheta = as.numeric(etheta)
  elogtheta = as.numeric(elogtheta)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  der1_logdga <- log(phi) + 1 - digamma(phi) + elogtheta - etheta
  
  
  g3 <- t( der1_logdga*nu2 ) %*% X2
  
  return(-g3)
}



PG_EM <- function(claims, X1,X2,  exposure = rep(1, length(claims)),  max_itr = 1000, tol = 1e-3){
  
  
  
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  beta1 <- beta1 -2*abs(beta1)
  
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  
  ll <-  sum(log_dnb(claims,mu,phi))
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    etheta <- Etheta(claims, mu, phi)
    elogtheta <- Elogtheta(claims,mu, phi)
    
    
    
    out_mean <- optim(c(beta1),
                      Q_P,
                      gr = Q_P_deriv,
                      claims = claims, 
                      exposure = exposure, 
                      X1 = X1, 
                      etheta = etheta, 
                      method = 'L-BFGS',
                      control = list(maxit = 20),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PG_beta_2,
                               gr = Q_PG_beta_2_deriv,
                               claims = claims, 
                               exposure =exposure, 
                               X2 = X2, 
                               etheta = etheta,
                               elogtheta = elogtheta,
                               method = 'L-BFGS',
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    mu <- nu1 *exp(log(exposure))
    phi <- nu2 
    
    new_ll <-  sum(log_dnb(claims,mu,phi))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    if(itr %% 10 == 0){
      print(paste('At iteration',itr,'log likelihood =',ll))
    }
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_deriv(beta1,  claims, exposure, X1,  etheta))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PG_beta_2_deriv(beta2,  claims, exposure, X2,  etheta, elogtheta))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2, beta1_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi))
  
}



PG_EM_psi <- function(claims, X1,X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)),  max_itr = 1000, tol = 1e-3){
  
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(-1, rep(0, ncol(X2) - 1)) 
  psi <- rep(1,p)*0.001
  lower_beta1 <- rep(-10, ncol(X1))
  lower_psi <- rep(1e-3, p)
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <-  sum(log_dnb(claims,mu,phi))
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    etheta <- Etheta(claims, mu, phi)
    elogtheta <- Elogtheta(claims,mu, phi)
    
    
    
    out_mean <- optim(c(beta1, psi),
                      Q_P_psi,
                      gr = Q_P_psi_deriv,
                      locs = locs, 
                      claims = claims, 
                      exposure = exposure, 
                      X1 = X1, 
                      agg_claims = agg_claims, 
                      years = years, 
                      etheta = etheta, 
                      w = w,
                      additive = additive,
                      method = 'L-BFGS-B',
                      lower = c(lower_beta1, lower_psi),
                      control = list(maxit = 20),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par[1:ncol(X1)]
    psi <- out_mean$par[(ncol(X1)+1):length(out_mean$par)]
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PG_beta_2,
                               gr = Q_PG_beta_2_deriv,
                               claims = claims, 
                               exposure =exposure, 
                               X2 = X2, 
                               etheta = etheta,
                               elogtheta = elogtheta,
                               method = 'L-BFGS',
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
    }
    phi <- nu2 
    
    new_ll <-  sum(log_dnb(claims,mu,phi))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    if(itr %% 10 == 0){
      print(paste('At iteration',itr,'log likelihood =',ll))
    }
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_psi_deriv(c(beta1, psi), locs, claims, exposure, X1, agg_claims, years, etheta, w, additive))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PG_beta_2_deriv(beta2, claims, exposure, X2,  etheta, elogtheta))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  
  
  return(list(beta1 = beta1, beta2 = beta2, psi = psi, mean_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, se =se1))
  
  
  
  
  
}


PG_EM_a <- function(claims, X1,X2, locs, years, agg_claims, additive, exposure = rep(1, length(claims)),  max_itr = 1000, lambda = 0, tol = 1e-3){
  
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  a  <- rep(1,p*(p+1)/2)*0.001
  A <- get_W_from_array(a,p)
  
  lower_beta1 <- rep(-10, ncol(X1))
  lower_a <- rep(1e-7, p*(p+1)/2)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <-  sum(log_dnb(claims,mu,phi))
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    etheta <- Etheta(claims, mu, phi)
    elogtheta <- Elogtheta(claims,mu, phi)
    
    
    
    out_mean <- optim(par = c(beta1, a),
                      fn = Q_P_a,
                      gr = Q_P_a_deriv,
                      locs = locs, 
                      claims = claims, 
                      exposure = exposure, 
                      X1 = X1, 
                      agg_claims = agg_claims, 
                      years = years, 
                      etheta = etheta, 
                      additive = additive,
                      p = p,
                      lambda = lambda,
                      method = 'L-BFGS-B',
                      lower = c(lower_beta1, lower_a),
                      control = list(maxit = 20),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par[1:ncol(X1)]
    a <- out_mean$par[(ncol(X1)+1):length(out_mean$par)]
    A <- get_W_from_array(a, p)
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PG_beta_2,
                               gr = Q_PG_beta_2_deriv,
                               claims = claims, 
                               exposure =exposure, 
                               X2 = X2, 
                               etheta = etheta,
                               elogtheta = elogtheta,
                               method = 'L-BFGS',
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
    }
    phi <- nu2 
    
    new_ll <-  sum(log_dnb(claims,mu,phi))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    if(itr %% 10 == 0){
      print(paste('At iteration',itr,'log likelihood =',ll))
    }
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_a_deriv(c(beta1, a), locs, claims, exposure, X1, agg_claims, years, etheta, additive, p, lambda))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PG_beta_2_deriv(beta2, claims, exposure, X2, etheta, elogtheta))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  
  
  return(list(beta1 = beta1, beta2 = beta2, a = a, A=A, mean_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, se = se1))
  
  
  
  
  
}


# PIG-----

Q_PIG_beta_2 <- function(param, claims, exposure, X2, etheta, einvtheta){
  
  etheta = as.numeric(etheta)
  einvtheta = as.numeric(einvtheta)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  #print(any(is.infinite(einvtheta)))
  
  ll = sum(log(phi) + phi^2 - 0.5*phi^2*(einvtheta + etheta))
  
  return(-ll)
  
}

Q_PIG_beta_2_deriv <- function(param, claims, exposure, X2, etheta, einvtheta){
  
  etheta = as.numeric(etheta)
  einvtheta = as.numeric(einvtheta)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  der1_logdgauss <- 1/phi + 2*phi - phi*einvtheta - phi*etheta
  g3 <- t(der1_logdgauss* nu2)%*% X2
  
  return(-g3)
}

PIG_EM <- function(claims, X1,X2, exposure = rep(1, length(claims)),  max_itr = 1000, z= "", tol = 1e-3){
  
  
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta1 <- rep(-1, ncol(X1))
  
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  
  ll <-  sum(log_dPIG(claims,mu,phi), na.rm = TRUE)
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    
    if(all(is.character(z))){
      # E-step
      etheta <- Etheta_gig(claims, phi^2 + 2*mu, phi^2)
      einvtheta <- Einvtheta_gig(claims, phi^2 + 2*mu, phi^2)
      
      
    }else{
      etheta <- z
      einvtheta <- 1/z
    }
    
    
    
    bool <- !(is.infinite(etheta) | is.infinite(einvtheta) | is.nan(etheta) | is.nan(einvtheta))
    
    etheta_tmp <- etheta[bool]
    einvtheta_tmp <- einvtheta[bool]
    X1_tmp <- X1[bool,,drop = T]
    X2_tmp <- X2[bool,,drop = T]
    exposure_tmp <- exposure[bool]
    claims_tmp <- claims[bool]
    
    
    #print(any(is.infinite(einvtheta)))
    # M-step
    out_mean <- optim(c(beta1),
                      Q_P,
                      gr = Q_P_deriv,
                      claims = claims_tmp, 
                      exposure = exposure_tmp, 
                      X1 = X1_tmp, 
                      etheta = etheta_tmp, 
                      method = 'L-BFGS',
                      control = list(maxit = 2),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par
    H_mean_par <- out_mean$hessian
    #print(beta2)
    out_beta2_deriv_em = optim(beta2,
                               Q_PIG_beta_2,
                               gr = Q_PIG_beta_2_deriv,
                               claims = claims_tmp, 
                               exposure =exposure_tmp, 
                               X2 = X2_tmp, 
                               etheta = etheta_tmp,
                               einvtheta = einvtheta_tmp,
                               method = 'L-BFGS',
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    #print(out_beta2_deriv_em$message)
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    mu <- nu1 *exp(log(exposure))
    phi <- nu2 
    
    tmp_ll <- log_dPIG(claims,mu,phi)
    
    new_ll <-  sum(log_dPIG(claims,mu,phi))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    if(itr %% 10 == 0){
      print(beta1)
      print(beta2)
      print(paste('At iteration',itr,'log likelihood =',ll))
    }
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_deriv(beta1, claims_tmp, exposure_tmp, X1_tmp, etheta_tmp))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PIG_beta_2_deriv(beta2, claims_tmp, exposure_tmp, X2_tmp,  etheta_tmp, einvtheta_tmp))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2, beta1_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, out_beta2 = out_beta2_deriv_em))
  
}


PIG_EM_psi <- function(claims, X1,X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)),  max_itr = 1000, z = "", tol = 1e-3){
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(1, rep(0, ncol(X2) - 1)) 
  psi <- rep(0,p)
  lower_beta1 <- rep(-10, ncol(X1))
  lower_psi <- rep(1e-8, p)
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  # ll_tmp <- log_dPIG(claims,mu,phi)
  # print(any(is.infinite(ll_tmp)))
  # print(claims[is.infinite(ll_tmp)])
  # print(mu[is.infinite(ll_tmp)])
  # print(phi[is.infinite(ll_tmp)])
  
  
  ll_tmp <- log_dPIG(claims,mu,phi)
  ll <-  sum(ll_tmp[!is.infinite(ll_tmp)], na.rm = T)
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    if(all(is.character(z))){
      etheta <- Etheta_gig(claims, phi^2 + 2*mu, phi^2)
      einvtheta <- Einvtheta_gig(claims, phi^2 + 2*mu, phi^2)
    }else{
      etheta <- z
      einvtheta <- 1/z
    }
    
    
    bool <- !(is.infinite(etheta) | is.infinite(einvtheta) | is.nan(etheta) | is.nan(einvtheta))
    
    etheta_tmp <- etheta[bool]
    einvtheta_tmp <- einvtheta[bool]
    X1_tmp <- X1[bool,,drop = T]
    X2_tmp <- X2[bool,,drop = T]
    exposure_tmp <- exposure[bool]
    years_tmp <- years[bool]
    claims_tmp <- claims[bool]
    locs_tmp <- locs[bool]
    
    #print(beta1)
    out_mean <- optim(c(beta1, psi),
                      Q_P_psi,
                      gr = Q_P_psi_deriv,
                      locs = locs_tmp, 
                      claims = claims_tmp, 
                      exposure = exposure_tmp, 
                      X1 = X1_tmp, 
                      agg_claims = agg_claims, 
                      years = years_tmp, 
                      etheta = etheta_tmp, 
                      w = w,
                      additive = additive,
                      method = 'L-BFGS-B',
                      lower = c(lower_beta1, lower_psi),
                      control = list(maxit = 2),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par[1:ncol(X1)]
    psi <- out_mean$par[(ncol(X1)+1):length(out_mean$par)]
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PIG_beta_2,
                               gr = Q_PIG_beta_2_deriv,
                               claims = claims_tmp, 
                               exposure =exposure_tmp, 
                               X2 = X2_tmp, 
                               etheta = etheta_tmp,
                               einvtheta = einvtheta_tmp,
                               method = 'L-BFGS-B',
                               lower = rep(-1, ncol(X2)),
                               upper = rep(5, ncol(X2)),
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
    }
    phi <- nu2 
    
    new_ll_tmp <- log_dPIG(claims,mu,phi)
    new_ll <-  sum(new_ll_tmp[!is.infinite(new_ll_tmp)], na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    if(itr %% 10 == 0){
      print(paste('At iteration',itr,'log likelihood =',ll))
    }
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_psi_deriv(c(beta1, psi), locs, claims, exposure, X1, agg_claims, years, etheta, w, additive))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PIG_beta_2_deriv(beta2, claims, exposure, X2, etheta, einvtheta))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  
  
  return(list(beta1 = beta1, beta2 = beta2, psi = psi, mean_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, se =se1))
  
  
  
  
  
}


PIG_EM_a <- function(claims, X1,X2, locs, years, agg_claims, additive, exposure = rep(1, length(claims)),  max_itr = 1000, lambda = 0, z = "", tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  a  <- rep(1,p*(p+1)/2)*0.00001
  A <- get_W_from_array(a,p)
  
  lower_beta1 <- rep(-10, ncol(X1))
  lower_beta2 <- rep(-0.5, ncol(X2))
  lower_a <- rep(1e-7, p*(p+1)/2)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll_tmp <- log_dPIG(claims,mu,phi)
  ll <-  sum(ll_tmp[!is.infinite(ll_tmp)], na.rm = TRUE)
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    if(all(is.character(z))){
      etheta <- Etheta_gig(claims, phi^2 + 2*mu, phi^2)
      einvtheta <- Einvtheta_gig(claims, phi^2 + 2*mu, phi^2)
    }else{
      etheta <- z
      einvtheta <- 1/z
    }
    
    bool <- !(is.infinite(etheta) | is.infinite(einvtheta) | is.nan(etheta) | is.nan(einvtheta))
    
    etheta_tmp <- etheta[bool]
    einvtheta_tmp <- einvtheta[bool]
    X1_tmp <- X1[bool,,drop = T]
    X2_tmp <- X2[bool,,drop = T]
    exposure_tmp <- exposure[bool]
    years_tmp <- years[bool]
    claims_tmp <- claims[bool]
    locs_tmp <- locs[bool]
    
    
    
    out_mean <- optim(par = c(beta1, a),
                      fn = Q_P_a,
                      gr = Q_P_a_deriv,
                      locs = locs_tmp, 
                      claims = claims_tmp, 
                      exposure = exposure_tmp, 
                      X1 = X1_tmp, 
                      agg_claims = agg_claims, 
                      years = years_tmp, 
                      etheta = etheta_tmp, 
                      additive = additive,
                      p = p,
                      lambda = lambda,
                      method = 'L-BFGS-B',
                      lower = c(lower_beta1, lower_a),
                      control = list(maxit = 2),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par[1:ncol(X1)]
    a <- out_mean$par[(ncol(X1)+1):length(out_mean$par)]
    A <- get_W_from_array(a, p)
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PIG_beta_2,
                               gr = Q_PIG_beta_2_deriv,
                               claims = claims_tmp, 
                               exposure =exposure_tmp, 
                               X2 = X2_tmp, 
                               etheta = etheta_tmp,
                               einvtheta = einvtheta_tmp,
                               method = 'L-BFGS-B',
                               lower = lower_beta2,
                               upper = rep(5, ncol(X2)),
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(out_beta2_deriv_em$message)
    
    
    
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
    }
    phi <- nu2 
    
    new_ll_tmp <- log_dPIG(claims,mu,phi)
    new_ll <-  sum(new_ll_tmp[!is.infinite(new_ll_tmp)], na.rm = TRUE)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    if(itr %% 10 == 0){
      print(paste('At iteration',itr,'log likelihood =',ll))
    }
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_a_deriv(c(beta1, a), locs, claims, exposure, X1, agg_claims, years, etheta, additive, p, lambda))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PIG_beta_2_deriv(beta2, claims, exposure, X2, etheta, einvtheta))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  
  
  return(list(beta1 = beta1, beta2 = beta2, a = a, A=A, mean_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, se = se1))
  
  
  
  
  
}




# ML maximization


PIG_ML <- function(claims, X1,X2, exposure = rep(1, length(claims)),  max_itr = 1000, z= "", tol = 1e-3){
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta1 <- rep(-1, ncol(X1))
  
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  objective_pig <- function(param, X1, X2, claims, exposure){
    
    d1 <- ncol(X1)
    d2 <- ncol(X2)
    
    
    beta1 <- param[1:d1]
    #print(length(beta1))
    beta2 <- param[(d1+1):(d2+d1)]
    
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    mu <- (nu1 )*exp(log(exposure))
    phi <- nu2
    
    return(-sum(log_dPIG(claims, mu, phi), na.rm = T) )
    
    
    
  }
  
  objective_pig_deriv <- function(param, X1, X2, claims, exposure){
    
    d1 <- ncol(X1)
    d2 <- ncol(X2)
    
    
    beta1 <- param[1:d1]
    #print(length(beta1))
    beta2 <- param[(d1+1):(d2+d1)]
    
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    mu <- (nu1 )*exp(log(exposure))
    phi <- nu2
    
    
    nu <- claims-0.5
    delta <- sqrt(phi^2 + 2*mu)
    const <- (nu/(phi*delta) - besselK(phi*delta,nu+1,T)/besselK(phi*delta,nu,T))
    
    fmu <- as.matrix(claims/mu - nu/delta^2 + const*phi/delta )
    fphi <- as.matrix(1/phi + 2*phi + nu*(1/phi - phi/delta^2) +  const*(delta + phi^2/delta))
    
    
    bool_mu <- !(is.nan(fmu) | is.na(fmu) | is.infinite(fmu))
    bool_phi <- !(is.nan(fmu) | is.na(fmu) | is.infinite(fmu))
    
    
    #beta 1 deriv
    g1 <- t( fmu[bool_mu] * nu1[bool_mu]*exposure[bool_mu]) %*% X1[bool_mu,,drop = F]
    
    # beta2 deriv
    g3 <- t( fphi[bool_phi]) %*% X2[bool_phi,,drop =F]
    
    grad <- -c(g1, g3)
    return(grad)
    
    
    
  }
  
  out <- optim(par = c(beta1, beta2),
               fn = objective_pig,
               gr = objective_pig_deriv, 
               X1 = X1, X2 = X2, 
               claims = claims, 
               exposure = exposure,
               method = 'L-BFGS')
  
  
  beta1 <- out$par[1:ncol(X1)]
  beta2 <- out$par[(ncol(X1)+1):length(out$par)]
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  return(list(beta1 = beta1, beta2 = beta2, mu = mu, phi = phi, out = out))
  
  
  
}


PIG_ML_psi <- function(claims, X1,X2, additive, agg_claims, years, locs, w, exposure = rep(1, length(claims)),  max_itr = 1000, z= "", tol = 1e-3){
  
  # Set initial parameters
  p <- length(unique(locs))
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta1 <- beta1 - abs(beta1)
  upper_beta1 <- beta1 + abs(beta1)
  psi <- rep(1e-6,p)
  lower_psi <- rep(1e-8,p)
  upper_psi <- rep(0.5,p)
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exposure
  }else{
    mu <- (nu1*(1+se1$spatial_effect))*exposure
  }
  phi <- nu2
  
  objective_pig_add_psi <- function(param, X1, X2, p, claims, agg_claims, locs, years, w,additive,exposure){
    
    d1 <- ncol(X1)
    d2 <- ncol(X2)
    
    beta1 <- param[1:d1]
    #print(length(beta1))
    beta2 <- param[(d1+1):(d2+d1)]
    #print(length(beta2))
    psi1 <- param[(d1+d2+1):length(param)]
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exposure
    }else{
      mu <- (nu1*(1+se1$spatial_effect))*exposure
    }
    phi <- nu2
    
    return(-sum(log_dPIG(claims, mu, phi), na.rm = T) )
    
    
    
  }
  
  objective_pig_add_psi_deriv <- function(param, X1, X2, p, claims, agg_claims, locs, years, w, additive, exposure){
    
    d1 <- ncol(X1)
    d2 <- ncol(X2)
    
    
    beta1 <- param[1:d1]
    #print(length(beta1))
    beta2 <- param[(d1+1):(d2+d1)]
    #print(length(beta2))
    psi1 <- param[(d1+d2+1):length(param)]
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi1, agg_claims, years)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exposure
    }else{
      mu <- (nu1*(1+se1$spatial_effect))*exposure
    }
    phi <- nu2
    
    
    nu <- claims-0.5
    delta <- sqrt(phi^2 + 2*mu)
    const <- (nu/(phi*delta) - besselK(phi*delta,nu+1,T)/besselK(phi*delta,nu,T))
    
    fmu <- as.matrix(claims/mu - nu/delta^2 + const*phi/delta )
    fphi <- as.matrix(1/phi + 2*phi + nu*(1/phi - phi/delta^2) +  const*(delta + phi^2/delta))
    
    bool_mu <- !(is.nan(fmu) | is.na(fmu) | is.infinite(fmu))
    bool_phi <- !(is.nan(fmu) | is.na(fmu) | is.infinite(fmu))
    
    #beta 1 deriv
    
    if(additive){
      g1 <- t( fmu[bool_mu] * nu1[bool_mu]*exposure[bool_mu]) %*% X1[bool_mu,,drop = F]
    }else{
      g1 <- t( fmu[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu])*exposure[bool_mu]) %*% X1[bool_mu,,drop = F]
      
    }
    
    
    # spatial derive
    if(additive){
      g22 <- as.numeric(fmu[bool_mu])*se1$agg_effect[bool_mu]*exposure[bool_mu]
    }else{
      g22 <- as.numeric(fmu[bool_mu])*se1$agg_effect[bool_mu]*exposure[bool_mu]*nu1[bool_mu]
      
    }
    
    g22 <- tapply(g22,locs[bool_mu],sum)
    
    # beta2 deriv
    g3 <- t( fphi[bool_phi]) %*% X2[bool_phi,,drop = F]
    
    grad <- -c(g1, g3, g22)
    return(grad)
    
    
    
  }
  
  
  
  out <- optim(par = c(beta1, beta2, psi),
               fn = objective_pig_add_psi,
               gr = objective_pig_add_psi_deriv, 
               X1 = X1, X2 = X2, 
               claims = claims, 
               p = p,
               agg_claims = agg_claims,
               locs = locs,
               years=years,
               w= w,
               additive = additive,
               exposure = exposure,
               lower = c(lower_beta1, c(-0.1, rep(0, ncol(X2) - 1)), lower_psi ),
               upper = c(upper_beta1, c(3, rep(0, ncol(X2) - 1)), upper_psi ),
               method = 'L-BFGS-B')
  
  print("finalizing")
  beta1 <- out$par[1:ncol(X1)]
  beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ncol(X2))]
  psi <- out$par[(ncol(X1)+ncol(X2)+1):length(out$par)]
  print(beta2)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exposure
  }else{
    mu <- (nu1*(1+se1$spatial_effect))*exposure
  }
  phi <- nu2
  
  return(list(beta1 = beta1, beta2 = beta2, mu = mu, phi = phi, out = out, psi = psi))
  
  
  
}

PIG_ML_a <- function(claims, X1,X2, additive, agg_claims, years, locs, exposure = rep(1, length(claims)),  max_itr = 1000, z= "", tol = 1e-3, lambda = 0){
  
  p <- length(unique(locs))
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta1 <- beta1 - abs(beta1)
  lower_a <- rep(1e-8, p*(p+1)/2)
  a <- rep(1e-5, p*(p+1)/2)
  A <- get_W_from_array(a,p)
  
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exposure
  }else{
    mu <- (nu1*(1+se1$spatial_effect))*exposure
  }
  phi <- nu2
  
  objective_pig_a <- function(param, X1, X2, p, claims, agg_claims, locs, years,additive, exposure, lambda){
    
    d1 <- ncol(X1)
    d2 <- ncol(X2)
    
    beta1 <- param[1:d1]
    #print(length(beta1))
    beta2 <- param[(d1+1):(d2+d1)]
    #print(length(beta2))
    a1 <- param[(d1+d2+1):length(param)]
    #print(length(a1))
    A1 <- get_W_from_array(a1, p)
    
    
    se1 <- get_spatial_aggregate_new(locs, A1, 1, agg_claims, years)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exposure
    }else{
      mu <- (nu1*(1+se1$spatial_effect))*exposure
    }
    phi <- nu2
    
    
    return(-sum(log_dPIG(claims, mu, phi), na.rm = T) )
    
    
    
  }
  
  objective_pig_a_deriv <- function(param, X1, X2, p, claims, agg_claims, locs, years, additive, exposure, lambda){
    
    d1 <- ncol(X1)
    d2 <- ncol(X2)
    
    
    beta1 <- param[1:d1]
    #print(length(beta1))
    beta2 <- param[(d1+1):(d2+d1)]
    #print(length(beta2))
    a1 <- param[(d1+d2+1):length(param)]
    #print(length(a1))
    A1 <- get_W_from_array(a1, p)
    
    
    se1 <- get_spatial_aggregate_new(locs, A1, 1, agg_claims, years)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exposure
    }else{
      mu <- (nu1*(1+se1$spatial_effect))*exposure
    }
    phi <- nu2
    
    
    nu <- claims-0.5
    delta <- sqrt(phi^2 + 2*mu)
    const <- (nu/(phi*delta) - besselK(phi*delta,nu+1,T)/besselK(phi*delta,nu,T))
    
    fmu <- as.matrix(claims/mu - nu/delta^2 + const*phi/delta )
    fphi <- as.matrix(1/phi + 2*phi + nu*(1/phi - phi/delta^2) +  const*(delta + phi^2/delta))
    
    bool_mu <- !(is.nan(fmu) | is.na(fmu) | is.infinite(fmu))
    bool_phi <- !(is.nan(fmu) | is.na(fmu) | is.infinite(fmu))
    
    
    #beta 1 deriv
    
    if(additive){
      g1 <- t( fmu[bool_mu] * nu1[bool_mu]*exposure[bool_mu]) %*% X1[bool_mu,,drop = F]
    }else{
      g1 <- t( fmu[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu])*exposure[bool_mu]) %*% X1[bool_mu,,drop = F]
      
    }
    
    
    
    
    # spatial derive
    if(additive){
      g2 <- as.numeric(fmu[bool_mu])*t(agg_claims[, years[bool_mu]])*exposure[bool_mu]
    }else{
      g2 <- as.numeric(fmu[bool_mu])*t(agg_claims[, years[bool_mu]])*exposure[bool_mu]*as.numeric(nu1[bool_mu])
    }
    
    
    g2 <- by(g2,locs[bool_mu], FUN=colSums)
    G2 = matrix(unlist(g2),nrow = nrow(A1), ncol = ncol(A1), byrow = T)
    g22 = G2[upper.tri(G2, diag = T)]
    diag(G2) = 0  # make sure we do not double count the diagonal when we add
    g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
    
    # beta2 deriv
    g3 <- t( fphi[bool_phi]) %*% X2[bool_phi,,drop = F]
    
    grad <- -c(g1, g3, g22)
    return(grad)
    
    
    
  }
  
  
  
  out <- optim(par = c(beta1, beta2, a),
               fn = objective_pig_a,
               gr = objective_pig_a_deriv, 
               X1 = X1, X2 = X2, 
               claims = claims, 
               p = p,
               agg_claims = agg_claims,
               locs = locs,
               years = years,
               additive = additive,
               exposure = exposure,
               lambda = lambda,
               lower = c(lower_beta1, c(-0.5, rep(0, ncol(X2) - 1)), lower_a ),
               method = 'L-BFGS-B')
  
  
  beta1 <- out$par[1:ncol(X1)]
  beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ncol(X2))]
  a <- out$par[(ncol(X1)+ncol(X2)+1):length(out$par)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exposure
  }else{
    mu <- (nu1*(1+se1$spatial_effect))*exposure
  }
  phi <- nu2
  
  return(list(beta1 = beta1, beta2 = beta2, mu = mu, phi = phi, out = out, a = a, A =A))
  
  
  
}



# PLN ----
# functions to perform integrations for pln
dlno<-function(z,phi) {
  (1/(sqrt(2*pi)*phi*z) ) *exp(  -((log(z)+((phi^2)/2))^2)  /(2*(phi^2)) )
}


jpdf<-function(y,mu,phi){
  fun=function(z) {
    dpois(y, mu*z)*dlno(z,phi)
  }
  
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  return(out)
}

jpdf <- Vectorize(jpdf, vectorize.args = c("y", "mu", "phi"))

pm2<-function(y,mu,phi, numer){
  fun=function(z) {
    (log(z)^2)* dpois(y, mu*z)* dlno(z,phi)
  }
  
  
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/numer)
}

pm2 <- Vectorize(pm2, vectorize.args = c("y", "mu", "phi","numer"))

pm1<-function(y,mu,phi, numer ){
  fun=function(z) {
    z * dpois(y, mu*z)* dlno(z,phi)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/numer)
}
pm1 <- Vectorize(pm1, vectorize.args = c("y", "mu", "phi","numer"))




log_dpln <- function(y,mu,phi){
  fun=function(z) {
    dpois(y, mu*z)* dlno(z,phi)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  if(!is.na(out)){
    if(out <= 1e-8){
      out <- 1e-8
    } 
  }
  return(log(out))
  
}
log_dpln <- Vectorize(log_dpln, vectorize.args = c("y", "mu", "phi"))



Q_PLN_beta_2 <- function(param, claims, exposure, X2, elogthetasq){
  
  elogthetasq = as.numeric(elogthetasq)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2
  
  
  ll = sum(-log(phi) - elogthetasq/(2*phi^2)  - phi^2/8)
  return(-ll)
  
}

Q_PLN_beta_2_deriv <- function(param, claims, exposure, X2, elogthetasq){
  
  elogthetasq = as.numeric(elogthetasq)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2
  
  der1_logdig <- -1/phi + elogthetasq/(phi ^3) - 2*phi/8
  
  
  g3 <- t( der1_logdig*nu2 ) %*% X2
  
  return(-g3)
}


PLN_EM <- function(claims, X1,X2, exposure = rep(1, length(claims)),  max_itr = 1000, tol = 1e-3){
  
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0.05, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- rep(-5, ncol(X1))
  upper_beta2 <- rep(0.5, ncol(X1))
  print(beta1)
  
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  ll_tmp <-log_dpln(claims,mu,phi)
  ll <-  sum(log_dpln(claims,mu,phi), na.rm = TRUE)
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    print(itr)
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    
    # E-step
    numer <- jpdf(claims, mu,phi)
    elogthetasq <- pm2(claims, mu,phi,numer)
    etheta <- pm1(claims, mu,phi,numer)
    
    
    
    
    bool <- !(is.infinite(etheta) | is.infinite(elogthetasq) | is.nan(etheta) | is.nan(elogthetasq) | is.na(etheta) | is.na(elogthetasq))
    
    etheta_tmp <- etheta[bool]
    elogthetasq_tmp <- elogthetasq[bool]
    X1_tmp <- X1[bool,,drop = T]
    X2_tmp <- X2[bool,,drop = T]
    exposure_tmp <- exposure[bool]
    claims_tmp <- claims[bool]
    
    
    # M-step
    out_mean <- optim(c(beta1),
                      Q_P,
                      gr = Q_P_deriv,
                      claims = claims_tmp, 
                      exposure = exposure_tmp, 
                      X1 = X1_tmp, 
                      etheta = etheta_tmp, 
                      method = 'L-BFGS',
                      control = list(maxit = 2),
                      hessian = calc_hessian)
    print(out_mean$message)
    
    beta1 <- out_mean$par
    H_mean_par <- out_mean$hessian
    #print(beta2)
    out_beta2_deriv_em = optim(beta2,
                               Q_PLN_beta_2,
                               gr = Q_PLN_beta_2_deriv,
                               claims = claims_tmp, 
                               exposure =exposure_tmp, 
                               X2 = X2_tmp,
                               elogthetasq = elogthetasq_tmp,
                               method = 'L-BFGS-B',
                               lower = lower_beta2,
                               upper = upper_beta2,
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    print(out_beta2_deriv_em$message)
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    mu <- nu1 *exp(log(exposure))
    phi <- nu2 
    
    tmp_ll <- log_dpln(claims,mu,phi)
    
    new_ll <-  sum(tmp_ll, na.rm = TRUE)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 == 0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',ll))
    #}
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_deriv(beta1, claims_tmp, exposure_tmp, X1_tmp, etheta_tmp))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PLN_beta_2_deriv(beta2, claims_tmp, exposure_tmp, X2_tmp, elogthetasq_tmp))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2, beta1_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, out_beta2 = out_beta2_deriv_em))
  
}


PLN_EM_psi <- function(claims, X1,X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)),  max_itr = 1000, z = "", tol = 1e-3){
  
  p <- length(unique(locs))
  
  # Set initial parameters
  if(additive){
    out <- Poisson_additive_psi(claims, X1, locs, years, agg_claims, w, exposure, max_itr = 100)
  }else{
    out <- Poisson_multiplicative_psi(claims, X1, locs, years, agg_claims, w, exposure, max_itr = 100)
  }
  beta1 <- out$beta1
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  psi <- out$psi
  
  lower_psi <- rep(1e-8, p)
  lower_beta1 <- beta1 - abs(beta1)
  upper_psi <- rep(0.5, p)
  upper_beta1 <- beta1 + abs(beta1)
  
  
  
  lower_beta2 <- rep(-5, ncol(X2)) 
  upper_beta2 <- rep(0.5, ncol(X2)) 
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  print(beta1)
  
  
  
  
  ll_tmp <- log_dpln(claims,mu,phi)
  ll <-  sum(ll_tmp[!is.infinite(ll_tmp)], na.rm = T)
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    if(all(is.character(z))){
      numer <- jpdf(claims, mu,phi)
      elogthetasq <- pm2(claims, mu,phi,numer)
      etheta <- pm1(claims, mu,phi,numer)
    }else{
      etheta <- z
      elogthetasq <- log(z)^2
    }
    
    
    
    
    bool <- !(is.infinite(etheta) | is.infinite(elogthetasq) | is.nan(etheta) | is.nan(elogthetasq)| is.na(etheta) | is.na(elogthetasq))
    
    etheta_tmp <- etheta[bool]
    elogthetasq_tmp <- elogthetasq[bool]
    X1_tmp <- X1[bool,,drop = T]
    X2_tmp <- X2[bool,,drop = T]
    exposure_tmp <- exposure[bool]
    years_tmp <- years[bool]
    claims_tmp <- claims[bool]
    locs_tmp <- locs[bool]
    
    #print(beta1)
    out_mean <- optim(c(beta1, psi),
                      Q_P_psi,
                      gr = Q_P_psi_deriv,
                      locs = locs_tmp, 
                      claims = claims_tmp, 
                      exposure = exposure_tmp, 
                      X1 = X1_tmp, 
                      agg_claims = agg_claims, 
                      years = years_tmp, 
                      etheta = etheta_tmp, 
                      w = w,
                      additive = additive,
                      method = 'L-BFGS-B',
                      lower = c(lower_beta1, lower_psi),
                      upper = c(upper_beta1, upper_psi),
                      control = list(maxit = 2),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par[1:ncol(X1)]
    psi <- out_mean$par[(ncol(X1)+1):length(out_mean$par)]
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PLN_beta_2,
                               gr = Q_PLN_beta_2_deriv,
                               claims = claims_tmp, 
                               exposure =exposure_tmp, 
                               X2 = X2_tmp, 
                               elogthetasq = elogthetasq_tmp,
                               method = 'L-BFGS-B',
                               lower = lower_beta2,
                               upper = upper_beta2,
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(psi)
    
    
    
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
    }
    phi <- nu2 
    
    new_ll_tmp <- log_dpln(claims,mu,phi)
    new_ll <-  sum(new_ll_tmp[!is.infinite(new_ll_tmp)], na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 == 0){
    print(beta1)
    print(psi)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',ll))
    #}
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_psi_deriv(c(beta1, psi), locs, claims, exposure, X1, agg_claims, years, etheta, w, additive))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PLN_beta_2_deriv(beta2, claims, exposure, X2, elogthetasq_tmp))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  
  
  
  return(list(beta1 = beta1, beta2 = beta2, psi = psi, mean_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, se =se1))
  
  
  
  
  
}


PLN_EM_a <- function(claims, X1,X2, locs, years, agg_claims, additive, exposure = rep(1, length(claims)),  max_itr = 1000, lambda = 0, z = "", tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  a  <- rep(1,p*(p+1)/2)*0.00001
  A <- get_W_from_array(a,p)
  print(beta1)
  
  lower_beta1 <- rep(-10, ncol(X1))
  lower_beta2 <- rep(-5, ncol(X2))
  lower_a <- rep(1e-3, p*(p+1)/2)
  
  upper_beta1 <- beta1 +abs(beta1)
  upper_beta2 <-  rep(0.5, ncol(X2))
  upper_a <- rep(0.5, p*(p+1)/2)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% beta1))
  nu2 <- as.numeric(exp(X2 %*% beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll_tmp <- log_dpln(claims,mu,phi)
  ll <-  sum(ll_tmp[!is.infinite(ll_tmp)], na.rm = TRUE)
  print(paste("ll at start = ", ll))
  has_converged <- FALSE
  last_iteration <- FALSE
  itr <- 1
  
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    # E-step
    if(all(is.character(z))){
      numer <- jpdf(claims, mu,phi)
      elogthetasq <- pm2(claims, mu,phi,numer)
      etheta <- pm1(claims, mu,phi,numer)
    }else{
      etheta <- z
      elogthetasq <- log(z)^2
    }
    
    
    bool <- !(is.infinite(etheta) | is.infinite(elogthetasq) | is.nan(etheta) | is.nan(elogthetasq)| is.na(etheta) | is.na(elogthetasq))
    
    etheta_tmp <- etheta[bool]
    elogthetasq_tmp <- elogthetasq[bool]
    X1_tmp <- X1[bool,,drop = T]
    X2_tmp <- X2[bool,,drop = T]
    exposure_tmp <- exposure[bool]
    years_tmp <- years[bool]
    claims_tmp <- claims[bool]
    locs_tmp <- locs[bool]
    
    
    
    out_mean <- optim(par = c(beta1, a),
                      fn = Q_P_a,
                      gr = Q_P_a_deriv,
                      locs = locs_tmp, 
                      claims = claims_tmp, 
                      exposure = exposure_tmp, 
                      X1 = X1_tmp, 
                      agg_claims = agg_claims, 
                      years = years_tmp, 
                      etheta = etheta_tmp, 
                      additive = additive,
                      p = p,
                      lambda = lambda,
                      method = 'L-BFGS-B',
                      lower = c(lower_beta1, lower_a),
                      upper = c(upper_beta1, upper_a),
                      control = list(maxit = 2),
                      hessian = calc_hessian)
    #print(out_mean$message)
    
    beta1 <- out_mean$par[1:ncol(X1)]
    a <- out_mean$par[(ncol(X1)+1):length(out_mean$par)]
    A <- get_W_from_array(a, p)
    H_mean_par <- out_mean$hessian
    
    
    
    out_beta2_deriv_em = optim(beta2,
                               Q_PLN_beta_2,
                               gr = Q_PLN_beta_2_deriv,
                               claims = claims_tmp, 
                               exposure =exposure_tmp, 
                               X2 = X2_tmp, 
                               elogthetasq = elogthetasq_tmp,
                               method = 'L-BFGS-B',
                               lower = lower_beta2,
                               upper = upper_beta2,
                               control = list(maxit = 2),
                               hessian = calc_hessian)
    
    beta2 <- out_beta2_deriv_em$par
    H_beta2 <- out_beta2_deriv_em$hessian
    
    #print(out_beta2_deriv_em$message)
    
    
    
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 *(1+se1$spatial_effect))*exp(log(exposure))
    }
    phi <- nu2 
    
    new_ll_tmp <- log_dpln(claims,mu,phi)
    new_ll <-  sum(new_ll_tmp[!is.infinite(new_ll_tmp)], na.rm = TRUE)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 == 0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',ll))
    #}
    
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    itr <- itr + 1
    
    
  }
  
  
  print("finalising")
  # The 
  mean_s <- as.numeric(Q_P_a_deriv(c(beta1, a), locs, claims, exposure, X1, agg_claims, years, etheta, additive, p, lambda))
  mean_H <- H_mean_par + outer(mean_s, mean_s)
  
  
  beta2_s <- as.numeric(Q_PLN_beta_2_deriv(beta2, claims, exposure, X2, elogthetasq))
  beta2_H <- H_beta2 + outer(beta2_s, beta2_s)
  
  return(list(beta1 = beta1, beta2 = beta2, a = a, A=A, mean_H = mean_H, beta2_H = beta2_H, mu = mu, phi = phi, se = se1))
  
}



# ZIP models -------------

# Poisson ZIP -------


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


ml_zip_psi_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, w){
  
  prop <- param[1]
  beta <- param[2:(ncol(X_mat)+1)]
  psi <- param[(ncol(X_mat)+2):length(param)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dZIP(claims, mu, prop, 1))
  
  
  
  return(-ll)
  
}

ml_zip_psi_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, w){
  
  prop <- param[1]
  beta <- param[2:(ncol(X_mat)+1)]
  psi <- param[(ncol(X_mat)+2):length(param)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  fphi_deriv <- log_dZIP_der(claims, mu, prop, 1)$fphi
  mu_deriv <- log_dZIP_der(claims, mu,prop, 1)$fmu
  
  # prop
  prop_deriv <- sum(fphi_deriv)
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
  }else{
    beta_deriv <- t(mu_deriv * nu1*(1+se1$spatial_effect) * exp(log(exposure)) ) %*% X_mat
  }
  
  
  # psi
  agg_effect <- colSums(((w %*% agg_claims[,years, drop = F])))
  
  if(additive){
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect
    
  }else{
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect *nu1
  }
  
  psi_deriv <- tapply(psi_deriv,locs,sum)
  
  return(-c(prop_deriv, beta_deriv, psi_deriv))
  
}



ml_zip_a_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, lambda, p){
  
  prop <- param[1]
  beta <- param[2:(ncol(X_mat)+1)]
  a <- param[(ncol(X_mat)+2):length(param)]
  A <- get_W_from_array(a,p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dZIP(claims, mu, prop, 1)) - sum(a*lambda)
  
  
  
  return(-ll)
  
}

ml_zip_a_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, lambda, p){
  
  prop <- param[1]
  beta <- param[2:(ncol(X_mat)+1)]
  a <- param[(ncol(X_mat)+2):length(param)]
  A <- get_W_from_array(a,p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  fphi_deriv <- log_dZIP_der(claims, mu, prop, 1)$fphi
  mu_deriv <- log_dZIP_der(claims, mu,prop, 1)$fmu
  
  # prop
  prop_deriv <- sum(fphi_deriv)
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
  }else{
    beta_deriv <- t(mu_deriv * nu1*(1+se1$spatial_effect) * exp(log(exposure)) ) %*% X_mat
  }
  
  
  # A
  if(additive){
    g2 <- mu_deriv*t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <- mu_deriv*t(agg_claims[, years]) * exp(log(exposure))*as.numeric(nu1)
  }
  
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  return(-c(prop_deriv, beta_deriv, g22))
  
}




Q_ZIP_beta <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, psi, w, prop, additive){
  
  beta <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dZIP(claims, mu, prop, etheta))
  
  
  
  return(-ll)
  
}


Q_ZIP_beta_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, psi, w, prop, additive){
  
  beta <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  mu_deriv <- log_dZIP_der(claims, mu,prop, etheta)$fmu
  
  
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
  }else{
    beta_deriv <- t(mu_deriv * nu1*se1$spatial_effect * exp(log(exposure)) ) %*% X_mat
  }
  
  
  
  return(-beta_deriv)
  
  
}


Q_ZIP_prop <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, psi, w, beta, additive){
  
  prop <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dZIP(claims, mu, prop, etheta))
  
  
  
  return(-ll)
  
}


Q_ZIP_prop_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, psi, w, beta, additive){
  
  prop <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  fphi_deriv <- log_dZIP_der(claims, mu, prop, etheta)$fphi
  
  
  return(-sum(fphi_deriv))
  
  
}



Q_ZIP_psi <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, beta, w, prop, additive){
  
  psi <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dZIP(claims, mu, prop, etheta))
  
  
  
  return(-ll)
  
}


Q_ZIP_psi_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, beta, w, prop, additive){
  
  psi <- param
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  agg_effect <- colSums(((w %*% agg_claims[,years, drop = F])))
  mu_deriv <- log_dZIP_der(claims, mu,prop, etheta)$fmu
  
  
  if(additive){
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect
    
  }else{
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect *nu1
  }
  
  psi_deriv <- tapply(psi_deriv,locs,sum)
  
  return(-psi_deriv)
  
  
}



Q_ZIP_a <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, beta, prop, additive, lambda = 0){
  
  a1 <- param
  A1 <- get_W_from_array(a1,p)
  
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, A1, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dZIP(claims, mu, prop, etheta)) - sum(lambda*a1)
  
  
  
  return(-ll)
  
}


Q_ZIP_a_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, beta,  prop, additive, lambda = 0){
  
  a1 <- param
  A1 <- get_W_from_array(a1,p)
  
  etheta <- as.numeric(etheta)
  
  se1 <- get_spatial_aggregate_new(locs, A1, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% as.matrix(beta)))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  mu_deriv <- log_dZIP_der(claims, mu,prop, etheta)$fmu
  
  
  if(additive){
    g2 <- mu_deriv*t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <- mu_deriv*t(agg_claims[, years]) * exp(log(exposure))*as.numeric(nu1)
  }
  
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A1), ncol = ncol(A1), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  return(-g22)
  
  
}



ZIP_psi <- function(claims, X1, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), max_itr = 1000, z = ""){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  psi <- rep(1,p)*0.001
  prop <- 0.5
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  
  lik <- sum(log_dZIP(claims, mu,prop, 1))
  print(lik)
  for(i in 1:max_itr){
    
    # update prop
    out_prop = optim(prop,
                     fn = Q_ZIP_prop,
                     gr = Q_ZIP_prop_deriv,
                     locs = locs,
                     claims = claims, 
                     exposure = exposure, 
                     X_mat = X1, 
                     agg_claims = agg_claims, 
                     years = years,
                     etheta = 1, 
                     psi = psi, 
                     w = w,
                     beta = beta, 
                     additive = additive,
                     method = 'L-BFGS-B',
                     control = list(maxit = 10),
                     lower = c(1e-3),
                     upper = c(0.99),
                     hessian = T)
    
    prop <- out_prop$par
    # update beta
    out_beta = optim(beta,
                     Q_ZIP_beta,
                     gr = Q_ZIP_beta_deriv,
                     locs = locs,
                     claims = claims, 
                     exposure = exposure, 
                     X_mat = X1, 
                     agg_claims = agg_claims, 
                     years = years,
                     etheta = 1, 
                     psi = psi, 
                     w = w,
                     prop = prop, 
                     additive = additive,
                     method = 'L-BFGS-B',
                     control = list(maxit = 10),
                     hessian = T)
    beta <- out_beta$par
    # update psi
    out_psi = optim(psi,
                    Q_ZIP_psi,
                    gr = Q_ZIP_psi_deriv,
                    locs = locs,
                    claims = claims, 
                    exposure = exposure, 
                    X_mat = X1, 
                    agg_claims = agg_claims, 
                    years = years,
                    etheta = 1, 
                    prop = prop, 
                    w = w,
                    beta = beta, 
                    additive = additive,
                    method = 'L-BFGS-B',
                    control = list(maxit = 10),
                    lower = rep(1e-4, p),
                    hessian = T)
    
    psi <- out_psi$par
    
    
    # finalize iteration
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    new_lik <- sum(log_dZIP(claims, mu,prop, 1))
    print(new_lik)
    
    
    lik <- new_lik
    
  }
  
  return(list(beta = beta, psi = psi, prop = prop, mu = mu))
}


ZIP_a <- function(claims, X1, locs, years, agg_claims, additive, exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0, z = ""){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  if(additive){
    a1 <- rep(0.1, p*(p+1)/2)
    A1 <-  get_W_from_array(a1, p)
  }else{
    A1 <- diag(p)
    a1 <- A1[upper.tri(A1,diag = TRUE)]
  }
  
  prop <- 0.7
  
  se1 <- get_spatial_aggregate_new(locs, A1, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * se1$spatial_effect)*exp(log(exposure))
  }
  
  
  
  lik <- sum(log_dZIP(claims, mu,prop, 1))
  print(lik)
  for(i in 1:max_itr){
    
    # update prop
    out_prop = optim(prop,
                     fn = Q_ZIP_prop,
                     gr = Q_ZIP_prop_deriv,
                     locs = locs,
                     claims = claims, 
                     exposure = exposure, 
                     X_mat = X1, 
                     agg_claims = agg_claims, 
                     years = years,
                     etheta = 1, 
                     psi = 1, 
                     w = A1,
                     beta = beta, 
                     additive = additive,
                     method = 'L-BFGS-B',
                     control = list(maxit = 10),
                     lower = c(1e-3),
                     upper = c(0.99),
                     hessian = T)
    
    prop <- out_prop$par
    # update beta
    out_beta = optim(beta,
                     Q_ZIP_beta,
                     gr = Q_ZIP_beta_deriv,
                     locs = locs,
                     claims = claims, 
                     exposure = exposure, 
                     X_mat = X1, 
                     agg_claims = agg_claims, 
                     years = years,
                     etheta = 1, 
                     psi = 1, 
                     w = A1,
                     prop = prop, 
                     additive = additive,
                     method = 'L-BFGS-B',
                     control = list(maxit = 100),
                     hessian = T)
    beta <- out_beta$par
    # update a
    lower <- diag(p)[upper.tri(diag(p), T)]
    lower[lower == 0] <- 1e-8
    out_A = optim(a1,
                  Q_ZIP_a,
                  gr = Q_ZIP_a_deriv,
                  locs = locs,
                  claims = claims, 
                  exposure = exposure, 
                  X_mat = X1, 
                  agg_claims = agg_claims, 
                  years = years,
                  etheta = 1, 
                  prop = prop, 
                  lambda = lambda,
                  beta = beta, 
                  additive = additive,
                  method = 'L-BFGS-B',
                  control = list(maxit = 100),
                  lower = lower,
                  hessian = T)
    
    a1 <- out_A$par
    A1 <- get_W_from_array(a1,p) 
    
    
    # finalize iteration
    se1 <- get_spatial_aggregate_new(locs, A1, 1, agg_claims, years)
    nu1 <- exp(X1 %*% as.matrix(beta))
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * se1$spatial_effect)*exp(log(exposure))
    }
    
    new_lik <- sum(log_dZIP(claims, mu,prop, 1))
    print(new_lik)
    
    
    lik <- new_lik
    
  }
  
  return(list(beta = beta, a = a1, prop = prop, mu = mu, A = A1))
}



ZIP_all_in_one <- function(claims, X1,  exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  ml_zip_obj <- function(param, claims, exposure, X_mat){
    
    prop <- param[1]
    beta <- param[2:(ncol(X_mat)+1)]
    
    
    nu1 <- as.numeric(exp(X_mat %*% beta))
    
    
    mu <- nu1*exp(log(exposure))
    
    
    ll <- sum(log_dZIP(claims, mu, prop, 1))
    
    
    
    return(-ll)
    
  }
  
  ml_zip_obj_deriv <- function(param,  claims, exposure, X_mat){
    
    prop <- param[1]
    beta <- param[2:(ncol(X_mat)+1)]
    
    nu1 <- as.numeric(exp(X_mat %*% beta))
    
    
    mu <- nu1*exp(log(exposure))
    
    fphi_deriv <- log_dZIP_der(claims, mu, prop, 1)$fphi
    mu_deriv <- log_dZIP_der(claims, mu,prop, 1)$fmu
    
    # prop
    prop_deriv <- sum(fphi_deriv)
    
    # beta
    
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
    
    return(-c(prop_deriv, beta_deriv))
    
  }
  
  
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  prop <- 0.5
  
  
  out = optim(c(prop, beta),
              ml_zip_obj,
              gr = ml_zip_obj_deriv,
              claims = claims, 
              exposure = exposure, 
              X_mat = X1, 
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(1e-3,rep(-20, ncol(X1))),
              upper = c(0.99,rep(Inf, ncol(X1))),
              hessian = T)
  
  
  prop <- out$par[1]
  beta <- out$par[2:(ncol(X1)+1)]
  
  
  
  nu1 <- exp(X1 %*% as.matrix(beta))
  
  mu <- (nu1 )*exp(log(exposure))
  
  
  
  return(list(beta1 = beta, prop = prop, mu = mu, out=out))
}



ZIP_psi_all_in_one <- function(claims, X1, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1, family = 'poisson')$coefficients
  
  psi <- rep(1,p)*0.001
  
  prop <- 0.5
  
  
  out = optim(c(prop, beta, psi),
              ml_zip_psi_obj,
              gr = ml_zip_psi_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat = X1, 
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              w = w,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(1e-3,rep(-20, ncol(X1)) ,rep(1e-4, p)),
              upper = c(0.99,rep(Inf, ncol(X1)), rep(1, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta <- out$par[2:(ncol(X1)+1)]
  psi <- out$par[(ncol(X1)+2):length(out$par)]
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, psi = psi, prop = prop, mu = mu, out=out))
}


ZIP_a_all_in_one <- function(claims, X1, locs, years, agg_claims, additive, exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  
  a <- rep(0.1, p*(p+1)/2)
  A <-  get_W_from_array(a, p)
  
  
  prop <- 0.7
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  lower <- rep(1e-8, p*(p+1)/2) 
  
  
  
  out = optim(par = c(prop, beta, a),
              fn = ml_zip_a_obj,
              gr = ml_zip_a_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat = X1, 
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              lambda = lambda,
              p = p,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(1e-3,rep(-20, ncol(X1)) ,lower),
              upper = c(0.99,rep(Inf, ncol(X1)), rep(1, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta <- out$par[2:(ncol(X1)+1)]
  a <- out$par[(ncol(X1)+2):length(out$par)]
  A <- get_W_from_array(a,p) 
  
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, A = A,a=a, prop = prop, mu = mu, out=out))
}

# ZINP -----

log_zinb <- function(x,mu,phi, prop){
  lab <- (x == 0)*1
  #obj_nb <- (gamma(phi + x)/(gamma(x+1)*gamma(phi)))*(phi^phi)*(mu^x)/((mu+phi)^(phi+x))
  
  #obj_nb[is.nan(obj_nb)] <- 1e-8
  #obj_nb[is.infinite(obj_nb)] <- 1e-8
  
  obj_nb <- exp(dnbinom(x, size = phi, mu = mu, log = TRUE))+1e-8
  
  ll <- lab*log(prop + (1-prop)*obj_nb) + (1-lab)*log(((1-prop)*obj_nb))
  #print(min(ll))
  
  
  return(ll)
}

log_dnb_der <- function(x,mu,phi){
  
  fmu <- (x/mu)- (phi+x)/(mu+phi)
  fphi <- 1 + log(phi) - (x + phi)/(mu + phi) - log(mu + phi) + 
    digamma(x + phi) - digamma(phi)
  
  return(list(fmu= fmu,fphi = fphi))
}


log_zinb_der <- function(x,mu,phi, prop){
  lab <- (x == 0)*1
  
  # derivative of log density db
  fmu_log <- (x/mu)- (phi+x)/(mu+phi)
  fphi_log <- 1 + log(phi) - (x + phi)/(mu + phi) - log(mu + phi) + 
    digamma(x + phi) - digamma(phi)
  
  # derivative of density db at x = 0
  fmu <- -(phi^phi)*phi*(mu+phi)^(-phi-1)
  fphi <- (phi^phi)*(phi + mu)^(-phi-1)*((phi+mu)*log(phi) - (phi+mu)*log(phi+mu) + mu)
  
  
  obj_nb <- exp(dnbinom(x, size = phi, mu = mu, log = TRUE))+1e-8
  
  
  const <- 1/(prop + (1-prop)*obj_nb)
  
  
  fmu_return <- lab*const*(1-prop)*fmu + (1-lab)*fmu_log
  fphi_return <- lab*const*(1-prop)*fphi + (1-lab)*fphi_log
  fprop <- lab*const*(1 - obj_nb) - (1-lab)/(1-prop)
  
  
  return(list(fmu = fmu_return, fphi = fphi_return, fprop = fprop))
}

ml_zinb_obj <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  mu <- nu1*exposure
  
  
  
  ll <- sum(log_zinb(claims, mu, phi, prop))
  
  return(-ll)
  
}

ml_zinb_obj_deriv <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  mu <- nu1*exposure
  
  
  derivs <- log_zinb_der(claims, mu, phi, prop)
  
  
  bool_fphi <- !(is.infinite(derivs$fphi) | is.na(derivs$fphi))
  derivs$fphi <- derivs$fphi[bool_fphi]
  
  bool_fmu <- !(is.infinite(derivs$fmu) | is.na(derivs$fmu))
  derivs$fmu <- derivs$fmu[bool_fmu]
  
  bool_fprop <- !(is.infinite(derivs$fprop) | is.na(derivs$fprop))
  derivs$fprop <- derivs$fprop[bool_fprop]
  
  fphi_deriv <- derivs$fphi
  mu_deriv <- derivs$fmu
  prop_deriv <- derivs$fprop
  
  # prop
  prop_deriv <- sum(prop_deriv[bool_fprop])
  
  # phi
  phi_deriv <- t( fphi_deriv*nu2[bool_fphi] ) %*% X_mat2[bool_fphi,,drop = F]
  
  # beta
  
  beta_deriv <- t(mu_deriv * nu1[bool_fmu] * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = F]
  
  
  
  return(-c(prop_deriv, phi_deriv, beta_deriv))
  
}

ml_zinb_psi_obj <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, w){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  psi <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll_tmp <- log_zinb(claims, mu, phi, prop)
  
  ll <- sum(ll_tmp, na.rm = T)
  
  
  return(-ll)
}

ml_zinb_psi_obj_deriv <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, w){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  psi <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  derivs <- log_zinb_der(claims, mu, phi, prop)
  
  fphi_deriv <- derivs$fphi
  mu_deriv <- derivs$fmu
  prop_deriv <- derivs$fprop
  
  # prop
  prop_deriv <- sum(prop_deriv)
  
  # phi
  phi_deriv <- t( fphi_deriv*nu2 ) %*% X_mat2
  
  # beta
  
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat1
  }else{
    beta_deriv <- t(mu_deriv * nu1*(1+se1$spatial_effect) * exp(log(exposure)) ) %*% X_mat1
  }
  
  
  if(additive){
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect
    
  }else{
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect *nu1
  }
  #print(bool_mu)
  psi_deriv <- tapply(psi_deriv,locs,sum)
  
  
  
  
  
  
  
  return(-c(prop_deriv, phi_deriv, beta_deriv, psi_deriv))
  
}

ml_zinb_a_obj <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, lambda, p){
  
  p <- length(unique(locs))
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  a <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  A <- get_W_from_array(a, p)
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_zinb(claims, mu, phi, prop)) - sum(a*lambda)
  
  return(-ll)
  
}

ml_zinb_a_obj_deriv <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, lambda, p){
  
  p <- length(unique(locs))
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  a <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  A <- get_W_from_array(a, p)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  derivs <- log_zinb_der(claims, mu, phi, prop)
  
  fphi_deriv <- derivs$fphi
  
  mu_deriv <- derivs$fmu
  
  prop_deriv <- derivs$fprop
  
  # prop
  prop_deriv <- sum(prop_deriv)
  
  # phi
  phi_deriv <- t( fphi_deriv*nu2 ) %*% X_mat2
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat1
  }else{
    beta_deriv <- t(mu_deriv * nu1*(1+se1$spatial_effect) * exp(log(exposure)) ) %*% X_mat1
  }
  
  
  
  # A
  if(additive){
    g2 <- mu_deriv[,1]*t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <- mu_deriv[,1]*t(agg_claims[, years]) * exp(log(exposure))*as.numeric(nu1)
  }
  
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  
  return(-c(prop_deriv, phi_deriv, beta_deriv, g22))
  
}


ZINB_all_in_one <- function(claims, X1, X2, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  start <- ZIP_all_in_one(claims, X1, exposure, max_itr = 100)
  beta <- start$beta
  prop <- start$prop
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  
  lower_beta1 <- beta - 2*abs(beta)
  upper_beta1 <- beta + 2*abs(beta)
  
  
  
  out = optim(c(prop, beta2, beta),
              ml_zinb_obj,
              gr = ml_zinb_obj_deriv,
              claims = claims, 
              exposure = exposure, 
              X_mat1 = X1, 
              X_mat2 = X2,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(0.01, rep(-2, ncol(X2)), lower_beta1),
              upper = c(0.99, rep(5, ncol(X2)), upper_beta1 ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  
  
  
  
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  mu <- nu1*exposure
  
  
  return(list(beta1 = beta,beta2 = beta2, prop = prop, mu = mu, out=out))
}

ZINB_psi_all_in_one <- function(claims, X1, X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  start <- ZIP_psi_all_in_one(claims, X1, locs, years, agg_claims, w, additive, exposure, max_itr = 100)
  beta <- start$beta
  psi <- start$psi
  prop <- start$prop
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  
  lower_beta1 <- beta - 2*abs(beta)
  upper_beta1 <- beta + 2*abs(beta)
  
  print(beta)
  print(beta2)
  print(prop)
  print(psi)
  
  
  
  out = optim(c(prop, beta2, beta, psi),
              ml_zinb_psi_obj,
              gr = ml_zinb_psi_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat1 = X1, 
              X_mat2 = X2,
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              w =w,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(1e-3, rep(-2, ncol(X2)),lower_beta1, rep(1e-8, p)),
              upper = c(0.99,  rep(2, ncol(X2)),upper_beta1, rep(1, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  psi <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta,beta2 = beta2, psi = psi, prop = prop, mu = mu, out=out))
}

ZINB_a_all_in_one <- function(claims, X1, X2, locs, years, agg_claims, additive, exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  start <- ZIP_a_all_in_one(claims, X1, locs, years, agg_claims, additive, exposure, max_itr = 100, lambda = lambda)
  beta <- start$beta
  a <- start$a
  prop <- start$prop
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  
  
  lower_beta1 <- beta - 2*abs(beta)
  upper_beta1 <- beta + 2*abs(beta)
  
  out = optim(par = c(prop, beta2, beta, a),
              fn = ml_zinb_a_obj,
              gr = ml_zinb_a_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat1 = X1, 
              X_mat2 = X2,
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              lambda = lambda,
              p = p,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(1e-3, rep(-2, ncol(X2)), lower_beta1 ,rep(1e-8, p*(p+1))),
              upper = c(0.99, rep(4, ncol(X2)), upper_beta1, rep(0.5, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  a <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
  A <- get_W_from_array(a, p)
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, beta2 = beta2, a = a, A=A, prop = prop, mu = mu, out=out))
}

# ZIPIG ----------------

# Poisson inverse Gaussian
log_dPIG <- function(x,mu,phi){
  nu <- x-0.5
  delta <- sqrt(phi^2 + 2*mu)
  res <-  0.5*log(2/pi) + log(phi) + phi^2 + 
    log(besselK(phi*delta,nu,expon.scaled = T))-phi*delta + 
    nu*log(phi/delta) + x*log(mu) - lgamma(x+1)
  
  if(any(is.infinite(res))){
    max_val <- max(res[!is.infinite(res)])
    res[is.infinite(res)] <- max_val*100 
  }
  return(res)
}

log_dPIG_der <- function(x,mu,phi){
  nu <- x-0.5
  delta <- sqrt(phi^2 + 2*mu)
  const <- (nu/(phi*delta) - besselK(phi*delta,nu+1,T)/besselK(phi*delta,nu,T))
  
  if(any(is.infinite(const))){
    max_val <- max(const[!is.infinite(const)])
    const[is.infinite(const)] <- max_val*100 
  }
  
  fmu <- x/mu - nu/delta^2 + const*phi/delta 
  fphi <- 1/phi + 2*phi + nu*(1/phi - phi/delta^2) +  const*(delta + phi^2/delta)
  
  return(list(fmu = fmu, fphi = fphi))
}



log_zipig <- function(x,mu,phi, prop){
  lab <- (x == 0)*1
  obj_nb <- exp(log_dPIG(x, mu, phi))
  ll_tmp <- lab*log(prop + (1-prop)*obj_nb) + (1-lab)*log(((1-prop)*obj_nb))
  
  
  return(ll_tmp)
}


K_05_deriv <- function(x){
  
  const1 <- -0.5*sqrt(pi/2)*exp(-x)*x^(-1.5)*(2*x+1)
  
  return(const1)
  
}




log_zipig_der <- function(x,mu,phi, prop){
  lab <- (x == 0)*1
  
  # derivative of log density db
  log_deriv <- log_dPIG_der(x, mu, phi)
  fmu_log <- log_deriv$fmu
  fphi_log <- log_deriv$fphi
  
  # derivative of density db at x = 0
  delta <- sqrt(phi ^2 + 2*mu)
  delta_deriv_mu <- (phi ^2 + 2*mu)^(-0.5)
  delta_deriv_phi <- (phi ^2 + 2*mu)^(-0.5)*phi
  
  # mu part
  const_mu <- sqrt(2/pi)*phi*exp(phi^2)*phi^(-0.5)
  
  fmu1 <- const_mu*K_05_deriv(phi*delta)*phi*delta_deriv_mu*delta^(0.5)
  fmu2 <- const_mu*besselK(phi*delta,-0.5)*0.5*delta^(-0.5)*delta_deriv_mu
  fmu <- fmu1 + fmu2
  #print(const_mu)
  #print(fmu2)
  
  # phi part
  const_phi <-  sqrt(2/pi)
  fphi1 <- const_phi*(0.5)*phi^(-0.5)*delta^0.5*besselK(phi*delta,-0.5)*exp(phi^2)
  fphi2 <- const_phi*phi^(0.5)*0.5*delta^(-0.5)*delta_deriv_phi*besselK(phi*delta,-0.5)*exp(phi^2)
  fphi3 <- const_phi*phi^(0.5)*delta^0.5*K_05_deriv(phi*delta)*(delta + phi*delta_deriv_phi)*exp(phi^2)
  fphi4 <- const_phi*phi^(0.5)*delta^0.5*besselK(phi*delta,-0.5)*exp(phi^2)*2*phi
  fphi <- fphi1 + fphi2 + fphi3+fphi4
  
  
  obj_pig <- exp(log_dPIG(x, mu, phi))
  const <- 1/(prop + (1-prop)*obj_pig)
  
  
  fmu_return <- lab*const*(1-prop)*fmu + (1-lab)*fmu_log
  fphi_return <- lab*const*(1-prop)*fphi + (1-lab)*fphi_log
  fprop <- lab*const*(1 - obj_pig) - (1-lab)/(1-prop)
  
  
  return(list(fmu = fmu_return, fphi = fphi_return, fprop = fprop))
}


ml_zipig_obj <- function(param, claims, exposure, X_mat1, X_mat2){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  
  
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  
  mu <- nu1*exp(log(exposure))
  
  ll_tmp <- log_zipig(claims, mu, phi, prop)
  ll_tmp <- ll_tmp[!(is.na(ll_tmp) | is.infinite(ll_tmp)) ]
  ll <- sum(ll_tmp)
  #print(ll)
  
  
  return(-ll)
  
}

ml_zipig_obj_deriv <- function(param, claims, exposure, X_mat1, X_mat2){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  
  
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  
  derivs <- log_zipig_der(claims, mu, phi, prop)
  
  bool_fphi <- !(is.infinite(derivs$fphi) | is.na(derivs$fphi))
  derivs$fphi <- derivs$fphi[bool_fphi]
  
  bool_fmu <- !(is.infinite(derivs$fmu) | is.na(derivs$fmu))
  derivs$fmu <- derivs$fmu[bool_fmu]
  
  bool_fprop <- !(is.infinite(derivs$fprop) | is.na(derivs$fprop))
  derivs$fprop <- derivs$fprop[bool_fprop]
  
  fphi_deriv <- derivs$fphi
  mu_deriv <- derivs$fmu
  prop_deriv <- derivs$fprop
  
  # prop
  prop_deriv <- sum(prop_deriv[bool_fprop])
  
  # phi
  phi_deriv <- t( fphi_deriv*nu2[bool_fphi] ) %*% X_mat2[bool_fphi,,drop = t]
  
  # beta
  beta_deriv <- t(mu_deriv * nu1[bool_fmu] * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = t]
  
  grad <- -c(prop_deriv, phi_deriv, beta_deriv)
  #print(grad)
  return(grad)
  
}





ml_zipig_psi_obj <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, w){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  psi <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  ll_tmp <- log_zipig(claims, mu+1e-8, phi, prop)
  ll_tmp <- ll_tmp[!(is.na(ll_tmp) | is.infinite(ll_tmp) | is.nan(ll_tmp)) ]
  ll <- sum(ll_tmp)
  
  return(-ll)
  
}

ml_zipig_psi_obj_deriv <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, w){
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  psi <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure)) + 1e-8
  }else{
    mu <- (nu1 * se1$spatial_effect)*exp(log(exposure)) + 1e-8
  }
  
  
  
  derivs <- log_zipig_der(claims, mu, phi, prop)
  
  bool_fphi <- !(is.infinite(derivs$fphi) | is.na(derivs$fphi)| is.nan(derivs$fphi))
  derivs$fphi[!bool_fphi] <- 0
  
  bool_fmu <- !(is.infinite(derivs$fmu) | is.na(derivs$fmu)| is.nan(derivs$fmu))
  #derivs$fmu <- derivs$fmu[bool_fmu]
  
  bool_fprop <- !(is.infinite(derivs$fprop) | is.na(derivs$fprop)| is.nan(derivs$fprop))
  derivs$fprop <- derivs$fprop[bool_fprop]
  
  
  fphi_deriv <- derivs$fphi
  mu_deriv <- derivs$fmu
  prop_deriv <- derivs$fprop
  
  
  # prop
  prop_deriv <- sum(prop_deriv[bool_fprop], na.rm = T)
  
  # phi
  phi_deriv <- t( fphi_deriv*nu2 ) %*% X_mat2
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv[bool_fmu] * nu1[bool_fmu] * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = t]
  }else{
    beta_deriv <- t(mu_deriv[bool_fmu] * nu1[bool_fmu] * (1+se1$spatial_effect[bool_fmu]) * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = t]
  }
  
  # psi
  agg_effect <- colSums(((w %*% agg_claims[,years, drop = F])))
  
  if(additive){
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect
    
  }else{
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect *nu1
  }
  
  psi_deriv <- tapply(psi_deriv,locs,sum, na.rm = T)
  
  return(-c(prop_deriv, phi_deriv, beta_deriv, psi_deriv))
  
  
  
  
}



ml_zipig_a_obj <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, lambda, p){
  
  p <- length(unique(locs))
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  a <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll_tmp <- log_zipig(claims, mu, phi, prop)
  ll_tmp <- ll_tmp[!(is.na(ll_tmp) | is.infinite(ll_tmp)) ]
  ll <- sum(ll_tmp) - sum(a*lambda)
  
  
  return(-ll)
  
}

ml_zipig_a_obj_deriv <- function(param, locs, claims, exposure, X_mat1, X_mat2, agg_claims, years, additive, lambda, p){
  
  p <- length(unique(locs))
  
  prop <- param[1]
  log_phi <- param[2:(ncol(X_mat2)+1)]
  beta <- param[(ncol(X_mat2)+2):(ncol(X_mat1)+ ncol(X_mat2)+1)]
  a <- param[(ncol(X_mat1)+ ncol(X_mat2)+2):length(param)]
  A <- get_W_from_array(a, p)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat1 %*% beta))
  nu2 <- exp(X_mat2 %*% as.matrix(log_phi))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  derivs <- log_zipig_der(claims, mu, phi, prop)
  bool_fphi <- !(is.infinite(derivs$fphi) | is.na(derivs$fphi))
  derivs$fphi <- derivs$fphi[bool_fphi]
  bool_fmu <- !(is.infinite(derivs$fmu) | is.na(derivs$fmu))
  derivs$fmu <- derivs$fmu[bool_fmu]
  bool_fprop <- !(is.infinite(derivs$fprop) | is.na(derivs$fprop))
  derivs$fprop <- derivs$fprop[bool_fprop]
  
  fphi_deriv <- derivs$fphi
  mu_deriv <- derivs$fmu
  prop_deriv <- derivs$fprop
  
  # prop
  prop_deriv <- sum(prop_deriv[bool_fprop])
  
  # phi
  phi_deriv <- t( fphi_deriv*nu2[bool_fphi] ) %*% X_mat2[bool_fphi,,drop = t]
  
  # beta
  beta_deriv <- t(mu_deriv * nu1[bool_fmu] * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = t]
  
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv * nu1[bool_fmu] * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = t]
  }else{
    beta_deriv <- t(mu_deriv * nu1[bool_fmu] * (1+se1$spatial_effect[bool_fmu]) * exp(log(exposure[bool_fmu])) ) %*% X_mat1[bool_fmu,,drop = t]
  }
  
  
  
  # A
  if(additive){
    g2 <- mu_deriv*t(agg_claims[, years[bool_fmu]]) * exp(log(exposure[bool_fmu]))
  }else{
    g2 <- mu_deriv*t(agg_claims[, years[bool_fmu]]) * exp(log(exposure[bool_fmu]))*as.numeric(nu1[bool_fmu])
  }
  
  
  g2 <- by(g2,locs[bool_fmu], FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  
  return(-c(prop_deriv, phi_deriv, beta_deriv, g22))
  
}




ZIPIG_all_in_one <- function(claims, X1, X2, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  # Set initial parameters
  out <- ZIP_all_in_one(claims,X1, exposure, max_itr)
  beta <- out$beta
  prop <- out$prop
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  
  print(beta)
  out = optim(c(prop, beta2, beta),
              ml_zipig_obj,
              gr = ml_zipig_obj_deriv,
              claims = claims, 
              exposure = exposure, 
              X_mat1 = X1, 
              X_mat2 = X2,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(0.1,rep(-1, ncol(X2)), rep(-5, ncol(X1))),
              upper = c(0.9, rep(2, ncol(X2)), rep(5, ncol(X1)) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  
  
  
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  
  mu <- nu1*exp(log(exposure))
  
  
  
  
  return(list(beta1 = beta, beta2 = beta2, prop = prop, mu = mu, out=out))
}


ZIPIG_psi_all_in_one <- function(claims, X1, X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  p <- length(unique(locs))
  
  
  
  out <- ZIP_psi_all_in_one(claims, X1, locs, years, agg_claims, w, additive, exposure, 100)
  
  
  # Set initial parameters
  beta <- out$beta
  prop <- out$prop
  psi <- out$psi  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  print(beta)
  print(prop)
  
  
  out = optim(c(prop, beta2, beta, psi),
              ml_zipig_psi_obj,
              gr = ml_zipig_psi_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat1 = X1, 
              X_mat2 = X2,
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              w = w,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(0.1,rep(-0.5, ncol(X2)), rep(-20, ncol(X1)) ,rep(1e-8, p)),
              upper = c(0.9, rep(1, ncol(X2)), rep(Inf, ncol(X1)), rep(1, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  psi <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * se1$spatial_effect)*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta,beta2 = beta2, psi = psi, prop = prop, mu = mu, out=out))
}

ZIPIG_a_all_in_one <- function(claims, X1, X2, locs, years, agg_claims, additive, exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0){
  
  
  p <- length(unique(locs))
  print(p)
  
  # Set initial parameters
  start <- ZIP_a_all_in_one(claims, X1, locs, years, agg_claims, FALSE, exposure, max_itr = 100,lambda = lambda)
  beta <- start$beta
  
  a <- rep(1e-3, p*(p+1)/2)
  A <-  get_W_from_array(a, p)
  
  prop <- start$prop
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  
  print(beta)
  
  lower <- rep(1e-8, p*(p+1)/2) 
  
  lower_beta1 <- beta -abs(beta)
  upper_beta1 <- beta +abs(beta)
  
  
  out = optim(par = c(prop, beta2, beta, a),
              fn = ml_zipig_a_obj,
              gr = ml_zipig_a_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat1 = X1, 
              X_mat2 = X2,
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              method = 'L-BFGS-B',
              lambda = lambda,
              p = p,
              control = list(maxit = max_itr),
              lower = c(1e-3,rep(-0.1, ncol(X2)), lower_beta1 , lower),
              upper = c(0.99, rep(3, ncol(X2)), upper_beta1, rep(1, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  a <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
  A <- get_W_from_array(a, p)
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, beta2 = beta2, a = a, A=A, prop = prop, mu = mu, out=out))
}



#ZIPLN ------

# functions to perform integrations for ZIpln
dlno<-function(z,phi) {
  
  (1/(sqrt(2*pi)*phi*z) ) *exp(  -((log(z)+((phi^2)/2))^2)  /(2*(phi^2)) )
  
}

dzip <- function(y, mu, prop){
  
  lab <- (y == 0)*1
  
  return(lab*(prop + (1-prop)*dpois(y, mu)) + (1-lab)*(1-prop)*dpois(y, mu) )
}


zijpdf<-function(y,mu,phi,prop){
  fun=function(z) {
    dzip(y, mu*z, prop)*dlno(z,phi)
  }
  
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  
  return(out)
}

zijpdf <- Vectorize(zijpdf, vectorize.args = c("y", "mu", "phi", "prop"))

zipm2<-function(y,mu,phi, prop, numer){
  fun=function(z) {
    (log(z)^2)* dzip(y, mu*z, prop)* dlno(z,phi)
  }
  
  
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(x)
}

zipm2 <- Vectorize(zipm2, vectorize.args = c("y", "mu", "phi", "prop", "numer"))

zipm1<-function(y,mu,phi, prop, numer){
  fun=function(z) {
    z * dzip(y, mu*z, prop)* dlno(z,phi)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(x)
}

zipm1 <- Vectorize(zipm1, vectorize.args = c("y", "mu", "phi", "prop","numer"))


log_dzipln <- function(y,mu,phi,prop){
  fun=function(z) {
    dzip(y, mu*z,prop)* dlno(z,phi)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  if(!is.na(out)){
    if(out <= 1e-8){
      out <- 1e-8
    }
  }
  
  return(log(out))
  
}
log_dzipln <- Vectorize(log_dzipln, vectorize.args = c("y", "mu", "phi", "prop"))



Q_zipln_numerator <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev){
  fun=function(z) {
    (log_dZIP(y, mu, prop, z) +  log(dlno(z,phi)))* dzip(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out)
  
}

Q_zipln_numerator <- Vectorize(Q_zipln_numerator, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev"))


Q_zipln <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev, denominator){
  
  numerator <- Q_zipln_numerator(y,mu, phi, prop, mu_prev, phi_prev, prop_prev)
  return(numerator/denominator)
  
}


prop_deriv_pln <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    (log_dZIP_der(y, mu, prop, z)$fphi)* dzip(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
prop_deriv_pln <- Vectorize(prop_deriv_pln, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))

mu_deriv_pln <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    (log_dZIP_der(y, mu, prop, z)$fmu)* dzip(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  
  return(out/denominator)
}
mu_deriv_pln <- Vectorize(mu_deriv_pln, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))

phi_deriv_pln <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    
    (-1/phi + (log(z)^2)/(phi ^3) - (phi)/4)* dzip(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
phi_deriv_pln <- Vectorize(phi_deriv_pln, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))



Q_ZIPLN_no_spatial <- function(param, param_prev,  claims, exposure, X1,X2){
  
  # Previous 
  prop_prev <- param_prev[1]
  beta1_prev <- param_prev[2:(ncol(X1)+1)]
  beta2_prev <- param_prev[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  
  
  # current estimate
  prop <- param[1]
  beta1 <- param[2:(ncol(X1)+1)]
  beta2 <- param[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  
  denominator <- exp(log_dzipln(claims, mu, phi, prop))
  ll <- Q_zipln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  
  ll <- log_dzipln(claims,mu,phi, prop)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_ZIPLN_no_spatial_deriv <- function(param, param_prev,  claims, exposure, X1,X2){
  
  # Previous 
  prop_prev <- param_prev[1]
  beta1_prev <- param_prev[2:(ncol(X1)+1)]
  beta2_prev <- param_prev[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  
  
  # current estimate
  prop <- param[1]
  beta1 <- param[2:(ncol(X1)+1)]
  beta2 <- param[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  denominator <- exp(log_dzipln(claims, mu, phi, prop))
  
  #print("prop deriv")
  # prop deriv 
  
  prop_deriv <- prop_deriv_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_prop <- !(is.na(prop_deriv) | is.infinite(prop_deriv) |is.nan(prop_deriv))
  prop_deriv<- sum(prop_deriv[bool_prop])
  
  
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_pln(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv) |is.nan(mu_deriv))
  
  
  
  
  beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv) |is.nan(phi_deriv))
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop=FALSE]
  
  
  
  return(-c(prop_deriv, beta1_deriv, beta2_deriv))
  
  
}


ZIPLN_no_spatial <- function(claims, X1,X2,exposure = rep(1, length(claims)), 
                             max_itr = 1000, z = "", tol = 1e-3){
  
  
  
  # Set initial parameters and bounds
  out <- ZIP_all_in_one(claims, X1, exposure, max_itr)
  beta1 <- out$beta
  prop <- out$prop
  print(beta1)
  print(prop)
  
  lower_beta1 <- beta1 - abs(beta1)
  upper_beta1 <- beta1 + abs(beta1)
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-5, rep(0, ncol(X2)-1))
  upper_beta2 <- c(0.5, rep(0, ncol(X2)-1))
  
  lower_prop <- 1e-3
  upper_prop <- 0.999
  
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2%*% as.matrix(beta2))
  mu <- nu1*exp(log(exposure))
  phi <- nu2
  
  
  ll <-  sum(log_dzipln(claims,mu,phi, prop))
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    
    out    <- optim(c(prop, beta1, beta2),
                    Q_ZIPLN_no_spatial,
                    gr = Q_ZIPLN_no_spatial_deriv,
                    param_prev = c(prop, beta1, beta2),
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    method = 'L-BFGS-B',
                    lower = c(lower_prop, lower_beta1, lower_beta2),
                    upper = c(upper_prop, upper_beta1, upper_beta2),
                    control = list(maxit = 2, trace = 0), hessian = calc_hessian)
    
    prop <- out$par[1]
    beta1 <- out$par[2:(ncol(X1)+1)]
    beta2 <- out$par[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
    
    
    # Create new estimates
    nu1 <- as.numeric(exp(X1 %*% beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    mu <- nu1*exp(log(exposure))
    
    # finalize
    new_ll <-  sum(log_dzipln(claims,mu,phi, prop))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(beta1)
    print(beta2)
    print(paste('At iteration', itr, 'log likelihood =', new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_ZIPLN_no_spatial_deriv(out$par, out$par,  claims, exposure, X1, X2))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}



Q_ZIPLN_psi <- function(param, param_prev, locs, claims, exposure, X1,X2, agg_claims, years,  additive, w){
  
  # Previous 
  prop_prev <- param_prev[1]
  beta1_prev <- param_prev[2:(ncol(X1)+1)]
  beta2_prev <- param_prev[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+2):length(param_prev)]
  
  
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years )
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  denominator <- exp(log_dzipln(claims, mu_prev, phi_prev, prop_prev))
  
  
  # current estimate
  prop <- param[1]
  beta1 <- param[2:(ncol(X1)+1)]
  beta2 <- param[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  psi <- param[(ncol(X1)+ ncol(X2)+2):length(param)]
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years )
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  
  ll <- Q_zipln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_ZIPLN_psi_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, agg_claims, years,  additive, w){
  
  # Previous 
  prop_prev <- param_prev[1]
  beta1_prev <- param_prev[2:(ncol(X1)+1)]
  beta2_prev <- param_prev[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+2):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years )
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  denominator <- exp(log_dzipln(claims, mu_prev, phi_prev, prop_prev))
  
  
  # current estimate
  prop <- param[1]
  beta1 <- param[2:(ncol(X1)+1)]
  beta2 <- param[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  psi <- param[(ncol(X1)+ ncol(X2)+2):length(param)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years )
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  #print("prop deriv")
  # prop deriv 
  
  prop_deriv <- prop_deriv_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_prop <- !(is.na(prop_deriv) | is.infinite(prop_deriv) |is.nan(prop_deriv))
  prop_deriv<- sum(prop_deriv[bool_prop])
  
  
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_pln(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv) |is.nan(mu_deriv))
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv) |is.nan(phi_deriv))
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop=FALSE]
  
  
  
  
  # psi deriv
  # print("psi 1 deriv")
  if(additive){
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu]
    
  }else{
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu] *nu1[bool_mu]
  }
  
  psi_deriv <- tapply(psi_deriv,locs[bool_mu],sum)
  
  
  
  
  
  return(-c(prop_deriv, beta1_deriv, beta2_deriv, psi_deriv))
  
  
}


ZIPLN_psi <- function(claims, X1,X2, locs, years, agg_claims, w, exposure = rep(1, length(claims)), additive, 
                      max_itr = 1000, z = "", tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  
  # Set initial parameters and bounds
  
  start <- ZIP_psi_all_in_one(claims, X1, locs, years, agg_claims, w, additive, exposure = exposure, max_itr = 100)
  
  beta1 <- start$beta
  prop <- start$prop
  psi <- start$psi
  print(beta1)
  print(prop)
  
  lower_beta1 <- beta1 -2*abs(beta1)
  upper_beta1 <- beta1 +2*abs(beta1)
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-5, rep(5, ncol(X2)-1))
  upper_beta2 <- c(0.5, rep(5, ncol(X2)-1))
  
  
  lower_psi <- rep(1e-4, p)
  upper_psi <- rep(0.5, p)
  
  
  lower_prop <- 1e-3
  upper_prop <- 0.999
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2%*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  phi <- nu2
  
  
  ll <-  sum(log_dzipln(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    out    <- optim(c(prop, beta1, beta2, psi),
                    Q_ZIPLN_psi,
                    gr = Q_ZIPLN_psi_deriv,
                    param_prev = c(prop, beta1, beta2, psi),
                    locs = locs,
                    claims = claims,
                    exposure =exposure,
                    X1 = X1,
                    X2 = X2,
                    agg_claims = agg_claims,
                    years = years,
                    additive = additive,
                    w = w,
                    method = 'L-BFGS-B',
                    lower = c(lower_prop, lower_beta1, lower_beta2, lower_psi),
                    upper = c(upper_prop, upper_beta1, upper_beta2, upper_psi),
                    control = list(maxit = 1, trace = 0), hessian = calc_hessian)
    
    prop <- out$par[1]
    beta1 <- out$par[2:(ncol(X1)+1)]
    beta2 <- out$par[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
    psi <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
    
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dzipln(claims,mu,phi, prop))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_ZIPLN_psi_deriv(out$par, out$par, locs, claims, exposure, X1, X2, agg_claims, years,  additive, w))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, psi = psi, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}




Q_ZIPLN_a <- function(param, param_prev, locs, claims, exposure, X1,X2, agg_claims, years,  additive, p, lambda){
  
  
  
  # Previous 
  prop_prev <- param_prev[1]
  beta1_prev <- param_prev[2:(ncol(X1)+1)]
  beta2_prev <- param_prev[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+2):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years )
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  denominator <- exp(log_dzipln(claims, mu_prev, phi_prev, prop_prev))
  
  
  # current estimate
  prop <- param[1]
  beta1 <- param[2:(ncol(X1)+1)]
  beta2 <- param[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  a <- param[(ncol(X1)+ ncol(X2)+2):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years )
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll <- Q_zipln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator) - sum(a*lambda)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_ZIPLN_a_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, agg_claims, years,  additive, p, lambda){
  
  # Previous 
  prop_prev <- param_prev[1]
  beta1_prev <- param_prev[2:(ncol(X1)+1)]
  beta2_prev <- param_prev[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+2):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years )
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  denominator <- exp(log_dzipln(claims, mu_prev, phi_prev, prop_prev))
  
  
  # current estimate
  prop <- param[1]
  beta1 <- param[2:(ncol(X1)+1)]
  beta2 <- param[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
  a <- param[(ncol(X1)+ ncol(X2)+2):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years )
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  #print("prop deriv")
  # prop deriv 
  
  prop_deriv <- prop_deriv_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_prop <- !(is.na(prop_deriv) | is.infinite(prop_deriv) |is.nan(prop_deriv))
  prop_deriv<- sum(prop_deriv[bool_prop])
  
  
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_pln(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv) |is.nan(mu_deriv))
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv) |is.nan(phi_deriv))
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop=FALSE]
  
  
  
  
  # A
  if(additive){
    g2 <- as.numeric(mu_deriv[bool_mu])*t(agg_claims[, years[bool_mu]]) * exp(log(exposure[bool_mu]))
  }else{
    g2 <- as.numeric(mu_deriv[bool_mu])*t(agg_claims[, years[bool_mu]]) * exp(log(exposure[bool_mu]))*as.numeric(nu1[bool_mu])
  }
  
  
  g2 <- by(g2,locs[bool_mu], FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  
  
  
  return(-c(prop_deriv, beta1_deriv, beta2_deriv, g22))
  
  
}




ZIPLN_a <- function(claims, X1,X2, locs, years, agg_claims, exposure, additive, 
                    max_itr = 1000, z = "", lambda = 0, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  start <- ZIP_a_all_in_one(claims, X1, locs, years, agg_claims, additive, exposure = exposure, max_itr = 100, lambda = lambda)
  
  beta1 <- start$beta
  prop <- start$prop
  a <- start$a
  print(beta1)
  print(prop)
  
  
  lower_beta1 <- beta1 -2*abs(beta1)
  upper_beta1 <- beta1 +2*abs(beta1)
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-2, rep(5, ncol(X2)-1))
  upper_beta2 <- c(0.3, rep(5, ncol(X2)-1))
  
  a <- rep(1,p*(p+1)/2)*0.001
  A <- get_W_from_array(a, p)
  lower_a <- rep(1e-4, p*(p+1)/2)
  upper_a <- rep(0.5, p*(p+1)/2)
  
  lower_prop <- 0.01
  upper_prop <- 0.99
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2%*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  phi <- nu2
  
  
  ll <-  sum(log_dzipln(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    out    <- optim(par = c(prop, beta1, beta2, a),
                    fn = Q_ZIPLN_a,
                    gr = Q_ZIPLN_a_deriv,
                    param_prev = c(prop, beta1, beta2, a),
                    locs = locs,
                    claims = claims,
                    exposure =exposure,
                    X1 = X1,
                    X2 = X2,
                    agg_claims = agg_claims,
                    years = years,
                    additive = additive,
                    p = p,
                    lambda = lambda,
                    method = 'L-BFGS-B',
                    lower = c(lower_prop, lower_beta1, lower_beta2, lower_a),
                    upper = c(upper_prop, upper_beta1, upper_beta2, upper_a),
                    control = list(maxit = 1, trace = 0), hessian = calc_hessian)
    print(out$message)
    
    prop <- out$par[1]
    beta1 <- out$par[2:(ncol(X1)+1)]
    beta2 <- out$par[(ncol(X1)+2):(ncol(X1)+ ncol(X2)+1)]
    a <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
    A <- get_W_from_array(a, p)
    
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dzipln(claims,mu,phi, prop), na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_ZIPLN_a_deriv(out$par, out$par, locs, claims, exposure, X1, X2, agg_claims, years,  additive, p, lambda))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, A = A, a=a, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}











# Hurdle Models ------

# Poisson Hurdle------

log_dHP <- function(x,mu,prop, z){
  lab <- (x==0)*1
  lab*log(prop) + (1-lab)*(log((1-prop)) + dpois(x,mu*z, log = TRUE) - log(1- dpois(0,mu*z, log = FALSE)))
}

log_dHP_der <- function(x,mu,prop, z){
  lab <- (x == 0)*1
  
  
  const <- (1/(1- dpois(0,mu*z, log = FALSE)))
  fphi <- lab*(1/prop) - (1-lab)/(1-prop)  # phi is prop...change name?
  fmu <- (1-lab)*(x/mu - z) - (1-lab)*const*exp(-mu*z)*z
  
  return(list(fmu = fmu,fphi = fphi))
}


dhurdle <- function(y, mu, prop){
  
  lab <- (y == 0)*1
  
  return(lab*prop + (1-lab)*(1-prop)*dpois(y, mu)/(1-dpois(0, mu)) )
}


ml_hp_obj <- function(param, claims, exposure, X_mat, additive, prop){
  
  beta <- param[1:(ncol(X_mat))]
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  mu <- nu1*exp(log(exposure))
  
  ll <- sum(log_dHP(claims, mu, prop, 1))
  
  return(-ll)
  
}

ml_hp_obj_deriv <- function(param, locs, claims, exposure, X_mat, years, additive, prop){
  
  beta <- param[1:(ncol(X_mat))]
  
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  mu <- nu1*exp(log(exposure))
  
  mu_deriv <- log_dHP_der(claims, mu,prop, 1)$fmu
  
  
  beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
  
  
  return(-c(beta_deriv))
  
}



ml_hp_psi_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, prop,w){
  
  
  
  
  
  beta <- param[1:(ncol(X_mat))]
  psi <- param[(ncol(X_mat)+1):length(param)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  #print(max(se1$spatial_effect))
  #print(min(mu))
  
  ll <- sum(log_dHP(claims, mu, prop, 1))
  
  
  
  return(-ll)
  
}

ml_hp_psi_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, prop,w){
  
  beta <- param[1:(ncol(X_mat))]
  psi <- param[(ncol(X_mat)+1):length(param)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  mu_deriv <- log_dHP_der(claims, mu,prop, 1)$fmu
  
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
  }else{
    beta_deriv <- t(mu_deriv * nu1*(1+se1$spatial_effect) * exp(log(exposure)) ) %*% X_mat
  }
  
  
  # psi
  agg_effect <- colSums(((w %*% agg_claims[,years, drop = F])))
  
  if(additive){
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect
    
  }else{
    psi_deriv <- mu_deriv*exp(log(exposure)) * se1$agg_effect *nu1
  }
  
  psi_deriv <- tapply(psi_deriv,locs,sum)
  
  return(-c(beta_deriv, psi_deriv))
  
}



ml_hp_a_obj <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, lambda, prop, p){
  
  beta <- param[1:(ncol(X_mat))]
  a <- param[(ncol(X_mat)+1):length(param)]
  A <- get_W_from_array(a,p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <- sum(log_dHP(claims, mu, prop, 1)) - sum(a*lambda)
  
  
  
  return(-ll)
  
}

ml_hp_a_obj_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, additive, lambda, prop,p){
  
  
  beta <- param[1:(ncol(X_mat))]
  a <- param[(ncol(X_mat)+1):length(param)]
  A <- get_W_from_array(a,p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X_mat %*% beta))
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  mu_deriv <- log_dHP_der(claims, mu,prop, 1)$fmu
  
  
  # beta
  if(additive){
    beta_deriv <- t(mu_deriv * nu1 * exp(log(exposure)) ) %*% X_mat
  }else{
    beta_deriv <- t(mu_deriv * nu1*(1+se1$spatial_effect) * exp(log(exposure)) ) %*% X_mat
  }
  
  
  # A
  if(additive){
    g2 <- mu_deriv*t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <- mu_deriv*t(agg_claims[, years]) * exp(log(exposure))*as.numeric(nu1)
  }
  
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  return(-c( beta_deriv, g22))
  
}

HP_all_in_one <- function(claims, X1, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  prop <- mean(claims == 0)
  
  
  out = optim(c(beta),
              fn = ml_hp_obj,
              gr = ml_hp_obj_deriv,
              claims = claims, 
              exposure = exposure, 
              X_mat = X1, 
              additive = additive,
              prop = prop,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(rep(-20, ncol(X1)) ),
              upper = c(rep(10, ncol(X1)) ),
              hessian = T)
  
  
  beta <- out$par[1:(ncol(X1))]
  
  nu1 <- exp(X1 %*% as.matrix(beta))
  
  mu <- nu1*exp(log(exposure))
  
  
  
  return(list(beta1 = beta, prop = prop, mu = mu, out=out))
}

HP_psi_all_in_one <- function(claims, X1, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), max_itr = 1000){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)) , family = 'poisson')$coefficients
  
  psi <- rep(1,p)*0.001
  
  prop <- mean(claims == 0)
  
  
  out = optim(c(beta, psi),
              fn = ml_hp_psi_obj,
              gr = ml_hp_psi_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat = X1, 
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              prop = prop,
              w = w,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(rep(-20, ncol(X1)) ,rep(1e-4, p)),
              upper = c(rep(10, ncol(X1)), rep(1, p*(p+1)/2) ),
              hessian = T)
  
  
  beta <- out$par[1:(ncol(X1))]
  psi <- out$par[(ncol(X1)+1):length(out$par)]
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, psi = psi, prop = prop, mu = mu, out=out))
}

HP_a_all_in_one <- function(claims, X1, locs, years, agg_claims, additive, exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  
  a <- rep(0.1, p*(p+1)/2)
  A <- get_W_from_array(a, p)
  
  prop <- mean(claims == 0)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  
  lower <- rep(1e-8, p*(p+1)/2) 
  
  
  
  out = optim(par = c( beta, a),
              fn = ml_hp_a_obj,
              gr = ml_hp_a_obj_deriv,
              locs = locs,
              claims = claims, 
              exposure = exposure, 
              X_mat = X1, 
              agg_claims = agg_claims, 
              years = years,
              additive = additive,
              lambda = lambda,
              prop = prop,
              p = p,
              method = 'L-BFGS-B',
              control = list(maxit = max_itr),
              lower = c(rep(-20, ncol(X1)) ,lower),
              upper = c(rep(Inf, ncol(X1)), rep(10, p*(p+1)/2) ),
              hessian = T)
  
  
  
  beta <- out$par[1:(ncol(X1))]
  a <- out$par[(ncol(X1)+1):length(out$par)]
  A <- get_W_from_array(a,p) 
  
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- exp(X1 %*% as.matrix(beta))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, A = A,a=a, prop = prop, mu = mu, out=out))
}

# Poisson Gamma Hurdle -----


log_dhpg <- function(y,mu,phi,prop){
  fun=function(z) {
    dhurdle(y, mu*z,prop)* dgamma(z,shape = phi, rate = phi)
  }
  out <- tryCatch({
    lower <- qgamma(0.0001,shape = phi, rate = phi)
    upper <- qgamma(0.9999,shape = phi, rate = phi)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  
  if(!is.na(out)){
    if(out <=1e-6){
      out <- 1e-6
    }
  }
  
  
  return(log(out))
  
}
log_dhpg <- Vectorize(log_dhpg, vectorize.args = c("y", "mu", "phi", "prop"))

Q_hpg_numerator <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev){
  fun=function(z) {
    (log_dHP(y, mu, prop, z) +  dgamma(z, shape = phi, rate = phi,log = TRUE))*dhurdle(y, mu_prev*z, prop_prev)* dgamma(z,shape = phi_prev, rate = phi_prev)
  }
  out <- tryCatch({
    lower <- qgamma(0.0001,shape = phi_prev, rate = phi_prev)
    upper <- qgamma(0.9999,shape = phi_prev, rate = phi_prev)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  
  return(out)
  
}

Q_hpg_numerator <- Vectorize(Q_hpg_numerator, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev"))

Q_hpg <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev, denominator){
  
  numerator <- Q_hpg_numerator(y,mu, phi, prop, mu_prev, phi_prev, prop_prev)
  return(numerator/denominator)
  
}


mu_deriv_hurdle_pg <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    (log_dHP_der(y, mu, prop, z)$fmu)* dhurdle(y, mu_prev*z, prop_prev)* dgamma(z, shape = phi_prev, rate = phi_prev)
  }
  out <- tryCatch({
    lower <- qgamma(0.0001,shape = phi_prev, rate = phi_prev)
    upper <- qgamma(0.9999,shape = phi_prev, rate = phi_prev)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
mu_deriv_hurdle_pg <- Vectorize(mu_deriv_hurdle_pg, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))


phi_deriv_hurdle_pg <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    
    (log(phi) + 1 - digamma(phi) + log(z) - z)* dhurdle(y, mu_prev*z, prop_prev)* dgamma(z, shape = phi_prev, rate = phi_prev)
  }
  out <- tryCatch({
    lower <- qgamma(0.0001,shape = phi_prev, rate = phi_prev)
    upper <- qgamma(0.9999,shape = phi_prev, rate = phi_prev)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
phi_deriv_hurdle_pg <- Vectorize(phi_deriv_hurdle_pg, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))


Q_HPG_no_spatial <- function(param, param_prev, claims, exposure, X1,X2, prop){
  
  prop_prev <- prop
  
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  
  
  # current estimate
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  denominator <- exp(log_dhpg(claims, mu, phi, prop))
  
  
  ll <- Q_hpg(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  bool <- !(is.na(ll) | is.infinite(ll))
  ll <- ll[bool]
  
  
  
  
  return(-sum(ll))
  
}

Q_HPG_no_spatial_deriv <- function(param, param_prev,  claims, exposure, X1,X2, prop){
  
  
  
  
  
  prop_prev <- prop
  # Previous 
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  
  
  # current estimate
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  denominator <- exp(log_dhpg(claims, mu, phi, prop))
  
  # 
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_hurdle_pg(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv))
  
  beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop = FALSE]
  
  
  # beta 2 deriv
  # beta 2 deriv
  
  phi_deriv <- phi_deriv_hurdle_pg(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv))
  
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = FALSE]
  
  return(-c(beta1_deriv, beta2_deriv))
  
  
}
HPG_no_spatial <- function(claims, X1,X2, exposure = rep(1, length(claims)), 
                           max_itr = 1000, z = "", tol = 1e-3){
  
  # Set initial parameters and bounds
  
  out <- HP_all_in_one(claims, X1, exposure, max_itr )
  beta1 <- out$beta
  
  
  lower_beta1 <- beta1 - 2*abs(beta1)
  upper_beta1 <- beta1 + 2*abs(beta1)
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-0.2, rep(5, ncol(X2)-1))
  upper_beta2 <- c(2, rep(5, ncol(X2)-1))
  
  
  prop <- mean(claims == 0)
  #lower_prop <- 1e-3
  #upper_prop <- 0.999
  
  
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2%*% as.matrix(beta2))
  mu <- nu1*exp(log(exposure))
  phi <- nu2
  
  
  ll <-  sum(log_dhpg(claims,mu,phi, prop),na.rm = TRUE)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F 
    }
    
    print(beta1)
    print(beta2)
    print(prop)
    out  <- optim(par = c(beta1, beta2),
                  fn = Q_HPG_no_spatial,
                  gr  = Q_HPG_no_spatial_deriv,
                  param_prev = c(beta1, beta2),
                  claims = claims,
                  exposure = exposure,
                  X1 = X1,
                  X2 = X2,
                  prop = prop,
                  method = 'L-BFGS-B',
                  lower = c(lower_beta1, lower_beta2),
                  upper = c(upper_beta1, upper_beta2),
                  control = list(maxit = 2, trace = 0, factr = 1), hessian = calc_hessian)
    
    print(out$message)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    
    
    # Create new estimates
    nu1 <- as.numeric(exp(X1 %*% beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    mu <- nu1*exp(log(exposure))
    
    # finalize
    new_ll <- sum(log_dhpg(claims,mu,phi, prop),na.rm = TRUE)
    vre <- abs(new_ll - ll)
    ll <- new_ll

    print(paste('At iteration',itr,'log likelihood =',new_ll))
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    itr <- itr + 1
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPG_no_spatial_deriv(out$par, out$par, claims, exposure, X1, X2, prop))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}


Q_HPG_psi<- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, w, agg_claims, additive){
  prop_prev <- prop
  
  # Previous 
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi <- param[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpg(claims, mu, phi, prop))
  
  
  ll <- Q_hpg(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_HPG_psi_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, w, agg_claims, additive){
  prop_prev <- prop
  # Previous 
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi <- param[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpg(claims, mu, phi, prop))
  
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_hurdle_pg(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv))
  
  phi_deriv <- phi_deriv_hurdle_pg(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv))
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop = F]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop = F]
  }
  
  
  # beta 2 deriv
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = F]
  
  
  if(additive){
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu]
    
  }else{
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu] *nu1[bool_mu]
  }
  
  psi_deriv <- tapply(psi_deriv,locs[bool_mu],sum)
  
  
  
  return(-c(beta1_deriv, beta2_deriv, psi_deriv))
  
  
}


HPG_psi <- function(claims, X1,X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), 
                    max_itr = 1000, z = "", max_int = Inf, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  out <- HP_psi_all_in_one(claims, X1, locs, years, agg_claims, w, additive,  exposure, max_itr )
  beta1 <- out$beta
  
  
  lower_beta1 <- beta1 - abs(beta1)
  upper_beta1 <- beta1 + abs(beta1)
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-0.5, rep(5, ncol(X2)-1))
  upper_beta2 <- c(2, rep(5, ncol(X2)-1))
  
  
  psi <- out$psi
  lower_psi <- rep(1e-6, p)
  upper_psi <- rep(0.5, p)
  
  prop <- mean(claims == 0)
  #lower_prop <- 1e-3
  #upper_prop <- 0.999
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll <-  sum(log_dhpg(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    print(beta1)
    print(beta2)
    print(prop)
    out    <- optim(par = c(beta1, beta2, psi),
                    fn = Q_HPG_psi,
                    gr = Q_HPG_psi_deriv,
                    param_prev = c(beta1, beta2, psi),
                    locs = locs,
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    years = years,
                    prop = prop,
                    w = w,
                    agg_claims = agg_claims,
                    additive = additive,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2, lower_psi),
                    upper = c(upper_beta1, upper_beta2, upper_psi),
                    control = list(maxit = 1, trace = 0, factr = 1e-3), hessian = calc_hessian)
    print(out$message)
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    psi <- out$par[(ncol(X1)+ ncol(X2)+1):length(out$par)]
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dhpg(claims,mu,phi, prop), na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    
    print(paste('At iteration',itr,'log likelihood =',new_ll))

    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPG_psi_deriv(out$par, out$par, locs, claims, exposure, X1, X2,years, prop, w, agg_claims, additive))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, psi = psi, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}



Q_HPG_a<- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop,  agg_claims, additive, p, lambda){
  prop_prev <- prop
  # Previous 

  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate

  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a <- param[(ncol(X1)+ ncol(X2)+1):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpg(claims, mu, phi, prop))
  
  
  ll <- Q_hpg(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator) - sum(a*lambda)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_HPG_a_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, agg_claims, additive, p, lambda){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a <- param[(ncol(X1)+ ncol(X2)+1):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpg(claims, mu, phi, prop))
  
  # print("prop deriv")
  
  
  
  mu_deriv <- mu_deriv_hurdle_pg(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv))
  
  phi_deriv <- phi_deriv_hurdle_pg(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv))
  
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop = F]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop = F]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop=FALSE]
  
  
  # A
  if(additive){
    g2 <- as.numeric(mu_deriv[bool_phi])*t(agg_claims[, years[bool_phi]]) * exp(log(exposure[bool_phi]))
  }else{
    g2 <- as.numeric(mu_deriv[bool_phi])*t(agg_claims[, years[bool_phi]]) * exp(log(exposure[bool_phi]))*as.numeric(nu1[bool_phi])
  }
  
  
  g2 <- by(g2,locs[bool_phi], FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  
  
  
  
  return(-c(beta1_deriv, beta2_deriv, g22))
  
  
}



HPG_a <- function(claims, X1,X2, locs, years, agg_claims,  additive, exposure = rep(1, length(claims)), 
                  max_itr = 1000, z = "", lambda = 0, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  
  # Set initial parameters and bounds
  out <- HP_a_all_in_one(claims, X1, locs, years, agg_claims, additive,  exposure, max_itr )
  beta1 <- out$beta
  
  
  lower_beta1 <- beta1 - 2*abs(beta1)
  upper_beta1 <- beta1 + 2*abs(beta1)
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-1, rep(5, ncol(X2)-1))
  upper_beta2 <- c(3, rep(5, ncol(X2)-1))
  
  
  a <- out$a
  A <- get_W_from_array(a, p)
  lower_a <- rep(1e-6, p*(p+1)/2)
  upper_a <- rep(0.5, p*(p+1)/2)
  
  prop <- mean(claims == 0)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll <-  sum(log_dhpg(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    print(beta1)
    print(beta2)
    print(prop)
    out    <- optim(par = c(beta1, beta2, a),
                    fn = Q_HPG_a,
                    gr = Q_HPG_a_deriv,
                    param_prev = c(beta1, beta2, a),
                    locs = locs,
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    years = years,
                    prop = prop,
                    agg_claims = agg_claims,
                    additive = additive,
                    p = p,
                    lambda = lambda,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2, lower_a),
                    upper = c(upper_beta1, upper_beta2, upper_a),
                    control = list(maxit = 1, trace = 0, factr = 1e-1), hessian = calc_hessian)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    a <- out$par[(ncol(X1)+ ncol(X2)+1):length(out$par)]
    A <- get_W_from_array(a,p)
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dhpg(claims,mu,phi, prop), na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPG_a_deriv(out$par, out$par, locs, claims, exposure, X1, X2,years, prop, agg_claims, additive, p, lambda))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, a = a, A=A, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}




# Poisson Inverse Gaussian Hurdle -----


log_dhpig <- function(y,mu,phi,prop){
  fun=function(z) {
    dhurdle(y, mu*z,prop)* actuar::dinvgauss(z,mean = 1, shape = phi^2)
  }
  out <- tryCatch({
    lower <- actuar::qinvgauss(0.0001,mean = 1, shape = phi^2)
    upper <- actuar::qinvgauss(0.9999,mean = 1, shape = phi^2)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  
  if(!is.na(out)){
    if(out <=1e-8){
      out <- 1e-8
    }
  }
  
  
  return(log(out))
  
}
log_dhpig <- Vectorize(log_dhpig, vectorize.args = c("y", "mu", "phi", "prop"))




Q_hpig_numerator <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev){
  fun=function(z) {
    (log_dHP(y, mu, prop, z) +  actuar::dinvgauss(z,mean = 1, shape = phi^2,log = TRUE))*dhurdle(y, mu_prev*z, prop_prev)* actuar::dinvgauss(z,mean = 1, shape = phi_prev^2)
  }
  out <- tryCatch({
    lower <- actuar::qinvgauss(0.0001,mean = 1, shape = phi_prev^2)
    upper <- actuar::qinvgauss(0.9999,mean = 1, shape = phi_prev^2)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  
  return(out)
  
}

Q_hpig_numerator <- Vectorize(Q_hpig_numerator, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev"))


Q_hpig <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev, denominator){
  
  numerator <- Q_hpig_numerator(y,mu, phi, prop, mu_prev, phi_prev, prop_prev)
  return(numerator/denominator)
  
}



mu_deriv_hurdle_pig <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    (log_dHP_der(y, mu, prop, z)$fmu)* dhurdle(y, mu_prev*z, prop_prev)* actuar::dinvgauss(z,mean = 1, shape = phi_prev^2)
  }
  out <- tryCatch({
    lower <- actuar::qinvgauss(0.0001,mean = 1, shape = phi_prev^2)
    upper <- actuar::qinvgauss(0.9999,mean = 1, shape = phi_prev^2)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
mu_deriv_hurdle_pig <- Vectorize(mu_deriv_hurdle_pig, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))


phi_deriv_hurdle_pig <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    
    (1/phi + 2*phi - phi/z - phi*z)* dhurdle(y, mu_prev*z, prop_prev)* actuar::dinvgauss(z,mean = 1, shape = phi_prev^2)
  }
  out <- tryCatch({
    lower <- actuar::qinvgauss(0.0001,mean = 1, shape = phi_prev^2)
    upper <- actuar::qinvgauss(0.9999,mean = 1, shape = phi_prev^2)
    integrate(fun, lower, upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
phi_deriv_hurdle_pig <- Vectorize(phi_deriv_hurdle_pig, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))


Q_HPIG_no_spatial <- function(param, param_prev,  claims, exposure, X1,X2,prop){
  
  prop_prev <- prop
  
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  denominator <- exp(log_dhpig(claims, mu, phi, prop))
  
  
  ll <- Q_hpig(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  bool <- !(is.na(ll) | is.infinite(ll) | is.nan(ll))
  print(paste("Nr removed ", sum(!bool)))
  ll <- ll[bool]
  
  return(-sum(ll))
  
}

Q_HPIG_no_spatial_deriv <- function(param, param_prev, claims, exposure, X1,X2, prop){
  
  
  
  
  
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  denominator <- exp(log_dhpig(claims, mu, phi, prop))
  mu_deriv <- mu_deriv_hurdle_pig(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv))
  print(paste("Nr removed mu ", sum(!bool_mu)))
  
  
  
  beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop = FALSE]
  
  
  # beta 2 deriv
  phi_deriv <- phi_deriv_hurdle_pig(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv))
  print(paste("Nr removed phi ", sum(!bool_phi)))
  
  
  
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = FALSE]
  
  
  
  
  return(-c(beta1_deriv, beta2_deriv))
  
  
}
HPIG_no_spatial <- function(claims, X1,X2, exposure = rep(1, length(claims)), 
                            max_itr = 1000, z = "", max_int = Inf, tol = 1e-3){
  
  
  # Set initial parameters and bounds
  
  out <- HP_all_in_one(claims, X1, exposure, max_itr )
  beta1 <- out$beta
  
  
  lower_beta1 <- beta1 - abs(beta1)
  upper_beta1 <- beta1 + abs(beta1)
  
  #print(lower_beta1)
  #print(upper_beta1)
  
  beta2 <- c(1, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-1, rep(5, ncol(X2)-1))
  upper_beta2 <- c(3, rep(5, ncol(X2)-1))
  
  
  
  
  prop <- mean(claims == 0)
  
  
  
  
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2%*% as.matrix(beta2))
  mu <- nu1*exp(log(exposure))
  phi <- nu2
  
  
  ll <-  sum(log_dhpig(claims,mu,phi, prop),na.rm = TRUE)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F 
    }
    
    print(beta1)
    print(beta2)
    print(prop)
    out  <- optim(par = c(beta1, beta2),
                  fn = Q_HPIG_no_spatial,
                  gr = Q_HPIG_no_spatial_deriv,
                  param_prev = c(beta1, beta2),
                  claims = claims,
                  exposure = exposure,
                  X1 = X1,
                  X2 = X2,
                  prop = prop,
                  method = 'L-BFGS-B',
                  lower = c(lower_beta1, lower_beta2),
                  upper = c(upper_beta1, upper_beta2),
                  control = list(maxit = 2, trace = 0, factr = 1e-5), hessian = calc_hessian)
    print(out$message)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    
    
    # Create new estimates
    nu1 <- as.numeric(exp(X1 %*% beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    mu <- nu1*exp(log(exposure))
    
    # finalize
    new_ll <- sum(log_dhpig(claims,mu,phi, prop),na.rm = TRUE)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    itr <- itr + 1
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPIG_no_spatial_deriv(out$par, out$par, claims, exposure, X1, X2, prop))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}


Q_HPIG_psi<- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, w, agg_claims, additive){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi <- param[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpig(claims, mu, phi, prop))
  
  
  ll <- Q_hpig(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  bool <- !(is.na(ll) | is.infinite(ll) | is.nan(ll))
  if(any(!bool)){
    print(paste("Nr removed ", sum(!bool)))
  }
  ll <- ll[bool]
  
  return(-sum(ll))
  
}

Q_HPIG_psi_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, w, agg_claims, additive){
  prop_prev <- prop
  # Previous 
  
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi <- param[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpig(claims, mu, phi, prop))
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_hurdle_pig(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv))
  if(any(!bool_mu)){
    print("NA in mu deriv")
  }
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_hurdle_pig(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv))
  
  if(any(!bool_phi)){
    
    print("NA in mu deriv")
  }
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = F]
  
  
  if(additive){
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu]
    
  }else{
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu] *nu1[bool_mu]
  }
  
  psi_deriv <- tapply(psi_deriv,locs[bool_mu],sum)
  
  
  
  
  return(-c(beta1_deriv, beta2_deriv, psi_deriv))
  
  
}


HPIG_psi <- function(claims, X1,X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), 
                     max_itr = 1000, z = "", max_int = Inf, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  out <- HP_psi_all_in_one(claims, X1, locs, years, agg_claims, w, additive,  exposure, max_itr )
  beta1 <- out$beta
  
  lower_beta1 <- beta1 - abs(beta1)
  upper_beta1 <- beta1 + abs(beta1)
  
  beta2 <- c(0.1, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-1, rep(5, ncol(X2)-1))
  upper_beta2 <- c(3, rep(5, ncol(X2)-1))
  
  psi <- out$psi
  lower_psi <- rep(1e-6, p)
  upper_psi <- rep(0.5, p)
  
  prop <- mean(claims == 0)
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll <-  sum(log_dhpig(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    print(beta1)
    print(beta2)
    print(prop)
    out    <- optim(par = c(beta1, beta2, psi),
                    fn = Q_HPIG_psi,
                    gr = Q_HPIG_psi_deriv,
                    param_prev = c(beta1, beta2, psi),
                    locs = locs,
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    years = years,
                    prop = prop,
                    w = w,
                    agg_claims = agg_claims,
                    additive = additive,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2, lower_psi),
                    upper = c(upper_beta1, upper_beta2, upper_psi),
                    control = list(maxit = 1, trace = 0, factr = 1e-3), hessian = calc_hessian)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    psi <- out$par[(ncol(X1)+ ncol(X2)+1):length(out$par)]
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dhpig(claims,mu,phi, prop), na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPIG_psi_deriv(out$par, out$par, locs, claims, exposure, X1, X2,years, prop, w, agg_claims, additive))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, psi = psi, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}



Q_HPIG_a<- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop,  agg_claims, additive, p, lambda){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a <- param[(ncol(X1)+ ncol(X2)+1):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpig(claims, mu, phi, prop))
  
  
  ll <- Q_hpig(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  bool <- !(is.na(ll) | is.infinite(ll) | is.nan(ll))
  if(any(!bool)){
    print(paste("Nr removed ", sum(!bool)))
  }
  ll <- ll[bool]
  
  return(-sum(ll))
  
}

Q_HPIG_a_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, agg_claims, additive, p, lambda){
  
  prop_prev <- prop
  
  # Previous 
  
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a <- param[(ncol(X1)+ ncol(X2)+1):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpig(claims, mu, phi, prop))
  
  
  
  # beta 1 deriv
  mu_deriv <- mu_deriv_hurdle_pig(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv))
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_hurdle_pig(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv))
  
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = F]
  
  
  # A
  if(additive){
    g2 <- as.numeric(mu_deriv[bool_mu])*t(agg_claims[, years[bool_mu]]) * exp(log(exposure[bool_mu]))
  }else{
    g2 <- as.numeric(mu_deriv[bool_mu])*t(agg_claims[, years[bool_mu]]) * exp(log(exposure))*as.numeric(nu1[bool_mu])
  }
  
  
  g2 <- by(g2,locs[bool_mu], FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  return(-c(beta1_deriv, beta2_deriv, g22))
  
  
}



HPIG_a <- function(claims, X1,X2, locs, years, agg_claims,  additive, exposure = rep(1, length(claims)), 
                   max_itr = 1000, z = "", max_int = Inf, lambda = 0, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  
  # Set initial parameters and bounds
  out <- HP_a_all_in_one(claims, X1, locs, years, agg_claims, additive,  exposure, max_itr )
  beta1 <- out$beta
  
  
  lower_beta1 <- beta1 - abs(beta1)
  upper_beta1 <- beta1 + abs(beta1)
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-0.3, rep(5, ncol(X2)-1))
  upper_beta2 <- c(2, rep(5, ncol(X2)-1))
  
  
  a <- out$a
  A <- get_W_from_array(a, p)
  lower_a <- rep(1e-6, p*(p+1)/2)
  upper_a <- rep(0.5, p*(p+1)/2)
  
  prop <- mean(claims == 0)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  ll <-  sum(log_dhpig(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    print(beta1)
    print(beta2)
    print(prop)
    out    <- optim(par = c(beta1, beta2, a),
                    fn = Q_HPIG_a,
                    gr = Q_HPIG_a_deriv,
                    param_prev = c(beta1, beta2, a),
                    locs = locs,
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    years = years,
                    prop = prop,
                    agg_claims = agg_claims,
                    additive = additive,
                    p = p,
                    lambda = lambda,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2, lower_a),
                    upper = c(upper_beta1, upper_beta2, upper_a),
                    control = list(maxit = 1, trace =3, factr = 1), hessian = calc_hessian)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    a <- out$par[(ncol(X1)+ ncol(X2)+1):length(out$par)]
    A <- get_W_from_array(a,p)
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dhpig(claims,mu,phi, prop), na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPIG_a_deriv(out$par, out$par, locs, claims, exposure, X1, X2,years, prop, agg_claims, additive, p, lambda))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, a = a, A=A, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}


# Poisson Log-Normal Hurdle ------

log_dhpln <- function(y,mu,phi,prop){
  fun=function(z) {
    dhurdle(y, mu*z,prop)* dlno(z,phi)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun,lower,upper)$value
  }, error = function(e) NA)
  
  if(!is.na(out)){
    if(out <= 1e-8){
      out <- 1e-8
    }
  }
  
  return(log(out))
  
}
log_dhpln <- Vectorize(log_dhpln, vectorize.args = c("y", "mu", "phi", "prop"))




Q_hpln_numerator <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev){
  fun=function(z) {
    (log_dHP(y, mu, prop, z) +  log(dlno(z,phi)))* dhurdle(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun,lower,upper)$value
  }, error = function(e) NA)
  
  
  return(out)
  
}

Q_hpln_numerator <- Vectorize(Q_hpln_numerator, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev"))


Q_hpln <- function(y,mu,phi,prop,mu_prev, phi_prev,prop_prev, denominator){
  
  numerator <- Q_zipln_numerator(y,mu, phi, prop, mu_prev, phi_prev, prop_prev)
  return(numerator/denominator)
  
}



mu_deriv_hurdle_pln <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    (log_dHP_der(y, mu, prop, z)$fmu)* dhurdle(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun,lower,upper)$value
  }, error = function(e) NA)
  
  
  return(out/denominator)
}
mu_deriv_hurdle_pln <- Vectorize(mu_deriv_hurdle_pln, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))

phi_deriv_hurdle_pln <- function(y,mu,phi,prop, mu_prev, phi_prev, prop_prev, denominator){
  fun=function(z) {
    
    (-1/phi + (log(z)^2)/(phi ^3) - (phi)/4)* dhurdle(y, mu_prev*z, prop_prev)* dlno(z,phi_prev)
  }
  out <- tryCatch({
    lower <- qlnorm(0.0001,meanlog = -phi^2/2, sdlog = phi)
    upper <- qlnorm(0.9999,meanlog = -phi^2/2, sdlog = phi)
    integrate(fun,lower,upper)$value
  }, error = function(e) NA)
  
  return(out/denominator)
}
phi_deriv_hurdle_pln <- Vectorize(phi_deriv_hurdle_pln, vectorize.args = c("y", "mu", "phi", "prop", "mu_prev", "phi_prev", "prop_prev", "denominator"))


Q_HPLN_no_spatial <- function(param, param_prev, claims, exposure, X1,X2,  prop){
  
  prop_prev <- prop
  
  # Previous 
  
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  
  denominator <- exp(log_dhpln(claims, mu_prev, phi_prev, prop_prev))
  
  
  # current estimate
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  
  
  ll <- Q_hpln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_HPLN_no_spatial_deriv <- function(param, param_prev, claims, exposure, X1,X2, prop){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  mu_prev <- nu1_prev*exp(log(exposure))
  denominator <- exp(log_dhpln(claims, mu_prev, phi_prev, prop_prev))
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  
  
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  mu <- nu1*exp(log(exposure))
  
  
  # print("beta 1 deriv")
  # beta 1 deriv
  mu_deriv <- mu_deriv_hurdle_pln(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv) | is.nan(mu_deriv))
  
  
  beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop =FALSE]
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_hurdle_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv) | is.nan(phi_deriv))
  
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = FALSE]
  
  
  return(-c(beta1_deriv, beta2_deriv))
  
  
}
HPLN_no_spatial <- function(claims, X1,X2, exposure = rep(1, length(claims)), 
                            max_itr = 1000, z = "", tol = 1e-3){
  
  
  # Set initial parameters and bounds
  start <- HP_all_in_one(claims, X1, exposure = exposure, max_itr = 1000)
  beta1 <- start$beta
  
  lower_beta1 <- beta1 - 2*abs(beta1)
  upper_beta1 <- beta1 + 2*abs(beta1)
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-5, rep(5, ncol(X2)-1))
  upper_beta2 <- c(0.5, rep(5, ncol(X2)-1))
  
  
  prop <- mean(claims == 0)
  
  
  
  nu1 <- exp(X1 %*% as.matrix(beta1))
  nu2 <- exp(X2%*% as.matrix(beta2))
  mu <- nu1*exp(log(exposure))
  phi <- nu2
  
  
  ll <-  sum(log_dhpln(claims,mu,phi, prop))
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    
    out    <- optim(par = c(beta1, beta2),
                    fn = Q_HPLN_no_spatial,
                    gr = Q_HPLN_no_spatial_deriv,
                    param_prev = c(beta1, beta2),
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    prop = prop,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2),
                    upper = c(upper_beta1, upper_beta2),
                    control = list(maxit = 1, trace = 0), hessian = calc_hessian)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    
    
    # Create new estimates
    nu1 <- as.numeric(exp(X1 %*% beta1))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    mu <- nu1*exp(log(exposure))
    
    # finalize
    new_ll <-  sum(log_dhpln(claims,mu,phi, prop))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPLN_no_spatial_deriv(out$par, out$par, claims, exposure, X1, X2, prop))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}


Q_HPLN_psi<- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, w, agg_claims, additive){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi <- param[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  denominator <- exp(log_dhpln(claims, mu, phi, prop))
  
  ll <- Q_hpln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    # print(mu[is.na(ll)])
    # print(phi[is.na(ll)])
    # print(claims[is.na(ll)])
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_HPLN_psi_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, w, agg_claims, additive){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1_prev <- get_spatial_aggregate_new(locs, w, psi_prev, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  psi <- param[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpln(claims, mu, phi, prop))
  
  # print("beta 1 deriv")
  # beta 1 deriv
  mu_deriv <- mu_deriv_hurdle_pln(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv) | is.nan(mu_deriv))
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_hurdle_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv) | is.nan(phi_deriv))
  
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = FALSE]
  
  
  
  if(additive){
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu]
    
  }else{
    psi_deriv <- mu_deriv[bool_mu]*exp(log(exposure[bool_mu])) * se1$agg_effect[bool_mu] *nu1[bool_mu]
  }
  
  psi_deriv <- tapply(psi_deriv,locs[bool_mu],sum)
  
  
  
  
  return(-c(beta1_deriv, beta2_deriv, psi_deriv))
  
  
}


HPLN_psi <- function(claims, X1,X2, locs, years, agg_claims, w, additive, exposure = rep(1, length(claims)), 
                     max_itr = 1000, z = "", max_int = Inf, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  
  # Set initial parameters and bounds
  start <- HP_psi_all_in_one(claims, X1, locs, years,agg_claims,w,additive,exposure, 100)
  beta1 <- start$beta
  psi <- start$psi
  
  print(beta1)
  
  lower_beta1 <- beta1 - 2*abs(beta1)
  upper_beta1 <- beta1 + 2*abs(beta1)
  
  
  prop <- mean(claims == 0)
  
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-5, rep(5, ncol(X2)-1))
  upper_beta2 <- c(0.5, rep(5, ncol(X2)-1))
  
  psi <- rep(0.001, p)
  lower_psi <- rep(1e-6, p)
  upper_psi <- rep(1, p)
  
  
  
  
  
  
  se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll <-  sum(log_dhpln(claims,mu,phi, prop))
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    print(beta1)
    print(prop)
    out    <- optim(par = c(beta1, beta2, psi),
                    fn = Q_HPLN_psi,
                    gr = Q_HPLN_psi_deriv,
                    param_prev = c(beta1, beta2, psi),
                    locs = locs,
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    years = years,
                    prop = prop,
                    w = w,
                    agg_claims = agg_claims,
                    additive = additive,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2, lower_psi),
                    upper = c(upper_beta1, upper_beta2, upper_psi),
                    control = list(maxit = 1, trace = 0, factr = 1e-3), hessian = calc_hessian)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    psi <- out$par[(ncol(X1)+ ncol(X2)+1):length(out$par)]
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, w, psi, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dhpln(claims,mu,phi, prop))
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPLN_psi_deriv(out$par, out$par, locs, claims, exposure, X1, X2,years, prop, w, agg_claims, additive))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, psi = psi, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}




Q_HPLN_a<- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, agg_claims, additive, p, lambda){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a <- param[(ncol(X1)+ ncol(X2)+1):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpln(claims, mu, phi, prop))
  
  
  ll <- Q_hpln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev,  denominator) - sum(a*lambda)
  
  if(any(is.na(ll))){
    print("NA in ll")
    print(sum(is.na(ll)))
    ll[is.na(ll)] <- 0
  }
  
  return(-sum(ll))
  
}

Q_HPLN_a_deriv <- function(param, param_prev, locs, claims, exposure, X1,X2, years, prop, agg_claims, additive, p, lambda){
  prop_prev <- prop
  # Previous 
  #prop_prev <- param_prev[1]
  beta1_prev <- param_prev[1:(ncol(X1)+0)]
  beta2_prev <- param_prev[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a_prev <- param_prev[(ncol(X1)+ ncol(X2)+1):length(param_prev)]
  A_prev <- get_W_from_array(a_prev, p)
  
  se1_prev <- get_spatial_aggregate_new(locs, A_prev, 1, agg_claims, years)
  nu1_prev <- as.numeric(exp(X1 %*% as.matrix(beta1_prev)))
  nu2_prev <- exp(X2 %*% as.matrix(beta2_prev))
  phi_prev <- nu2_prev
  
  if(additive){
    mu_prev <- (nu1_prev + se1_prev$spatial_effect)*exp(log(exposure))
  }else{
    mu_prev <- (nu1_prev * (1+se1_prev$spatial_effect))*exp(log(exposure))
  }
  
  
  
  
  # current estimate
  #prop <- param[1]
  beta1 <- param[1:(ncol(X1)+0)]
  beta2 <- param[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
  a <- param[(ncol(X1)+ ncol(X2)+1):length(param)]
  A <- get_W_from_array(a, p)
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  denominator <- exp(log_dhpln(claims, mu, phi, prop))
  
  mu_deriv <- mu_deriv_hurdle_pln(claims,mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_mu <- !(is.na(mu_deriv) | is.infinite(mu_deriv) | is.nan(mu_deriv))
  
  
  if(additive){
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu] * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }else{
    beta1_deriv <- t(mu_deriv[bool_mu] * nu1[bool_mu]*(1+se1$spatial_effect[bool_mu]) * exp(log(exposure[bool_mu])) ) %*% X1[bool_mu,,drop=FALSE]
  }
  
  
  # beta 2 deriv
  #print("beta 2 deriv")
  phi_deriv <- phi_deriv_hurdle_pln(claims, mu, phi, prop, mu_prev, phi_prev, prop_prev, denominator)
  bool_phi <- !(is.na(phi_deriv) | is.infinite(phi_deriv) | is.nan(phi_deriv))
  
  beta2_deriv <-  t( phi_deriv[bool_phi]*nu2[bool_phi] ) %*% X2[bool_phi,,drop = FALSE]
  
  
  # A
  if(additive){
    g2 <- as.numeric(mu_deriv[bool_mu])*t(agg_claims[, years[bool_mu]]) * exp(log(exposure[bool_mu]))
  }else{
    g2 <- as.numeric(mu_deriv[bool_mu])*t(agg_claims[, years[bool_mu]]) * exp(log(exposure[bool_mu]))*as.numeric(nu1[bool_mu])
  }
  
  
  g2 <- by(g2,locs[bool_mu], FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) <- 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  
  
  
  
  return(-c(beta1_deriv, beta2_deriv, g22))
  
  
}



HPLN_a <- function(claims, X1,X2, locs, years, agg_claims,  additive, exposure = rep(1, length(claims)), 
                   max_itr = 1000, z = "", lambda = 0, tol = 1e-3){
  
  
  p <- length(unique(locs))
  
  
  # Set initial parameters and bounds
  start <- HP_a_all_in_one(claims, X1, locs, years,agg_claims,additive, exposure = exposure, max_itr = 100, lambda = lambda)
  
  beta1 <- start$beta
  #a <- start$a
  
  lower_beta1 <- beta1 - 2*abs(beta1)
  upper_beta1 <- beta1 + 2*abs(beta1)
  
  
  
  print(beta1)
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  lower_beta2 <- c(-5, rep(5, ncol(X2)-1))
  upper_beta2 <- c(0.5, rep(5, ncol(X2)-1))
  
  a <-rep(1e-5, p*(p+1)/2)
  A <- get_W_from_array(a, p)
  lower_a <- rep(1e-6, p*(p+1)/2)
  upper_a <- rep(0.5, p*(p+1)/2)
  print(A)
  
  prop <- mean(claims == 0)
  
  
  print("beta1")
  print(beta1)
  print(lower_beta1)
  print(upper_beta1)
  
  print("beta2")
  print(beta2)
  print(lower_beta2)
  print(upper_beta2)
  
  
  se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
  nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  phi <- nu2
  
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  ll <-  sum(log_dhpln(claims,mu,phi, prop), na.rm = T)
  print(paste("ll at start = ", ll))
  itr <- 1
  has_converged <- FALSE
  last_iteration <- FALSE
  while((itr <=max_itr) & (!has_converged)) {
    
    if(last_iteration || itr == max_itr){
      calc_hessian <- T
    }else{
      calc_hessian <- F
    }
    
    out    <- optim(par = c(beta1, beta2, a),
                    fn = Q_HPLN_a,
                    gr = Q_HPLN_a_deriv,
                    param_prev = c(beta1, beta2, a),
                    locs = locs,
                    claims = claims,
                    exposure = exposure,
                    X1 = X1,
                    X2 = X2,
                    years = years,
                    prop = prop,
                    agg_claims = agg_claims,
                    additive = additive,
                    p = p,
                    lambda = lambda,
                    method = 'L-BFGS-B',
                    lower = c(lower_beta1, lower_beta2, lower_a),
                    upper = c(upper_beta1, upper_beta2, upper_a),
                    control = list(maxit = 1, trace = 3), hessian = calc_hessian)
    
    #prop <- out$par[1]
    beta1 <- out$par[1:(ncol(X1)+0)]
    beta2 <- out$par[(ncol(X1)+1):(ncol(X1)+ ncol(X2)+0)]
    a <- out$par[(ncol(X1)+ ncol(X2)+1):length(out$par)]
    A <- get_W_from_array(a,p)
    
    # Create new estimates
    se1 <- get_spatial_aggregate_new(locs, A, 1, agg_claims, years)
    nu1 <- as.numeric(exp(X1 %*% as.matrix(beta1)))
    nu2 <- exp(X2 %*% as.matrix(beta2))
    phi <- nu2
    
    if(additive){
      mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
    }else{
      mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
    }
    
    # finalize
    new_ll <-  sum(log_dhpln(claims,mu,phi, prop), na.rm = T)
    vre <- abs(new_ll - ll)
    ll <- new_ll
    #if(itr %% 10 ==0){
    print(beta1)
    print(beta2)
    print(paste('At iteration',itr,'log likelihood =',new_ll))
    #}
    
    if(last_iteration){
      has_converged <- TRUE
    }
    
    if(vre <= tol){
      print("Stopping converged")
      last_iteration <-  TRUE
    }
    
    
    
    itr <- itr + 1
    
  }
  
  
  print("finalising")
  
  H_par <- out$hessian
  
  par_s <- as.numeric(Q_HPLN_a_deriv(out$par, out$par, locs, claims, exposure, X1, X2,years, prop, agg_claims, additive, p, lambda))
  par_H <- H_par + outer(par_s, par_s)
  
  
  return(list(beta1 = beta1, beta2 = beta2,  prop = prop, a = a, A = A, mu = mu, phi = phi, par_H = par_H, out = out, par_s = par_s, H_par = H_par))
}





# Generate Data -----


generate_no_spatial <- function(beta2_true, type = NULL, dist = "ordinary", n = 100){
  n <- n
  seed <- 10
  years <- 8  # number of time points
  area <- 5  # number of regions
  set.seed(seed)
  exposure <-   rpois(n, 10)#rep(1,n)#
  exposure <- rep(exposure, years)
  
  # store claim data
  locs <- sample(1:area,n,replace = T)
  claims <- matrix(0,nrow = n,ncol=years+1)
  mu <- matrix(0,nrow=n,ncol=years)
  phi <- matrix(0,nrow=n,ncol=years)
  z <- matrix(0,nrow=n,ncol=years)
  years_vec = rep(1:years,each=n)
  locs_vec = c()
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),rnorm(n*years,0.1,sd = 0.1),rnorm(n*years,-0.3,sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1,-1,1)
  
  
  
  
  prop_true <- 0.7
  p <- length(unique(locs))
  # Generate  data
  for(t in 1:years){
    agg_sds <- tapply(claims[,t],locs,sd)
    
    mu[,t] <- exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true)*exposure[((t-1)*n +1): (t*n)]
    
    if(is.null(beta2_true)){
      print("no beta")
      z[,t] <- 1
    }else{
      phi[,t] <- exp(X_mat2[((t-1)*n +1): (t*n),, drop = FALSE] %*% beta2_true)
      if(type == "PG"){
        z[,t] <- rgamma(n,phi[,t],phi[,t]) 
      }else if (type == "PIG"){
        z[,t] <- actuar::rinvgauss(n,mean = 1, shape = phi[,t]^2)
      }else if (type == "PLN"){
        z[,t] <- rlnorm(n,meanlog = -phi[,t]^2/2, sdlog = phi[,t])
      }else{
        stop()
      }
    }
    u <- runif(n)
    
    if(dist == 'ordinary'){
      claims[,t+1] <- rpois(n, lambda = z[,t]*mu[, t])
    }else if(dist == 'zip'){

      claims[,t+1] <- (u>=prop_true)* rpois(n, mu[,t]* z[,t])

      
    }else if(dist == 'hurdle'){
      # 
      
      claims[,t+1] <- ((u>=prop_true)*qpois(runif(n, dpois(0,mu[,t]* z[,t])+1e-8, 1), mu[,t]* z[,t]))
    }else{
      stop("ekkert")
    }
    
    locs_vec <- c(locs_vec,locs)
  }
  
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  
  agg_claims <- tapply(claims,list(locs_vec,years_vec),sum)/tapply(exposure,list(locs_vec,years_vec),sum)
  agg_claims <- cbind(rep(0.1, area), agg_claims)
  agg_claims
  
  
  return(list(claims = claims, z = z, locs = locs_vec, years = years_vec, beta = beta1_true,
              X1=X_mat1, X2=X_mat2, prop = prop_true, agg_claims = agg_claims, exposure = exposure ) )
  
}

generate_psi_add <- function(beta2_true, type = NULL, dist = 'ordinary', n = 100){
  n <- n
  seed <- 10
  years <- 15  # number of time points
  area <- 5  # number of regions
  set.seed(seed)
  exposure <- rpois(n, 10)
  exposure[exposure == 0] <- 1
  exposure <- rep(exposure, years)
  
  # store claim data
  locs <- sample(1:area,n,replace = T)
  claims <- matrix(0,nrow = n,ncol=years+1)
  mu <- matrix(0,nrow=n,ncol=years)
  phi <- matrix(0,nrow=n,ncol=years)
  z <- matrix(0,nrow=n,ncol=years)
  years_vec = rep(1:years,each=n)
  locs_vec = c()
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent< 0.7] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),rnorm(n*years,0.1,sd = 0.1),rnorm(n*years,-0.3,sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1,-1,1)
  
  tmp <- runif(area,0,2)
  w <- outer(tmp, tmp)
  w <- 1/(1+w)
  diag(w) <- 1
  
  
  psi1_true <- runif(area, 0.2, 0.5)
  
  prop_true <- 0.7
  p <- length(unique(locs))
  # Generate  data
  for(t in 1:years){
    
    mu[,t] <- (exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) + (psi1_true * (w %*% y_latent[,t]))[locs])*exposure[((t-1)*n +1): (t*n)]
    
    if(is.null(beta2_true)){
      z[,t] <- 1
    }else{
      phi[,t] <- exp(X_mat2[((t-1)*n +1): (t*n),, drop = FALSE] %*% beta2_true)
      if(type == "PG"){
        z[,t] <- rgamma(n, shape = phi[,t], rate = phi[,t]) 
      }else if (type == "PIG"){
        z[,t] <- actuar::rinvgauss(n,mean = 1, shape = phi[,t]^2)
      }else if (type == "PLN"){
        z[,t] <- rlnorm(n,meanlog = -phi[,t]^2/2, sdlog = phi[,t])
      }
    }
    
    
    
    u <- runif(n)
    if(dist == 'ordinary'){
      claims[,t+1] <- rpois(n, lambda = z[,t]*mu[, t])
    }else if(dist == 'zip'){

      claims[,t+1] <- (u>=prop_true)* rpois(n, mu[,t]* z[,t])

      
    }else if(dist == 'hurdle'){
      claims[,t+1] <- ((u>=prop_true)*qpois(runif(n, dpois(0,mu[,t]* z[,t])+1e-5, 1), mu[,t]* z[,t]))
    }else{
      stop("ekkert")
    }
    locs_vec <- c(locs_vec,locs)
  }
  
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  
  agg_claims <- tapply(claims,list(locs_vec,years_vec),sum)/tapply(exposure,list(locs_vec,years_vec),sum)
  agg_claims <- cbind(rep(0.1, area), agg_claims)
  agg_claims
  
  
  
  return(list(agg_claims = agg_claims , claims = claims, z = z, locs = locs_vec, years = years_vec, 
              X1=X_mat1, X2=X_mat2, w = w, beta1 = beta1_true, 
              beta2_true = beta2_true, psi = psi1_true, prop = prop_true,
              y_latent = y_latent, exposure = exposure) )
  
}

generate_psi_mult <- function(beta2_true, type = NULL, dist = 'ordinary', n = 100){
  n <- n
  seed <- 10
  years <- 8  # number of time points
  area <- 5  # number of regions
  set.seed(seed)
  exposure <- rpois(n, 10)
  exposure <- rep(exposure, years)
  
  # store claim data
  locs <- sample(1:area,n,replace = T)
  claims <- matrix(0,nrow = n,ncol=years+1)
  mu <- matrix(0,nrow=n,ncol=years)
  phi <- matrix(0,nrow=n,ncol=years)
  z <- matrix(0,nrow=n,ncol=years)
  years_vec = rep(1:years,each=n)
  locs_vec = c()
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent< 0.7] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),rnorm(n*years,0.1,sd = 0.1),rnorm(n*years,-0.3,sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1,-1,1)
  
  tmp <- runif(area,0,2)
  w <- outer(tmp, tmp)
  w <- 1/(1+w)
  diag(w) <- 1
  
  
  psi1_true <- runif(area, 0.01, 0.1)
  
  prop_true <- 0.7
  p <- length(unique(locs))
  # Generate  data
  for(t in 1:years){
    mu[,t] <- (exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) * (1+ (psi1_true * (w %*% y_latent[,t]))[locs]))*exposure[((t-1)*n +1): (t*n)]
    
    if(is.null(beta2_true)){
      z[,t] <- 1
    }else{
      phi[,t] <- exp(X_mat2[((t-1)*n +1): (t*n),, drop = FALSE] %*% beta2_true)
      if(type == "PG"){
        z[,t] <- rgamma(n,phi[,t],phi[,t]) 
      }else if (type == "PIG"){
        z[,t] <- actuar::rinvgauss(n,mean = 1, shape = phi[,t]^2)
      }else if (type == "PLN"){
        z[,t] <- rlnorm(n,meanlog = -phi[,t]^2/2, sdlog = phi[,t])
      }
    }
    
    
    
    u <- runif(n)
    if(dist == 'ordinary'){
      claims[,t+1] <- rpois(n, lambda = z[,t]*mu[, t])
    }else if(dist == 'zip'){
      claims[,t+1] <- (u>=prop_true)*rpois(n, mu[,t]* z[,t])
    }else if(dist == 'hurdle'){
      claims[,t+1] <- ((u>=prop_true)*qpois(runif(n, dpois(0,mu[,t]* z[,t])+1e-5, 1), mu[,t]* z[,t]))
    }else{
      stop("ekkert")
    }
    locs_vec <- c(locs_vec,locs)
  }
  
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  
  agg_claims <- tapply(claims,list(locs_vec,years_vec),sum)/tapply(exposure,list(locs_vec,years_vec),sum)
  agg_claims <- cbind(rep(0.1, area), agg_claims)
  
  return(list(agg_claims = agg_claims , claims = claims, z = z, locs = locs_vec, years = years_vec, beta1 = beta1_true, psi =psi1_true, 
              X1=X_mat1, X2=X_mat2, w = w, y_latent = y_latent, exposure = exposure, psi = psi1_true, prop = prop_true ) )
  
}

generate_a_add <- function(beta2_true, type = NULL, dist = 'ordinary', n = 100){
  n <- n
  seed <- 10
  years <- 8  # number of time points
  area <- 5  # number of regions
  set.seed(seed)
  
  # store claim data
  locs <- sample(1:area,n,replace = T)
  claims <- matrix(0,nrow = n,ncol=years+1)
  mu <- matrix(0,nrow=n,ncol=years)
  phi <- matrix(0,nrow=n,ncol=years)
  z <- matrix(0,nrow=n,ncol=years)
  years_vec = rep(1:years,each=n)
  locs_vec = c()
  
  exposure <- rpois(n, 10)
  exposure <- rep(exposure, years)
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent< 0.7] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),rnorm(n*years,0.1,sd = 0.1),rnorm(n*years,-0.3,sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1,-1,1)
  
  A1_true <- outer(0.2*runif(area,0,2),0.2*runif(area,0,2))
  mask_true = matrix(rbinom(area*area,size = 1, 0.3), nrow = area, ncol = area)
  diag(mask_true) = T
  A1_true = A1_true*mask_true
  A1_true = A1_true + t(A1_true)
  diag(A1_true) = 0.1
  
  
  prop_true <- 0.7
  p <- length(unique(locs))
  # Generate  data
  for(t in 1:years){
    
    mu[,t] <- (exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) + (1 * (A1_true %*% y_latent[,t]))[locs])*exposure[((t-1)*n +1): (t*n)]
    
    if(is.null(beta2_true)){
      z[,t] <- 1
    }else{
      phi[,t] <- exp(X_mat2[((t-1)*n +1): (t*n),, drop = FALSE] %*% beta2_true)
      if(type == "PG"){
        z[,t] <- rgamma(n,phi[,t],phi[,t]) 
      }else if (type == "PIG"){
        z[,t] <- actuar::rinvgauss(n,mean = 1, shape = phi[,t]^2)
      }else if (type == "PLN"){
        z[,t] <- rlnorm(n,meanlog = -phi[,t]^2/2, sdlog = phi[,t])
      }
    }
    
    
    u <- runif(n)
    if(dist == 'ordinary'){
      claims[,t+1] <- rpois(n, lambda = z[,t]*mu[, t])
    }else if(dist == 'zip'){
      claims[,t+1] <- (u>=prop_true)*rpois(n, mu[,t]* z[,t])
    }else if(dist == 'hurdle'){
      claims[,t+1] <- ((u>=prop_true)*qpois(runif(n, dpois(0,mu[,t]* z[,t])+1e-5, 1), mu[,t]* z[,t]))
    }else{
      stop("ekkert")
    }
    locs_vec <- c(locs_vec,locs)
  }
  
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  
  agg_claims <- tapply(claims,list(locs_vec,years_vec),sum)/tapply(exposure,list(locs_vec,years_vec),sum)
  
  
  a <- A1_true[upper.tri(A1_true, diag = T)]
  
  return(list(agg_claims = agg_claims , claims = claims, z = z, locs = locs_vec, years = years_vec, beta1 = beta1_true, beta2 = beta2_true,
              X1=X_mat1, X2=X_mat2, A = A1_true, prop = prop_true, y_latent = y_latent, exposure = exposure, a = a ) )
  
}

generate_a_mult <- function(beta2_true, type = NULL, dist = 'ordinary', n = 100){
  n <- n
  seed <- 10
  years <- 8  # number of time points
  area <- 5  # number of regions
  set.seed(seed)
  
  # store claim data
  locs <- sample(1:area,n,replace = T)
  claims <- matrix(0,nrow = n,ncol=years+1)
  mu <- matrix(0,nrow=n,ncol=years)
  phi <- matrix(0,nrow=n,ncol=years)
  z <- matrix(0,nrow=n,ncol=years)
  years_vec = rep(1:years,each=n)
  locs_vec = c()
  
  exposure <-  rpois(n, 10)
  exposure <- rep(exposure, years)
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent< 0.7] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),rnorm(n*years,0.1,sd = 0.1),rnorm(n*years,-0.3,sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1,-1,1)
  
  A1_true <- outer(0.3*runif(area,0,2),0.3*runif(area,0,2))
  mask_true = matrix(rbinom(area*area,size = 1, 0.5), nrow = area, ncol = area)
  diag(mask_true) = T
  A1_true = A1_true*mask_true
  A1_true = A1_true + t(A1_true)
  diag(A1_true) = 0.1
  
  prop_true <- 0.3
  p <- length(unique(locs))
  # Generate  data
  for(t in 1:years){
    
    
    mu[,t] <- (exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) * (1 + (A1_true %*% y_latent[,t]))[locs] ) *exposure[((t-1)*n +1): (t*n)]
    
    if(is.null(beta2_true)){
      z[,t] <- 1
    }else{
      phi[,t] <- exp(X_mat2[((t-1)*n +1): (t*n),, drop = FALSE] %*% beta2_true)
      if(type == "PG"){
        z[,t] <- rgamma(n,phi[,t],phi[,t]) 
      }else if (type == "PIG"){
        z[,t] <- actuar::rinvgauss(n,mean = 1, shape = phi[,t]^2)
      }else if (type == "PLN"){
        z[,t] <- rlnorm(n,meanlog = -phi[,t]^2/2, sdlog = phi[,t])
      }
    }
    
    
    
    u <- runif(n)
    if(dist == 'ordinary'){
      claims[,t+1] <- rpois(n, lambda = z[,t]*mu[, t])
    }else if(dist == 'zip'){
      claims[,t+1] <- (u>=prop_true)*rpois(n, mu[,t]* z[,t])
    }else if(dist == 'hurdle'){
      claims[,t+1] <- ((u>=prop_true)*qpois(runif(n, dpois(0,mu[,t]* z[,t])+1e-5, 1), mu[,t]* z[,t]))
    }else{
      stop("ekkert")
    }
    locs_vec <- c(locs_vec,locs)
  }
  
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  
  agg_claims <- tapply(claims,list(locs_vec,years_vec),sum)/tapply(exposure,list(locs_vec,years_vec),sum)
  
  
  a <- A1_true[upper.tri(A1_true, diag = T)]
  
  return(list(agg_claims = agg_claims , claims = claims, z = z, locs = locs_vec, years = years_vec, beta1 = beta1_true, beta2 = beta2_true,
              X1=X_mat1, X2=X_mat2, A = A1_true, prop = prop_true, y_latent = y_latent, exposure = exposure, a=a ) )
  
}










