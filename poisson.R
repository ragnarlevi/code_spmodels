#### This file contains Possion models. That is - no mixing and no EM needed.

set.seed(123)
source("utils.R")
source("simulate_data.R")


##########
# Define a Poisson model with no mixing - base
######
# Objective of Poisson

Q_P <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p){
  
  etheta <- as.numeric(etheta)
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):length(param)]
    a <- 1
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):length(param)]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  ll <- sum(claims*log(mu) -mu*etheta) - sum(a*lambda)*(model_type == "learn_graph")
  
  return(-ll)
  
}


# deriviative of the Poisson
Q_P_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p){
  

  etheta <- as.numeric(etheta)
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):length(param)]
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):length(param)]
    A <- get_W_from_array(a, p)
    psi <- NA
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
      g22 = g22 + t(G2)[upper.tri(G2, diag = T)] - 0
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
    
    grad <- c(g1, g22)
  }
  
  
  
  return(-grad)
  
}



# Poisson model
Poisson <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure,  lambda = 0, nr_em = 100, max_itr = 1000){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  
  
  if(model_type == "learn_graph"){
    a <- rep(1,p*(p+1)/2)*0.001
    A <- get_W_from_array(a, p)
    lower_a <- rep(1e-4, p*(p+1)/2)
    upper_a <- rep(0.5, p*(p+1)/2)
  }else if(model_type == "learn_psi"){
    psi <- rep(1,p)
  }
  
  
  
  
  # add initial param depending on model type
  if(model_type == "ordinary"){
    theta0 <- c(beta1)
    lower <- c(rep(-20,ncol(X1)))
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
  }else if(model_type == "learn_graph"){
    theta0 <- c(beta1, a)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2))
    psi <- NA
  }else{
    stop("model type not known")
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
                 method = 'L-BFGS-B',
                 control = list(maxit = max_itr),
                 lower = lower,
                 hessian = T)
    

  
  
  if(model_type == "ordinary"){
    beta1 <- out$par[1:ncol(X1)]
  }else if(model_type == "learn_psi"){
    beta1 <- out$par[1:ncol(X1)]
    psi <- out$par[(ncol(X1)+1):length(out$par)]
  }else if(model_type == "learn_graph"){
    beta1 <- out$par[1:ncol(X1)]
    a <- out$par[(ncol(X1)+1):length(out$par)]
    A <- get_W_from_array(a, p)
  }
  
  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], beta2 = 1 , H = out$hessian, H_beta2 = NA, mu = mu, optim_obj = out, model_type = model_type))
  
  
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
Q_PG_beta_2 <- function(param,  claims, exposure, X2,  etheta, eother){
  
  etheta = as.numeric(etheta)
  elogtheta = as.numeric(eother)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  
  ll = sum(phi*log(phi) -lgamma(phi) + (phi-1)*elogtheta - phi*etheta )
  return(-ll)
  
}

Q_PG_beta_2_deriv <- function(param,  claims, exposure, X2,  etheta, eother){
  
  etheta = as.numeric(etheta)
  elogtheta = as.numeric(eother)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  der1_logdga <- log(phi) + 1 - digamma(phi) + elogtheta - etheta
  
  
  g3 <- t( der1_logdga*nu2 ) %*% X2
  
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
Q_PIG_beta_2 <- function(param, claims, exposure, X2, etheta, eother){
  
  etheta = as.numeric(etheta)
  einvtheta = as.numeric(eother)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  
  ll = sum(log(phi) + phi^2 - 0.5*phi^2*(einvtheta + etheta))
  
  return(-ll)
  
}

Q_PIG_beta_2_deriv <- function(param, claims, exposure, X2, etheta, eother){
  
  etheta = as.numeric(etheta)
  einvtheta = as.numeric(eother)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2 
  
  der1_logdgauss <- 1/phi + 2*phi - phi*einvtheta - phi*etheta
  g3 <- t(der1_logdgauss* nu2)%*% X2
  
  return(-g3)
}







############# Poisson LN ##########




# Q function for the mixing log normal
Q_PLN_beta_2 <- function(param, claims, exposure, X2, etheta = 1, eother){
  
  elogthetasq = as.numeric(eother)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2
  
  
  ll = sum(-log(phi) - elogthetasq/(2*phi^2)  - phi^2/8)
  return(-ll)
  
}

Q_PLN_beta_2_deriv <- function(param, claims, exposure, X2, etheta = 1, eother){
  
  elogthetasq = as.numeric(eother)
  
  nu2 <- exp(X2 %*% as.matrix(param))
  phi <- nu2
  
  der1_logdig <- -1/phi + elogthetasq/(phi ^3) - 2*phi/8
  
  
  g3 <- t( der1_logdig*nu2 ) %*% X2
  
  return(-g3)
}



# Define function to calculate E-step numerically
expected_ord_w <- function(y, phi, mu, z_density, w_fun) {
  

    integrand_num <- function(z) {
      f_y_given_z <- dpois(y, z * mu)
      f_y_given_z * w_fun(z) * z_density(z, phi)
    }
    integrand_denom <- function(z) {
      f_y_given_z <- dpois(y, z * mu)
      f_y_given_z * z_density(z, phi)
    }
  
  
  num <- integrate(integrand_num, lower = 1e-6, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  denom <- integrate(integrand_denom, lower = 1e-6, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  if (denom == 0) return(0)
  return(num / denom)
}
expected_ord_w <- Vectorize(expected_ord_w, vectorize.args = c("y", "mu"))



# Poisson mixture model 
Poisson_mixed <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure, lambda = 0, nr_em = 100, max_itr = 1000, mixing_var = "gamma", z= ""){
  

  p <- length(unique(locs))
  
  # Set initial parameters
  out_poisson <- Poisson(claims, X1, locs, years, agg_claims, A, additive, model_type, lambda = lambda, exposure, nr_em = 100, max_itr = 1000)
  beta1 <- out_poisson$beta1
  
  
  if(model_type == "learn_graph"){
  a <- out_poisson$a
  A <- get_W_from_array(a, p)
  lower_a <- rep(1e-4, p*(p+1)/2)
  upper_a <- rep(0.5, p*(p+1)/2)
  }else if(model_type == "learn_psi"){
    psi <- out_poisson$psi
  }
  
  
  
  
  # add initial param depending on model type
  if(model_type == "ordinary"){
    theta0 <- c(beta1)
    lower <- c(rep(-20,ncol(X1)))
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
  }else if(model_type == "learn_graph"){
    theta0 <- c(beta1, a)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2))
    psi <- NA
  }else{
    stop("model type not known")
  }
  
  # set initial for the mixing
  if(mixing_var == "gamma"){
    beta20 <- 1
  }else if(mixing_var == "ig"){
    beta20 <- 0
  }else if(mixing_var == "ln"){
    beta20 <- 0
    phi <- exp(beta20)
    z_density <- ln_density
  }else{
    stop("mixing type not known")
  }
  X2 <-  as.matrix(rep(1, nrow(X1)))


  for(i in 1:nr_em){
    print(i)
    
    se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
    nu <- as.numeric(exp(X1 %*% beta1))
    mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
    phi <- exp(beta20)
    
    if(mixing_var == "gamma"){
        etheta <- Etheta(claims, mu, phi)
        eother <- Elogtheta(claims,mu, phi)
        Q_beta_2 <- Q_PG_beta_2
        Q_beta_2_deriv <- Q_PG_beta_2_deriv
      
    }else if(mixing_var == "ig"){
      etheta <- Etheta_gig(claims, phi^2 + 2*mu, phi^2)
      eother <- Einvtheta_gig(claims, phi^2 + 2*mu, phi^2)
      
      Q_beta_2 <- Q_PIG_beta_2
      Q_beta_2_deriv <- Q_PIG_beta_2_deriv
    }else if(mixing_var == "ln"){
      
      numer <- jpdf(claims, mu,phi)
      eother <- expected_ord_w(y = claims, phi = phi, mu = mu, z_density = z_density, w_fun = function(x) log(x)^2)
      etheta <- expected_ord_w(y = claims, phi = phi, mu = mu, z_density = z_density, w_fun = identity)
      
      Q_beta_2 <- Q_PLN_beta_2
      Q_beta_2_deriv <- Q_PLN_beta_2_deriv
    }
    
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
                 method = 'L-BFGS-B',
                 control = list(maxit = ifelse(mixing_var == "none", 1000, 2)),
                 lower = lower,
                 hessian = FALSE)
  
    
    ## mixing part
    
    out_beta2 = optim(beta20,
                       Q_beta_2,
                       gr = Q_beta_2_deriv,
                       claims = claims, 
                       exposure =exposure, 
                       X2 = X2, 
                       etheta = etheta,
                       eother = eother,
                       method = 'L-BFGS',
                       control = list(maxit = 2),
                       hessian = FALSE)
    
    beta2 <- out_beta2$par
    H_beta2 <- out_beta2$hessian
    

    
    if(model_type == "ordinary"){
      beta1 <- out$par[1:ncol(X1)]
      theta <- beta1
    }else if(model_type == "learn_psi"){
      beta1 <- out$par[1:ncol(X1)]
      psi <- out$par[(ncol(X1)+1):length(out$par)]
      theta <- c(beta1, psi)
    }else if(model_type == "learn_graph"){
      beta1 <- out$par[1:ncol(X1)]
      a <- out$par[(ncol(X1)+1):length(out$par)]
      A <- get_W_from_array(a, p)
      theta <- c(beta1, a)
    }
    
    # check convergence
    
    if(sum(abs(theta-theta0))/sum(abs(theta0)) < 1e-5 &  sum(abs(beta2-beta20))/sum(abs(beta20)) < 1e-5 ){
      break
    }else{
      theta0 <- theta
      beta20 <- beta2
    }
  
    
    
  }



  if(model_type == "ordinary"){
    beta1 <- out$par[1:ncol(X1)]
  }else if(model_type == "learn_psi"){
    beta1 <- out$par[1:ncol(X1)]
    psi <- out$par[(ncol(X1)+1):length(out$par)]
  }else if(model_type == "learn_graph"){
    beta1 <- out$par[1:ncol(X1)]
    a <- out$par[(ncol(X1)+1):length(out$par)]
    A <- get_W_from_array(a, p)
  }
  
  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], beta2 = beta2 , H = out$hessian, H_beta2 = H_beta2, mu = mu, optim_obj = out, model_type = model_type))
  
  
}

##### Test some data set ###########

sim <- simulate_claims(100, 20, "graph", TRUE, mixing = "ig", model_type = "zip")

out_poisson <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, TRUE, "learn_graph", lambda = 0, 
              exposure = sim$exposure, max_itr = 2, mixing_var = "ln", nr_em = 5)




######### Simulate Data ################

# A_est_p <- list()
# beta_est_p <- list()
# 
# A_est_ig <- list()
# beta_est_ig <- list()
# 
# ts <- c(10, 20, 50, 70, 100, 200, 500, 1000)
# for(t in ts){
#   print(t)
#   sim <- simulate_claims(100, t, FALSE, FALSE, mixing = "ig")
#   
#   out_none <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, FALSE, "learn_graph", lambda = 0, exposure = sim$exposure, max_itr = 100)
#   A_est_p[[t]] <- out_none$a
#   beta_est_p[[t]] <- out_none$beta1
#   
#   out_ig <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, FALSE, "learn_graph", lambda = 0, 
#                             exposure = sim$exposure, max_itr = 100, mixing_var = "ln")
#   A_est_ig[[t]] <- out_ig$a
#   beta_est_ig[[t]] <- out_ig$beta1
#   
#   
#   
# }
# 
# a_none_error <- lapply(A_est_p, function(x, y){sum(abs(x-y))}, y = 15*sim$A[upper.tri(sim$A, T)])
# a_ig_error <- lapply(A_est_ig, function(x, y){sum(abs(x-y))}, y = 15*sim$A[upper.tri(sim$A, T)])
# 
# 
# b_none_error <- lapply(beta_est_p, function(x, y){sum(abs(x-y))}, y = sim$beta1)
# b_ig_error <- lapply(beta_est_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)
# 
# 
# library(tidyverse)
# ggplot(data.frame(t = ts,
#                    poissson_error = unlist(a_none_error[ts]),
#                    ig_error = unlist(a_ig_error[ts])
#                      )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "IG"))
# 
# 
# 
# ggplot(data.frame(t = ts,
#                   poissson_error = unlist(b_none_error[ts]),
#                   ig_error = unlist(b_ig_error[ts])
# )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "IG"))
# 
# 
# 
# 
# 
# A_est_p_l <- list()
# beta_est_p_l <- list()
# 
# A_est_ig_l <- list()
# beta_est_ig_l <- list()
# 
# lambdas <- c(0, 0.01, 0.02, 0.025, 0.05, 0.1, 0.5, 1)*100*250
# for(lambda in lambdas){
#   print(t)
#   sim <- simulate_claims(100, 250, FALSE, FALSE, mixing = "ig")
#   
#   out_none <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, FALSE, 
#                       "learn_graph", lambda = 25000, exposure = sim$exposure, max_itr = 100)
#   A_est_p_l[[as.character(lambda)]] <- out_none$a
#   beta_est_p_l[[as.character(lambda)]] <- out_none$beta1
#   
#   out_ig <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, FALSE, 
#                           "learn_graph", lambda = lambda, 
#                           exposure = sim$exposure, max_itr = 100, mixing_var = "ig")
#   A_est_ig_l[[as.character(lambda)]] <- out_ig$a
#   beta_est_ig_l[[as.character(lambda)]] <- out_ig$beta1
#   
#   
#   
# }
# 
# a_none_nr_edgegs <- unlist(lapply(A_est_p, function(x){sum(abs(x>1e-3))   }  ))
# a_ig_nr_edgegs <- unlist(lapply(A_est_ig_l, function(x){sum(abs(x>1e-3))   }  ))
# 
# 
# 
# 
# library(tidyverse)
# ggplot(data.frame(lambda = lambdas,
#                   poissson_error = unlist(a_none_error[ts]),
#                   ig_error = unlist(a_ig_error[ts])
# )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "IG"))
# 
# 
# 
# ggplot(data.frame(t = ts,
#                   poissson_error = unlist(b_none_error[ts]),
#                   ig_error = unlist(b_ig_error[ts])
# )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "IG"))
# 
# 
# 
# save.image("Ordinary.RData")
# 











