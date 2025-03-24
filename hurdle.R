# Source code for hurdle model (no mixing)
source("utils.R")

# Objective of hurdle

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


Q_hurdle <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, prop){
  
  
  etheta <- as.numeric(etheta)
  
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param))]
    a <- 1
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param))]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  ll <- sum(log_dHP(claims, mu, prop, etheta))
  
  return(-ll)
  
}


Q_hurdle_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p, prop){
  
  etheta <- as.numeric(etheta)
  
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param))]
    a <- 1
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param))]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  
  derivs <-  log_dHP_der(claims, mu, prop, 1)
  fmu <- as.matrix(derivs$fmu)
  fphi_deriv <- derivs$fphi
  

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
    
    grad <- c(g1, g22)
  }
  
  
  return(-grad)  # adding prop deriv
  
  
}









hurdle <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure,  lambda = 0, nr_em = 100, max_itr = 1000){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  beta1 <- glm(claims ~ -1 + X1 + offset(log(exposure)), family = 'poisson')$coefficients
  prop <- mean(claims == 0)
  
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
    theta0 <- beta1
    lower <- c(rep(-20,ncol(X1)))
    upper <-  c(rep(Inf,ncol(X1)))
    psi <- NA
    A <- NA
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
    upper <-  c(rep(Inf,ncol(X1)+p))
  }else if(model_type == "learn_graph"){
    theta0 <- c(beta1, a)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2))
    upper <-  c(rep(Inf,ncol(X1)+p*(p+1)/2))
    psi <- NA
  }else{
    stop("model type not known")
  }
  
  
  
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X1 %*% beta1))
  mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
  
  
  etheta <- 1
  
  ### hurdle poisson part
  out <- optim(par = theta0,
               fn = Q_hurdle,
               gr = Q_hurdle_deriv,
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
               prop = prop,
               method = 'L-BFGS-B',
               control = list(maxit = max_itr),
               lower = lower,
               upper = upper,
               hessian = T)
  
  
  if(model_type == "ordinary"){
    beta1 <- out$par[1:ncol(X1)]
  }else if(model_type == "learn_psi"){
    beta1 <- out$par[1:ncol(X1)]
    psi <- out$par[(ncol(X1)+1):(length(out$par))]
  }else if(model_type == "learn_graph"){
    beta1 <- out$par[1:ncol(X1)]
    a <- out$par[(ncol(X1)+1):(length(out$par))]
    A <- get_W_from_array(a, p)
  }

  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], 
              beta2 = 1 , H = out$hessian, H_beta2 = NA, mu = mu, optim_obj = out, 
              model_type = model_type, prop = prop))
  
  
}



source("simulate_data.R")

data_sim <- simulate_hurdle_claims(100, 200, "psi", TRUE, mixing = "none")

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
mixing <- "none"


# Initialize parameters
out_hurdle <- hurdle(data_sim$claims, data_sim$X, data_sim$locs, data_sim$years,  data_sim$agg_claims, 
               data_sim$A, additive, model_type, exposure = data_sim$exposure,  lambda = 0, max_itr = 300)

out_hurdle$optim_obj$message

data_sim$psi
out_hurdle$psi

15*data_sim$a
out_hurdle$a



