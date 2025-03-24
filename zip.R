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


 

Q_ZIP_deriv <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, additive, model_type, lambda, p){
  
  etheta <- as.numeric(etheta)
  
  prop <- param[length(param)]
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param)-1)]
    a <- 1
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param)-1)]
    A <- get_W_from_array(a, p)
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
    
    grad <- c(g1, g22)
  }
  
  
  return(-c(grad, g3))  # adding prop deriv
  
  
}



Q_ZIP <- function(param, locs, claims, exposure, X_mat, agg_claims, years, etheta, A, 
                  additive, model_type, lambda, p, mu_old){
  
  
  etheta <- as.numeric(etheta)
  
  prop <- param[length(param)]
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param)-1)]
    a <- 1
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param)-1)]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  ll <- sum(log_dZIP(claims, mu, prop, etheta))
  
  return(-ll)
  
}



 

zip <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure,  lambda = 0, nr_em = 100, max_itr = 1000){
  
  
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
    theta0 <- c(beta1, 0.3)
    lower <- c(rep(-20,ncol(X1)), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)), 1-1e-3)
    psi <- NA
    A <- NA
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi, 0.3)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)+p), 1-1e-3)
  }else if(model_type == "learn_graph"){
    theta0 <- c(beta1, a, 0.3)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2), 1e-3)
    upper <-  c(rep(Inf,ncol(X1)+p*(p+1)/2), 1-1e-3)
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
               method = 'L-BFGS-B',
               control = list(maxit = max_itr),
               lower = lower,
               upper = upper,
               hessian = T)
  

  if(model_type == "ordinary"){
    beta1 <- out$par[1:ncol(X1)]
  }else if(model_type == "learn_psi"){
    beta1 <- out$par[1:ncol(X1)]
    psi <- out$par[(ncol(X1)+1):(length(out$par)-1)]
  }else if(model_type == "learn_graph"){
    beta1 <- out$par[1:ncol(X1)]
    a <- out$par[(ncol(X1)+1):(length(out$par)-1)]
    A <- get_W_from_array(a, p)
  }
  
  prop <- out$par[length(out$par)]
  
  return(list(beta1 = beta1, psi = psi, a = A[upper.tri(A, diag = TRUE)], 
              beta2 = 1 , H = out$hessian, H_beta2 = NA, mu = mu, optim_obj = out, 
              model_type = model_type, prop = prop))
  
  
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



Q_ZIP_mixed <- function(param, locs, claims, exposure, X_mat, agg_claims, years, A, additive, model_type, lambda, p, mu_old, phi, prop){
  

  
  prop <- param[length(param)]
  
  if(model_type == "ordinary"){
    beta1 <- param[1:ncol(X_mat)]
  }else if(model_type == "learn_psi"){
    beta1 <- param[1:ncol(X_mat)]
    psi <- param[(ncol(X_mat)+1):(length(param)-1)]
    a <- 1
  }else if(model_type == "learn_graph"){
    beta1 <- param[1:ncol(X_mat)]
    a <- param[(ncol(X_mat)+1):(length(param)-1)]
    A <- get_W_from_array(a, p)
    psi <- NA
  }
  
  se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
  nu <- as.numeric(exp(X_mat %*% beta1))
  
  mu <- get_mu(nu, se$spatial_effect, exposure, additive)
  
  Q_val <- E_fzy(y, mu,mu_old,prop, phi)
  
  ll <- sum(Q_val)
  
  return(-ll)
  
}



# Poisson mixture model 
zip_mixed <- function(claims, X1, locs, years, agg_claims, A, additive, model_type, exposure, lambda = 0, nr_em = 100, max_itr = 1000, mixing_var = "gamma", z= ""){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  out_zip <- zip(claims, X1, locs, years, agg_claims, A, additive, model_type, lambda = lambda, exposure, nr_em = 100, max_itr = 300)
  beta1 <- out_zip$beta1
  prop <- out_zip$prop
  
  
  if(model_type == "learn_graph"){
    a <- rep(1,p*(p+1)/2)*0.001
    A <- get_W_from_array(a, p)
    lower_a <- rep(1e-4, p*(p+1)/2)
    upper_a <- rep(0.5, p*(p+1)/2)
  }else if(model_type == "learn_psi"){
    psi <- rep(1,p)
  }
  
  
  
  
  # add gradients depending on model type
  if(model_type == "ordinary"){
    theta0 <- c(beta1, prop)
    lower <- c(rep(-20,ncol(X1)))
  }else if(model_type == "learn_psi"){
    theta0 <- c(beta1, psi, prop)
    lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
  }else if(model_type == "learn_graph"){
    theta0 <- c(beta1, a, prop)
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
  }else{
    stop("mixing type not known")
  }
  
  
  X2 <-  as.matrix(rep(1, nrow(X1)))
  
  
  for(i in 1:nr_em){
    
    
    se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
    nu <- as.numeric(exp(X1 %*% beta1))
    mu <-   get_mu(nu, se$spatial_effect, exposure, additive)
    phi <- exp(beta20)
    

 
    
    Q_ZIP_mixed(theta0, locs, claims, exposure, X1, agg_claims, years, A, additive, model_type, lambda, p, mu, phi, prop)
    
    ### zip part
    out <- optim(par = theta0,
                 fn = Q_ZIP_mixed,
                 #gr = Q_P_deriv,
                 locs = locs,
                 claims = claims, 
                 exposure = exposure, 
                 X_mat = X1, 
                 agg_claims = agg_claims, 
                 years = years,
                 A = A,
                 additive = additive,
                 model_type = model_type,
                 lambda = lambda,
                 p = p,
                 mu_old = mu,
                 phi = phi,
                 prop = prop,
                 method = 'L-BFGS-B',
                 control = list(maxit = ifelse(mixing_var == "none", 1000, 1)),
                 lower = lower,
                 hessian = T)
    
    
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
                      hessian = TRUE)
    
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









# ######### Simulate Data ################
# simulate_zip_claims <- function(n, years, spatial_type, additive, mixing, area = 5){
#   set.seed(123)
#   
#   
#   
#   # store claim data
#   locs <- sample(1:area,n,replace = T)
#   claims <- matrix(0,nrow = n,ncol=years+1)
#   mu <- matrix(0,nrow=n,ncol=years)
#   phi <- matrix(0,nrow=n,ncol=years)
#   z <- matrix(0,nrow=n,ncol=years)
#   years_vec = rep(1:years,each=n)
#   locs_vec = c()
#   
#   exposure <- 1#1rpois(n, 10)
#   exposure <- rep(1, years*n)
#   
#   y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
#   y_latent[y_latent< 0.4] <- 0
#   y_latent[y_latent != 0] <- 1
#   
#   
#   # Generate known parameters
#   X_mat1 <- cbind(rep(1, n*years),rnorm(n*years,0.1,sd = 0.1),rnorm(n*years,-0.3,sd = 0.1))
#   X_mat2 <- as.matrix(rep(1, n*years))
#   beta1_true <- c(1,-1,1)
#   
#   A1_true <- outer(0.2*runif(area,0,2),0.2*runif(area,0,2))
#   mask_true = matrix(rbinom(area*area,size = 1, 0.3), nrow = area, ncol = area)
#   diag(mask_true) = T
#   A1_true = A1_true*mask_true
#   A1_true = A1_true + t(A1_true)
#   diag(A1_true) <- 0.1
#   
#   psi_true <- runif(area, 1, 3)
#   
#   prop_true <- 0.7
#   p <- length(unique(locs))
#   beta2_true <- -0.5
#   # Generate  data
#   for(t in 1:years){
#     
#     if(spatial_type == "psi"){
#       if(additive){
#         mu[,t] <- (exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) + psi_true[locs]*(1 * (A1_true %*% y_latent[,t]))[locs])*exposure[((t-1)*n +1): (t*n)]
#       }else{
#         mu[,t] <- (1 + psi_true[locs]*(1 * (A1_true %*% y_latent[,t]))[locs] )*exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) *exposure[((t-1)*n +1): (t*n)]
#       }
#       
#     }else if(spatial_type == "graph") {
#       if(additive){
#         mu[,t] <- (exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) + (15 * (A1_true %*% y_latent[,t]))[locs])*exposure[((t-1)*n +1): (t*n)]
#       }else{
#         mu[,t] <- (1+(15 * (A1_true %*% y_latent[,t]))[locs])*exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true) *exposure[((t-1)*n +1): (t*n)]
#       }
#     }else if(spatial_type == "none") {
#         mu[,t] <- exp(X_mat1[((t-1)*n +1): (t*n),] %*% beta1_true)*exposure[((t-1)*n +1): (t*n)]
#     }
#     
#     
#     
#     
#     
#     
#     if(mixing == "none"){
#       z[,t] <- 1
#     }else if(mixing == "gamma"){
#       z[,t] <- rgamma(n,exp(beta2_true),exp(beta2_true))
#     }else if(mixing == "ig"){
#       z[,t] <- actuar::rinvgauss(n,mean = 1, shape = exp(beta2_true)^2)
#     }else if(mixing == "ln"){
#       z[,t] <- rlnorm(n,meanlog = -exp(beta2_true)^2/2, sdlog = exp(beta2_true))
#     }
#     
#     u <- runif(n)
#     claims[,t+1] <- (u>=prop_true)* rpois(n, mu[,t]* z[,t])
#     
#     locs_vec <- c(locs_vec,locs)
#   }
#   
#   
#   claims <- as.vector(claims[,-1])
#   z <- as.vector(z)
#   mu <- as.vector(mu)
#   phi <- as.vector(phi)
#   
#   
#   return(list(X=X_mat1, claims = claims, z = z, mu = mu, phi = phi, A = A1_true, psi = psi_true, beta2 = beta2_true, beta1 = beta1_true, agg_claims = y_latent,
#               years = years_vec, locs = locs_vec, exposure = exposure, a = A1_true[upper.tri(A1_true, TRUE)] ))
# }
# 
# 
# data_sim <- simulate_zip_claims(100, 20, "graph", FALSE, mixing = "none")
# 
# 
# out_zip <- zip(data_sim$claims, data_sim$X, data_sim$locs, data_sim$years,  data_sim$agg_claims, 
#                data_sim$A, TRUE, "learn_psi", lambda = 0, exposure = data_sim$exposure, max_itr = 300)
# out_zip$optim_obj$message
# 
# out_zip$psi
# 
# 
# out_zip$a
# 15*sim$a
# 
# #out_zip$psi
# #sim$psi
# 
# 
# 
# # 
# # sim <- simulate_zip_claims(100, 20, "graph", FALSE, mixing = "gamma")
# # claims <- sim$claims
# # X1 <- sim$X
# # years <- sim$years
# # locs <- sim$locs
# # agg_claims <- sim$agg_claims
# # A <- sim$A
# # exposure <- sim$exposure
# # model_type <- "learn_psi"
# # X_mat <- X1
# # additive <- TRUE
# # 
# 
# 
# 
# model_type = "learn_psi"
# mixing_var = "gamma"
# # Set initial parameters
# out_zip <- zip(claims, X, locs, years, agg_claims, A, additive, model_type, lambda = lambda, exposure, nr_em = 100, max_itr = 300)
# beta1 <- out_zip$beta1
# prop <- out_zip$prop
# out_zip$psi
# 
# p <- length(unique(locs))
# 
# 
# if(model_type == "learn_graph"){
#   a <- rep(1,p*(p+1)/2)*0.001
#   A <- get_W_from_array(a, p)
#   lower_a <- rep(1e-4, p*(p+1)/2)
#   upper_a <- rep(0.5, p*(p+1)/2)
# }else if(model_type == "learn_psi"){
#   psi <- rep(1,p)
# }
# 
# 
# 
# 
# # add gradients depending on model type
# if(model_type == "ordinary"){
#   theta0 <- c(beta1, prop)
#   lower <- c(rep(-20,ncol(X1)))
# }else if(model_type == "learn_psi"){
#   theta0 <- c(beta1, psi, prop)
#   lower <- c(rep(-20,ncol(X1)), rep(1e-8, p))
# }else if(model_type == "learn_graph"){
#   theta0 <- c(beta1, a, prop)
#   lower <- c(rep(-20,ncol(X1)), rep(1e-8, p*(p+1)/2))
#   psi <- NA
# }else{
#   stop("model type not known")
# }
# 
# # set initial for the mixing
# if(mixing_var == "gamma"){
#   beta20 <- 1
# }else if(mixing_var == "ig"){
#   beta20 <- 0
# }else if(mixing_var == "ln"){
#   beta20 <- 0
# }else{
#   stop("mixing type not known")
# }
# 
# 
# X2 <-  as.matrix(rep(1, nrow(X1)))
# 
# prop <- 0.7
# for(i in 1:10){
#   
#   
#   se <- get_spatial_aggregate(locs, A, psi, agg_claims, years, model_type)
#   nu <- as.numeric(exp(X1 %*% beta1))
#   mu_old <-   get_mu(nu, se$spatial_effect, exposure, additive)
#   phi <- exp(-0.5)
#   
#   
#   # do gradient
#   
#   log_dZIP_der_mu(claims, mu, mu_old, prop)
#   
#   int_fun <- function(z, mu, mu_old, y) grad_fun(z, y, mu, prop)*exp(log_dZIP(y,mu_old,prop, z))*dgamma(z,shape = phi, rate = phi )
#   
#   y_max <- max(y)
#   
#   out <- mapply(function(y, mu, mu_old){pracma::integral(int_fun, 1e-3, y_max*10, mu = mu, mu_old = mu_old, y = y )}, 
#                 y = y, mu = mu, mu_old = mu_old)
#   
#   
#   denominator <- fyzz(y,mu,prop, phi)
#   
#   
#   
#   
#   g1 <- t( out * nu * exp(log(exposure)) ) %*% X_mat
#   # psi
#   g2 <- out[,1]*exp(log(exposure)) * se$agg_effect
#   g2 <- tapply(g2,locs,sum)
#   
# 
#   
# }
# 
# 
# 
# 
# 
# 
# 
# 
