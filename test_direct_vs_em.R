
# Zero inflated Poisson with Gamma mixing - direct ML maximization

source("utils.R")
source("zip.R")
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
  
  fmu_return[is.na(fmu_return)] <-0
  fphi_return[is.na(fphi_return)] <-0
  fprop[is.na(fprop)] <-0
  
  
  
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
  
  
  se1 <- get_spatial_aggregate(locs, w, psi, agg_claims, years, "learn_psi")
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
  
  
  se1 <- get_spatial_aggregate(locs, w, psi, agg_claims, years, "learn_psi")
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
  
  
  
  
  se1 <- get_spatial_aggregate(locs, A, 1, agg_claims, years, "learn_graph")
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
  
  
  se1 <- get_spatial_aggregate(locs, A, 1, agg_claims, years, "learn_graph")
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
  
  
  return(list(beta1 = beta,beta2 = beta2, prop = prop, mu = mu, out=out, Hessian = out$hessian))
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
  
  
  
  
  se1 <- get_spatial_aggregate(locs, w, psi, agg_claims, years, "learn_psi")
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta,beta2 = beta2, psi = psi, prop = prop, mu = mu, out=out, Hessian = out$hessian))
}

ZINB_a_all_in_one <- function(claims, X1, locs, years, agg_claims, additive,  X2 = matrix(1, nrow = length(claims)), exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0){
  
  
  p <- length(unique(locs))
  
  # Set initial parameters
  out_zip <- zip(claims, X1, locs, years,  agg_claims,
                 A, additive, "learn_graph", lambda = lambda, exposure = exposure, max_itr = 300)
  
  
  beta <- out_zip$beta1
  a <- out_zip$a
  prop <- out_zip$prop
  
  
  beta2 <- c(0, rep(0, ncol(X2) - 1)) 
  
  
  lower_beta1 <- rep(-20, ncol(X1))
  upper_beta1 <- rep(Inf, ncol(X1))
  
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
              upper = c(0.99, rep(Inf, ncol(X2)), upper_beta1, rep(Inf, p*(p+1)/2) ),
              hessian = T)
  
  
  prop <- out$par[1]
  beta2 <- out$par[2:(ncol(X2)+1)]
  beta <- out$par[(ncol(X2)+2):(ncol(X1)+ ncol(X2)+1)]
  a <- out$par[(ncol(X1)+ ncol(X2)+2):length(out$par)]
  A <- get_W_from_array(a, p)
  
  
  
  
  se1 <- get_spatial_aggregate(locs, A, 1, agg_claims, years, "learn_graph")
  nu1 <- exp(X1 %*% as.matrix(beta))
  nu2 <- exp(X2 %*% as.matrix(beta2))
  if(additive){
    mu <- (nu1 + se1$spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu1 * (1+se1$spatial_effect))*exp(log(exposure))
  }
  
  
  return(list(beta1 = beta, beta2 = beta2, a = a, A=A, prop = prop, mu = mu, out=out, Hessian = out$hessian))
}

# source("simulate_data.R")
# source("zip.R")
# n <- 50
# t <- 50
# data_sim <- simulate_claims(n, t, "graph", FALSE, mixing = "gamma", model_type = "zip",  exposure_lambda = 0, area = 10, seed = 1)
# 
# # Extract variables from simulation
# claims <- data_sim$claims
# X <- data_sim$X
# years <- data_sim$years
# locs <- data_sim$locs
# agg_claims <- data_sim$agg_claims
# A <- data_sim$A
# exposure <- data_sim$exposure
# model_type <- "learn_graph"
# additive <- FALSE
# mixing <- "gamma"
# 
# 
# source("utils.R")
# source("zip_mixed.R")
# out_zip <- zip(claims, X, locs, years,  agg_claims,
#                A, additive, model_type, lambda = 0, exposure = exposure, max_itr = 300)
# 
# out_mzip_direct <- ZINB_a_all_in_one(claims, X, locs, years, agg_claims, additive,  X2 = matrix(1, nrow = length(claims)), exposure = rep(1, length(claims)), max_itr = 1000, lambda = 0)
# 
# out <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, model_type, additive, mixing,  Emethod = "integration",
#               n_iter = 10, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
#               optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
#               batch_size = 100, verbose = 2, do_optim = FALSE, calc_se = FALSE)



source("simulate_data.R")
do_zip_t_direct <- function(density, id){
  
  
  A_est_ig <- list()
  A_ig_info <- list()
  beta_est_ig <- list()
  beta2_est_ig <- list()
  H_ig_p <- list()
  var_ig_p <- list()
  time_ig <- list()
  
  additive <- TRUE
  nr_regions <- 10
  sim <- simulate_claims(50, 5000, spatial_type = "graph",additive =  additive, area = nr_regions, 
                         model_type = "zip", mixing = "gamma", density = density,  seed = id)
  A <- sim$A
  
  
  beta_est_known_ig <- list()
  beta2_est_known_ig <- list()
  
  ts <- c(10, 20, 50, 100, 300, 500, 1000)#
  
  for(t in ts){
    
    
    print(t)
    # subset claims
    claims <- sim$claims[sim$years <= t]
    locs <- sim$locs[sim$years <= t]
    agg_claims <- sim$agg_claims[, 1:t]
    X <- sim$X[sim$years <= t, ]
    years <- sim$years[sim$years <= t]
    exposure <- sim$exposure[sim$years <= t]
    
    
    
    

    
    # Mixed Poisson
    a1 <- Sys.time()
    
    out_ig <-  ZINB_a_all_in_one(claims, X, locs, years, agg_claims, additive,  
                                 X2 = matrix(1, nrow = length(claims)), 
                                 exposure = exposure, max_itr = 1000, lambda = 0)
    
    
    time_ig[[as.character(t)]] <-  Sys.time() - a1
    
    A_est_ig[[as.character(t)]] <- out_ig$a
    A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
    beta_est_ig[[as.character(t)]] <- out_ig$beta1
    H_ig_p[[as.character(t)]] <- out_ig$Hessian
    beta2_est_ig[[as.character(t)]] <- out_ig$beta2
  
    save(
      list = ls(envir = environment()),      # all names in this env
      file =  paste0("sim_data/zip_direct_gamma_per_t",round(density, 2), "id", id , ".RData" ),               # where to write
      envir = environment()                  # whose vars to save
    )
    
  }
  
}


for(id in 1:10){
  do_zip_t_direct(0.4, id)
  do_zip_t_direct(0.8, id)
}


source("zip_mixed.R")
do_zip_t <- function(density, id, calc_se){
  
  A_est_p <- list()
  A_p_info <- list()
  beta_est_p <- list()
  H_p <- list()
  time_p <- list()
  
  A_est_ig <- list()
  A_ig_info <- list()
  beta_est_ig <- list()
  beta2_est_ig <- list()
  H_ig_p <- list()
  var_ig_p <- list()
  time_ig <- list()
  
  additive <- TRUE
  nr_regions <- 10
  sim <- simulate_claims(50, 5000, spatial_type = "graph",additive =  additive, area = nr_regions, 
                         model_type = "zip", mixing = "gamma", density = density,  seed = id)
  A <- sim$A
  
  
  beta_est_known_ig <- list()
  beta2_est_known_ig <- list()
  
  ts <- c(10, 20, 50, 100, 300, 500, 1000)#
  
  for(t in ts){
    
    print(t)
    # subset claims
    claims <- sim$claims[sim$years <= t]
    locs <- sim$locs[sim$years <= t]
    agg_claims <- sim$agg_claims[, 1:t]
    X <- sim$X[sim$years <= t, ]
    years <- sim$years[sim$years <= t]
    exposure <- sim$exposure[sim$years <= t]
    
    
    
    
    # Poisson
    a1 <- Sys.time()
    
    out_none <- zip(claims, X, locs, years,  agg_claims, 
                    A = A, additive, model_type = "learn_graph", lambda = 0, exposure = exposure, max_itr = 300, a_known = FALSE)
    
    time_p[[as.character(t)]] <-  Sys.time() - a1
    A_p_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
    
    A_est_p[[as.character(t)]] <- out_none$a
    beta_est_p[[as.character(t)]] <- out_none$beta1
    H_p[[as.character(t)]] <- out_none$optim_obj$hessian
    
    # Mixed Poisson
    a1 <- Sys.time()
    
    out_ig <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, 
                         model_type = "learn_graph", additive = additive, mixing = "gamma",  Emethod = "integration",
                         n_iter = 80, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                         optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                         batch_size = 100, param_tol = 1e-4, Q_tol = 1e-6, verbose = 2, do_optim = FALSE, a_known = FALSE, calc_se = calc_se)
    
    
    time_ig[[as.character(t)]] <-  Sys.time() - a1
    
    A_est_ig[[as.character(t)]] <- out_ig$a
    A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
    beta_est_ig[[as.character(t)]] <- out_ig$beta1
    H_ig_p[[as.character(t)]] <- out_ig$Hessian
    var_ig_p[[as.character(t)]] <- out_ig$var_loglik
    beta2_est_ig[[as.character(t)]] <- out_ig$beta2
    
    
    # Mixed Poisson known A
    out_ig_known <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, 
                               model_type = "learn_graph", additive = additive, mixing = "gamma",  Emethod = "integration",
                               n_iter = 80, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                               optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                               batch_size = 100, param_tol = 1e-5, Q_tol = 1e-9, verbose = 2, do_optim = FALSE, a_known = TRUE, calc_se =FALSE)
    
    
    
    beta_est_known_ig[[as.character(t)]] <- out_ig_known$beta1
    beta2_est_known_ig[[as.character(t)]] <- out_ig_known$beta2
    
    
    save(
      list = ls(envir = environment()),      # all names in this env
      file =  paste0("sim_data/zip_gamma_per_t",round(density, 2), "id", id , ".RData" ),               # where to write
      envir = environment()                  # whose vars to save
    )
    
  }
  
}


for(id in 1:10){
  print(id)
  do_zip_t(0.4, id, FALSE)
  do_zip_t(0.8, id, FALSE)
}







