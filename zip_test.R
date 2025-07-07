# Script to test ZIP models both learn A and learn Psi


source("zip_mixed.R")
source("simulate_data.R")
#######################################
# The goal of this Document is to perform the tests above for the
# ordinary poisson models.



# Learn A with time

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
                           model_type = "zip", mixing = "ig", density = density,  seed = id)
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
                           model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
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
                                 model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
                                 n_iter = 80, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                 optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                 batch_size = 100, param_tol = 1e-5, Q_tol = 1e-9, verbose = 2, do_optim = FALSE, a_known = TRUE, calc_se =FALSE)
  
  
  
      beta_est_known_ig[[as.character(t)]] <- out_ig_known$beta1
      beta2_est_known_ig[[as.character(t)]] <- out_ig_known$beta2
  
      
      save(
        list = ls(envir = environment()),      # all names in this env
        file =  paste0("sim_data/zip_ln_per_t",round(density, 2), "id", id , ".RData" ),               # where to write
        envir = environment()                  # whose vars to save
      )
      
    }
    
  }


for(id in 1:10){
  print(id)
  if(id == 1){
    do_zip_t(0.2, id, TRUE)
    do_zip_t(0.4, id, TRUE)
    do_zip_t(0.8, id, TRUE)
  }else{
    do_zip_t(0.2, id, FALSE)
    do_zip_t(0.4, id, FALSE)
    do_zip_t(0.8, id, FALSE)
  }
}


# Learn A for different lambdas
do_zip_lambda <- function(density, id, calc_se = FALSE){
    
    A_est_p <- list()
    A_p_info <- list()
    beta_est_p <- list()
    
    A_est_ig <- list()
    A_ig_info <- list()
    beta_est_ig <- list()
    time_ig <- list()
    
    beta2_est <- list()
    
    
    additive <- FALSE
    nr_regions <- 10
    sim <- simulate_claims(50, 5000, spatial_type = "graph",additive =  additive, area = nr_regions, 
                           model_type = "zip", mixing = "ig", density = density, seed = id)
    A <- sim$A
    
    
    lambdas <- c(0, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7, 10, 25, 50)
    additive <- TRUE
    for(lambda in lambdas){
      for(t in c(10, 50, 100, 200)){
        
        
        # subset claims
        claims <- sim$claims[sim$years <= t]
        locs <- sim$locs[sim$years <= t]
        agg_claims <- sim$agg_claims[, 1:t]
        X <- sim$X[sim$years <= t, ]
        years <- sim$years[sim$years <= t]
        exposure <- sim$exposure[sim$years <= t]
        
        # Poisson
        out_none <- zip(claims, X, locs, years,  agg_claims,
                        A = A, additive, model_type = "learn_graph", lambda = lambda, exposure = exposure, max_itr = 300, a_known = FALSE)

        A_est_p[[paste0(t, "-", lambda)]] <- out_none$a
        beta_est_p[[paste0(t, "-", lambda)]] <- out_none$beta1
        A_p_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))

        # Mixed Poisson
        a1 <- Sys.time()
        out_ig <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, 
                             model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
                             n_iter = 80, lambda = lambda, optimizer_beta = "gd", optimizer_psi = "gd",
                             optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                             batch_size = 100, param_tol = 1e-4, Q_tol = 1e-6, verbose = 2, do_optim = FALSE, a_known = FALSE, calc_se = calc_se)

        
        time_ig[[paste0(t, "-", lambda)]] <-  Sys.time() - a1
        
        A_est_ig[[paste0(t, "-", lambda)]] <- out_ig$a
        A_ig_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
        beta_est_ig[[paste0(t, "-", lambda)]] <- out_ig$beta1
        beta2_est[[paste0(t, "-", lambda)]] <- out_ig$beta2
        
        
        save(
          list = ls(envir = environment()),      # all names in this env
          file =  paste0("sim_data/2zip_ig_per_t_lambda",round(density, 2), "id", id , ".RData" ),               # where to write
          envir = environment()                  # whose vars to save
        )
        
      }
      
    }  
  }
  
for(id in 1:10){
    print(id)
    do_zip_lambda(0.2, id, FALSE)
    do_zip_lambda(0.4, id, FALSE)
    do_zip_lambda(0.8, id, FALSE)
   
  }
  
  
  
# Learn psi
  do_zip_psi_t <- function(density, id, calc_se){
    
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
    sim <- simulate_claims(50, 5000, spatial_type = "psi",additive =  additive, area = nr_regions, 
                           model_type = "zip", mixing = "ig", density = density,  seed = id, agg_spar = 0.4)
    A <- sim$A
    print(sim$psi)
    
    
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
                      A = A, additive, model_type = "learn_psi", lambda = 0, exposure = exposure, max_itr = 300, a_known = FALSE)
      
      time_p[[as.character(t)]] <-  Sys.time() - a1
      
      A_est_p[[as.character(t)]] <- out_none$psi
      beta_est_p[[as.character(t)]] <- out_none$beta1
      H_p[[as.character(t)]] <- out_none$optim_obj$hessian
      
      # Mixed Poisson
      a1 <- Sys.time()
      
      out_ig <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, 
                           model_type = "learn_psi", additive = additive, mixing = "ig",  Emethod = "integration",
                           n_iter = 80, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                           optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                           batch_size = 100, param_tol = 1e-4, Q_tol = 1e-6, verbose = 2, do_optim = FALSE, a_known = FALSE, calc_se = calc_se,
                           control_list = list(beta_phi = list(lr = 3/nrow(X)),
                                               psi = list(lr = 1/t)))
      
      
      time_ig[[as.character(t)]] <-  Sys.time() - a1
      
      A_est_ig[[as.character(t)]] <- out_ig$psi
      beta_est_ig[[as.character(t)]] <- out_ig$beta1
      H_ig_p[[as.character(t)]] <- out_ig$Hessian
      var_ig_p[[as.character(t)]] <- out_ig$var_loglik
      beta2_est_ig[[as.character(t)]] <- out_ig$beta2
      
      

      
      
      save(
        list = ls(envir = environment()),      # all names in this env
        file =  paste0("sim_data/5zip_psi_ig_per_t",round(density, 2), "id", id , ".RData" ),               # where to write
        envir = environment()                  # whose vars to save
      )
      
    }
    
  }
  
  for(id in 1:10){
    print(id)
    if(id %in% c(1,2)){
      do_zip_psi_t(0.4, id, FALSE)
      do_zip_psi_t(0.8, id, FALSE)
    }else{
      do_zip_psi_t(0.4, id, FALSE)
      do_zip_psi_t(0.8, id, FALSE)
    }
    
  }
