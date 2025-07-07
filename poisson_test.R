# We need to test
# 1. Do we uncover true A.
# 2. Assume A known then assume unkown, does beta change?
# 3. Look at how fit behaves with n and t. YES
# 4. Computational complexity with number of parameter and data points.
# 5. How do standard errors converge. Plot standard errors of estimates
# 6. Conway Maxwell Poisson (COMP) distributions



source("poisson.R")
source("simulate_data.R")
#######################################
# The goal of this Document is to perform the tests above for the
# ordinary poisson models.




##########
######### Test behaviour with t ############
##########

do_poisson_t <- function(density, id){
  
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
  
  additive <- FALSE
  nr_regions <- 10
  sim <- simulate_claims(50, 5000, spatial_type = "graph", additive =  additive, area = nr_regions, 
                         model_type = "poisson", mixing = "ln", density = density, seed = id)
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
    out_none <- Poisson(claims, X, locs, years, agg_claims, NA, additive = additive, model_type = "learn_graph", 
                        lambda = 0, exposure = exposure, max_itr = 100)
    
    time_p[[as.character(t)]] <-  Sys.time() - a1
    A_p_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
    
    A_est_p[[as.character(t)]] <- out_none$a
    beta_est_p[[as.character(t)]] <- out_none$beta1
    H_p[[as.character(t)]] <- out_none$optim_obj$hessian
    
    # Mixed Poisson
    a1 <- Sys.time()
    out_ig <- Poisson_mixed(claims, X, locs, years, agg_claims, A = NA, additive = additive, model_type = "learn_graph", 
                            lambda = 0,  Q_tol = 0,
                            exposure = exposure, mixing_var = "ln", nr_em = 100, verbose = 2)
    
    time_ig[[as.character(t)]] <-  Sys.time() - a1
    
    A_est_ig[[as.character(t)]] <- out_ig$a
    A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
    beta_est_ig[[as.character(t)]] <- out_ig$beta1
    H_ig_p[[as.character(t)]] <- out_ig$Hessian
    var_ig_p[[as.character(t)]] <- out_ig$var_loglik
    beta2_est_ig[[as.character(t)]] <- out_ig$beta2
    
    
    # Mixed Poisson known A
    out_ig_known <- Poisson_mixed(claims, X, locs, years, agg_claims, A = A, additive = additive, model_type = "learn_graph", 
                                  lambda = 0, 
                                  exposure = exposure, max_itr = 5, mixing_var = "ln", nr_em = 100, verbose = 2, a_known = TRUE)
    
    
    
    beta_est_known_ig[[as.character(t)]] <- out_ig_known$beta1
    beta2_est_known_ig[[as.character(t)]] <- out_ig_known$beta2
    

    save(
      list = ls(envir = environment()),      # all names in this env
      file =  paste0("Poisson_ln_per_t",round(density, 2), "id", id , ".RData" ),               # where to write
      envir = environment()                  # whose vars to save
    )
    
  }
  
}

# Learn A with time
for(id in 1:10){
  print(id)
  do_poisson_t(0.4, id)
  do_poisson_t(0.8, id)
  
}


##########
######### Test behaviour with lambda ############
##########


# Learn A for different lambdas
do_poisson_lambda <- function(density, id){
  
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
                         model_type = "poisson", mixing = "ln", density = density, seed = id)
  A <- sim$A
  
  
  lambdas <- c(0, 0.1, 0.25, 0.5, 0.75, 1, 2.5, 5, 7, 10, 25, 50, 75, 100)
  additive <- FALSE
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
      out_none <- Poisson(claims, X, locs, years, agg_claims, NA, additive = additive, "learn_graph", 
                          lambda = lambda, exposure = exposure, max_itr = 100)
      A_p_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(A, get_W_from_array(out_none$a, 10))
      
      A_est_p[[paste0(t, "-", lambda)]] <- out_none$a
      beta_est_p[[paste0(t, "-", lambda)]] <- out_none$beta1
      
      # Mixed Poisson
      a1 <- Sys.time()
      out_ig <- Poisson_mixed(claims, X, locs, years, agg_claims, A = NA, additive = additive, model_type = "learn_graph", 
                              lambda = lambda, 
                              exposure = exposure, max_itr = 10, mixing_var = "ln", nr_em = 50, verbose = 2)
      time_ig[[paste0(t, "-", lambda)]] <-  Sys.time() - a1
      
      A_est_ig[[paste0(t, "-", lambda)]] <- out_ig$a
      A_ig_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
      beta_est_ig[[paste0(t, "-", lambda)]] <- out_ig$beta1
      beta2_est[[paste0(t, "-", lambda)]] <- out_ig$beta2
      
      
      save(
        list = ls(envir = environment()),      # all names in this env
        file =  paste0("Poisson_ln_per_t_lambda",round(density, 2), "id", id , ".RData" ),               # where to write
        envir = environment()                  # whose vars to save
      )
      
    }
    
  }  
}

for(id in 1:10){
  print(id)
  do_poisson_lambda(0.2, id)
  do_poisson_lambda(0.4, id)
  do_poisson_lambda(0.8, id)
  
}


