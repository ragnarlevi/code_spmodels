# We need to test
# 1. Do we uncover true A.
# Assume A known then assume unkown, does beta change?
# 2. Look at how fit behaves with n and t.
# 3. Computational complexity with number of parameter and data points.
# 4. How do standard errors converge. Plot standard errors of estimates
# 5. Conway Maxwell Poisson (COMP) distributions


# Note as pi is bigger from the initial zip estimated, beta_phi will be higher
#


source("zip_mixed.R")
source("simulate_data.R")
#######################################
# The goal of this Document is to perform the tests above for the
# ordinary poisson models.




A_est_p <- list()
A_p_info <- list()
beta_est_p <- list()
H_p <- list()
time_p <- list()

A_est_ig <- list()
A_ig_info <- list()
beta_est_ig <- list()
H_ig_p <- list()
var_ig_p <- list()
time_ig <- list()


A_est_ig_lfbgs <- list()
A_ig_info_lfbgs <- list()
beta_est_ig_lfbgs <- list()
H_ig_p_lfbgs <- list()
var_ig_p_lfbgs <- list()
time_ig_lfbgs <- list()


ts <- c(10, 20, 50, 100, 300, 500, 1000, 2000)
additive <- FALSE
for(t in ts){
  print(t)
  sim <- simulate_claims(100, t, spatial_type = "graph",additive =  additive, area = 10, model_type = "zip", mixing = "ig")
  
  # Poisson
  a1 <- Sys.time()
  

  out_none <- zip(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, additive = additive, "learn_graph", 
                      lambda = 0, exposure = sim$exposure, max_itr = 100)
  time_p[[as.character(t)]] <-  Sys.time() - a1
  A_p_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
  
  A_est_p[[as.character(t)]] <- out_none$a
  beta_est_p[[as.character(t)]] <- out_none$beta1
  H_p[[as.character(t)]] <- out_none$optim_obj$hessian
  
  
  
  # Mixed Poisson
  a1 <- Sys.time()
  source("zip_mixed.R")

  out_ig <- zip_mixed (sim$claims, sim$X, sim$years, sim$locs, sim$agg_claims, NA, sim$exposure, 
                       model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
                      n_iter = 50, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                      optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                      batch_size = 100, param_tol = 1e-5, Q_tol = 1e-5, verbose = 2, do_optim = FALSE)
                                                                                                                                                                                          
  time_ig[[as.character(t)]] <-  Sys.time() - a1
  
  A_est_ig[[as.character(t)]] <- out_ig$a
  A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
  beta_est_ig[[as.character(t)]] <- out_ig$beta1
  H_ig_p[[as.character(t)]] <- out_ig$Hessian
  
  
 #  Mixed Poisson
  a1 <- Sys.time()
  out_ig_lfbgs <- zip_mixed (sim$claims, sim$X, sim$years, sim$locs, sim$agg_claims, NA, sim$exposure,
                       model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
                       n_iter = 50, lambda = 0, param_tol = 1e-5, Q_tol = 1e-5, verbose = 2,do_optim = TRUE)

  time_ig_lfbgs[[as.character(t)]] <-  Sys.time() - a1

  A_est_ig_lfbgs[[as.character(t)]] <- out_ig_lfbgs$a
  A_ig_info_lfbgs[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig_lfbgs$a, 10))
  beta_est_ig_lfbgs[[as.character(t)]] <- out_ig_lfbgs$beta1
  H_ig_p_lfbgs[[as.character(t)]] <- out_ig_lfbgs$Hessian
  
  
  
  
  save.image("zip_ig_per_t2.RData")
  
}



# load("Poisson_ln_per_t.RData")
# 
# 
# a_none_error <- lapply(A_est_p, function(x, y){sum(abs(x-y))}, y = sim$A[upper.tri(sim$A, T)])
# a_ig_error <- lapply(A_est_ig, function(x, y){sum(abs(x-y))}, y =  sim$A[upper.tri(sim$A, T)])
# 
# 
# b_none_error <- lapply(beta_est_p, function(x, y){sum(abs(x-y))}, y = sim$beta1)
# b_ig_error <- lapply(beta_est_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)
# 
# 
# 
# library(tidyverse)
# ggplot(data.frame(t = ts,
#                   poissson_error = unlist(a_none_error),
#                   ig_error = unlist(a_ig_error)
# )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "LN"))
# 
# 
# 
# ggplot(data.frame(t = ts,
#                   poissson_error = unlist(b_none_error),
#                   ig_error = unlist(b_ig_error)
# )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "LN"))
# 
# 
# A_ig_info
# 
# 
# 
# ggplot(data.frame(t = ts,
#                   poissson_error = unlist(lapply(A_p_info, function(x){x$f1_score})),
#                   ig_error = unlist(lapply(A_ig_info, function(x){x$f1_score}))
# )) +
#   geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
#   geom_line(aes(x = t, y = ig_error, color = "LN"))
# 







##########
############### Test behaviour with lambda ###########
##########

source("zip_mixed.R")
source("simulate_data.R")


A_est_p <- list()
A_p_info <- list()
beta_est_p <- list()

A_est_ig <- list()
A_ig_info <- list()
beta_est_ig <- list()


lambdas <- c(0, 1, 2,  5, 7, 10, 15, 25, 50,)
additive <- FALSE
for(lambda in lambdas){
  for(t in c(50, 100, 200)){
    
    sim <- simulate_claims(100, t, spatial_type = "graph",additive =  additive, area = 10, model_type = "zip", mixing = "ig")
    
    # zip
    out_none <- zip(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, additive = additive, "learn_graph", 
                    lambda = lambda, exposure = sim$exposure, max_itr = 100)
    A_p_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
    
    A_est_p[[paste0(t, "-", lambda)]] <- out_none$a
    beta_est_p[[paste0(t, "-", lambda)]] <- out_none$beta1
    
    # Mixed zip
    out_ig <- zip_mixed (sim$claims, sim$X, sim$years, sim$locs, sim$agg_claims, NA, sim$exposure, 
                         model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
                         n_iter = 50, lambda = lambda, optimizer_beta = "gd", optimizer_psi = "gd",
                         optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                         batch_size = 100, param_tol = 1e-5, Q_tol = 1e-5, verbose = 2, do_optim = FALSE)
    
    
    A_est_ig[[paste0(t, "-", lambda)]] <- out_ig$a
    A_ig_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
    beta_est_ig[[paste0(t, "-", lambda)]] <- out_ig$beta1
    
    
    
    save.image("zip_ig_per_t_lambda.RData")
    
  }
  
  
  
}



