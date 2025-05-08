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



do_zip_t <- function(density, id){
  
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
                    A = A, additive, model_type = "learn_graph", lambda = 0, exposure = exposure, max_itr = 300, a_known = TRUE)
    
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
                         batch_size = 100, param_tol = 1e-5, Q_tol = 1e-9, verbose = 2, do_optim = FALSE, a_known = FALSE)
    
    
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
                               batch_size = 100, param_tol = 1e-5, Q_tol = 1e-9, verbose = 2, do_optim = FALSE, a_known = TRUE)



    beta_est_known_ig[[as.character(t)]] <- out_ig_known$beta1
    beta2_est_known_ig[[as.character(t)]] <- out_ig_known$beta2

    
    save(
      list = ls(envir = environment()),      # all names in this env
      file =  paste0("zip_ln_per_t",round(density, 2), "id", id , ".RData" ),               # where to write
      envir = environment()                  # whose vars to save
    )
    
  }
  
}


for(id in 1:10){
  print(id)
  do_zip_t(0.4, id)
  do_zip_t(0.8, id)
}



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
  out_ig <- zip_mixed (sim$claims, sim$X, sim$years, sim$locs, sim$agg_claims, NA, sim$exposure, 
                       model_type = "learn_graph", additive = additive, mixing = "ig",  Emethod = "integration",
                      n_iter = 80, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                      optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                      batch_size = 100, param_tol = 1e-5, Q_tol = 1e-5, verbose = 2, do_optim = FALSE)
                                                                                                                                                                                          
  time_ig[[as.character(t)]] <-  Sys.time() - a1
  
  A_est_ig[[as.character(t)]] <- out_ig$a
  A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
  beta_est_ig[[as.character(t)]] <- out_ig$beta1
  H_ig_p[[as.character(t)]] <- out_ig$Hessian
  var_ig_p[[as.character(t)]]<- out_ig$var_loglik

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
  var_ig_p_lfbgs[[as.character(t)]]<- out_ig_lfbgs$var_loglik
  
  
  
  
  save.image("zip_ig_per_t2.RData")
  
}


library(tidyverse)
load("zip_ig_per_t2.RData")
a_none_error <- lapply(A_est_p, function(x, y){sum(abs(x-y))}, y = sim$A[upper.tri(sim$A, T)])
a_ig_error <- lapply(A_est_ig, function(x, y){sum(abs(x-y))}, y =  sim$A[upper.tri(sim$A, T)])
a_ig_bfgs_error <- lapply(A_est_ig_lfbgs, function(x, y){sum(abs(x-y))}, y =  sim$A[upper.tri(sim$A, T)])


b_none_error <- lapply(beta_est_p, function(x, y){sum(abs(x-y))}, y = sim$beta1)
b_ig_error <- lapply(beta_est_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)
b_ig_bfgs_error <- lapply(beta_est_ig_lfbgs, function(x, y){sum(abs(x-y))}, y = sim$beta1)

se_error_ig <- do.call(rbind,
                       mapply(
                         function(x, idx) {
                           # idx is a single integer (1, 2, 3, …)
                           mat <- x - 0.5 * diag(nrow(x)) * (idx == 1)
                           diag( solve(mat) )
                         },
                         H_ig_p,
                         seq_along(H_p),
                         SIMPLIFY = FALSE
                       )
)

se_error_p <- do.call(rbind,
                      mapply(
                        function(x, idx) {
                          # idx is a single integer (1, 2, 3, …)
                          mat <- x + 0.5 * diag(nrow(x)) * (idx == 1)
                          diag( solve(mat) )
                        },
                        H_p,
                        seq_along(H_p),
                        SIMPLIFY = FALSE
                      )
)

ggplot() + 
  geom_line(aes( x = ts, y = se_error_p[,1], color = "Poisson"))+
  geom_line(aes( x = ts, y = se_error_ig[,1], color = "LN"))




est <- sapply(beta_est_ig, function(x) x[1])
se_err <- se_error_p[,1]
ggplot() + geom_line(aes( x = log( ts[1:7]), y = est)) +
  geom_ribbon(aes(x = log( ts[1:7]), ymin = est - 2*se_err,  ymax = est + 2*se_err ), alpha = 0.5)


est <- sapply(A_est_p, function(x) x[2])
se_err <- se_error_p[,3+2]
ggplot() + geom_line(aes( x = log(ts), y = est)) +
  geom_ribbon(aes(x = log(ts), ymin = est - 2*se_err,  ymax = est + 2*se_err ), alpha = 0.5)


# 
# 
# 

ggplot(data.frame(t = ts,
                  poissson_error = unlist(a_none_error),
                  ig_error = unlist(a_ig_error),
                  bfgs_error = unlist(a_ig_bfgs_error)
)) +
  geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = t, y = ig_error, color = "LN"))+
  geom_line(aes(x = t, y = bfgs_error, color = "LN BFGS"))
# 
# 
# 
ggplot(data.frame(t = ts,
                  poissson_error = unlist(b_none_error),
                  ig_error = unlist(b_ig_error),
                  ig_known_error = unlist(b_ig_known_error)
)) +
  geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = t, y = ig_error, color = "LN"))+
  geom_line(aes(x = t, y = ig_known_error, color = "LN, A known"))




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



