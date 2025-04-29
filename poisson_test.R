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



beta_est_known_ig <- list()

ts <- c(10, 20, 50, 100, 300, 500, 1000, 2000)
additive <- FALSE
for(t in ts){
  
  print(t)
  sim <- simulate_claims(100, t, spatial_type = "graph",additive =  additive, area = 10, model_type = "poisson", mixing = "ln")
  
  # Poisson
  a1 <- Sys.time()
  out_none <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, additive = additive, "learn_graph", 
                      lambda = 0, exposure = sim$exposure, max_itr = 100)
  time_p[[as.character(t)]] <-  Sys.time() - a1
  A_p_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
  
  A_est_p[[as.character(t)]] <- out_none$a
  beta_est_p[[as.character(t)]] <- out_none$beta1
  H_p[[as.character(t)]] <- out_none$optim_obj$hessian

  # Mixed Poisson
  a1 <- Sys.time()
  out_ig <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, A = NA, additive = additive, model_type = "learn_graph", 
                          lambda = 0, 
                          exposure = sim$exposure, max_itr = 10, mixing_var = "ln", nr_em = 50, verbose = 2)
  time_ig[[as.character(t)]] <-  Sys.time() - a1
  
  A_est_ig[[as.character(t)]] <- out_ig$a
  A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
  beta_est_ig[[as.character(t)]] <- out_ig$beta1
  H_ig_p[[as.character(t)]] <- out_ig$Hessian
  var_ig_p[[as.character(t)]] <- out_ig$var_der_log_lik
  
  
  # Mixed Poisson known A
  out_ig_known <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, A = sim$A, additive = additive, model_type = "learn_graph", 
                                lambda = 0, 
                                exposure = sim$exposure, max_itr = 10, mixing_var = "ln", nr_em = 50, verbose = 2, a_known = TRUE)
  
  
  
  beta_est_known_ig[[as.character(t)]] <- out_ig_known$beta1
  

  
  save.image("Poisson_ln_per_t.RData")

}





library(tidyverse)
# 
# load("Poisson_ln_per_t.RData")
# 
# 
a_none_error <- lapply(A_est_p, function(x, y){sum(abs(x-y))}, y = sim$A[upper.tri(sim$A, T)])
a_ig_error <- lapply(A_est_ig, function(x, y){sum(abs(x-y))}, y =  sim$A[upper.tri(sim$A, T)])


b_none_error <- lapply(beta_est_p, function(x, y){sum(abs(x-y))}, y = sim$beta1)
b_ig_error <- lapply(beta_est_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)
b_ig_known_error <- lapply(beta_est_known_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)





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
ggplot() + geom_line(aes( x = log(ts), y = est)) +
  geom_ribbon(aes(x = log(ts), ymin = est - 2*se_err,  ymax = est + 2*se_err ), alpha = 0.5)


est <- sapply(A_est_p, function(x) x[2])
se_err <- se_error_p[,3+2]
ggplot() + geom_line(aes( x = log(ts), y = est)) +
  geom_ribbon(aes(x = log(ts), ymin = est - 2*se_err,  ymax = est + 2*se_err ), alpha = 0.5)


# 
# 
# 

ggplot(data.frame(t = ts,
                   poissson_error = unlist(a_none_error),
                   ig_error = unlist(a_ig_error)
                     )) +
  geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = t, y = ig_error, color = "LN"))
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
# 
# 
# 
# 
# 
# A_ig_info[["10"]]$relative_error
# 
# 
# 



##########
######### Test behaviour with x ############
##########

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



beta_est_known_ig <- list()

xs <- c(10, 20, 50, 100, 300, 500, 1000)
additive <- FALSE
for(t in xs){
  
  print(t)
  sim <- simulate_claims(x, 100, spatial_type = "graph",additive =  additive, area = 10, model_type = "poisson", mixing = "ln")
  
  # Poisson
  a1 <- Sys.time()
  out_none <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, additive = additive, "learn_graph", 
                      lambda = 0, exposure = sim$exposure, max_itr = 100)
  time_p[[as.character(t)]] <-  Sys.time() - a1
  A_p_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
  
  A_est_p[[as.character(t)]] <- out_none$a
  beta_est_p[[as.character(t)]] <- out_none$beta1
  H_p[[as.character(t)]] <- out_none$optim_obj$hessian
  
  # Mixed Poisson
  a1 <- Sys.time()
  out_ig <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, A = NA, additive = additive, model_type = "learn_graph", 
                          lambda = 0, 
                          exposure = sim$exposure, max_itr = 10, mixing_var = "ln", nr_em = 50, verbose = 2)
  time_ig[[as.character(t)]] <-  Sys.time() - a1
  
  A_est_ig[[as.character(t)]] <- out_ig$a
  A_ig_info[[as.character(t)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
  beta_est_ig[[as.character(t)]] <- out_ig$beta1
  H_ig_p[[as.character(t)]] <- out_ig$Hessian
  var_ig_p[[as.character(t)]] <- out_ig$var_der_log_lik
  
  
  # Mixed Poisson known A
  out_ig_known <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, A = sim$A, additive = additive, model_type = "learn_graph", 
                                lambda = 0, 
                                exposure = sim$exposure, max_itr = 10, mixing_var = "ln", nr_em = 50, verbose = 2, a_known = TRUE)
  
  
  
  beta_est_known_ig[[as.character(t)]] <- out_ig_known$beta1
  
  
  
  save.image("Poisson_ln_per_x.RData")
  
}


##########
############### Test behaviour with lambda ###########
##########

source("poisson.R")
source("simulate_data.R")


A_est_p <- list()
A_p_info <- list()
beta_est_p <- list()

A_est_ig <- list()
A_ig_info <- list()
beta_est_ig <- list()


lambdas <- c(0, 1, 7, 5, 10, 25, 50, 75, 100, 125, 500, 1000)
additive <- FALSE
for(lambda in lambdas){
  for(t in c(50, 100, 200)){
    
    sim <- simulate_claims(100, t, spatial_type = "graph",additive =  additive, area = 10, model_type = "poisson", mixing = "ln")
    
    # Poisson
    out_none <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, NA, additive = additive, "learn_graph", 
                        lambda = lambda, exposure = sim$exposure, max_itr = 100)
    A_p_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_none$a, 10))
    
    A_est_p[[paste0(t, "-", lambda)]] <- out_none$a
    beta_est_p[[paste0(t, "-", lambda)]] <- out_none$beta1
    
    # Mixed Poisson
    out_ig <- Poisson_mixed(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, A = NA, additive = additive, model_type = "learn_graph", 
                            lambda = lambda, 
                            exposure = sim$exposure, max_itr = 10, mixing_var = "ln", nr_em = 50, verbose = 2)
    
    
    A_est_ig[[paste0(t, "-", lambda)]] <- out_ig$a
    A_ig_info[[paste0(t, "-", lambda)]] <- evaluate_matrix_recovery(sim$A, get_W_from_array(out_ig$a, 10))
    beta_est_ig[[paste0(t, "-", lambda)]] <- out_ig$beta1
    
    

    save.image("Poisson_ln_per_t_lambda.RData")
    
  }
  
  
  
}





load("Poisson_ln_per_t_lambda.RData")

A_p_info_list <- list()
A_est_p_list <- list()
beta_est_p_list <- list()

A_ig_info_list <- list()
A_est_ig_list <- list()
beta_est_ig_list <- list()
for(t in c(50, 100, 200)){
  A_p_info_list[[as.character(t)]] <- list()
  A_est_p_list[[as.character(t)]]<- list()
  A_p_info_list[[as.character(t)]]<- list()
  
  A_ig_info_list[[as.character(t)]]<- list()
  A_est_ig_list[[as.character(t)]]<- list()
  beta_est_ig_list[[as.character(t)]]<- list()
  
  for(lambda in lambdas){

    A_p_info_list[[as.character(t)]][[as.character(lambda)]] <- A_p_info[[paste0(t, "-", lambda)]]
    A_est_p_list[[as.character(t)]][[as.character(lambda)]] <- A_est_p[[paste0(t, "-", lambda)]]
    beta_est_p_list[[as.character(t)]][[as.character(lambda)]] <- beta_est_p[[paste0(t, "-", lambda)]]
    A_ig_info_list[[as.character(t)]][[as.character(lambda)]] <- A_ig_info[[paste0(t, "-", lambda)]]
    A_est_ig_list[[as.character(t)]][[as.character(lambda)]] <- A_est_ig[[paste0(t, "-", lambda)]]
    beta_est_ig_list[[as.character(t)]][[as.character(lambda)]] <- beta_est_ig[[paste0(t, "-", lambda)]]

    
  }
  
  
}





ggplot(data.frame(lambda = lambdas,
                  poissson_error = unlist(lapply(A_p_info_list[["50"]], function(x){x$f1_score})),
                  ig_error = unlist(lapply(A_ig_info_list[["50"]], function(x){x$f1_score}))
)) +
  geom_line(aes(x = lambda, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = lambda, y = ig_error, color = "LN"))



ggplot(data.frame(lambda = lambdas,
                  poissson_error = unlist(lapply(A_p_info_list[["100"]], function(x){x$precision})),
                  ig_error = unlist(lapply(A_ig_info_list[["100"]], function(x){x$precision}))
)) +
  geom_line(aes(x = lambda, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = lambda, y = ig_error, color = "LN"))



ggplot(data.frame(lambda = lambdas,
                  poissson_error = unlist(lapply(A_p_info_list[["100"]], function(x){x$recall})),
                  ig_error = unlist(lapply(A_ig_info_list[["100"]], function(x){x$recall}))
)) +
  geom_line(aes(x = lambda, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = lambda, y = ig_error, color = "LN"))






A_est <- get_W_from_array(A_est_ig_list[["50"]][["10"]], 10)
A_est[abs(A_est) <= 1e-3 ]<- 0
A_true <- sim$A

graph_plot_out <- plot_graphs(get_W_from_array(A_est_ig_list[["50"]][["10"]], 10), sim$A)


graph_plot_out$side_by_side

graph_plot_out$true_false_graph

