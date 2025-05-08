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
  sim <- simulate_claims(50, 5000, spatial_type = "graph",additive =  additive, area = nr_regions, 
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

# do_poisson_t(0.4, 1)
# do_poisson_t(0.8, 1)

do_poisson_t(0.4, 2)
do_poisson_t(0.8, 2)

do_poisson_t(0.4, 3)
do_poisson_t(0.8, 3)

do_poisson_t(0.4, 4)
do_poisson_t(0.8, 4)

do_poisson_t(0.4, 5)
do_poisson_t(0.8, 5)

do_poisson_t(0.4, 6)
do_poisson_t(0.8, 6)

do_poisson_t(0.4, 7)
do_poisson_t(0.8, 7)

do_poisson_t(0.4, 8)
do_poisson_t(0.8, 8)

do_poisson_t(0.4, 9)
do_poisson_t(0.8, 9)

do_poisson_t(0.4, 10)
do_poisson_t(0.8, 10)


##########
######### Test behaviour with lambda ############
##########


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
  
  
  lambdas <- c(0, 0.5, 1, 2.5, 5, 7, 10, 25, 50, 75, 100)
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


do_poisson_lambda(0.4, 1)
do_poisson_lambda(0.8, 1)

do_poisson_lambda(0.4, 2)
do_poisson_lambda(0.8, 2)

do_poisson_lambda(0.4, 3)
do_poisson_lambda(0.8, 3)

do_poisson_lambda(0.4, 4)
do_poisson_lambda(0.8, 4)

do_poisson_lambda(0.4, 5)
do_poisson_lambda(0.8, 5)

do_poisson_lambda(0.4, 6)
do_poisson_lambda(0.8, 6)

do_poisson_lambda(0.4, 7)
do_poisson_lambda(0.8, 7)

do_poisson_lambda(0.4, 8)
do_poisson_lambda(0.8, 8)

do_poisson_lambda(0.4, 9)
do_poisson_lambda(0.8, 9)

do_poisson_lambda(0.4, 10)
do_poisson_lambda(0.8, 10)





tmp <- new.env()
loaded_objs <- load("Poisson_ln_per_t0.4id1.RData", envir = tmp)
my_list <- mget(loaded_objs, envir = tmp)







plot_per_t <- function(ids, dens){
  betas_p <- list()
  betas_mp <- list()
  betas_mp_known <- list()
  
  betas2_mp <- list()
  betas2_mp_known <- list()
  
  a_p <- list()
  a_mp <- list()
  a_mp_known <- list()
  for(id in ids){
    
    tmp <- new.env()
    loaded_objs <- load(paste0("Poisson_ln_per_t", dens, "id", id, ".RData") , envir = tmp)
    my_list <- mget(loaded_objs, envir = tmp)
    
    
    betas_p[[as.character(id)]] <- bind_rows(my_list$beta_est_p)
    betas_p[[as.character(id)]]$id <- id
    betas_p[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_p))
    
    betas_mp[[as.character(id)]] <- bind_rows(my_list$beta_est_ig)
    betas_mp[[as.character(id)]]$id <- id
    betas_mp[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
    
    betas_mp_known[[as.character(id)]] <- bind_rows(my_list$beta_est_known_ig)
    betas_mp_known[[as.character(id)]]$id <- id
    betas_mp_known[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_known_ig))
    
    a_p[[as.character(id)]] <- as.data.frame(t(bind_rows(my_list$A_est_p)))
    a_p[[as.character(id)]]$id <- id
    a_p[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_p))
    
    a_mp[[as.character(id)]] <- as.data.frame(t(bind_rows(my_list$A_est_ig)))
    a_mp[[as.character(id)]]$id <- id
    a_mp[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
    
    betas2_mp[[as.character(id)]] <- as.data.frame(t(data.frame(my_list$beta2_est_ig)))
    betas2_mp[[as.character(id)]]$id  <- id
    betas2_mp[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
    
    betas2_mp_known[[as.character(id)]]<- as.data.frame(t(data.frame(my_list$beta2_est_known_ig)))
    betas2_mp_known[[as.character(id)]]$id  <- id
    betas2_mp_known[[as.character(id)]]$t <- as.numeric(names(my_list$beta_est_ig))
    
    
    a_true <- my_list$sim$a
    beta_true <- my_list$sim$beta1
    beta2_true <- my_list$sim$beta2
    
  }
  
  betas_p <- bind_rows(betas_p)
  betas_mp <- bind_rows(betas_mp)
  betas_mp_known <- bind_rows(betas_mp_known)
  a_p <- bind_rows(a_p)
  a_mp <- bind_rows(a_mp)
  
  betas2_mp <- bind_rows(betas2_mp)
  betas2_mp_known <- bind_rows(betas2_mp_known)
  
  # error 
  
  beta_true_matrix <- matrix(beta_true,
                             nrow = nrow(betas_p),
                             ncol = length(beta_true),
                             byrow = TRUE)
  
  betas_p$error <- rowSums(abs(betas_p[, c("X11", "X12", "X13")] - beta_true_matrix))
  betas_mp$error <- rowSums(abs(betas_mp[, c("X11", "X12", "X13")] - beta_true_matrix))
  betas_mp_known$error <- rowSums(abs(betas_mp_known[, c("X11", "X12", "X13")] - beta_true_matrix))
  
  
  a_true_matrix <- matrix(a_true,
                          nrow = nrow(a_p),
                          ncol = length(a_true),
                          byrow = TRUE)
  a_p$error <- rowSums(abs(a_p[ ,1:(10*11/2)] - a_true_matrix))
  a_mp$error <- rowSums(abs(a_mp[ ,1:(10*11/2)] - a_true_matrix))
  
  
  ret <- list()
  

  df1 <- betas_p %>% mutate(model = "Poisson")
  df2 <- betas_mp %>% mutate(model = "Mixed Poisson")
  df3 <- betas_mp_known %>% mutate(model = "Mixed Poisson A known")
  
  # bind into one
  betas_all <- bind_rows(df1, df2, df3)
  
  # dodged box‐plot
  ret$beta1_plot <- ggplot(betas_all, aes(x = factor(t), y = X11, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = 1, group = 1),
      inherit.aes = FALSE,
      color       = "black",
      size        = 1,
      alpha = 0.5
    ) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) + 
    labs(
      x     = "Time (t)",
      y     = expression(beta[1] ~ " estimate"),
      fill  = "Model spec",
      title = expression(beta[1] ~ " over Time, by Model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  # dodged box‐plot
  ret$beta_error_plot <- ggplot(betas_all, aes(x = factor(t), y = error, fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) + 
    labs(
      x     = "Time (t)",
      y     = expression(l[1] ~ " error"),
      fill  = "Model spec",
      title = expression(l[1] ~ " error over Time, by Model")
    ) +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  
  df1 <- a_p %>% mutate(model = "Poisson")
  df2 <- a_mp %>% mutate(model = "Mixed Poisson")
  
  # bind into one
  as_all <- bind_rows(df1, df2)
  
  # dodged box‐plot
  ret$edge_plot <- ggplot(as_all, aes(x = factor(t), y = V3, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = a_true[3], group = 1),
      inherit.aes = FALSE,
      color       = "black",
      size        = 1,
      alpha = 0.5
    ) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) + 
    labs(
      x     = "Time (t)",
      y     = latex2exp::TeX("$a_{1,1}$ estimate"),
      fill  = "Model spec",
      title = latex2exp::TeX("$a_{1,1}$ over time, by model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  
  # dodged box‐plot
  ret$edge_error_plot <- ggplot(as_all, aes(x = factor(t), y = log(error), fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) + 
    labs(
      x     = "Time (t)",
      y     = latex2exp::TeX("$a_{1,1}$ estimate"),
      fill  = "Model spec",
      title = latex2exp::TeX("$a_{1,1}$ over time, by model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  
  
  
  
  
  
  
  
  df1 <- betas2_mp %>% mutate(model = "Mixed Poisson")
  df2 <- betas2_mp_known %>% mutate(model = "Mixed Poisson known A")
  
  # bind into one
  beta2_all <- bind_rows(df1, df2)
  
  # dodged box‐plot
  ret$beta2_plot <- ggplot(beta2_all, aes(x = factor(t), y = V1, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = 0, group = 1),
      inherit.aes = FALSE,
      color       = "black",
      size        = 1,
      alpha = 0.5
    ) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) + 
    labs(
      x     = "Time (t)",
      y     = latex2exp::TeX("$a_{1,1}$ estimate"),
      fill  = "Model spec",
      title = latex2exp::TeX("$a_{1,1}$ over time, by model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  
return(ret)
}

ids <- 2:10
dens <- 0.4

out_plot_per_t <- plot_per_t(ids, dens)


library(tidyverse)
# 
# load("Poisson_ln_per_t0.4id1.RData", , envir = environment())
# 
# 
a_none_error <- lapply(A_est_p, function(x, y){sum(abs(x-y))}, y = sim$A[upper.tri(sim$A, T)])
a_ig_error <- lapply(A_est_ig, function(x, y){sum(abs(x-y))}, y =  sim$A[upper.tri(sim$A, T)])


b_none_error <- lapply(beta_est_p, function(x, y){sum(abs(x-y))}, y = sim$beta1)
b_ig_error <- lapply(beta_est_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)
b_ig_known_error <- lapply(beta_est_known_ig, function(x, y){sum(abs(x-y))}, y = sim$beta1)





se_error_ig <- do.call(rbind,
                       mapply(
                         function(x, idx, y) {
                           # idx is a single integer (1, 2, 3, …)
                           mat <- x-y + 1 * diag(nrow(x)) * (idx == 1)
                           diag( solve(mat) )
                         },
                         H_ig_p,
                         seq_along(H_p),
                         var_ig_p,
                         SIMPLIFY = FALSE
                       )
)

se_error_p <- do.call(rbind,
                      mapply(
                        function(x, idx) {
                          # idx is a single integer (1, 2, 3, …)
                          mat <- x + 1 * diag(nrow(x)) * (idx == 1)
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
se_err <- se_error_ig[,2]
ggplot() + geom_line(aes( x = log(ts), y = est)) +
  geom_ribbon(aes(x = log(ts), ymin = est - 2*se_err,  ymax = est + 2*se_err ), alpha = 0.5)


est <- sapply(A_est_p, function(x) x[1])
se_err <- se_error_p[,1 + 3+1]
ggplot() + geom_line(aes( x = log(ts), y = est)) +
  geom_ribbon(aes(x = log(ts), ymin = est - 2*se_err,  ymax = est + 2*se_err ), alpha = 0.5)


# 
# 
# 



ggplot(data.frame(t = log(ts),
                   poissson_error = unlist(a_none_error),
                   ig_error = unlist(a_ig_error)
                     )) +
  geom_line(aes(x = t, y = poissson_error, color = "Poisson"))+
  geom_line(aes(x = t, y = ig_error, color = "LN"))
# 
# 
# 
ggplot(data.frame(t = log(ts),
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

