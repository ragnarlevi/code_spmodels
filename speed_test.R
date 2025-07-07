source("zip_mixed.R")
source("simulate_data.R")
#######################################
# The goal of this Document is to analyze the time complexity of each iteration in the EM algorithm


do_zip_speed <- function(density, id, calc_se){
  

  time_ig <- list()
  nr_data_points <- list()
  time_to_init <- list() 
  iteration_times <- list()

  additive <- TRUE
  nr_regions <- c(10, 20, 50, 70, 100, 200)
  
  for(t in c(50, 100, 500, 1000)){
    for(nr_region in nr_regions){
      

      print(paste(t, "-", nr_region))
      
      # Order regions so we do not have an infinite loop in the data_generation process
      sim <- simulate_claims(200, t, spatial_type = "graph", additive =  additive, area = nr_region, 
                             model_type = "zip", mixing = "gamma", density = density,  seed = id, , order_regions = TRUE)
      print("sim over")
      
      # subset claims
      claims <- sim$claims
      locs <- sim$locs
      agg_claims <- sim$agg_claims
      X <- sim$X
      years <- sim$years
      exposure <- sim$exposure
      
      
      # Mixed Poisson
      a1 <- Sys.time()
      
      out_ig <- zip_mixed (claims, X, years, locs, agg_claims, A, exposure, 
                           model_type = "learn_graph", additive = additive, mixing = "gamma",  Emethod = "integration",
                           n_iter = 5, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                           optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                           batch_size = 100, param_tol = 0, Q_tol =0, verbose = 1, do_optim = FALSE, a_known = FALSE, calc_se = FALSE, 
                           max_itr_poisson = 10)
      
      
      time_ig[[paste0(t,"-" ,nr_region)]] <-  Sys.time() - a1
      nr_data_points[[paste0(t,"-" ,nr_region)]] <-  nrow(X)
      time_to_init[[paste0(t,"-" ,nr_region)]] <- out_ig$time_to_init
      iteration_times[[paste0(t,"-" ,nr_region)]] <- out_ig$iteration_times
      
      save(
        list = ls(envir = environment()),      # all names in this env
        file =  paste0("sim_data/2zip_speed",round(density, 2), "id", id , ".RData" ),               # where to write
        envir = environment()                  # whose vars to save
      )
      
    } 
  }
  
  
  
}


do_zip_speed(0.4, 1, FALSE)
