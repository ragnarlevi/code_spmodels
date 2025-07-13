ReadMe
================

# code_spmodels

The code requires the package [tidyverse](https://www.tidyverse.org/) to
be installed: `install.packages("tidyverse")`.

To run the case study on a generate toy data set, please run
`toy_data_analysis.R`.

The experiments in the paper were run by running `poisson_test.R`,
`zip_test.R`, `hurdle_test.R`, `speed_test.R` and `test_direct_vs_em.R`,
but they do take quite some time to run. To generate the plots, run
`visualize.R` which reads the data from the simulation experiments from
the folder `data_sim`.

Example code can be seen below.


## Run the models

``` r
source("utils.R")  # gradients and helper functions
source("poisson.R")  # All Poisson models models
source("zip_mixed.R")  # Zip and mixed zip models
source("hurdle_mixed.R") # hurdle and mixed hurdle models
source("simulate_data.R")  # to gernerate data


#####################
# Test Poisson
#####################
# Poisson is purely L-BFGS-B
sim <- simulate_claims(50, 100, spatial_type = "graph", additive =  FALSE, area = 10, 
                           model_type = "poisson", mixing = "gamma", density = 0.4,  seed = 1)

out_poisson <- Poisson(sim$claims, sim$X, sim$locs, sim$years, sim$agg_claims, sim$A, additive = FALSE, model_type = "learn_graph", 
                       lambda = 0, 
                       exposure = sim$exposure, 
                       max_itr = 500)

out_poisson_mixed <- Poisson_mixed(claims = sim$claims, 
                           X = sim$X, 
                           locs = sim$locs, 
                           years = sim$years, 
                           agg_claims = sim$agg_claims,
                           A = sim$A,   # can be NA if A is learned
                           additive = FALSE, 
                           model_type = "learn_graph", 
                           lambda = 0,
                           exposure = sim$exposure,  
                           mixing_var = "gamma", 
                            Q_tol = 0,
                           nr_em = 60,
                           verbose = 0  # increase number for info
                           )
```

    ## [1] "Breaking because parameters have stopped changing"

``` r
# Louis method error, first few parameters
sqrt(diag(solve(out_poisson_mixed$Hessian-out_poisson_mixed$var_loglik)))[1:5]
```

    ##                   X11        X12        X13            
    ## 0.03127191 0.09569682 0.13410355 0.13721541 0.36369486


``` r
# Beta error
cat("Poisson model beta1 error:", sum(abs(out_poisson$beta1 - sim$beta1)), "\n")
```

    ## Poisson model beta1 error: 0.3760549

``` r
cat("Poisson mixed model beta1 error:", sum(abs(out_poisson_mixed$beta1 - sim$beta1)), "\n")
```

    ## Poisson mixed model beta1 error: 0.3562951



``` r
# Alpha Error
cat("Poisson model a error:", sum(abs(out_poisson$a - sim$a)), "\n")
```

    ## Poisson model a error: 8.392337

``` r
cat("Poisson mixed model a error:", sum(abs(out_poisson_mixed$a - sim$a)), "\n")
```

    ## Poisson mixed model a error: 7.922855



``` r
cat("Poisson mixed model beta phi error:", sum(abs(out_poisson_mixed$beta2 - sim$beta2)), "\n")
```

    ## Poisson mixed model beta phi error: 0.01698877

``` r
#####################
# Test Zip
#####################
source("zip_mixed.R")  # Zip and mixed zip models
zip_sim <- simulate_claims(50, 100, spatial_type = "graph", additive =  TRUE, area = 10, 
                           model_type = "zip", mixing = "ln", density = 0.4,  seed = 1)

# not mixed
out_zip <- zip(claims = zip_sim$claims, X1 = zip_sim$X, locs = zip_sim$locs, 
               years = zip_sim$years,  agg_claims = zip_sim$agg_claims, 
                 A = zip_sim$A, additive = TRUE, 
                model_type = "learn_graph", 
               exposure = zip_sim$exposure, 
                 lambda = 0, max_itr = 500)
# mixed
out_zip_mixed <- zip_mixed (claims = zip_sim$claims, X = zip_sim$X, years = zip_sim$years, locs = zip_sim$locs, agg_claims = zip_sim$agg_claims, 
                            A = zip_sim$A,  # can be NA in this case (is ignoted when A is learned)
                            exposure = zip_sim$exposure, model_type = "learn_graph", additive = TRUE, mixing = "ln",  
                            Emethod = "integration",
                            n_iter = 100, 
                            lambda = 0, 
                            optimizer_beta = "gd", optimizer_psi = "gd",
                            optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd",  # use gradient descent for each parameter
                            sgd = FALSE,  # can allow stochastic gradient but EM not guaranteed to take the correct step
                            batch_size = 100, # only used if sgd is true
                            verbose = 0, 
                            do_optim = FALSE,   # generally slower
                            calc_se = TRUE,
                            beta2_start = 0,   # takes some time
                            control_list = list(a = list(lr = 0.001),  # set learning rate, too high can lead to a decrease in lkelihood
                                                  beta_phi  = list(lr = 0.0001)),
                            max_itr_poisson = 500  # max number of iterations to start Poisson
                            )
```

    ## 
    ## Starting EM updates:
    ## [1] "Breaking because log-likelihood has stopped changing"

``` r
# Louis method error
sqrt(diag(solve(out_zip_mixed$Hessian-out_zip_mixed$var_loglik)))[1:5]
```

    ##                                 X11         X12         X13 
    ## 0.044667292 0.009032701 0.221433961 0.618717206 0.600556352



``` r
cat("Poisson model beta1 error:", sum(abs(out_zip$beta1 - zip_sim$beta1)), "\n")
```

    ## Poisson model beta1 error: 0.7673272

``` r
cat("Poisson mixed model beta1 error:", sum(abs(out_zip_mixed$beta1 - zip_sim$beta1)), "\n")
```

    ## Poisson mixed model beta1 error: 0.6042119



``` r
cat("Poisson model a error:", sum(abs(out_zip$a - zip_sim$a)), "\n")
```

    ## Poisson model a error: 17.95661

``` r
cat("Poisson mixed model a error:", sum(abs(out_zip_mixed$a - zip_sim$a)), "\n")
```

    ## Poisson mixed model a error: 16.76104



``` r
cat("Poisson mixed model beta2 error:", sum(abs(out_zip_mixed$beta2 - zip_sim$beta2)), "\n")
```

    ## Poisson mixed model beta2 error: 0.121224

``` r
#####################
# Test Zip psi
#####################

# We can also fit a model where A is known
zip_sim_psi <- simulate_claims(50, 100, spatial_type = "psi", additive =  TRUE, area = 10, 
                           model_type = "zip", mixing = "ln", density = 0.4,  seed = 1)

# not mixed
out_zip_psi <- zip(claims = zip_sim_psi$claims, X1 = zip_sim_psi$X, locs = zip_sim_psi$locs, 
               years = zip_sim_psi$years,  agg_claims = zip_sim_psi$agg_claims, 
                 A = zip_sim_psi$A, additive = TRUE, 
                model_type = "learn_psi", 
               exposure = zip_sim_psi$exposure, 
                 lambda = 0, max_itr = 500)
# mixed
out_zip_psi_mixed <- zip_mixed (claims = zip_sim_psi$claims, X = zip_sim_psi$X, years = zip_sim_psi$years, 
                                locs = zip_sim_psi$locs, 
                                agg_claims = zip_sim_psi$agg_claims, 
                                A = zip_sim_psi$A,   # has to be given when psi is learned
                                exposure = zip_sim_psi$exposure, model_type = "learn_psi", additive = TRUE, mixing = "ln",  
                                Emethod = "integration",
                                n_iter = 30, 
                                lambda = 0, # not used when psi is learned
                                verbose = 0, # set to 2 for convergence info
                                do_optim = FALSE,   # generally slower to use the L_BFGS within the EM.
                                calc_se = TRUE,
                                beta2_start = 0.1,   
                                control_list = list(psi = list(lr = 0.0001),  # set learning rate, too high can lead to a decrease in likelihood
                                                      beta_phi  = list(lr = 0.00001),
                                                    beta = list(lr = 0.0001),
                                                    pi = list(lr = 0.00001)),
                                max_itr_poisson = 100  # max number of iterations to start Poisson
                            )
```

    ## 
    ## Starting EM updates:

``` r
# Louis method error
sqrt(diag(solve(out_zip_psi_mixed$Hessian-out_zip_psi_mixed$var_loglik)))
```

    ##                              X11        X12        X13                                                                                                               
    ## 0.04225049 0.00898712 0.22837113 0.69710626 0.69823973 0.65033708 1.20946297 0.83040650 0.69251232 0.96627441 0.33314975 0.57568987 0.89359822 0.48720342 0.39350259



``` r
cat("ZIP psi model beta1 error:", sum(abs(out_zip_psi$beta1 - zip_sim_psi$beta1)), "\n")
```

    ## ZIP psi model beta1 error: 1.508419

``` r
cat("ZIP psi mixed model beta1 error:", sum(abs(out_zip_psi_mixed$beta1 - zip_sim_psi$beta1)), "\n")
```

    ## ZIP psi mixed model beta1 error: 1.407495


``` r
cat("ZIP psi model a error:", sum(abs(out_zip_psi$psi - zip_sim_psi$psi)), "\n")
```

    ## ZIP psi model a error: 5.169256

``` r
cat("ZIP psi mixed model a error:", sum(abs(out_zip_psi_mixed$psi - zip_sim_psi$psi)), "\n")
```

    ## ZIP psi mixed model a error: 5.167184



``` r
cat("ZIP psi mixed model beta2 error:", sum(abs(out_zip_psi_mixed$beta2 - zip_sim_psi$beta2)), "\n")
```

    ## ZIP psi mixed model beta2 error: 0.07877959

``` r
#####################
# Test Zip psi
#####################

# Test mixed hurdle
source("hurdle_mixed.R")  # Zip and mixed zip models
hurdle_sim <- simulate_claims(50, 100, spatial_type = "graph", additive =  TRUE, area = 10, 
                           model_type = "hurdle", mixing = "ig", density = 0.4,  seed = 1)

out_hurdle <- hurdle(hurdle_sim$claims, hurdle_sim$X, hurdle_sim$locs, hurdle_sim$years,  hurdle_sim$agg_claims, 
                 hurdle_sim$A, TRUE, "learn_graph", lambda = 0, exposure = hurdle_sim$exposure, max_itr = 500)

out_hurdle_mixed <- hurdle_mixed (claims = hurdle_sim$claims, X = hurdle_sim$X, years = hurdle_sim$years, locs = hurdle_sim$locs, agg_claims = hurdle_sim$agg_claims, 
                            A = hurdle_sim$A, exposure = hurdle_sim$exposure, model_type = "learn_graph", additive = TRUE, mixing = "ig",  
                            Emethod = "integration",
                            n_iter = 60, 
                            lambda = 0, 
                            optimizer_beta = "gd", optimizer_psi = "gd",
                            optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd",  # use gradient descent for each parameter
                            sgd = FALSE,  # can allow stochastic gradient but EM not guaranteed to take the correct step
                            batch_size = 100, 
                            verbose = 0, 
                            do_optim = FALSE,   # generally slower
                            calc_se = TRUE,
                            beta2_start = 0.2,   # takes some time
                            control_list = list(psi = list(a = 0.000001),  # set learning rate, too high can lead to a decrease in likelihood
                                                  beta_phi  = list(lr = 0.0001),
                                                beta = list(lr = 0.1)),
                            max_itr_poisson = 100  # number of iterations to start Poisson
                            )
```



``` r
cat("Hurdle model beta1 error:", sum(abs(out_hurdle$beta1 - hurdle_sim$beta1)), "\n")
```

    ## Hurdle model beta1 error: 1.297332

``` r
cat("Hurdle mixed model beta1 error:", sum(abs(out_hurdle_mixed$beta1 - hurdle_sim$beta1)), "\n")
```

    ## Hurdle mixed model beta1 error: 0.8678489

``` r
# Alpha error
cat("\n# Alpha Error\n")
```

    ## 
    ## # Alpha Error

``` r
cat("Hurdle model a error:", sum(abs(out_hurdle$a - hurdle_sim$a)), "\n")
```

    ## Hurdle model a error: 20.84158

``` r
cat("Hurdle mixed model a error:", sum(abs(out_hurdle_mixed$a - hurdle_sim$a)), "\n")
```

    ## Hurdle mixed model a error: 8.036557

``` r
# Beta_phi error
cat("\n# Beta_phi Error\n")
```

    ## 
    ## # Beta_phi Error

``` r
cat("Hurdle mixed model beta2 error:", sum(abs(out_hurdle_mixed$beta2 - hurdle_sim$beta2)), "\n")
```

    ## Hurdle mixed model beta2 error: 0.04425221


