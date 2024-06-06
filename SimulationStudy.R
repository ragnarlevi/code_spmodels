
library(tidyverse)


# remove everything from environment
rm(list = ls())
get_errors <- function(beta1_list, a_list, beta_true_list, a_true_list, betas, ns, is_matrix = TRUE){
  
  error_a <- list()
  error_beta1 <- list()
  beta1_ord <- bind_rows(beta1_list)
  a_ord <- bind_rows(a_list)
  
  for(ni in ns){
    print(ni)
    for(betai in betas){
      
      tmp_beta <- bind_rows(beta1_list) %>% filter(n == ni, beta==betai)
      tmp_a <- bind_rows(a_list) %>% filter(n == ni, beta==betai)
      
      beta1_tmp <- beta_true_list[[paste(betai, ni)]]
      error_beta1[[paste(ni,betai)]] <- data.frame(n = ni, beta = betai,
                                                   Poisson = sum(abs(tmp_beta$Poisson - beta1_tmp )),
                                                   PG = sum(abs(tmp_beta$PG - beta1_tmp )),
                                                   PIG = sum(abs(tmp_beta$PIG - beta1_tmp )),
                                                   PLN = sum(abs(tmp_beta$PLN - beta1_tmp )))
      a_tmp <- a_true_list[[paste(betai, ni)]]
      
      if(is_matrix){
        error_a[[paste(ni,betai)]] <- data.frame(n = ni, beta = betai,
                                                 Poisson = sum(abs(tmp_a$Poisson - a_tmp[upper.tri(a_tmp, diag = T)])),
                                                 PG = sum(abs(tmp_a$PG - a_tmp[upper.tri(a_tmp, diag = T)] )),
                                                 PIG = sum(abs(tmp_a$PIG - a_tmp[upper.tri(a_tmp, diag = T)] )),
                                                 PLN = sum(abs(tmp_a$PLN - a_tmp[upper.tri(a_tmp, diag = T)] )))
      }else{
        error_a[[paste(ni,betai)]] <- data.frame(n = ni, beta = betai,
                                                 Poisson = sum(abs(tmp_a$Poisson - a_tmp)),
                                                 PG = sum(abs(tmp_a$PG - a_tmp )),
                                                 PIG = sum(abs(tmp_a$PIG - a_tmp )),
                                                 PLN = sum(abs(tmp_a$PLN - a_tmp )))
      }
      
      
      
      
    }
    
    
  }
  
  
  return(list(error_a = error_a, error_beta1 = error_beta1))
}


# Params and load all functions -----

source("utils_final.R")

ns <- c(100, 500, 1000, 5000)


# Ordinary test -----

beta1_ord <- list()
beta2_ord <- list()
a_ord  <- list()
liks_ord <- list()
beta1_true_ord <- list()
a_true_ord <- list()
betas_ord <- c(-0.1, 0.5, 2)
for(n in ns){
  print(n)
  for(beta in betas_ord){
    out <- generate_a_add(beta2_true = c(beta), type = 'PIG', dist = 'ordinary', n = n)
    agg_claims <- out$agg_claims
    claims <- out$claims
    z <- out$z
    locs <- out$locs
    years <- out$years
    X1 <- out$X1
    X2 <- out$X2
    w <- out$w
    exposure <- out$exposure
    y_latent <- out$y_latent
    a_true_ord[[paste(beta,n)]] <- out$A
    beta1_true_ord[[paste(beta,n)]] <- out$beta1
    
    
    # Poisson
    
    p_add_a <- Poisson_additive_a(claims, X1, locs, years, y_latent, exposure, max_itr = 100)
    p_lik <- sum(dpois(claims, lambda = p_add_a$mu, log = TRUE ))
    
    pg_add_a <- PG_EM_a(claims, X1, X2, locs, years, y_latent, TRUE, exposure, max_itr = 100)
    pg_lik <- sum(log_dnb(claims,pg_add_a$mu,pg_add_a$phi))
    
    pig_add_a <- PIG_EM_a(claims, X1, X2, locs, years, y_latent, TRUE, exposure, max_itr = 100)
    pig_lik <- sum(log_dPIG(claims,pig_add_a$mu,pig_add_a$phi))
    
    pln_add_a <- PLN_EM_a(claims, X1, X2, locs, years, y_latent, TRUE, exposure, max_itr = 100, tol = 1)
    pln_lik <- sum(log_dpln(claims,pln_add_a$mu,pln_add_a$phi))
    
    beta1_ord[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                             Poisson = p_add_a$beta1,
                                             PG = pg_add_a$beta1,
                                             PIG = pig_add_a$beta1,
                                             PLN = pln_add_a$beta1)
    
    beta2_ord[[paste(beta,n)]] <- data.frame(n = n, beta = beta,
                                             Poisson = NA,
                                             PG = pg_add_a$beta2,
                                             PIG = pig_add_a$beta2,
                                             PLN = pln_add_a$beta2)
    
    a_ord[[paste(beta,n)]] <- data.frame(n = n, beta = beta,
                                         Poisson = p_add_a$a,
                                         PG = pg_add_a$a,
                                         PIG = pig_add_a$a,
                                         PLN = pln_add_a$a)
    
    liks_ord[[paste(beta,n)]] <- data.frame(n = n, beta = beta,
                                            Poisson = p_lik,
                                            PG = pg_lik,
                                            PIG = pig_lik,
                                            PLN = pln_lik)
    
  }
  
  
}




errors_ord <- get_errors(beta1_ord, a_ord, beta1_true_ord, a_true_ord, betas_ord, ns)

error_a <- errors_ord$error_a
error_beta1 <- errors_ord$error_beta1



ggplot(gather(bind_rows(liks_ord),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Log-likelihood") +
  theme_minimal(base_size = 18)


ggplot(gather(bind_rows(error_a),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("A error") +
  theme_minimal(base_size = 18)




ggplot(gather(bind_rows(error_beta1),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("beta error") +
  theme_minimal(base_size = 18)




# ZIP ------

beta1_zip <- list()
beta2_zip <- list()
a_zip  <- list()
liks_zip <- list()
beta1_true_zip <- list()
a_true_zip <- list()

betas_zip <- c(-0.3, 0.2, 2)
for(n in ns){
  print(n)
  for(beta in betas_zip){
    print(beta)
    out <- generate_a_mult(beta2_true = c(beta), type = 'PG', dist = 'zip', n = n)
    agg_claims <- out$agg_claims
    claims <- out$claims
    z <- out$z
    locs <- out$locs
    years <- out$years
    X1 <- out$X1
    X2 <- out$X2
    w <- out$w
    exposure <- out$exposure
    y_latent <- out$y_latent
    a_true_zip[[paste(beta,n)]] <- out$A
    beta1_true_zip[[paste(beta,n)]] <- out$beta1
    
    
    # Poisson
    
    print("poisson")
    p_mult_a <- ZIP_a_all_in_one(claims, X1, locs, years, y_latent, FALSE, exposure, max_itr = 100)
    p_lik <- sum(log(dzip(claims,  p_mult_a$mu, p_mult_a$prop )))
    print("pg")
    pg_mult_a <- ZINB_a_all_in_one(claims, X1, X2, locs, years, y_latent, FALSE, exposure, max_itr = 100)
    pg_lik <- sum(log_zinb(claims,pg_mult_a$mu,exp(pg_mult_a$beta2), pg_mult_a$prop))
    
    print("pig")
    pig_mult_a <- ZIPIG_a_all_in_one(claims, X1, X2, locs, years, y_latent, FALSE, exposure, max_itr = 100)
    pig_lik <- sum(log_zipig(claims,pig_mult_a$mu,exp(pig_mult_a$beta2),pig_mult_a$prop))
    print("pln")
    pln_mult_a <- ZIPLN_a(claims, X1, X2, locs, years, y_latent, exposure, FALSE, max_itr = 30, tol = 1)
    pln_lik <- sum(log_dzipln(claims,pln_mult_a$mu, exp(pln_mult_a$beta2), pln_mult_a$prop))
    
    beta1_zip[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                             Poisson = p_mult_a$beta,
                                             PG = pg_mult_a$beta1,
                                             PIG = pig_mult_a$beta1,
                                             PLN = pln_mult_a$beta1)
    
    beta2_zip[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                             Poisson = NA,
                                             PG = pg_mult_a$beta2,
                                             PIG = pig_mult_a$beta2,
                                             PLN = pln_mult_a$beta2)
    
    a_zip[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                         Poisson = p_mult_a$a,
                                         PG = pg_mult_a$a,
                                         PIG = pig_mult_a$a,
                                         PLN = pln_mult_a$a)
    
    liks_zip[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                            Poisson = p_lik,
                                            PG = pg_lik,
                                            PIG = pig_lik,
                                            PLN = pln_lik)
    
  }
  
  
}



errors_zip <- get_errors(beta1_zip, a_zip, beta1_true_zip, a_true_zip, betas_zip, ns)


error_a <- errors_zip$error_a
error_beta1 <- errors_zip$error_beta1



ggplot(gather(bind_rows(liks_zip),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Log-likelihood") +
  theme_minimal(base_size = 18)


ggplot(gather(bind_rows(error_a),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("A error") +
  theme_minimal(base_size = 18)




ggplot(gather(bind_rows(error_beta1),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("beta error") +
  theme_minimal(base_size = 18)

# Hurdle -----


beta1_hurdle <- list()
beta2_hurdle <- list()
a_hurdle  <- list()
liks_hurdle <- list()
beta1_true_hurdle <- list()
a_true_hurdle <- list()
beta1_true_hurdle <- list()

betas_hurdle <- c(-1, -0.1, 0)
for(n in c(5000)){
  print(n)
  for(beta in betas_hurdle){
    print(n)
    print(beta)
    out <- generate_psi_mult(beta2_true = c(beta), type = 'PLN', dist = 'hurdle', n = n)
    agg_claims <- out$agg_claims
    claims <- out$claims
    z <- out$z
    locs <- out$locs
    years <- out$years
    X1 <- out$X1
    X2 <- out$X2
    w <- out$w
    exposure <- out$exposure
    y_latent <- out$y_latent
    a_true_hurdle[[paste(beta,n)]] <- out$psi
    beta1_true_hurdle[[paste(beta,n)]] <- out$beta1
    
    
    # Poisson
    
    print("Poisson")
    p_mult_a <- HP_psi_all_in_one(claims, X1, locs, years, y_latent, w, FALSE, exposure, max_itr = 100)
    p_lik <- sum(log(dhurdle(claims,  p_mult_a$mu, p_mult_a$prop )+1e-8))
    
    
    print("PG")
    pg_mult_a <- HPG_psi(claims, X1, X2, locs, years, y_latent, w, FALSE, exposure, max_itr = 20, tol = 0.1)
    pg_lik <- sum(log_dhpg(claims,pg_mult_a$mu,exp(pg_mult_a$beta2), pg_mult_a$prop))
    
    print("PIG")
    pig_mult_a <- HPIG_psi(claims, X1, X2, locs, years, y_latent, w, FALSE, exposure, max_itr = 20, tol = .1)
    pig_lik <- sum(log_dhpig(claims,pig_mult_a$mu,exp(pig_mult_a$beta2),pig_mult_a$prop))
    
    print("PLN")
    pln_mult_a <- HPLN_psi(claims, X1, X2, locs, years, y_latent, w, FALSE, exposure, max_itr = 20, tol = .1)
    pln_lik <- sum(log_dhpln(claims,pln_mult_a$mu, exp(pln_mult_a$beta2), pln_mult_a$prop))
    
    
    beta1_hurdle[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                                Poisson = p_mult_a$beta,
                                                PG = pg_mult_a$beta1,
                                                PIG = pig_mult_a$beta1,
                                                PLN = pln_mult_a$beta1)
    
    beta2_hurdle[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                                Poisson = NA,
                                                PG = pg_mult_a$beta2,
                                                PIG = pig_mult_a$beta2,
                                                PLN = pln_mult_a$beta2)
    
    a_hurdle[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                            Poisson = p_mult_a$psi,
                                            PG = pg_mult_a$psi,
                                            PIG = pig_mult_a$psi,
                                            PLN = pln_mult_a$psi)
    
    liks_hurdle[[paste(beta,n)]] <- data.frame(n = n, beta = beta, 
                                               Poisson = p_lik,
                                               PG = pg_lik,
                                               PIG = pig_lik,
                                               PLN = pln_lik)
    
  }
  
  
}



errors_hurdle <- get_errors(beta1_hurdle, a_hurdle, beta1_true_hurdle, a_true_hurdle, betas_hurdle, ns)


error_a <- errors_hurdle$error_a
error_beta1 <- errors_hurdle$error_beta1



ggplot(gather(bind_rows(liks_hurdle),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Log-likelihood") +
  theme_minimal(base_size = 18)


ggplot(gather(bind_rows(error_a),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Psi error") +
  theme_minimal(base_size = 18)


ggplot(gather(bind_rows(error_beta1),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("beta error") +
  theme_minimal(base_size = 18)



# EM vs ML PIG-----

beta1_ml_psi_add <- list()
beta1_ml_psi_mult <- list()
beta1_ml_a_add <- list()
beta1_ml_a_mult <- list()

beta2_ml_psi_add <- list()
beta2_ml_psi_mult <- list()
beta2_ml_a_add <- list()
beta2_ml_a_mult <- list()

lik_ml_psi_add <- list()
lik_ml_psi_mult <- list()
lik_ml_a_add <- list()
lik_ml_a_mult <- list()

psi_ml_psi_add <- list()
psi_ml_psi_mult <- list()
a_ml_a_add <- list()
a_ml_a_mult <- list()


betas_ord <- c(-0.1, 0.5, 2)
ns <- c(100, 500, 1000, 5000)
for(n in ns){
  print(n)
  for(beta in betas_ord){
    print(beta)
    # PSI add
    print("PSI add")
    out_psi_add <- generate_psi_add(beta2_true =beta, type = 'PIG', dist = 'ordinary', n = n)
    agg_claims <- out_psi_add$agg_claims
    claims <- out_psi_add$claims
    z <- out_psi_add$z
    locs <- out_psi_add$locs
    years <- out_psi_add$years
    X1 <- out_psi_add$X1
    X2 <- out_psi_add$X2
    w <- out_psi_add$w
    exposure <- out_psi_add$exposure
    y_latent <- out_psi_add$y_latent
    
    
    
    out_PIG_ML_psi_add <- PIG_ML_psi(claims, X1,X2, TRUE, y_latent, years, locs, w, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    out_PIG_EM_psi_add <- PIG_EM_psi(claims, X1,X2, locs, years, y_latent, w, TRUE, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    beta1_ml_psi_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                    ML = out_PIG_ML_psi_add$beta1, 
                                                    EM = out_PIG_EM_psi_add$beta1, 
                                                    par = out_psi_add$beta1)
    
    beta2_ml_psi_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                    ML = out_PIG_ML_psi_add$beta2, 
                                                    EM = out_PIG_EM_psi_add$beta2, 
                                                    par = out_psi_add$beta2)
    psi_ml_psi_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                  ML = out_PIG_ML_psi_add$psi, 
                                                  EM = out_PIG_EM_psi_add$psi, 
                                                  par = out_psi_add$psi)
    lik_ml_psi_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n,
                                                  ML = sum(log_dPIG(claims, out_PIG_ML_psi_add$mu,out_PIG_ML_psi_add$phi)),
                                                  EM =sum(log_dPIG(claims, out_PIG_EM_psi_add$mu,out_PIG_EM_psi_add$phi)))
    
    
    # PSI mult
    print("PSI mult")
    out_psi_mult <- generate_psi_add(beta2_true =beta, type = 'PIG', dist = 'ordinary', n = n)
    agg_claims <- out_psi_mult$agg_claims
    claims <- out_psi_mult$claims
    z <- out_psi_mult$z
    locs <- out_psi_mult$locs
    years <- out_psi_mult$years
    X1 <- out_psi_mult$X1
    X2 <- out_psi_mult$X2
    w <- out_psi_mult$w
    exposure <- out_psi_mult$exposure
    y_latent <- out_psi_mult$y_latent
    
    
    out_PIG_ML_psi_mult <- PIG_ML_psi(claims, X1,X2, FALSE, y_latent, years, locs, w, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    out_PIG_EM_psi_mult <- PIG_EM_psi(claims, X1,X2, locs, years, y_latent, w, FALSE, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    beta1_ml_psi_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                     ML = out_PIG_ML_psi_mult$beta1, 
                                                     EM = out_PIG_EM_psi_mult$beta1, 
                                                     par = out_psi_mult$beta1)
    
    beta2_ml_psi_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                     ML = out_PIG_ML_psi_mult$beta2, 
                                                     EM = out_PIG_EM_psi_mult$beta2, 
                                                     par = out_psi_mult$beta2)
    psi_ml_psi_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                   ML = out_PIG_ML_psi_mult$psi, 
                                                   EM = out_PIG_EM_psi_mult$psi, 
                                                   par = out_psi_mult$psi)
    lik_ml_psi_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n,
                                                   ML = sum(log_dPIG(claims, out_PIG_ML_psi_mult$mu,out_PIG_ML_psi_mult$phi)),
                                                   EM =sum(log_dPIG(claims, out_PIG_EM_psi_mult$mu,out_PIG_EM_psi_mult$phi)))
    
    
    # a add
    print("a add")
    out_a_add <- generate_a_add(beta2_true =beta, type = 'PIG', dist = 'ordinary', n = n)
    agg_claims <- out_a_add$agg_claims
    claims <- out_a_add$claims
    z <- out_a_add$z
    locs <- out_a_add$locs
    years <- out_a_add$years
    X1 <- out_a_add$X1
    X2 <- out_a_add$X2
    w <- out_a_add$w
    exposure <- out_a_add$exposure
    y_latent <- out_a_add$y_latent
    
    
    out_PIG_ML_a_add <- PIG_ML_a(claims, X1,X2, TRUE, y_latent, years, locs, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    out_PIG_EM_a_add <- PIG_EM_a(claims, X1,X2, locs, years, y_latent, TRUE, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    beta1_ml_a_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                  ML = out_PIG_ML_a_add$beta1, 
                                                  EM = out_PIG_EM_a_add$beta1, 
                                                  par = out_a_add$beta1)
    
    beta2_ml_a_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                  ML = out_PIG_ML_a_add$beta2, 
                                                  EM = out_PIG_EM_a_add$beta2, 
                                                  par = out_a_add$beta2)
    a_ml_a_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                              ML = out_PIG_ML_a_add$a, 
                                              EM = out_PIG_EM_a_add$a, 
                                              par = out_a_add$a)
    lik_ml_a_add[[paste(n,beta)]] <- data.frame(beta = beta, n = n,
                                                ML = sum(log_dPIG(claims, out_PIG_ML_a_add$mu,out_PIG_ML_a_add$phi)),
                                                EM =sum(log_dPIG(claims, out_PIG_EM_a_add$mu,out_PIG_EM_a_add$phi)))
    
    
    # a mult
    print("a mult")
    out_a_mult <- generate_a_add(beta2_true =beta, type = 'PIG', dist = 'ordinary', n = n)
    agg_claims <- out_a_mult$agg_claims
    claims <- out_a_mult$claims
    z <- out_a_mult$z
    locs <- out_a_mult$locs
    years <- out_a_mult$years
    X1 <- out_a_mult$X1
    X2 <- out_a_mult$X2
    w <- out_a_mult$w
    exposure <- out_a_mult$exposure
    y_latent <- out_a_mult$y_latent
    
    
    out_PIG_ML_a_mult <- PIG_ML_a(claims, X1,X2, FALSE, y_latent, years, locs, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    out_PIG_EM_a_mult <- PIG_EM_a(claims, X1,X2, locs, years, y_latent, FALSE, exposure = exposure,  max_itr = 1000, z= "", tol = 1e-3)
    beta1_ml_a_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                   ML = out_PIG_ML_a_mult$beta1, 
                                                   EM = out_PIG_EM_a_mult$beta1, 
                                                   par = out_a_mult$beta1)
    
    beta2_ml_a_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                                   ML = out_PIG_ML_a_mult$beta2, 
                                                   EM = out_PIG_EM_a_mult$beta2, 
                                                   par = out_a_mult$beta2)
    a_ml_a_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n, 
                                               ML = out_PIG_ML_a_mult$a, 
                                               EM = out_PIG_EM_a_mult$a, 
                                               par = out_a_mult$a)
    lik_ml_a_mult[[paste(n,beta)]] <- data.frame(beta = beta, n = n,
                                                 ML = sum(log_dPIG(claims, out_PIG_ML_a_mult$mu,out_PIG_ML_a_mult$phi)),
                                                 EM =sum(log_dPIG(claims, out_PIG_EM_a_mult$mu,out_PIG_EM_a_mult$phi)))
    
  }
}







ggplot(gather(bind_rows(lik_ml_a_mult),"Likelihood", "value",-n,-beta)) + geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Log-likelihood") +
  theme_minimal(base_size = 18)

bind_rows(beta2_ml_a_mult) %>% 
  mutate(ML_error = abs(ML-par),
         EM_error = abs(EM-par)) %>% 
  select(beta, n, ML_error, EM_error) %>% 
  gather("Likelihood", "value",-n,-beta) %>% 
  ggplot()+ geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Dispersion error") +
  theme_minimal(base_size = 18)



bind_rows(a_ml_a_mult) %>% 
  group_by(beta, n) %>% 
  mutate(ML_error = sum(abs(ML-par)),
         EM_error = sum(abs(EM-par))) %>% 
  select(beta, n, ML_error, EM_error) %>% 
  gather("Likelihood", "value",-n,-beta) %>% 
  ggplot()+ geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("beta error") +
  theme_minimal(base_size = 18)


bind_rows(beta1_ml_a_mult) %>% 
  group_by(beta, n) %>% 
  mutate(ML_error = sum(abs(ML-par)),
         EM_error = sum(abs(EM-par))) %>% 
  select(beta, n, ML_error, EM_error) %>% 
  gather("Likelihood", "value",-n,-beta) %>% 
  ggplot()+ geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("A error") +
  theme_minimal(base_size = 18)



bind_rows(a_ml_a_mult) %>% 
  group_by(beta, n) %>% 
  mutate(ML_error = sum(abs(ML-par)),
         EM_error = sum(abs(EM-par))) %>% 
  select(beta, n, ML_error, EM_error) %>% 
  gather("Likelihood", "value",-n,-beta) %>% 
  ggplot()+ geom_line(aes(x = log(n), y = value, color = Likelihood)) +
  facet_wrap(~beta) +
  ylab("Likelihood") +
  theme_minimal(base_size = 18)







