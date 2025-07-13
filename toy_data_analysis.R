

# Get data -----
library(tidyverse)
set.seed(123)



# ---------------- Step 2 Data cleaning ----

# Toy data 
data_model <- read.csv("toy_data.csv") #read.csv("BVDatabase.csv", sep = ";")
data_model %>% group_by(distr_code) %>% summarise(TExposure = sum(TExposure))

unique(data_model[, c("distr")])

# Spatial adjacency matrix - used if we learn spatial weights/psi
A <- read_rds("adj_w.rds")
# Adjacency matrix is alphabetric w.r.t distr, also see https://en.wikipedia.org/wiki/Regions_of_Greece
# North aegan and south Agan not in adjacency matrox

A <- A[-c(8,10), -c(8,10)]

# Map districtics to districts number
strings <- sort(unique(data_model$distr))
ds <- 1:length(strings)
names(ds) <- strings
data_model$distr_num <- ds[data_model$distr]

# Map year to time starting from 1
data_model$t <- data_model$Year - min(data_model$Year)+1


# policies with low exposure
nrow(data_model)
data_model <- data_model  %>% filter( TExposure>0.5)
nrow(data_model)


# Claim Count - all claims
data_model$ClaimCnt = data_model$TClaimsNo + data_model$ClaimNb_Flood + data_model$ClaimNb_Other + data_model$TClaimsNo



# Bin continuous variables
data_model$builtYear_cut <- cut(data_model$builtYear, breaks = c(-Inf, 1980, 2000, Inf), right = F )
data_model$nrOfFloors_cut <-  pmax(pmin(data_model$nrOfFloors, 3),0)

data_model = data_model %>% filter(Construction_Material_C %in% c("C1", "C2"))
data_model = data_model %>% filter(FloorsType_C %in% c("C1", "C2"))

# Generate risk groups - This is the data we will model
per_risks <- data_model %>% group_by(t,Year ,distr_num, builtYear_cut, nrOfFloors_cut, Construction_Material_C, FloorsType_C) %>% 
  summarise(ClaimCnt = sum(ClaimCnt),
            TExposure = sum(TExposure),
            TIncurred = sum(TIncurred),
            Freq = sum(ClaimCnt)/ sum(TExposure))


# --------- Step 2 create spatial data process -----------





# Claim Count -  per year and per locations
per_year_loc <- data_model %>% group_by(t, distr_num) %>% 
  summarise(ClaimCnt = sum(ClaimCnt),
            TExposure = sum(TExposure),
            TIncurred = sum(TIncurred),
            Freq = sum(ClaimCnt)/ sum(TExposure))

# Global frequency
freq_global <- sum(data_model$ClaimCnt)/sum(data_model$TExposure)

per_year_loc$isactive <- (per_year_loc$ClaimCnt > qpois(0.3, freq_global*per_year_loc$TExposure))*1
# Claim Cnt spatial
agg_claims <- spread(per_year_loc[, c("t", "distr_num", "isactive")], t, isactive) 
agg_claims <- agg_claims[,2:ncol(agg_claims)]
names(agg_claims) <- 1:ncol(agg_claims)
agg_claims <- cbind(agg_claims[,1], agg_claims)
agg_claims[is.na(agg_claims)] <- 0


# ---------------- Step 2 Split train test ----
set.seed(123)
train_ind <- sample(1:nrow(per_risks), 0.8*nrow(per_risks) )

data_model_train <- per_risks[train_ind, ]

# check if we have all time/location
length(unique(data_model_train$t)) == length(unique(per_risks$t))
length(unique(data_model_train$distr_num)) == length(unique(per_risks$distr_num))

data_model_train <-  per_risks[-train_ind, ]


#data_model_train <- data_model[train_ind,]
data_model_test <- per_risks[-train_ind, ]


X_mat_train <-  model.matrix(ClaimCnt~ as.factor(builtYear_cut)+as.factor(nrOfFloors_cut)+ 
                               as.factor(Construction_Material_C)+ 
                               as.factor(FloorsType_C), data = data_model_train)



years_train <- data_model_train$t # years go from 1 to Tend
years_test <- data_model_test$t
Tend <- max(years_train)

locs_train <- data_model_train$distr_num
Y_claims_train <- data_model_train$ClaimCnt
offset_train <- data_model_train$TExposure

locs_test <- data_model_test$distr_num
Y_claims_test <- data_model_test$ClaimCnt
offset_test <- data_model_test$TExposure



#--------- Model -------------

#--------- Ordinary Poisson-------------



file_name <- "toy_data_analysis.RData"
source("poisson.R")
poisson_models <- list()
for(mixing in c("gamma", "ig", "ln")){
  
  if(mixing == "gamma"){
    beta2_start <- 0.5
  }else if(mixing == "ig"){
    beta2_start <- 0.5
  }else if(mixing == "ln"){
    beta2_start <- 0
  }
  
  print(mixing)
  
  poisson_models[[mixing]] <- list()
  
  print(paste0(mixing, " - ","psi-add" ))
  poisson_models[[mixing]][["psi-add"]]  <- Poisson_mixed(Y_claims_train, X_mat_train, locs_train, years_train, agg_claims, A, TRUE, "learn_psi", lambda = 0,
                                                          exposure = offset_train, max_itr = 200, mixing_var = mixing, nr_em = 200, Q_tol = 0,
                                                          verbose = 1, sgd = FALSE, batch_size = 100, param_tol = 0, a_known = FALSE, beta2_start = beta2_start)
  save.image(file_name)
  print(paste0(mixing, " - ","psi-mult" ))
  poisson_models[[mixing]][["psi-mult"]] <- Poisson_mixed(Y_claims_train, X_mat_train, locs_train, years_train, agg_claims, A, FALSE, "learn_psi", lambda = 0,
                                                          exposure = offset_train, max_itr = 200, mixing_var = mixing, nr_em = 200, Q_tol = 0,
                                                          verbose = 1, sgd = FALSE, batch_size = 100, param_tol = 0, a_known = FALSE,beta2_start = beta2_start)
  
  save.image(file_name)
  print(paste0(mixing, " - ","a-add" ))
  poisson_models[[mixing]][["a-add"]] <- Poisson_mixed(Y_claims_train, X_mat_train, locs_train, years_train, agg_claims, A, TRUE, "learn_graph", lambda = 0,
                                                       exposure = offset_train, max_itr = 200, mixing_var = mixing, nr_em = 200, Q_tol = 0,
                                                       verbose = 1, sgd = FALSE, batch_size = 200, param_tol = 0, a_known = FALSE,beta2_start = beta2_start)
  save.image(file_name)
  print(paste0(mixing, " - ","a-add" ))
  poisson_models[[mixing]][["a-mult"]] <- Poisson_mixed(Y_claims_train, X_mat_train, locs_train, years_train, agg_claims, A, FALSE, "learn_graph", lambda = 0,
                                                        exposure = offset_train, max_itr = 200, mixing_var = mixing, nr_em = 200, Q_tol = 0,
                                                        verbose = 1, sgd = FALSE, batch_size = 200, param_tol = 0, a_known = FALSE,beta2_start = beta2_start)
  save.image(file_name)
}




#--------- ZIP Poisson-------------
source("zip_mixed.R")
zip_models <- list()
for(mixing in c("gamma", "ig", "ln")){
  
  if(mixing == "gamma"){
    beta2_start <- 0.5
  }else if(mixing == "ig"){
    beta2_start <- 0.5
  }else if(mixing == "ln"){
    beta2_start <- 0
  }
  

    zip_models[[mixing]] <- list()
    
    print(paste0(mixing, " - ","psi-add" ))
    zip_models[[mixing]][["psi-add"]] <- zip_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A,
                                                   exposure = offset_train, control_list = list(psi = list(lr = 0.00001),
                                                                                                beta_phi  = list(lr = 0.0001)),
                                                   model_type = "learn_psi", additive = TRUE, mixing = mixing,  Emethod = "integration",
                                                   n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                   optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                   batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
    save.image(file_name)
    print(paste0(mixing, " - ","psi-mult" ))
    zip_models[[mixing]][["psi-mult"]] <- zip_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A,
                                                    exposure = offset_train, control_list = list(psi = list(lr = 0.00001),
                                                                                                 beta_phi  = list(lr = 0.0001)),
                                                    model_type = "learn_psi", additive = FALSE, mixing = mixing,  Emethod = "integration",
                                                    n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                    optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                    batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
    
  
  
  save.image(file_name)
  print(paste0(mixing, " - ","a-add" ))
  zip_models[[mixing]][["a-add"]] <- zip_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A, 
                                               exposure = offset_train, control_list = list(a = list(lr = 0.000001),
                                                                                            beta_phi  = list(lr = 0.0001)),
                                               model_type = "learn_graph", additive = TRUE, mixing = mixing,  Emethod = "integration",
                                               n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                               optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                               batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
  save.image(file_name)
  print(paste0(mixing, " - ","a-add" ))
  zip_models[[mixing]][["a-mult"]] <- zip_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A,
                                                exposure = offset_train, control_list = list(a = list(lr = 0.000001),
                                                                                             beta_phi  = list(lr = 0.0001)),
                                                model_type = "learn_graph", additive = FALSE, mixing = mixing,  Emethod = "integration",
                                                n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
  save.image(file_name)
}



#--------- Hurdle Poisson-------------
source("hurdle_mixed.R")
hurdle_models <- list()
for(mixing in c("gamma", "ig", "ln")){
  
  if(mixing == "gamma"){
    beta2_start <- 0.5
  }else if(mixing == "ig"){
    beta2_start <- 0.5
  }else if(mixing == "ln"){
    beta2_start <- 0
  }
  
  print(mixing)
  
  hurdle_models[[mixing]] <- list()
  
  print(paste0(mixing, " - ","psi-add" ))
  hurdle_models[[mixing]][["psi-add"]] <- hurdle_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A, 
                                                       exposure = offset_train,   control_list = list(psi = list(lr = 0.00001),
                                                                                                      beta_phi  = list(lr = 0.0001)),
                                                       model_type = "learn_psi", additive = TRUE, mixing = mixing,  Emethod = "integration",
                                                       n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                       optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                       batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
  save.image(file_name)
  print(paste0(mixing, " - ","psi-mult" ))
  hurdle_models[[mixing]][["psi-mult"]] <- hurdle_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A, 
                                                        exposure = offset_train,   control_list = list(psi = list(lr = 0.00001),
                                                                                                       beta_phi  = list(lr = 0.0001)),
                                                        model_type = "learn_psi", additive = FALSE, mixing = mixing,  Emethod = "integration",
                                                        n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                        optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                        batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
  
  save.image(file_name)
  print(paste0(mixing, " - ","a-add" ))
  hurdle_models[[mixing]][["a-add"]] <- hurdle_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A, 
                                                     exposure = offset_train, control_list = list(a = list(lr = 0.000001),
                                                                                                  beta_phi  = list(lr = 0.0001)),
                                                     model_type = "learn_graph", additive = TRUE, mixing = mixing,  Emethod = "integration",
                                                     n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                     optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                     batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
  save.image(file_name)
  print(paste0(mixing, " - ","a-add" ))
  hurdle_models[[mixing]][["a-mult"]] <- hurdle_mixed(claims = Y_claims_train, X = X_mat_train, years = years_train, locs = locs_train, agg_claims = agg_claims, A = A, 
                                                      exposure = offset_train, control_list = list(a = list(lr = 0.000001),
                                                                                                   beta_phi  = list(lr = 0.0001)),
                                                      model_type = "learn_graph", additive = FALSE, mixing = mixing,  Emethod = "integration",
                                                      n_iter = 200, lambda = 0, optimizer_beta = "gd", optimizer_psi = "gd",
                                                      optimizer_a = "gd", optimizer_pi = "gd", optimizer_beta_phi = "gd", sgd = FALSE,
                                                      batch_size = 100, verbose = 1, do_optim = FALSE, calc_se = TRUE, a_known = FALSE, beta2_start = beta2_start)
  save.image(file_name)
}


