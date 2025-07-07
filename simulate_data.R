
#########################
# 1. Simulate Data
#########################
source("utils.R")


# Function to generate a sparse matrix with positive numbers
generateSparseMatrix <- function(nrows, ncols, density, minVal = 0.1, maxVal = 1) {
  set.seed(123)
  library(Matrix)
  # Ensure that density is between 0 and 1
  if (density < 0 || density > 1) {
    stop("Density must be between 0 and 1.")
  }
  
  # Use rsparsematrix to create the sparse matrix
  # The rand.x function generates positive numbers (from runif)
  sparse_mat <- rsparsematrix(nrow = nrows,
                              ncol = ncols,
                              density = density,
                              symmetric = TRUE,
                              rand.x = function(n) runif(n, min = minVal, max = maxVal))
  
  diag(sparse_mat) <- max(sparse_mat)*runif(nrows, 1.2, 1.5)
  return(sparse_mat)
}


simulate_claims <- function(n, years, spatial_type, additive, mixing, area = 5, model_type = "poisson", density = 0.4, 
                            exposure_lambda = 0, seed = 123, agg_spar = 0.4, order_regions = FALSE){
  
  A1_true <- as.matrix(generateSparseMatrix(area, area, density))
  set.seed(123)
  psi_true <- runif(area, 1, 3)
  set.seed(seed)
  
  # store claim data
  if(order_regions){
    locs <- rep(1, n)
    for(i in 1:(n)){
      locs[i] <- ((i - 1) %% area) + 1
    }
    
  }else{
    while(TRUE){
      locs <- sample(1:area, n, replace = TRUE)
      if(length(unique(locs)) == area){
        break
      }
    }  
  }

  
  claims <- matrix(0, nrow = n, ncol = years+1)
  mu <- matrix(0, nrow = n, ncol = years)
  phi <- matrix(0, nrow = n, ncol = years)
  z <- matrix(0, nrow = n, ncol = years)
  years_vec <- rep(1:years, each = n)
  locs_vec <- c()
  
  exposure <- rpois(years*n, lambda = exposure_lambda) + 1
  
  y_latent <- matrix(rnorm(years*area), nrow = area, ncol = years)
  y_latent[y_latent < agg_spar] <- 0
  y_latent[y_latent != 0] <- 1
  
  # Generate known parameters
  X_mat1 <- cbind(rep(1, n*years),
                  rnorm(n*years, 0.1, sd = 0.1),
                  rnorm(n*years, -0.3, sd = 0.1))
  X_mat2 <- as.matrix(rep(1, n*years))
  beta1_true <- c(1, -1, 1)
  
  
  

  if(model_type == "hurdle"){
    prop_true <- 0.4
  }else{
    prop_true <- 0.7
  }
  p <- length(unique(locs))
  
  if(mixing == "none"){
    beta2_true <- 0
  } else if(mixing == "gamma"){
    beta2_true <- 0.5
  } else if(mixing == "ig"){
    beta2_true <- 0.5
  } else if(mixing == "ln"){
    beta2_true <- 0
  }
  
  # Generate data over time
  for(t in 1:years){
    if(spatial_type == "psi"){
      if(additive){
        mu[, t] <- (exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) +
                      psi_true[locs] * ( (A1_true %*% y_latent[, t]) )[locs]) *
          exposure[((t-1)*n + 1):(t*n)]
      } else {
        mu[, t] <- (1 + psi_true[locs] * ( (A1_true %*% y_latent[, t]) )[locs]) *
          exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
          exposure[((t-1)*n + 1):(t*n)]
      }
    } else if(spatial_type == "graph") {
      if(additive){
        mu[, t] <- (exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) +
                      ( (A1_true %*% y_latent[, t]))[locs]) *
          exposure[((t-1)*n + 1):(t*n)]
      } else {
        mu[, t] <- (1 + ( (A1_true %*% y_latent[, t]))[locs]) *
          exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
          exposure[((t-1)*n + 1):(t*n)]
      }
    } else if(spatial_type == "none") {
      mu[, t] <- exp(X_mat1[((t-1)*n + 1):(t*n),] %*% beta1_true) *
        exposure[((t-1)*n + 1):(t*n)]
    }
    
    if(mixing == "none"){
      z[, t] <- 1
    } else if(mixing == "gamma"){
      z[, t] <- rgamma(n, exp(beta2_true), exp(beta2_true))
    } else if(mixing == "ig"){
      z[, t] <- actuar::rinvgauss(n, mean = 1, shape = exp(beta2_true)^2)
    } else if(mixing == "ln"){
      z[, t] <- rlnorm(n, meanlog = -exp(beta2_true)^2/2, sdlog = exp(beta2_true))
    }
    
    u <- runif(n)
    if(model_type == "poisson"){
      claims[, t+1] <- rpois(n, mu[, t] * z[, t])
    }else if(model_type == "zip"){
      claims[, t+1] <- (u >= prop_true) * rpois(n, mu[, t] * z[, t])
    } else if(model_type == "hurdle"){
      for(i in 1:n){
        # Decide if observation is a hurdle (y=0) or positive outcome (y>0)
        if (runif(1) <= prop_true) {
          claims[i,t+1] <- 0
        } else {
          lambda_val <- mu[i,t]* z[i,t]
          # Rejection sampling for truncated Poisson (ensure y > 0)
          repeat {
            y_sample <- rpois(1, lambda = lambda_val)
            if (y_sample > 0) break
          }
          claims[i,t+1] <- y_sample
        }
      }
    }

    locs_vec <- c(locs_vec, locs)
  }
  
  claims <- as.vector(claims[,-1])
  z <- as.vector(z)
  mu <- as.vector(mu)
  phi <- as.vector(phi)
  
  return(list(X = X_mat1, claims = claims, z = z, mu = mu, phi = phi, A = A1_true, psi = psi_true, 
              beta2 = beta2_true, beta1 = beta1_true, agg_claims = y_latent,
              years = years_vec, locs = locs_vec, exposure = exposure, 
              a = A1_true[upper.tri(A1_true, TRUE)] ))
}



evaluate_matrix_recovery <- function(A, A_est, threshold = 1e-3) {
  # Check if matrices have identical dimensions
  if (!all(dim(A) == dim(A_est))) {
    stop("A and A_est must have the same dimensions.")
  }
  
  # --- Method 1: Relative Reconstruction Error ---
  # Frobenius norm error between A and its estimate A_est
  relative_error <- norm(A - A_est, type = "F") / norm(A, type = "F")
  
  # --- Method 2: Link Recovery Metrics ---
  # In A, consider an entry to represent a "link" if it is nonzero
  true_links <- (A != 0)
  
  # For A_est, count an entry as a recovered link if its value exceeds the threshold.
  # (If negative values might occur and you only care about magnitude, you can use abs(A_est).)
  estimated_links <- (A_est > threshold)
  
  # Compute confusion counts
  TP <- sum(true_links & estimated_links)       # True links correctly uncovered
  FP <- sum(!true_links & estimated_links)      # False links (non-links predicted as links)
  FN <- sum(true_links & !estimated_links)      # Missed links
  TN <- sum(!true_links & !estimated_links)     # Correctly absent links
  
  # Calculate precision, recall, and F1 score
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  f1_score  <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall)
  else NA
  
  # Create a confusion matrix.
  # Rows are actual classes, and columns are predicted classes.
  confusion_matrix <- matrix(c(TP, FN, FP, TN), nrow = 2, byrow = TRUE)
  rownames(confusion_matrix) <- c("Actual Positive", "Actual Negative")
  colnames(confusion_matrix) <- c("Predicted Positive", "Predicted Negative")
  
  # Return all computed metrics as a list
  list(
    relative_error = relative_error,
    precision = precision,
    recall = recall,
    f1_score = f1_score,
    confusion_matrix = confusion_matrix
  )
}


plot_graphs <- function(A_est, A_true){
  A_est[abs(A_est) <= 1e-3 ]<- 0
  # suppose A_true and A_est are your two N×N 0/1 matrices
  df <- expand.grid(i = 1:nrow(A_true), j = 1:ncol(A_true)) %>%
    mutate(
      true = as.vector(A_true),
      est  = as.vector(A_est),
      type = case_when(
        abs(true) > 1e-3 & abs(est) > 1e-3 ~ "TP",
        abs(true) > 1e-3 & abs(est) <= 1e-3 ~ "FN",
        abs(true) <= 1e-3  & abs(est) > 1e-3 ~ "FP",
        TRUE                ~ "TN"
      )
    )
  
  true_false_graph <- ggplot(df, aes(x = j, y = i, fill = type)) +
    geom_tile(color = "grey80") +
    # flip the y-axis so row 1 appears at the top
    scale_y_reverse(breaks = 1:nrow(A_true)) +
    coord_equal() +
    scale_fill_manual(
      values = c(TP = "#1b9e77", FN = "#d95f02", FP = "#DE3163", TN = "white")
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Column (j)",
      y = "Row (i)",
      fill = "Edge type",
      title = "True vs. Estimated Adjacency"
    )
  
  # suppose A_true and A_est are your two N×N 0/1 matrices
  # turn them into a long data frame with a “which” indicator
  df_true <- expand.grid(i = 1:nrow(A_true), j = 1:ncol(A_true)) %>%
    mutate(type = "True", value = as.vector(A_true))
  
  df_est  <- expand.grid(i = 1:nrow(A_est), j = 1:ncol(A_est)) %>%
    mutate(type = "Estimated", value = as.vector(A_est))
  
  df_long <- bind_rows(df_true, df_est)
  
  # now plot with facet
  side_by_side <- ggplot(df_long, aes(x = j, y = i, fill = factor(value))) +
    geom_tile(color = "grey80") +
    scale_y_reverse(breaks = 1:nrow(A_true)) +
    coord_equal() +
    scale_fill_manual(
      values = c(`0` = "white", `1` = "#1b9e77"),
      labels = c("0" = "no edge", "1" = "edge")
    ) +
    facet_wrap(~type, ncol = 2) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Column (j)",
      y = "Row (i)",
      fill = "Presence",
      title = "True vs. Estimated Adjacency Matrices"
    )
  
  
  
  return(list(true_false_graph = true_false_graph, side_by_side = side_by_side))
  
}



plot_per_lambda <- function(model, ids, dens, A_compare = "50-5"){
  betas_p <- list()
  betas_mp <- list()
  betas_mp_known <- list()
  
  betas2_mp <- list()
  betas2_mp_known <- list()
  
  a_p <- list()
  a_mp <- list()
  a_mp_known <- list()
  
  a_p_info <- list()
  a_ig_info <- list()
  for(id in ids){
    
    tmp <- new.env()
    loaded_objs <- load(paste0("sim_data/", model, "_per_t_lambda", dens, "id", id, ".RData") , envir = tmp)
    my_list <- mget(loaded_objs, envir = tmp)
    
    
    betas_p[[as.character(id)]] <- bind_rows(my_list$beta_est_p)
    betas_p[[as.character(id)]]$id <- id
    mat <- do.call(rbind, strsplit(names(my_list$beta_est_p), "-", fixed = TRUE))
    betas_p[[as.character(id)]]$t <- as.numeric(mat[,1])
    betas_p[[as.character(id)]]$lambda <- as.numeric(mat[,2])
    
    betas_mp[[as.character(id)]] <- bind_rows(my_list$beta_est_ig)
    betas_mp[[as.character(id)]]$id <- id
    betas_mp[[as.character(id)]]$t <- as.numeric(mat[,1])
    betas_mp[[as.character(id)]]$lambda <- as.numeric(mat[,2])
    
    
    a_p[[as.character(id)]] <- as.data.frame(t(bind_rows(my_list$A_est_p)))
    a_p[[as.character(id)]]$id <- id
    a_p[[as.character(id)]]$t <- as.numeric(mat[,1])
    a_p[[as.character(id)]]$lambda <- as.numeric(mat[,2])
    
    a_mp[[as.character(id)]] <- as.data.frame(t(bind_rows(my_list$A_est_ig)))
    a_mp[[as.character(id)]]$id <- id
    a_mp[[as.character(id)]]$t <- as.numeric(mat[,1])
    a_mp[[as.character(id)]]$lambda <- as.numeric(mat[,2])
    
    
    a_p_info[[as.character(id)]] <- data.frame(t = mat[, 1], lambda = mat[, 2])
    a_p_info[[as.character(id)]]$f1 <- sapply(my_list$A_p_info, function(x) x$f1_score)
    a_p_info[[as.character(id)]]$recall <- sapply(my_list$A_p_info, function(x) x$recall)
    a_p_info[[as.character(id)]]$precision <- sapply(my_list$A_p_info, function(x) x$precision)
    a_p_info[[as.character(id)]]$relative_error <- sapply(my_list$A_p_info, function(x) x$relative_error)
    
    a_ig_info[[as.character(id)]] <- data.frame(t = mat[, 1], lambda = mat[, 2])
    a_ig_info[[as.character(id)]]$f1 <- sapply(my_list$A_ig_info, function(x) x$f1_score)
    a_ig_info[[as.character(id)]]$recall <- sapply(my_list$A_ig_info, function(x) x$recall)
    a_ig_info[[as.character(id)]]$precision <- sapply(my_list$A_ig_info, function(x) x$precision)
    a_ig_info[[as.character(id)]]$relative_error <- sapply(my_list$A_ig_info, function(x) x$relative_error)
    
    
    
    a_true <- my_list$sim$a
    beta_true <- my_list$sim$beta1
    beta2_true <- my_list$sim$beta2
    
  }
  
  betas_p <- bind_rows(betas_p)
  betas_mp <- bind_rows(betas_mp)
  a_p <- bind_rows(a_p)
  a_mp <- bind_rows(a_mp)
  a_p_info <- bind_rows(a_p_info)
  a_ig_info <- bind_rows(a_ig_info)
  
  
  
  # error 
  
  beta_true_matrix <- matrix(beta_true,
                             nrow = nrow(betas_p),
                             ncol = length(beta_true),
                             byrow = TRUE)
  
  betas_p$error <- rowSums(abs(betas_p[, c("X11", "X12", "X13")] - beta_true_matrix))
  betas_mp$error <- rowSums(abs(betas_mp[, c("X11", "X12", "X13")] - beta_true_matrix))
  
  
  a_true_matrix <- matrix(a_true,
                          nrow = nrow(a_p),
                          ncol = length(a_true),
                          byrow = TRUE)
  a_p$error <- rowSums(abs(a_p[ ,1:(10*11/2)] - a_true_matrix))
  a_mp$error <- rowSums(abs(a_mp[ ,1:(10*11/2)] - a_true_matrix))
  
  
  ret <- list()
  
  
  df1 <- betas_p %>% mutate(model = "Poisson")
  df2 <- betas_mp %>% mutate(model = "Mixed Poisson")
  
  # bind into one
  betas_all <- bind_rows(df1, df2)
  
  # dodged box‐plot
  ret$beta1_plot <- ggplot(betas_all %>% group_by(t, lambda, model) %>% summarise(y = mean(X11)), aes(x = t, y = y)) +
    # true parameter line + points
    geom_line(
      aes( color = as.factor(lambda))
    ) +
    geom_line(
      aes(x = t, y = 1, group = 1),
      inherit.aes = FALSE,
      color       = "black",
      size        = 1,
      alpha = 0.5
    ) +
    labs(
      x     = "Time (t)",
      y     = expression(beta[1] ~ " estimate"),
      fill  = "Model spec",
      title = expression(beta[1] ~ " over Time, by Model"),
      color  = "Lambda"
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    ) + 
    facet_wrap(~model)
  
  ret$beta1_plot
  
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
    ) + 
    facet_wrap(~ lambda)
  
  ret$beta_error_plot 
  
  df1 <- a_p %>% mutate(model = "Poisson")
  df2 <- a_mp %>% mutate(model = "Mixed Poisson")
  
  # bind into one
  as_all <- bind_rows(df1, df2)
  
  # dodged box‐plot
  ret$edge_plot <- ggplot(as_all %>% group_by(t, lambda, model) %>% summarise(y = mean(V8)), aes(x = T, y = y)) +
    # true parameter line + points
    geom_line(
      aes(x = t, y = a_true[8], group = 1),
      inherit.aes = FALSE,
      color       = "black",
      size        = 1,
      alpha = 1
    ) +
    geom_line(
      aes(x =t, color = as.factor(lambda))
    ) +
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
    ) + facet_wrap(~ model)
  
  ret$edge_plot
  
  # dodged box‐plot
  ret$edge_error_plot <- ggplot(as_all %>% filter(lambda <= 10), aes(x = factor(t), y = log(error), fill = model)) +
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
    ) + facet_wrap(~lambda)
  
  ret$edge_error_plot 
  
  
  
  a_p_info$model <- "Poisson"
  a_ig_info$model <- "Mixed Poisson"
  
  df_info <- a_ig_info# bind_rows(a_p_info, a_ig_info)
  
  df_info$t <- paste0("t = ",df_info$t)
  df_info$t <- factor( df_info$t, levels = c("t = 10","t = 50", "t = 100", "t = 200"), ordered  = TRUE)
  
  
  # dodged box‐plot
  ret$f1_plot <- ggplot(df_info, 
                        aes(x = as.numeric(lambda),    # make sure λ is numeric
                            y = f1, 
                            color = model, 
                            fill  = model,
                            group = model)) +
    # ribbon: mean ± SE
    stat_summary(fun.data = mean_se, 
                 geom     = "ribbon", 
                 alpha    = 0.2, 
                 color    = NA) +
    # mean line
    stat_summary(fun      = mean, 
                 geom     = "line", 
                 size     = 1) +
    labs(
      x     = expression(lambda),
      y     = "f1 value",
      color = "Model spec",
      fill  = "Model spec",
      title = "f1 value over λ, by Model, Time, and time points"
    ) +
    theme_minimal(base_size = 24) +
    theme(
      plot.title  = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 0.5)
    ) +
    facet_wrap(~ t, scales = "free_x")
  
  ret$f1_plot
  
  
  # dodged box‐plot
  ret$recall_plot <- ggplot(df_info, 
                            aes(x = as.numeric(lambda),    # make sure λ is numeric
                                y = recall, 
                                color = model, 
                                fill  = model,
                                group = model)) +
    # ribbon: mean ± SE
    stat_summary(fun.data = mean_se, 
                 geom     = "ribbon", 
                 alpha    = 0.2, 
                 color    = NA) +
    # mean line
    stat_summary(fun      = mean, 
                 geom     = "line", 
                 size     = 1) +
    labs(
      x     = expression(lambda),
      y     = "f1 value",
      color = "Model spec",
      fill  = "Model spec",
      title = "f1 value over λ, by Model, Time, and time points"
    ) +
    theme_minimal(base_size = 24) +
    theme(
      plot.title  = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 0.5)
    ) +
    facet_wrap(~ t, scales = "free_x")
  
  ret$recall_plot
  
  # dodged box‐plot
  ret$precision_plot <- ggplot(df_info, 
                               aes(x = as.numeric(lambda),    # make sure λ is numeric
                                   y = precision, 
                                   color = model, 
                                   fill  = model,
                                   group = model)) +
    # ribbon: mean ± SE
    stat_summary(fun.data = mean_se, 
                 geom     = "ribbon", 
                 alpha    = 0.2, 
                 color    = NA) +
    # mean line
    stat_summary(fun      = mean, 
                 geom     = "line", 
                 size     = 1) +
    labs(
      x     = expression(lambda),
      y     = "f1 value",
      color = "Model spec",
      fill  = "Model spec",
      title = "f1 value over λ, by Model, Time, and time points"
    ) +
    theme_minimal(base_size = 24) +
    theme(
      plot.title  = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 0.5)
    ) +
    facet_wrap(~ t, scales = "free_x")
  
  ret$precision_plot
  
  
  
  graph_plot_out <- plot_graphs(get_W_from_array(my_list$A_est_ig[[A_compare]], 10), my_list$sim$A)
  
  
  ret$side_by_side <- graph_plot_out$side_by_side
  
  ret$true_false_graph<- graph_plot_out$true_false_graph
  
  
  return(ret)
}


plot_per_t <- function(model, ids, dens){
  betas_p <- list()
  betas_mp <- list()
  betas_mp_known <- list()
  
  betas2_mp <- list()
  betas2_mp_known <- list()
  
  a_p <- list()
  a_mp <- list()
  a_mp_known <- list()
  
  time_p <- list()
  time_mp <- list()
  for(id in ids){
    
    tmp <- new.env()
    loaded_objs <- load(paste0("sim_data/", model, "_per_t", dens, "id", id, ".RData") , envir = tmp)
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
    
    
    time_p[[as.character(id)]] <- data.frame(time = sapply(my_list$time_p,  function(dt) as.numeric(dt, units = "mins")),
                                             id = id,
                                             t = as.numeric(names(my_list$beta_est_ig)))
    
    time_mp[[as.character(id)]] <- data.frame(time = sapply(my_list$time_ig,  function(dt) as.numeric(dt, units = "mins")),
                                             id = id,
                                             t = as.numeric(names(my_list$beta_est_ig)))
    
    a_true <- my_list$sim$a
    beta_true <- my_list$sim$beta1
    beta2_true <- my_list$sim$beta2
    
  }
  
  betas_p <- bind_rows(betas_p)
  betas_mp <- bind_rows(betas_mp)
  betas_mp_known <- bind_rows(betas_mp_known)
  a_p <- bind_rows(a_p)
  a_mp <- bind_rows(a_mp)
  time_mp <- bind_rows(time_mp)
  time_p <- bind_rows(time_p)
  
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
  ret$edge_plot <- ggplot(as_all, aes(x = factor(t), y = V8, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = a_true[8], group = 1),
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
  
  
  
  
  df1 <- time_mp %>% mutate(model = "Mixed Poisson")
  df2 <- time_p %>% mutate(model = "Poisson")
  times_df <- bind_rows(df1, df2)
  
  ret$times <- ggplot(times_df, aes(x = as.factor(t), y = log(time), fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2)+ 
    labs(
      x     = "Time points (t)",
      y     = latex2exp::TeX("$log$ complexity"),
      fill  = "Model spec",
      title = latex2exp::TeX("$a_{1,1}$ over time, by model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  ret$times
  
  df1 <- betas2_mp %>% mutate(model = "Mixed Poisson")
  df2 <- betas2_mp_known %>% mutate(model = "Mixed Poisson known A")
  
  # bind into one
  beta2_all <- bind_rows(df1, df2)
  
  # dodged box‐plot
  ret$beta2_plot <- ggplot(beta2_all, aes(x = factor(t), y = V1, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = beta2_true, group = 1),
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

plot_per_t_psi <- function(model, ids, dens){

  
  
  betas_p <- list()
  betas_mp <- list()
  betas_mp_known <- list()
  
  betas2_mp <- list()
  betas2_mp_known <- list()
  
  a_p <- list()
  a_mp <- list()
  a_mp_known <- list()
  
  time_p <- list()
  time_mp <- list()
  for(id in ids){
    
    tmp <- new.env()
    loaded_objs <- load(paste0("sim_data/", model, "_per_t", dens, "id", id, ".RData") , envir = tmp)
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
    
    time_p[[as.character(id)]] <- data.frame(time = sapply(my_list$time_p,  function(dt) as.numeric(dt, units = "mins")),
                                             id = id,
                                             t = as.numeric(names(my_list$beta_est_ig)))
    
    time_mp[[as.character(id)]] <- data.frame(time = sapply(my_list$time_ig,  function(dt) as.numeric(dt, units = "mins")),
                                              id = id,
                                              t = as.numeric(names(my_list$beta_est_ig)))
    
    a_true <- my_list$sim$psi
    beta_true <- my_list$sim$beta1
    beta2_true <- my_list$sim$beta2
    
  }
  
  betas_p <- bind_rows(betas_p)
  betas_mp <- bind_rows(betas_mp)
  a_p <- bind_rows(a_p)
  a_mp <- bind_rows(a_mp)
  time_mp <- bind_rows(time_mp)
  time_p <- bind_rows(time_p)
  
  betas2_mp <- bind_rows(betas2_mp)

  
  # error 
  
  beta_true_matrix <- matrix(beta_true,
                             nrow = nrow(betas_p),
                             ncol = length(beta_true),
                             byrow = TRUE)
  
  betas_p$error <- rowSums(abs(betas_p[, c("X11", "X12", "X13")] - beta_true_matrix))
  betas_mp$error <- rowSums(abs(betas_mp[, c("X11", "X12", "X13")] - beta_true_matrix))

  
  
  a_true_matrix <- matrix(a_true,
                          nrow = nrow(a_p),
                          ncol = length(a_true),
                          byrow = TRUE)
  a_p$error <- rowSums(abs(a_p[ ,1:10] - a_true_matrix))
  a_mp$error <- rowSums(abs(a_mp[ ,1:10] - a_true_matrix))
  
  
  ret <- list()
  
  
  df1 <- betas_p %>% mutate(model = "Poisson")
  df2 <- betas_mp %>% mutate(model = "Mixed Poisson")
  
  # bind into one
  betas_all <- bind_rows(df1, df2)
  
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
  ret$edge_plot <- ggplot(as_all, aes(x = factor(t), y = V1, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = a_true[1], group = 1),
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
      y     = latex2exp::TeX("$\\phi_{1}$ estimate"),
      fill  = "Model spec",
      title = latex2exp::TeX("$\\phi_{1}$ over time, by model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  
  
  # dodged box‐plot
  ret$edge_error_plot <- ggplot(as_all, aes(x = factor(t), y = error, fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2) + 
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  
  
  
  
  df1 <- time_mp %>% mutate(model = "Mixed Poisson")
  df2 <- time_p %>% mutate(model = "Poisson")
  times_df <- bind_rows(df1, df2)
  
  ret$times <- ggplot(times_df, aes(x = as.factor(t), y = log(time), fill = model)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", 
                 position = position_dodge(width = 0.8), 
                 width = 0.2)+ 
    labs(
      x     = "Time points (t)",
      y     = latex2exp::TeX("$log$ complexity"),
      fill  = "Model spec",
      title = latex2exp::TeX("$a_{1,1}$ over time, by model")
    )  +
    theme_minimal(base_size = 24) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      axis.text.x  = element_text(angle = 45, vjust = 0.5)
    )
  ret$times
  
  df1 <- betas2_mp %>% mutate(model = "Mixed Poisson")

  
  # bind into one
  beta2_all <- bind_rows(df1)
  
  # dodged box‐plot
  ret$beta2_plot <- ggplot(beta2_all, aes(x = factor(t), y = V1, fill = model)) +
    # true parameter line + points
    geom_line(
      aes(x = factor(t), y = beta2_true, group = 1),
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



plot_error <- function(my_list, idx1 , idx2, is_ordinary){
  
  if(is_ordinary){
    extra_id <- 0
  }else{
    extra_id <- 1
  }
  
  se_error_ig <- do.call(rbind,
                         mapply(
                           function(x, idx, y) {
                             # idx is a single integer (1, 2, 3, …)
                             mat <- x-y + 1 * diag(nrow(x)) * (idx == 1) + 1 * diag(nrow(x)) * (idx == 2) + 2 * diag(nrow(x)) * (idx == 3)
                             diag( solve(mat) )
                           },
                           my_list$H_ig_p,
                           seq_along(my_list$H_ig_p),
                           my_list$var_ig_p,
                           SIMPLIFY = FALSE
                         )
  )
  se_error_ig
  
  se_error_p <- do.call(rbind,
                        mapply(
                          function(x, idx) {
                            # idx is a single integer (1, 2, 3, …)
                            mat <- x  + .5 * diag(nrow(x)) * (idx == 1) + 2 * diag(nrow(x)) * (idx == 2) + 2 * diag(nrow(x)) * (idx == 3)
                            diag( solve(mat) )
                          },
                          my_list$H_p,
                          seq_along(my_list$H_p),
                          SIMPLIFY = FALSE
                        )
  )
  
  
  
  
  
  est1 <- sapply(my_list$beta_est_ig, function(x) x[1])
  se_err1 <- se_error_ig[,2+extra_id]
  est2 <- sapply(my_list$beta_est_ig, function(x) x[2])
  se_err2 <- se_error_ig[,3+extra_id]
  
  
  aest1 <- sapply(my_list$A_est_ig, function(x) x[idx1])
  ase_err1 <- se_error_ig[,5+ idx1-1 +extra_id]
  aest2 <- sapply(my_list$A_est_ig, function(x) x[idx2])
  ase_err2 <- se_error_ig[,5+idx2-1+extra_id]
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # 1. Build a tidy data‐frame for the betas
  beta_df <- tibble(
    ts        = log(my_list$ts),
    beta1     = est1,
    beta1_lo  = est1 - 2*se_err1,
    beta1_hi  = est1 + 2*se_err1,
    beta2     = est2,
    beta2_lo  = est2 - 2*se_err2,
    beta2_hi  = est2 + 2*se_err2
  ) %>%
    pivot_longer(
      cols      = -ts,
      names_to  = "name",
      values_to = "value"
    ) %>%
    separate(
      name,
      into = c("param","bound"),
      sep  = "_",
      fill = "right"
    ) %>%
    mutate(
      bound = ifelse(is.na(bound), "estimate", bound)
    ) %>%
    pivot_wider(
      names_from = bound,
      values_from = value
    )
  
  # 2. Plot β1 & β2
  beta_plot <- ggplot(beta_df, aes(x = ts, y = estimate, color = param, fill = param)) +
    geom_hline(yintercept = c( 1, -1), linetype="dashed", color="black", alpha=0.6) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha=0.2, color=NA) +
    geom_line(size=1) +
    scale_color_brewer("Parameter", palette="Set1",
                       labels=c(expression(beta[1]),expression(beta[2]))) +
    scale_fill_brewer("Parameter",  palette="Set1",
                      labels=c(expression(beta[1]),expression(beta[2]))) +
    labs(
      x     = expression(log(t)),
      y     = "Estimate ± 2 SE",
      title = expression(beta[1]~" & "~beta[2]~" over Time")
    ) +
    theme_minimal(base_size=36) +
    theme(
      plot.title      = element_text(hjust=0.5),
      legend.position = c(0.8,0.5),
      legend.background = element_rect(fill="white")
    )
  
  
  
  # 1. Tidy the a‐estimates
  a_df <- tibble(
    ts     = log(my_list$ts),
    a11    = aest1,
    a11_lo = aest1 - 2*ase_err1,
    a11_hi = aest1 + 2*ase_err1,
    a31    = aest2,
    a31_lo = aest2 - 2*ase_err2,
    a31_hi = aest2 + 2*ase_err2
  ) %>%
    pivot_longer(
      cols      = -ts,
      names_to  = "name",
      values_to = "value"
    ) %>%
    separate(
      name,
      into = c("param","bound"),
      sep  = "_",
      fill = "right"
    ) %>%
    mutate(
      bound = ifelse(is.na(bound), "estimate", bound)
    ) %>%
    pivot_wider(
      names_from  = bound,
      values_from = value
    )
  
  a_true <- my_list$sim$A[upper.tri(my_list$sim$A,  TRUE)]
  
  # 2. Plot a11 & a31
  a_plot <- ggplot(a_df, aes(x = ts, y = estimate, color = param, fill = param)) +
    geom_hline(yintercept = c(a_true[idx1], a_true[idx2]),
               linetype = "dashed",
               color    = "black",
               alpha    = 0.6) +
    geom_ribbon(aes(ymin = lo, ymax = hi),
                alpha = 0.2, color = NA) +
    geom_line(size = 1) 
  
  return(list(beta_plot = beta_plot, a_plot = a_plot))
}
