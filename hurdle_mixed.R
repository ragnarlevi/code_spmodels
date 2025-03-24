# Set seed for reproducibility
set.seed(123)
source("utils.R")
source("hurdle.R")
source("simulate_data.R")

#########################
# 1. Define Expected Gradient Functions
#########################



expected_hurdle_h <- function(y, beta, phi, pi, mu, z_density) {

  # For y == 0 the derivative (w.r.t. beta) is zero (due to the indicator 1(y > 0))
  if (y <= 0) {
    return(0)
  }
  
  # Define h(z): derivative with respect to mu
  h <- function(z) {
    y/mu - z - (z * exp(-z * mu)) / (1 - exp(-z * mu))
  }
  
  # f(y|z) for y > 0 under the hurdle Poisson model:
  #   f(y|z) = (1 - pi) * dpois(y, lambda = z * mu) / (1 - exp(-z * mu))
  f_y_given_z <- function(z) {
    (1 - pi) * dpois(y, lambda = z * mu) / (1 - exp(-z * mu))
  }
  
  
  # Numerator: integrate h(z)*f(y|z)*f(z) over z
  numerator_integrand <- function(z) {
    h(z) * f_y_given_z(z) * z_density(z)
  }
  
  # Denominator: integrate f(y|z)*f(z) over z
  denominator_integrand <- function(z) {
    f_y_given_z(z) * z_density(z)
  }
  
  num_result <- integrate(numerator_integrand, lower = 0, subdivisions = 1000,  upper = Inf, rel.tol = 1e-8)
  den_result <- integrate(denominator_integrand, lower = 0, subdivisions = 1000, upper = Inf, rel.tol = 1e-8)
  
  if (den_result$value == 0) {
    warning("Denominator evaluated to zero; check parameter values.")
    return(0)
  }
  
  grad_mu <- num_result$value / den_result$value  # Expected derivative w.r.t. mu
  return(grad_mu)
}
expected_hurdle_h <- Vectorize(expected_hurdle_h, vectorize.args = c("y", "mu"))



expected_hurdle_w <- function(y, phi, pi, mu, z_density, w_fun) {
  
  if (y == 0) {
    integrand_num <- function(z) {
      f_y_given_z <- pi 
      f_y_given_z * w_fun(z) * z_density(z)
    }
    integrand_denom <- function(z) {
      f_y_given_z <- pi 
      f_y_given_z * z_density(z)
    }
  } else {
    integrand_num <- function(z) {
      f_y_given_z <- (1 - pi) * dpois(y, z * mu) / (1-dpois(0, z * mu))
      f_y_given_z * w_fun(z) * z_density(z)
    }
    integrand_denom <- function(z) {
      f_y_given_z <- (1 - pi) * dpois(y, z * mu) / (1-dpois(0, z * mu))
      f_y_given_z * z_density(z)
    }
  }
  
  num <- integrate(integrand_num, lower = 1e-6, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  denom <- integrate(integrand_denom, lower = 1e-6, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  if (denom == 0) return(0)
  return(num / denom)
}
expected_hurdle_w <- Vectorize(expected_hurdle_w, vectorize.args = c("y", "mu"))




#########################
# 4. Test
#########################
source("simulate_data.R")
data_sim <- simulate_hurdle_claims(100, 50, "graph", TRUE, mixing = "gamma")

# Extract variables from simulation
claims <- data_sim$claims
X <- data_sim$X
years <- data_sim$years
locs <- data_sim$locs
agg_claims <- data_sim$agg_claims
A <- data_sim$A
exposure <- data_sim$exposure
model_type <- "learn_graph"
additive <- TRUE
mixing <- "gamma"

if(mixing == "gamma"){
  z_density <- gamma_density  
}else if(mixing == "ln"){
  z_density <- ln_density  
}else if(mixing == "ig"){
  z_density <- ig_density  
}



# Adam parameters for beta and pi
n_iter <- 5000
d <- length(data_sim$beta1)  # dimension of beta
beta1_mom <- 0.9
beta2_mom <- 0.999
epsilon <- 1e-8
momentum <- 0.07

# Get number of regions
nr_regions <- length(unique(data_sim$locs))
nr_edges <- nr_regions*(nr_regions+1)/2
beta_phi <- 1
phi <- exp(beta_phi)


lambda <- 0


# Initialize parameters
out_hurdle <- hurdle(data_sim$claims, data_sim$X, data_sim$locs, data_sim$years,  data_sim$agg_claims, 
               data_sim$A, additive, "learn_psi", lambda = 0, exposure = data_sim$exposure, max_itr = 300)


beta_est <- out_hurdle$beta1    # beta (d-dimensional)
pi_est <- out_hurdle$prop
psi_est <- out_hurdle$psi     # psi (vector of length nr_regions)
a_est <- out_hurdle$a

# Initialize controls
optimizer_beta <- "adam"
optimizer_psi <- "adam"
optimizer_a <- "adam"
optimizer_pi <- "adam"
optimizer_beta_phi <- "adam"
sgd <- FALSE
batch_size <- 500



controls <- list()
controls$beta <- get_controls(optimizer_beta, 0.1, d, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$psi <- get_controls(optimizer_psi, 0.01, nr_regions, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$pi <- get_controls(optimizer_pi, 0.001, 1, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$beta_phi <- get_controls(optimizer_beta_phi, 0.01, 1, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$a <- get_controls(optimizer_a, 0.01, nr_edges, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)







N <- length(claims)
# (Assuming get_spatial_aggregate() and get_mu() are defined in utils.R)
cat("\nStarting updates for beta, pi (Adam) and psi (Adagrad):\n")
for (iter in 1:n_iter) {
  
  if(model_type == "learn_graph"){
    A <- get_W_from_array(a_est, nr_regions)
    psi <- NA
  }else if(model_type == "learn_psi"){
    a_est <- 1
  }
  
  
  # If using SGD, sample a mini-batch
  if (sgd) {
    idx <- c()
    for(i in 1:nr_regions){
      idx <- c(idx, sample(which(locs == i), batch_size, replace = FALSE))
    }
    claims_batch   <- claims[idx]
    X_batch        <- X[idx, , drop = FALSE]
    exposure_batch <- exposure[idx]
    locs_batch     <- locs[idx]
    years_batch <- years[idx]
  } else {
    idx <- 1:N
    claims_batch   <- claims
    X_batch        <- X
    exposure_batch <- exposure
    locs_batch     <- locs
    years_batch <- years
  }
  
  
  phi <- exp(beta_phi)
  
  # Compute spatial aggregate effects and baseline nu
  se <- get_spatial_aggregate(locs_batch, A, psi_est, agg_claims, years_batch, model_type)
  nu <- as.numeric(exp(X_batch %*% beta_est))  # baseline mu from covariates
  
  # Get full mu by combining covariate effects (nu), spatial effects, and exposure
  mu <- get_mu(nu, se$spatial_effect, exposure_batch, additive)
  
  # E-steps
  # For the mean parameters: Compute expected derivative with respect to mu
  grad_mu <- expected_hurdle_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density)
  
  # For the mixing parameter:
  if(mixing == "gamma"){
    w1 <- expected_hurdle_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = identity)
    w2 <- expected_hurdle_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = log)
  }else if(mixing == "ig"){
    w1 <- expected_hurdle_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = identity)
    w3 <- expected_hurdle_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = function(x) 1/x)
  } else if(mixing == "ln"){
    w4 <- expected_hurdle_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = function(x) log(x)^2)
  }
  
  
  # Compute gradient for beta via chain rule: sum_i (grad_mu_i * nu_i * x_i)
  grad_beta_total <- get_grad_beta(additive)
  
  # Compute psi/a gradient
  if(model_type == "learn_graph"){
    grad_a_total <- get_grad_a(additive)
  }else if (model_type == "learn_psi"){
    grad_psi_total <- get_grad_psi(additive)
  }

  # phi grad
  grad_beta_phi_total <- get_grad_beta_phi()
  
  # scale gradient for sgd
  if(sgd){
    grad_beta_total <- grad_beta_total* (N / batch_size)  
    if(model_type == "learn_graph"){
      grad_a_total <- grad_a_total* (N / batch_size)
    }else if (model_type == "learn_psi"){
      grad_psi_total <- grad_psi_total* (N / batch_size)
    }
    grad_beta_phi_total <- grad_beta_phi_total* (N / batch_size)
  }
  
  
  # update beta
  beta_out <- update_gradient(optimizer_beta, beta_est, grad_beta_total, controls$beta)
  beta_est <- beta_out$param
  controls$beta <- beta_out$controls
  
  # update psi/a
  if(model_type == "learn_graph"){
    a_out <- update_gradient(optimizer_a, a_est, grad_a_total, controls$a)
    a_est <- a_out$param
    controls$a <- a_out$controls
    a_est <- pmax(a_est, 1e-6)
  }else if(model_type == "learn_psi"){
    psi_out <- update_gradient(optimizer_psi, psi_est, grad_psi_total, controls$psi)
    psi_est <- psi_out$param
    controls$psi <- psi_out$controls
    psi_est <- pmax(psi_est, 1e-6)
  }


  # update phi
  beta_phi_out <- update_gradient(optimizer_beta_phi, beta_phi, grad_beta_phi_total, controls$beta_phi)
  beta_phi <- beta_phi_out$param
  controls$beta_phi <- beta_phi_out$controls
  
  
  
  if(model_type == "learn_graph"){
    cat(sprintf("Iteration %d: beta = [%s], a = [%s], beta_phi = %f\n,",
                iter,
                paste(round(beta_est, 4), collapse = ", "),
                paste(round(a_est, 4), collapse = ", "),
                beta_phi))
  }else if (model_type == "learn_psi"){
    cat(sprintf("Iteration %d: beta = [%s], psi = [%s], beta_phi = %f\n,",
                iter,
                paste(round(beta_est, 4), collapse = ", "),
                paste(round(psi_est, 4), collapse = ", "),
                beta_phi))
  }
  
}


save.image("hurdle_test.R")

