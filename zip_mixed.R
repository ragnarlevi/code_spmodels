# Set seed for reproducibility
set.seed(123)
source("utils.R")
source("zip.R")
source("simulate_data.R")

#########################
# 1. Define Expected Gradient Functions
#########################



expected_h <- function(y, beta, phi, pi, mu, z_density, do_integral = TRUE) {
  
  # THE DO_INTEGRAL IS JUST TO TEST DERIVATIVES
  
  if(do_integral){
    if (y == 0) {
      integrand_num <- function(z) {
        den <- pi * exp(z * mu) + (1 - pi)
        h_val <- ifelse(is.finite(den), -((1 - pi) * z) / den, 0)
        lik <- pi + (1 - pi) * exp(-z * mu)
        val <- h_val * lik * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
      integrand_denom <- function(z) {
        lik <- pi + (1 - pi) * exp(-z * mu)
        val <- lik * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
    } else {
      integrand_num <- function(z) {
        h_val <- (y / mu) - z
        lik <- (1 - pi) * ((z * mu)^y * exp(-z * mu)) / gamma(y + 1)
        val <- h_val * lik * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
      integrand_denom <- function(z) {
        lik <- (1 - pi) * ((z * mu)^y * exp(-z * mu)) / gamma(y + 1)
        val <- lik * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
    }
    
    num <- integrate(integrand_num, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
    denom <- integrate(integrand_denom, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
    if (denom == 0) return(0)
    return(num / denom)
  }else{
    if (y == 0) {

        den <- pi * exp(1 * mu) + (1 - pi)
        h_val <- ifelse(is.finite(den), -((1 - pi) * 1) / den, 0)
        val <- h_val
        val[!is.finite(val)] <- 0
        return(val)
    } else {
        h_val <- (y / mu) - 1
        val <- h_val 
        val[!is.finite(val)] <- 0
        return(val)
      

    }
    
  }


}
expected_h <- Vectorize(expected_h, vectorize.args = c("y", "mu"))

expected_pi_grad <- function(y, beta, phi, pi, mu, z_density, do_integral = TRUE) {

  # THE DO_INTEGRAL IS JUST TO TEST DERIVATIVES
  if(do_integral){
    if (y == 0) {
      integrand_num <- function(z) {
        num <- (1 - exp(-z * mu))
        denom <- pi + (1 - pi) * exp(-z * mu)
        val <- (num / denom) * (pi + (1 - pi) * exp(-z * mu)) * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
      integrand_denom <- function(z) {
        lik <- pi + (1 - pi) * exp(-z * mu)
        val <- lik * z_density(z, phi)
        val[!is.finite(val)] <- 0
        return(val)
      }
      num <- integrate(integrand_num, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
      denom <- integrate(integrand_denom, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
      if (denom == 0) return(0)
      return(num / denom)
    } else {
      return(-1 / (1 - pi))
    } 
  }else{
    if (y == 0) {
        num <- (1 - exp(-1 * mu))
        denom <- pi + (1 - pi) * exp(-1 * mu)
        val <- (num / denom) 
        val[!is.finite(val)] <- 0
      return(val)
    } else {
      return(-1 / (1 - pi))
    } 
    
  }


}
expected_pi_grad <- Vectorize(expected_pi_grad, vectorize.args = c("y", "mu"))



expected_w <- function(y, phi, pi, mu, z_density, w_fun) {

  if (y == 0) {
    integrand_num <- function(z) {
      f_y_given_z <- pi + (1 - pi) * exp(-z * mu)
      f_y_given_z * w_fun(z) * z_density(z, phi)
    }
    integrand_denom <- function(z) {
      f_y_given_z <- pi + (1 - pi) * exp(-z * mu)
      f_y_given_z * z_density(z, phi)
    }
  } else {
    integrand_num <- function(z) {
      f_y_given_z <- (1 - pi) * ((z * mu)^y * exp(-z * mu)) / factorial(y)
      f_y_given_z * w_fun(z) * z_density(z, phi)
    }
    integrand_denom <- function(z) {
      f_y_given_z <- (1 - pi) * ((z * mu)^y * exp(-z * mu)) / factorial(y)
      f_y_given_z * z_density(z, phi)
    }
  }
  
  num <- integrate(integrand_num, lower = 1e-6, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  denom <- integrate(integrand_denom, lower = 1e-6, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  if (denom == 0) return(0)
  return(num / denom)
}
expected_w <- Vectorize(expected_w, vectorize.args = c("y", "mu"))




get_grad_beta <- function(additive){
  if(additive){
    return(as.numeric(t(grad_mu * nu * exposure_batch) %*% X_batch))
  }else{
    return(as.numeric(t(grad_mu * nu * exposure_batch * (1+se$spatial_effect)) %*% X_batch))
  }
}

get_grad_psi <- function(additive){
  if(additive){
    grad_psi <- grad_mu * exposure_batch * se$agg_effect
  }else{
    grad_psi <- grad_mu * exposure_batch * se$agg_effect *  nu
  }
  
  sapply(1:nr_regions, function(r) sum(grad_psi[locs_batch == r]))
}

get_grad_a <- function(additive){
  if(additive){
    g2 <- grad_mu * t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <-  grad_mu * t(agg_claims[, years]) * exp(log(exposure)) * as.numeric(nu)
  }
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nrow(A), ncol = ncol(A), byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) = 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  return(g22)
}



get_grad_pi <- function(){
  sum(grad_pi)
}


get_grad_beta_phi <- function(){
  if(mixing == "gamma"){
    grad <- sum( phi * ( log(phi) + 1 - digamma(phi) + w2 - w1 ) )
  }else if(mixing == "ig"){
    grad <- sum(1 + 2*phi^2 - phi^2*w3 - w1*phi^2)
  }else if(mixing == "ln"){
    grad <- sum(-1 + w4/phi^2 - phi^2/4)
  }
  return(grad)
}

#########################
# 2.1. Define gradient update function
#########################

get_controls <- function(optimizer, lr, len, ...){
  
  input <- list(...)
  
  if(optimizer == "gd"){
    return(list(lr = lr))
  }else if(optimizer == "adam") {
    m <- rep(0, len)
    v <- rep(0, len)
    return(list(m = m, v = v, lr = lr, beta1 = input$beta1, beta2 = input$beta2, epsilon = input$epsilon, iter = input$iter))
  } else if(optimizer == "adagrad") {
    G <- rep(0, len)
    return(list(G=G, lr = lr , epsilon = input$epsilon))
  } else if(optimizer == "nesterov") {
    v <- rep(0, len)
    return(list(lr = lr, v = v, momentum = input$momentum))
  }
  
  
}

update_gd <- function(param, grad, controls){
  list(param = param + controls$lr*grad, controls = controls)
}

# Adam update function (for scalar or vector parameters)
update_adam <- function(param, grad, controls) {
  controls$m <- controls$beta1 * controls$m + (1 - controls$beta1) * grad
  controls$v <- controls$beta2 * controls$v + (1 - controls$beta2) * (grad^2)
  m_hat <- controls$m / (1 - controls$beta1^controls$iter)
  v_hat <- controls$v / (1 - controls$beta2^controls$iter)
  
  param_new <- param + controls$lr * m_hat / (sqrt(v_hat) + controls$epsilon)
  
  controls$iter <- controls$iter+1
  list(param = param_new, controls=controls)
}

# Adagrad update function
update_adagrad <- function(param, grad, controls) {
  controls$G <- controls$G + grad^2
  param_new <- param + controls$lr * grad / (sqrt(controls$G) + controls$epsilon)
  list(param = param_new, controls=controls)
}

# Nesterov momentum update function
# Here we assume that "grad" is computed at the lookahead position.
update_nesterov <- function(param, grad, controls) {
  # v: velocity vector
  controls$v <- controls$momentum * controls$v + controls$lr * grad
  param_new <- param + controls$v
  list(param = param_new, controls=controls)
}

update_gradient <- function(optimizer, param, grad, controls){
  if(optimizer == "gd"){
    out = update_gd(param, grad, controls)
  }else if(optimizer == "adam") {
    out = update_adam(param, grad, controls)
  } else if(optimizer == "adagrad") {
    out = update_adagrad(param, grad, controls)
  } else if(optimizer == "nesterov") {
    out = update_nesterov(param, grad, controls)
  }
  
  return(out)
}



#########################
# 4. Test
#########################
data_sim <- simulate_claims(100, 20, "psi", TRUE, mixing = "none", model_type = "zip")

# Extract variables from simulation
claims <- data_sim$claims
X <- data_sim$X
years <- data_sim$years
locs <- data_sim$locs
agg_claims <- data_sim$agg_claims
A <- data_sim$A
exposure <- data_sim$exposure
model_type <- "learn_psi"
additive <- TRUE
mixing <- "ig"

if(mixing == "gamma"){
  z_density <- gamma_density  
}else if(mixing == "ln"){
  z_density <- ln_density  
}else if(mixing == "ig"){
  z_density <- ig_density  
}



# Adam parameters for beta and pi
learning_rate_beta <- 0.05
learning_rate_pi <- 0.01
n_iter <- 500
d <- length(data_sim$beta1)  # dimension of beta
beta1_mom <- 0.9
beta2_mom <- 0.999
epsilon <- 1e-8
momentum <- 0.07

# Get number of regions
nr_regions <- length(unique(data_sim$locs))
nr_edges <- nr_regions*(nr_regions+1)/2
beta_phi <- 0
phi <- exp(beta_phi)


lambda <- 0


# Initialize parameters
out_zip <- zip(data_sim$claims, data_sim$X, data_sim$locs, data_sim$years,  data_sim$agg_claims, 
               data_sim$A, additive, "learn_psi", lambda = 0, exposure = data_sim$exposure, max_itr = 300)


beta_est <- c(0,0,0) #out_zip$beta1    # beta (d-dimensional)
pi_est <- 0.5#out_zip$prop           # pi (scalar)
psi_est <- rep(0, nr_regions)# out_zip$psi     # psi (vector of length nr_regions)
a_est <- rep(0, nr_regions*(nr_regions + 1)/2)# out_zip$a

# Initialize controls
optimizer_beta <- "gd"
optimizer_psi <- "gd"
optimizer_a <- "gd"
optimizer_pi <- "gd"
optimizer_beta_phi <- "gd"
sgd <- FALSE
batch_size <- 100



controls <- list()
controls$beta <- get_controls(optimizer_beta, 0.001, d, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$psi <- get_controls(optimizer_psi, 0.1, nr_regions, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$pi <- get_controls(optimizer_pi, 0.0001, 1, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$beta_phi <- get_controls(optimizer_beta_phi, 0.001, 1, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)
controls$a <- get_controls(optimizer_a, 0.001, nr_edges, beta1 = beta1_mom, beta2 = beta2_mom, epsilon = epsilon, iter = 1)







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
  # Compute expected derivative with respect to mu
  grad_mu <- expected_h(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, TRUE)
  grad_pi <- expected_pi_grad(y = claims_batch, beta = beta_est, phi = phi, pi = pi_est, mu = mu, z_density = z_density, TRUE)

  
  if(mixing == "gamma"){
    w1 <- expected_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = identity)
    w2 <- expected_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = log)
  }else if(mixing == "ig"){
    w1 <- expected_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = identity)
    w3 <- expected_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = function(x) 1/x)
  } else if(mixing == "ln"){
    w4 <- expected_w(y = claims_batch, phi = phi, pi = pi_est, mu = mu, z_density = z_density, w_fun = function(x) log(x)^2)
  }
  

  
  
  # Compute gradient for beta via chain rule: sum_i (grad_mu_i * nu_i * x_i)
  grad_beta_total <- get_grad_beta(additive)
  
  # Compute psi/a gradient
  if(model_type == "learn_graph"){
    grad_a_total <- get_grad_a(additive)
  }else if (model_type == "learn_psi"){
    grad_psi_total <- get_grad_psi(additive)
  }
  
  # Compute gradient for pi (scalar)
  grad_pi_total <- get_grad_pi()
  
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
    grad_pi_total <- grad_pi_total* (N / batch_size)  
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

  
  # update pi
  pi_out <- update_gradient(optimizer_pi, pi_est, grad_pi_total, controls$pi)
  pi_est <- pi_out$param
  controls$pi <- pi_out$controls
  pi_est <- pmin(pmax(pi_est, 1e-6), 1-1e-6)
  
  
  # update phi
  beta_phi_out <- update_gradient(optimizer_beta_phi, beta_phi, grad_beta_phi_total, controls$beta_phi)
  beta_phi <- beta_phi_out$param
  controls$beta_phi <- beta_phi_out$controls
  
  print(sum(w3==0))

  
  if(model_type == "learn_graph"){
    cat(sprintf("Iteration %d: beta = [%s], pi = %f, a = [%s], beta_phi = %f\n,",
                iter,
                paste(round(beta_est, 4), collapse = ", "),
                pi_est,
                paste(round(a_est, 4), collapse = ", "),
                beta_phi))
  }else if (model_type == "learn_psi"){
    cat(sprintf("Iteration %d: beta = [%s], pi = %f, psi = [%s], beta_phi = %f\n,",
                iter,
                paste(round(beta_est, 4), collapse = ", "),
                pi_est,
                paste(round(psi_est, 4), collapse = ", "),
                beta_phi))
  }

}


