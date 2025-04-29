### This files contains utillity functions that are needed across all models

############## Mixing Densities ##############
gamma_density <- function(z, phi) dgamma(z, shape = phi, rate = phi)
ig_density <- function(z, phi) actuar::dinvgauss(z, mean = 1, shape = phi^2)
ln_density <- function(z, phi) dlnorm(z, meanlog = -phi^2/2, sdlog = phi)

rgamma_density <- function(n, phi) rgamma(n, shape = phi, rate = phi)
rig_density <- function(n, phi) actuar::rinvgauss(n, mean = 1, shape = phi^2)
rln_density <- function(n, phi) rlnorm(n, meanlog = -phi^2/2, sdlog = phi)



####################
# Functions for E-step
###################
# 1. Numerical Integration Version --------------------------------------------
expected_h_integration <- function(y, phi, pi, mu, z_density, h_func, likelihood_func, special_case, special_fun, 
                                   phi_old = phi, pi_old = pi, mu_old = mu) {
  
  if(special_case){
    return(special_fun(z, y, pi, mu))
  }
  
  # print(length(y))
  # print(length(pi))
  # print(length(pi_old))
  # print(length(mu))
  # print(length(mu_old))
  # print(length(phi))
  # print(length(phi_old))
  
  integrand_num <- function(z) {
    h_val <- h_func(z, y, pi, mu)
    lik   <- likelihood_func(z, y, pi_old, mu_old)
    val   <- h_val * lik * z_density(z, phi_old)
    val[!is.finite(val)] <- 0
    return(val)
  }
  integrand_denom <- function(z) {
    lik <- likelihood_func(z, y, pi_old, mu_old)
    val <- lik * z_density(z, phi_old)
    val[!is.finite(val)] <- 0
    return(val)
  }
  
  num   <- integrate(integrand_num, lower = 1e-3, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  denom <- integrate(integrand_denom, lower = 1e-3, upper = Inf, subdivisions = 1000, rel.tol = 1e-8)$value
  if (denom == 0) return(0)
  return(num / denom)
}

# 2. Monte Carlo EM (MCEM) Version --------------------------------------------
expected_h_MCEM <- function(y,  phi, pi, mu, h_func, likelihood_func, rzdensity, special_case, special_fun, S = 1000,
                            phi_old = phi, pi_old = pi, mu_old = mu) {
  
  if(special_case){
    return(special_fun(z, y, pi, mu))
  }
  # Draw S samples from the z prior:
  z_samples <- rzdensity(S, phi_old)
  
  # Compute h(z) and likelihood for each sample
  h_vals   <- sapply(z_samples, function(z) h_func(z, y, pi, mu))
  lik_vals <- sapply(z_samples, function(z) likelihood_func(z, y, pi_old, mu_old))
  
  num   <- sum(h_vals * lik_vals)
  denom <- sum(lik_vals)
  if (denom == 0) return(0)
  return(num / denom)
}

# 3. Gauss–Legendre Quadrature Version -----------------------------------------
expected_h_legendre <- function(y, phi, pi, mu, z_density, h_func, likelihood_func, special_case, special_fun, n_nodes = 50,
                                 phi_old = phi, pi_old = pi, mu_old = mu) {
  
  
  if(special_case){
    return(special_fun(z, y, pi, mu))
  }
  
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required for Gauss-Legendre quadrature. Please install it.")
  }
  # Get Gauss-Legendre nodes and weights on [-1,1]
  quad <- statmod::gauss.quad(n_nodes, kind = "legendre")
  # Transform nodes and weights to [0,1]
  t_nodes <- (quad$nodes + 1) / 2
  t_weights <- quad$weights / 2
  
  # Transformation: z = t/(1-t), dz/dt = 1/(1-t)^2.
  z_vals   <- t_nodes / (1 - t_nodes)
  jacobian <- 1 / (1 - t_nodes)^2
  
  integrand_num_vals   <- numeric(length(z_vals))
  integrand_denom_vals <- numeric(length(z_vals))
  
  for (i in seq_along(z_vals)) {
    z <- z_vals[i]
    h_val <- h_func(z, y, pi, mu)
    lik   <- likelihood_func(z, y, pi_old, mu_old)
    integrand_num_vals[i]   <- h_val * lik * z_density(z, phi_old) * jacobian[i]
    integrand_denom_vals[i] <- lik * z_density(z, phi_old) * jacobian[i]
  }
  
  num_quad   <- sum(t_weights * integrand_num_vals, na.rm = T)
  denom_quad <- sum(t_weights * integrand_denom_vals, na.rm = T)
  
  if (denom_quad == 0) return(0)
  return(num_quad / denom_quad)
}


# --- 4. Gauss–Hermite Quadrature Version ---
expected_h_hermite <- function(y,  phi, pi, mu, z_density, h_func, likelihood_func, special_case, special_fun, n_nodes = 50,
                               phi_old = phi, pi_old = pi, mu_old = mu) {
  
  if(special_case){
    return(special_fun(z, y, pi, mu))
  }
  
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required for Gauss-Hermite quadrature. Please install it.")
  }
  
  # Gauss–Hermite quadrature approximates:
  #   ∫_{-∞}^∞ f(x)e^{-x^2} dx
  # We change variable: z = exp(x), so dz = exp(x) dx.
  # Then the numerator becomes:
  #   I_num = ∫_{-∞}^∞ h(e^x)L(e^x)π(e^x)e^x dx.
  # To bring it into Gauss–Hermite form, define:
  #   f_num(x) = h(e^x)L(e^x)π(e^x)e^x * exp(x^2)
  # Similarly, for the denominator:
  #   f_denom(x) = L(e^x)π(e^x)e^x * exp(x^2)
  
  gh <- statmod::gauss.quad(n_nodes, kind = "hermite")
  nodes <- gh$nodes    # x values over (-∞, ∞)
  weights <- gh$weights
  
  f_num_vals <- sapply(nodes, function(x) {
    z <- exp(x)
    h_val <- h_func(z, y, pi, mu)
    lik <- likelihood_func(z, y, pi_old, mu_old)
    prior <- z_density(z, phi_old)
    return(h_val * lik * prior * exp(x) * exp(x^2))
  })
  
  f_denom_vals <- sapply(nodes, function(x) {
    z <- exp(x)
    lik <- likelihood_func(z, y, pi_old, mu_old)
    prior <- z_density(z, phi_old)
    return(lik * prior * exp(x) * exp(x^2))
  })
  
  num_quad <- sum(weights * f_num_vals)
  denom_quad <- sum(weights * f_denom_vals)
  
  if (denom_quad == 0) {
    
    return(0)
    }
  
  return(num_quad / denom_quad)
}

# Wrapper Function ---------------------------------------------------------
expected_h <- function(y, phi, pi, mu, z_density, h_func, likelihood_func, rzdensity, special_case = FALSE, special_fun = NA,  
                       method = c("integration", "MCEM", "GaussLegendre", "GaussHermite"),
                       S = 1000, n_nodes = 50, phi_old = phi, pi_old = pi, mu_old = mu) {
  
  method <- match.arg(method)
  
  
  if (method == "integration") {
    
    # return(expected_h_integration(y, phi, pi, mu,
    #                               z_density, h_func, likelihood_func,
    #                               special_case, special_fun,
    #                               phi_old = phi_old,
    #                               pi_old  = pi_old,
    #                               mu_old  = mu_old))
    
    return( tryCatch(
      {
        expected_h_integration(y, phi, pi, mu,
                               z_density, h_func, likelihood_func,
                               special_case, special_fun,
                               phi_old = phi_old,
                               pi_old  = pi_old,
                               mu_old  = mu_old)
      },
      error = function(e) {
        print(y)
        print(mu)
        print(pi)
        print(phi)

        warning("`integrate()` failed: ", paste0(y, " ", mu, " ", pi, " ", phi),
                "\n  falling back to Gauss–Legendre quadrature.")

        expected_h_legendre(y, phi, pi, mu,
                            z_density, h_func, likelihood_func,
                            special_case, special_fun,
                            n_nodes = n_nodes,
                            phi_old = phi_old,
                            pi_old  = pi_old,
                            mu_old  = mu_old)
      }
    ) )
  } else if (method == "MCEM") {
    return(expected_h_MCEM(y,  phi, pi, mu, h_func,  likelihood_func, rzdensity, special_case, special_fun, S,
                            phi_old = phi_old, pi_old = pi_old, mu_old = mu_old))
  } else if (method == "GaussLegendre") {
    return(expected_h_legendre(y,  phi, pi, mu, z_density, h_func, likelihood_func, special_case, special_fun, n_nodes,
                               phi_old = phi_old, pi_old = pi_old, mu_old = mu_old))
  }else if (method == "GaussHermite") {
    return(expected_h_hermite(y,  phi, pi, mu, z_density, h_func, likelihood_func, special_case, special_fun, n_nodes,
                              phi_old = phi_old, pi_old = pi_old, mu_old = mu_old))
  }
}

# Optionally, vectorize the wrapper over y and mu if needed:
expected_h <- Vectorize(expected_h, vectorize.args = c("y", "mu", "special_case", "mu_old"))







#########################
# Define general gradients
#########################

get_grad_beta <- function(additive, grad_mu, nu, exposure_batch, X_batch, spatial_effect, locs_batch){
  if(additive){
    return(as.numeric(t(grad_mu * nu * exposure_batch) %*% X_batch))
  }else{
    return(as.numeric(t(grad_mu * nu * exposure_batch * (1+spatial_effect)) %*% X_batch))
  }
}

get_grad_psi <- function(additive, grad_mu, agg_effect, nu, exposure_batch, locs_batch, nr_regions){
  if(additive){
    grad_psi <- grad_mu * exposure_batch * agg_effect
  }else{
    grad_psi <- grad_mu * exposure_batch * agg_effect *  nu
  }
  
  sapply(1:nr_regions, function(r) sum(grad_psi[locs_batch == r]))
}

get_grad_a <- function(additive, grad_mu, agg_claims, years,locs, exposure, nu, lambda, nr_regions){
  
  if(additive){
    g2 <- grad_mu * t(agg_claims[, years]) * exp(log(exposure))
  }else{
    g2 <-  grad_mu * t(agg_claims[, years]) * exp(log(exposure)) * as.numeric(nu)
  }
  
  g2 <- by(g2,locs, FUN=colSums)
  G2 <- matrix(unlist(g2),nrow = nr_regions, ncol = nr_regions, byrow = T)
  g22 <- G2[upper.tri(G2, diag = T)]
  diag(G2) = 0  # make sure we do not double count the diagonal when we add
  g22 <- g22 + t(G2)[upper.tri(G2, diag = T)] - lambda
  
  return(g22)
}



get_grad_pi <- function(grad_pi){
  sum(grad_pi)
}


get_grad_beta_phi <- function(phi, mixing, w1 = NA, w2 = NA, w3 = NA, w4 = NA){
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
# Define gradient update functions
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




########### Other Function ##############

# get default value if a is
`%||%` <- function(a, b) if (!is.null(a)) a else b


# Spatial aggregate function - used for loading models
# Assumes that the first column of agg_claims is claims at time t=0
# The function aggregats spatial effects for each locations.
# It returns a data frame, 

get_spatial_aggregate <- function(locs, w, psi, agg_claims, years, model_type){
  
  
  if(model_type == "ordinary"){
    spatial_effect <- 0
    agg_effect <- 0
  }
  else if(model_type == "learn_psi"){
    spatial_effect <- unname(rowSums(t(agg_claims[, years])*w[locs,]))*psi[locs]
    agg_effect <- unname(rowSums(t(agg_claims[, years])*w[locs,])) 
    
  }else if(model_type == "learn_graph"){
    spatial_effect <- unname(rowSums(t(agg_claims[, years])*w[locs,])) 
    agg_effect <- unname(rowSums(t(agg_claims[, years])*w[locs,])) 
    
  }
  
  return(list(spatial_effect = spatial_effect, agg_effect = agg_effect))
}




# Function to build matrix from vector
get_W_from_array <- function(a, p){
  
  A <- matrix(0, nrow =p, ncol = p)
  A[upper.tri(A, diag = T)] <- a
  A <- A +  t(A)
  diag(A) <- diag(A)*0.5
  
  return(A)
  
}


# Get mean function
get_mu <- function(nu, spatial_effect, exposure, additive){
  
  if(additive){
    mu <- (nu + spatial_effect)*exp(log(exposure))
  }else{
    mu <- (nu * (1 + spatial_effect))*exp(log(exposure))
  }
  
}




