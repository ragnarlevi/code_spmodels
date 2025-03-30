### This files contains utillity functions that are needed across all models

############## Mixing Densities ##############
gamma_density <- function(z, phi) dgamma(z, shape = phi, rate = phi)
ig_density <- function(z, phi) actuar::dinvgauss(z, mean = 1, shape = phi^2)
ln_density <- function(z, phi) dlnorm(z, meanlog = -phi^2/2, sdlog = phi)



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

get_grad_psi <- function(additive, grad_mu, agg_effect, nu, exposure_batch, locs_batch){
  if(additive){
    grad_psi <- grad_mu * exposure_batch * agg_effect
  }else{
    grad_psi <- grad_mu * exposure_batch * agg_effect *  nu
  }
  
  sapply(1:nr_regions, function(r) sum(grad_psi[locs_batch == r]))
}

get_grad_a <- function(additive, grad_mu, agg_claims, years,locs, exposure, nu, lambda){
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


get_grad_beta_phi <- function(phi, w1 = NA, w2 = NA, w3 = NA, w4 = NA){
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




