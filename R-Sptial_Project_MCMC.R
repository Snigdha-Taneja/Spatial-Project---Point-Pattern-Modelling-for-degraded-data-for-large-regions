library(spatstat)
library(Matrix)
library(invgamma)
library(RSpectra)
library(fields)
library(truncnorm)
library(SMUT)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
sourceCpp("matrix_inv.cpp")
set.seed(123)


data <- read.csv('aurea_variables.csv')
m <- 450
random_indices <- sample(1:nrow(data), m)
r <- 200
n <- length(data$latitude)

# Select the rows from the dataframe
samp_data <- data[random_indices, ]
rem_data <- data[-random_indices, ]


random_indices2 <- sample(1:nrow(samp_data), r)
samp_data_r <- samp_data[random_indices2, ]
rem_data_r <- samp_data[-random_indices2, ]
m_locations <- samp_data
r_locations <- samp_data_r

n_iterations <- 1000 # Number of MCMC iterations
burn_in <- 250 # Number of burn-in iterations
l = 7 #the number of regressors
r_knot = r


epsilon <- 1e-6

scaled_vars_m = scale(m_locations[,4:9])
scaled_vars_rem = scale(rem_data[,4:9])

X1 <- scaled_vars_m
X2 <- scaled_vars_rem
X1_f <- cbind(rep(1,m),X1)
X2_f <- cbind(rep(1,n-m),X2)


#epsi_inv <- R_Ir %*% R_r_inv %*% R_rI
#epsi <- numeric(2000)
#for (i in 1:2000)
#{
# epsi[i] <- rnorm(1,0, 1- diag(epsi_inv)[i])
#}

# Initialize parameters
beta <- rep(1, 7)  # Initialize beta
sigma2 <- 0.05   # Initialize sigma^2
phi <- 1  # Initialize phi
#w <- rnorm(r_knot)  # Initialize w
lambda <- rep(1, m)  # Initialize lambda

# Define prior parameters
beta0 <- rep(1, 7)  # Prior mean for beta
Sigma0 <- diag(7)  # Prior covariance for beta
a0 <- 1.5  # Shape parameter for inverse gamma prior
b0 <- 2  # Scale parameter for inverse gamma prior
phi0 <- 0.01 # Prior minimum for phi
phi1 <- 1  # Prior maximum for phi

# Initialize φ and other parameters
phi_current <- 0.01
phi_proposal_sd <- 0.05
likeli = numeric(length = n_iterations)


prior_phi <- function(phi, phi0, phi1) {
  if (phi >= phi0 && phi <= phi1) {
    return(1 / (phi1 - phi0))
  } else {
    return(0)
  }
}

likelihood <- function(lambda, phi, L1, R_rm, M1, X1, beta, sigma2, a) {
  e <- as.matrix(log(lambda)) - eigenMapMatMult(as.matrix(X1_f),as.matrix(t(beta)))
  eig <- eigs_sym(a,k = m)$values
  S_phi <- -sum(log(abs(eig)))/2
  #S_phi <- -log(abs((det(solve(diag(1,r)) + R_rm %*% solve(M1) %*% L1) * (det(M1)) * det(diag(1,r)))) / 2
  S_phi <- S_phi - eigenMapMatMult(eigenMapMatMult(t(e),mult), e)/(2 * sigma2)
  S_phi <- as.numeric(S_phi)
  return(S_phi)
}

# Create vectors to store posterior samples
beta_samples <- matrix(NA, nrow = n_iterations, ncol = 7)
sigma2_samples <- rep(NA, n_iterations)
phi_samples <- rep(NA, n_iterations)
w_samples <- matrix(NA, nrow = n_iterations, ncol = r)  # Adjust dimensions
phi_samples <- numeric(n_iterations+1)
phi_samples[1] <- phi_current
DIC <- rep(NA, n_iterations)

# Number of samples from the posterior distribution
n_samples <- nrow(beta_samples)

# Initialize storage for predicted values
predicted_lambda <- matrix(NA, nrow = n_samples, ncol = nrow(X2_f))



lambda_samples <- matrix(NA, nrow = n_iterations+1, ncol = m)
lambda_samples[1, ] <- lambda[1:m]  # Initialize with the initial values of λ

# Define the proposal distribution parameters (you can adjust them)
proposal_mean <- 0.15  
proposal_sd <- 2


# MCMC loop
for (iteration in 1:n_iterations) {
  print("iteration")
  print(iteration)
  # step (a)
  R_mr <- exp(-rdist(m_locations[,2:3],r_locations[,2:3]))/phi_current
  #since R_r is not singular
  
  R_r_reg <- exp(-rdist(r_locations[,2:3]))/phi + epsilon * diag(nrow(r_locations))
  R_r_inv <- solve(R_r_reg)
  
  L1 <- eigenMapMatMult(R_mr , R_r_inv)
  
  R_Imr <- exp(-rdist(rem_data[,2:3],r_locations[,2:3]))/phi_current
  L2 <- eigenMapMatMult(R_Imr , R_r_inv)
  
  R_rm <- exp(-rdist(r_locations[,2:3],m_locations[,2:3]))/phi_current
  
  L_rm = eigenMapMatMult(L1,R_rm)
  
  M1 <- diag(1,m) - diag(diag(L_rm))
  
  R_rIm <- exp(-rdist(r_locations[,2:3],rem_data[,2:3]))/phi_current
  M2 <- diag(1,n-m) - diag(diag(eigenMapMatMult(L2,R_rIm)))
  
  w <- rmvnorm(1, mu = rep(0,r), Sigma = R_r_reg)
  R_Ir <- exp(-rdist(data[,2:3],r_locations[,2:3]))/phi_current
  #w_curl <- R_Ir %*% R_r_inv %*% t(w)
  R_rI <- exp(-rdist(r_locations[,2:3],data[,2:3]))/phi_current
  
  mult = inverseMatrix(L_rm + M1)
  
  # Step (b): Sample beta
  Sigma_beta_inverse <- solve(Sigma0) + (1 / (sigma2)) * t(X1_f) %*% mult %*% as.matrix(X1_f)
  mu_beta <- inverseMatrix(Sigma_beta_inverse) %*% (solve(Sigma0) %*% beta0 + (1 / (sigma2)) * t(X1_f) %*% mult %*% as.matrix(log(lambda)))
  near_positive_definite_Sigma <- nearPD(solve(Sigma_beta_inverse))$mat
  beta <- rmvt(1, df = 1, mu = mu_beta, Sigma = near_positive_definite_Sigma)
  #beta <- rmvt(1, df = 1, mu = mu_beta, Sigma = solve(Sigma_beta_inverse)+0.5*diag(7))
  
  # Step (c): Sample sigma^2
  a_sigma2 <- a0 + (m / 2)
  e <- log(lambda) - as.matrix(X1_f) %*% as.matrix(t(beta))
  #b_sigma2 <- abs(b0+as.numeric((t(e) %*% mult %*% as.matrix(e)) / 2))
  b_sigma2 <- b0+as.numeric((t(e) %*% mult %*% as.matrix(e))) / 2
  print(b_sigma2)
  sigma2 <- rinvgamma(n = 1, shape = a_sigma2, rate = abs(b_sigma2))
  print(sigma2)
  
  # Step (d): Sample w
  Sigma_w_inverse <- (R_r_inv) + (1 / (sigma2)) * (eigenMapMatMult(eigenMapMatMult(t(L1),inverseMatrix(M1)),L1))
  #Sigma_w_reg  = Sigma_w_inverse + 1e-3 * diag(nrow(Sigma_w_inverse))
  #chol_decomp <- chol(Sigma_w_reg, pivot = TRUE,tol = 1e-10)
  mu_w <- inverseMatrix(Sigma_w_inverse) %*% ((1 / (sigma2 )) * eigenMapMatMult(eigenMapMatMult(t(L1),inverseMatrix(M1)),as.matrix((log(lambda) - eigenMapMatMult(as.matrix(X1_f),as.matrix(t(beta)))))))
  w_temp = rnorm(nrow(Sigma_w_inverse))
  w = mu_w + eigenMapMatMult(inverseMatrix(Sigma_w_inverse),as.matrix(w_temp))
  #w <- rmvnorm(1, mu = mu_w, Sigma = ginv(Sigma_w_inverse))
  
  # Store posterior samples
  if (iteration > burn_in) {
    beta_samples[iteration - burn_in, ] <- beta
    sigma2_samples[iteration - burn_in] <- sigma2
    w_samples[iteration - burn_in, ] <- w
  }
  
  # Step (e): Sample phi
  # You mentioned using a Metropolis-Hastings update for phi, which is a more complex step not covered here.
  #sampling metropolis
  
  phi_proposal <- rnorm(1, mean = phi_current, sd = phi_proposal_sd)
  
  R_mr_p <- exp(-rdist(m_locations[,2:3],r_locations[,2:3]))/phi_proposal
  #since R_r is not singular
  
  R_r_reg_p <- exp(-rdist(r_locations[,2:3]))/phi_proposal + epsilon * diag(nrow(r_locations))
  R_r_inv_p <- solve(R_r_reg_p)
  
  L1_p <- eigenMapMatMult(R_mr_p ,R_r_inv_p)
  
  R_rm_p <- exp(-rdist(r_locations[,2:3],m_locations[,2:3]))/phi_proposal
  
  L_rm_p = eigenMapMatMult(L1_p,R_rm_p)
  
  M1_p <- diag(1,m) - diag(diag(L_rm_p))
  
  mult_p = inverseMatrix(L_rm_p + M1_p)
  a <- inverseMatrix(mult)
  b <- inverseMatrix(mult_p)
  
  # Calculate the prior and likelihood ratios
  prior_ratio <- prior_phi(phi_proposal, phi0, phi1) / prior_phi(phi_current, phi0, phi1)
  lik1 <- likelihood(lambda, phi_current, L1, R_rm, M1, X1_f, beta, sigma2, a) 
  lik2 <- likelihood(lambda, phi_proposal, L1_p, R_rm_p, M1_p, X1_f, beta, sigma2, b)
  likelihood_ratio <- lik1/lik2
  
  # Calculate the acceptance probability
  accept_prob <- min(1, likelihood_ratio * prior_ratio)
  if (is.nan(accept_prob)) {
    accept_prob <- 1
  }
  if (runif(1) < accept_prob) {
    phi_current <- phi_proposal
  }
  phi_samples[iteration + 1] <- phi_current
  
  #Step (f)
  
  y = rpois(m,lambda_samples[iteration,])
  for (i in 1:m) {
    print("lambda")
    print(i)
    # Propose a new λ value
    lambda_proposal <- rtruncnorm(1,a = 0, b = Inf, mean = lambda[i], sd = proposal_sd)
    
    # Calculate the Log-Normal likelihood for the proposed value
    log_likelihood_candidate <- dlnorm(lambda_proposal, meanlog = as.numeric(as.matrix(t(X1_f[i, ])) %*% as.matrix(t(beta))+ (L1 %*% as.matrix(w))[i]), sdlog = sqrt(sigma2*abs(M1[i,i] )))
    log_likelihood_current <- dlnorm(lambda[i], meanlog = as.numeric(as.matrix(t(X1_f[i, ])) %*% as.matrix(t(beta)) + (L1 %*% as.matrix(w))[i]), sdlog = sqrt(sigma2*abs(M1[i,i] )))
    
    # Calculate the Poisson likelihood (replace with your actual Poisson likelihood)
    log_poisson_likelihood_candidate <- log(dpois(y[i], lambda_proposal))
    log_poisson_likelihood_current <- log(dpois(y[i], lambda[i]))
    
    # Calculate the total log-likelihood
    total_log_likelihood_candidate <- log_likelihood_candidate + log_poisson_likelihood_candidate
    total_log_likelihood_current <- log_likelihood_current + log_poisson_likelihood_current
    
    # Calculate the acceptance probability
    acceptance_ratio_l <- exp(total_log_likelihood_candidate - total_log_likelihood_current)
    
    # Accept or reject the proposal
    if (is.nan(acceptance_ratio_l)) {
      acceptance_ratio_l <- 1
    }
    if (runif(1) < acceptance_ratio_l) 
    {
      lambda[i] <- lambda_proposal
    }
    lambda_samples[iteration + 1, i] <- lambda[i]
  }
  
  for (i in 1:nrow(X2_f)) {
    print('prediction')
    print(i)
    prediction_i <- as.numeric(t(as.matrix(X2_f[i, ])) %*% t(as.matrix(beta)) + t(L2%*%w)[i]) + sqrt(sigma2*abs(M2[i,i] )) * rnorm(1,mean = 0 , sd = 1)
    # Store the predicted value in the matrix
    predicted_lambda[iteration, i] <- prediction_i
  }
  likeli[i] = lik1
}


#write.csv(beta_samples,'beta_samples_aurea_1.csv')
#write.csv(lambda_samples,'lambda_samples_aurea_1.csv')
#write.csv(sigma2_samples,'sigma2_samples_aurea_1.csv')
#write.csv(predicted_lambda,'predicted_lambda_aurea_1.csv')
#write.csv(w_samples,'w_samples_aurea_1.csv')
#write.csv(phi_samples,'phi_samples_aurea_1.csv')

#aic = numeric(length=4)
loglik = sum(likeli)
#k = m+r
aic[3] = -2*loglik +2*k


View(beta_samples)
