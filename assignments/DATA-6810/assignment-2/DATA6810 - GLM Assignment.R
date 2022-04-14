library(tidyverse)
library(mvtnorm)

# ------------------------------------------------------------------------------
# ------------------------------ Initialize Data -------------------------------
# ------------------------------------------------------------------------------
woodlark <- read_csv("woodlark.csv")
N = 71

# I will be using the 2008 data, as suggested in the assignment description
ggplot(woodlark) +
  geom_point(aes(x = seq(1, N, 1), y = Y08)) +
  theme_bw()

# Input and dependent variables
y <- woodlark$Y08
t <- seq(1, N)
x <- matrix(c(t, t^2), ncol = 2)
X <- matrix(c(rep(1, N), t, t^2), ncol = 3)


# ------------------------------------------------------------------------------
# ---------------------------- Part 1: MLE Estimates ---------------------------
# ------------------------------------------------------------------------------

# Fit GLM model with poisson family
poisson_glm_fit <- glm(
  formula = y ~ x,
  family = poisson(link = "log"),
  data.frame(y = y, t = x)
)

summary(poisson_glm_fit)

# MLE estimates of mean and covariance for beta
beta_hat <- matrix(summary(poisson_glm_fit)$coefficients[ , 1])
V_hat <- summary(poisson_glm_fit)$cov.unscaled

# ------------------------------------------------------------------------------
# ---------------------------- Part 2: MCMC-MH ---------------------------------
# ------------------------------------------------------------------------------

sigma_beta <- diag(c(1^2, 0.25^2, 0.25^2))

log_posterior <- function(beta, X, Y, sigma_beta){
  
  beta <- as.matrix(beta)
  
  lambda <- exp(X %*% beta)
  
  sum(X %*% beta * Y - exp(X %*% beta))  - 
    0.5 * t(beta) %*% solve(sigma_beta) %*% beta
  
}

Y <- matrix(y)

n_mcmc <- 1000

beta_sample <- array(NA, dim = c(3, 1, n_mcmc))
accept_sample <- array(NA, dim = c(n_mcmc))

beta_sample[ , , 1] <- beta_hat
accept_sample[1] <- 1

for(i in 2:n_mcmc){
  
  cat(sprintf("\rIteration %i/%i", i, n_mcmc))
  
  beta_proposed <- c(rmvnorm(1, beta_hat, V_hat))
  
  lp_curr <- log_posterior(beta_sample[ , , i-1], X, Y, sigma_beta)
  lp_prop <- log_posterior(beta_proposed, X, Y, sigma_beta)
  
  q_curr <- dmvnorm(beta_sample[ , , i-1], beta_hat, V_hat, log = TRUE)
  q_prop <- dmvnorm(beta_proposed, beta_hat, V_hat, log = TRUE)
  
  alpha <- min(1, exp(lp_prop - lp_curr + q_curr - q_prop))
  
  u <- runif(1, 0, 1)
  
  if(alpha > u){
    beta_sample[ , , i] <- beta_proposed
    accept_sample[i] <- 1
  } else {
    beta_sample[ , , i] <- beta_sample[ , , i-1]
    accept_sample[i] <- 0
  }
  
}

## Percent accepted samples
length(which(accept_sample == 1))/n_mcmc

# ------------------------------------------------------------------------------
# ------------------------ Simulation Study - Part 6 ---------------------------
# ------------------------------------------------------------------------------

burnin <- 500

plot(beta_sample[1, , burnin:n_mcmc], type = "l")
plot(density(beta_sample[1, , burnin:n_mcmc]))

plot(beta_sample[2, , burnin:n_mcmc], type = "l")
plot(density(beta_sample[2, , burnin:n_mcmc]))


plot(beta_sample[3, , burnin:n_mcmc], type = "l")
plot(density(beta_sample[3, , burnin:n_mcmc]))


beta_sample_tbl <- tibble(
  i = burnin:n_mcmc,
  beta0 = beta_sample[1, 1, burnin:n_mcmc],
  beta1 = beta_sample[2, 1, burnin:n_mcmc],
  beta2 = beta_sample[3, 1, burnin:n_mcmc]
)

ggplot(beta_sample_tbl) +
  geom_point(aes(x = beta0, y = beta1), alpha = 0.1) +
  geom_density_2d(aes(x = beta0, y = beta1)) +
  geom_point(
    data = tibble(beta0 = beta[1,], beta1 = beta[2,]),
    aes(x = beta0, y = beta1), colour = "red"
  ) +
  theme_bw()
