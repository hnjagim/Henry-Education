#Load the necessary packages
library(rugarch)
library(Rcpp)
library(RcppArmadillo)
library(BayesianTools)
library(quantmod)
library(bayesGARCH)
library(psych)
library(fBasics)


#Define the model
spec <- bayesGARCH::rhgarch(
  formula = ~ garch(1, 1) + rhgarch(1, 1),
  data = returns,
  start = c(0.01, 0.05, 0.9, 0.1),
  fixed.pars = list(v = 1, omega = 0.01, gamma = 0.1),
  constraints = list(beta1 = c(0, 1), beta2 = c(0, 1))
)




#Prepare the data
symbols <- c("^N225", "^FTSE", "^DJI")
getSymbols(symbols, src = "yahoo", from = "2021-01-01", to = "2023-04-12")
FTSE100 <- Ad(get("FTSE"))
DJI100 <- Ad(get("DJI"))

NIKKEI225_Close <- Ad(N225)
# Reset the index of the N225 data frame
N225_df <- data.frame(Date = index(N225), N225.Low = N225$N225.Adjusted)
row.names(N225_df) <- NULL
# Extract the Date and N225.Low columns into separate data frames
dates <- N225_df$Date
prices <- N225_df$N225.Adjusted
returns <- diff(log(prices))

#===============================================
#Summary Statistics
#===============================================
# Calculate daily returns
returns <- diff(prices)/prices[-length(prices)]
# Calculate log returns
log_returns <- diff(log(prices))
# Summary statistics for both types of returns
describe(cbind(returns, log_returns))
basicStats(returns)
#===============================================

#Compute realized volatility using one of the methods 
#provided in the package Realized, 
#e.g., using the bi-power variation estimator:
rv <- realized_bi_power(returns, bipower_coef = 1.5, robust = FALSE)

#Prepare the data for the Bayesian model
y <- returns[-1]
x <- cbind(1, rv)

#prior distributions for the model parameters
prior <- prior_list(mu = NULL,
                    alpha = c(1,1),
                    beta = c(1,1),
                    omega = NULL,
                    theta = NULL,
                    nu = c(2,2),
                    lambda = c(2,2),
                    psi = c(2,2),
                    kappa = c(2,2),
                    xi = c(2,2),
                    phi = NULL,
                    rho = NULL,
                    sigma = NULL,
                    tau = NULL,
                    delta = NULL,
                    gamma = NULL,
                    psi = NULL)
#likelihood function for the model
ll <- function(params, data) {
  alpha <- params["alpha"]
  beta <- params["beta"]
  nu <- params["nu"]
  lambda <- params["lambda"]
  psi <- params["psi"]
  kappa <- params["kappa"]
  xi <- params["xi"]
  mu <- params["mu"]
  theta <- params["theta"]
  phi <- params["phi"]
  rho <- params["rho"]
  sigma <- params["sigma"]
  tau <- params["tau"]
  delta <- params["delta"]
  gamma <- params["gamma"]
  psi <- params["psi"]
  
  spec <- ugarchspec(variance.model = list(model = "RHGARCH",
                                           garchOrder = c(kappa, xi),
                                           harchOrder = c(nu, lambda),
                                           mci = FALSE),
                     mean.model = list(armaOrder = c(theta, phi), include.mean = mu),
                     distribution.model = "std")
  fit <- ugarchfit(spec, data = data, solver = "hybrid", solver.control = list(trace = 0))
  
  loglik <- sum(dnorm(data, mean = fit@model$mu, sd = fit@sigma, log = TRUE))
  
  prior_alpha <- dgamma(alpha, shape = 2, rate = 2, log = TRUE)
  prior_beta <- dgamma(beta, shape = 2, rate = 2, log = TRUE)
  prior_nu <- dgamma(nu, shape = 2, rate = 2, log = TRUE)
  prior_lambda <- dgamma(lambda, shape = 2, rate = 2, log = TRUE)
  prior_psi <- dgamma(psi, shape = 2, rate = 2, log = TRUE)
  prior_kappa <- dgamma(kappa, shape = 2, rate = 2, log = TRUE)
  prior_xi <- dgamma(xi, shape = 2, rate = 2, log = TRUE)
  prior_mu <- dnorm(mu, mean = 0, sd = 10, log = TRUE)
  prior_theta <- dnorm(theta, mean = 0, sd = 10, log = TRUE)
  prior_phi <- dnorm(phi, mean = 0, sd = 10, log = TRUE)
  prior_rho <- dnorm(rho, mean = 0, sd = 10, log = TRUE)
  prior_sigma <- dgamma(sigma, shape = 2, rate = 2, log = TRUE)
  prior_tau <- dgamma(tau, shape = 2, rate = 2, log = TRUE)
  prior_delta <- dgamma(delta, shape = 2, rate = 2, log = TRUE)
  prior_gamma <- dgamma(gamma, shape = 2, rate = 2, log = TRUE)
  prior_psi <- dgamma(psi, shape = 2, rate = 2, log = TRUE)
  
  prior <- prior_alpha + prior_beta + prior_nu + prior_lambda + prior_psi + prior_kappa +
    prior_xi + prior_mu + prior_theta + prior_phi + prior_rho + prior_sigma + prior_tau + 
    prior_delta + prior_gamma + prior_psi
  
  return(list(loglik = loglik, prior = prior))
}

#Run the MCMC algorithm
model <- MCMCmodel(ll = ll, prior = prior, data = y, params = c("alpha", "beta", "nu", "lambda", "psi", "kappa",
                                                                "xi", "mu", "theta", "phi", "rho", "sigma",
                                                                "tau", "delta", "gamma", "psi"), 
                   iter = 10000, warmup = 2000, chains = 4, cores = 4, thin = 1)
fit <- runMCMC(model)

#Compute posterior quantiles
q <- c(0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99)
forecasts <- ugarchforecast(spec, data = y, n.ahead = 1, n.roll = 0, 
                            forecast.length = 1, always.use.starting.values = TRUE,
                            horizon = 1, VaR.alpha = q, VaR.cr = NULL, 
                            VaR.model = NULL, seed = NULL, rec.init = TRUE, 
                            solver = "hybrid", solver.control = list(trace = 0))
posterior_quantiles <- forecasts@VaR$q[q]

  
