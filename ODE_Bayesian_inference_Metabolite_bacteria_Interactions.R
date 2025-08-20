# Libraries needed
library(deSolve)
library(rstan)
library(tidyr)
library(ggplot2)
library(dplyr)
library(posterior)
library(stringr)

# Prepare for running stan optimally
options(mc.cores = parallel::detectCores())
rstan_options(threads_per_chain = 2)

#--------------------------------------------------------------------------------------
# Generating a bacterial community including metabolic interaction and dynamics 
# in silico using a set of coupled ordinary differential equations
#--------------------------------------------------------------------------------------
{# Simulate metabolic and bacterial dynamics using a set of ODEs 
  bacteria_metabolite_ode <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      N_met <- length(S_t0)
      N_bac <- length(y0_bacteria)
      B_total <- sum(state[(N_met + 1):(N_met + N_bac)])
      
      dydt <- numeric(N_met + N_bac)
      
      # Metabolite dynamics (produced or consumed)
      for (i in 1:N_met) {
        dydt[i] <- 0
        for (j in 1:N_bac) {
          dydt[i] <- dydt[i] - Cij[i, j] * state[i] * state[N_met + j]
        }
      }
      
      # Bacteria dynamics (Only consumed metabolites contribute to bacterial growth)
      for (j in 1:N_bac) {
        growth <- -d_bac
        for (i in 1:N_met) {
          growth <- growth + C0 * Cij[i, j] * state[i]
        }
        dydt[N_met + j] <- state[N_met + j] * growth * (k - B_total)
      }
      
      list(dydt)
    })
  }
  
  # Constants
  N_met <- 10   # Number of metabolites
  N_bac <- 3    # Number of bacteria
  Reps <- 6     # Number of "experimental" replicates
  Obs <- 4      # Number of "experimental" observations
  k <- 14       # Carrying capacity
  delta <-1e-6  # Small number to ensure numerical stability
  C0 <- 0.1     # Metabolism conversion coefficient
  d_bac <- 0.1  # Death rate per bacterium
  n_sims <- 10  # Number of independent simulations, each simulation uses a separate set of Cij values
  t0 <- 0       # Starting point time series
  tmax <- 4     # End point time series
  ts <- seq(1e-6, tmax, length.out = Obs) # Time steps for the ODEs
  set.seed(42)
  
  # Noise addition
  add_noise <- FALSE      # If TRUE adds random noise
  use_ar1_noise <- FALSE  # If TRUE adds autocorrelated noise across observations per replicate 
  rho_sim <- 0.8          # Parameter needed if use_ar1_noise is TRUE
  Sigma_met <- runif(N_met, min = 0.1, max = 0.5)   # Noise added to the metabolites
  Sigma_bac <- runif(N_bac, min = 0.01, max = 0.05) # Noise added to the bacteria
  Sigma_noise <- c(Sigma_met, Sigma_bac)
  
  # Determine the amount of metabolites per category, see N_met_sim
  # Make sure that the order of metabolites is consumed --> produced
  consumed_metabolites <- c(1:7)
  produced_metabolites <- c(8:10)
  
  # Initialize storage
  Cij_list <- vector("list", n_sims)
  S_t0_list_all <- vector("list", n_sims)
  y0_bacteria_list_all <- vector("list", n_sims)
  y_obs_all <- array(0, dim = c(n_sims, Reps, Obs, N_met + N_bac))
  
  for (sim in 1:n_sims) {
    # Generate true Cij values
    Cij_consumed <- matrix(runif(length(consumed_metabolites) * N_bac, 0, 1),
                           nrow = length(consumed_metabolites), ncol = N_bac)
    Cij_produced <- matrix(runif(length(produced_metabolites) * N_bac, -1, 0),
                           nrow = length(produced_metabolites), ncol = N_bac)
    Cij <- matrix(0, N_met, N_bac)
    Cij[consumed_metabolites, ] <- Cij_consumed 
    Cij[produced_metabolites, ] <- Cij_produced
    
    # Save true Cij values for current simulation
    Cij_list[[sim]] <- Cij
    
    # Simulate initial conditions
    S_t0_list <- replicate(Reps, {
      c(runif(length(consumed_metabolites), min = 0.5, max = 1),    
        runif(length(produced_metabolites), min = 1e-3, max = 1e-2))
    }, simplify = FALSE)
    
    y0_bacteria_list <- replicate(Reps,
                                  runif(N_bac, min = 0.05, max = 0.2),
                                  simplify = FALSE)
    
    # Save initial conditions
    S_t0_list_all[[sim]] <- S_t0_list
    y0_bacteria_list_all[[sim]] <- y0_bacteria_list
    
    # Solve the set of ODEs
    for (r in 1:Reps) {
      state_init <- c(S_t0_list[[r]], y0_bacteria_list[[r]])
      
      parameters <- list(
        Cij = Cij,
        k = k,
        S_t0 = S_t0_list[[r]],
        y0_bacteria = y0_bacteria_list[[r]],
        C0 = C0,
        consumed_metabolites = consumed_metabolites,
        produced_metabolites = produced_metabolites
      )
      
      out <- ode(y = state_init, times = ts, func = bacteria_metabolite_ode,
                 parms = parameters, method = "ode45")
      
      for (i in 1:(N_met + N_bac)) {
        true_vals <- out[, i + 1]
        
        if (!add_noise) {
          y_obs_all[sim, r, , i] <- true_vals + delta
          
        } else if (add_noise && !use_ar1_noise) {
          y_obs_all[sim, r, , i] <- rlnorm(n = Obs, meanlog = log(true_vals + delta), sdlog = Sigma_noise[i])
          
        } else if (add_noise && use_ar1_noise) {
          noise <- numeric(Obs)
          noise[1] <- rnorm(1, 0, Sigma_noise[i] / sqrt(1 - rho^2))
          for (t in 2:T_sim) {
            noise[t] <- rho * noise[t - 1] + rnorm(1, 0, Sigma_noise[i])
          }
          y_obs_all[sim, r, , i] <- exp(log(true_vals + delta) + noise)
        }
      }
    }}
  
  # Output:
  # Cij_list:             list of n_sims Cij matrices
  # S_t0_list_all:        list of n_sims lists, each containing 6 metabolite vectors corresponding to the intial values
  # y0_bacteria_list_all: list of n_sims lists, each containing 6 bacterial vectors corresponsing to the intial values
  # y_obs_all:            4D array [n_sims simulations, 6 replicates, Obs, N_met + N_bac]
  
  # Make sure to save the generated observations and variables before parameter 
  # estimation using the Bayesian inference approach
  
  save(Cij_list, S_t0_list_all, y0_bacteria_list_all, y_obs_all, file = "...")}

#--------------------------------------------------------------------------------------
#Bayesian inference framework
#--------------------------------------------------------------------------------------
{# The Bayesian inference framework is written in Rstan
  Stan_code <- {"
functions {
  // ODE system describing bacteria-metabolite interactions
  vector Bac_Met(real t,
               vector y,
               matrix Cij,
               real k,
               vector S_t0,
               int[] consumed_metabolites,
               real C0) {
    int N_met = rows(Cij);                                  
    int N_bac = cols(Cij);                                  
    vector[N_met + N_bac] dydt;                            
    real B_total = sum(y[(N_met + 1):(N_met + N_bac)]);    
   
    // Metabolite dynamics
    for (i in 1:N_met) {
      dydt[i] = 0;
        for (j in 1:N_bac) {
          dydt[i] -= Cij[i, j] * y[i] * y[N_met + j];
      }
    }

    // Bacterial growth dynamics
    for (j in 1:N_bac) {
      real growth = -0.1;
        for (i in consumed_metabolites) {
          growth += C0 * Cij[i, j] * y[i];
  }
  dydt[N_met + j] = y[N_met + j] * growth * (k - B_total);
}
    return dydt;
  }
}

data {
  int<lower=1> R;
  int<lower=1> T;
  int<lower=1> N_met;
  int<lower=1> N_bac;
  array[T] real ts;
  array[R, T] vector[N_met + N_bac] y_obs;
  vector[N_met] S_t0[R];
  vector[N_bac] y0_bacteria[R];
  real t0;
  real<lower=0> k;
  int<lower=0> N_produced;
  int<lower=0> N_consumed;
  int produced_metabolites[N_produced];
  int consumed_metabolites[N_consumed];
  int<lower=0, upper=1> use_ar1;
  real<lower=0> delta;
  real<lower=0, upper=1> C0;
}

parameters {
  // Cij matrix, the priors for the cij matrix are defined here since those are uniform and not in the model block
  // Make sure the order of the metabolites is consumed --> produced
  matrix<lower=0, upper=1>[N_consumed, N_bac] Cij_pos;  // Uniform prior between 0 and 1
  matrix<lower=-1, upper=0>[N_produced, N_bac] Cij_neg; // Uniform prior between -1 and 0
  
  // Residual noise
  real<lower=0> sigma_ar; 
  
  // Parameters for the inclusion of noise autocorrelation
  vector[R] rho_raw;
  real mu_rho;
  real<lower=0> sigma_rho;
}

transformed parameters {
  // Inclusion of noise autocorrelation
  vector[R] rho;
  if (use_ar1) {
    rho = 2 * inv_logit(mu_rho + sigma_rho * rho_raw) - 1;
  } else {
    rho = rep_vector(0.0, R);
  }

  // Construct full Cij matrix by combining positive and negative parts
  matrix[N_met, N_bac] Cij;
  for (i in 1:N_consumed) {
    for (j in 1:N_bac) {
      Cij[i, j] = Cij_pos[i, j];
    }
  }
  
  for (i in 1:N_produced) {
    int idx = N_consumed + i;
    for (j in 1:N_bac) {
      Cij[idx, j] = Cij_neg[i, j];
    }
  }
}

model {
  array[R, T] vector[N_met + N_bac] mu;
  vector[N_met + N_bac] y0;

  for (r in 1:R) {
    y0[1:N_met] = S_t0[r];
    y0[(N_met + 1):(N_met + N_bac)] = y0_bacteria[r];
    mu[r] = ode_rk45(Bac_Met,
                     y0,
                     t0,
                     ts,
                     Cij,
                     k,
                     S_t0[r],
                     consumed_metabolites,
                     C0);
  }

// Priors for the inclusion of noise autocorrelation
rho_raw ~ std_normal();
mu_rho ~ normal(0, 1);
sigma_rho ~ exponential(2);

// Prior for the residual noise
sigma_ar ~ normal(0, 0.2);


for (r in 1:R) {
    for (l in 1:(N_met + N_bac)) {
      vector[T] eps;
      for (t in 1:T) {
        eps[t] = log(y_obs[r, t, l] + delta) - log(mu[r, t, l] + delta);
      }

      if (use_ar1 == 1) {
        target += normal_lpdf(eps[1] | 0, sigma_ar / sqrt(1 - square(rho[r])));
        for (t in 2:T) {
          target += normal_lpdf(eps[t] | rho[r] * eps[t - 1], sigma_ar);
        }
      } else {
        for (t in 1:T) {
          target += normal_lpdf(eps[t] | 0, sigma_ar);
        }
      }
    }
  }
}

generated quantities {
  array[R, T] vector[N_met + N_bac] y_sim;
  vector[N_met + N_bac] y0_sim;
  array[R, T, N_met + N_bac] real log_lik;
  real rho_pop = 2 * inv_logit(mu_rho) - 1;

  for (r in 1:R) {
    y0_sim[1:N_met] = S_t0[r];
    y0_sim[(N_met + 1):(N_met + N_bac)] = y0_bacteria[r];
    y_sim[r] = ode_rk45(Bac_Met,
                        y0_sim,
                        t0,
                        ts,
                        Cij,
                        k,
                        S_t0[r],
                        consumed_metabolites,
                        C0);

    for (l in 1:(N_met + N_bac)) {
      for (t in 1:T) {
        real eps = log(y_obs[r, t, l] + delta) - log(y_sim[r, t, l] + delta);
        if (use_ar1 == 1) {
          if (t == 1) {
            log_lik[r, t, l] = normal_lpdf(eps | 0, sigma_ar / sqrt(1 - square(rho[r])));
          } else {
            real prev_eps = log(y_obs[r, t - 1, l] + delta) - log(y_sim[r, t - 1, l] + delta);
            log_lik[r, t, l] = normal_lpdf(eps | rho[r] * prev_eps, sigma_ar);
          }
        } else {
          log_lik[r, t, l] = normal_lpdf(eps | 0, sigma_ar);
        }
      }
    }
  }
}
"}

# Compiling and parameter estimation is written in R
# Make sure to have saved the input variable data, the datalist is not saved here
# Compile the model
Stan_code <- stan_model(model_code = Stan_code)

# Storage objects
fits_test <- vector("list", n_sims)
true_Cij_list <- vector("list", n_sims)
posterior_Cij_list <- vector("list", n_sims)
n_sims <- 1 # Number of simulations

# Parameter estimation
for (i in 1:n_sims) {
  S_t0_i <- do.call(rbind, S_t0_list_all[[i]])
  y0_bacteria_i <- do.call(rbind, y0_bacteria_list_all[[i]])
  y_obs_i <- y_obs_all[i, , , ]
  
  # Create a data list for Stan
  data_list_Stan <- list(
    N_met = ncol(S_t0_i),
    N_bac = ncol(y0_bacteria_i),
    T = dim(y_obs_i)[2],
    R = dim(y_obs_i)[1],
    ts = ts_sim,
    y_obs = y_obs_i,
    S_t0 = S_t0_i,
    y0_bacteria = y0_bacteria_i,
    t0 = 0,
    k = 14,
    N_produced = length(produced_metabolites),
    N_consumed = length(consumed_metabolites),
    produced_metabolites = produced_metabolites,
    consumed_metabolites = consumed_metabolites,
    C0 = 0.1,
    use_ar1 = 0,
    delta = delta_sim
  )
  
  # Fit Stan model
  fit <- sampling(Stan_code,
                  data = data_list_Stan,
                  iter = 4000,
                  chains = 4,
                  cores = 4,
                  warmup = 2000,
                  control = list(adapt_delta = 0.99))
  
  # Save the fit object
  fits_test[[i]] <- fit
  
  # Save true Cij values used for this simulation
  true_Cij_list[[i]] <- Cij_list[[i]]
  
  # Extract posterior median for Cij
  Cij_post <- rstan::extract(fit, pars = "Cij")$Cij
  posterior_Cij_list[[i]] <- apply(Cij_post, c(2,3), median)
}
saveRDS(fit, file = ".......rds")
saveRDS(true_Cij_list, file = "....rds")
saveRDS(posterior_Cij_list, file = "....rds")
}
