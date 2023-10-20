# Parameters
MC_sim_data <- 1e3  # Number of datasets to simulate
MC_draws    <- 1e4  # Number of the MCMC draws
burnin      <- 1e3
thining     <- 30

# Choose the DGP specification
indicator     <- 1

# DGP parameters
T_len_value   <- c(100)
phi_values    <- c(0.95)
beta_values   <- c(0, 0.025, 0.05, 0.075, 0.1, 0.2)

# Generate possible values
possible_values <- expand.grid(T_len_value, phi_values, beta_values)

# Priors
Prior_cont_ar_in$a_R            <- c(0.1, 0.5) 
Prior_cont_ar_in$b_R            <- as.character(1)
Prior_cont_ar_in$Prior_variance <- "r2"

## Simulate data ##
data_4 <- list(T_len    = possible_values[indicator, 1],
               alpha    = 0.6,
               beta     = possible_values[indicator, 3],
               theta    = -3*(1 - possible_values[indicator, 2]),
               phi      = possible_values[indicator, 2],
               Sigma    = matrix(c(0.02, -0.02, -0.02, 0.04), nrow = 2))

# Simulate the data
data_sim_4 <- f_con_ar_sim_data(mc = MC_sim_data, data = data_4, FUN = f_data)

## Estimation ##

res_ar_cont <- f_con_ar_sim_mc_function(Data             = data_sim_4, 
                                        FUN              = f_Bayesian_control_function_Cholesky_AR_prior_R2, 
                                        Prior            = Prior_cont_ar_in,
                                        mc               = MC_draws, 
                                        burnin           = burnin,
                                        start_true_ind   = F,
                                        thin             = thining,
                                        true_data        = data_4)


## Results ##

# Convergence level
level <- 700

# Obtain the results
res <- fun.create_result(data_in       = res_ar_cont, 
                         level         = level, 
                         prior_r2_type = Prior_cont_ar_in$Prior_variance)

res
