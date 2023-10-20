# Prepare the data
data <- fun.create_data_goyal("annualy")

# Load the prepared data
data <- readRDS(paste0(path, "data/2.real_data/tbl.goyal_annual_v3.rds"))

# Choose variables to analyze
vars_in <- c("b.m", "ep.sp", "tms", "dfy", "dg.sp")

# Chose a specific variable
indicator <- 1

# Save it
var_in <- sym(vars_in[indicator])

# Time span
data <- data %>% 
  dplyr::mutate(x_in = !!var_in) %>% 
  filter(between(date, left = as.Date("1953-01-01"), right = as.Date("2022-12-01")))

# Save the length of the data
T_len <- nrow(data)

# Prepare the data for the estimation
data.final <- list(x_latent   = data$x_in[-T_len],
                   y_latent   = data$x_in[-1],
                   x_observe  = data$x_in[-T_len],
                   y_observe  = data$ret.crsp.l[-1],
                   dg_observe = NULL)

# Settings
MC_draws <- 1e4
burnin   <- 1e3
thining  <- 30

# Priors
Prior_cont_ar_in$a_R            <- c(0.1, 0.5) 
Prior_cont_ar_in$b_R            <- as.character(1)
Prior_cont_ar_in$Prior_variance <- "r2"

# Estimation
res_ar_cont <- f_Bayesian_control_function_Cholesky_AR_prior_R2(
  Data           = data.final,
  Prior          = Prior_cont_ar_in,
  mc             = MC_draws,
  burnin         = burnin,
  start_true_ind = F,
  thin           = thining,
  true_data      = NULL)








