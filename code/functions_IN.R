#### Packages ####

# VARs estimation
library(vars) 

# Bayesian IV package
library(bayesm)

# Estimate the mode of the MCMC
library(hBayesDM)

# For the posterior analysis
library(mcmcOutput)

# resave objects
library(cgwtools)

# Fit the distribution
library(fitdistrplus)

# For the correlation plot
library(corrplot)

# For the nice tables
library(stargazer)

# Time Series package
library(tseries)

# Time Series package
library(forecast)

# time-series
library(astsa)

# For the winsorize function
library(DescTools)

# For plotting
library(ggplot2)

# For plotting several plots one on grid
library(gridExtra)
library(ggpubr)

# For nicer axes of figures
library(scales)      

# Reshape data for ggplot
library(reshape2)

# For working with dates
library(lubridate)  

# Dynamic regression
library(dyn)

# For rolling sum
library(RcppRoll)

# For the ADF test
library(urca)
  
# HDP interval
library(HDInterval) 

# For effective sample size
library(coda)

# Simulate from the inverse gamma distribution
library(invgamma)

# Simulate from the multivariate normal distribution
library(mvtnorm)

# Simulate from the beta prime distribution
library(extraDistr)

# Draw colorful output
library(crayon)

# Dplyr friendly functions
library(tidyr)

# For grammar
library(tidyverse) 

# For string
library(stringr) 

# For skewness
library(moments) 

# For RMD
library(markdown) 

# For RMD
library(knitr)

col_in          <- c("gold", "magenta", "cyan",  "steelblue", "tomato", "green")

shape_in <- c("OLS"  = 3, 
              "RBE"  = 4, 
              "BAY"  = 17)


#### Real Data ####

fun.create_data_goyal <- function(frequency_in){
  
  if(frequency_in == "annualy"){
    
    tbl.goyal <- read.csv2(paste0(path, "data/2.real_data/annually_2021.csv"), na.strings="NaN", stringsAsFactors=FALSE, sep = ";")
    
    # Change to tibble
    tbl.goyal <- as_tibble(apply(tbl.goyal,2,as.double), column_name = colnames(tbl.goyal))
    
    tbl.goyal <- tbl.goyal %>% 
      dplyr::rename(date             = Datum,
                    P                = Index,
                    D                = D12,
                    E                = E12,
                    ret.rf.n         = Rfree,
                    retx.crsp.n      = CRSP_SPvwx,
                    ret.crsp.n       = CRSP_SPvw
      ) %>% 
      mutate(date             = ymd(date, truncated = 2L),     # Date
             dp.sp            = log(D) - log(P),               # Dividend Price Ratio (dp)
             dy.sp            = log(D) - log(dplyr::lag(P)),   # Dividend Yield (dy) 
             ep.sp            = log(E) - log(P),               # Earnings Price Ratio (ep) 
             de.sp            = log(D) - log(E),               # Dividend Payout Ratio (de)
             
             ret.infl.n       = infl,                        # Inflation return adjustment
             tms              = lty - ret.rf.n,              # Term Spread (tms)   
             dfy              = BAA - AAA,                   # Default Yield Spread (dfy)
             dfr              = corpr - ltr,                 # Default Yield Spread (dfr)
             
             dg.sp            = log(D) - log(dplyr::lag(D)),     # Dividend growth
             
             ret.sp.tot       = (P + D) / dplyr::lag(P),          # Total-Return SP 500
             ret.sp.l         = log(ret.sp.tot),                  # Log-Return SP 500
             ret.sp.n         = ret.sp.tot - 1,                   # Net-Return SP 500
             
             ret.crsp.l       = log(1 + ret.crsp.n),         # CRSP log return incl D
             retx.crsp.l      = log(1 + retx.crsp.n),        # CRSP log return excl D
             dp.crsp          = log((1 + ret.crsp.n)/(1 + retx.crsp.n) - 1), # Create a annual dividend price ratio
             dg.crsp          = dp.crsp - dplyr::lag(dp.crsp) + retx.crsp.l, # Create a annual dividend growth
             
             ret.crsp.n.cpi   = ret.crsp.n - ret.infl.n,                  # Real CRSP return
             ret.sp.n.cpi     = ret.sp.n - ret.infl.n,                    # Real return
             
             ret.crsp.n.ex    = ret.crsp.n - ret.rf.n,              # CRSP excess return 
             ret.sp.n.ex      = ret.sp.n - ret.rf.n,                # Excess return 
             
             ret.crsp.n.ex.cpi   = ret.crsp.n.ex - ret.infl.n,     # Real CRSP excess return 
             ret.sp.n.ex.cpi     = ret.sp.n.ex - ret.infl.n   # Real excess return
      )   
    
    
  }
  else{
    return("Please choose annualy")
  }
  
  return(tbl.goyal)
  
}


#### Simulated Data ####

# Applies data to the function n-times
f_con_ar_sim_mc_function <- function(Data, FUN, ...){
  return(lapply(Data, function(x) FUN(Data = x$reg_Bayesian, ...))) 
}

# Simulates mc datasets
f_con_ar_sim_data <- function(mc, data, FUN){
  replicate(n = mc, expr = FUN(data), simplify = F)
}

# Creates the data)
f_data <- function(data){
  
  # Save values
  T_len    = data$T_len
  alpha    = data$alpha
  beta     = data$beta
  theta    = data$theta
  phi      = data$phi
  Sigma    = data$Sigma
  
  # Create an empty vectors
  x <- {}
  y <- {}
  
  # Matrix of the errors
  eps = rmvnorm(n = T_len, mean = rep(0,2), sigma = Sigma)
  
  # First observation for the x
  x[1] <- rnorm(n = 1, mean = theta / (1 - phi), sd = sqrt(Sigma[1,1] / (1 - phi^2)))
  y[1] <- NA
  
  
  for(i in 2:T_len){
    # Simulate X
    x[i] <- theta + phi *x[i-1] + eps[i, 1]
    # Simulate y
    y[i] <- alpha + beta *x[i-1] + eps[i, 2]
  }
  
  w          <- matrix(rep(1, T_len), ncol = 1)
  z          <- cbind(1,x)
  
  return(list(data_init  = list(x          = x, 
                                y          = y,
                                w          = w, 
                                error      = eps,
                                alpha      = alpha,
                                beta       = beta),
              reg_observ = list(x          = x[-T_len], 
                                y          = y[-1],
                                w          = w[-1], 
                                error      = eps[-1,]),
              reg_latent = list(z          = z[-T_len,], 
                                x          = x[-1], 
                                error      = eps[-1,]),
              # Not sure
              reg_Bayesian = list(y_observe          = y[-1],
                                  x_observe          = x[-T_len], 
                                  y_latent           = x[-1], 
                                  x_latent           = x[-T_len], 
                                  w                  = w[-1, , drop = FALSE],
                                  error              = eps[-1,])))
}


#### Priors ####

Prior_cont_ar_in     <- list(bias_cor           = F,
                             Prior_type         = "reference",
                             Prior_variance     = "sig2_x",
                             boundary_bias_cor  = 0.99,
                             B_mu               = 0.2,
                             B_0_22             = 10^8,
                             B_0_11             = 10^12,
                             
                             mu_0_alpha_2       = 0,
                             mu_0_beta          = 0,
                             
                             Sigma_0_alpha_2    = 10,
                             Sigma_0_beta       = 0.01,
                             
                             Sigma_0_psi        = 10,  # for the Cholesky model
                             mu_psi             = 0,   # for the Cholesky model
                             
                             a_R                = 0.5,
                             b_R                = 3,
                             
                             v_0_x              = 4,
                             Sigma_0_x          = 0.06,
                             v_0_v              = 2.5,
                             Sigma_0_v          = 0.03)


#### Estimation ####

# Main function
f_Bayesian_control_function_Cholesky_AR_prior_R2 <- function(Data, 
                                                             Prior,
                                                             mc, 
                                                             burnin,
                                                             start_true_ind,
                                                             true_data = NULL,
                                                             thin = 10,
                                                             save_warnings = FALSE){
  
  ## Data ##
  
  x_lat <- Data$x_latent
  y_lat <- Data$y_latent
  
  x_obs <- Data$x_observe
  y_obs <- Data$y_observe
  
  data_bias <- list(x = x_lat, y = y_lat)
  
  ## Frequentist ##
  
  # bias corrected #
  res_bias_corrected <- f_bias_correction(data = data_bias, type = "corrected")
  
  # ols #
  ols <- summary(lm(y_obs ~ x_obs))
  
  # ols corrected #
  ols_corrected <- summary(lm(y_obs ~ x_obs + res_bias_corrected$error_corrected))
  
  # If we have a simulated data we can write the correct OLS conditioning on the errors
  if(!is.null(Data$error[,1])){
    ols_true <- summary(lm(y_obs ~ x_obs + Data$error[,1]))    
  }else{
    ols_true <- NULL
  }
  
  
  # Compute needed values
  N <- length(y_obs)
  #p <- ncol(x_obs)
  p <- 1
  
  #############################
  ## Priors Control Function ##
  #############################
  
  # alpha_2
  mu_0_alpha_2    <- Prior$mu_0_alpha_2
  Sigma_0_alpha_2 <- Prior$Sigma_0_alpha_2
  
  # beta
  mu_0_beta       <- Prior$mu_0_beta
  
  # psi
  mu_0_psi        <- Prior$mu_psi
  Sigma_0_psi     <- Prior$Sigma_0_psi
  
  # sig2_x
  v_0_x           <- Prior$v_0_x
  Sigma_0_x       <- Prior$Sigma_0_x
  
  # sig2_v      
  v_0_v           <- Prior$v_0_v
  Sigma_0_v       <- Prior$Sigma_0_v
  
  # R^2
  Prior_a_R <- ifelse(length(Prior$a_R) > 1, TRUE, FALSE)
  Prior_b_R <- ifelse(length(Prior$b_R) > 1, TRUE, FALSE)
  
  if(Prior_a_R & Prior_b_R){
    
    a_R  <- rep(NA_real_, mc)
    b_R  <- rep(NA_real_, mc)
    
    # Save two values
    a_R_in <- Prior$a_R
    b_R_in <- Prior$b_R
    
  }else if(Prior_a_R){
    
    a_R  <- rep(NA_real_, mc)
    b_R  <- rep(as.double(Prior$b_R), mc)
    
    # Save two values
    a_R_1 <- Prior$a_R[1]
    a_R_2 <- Prior$a_R[2]
    
    
  }else{ 
    
    return("Please use the previous version of the function that treats b_R as fixed")
    
      }
  
  
  
  # Pre-calculate
  shape_x_n_con   <- v_0_x + N/2
  shape_v_n_con   <- v_0_v + N/2
  
  # Depends on the Jeffreys or not
  if(Prior$Prior_type == "jeffreys_3"){
    shape_x_n_con <- 1/2 + N/2
    Sigma_0_x     <- 0
  }
  
  
  #############################
  ##        Priors AR        ##
  #############################
  
  ## Priors ##
  
  # Set up values equal to the biased corrected value if it doesn't exceed 0.99
  if(res_bias_corrected$phi_corrected < Prior$boundary_bias_cor){
    phi_in    <- res_bias_corrected$phi_corrected
    gamma_in  <- res_bias_corrected$theta_corrected
    
  }else if(res_bias_corrected$res$coefficients[2,1] < Prior$boundary_bias_cor){
    phi_in    <- res_bias_corrected$res$coefficients[2,1]
    gamma_in  <- res_bias_corrected$res$coefficients[1,1]
    
  }else{ ###### !
    phi_in    <- Prior$boundary_bias_cor
    gamma_in  <- mean(c(x_lat[1],y_lat)) * (1 - phi_in)
  }
  
  if(isTRUE(Prior$bias_cor)){
    
    b_mu    <- gamma_in/(1 - phi_in)
    b_0     <- Prior$N_obs*(1 - phi_in)
    a_0     <- Prior$N_obs*phi_in
    
    # We fix b_0 to 1
    if(isTRUE(Prior$b_0_fixed)){
      b_0     <- 1
      a_0     <- phi_in / (1 - phi_in)  
    }
    
    
  }else{
    # We chose b_mu equals the b_mu in the data
    b_mu    <- mean(c(x_lat[1],y_lat))
  }
  
  B_mu      <- Prior$B_mu
  
  # Auxiliary priors
  b_0_mu_phi     <- matrix(c(0,0), nrow = 2)
  B_0            <- diag(c(Prior$B_0_11, Prior$B_0_22))
  B_0_inv        <- solve(diag(c(Prior$B_0_11, Prior$B_0_22)))
  
  #############################
  ##    Store the values     ##
  #############################
  
  alpha_2        <- rep(NA_real_, mc)
  beta           <- rep(NA_real_, mc)
  sig2_x         <- rep(NA_real_, mc)
  sig2_v         <- rep(NA_real_, mc)
  psi            <- rep(NA_real_, mc)
  sig2_x_tilde   <- rep(NA_real_, mc) # need just inside if the matrix fails
  
  sig2_eta       <- rep(NA_real_, mc)
  Z_eta          <- rep(NA_real_, mc)
  Sigma_0_beta   <- rep(NA_real_, mc)
  
  gamma          <- rep(NA_real_, mc)
  phi            <- rep(NA_real_, mc)
  mu             <- rep(NA_real_, mc)
  
  mu_2_beta_n    <- rep(NA_real_, mc)
  Sigma_2_beta_n <- rep(NA_real_, mc)
  
  
  warnings_message <- matrix(data = NA_real_, ncol = 5, nrow = 1) 
  
  acceptance_phi_gamma = 0
  acceptance_psi       = 0
  acceptance_sig2_x    = 0
  acceptance_sig2_v    = 0
  
  #############################
  ##      Precalculate       ##
  #############################
  
  mu_2_0     <- matrix(c(mu_0_alpha_2, mu_0_beta), ncol = 1)
  
  c_T_ar <- N/2
  
  X_1_ar <- cbind(1, x_lat)
  B_T_ar <- solve(crossprod(X_1_ar, X_1_ar) + B_0_inv)
  
  
  #############################
  ##   Draw initial values   ##
  #############################
  #set.seed(1)
  
  if(isTRUE(start_true_ind)){
    
    alpha_2[1]       <- true_data$alpha
    beta[1]          <- true_data$beta
    
    sig2_x[1]        <- true_data$Sigma[1,1]
    psi[1]           <- true_data$Sigma[1,2] / sig2_x[1]
    sig2_v[1]        <- true_data$Sigma[2,2] - true_data$Sigma[1,1] * psi[1]^2
    sig2_x_tilde[1]  <- sig2_x[1] * sig2_v[1] / ( sig2_v[1] + sig2_x[1] * psi[1]^2)
    
    phi[1]           <- true_data$phi
    gamma[1]         <- true_data$theta
    mu[1]            <- gamma[1]/(1 - phi[1])
    
    Z_eta[1]         <- 2   # We chose this value as if in the case when Z_eta is fixed
    sig2_eta[1]      <- 0.1 # 
    Sigma_0_beta[1]  <- sig2_eta[1] * (sig2_v[1] / sig2_x[1] + psi[1]^2) * (1 - phi[1]^2) # deduced
    
    # Sample a_R and b_R
    if(isTRUE(Prior_a_R) & isTRUE(Prior_b_R)){
      
      # We sample ind.beta_prior that basically tells us which value to take
      ind.beta_prior   <- sample(x = 1:length(Prior$a_R), size = 1)
      a_R[1]           <- Prior$a_R[ind.beta_prior]
      b_R[1]           <- Prior$b_R[ind.beta_prior]
      
    }else{
      
      # We are on the safe side if we have 1 or 2 values
      a_R[1]           <- sample(x    = Prior$a_R, 
                                 size = 1)
      
    }
    
    
    
  }else{
    
    alpha_2[1]       <- rnorm(n = 1, mean = mu_0_alpha_2, sd = sqrt(Sigma_0_alpha_2))
    sig2_v[1]        <- 1/rgamma(n = 1, shape = v_0_v, rate = Sigma_0_v)
    psi[1]           <- rnorm(n = 1, mean = mu_0_psi, sd = sqrt(Sigma_0_psi))
    sig2_x_tilde[1]  <- 1 
    
    if(Prior$Prior_type == "jeffreys_3"){
      sig2_x[1]        <- 1/rgamma(n = 1, shape = v_0_x, rate = Sigma_0_x + 0.1) 
    }else{
      sig2_x[1]        <- 1/rgamma(n = 1, shape = v_0_x, rate = Sigma_0_x) 
    }
    
    
    if(Prior$Prior_type == "beta"){
      # We do this because for small a_0, b_0 we sample 1 very often!!
      phi[1]        <- rbeta(n = 1, shape1 = a_0, shape2 = b_0)  
      while(phi[1] > 0.99){
        phi[1]        <- rbeta(n = 1, shape1 = a_0, shape2 = b_0)  
      }     
    }else if(Prior$Prior_type == "reference"){
      phi[1]        <- sin(pi*(runif(n = 1, min = 0, max = 1))/2)
    }else if(Prior$Prior_type == "jeffreys_3"){
      phi[1]        <- runif(1, min = 0.9, max = 1) 
    }else{
      print("You forgot to specify the prior!!")
    }
    
    gamma[1]         <- rnorm(n = 1, mean = b_mu*(1 - phi[1]), sd = sqrt(B_mu*(1 - phi[1])^2))
    mu[1]            <- gamma[1]/(1 - phi[1]) 
    
    Z_eta[1]         <- 2   # We chose this value as if in the case when Z_eta is fixed
    sig2_eta[1]      <- 0.1 # 
    Sigma_0_beta[1]  <- sig2_eta[1] * (sig2_v[1] / sig2_x[1] + psi[1]^2) * (1 - phi[1]^2) # deduced
    
    # Sample a_R and b_R
    if(Prior_a_R & Prior_b_R){
      
      # We sample ind.beta_prior that basically tells us which value to take
      ind.beta_prior   <- sample(x = 1:length(Prior$a_R), size = 1)
      a_R[1]           <- Prior$a_R[ind.beta_prior]
      b_R[1]           <- Prior$b_R[ind.beta_prior]
      
    }else{
      
      # We are on the safe side if we have 1 or 2 values
      a_R[1]           <- sample(x    = Prior$a_R, 
                                 size = 1)
      
    }
    
    beta[1]          <- rnorm(n = 1, mean = mu_0_beta, sd = sqrt(Sigma_0_beta[1]))
    
  }
  
  #############################
  ##        For loop         ##
  #############################
  #i=2
  cat(green("I start, hope it works!" ,Sys.time() ,"\n"))
  for(i in 2:mc){
    
    #############################
    ## Observational Equation  ##
    #############################
    
    # Create Prior for the variance
    Sigma_2_0_inv <- diag(c(1/Sigma_0_alpha_2, 1/Sigma_0_beta[i-1]))
    
    # Create the regressor matrix
    X_2_con <- cbind(1, x_obs)
    
    # Find the sig2_y
    sig2_y <- sig2_v[i-1] + sig2_x[i-1] * psi[i-1]^2
    
    # Find the parameters
    Sigma_2_n <- chol2inv(chol(crossprod(X_2_con, X_2_con)/sig2_y + Sigma_2_0_inv))
    mu_2_n    <- Sigma_2_n %*% (crossprod(X_2_con, y_obs)/sig2_y + Sigma_2_0_inv %*% mu_2_0)
    
    # Save if we have NA
    if(is.na(sig2_v[i-1])){
      sig2_v[i-1] <- sig2_v[i-2]
      print(paste0("Itteration:",i,". Problem with sig2_v"))
    }
    
    # Sampling
    draw <- rmvnorm(n     = 1, 
                    mean  = mu_2_n, 
                    sigma = Sigma_2_n)
    
    # Assign
    alpha_2[i] <- draw[1]
    beta[i]    <- draw[2]
    
    mu_2_beta_n[i]    <- mu_2_n[2]
    Sigma_2_beta_n[i] <- Sigma_2_n[2,2]
    
    #############################
    ##     Latent Equation     ##
    #############################
    
    
    # Find u_2
    u_2 <- y_obs - x_obs*beta[i] - alpha_2[i]
    
    # Find \tilde{u}_i
    u_tilde <- u_2*sig2_x[i-1] * psi[i-1] / (sig2_v[i-1] + sig2_x[i-1]*psi[i-1]^2)
    
    # Find \tilde{\sigma}^2_v
    sig2_x_tilde[i] <- sig2_x[i-1]*sig2_v[i-1]/(sig2_v[i-1] + sig2_x[i-1]*psi[i-1]^2)
    
    # Find b_t_ar
    b_t_ar <- B_T_ar %*% (crossprod(X_1_ar, y_lat - u_tilde) + B_0_inv %*% b_0_mu_phi)
    
    # If we have an NA we save why this has happened
    if(is.na(sig2_x_tilde[i])){
      sig2_x_tilde[i] <- sig2_x_tilde[i-1]
      print(paste0("Iteration:",i,". Problem with sig2_x_tilde. Prior_type:", Prior$Prior_type))
    }
    
    # Sample phi and gamma
    gamma_phi_prop <- rmvnorm(n     = 1,
                              mean  = b_t_ar, 
                              sigma = B_T_ar*sig2_x_tilde[i])
    
    # Just assign values
    gamma_prop <- gamma_phi_prop[1]
    phi_prop   <- gamma_phi_prop[2]
    mu_prop    <- gamma_prop/(1 - phi_prop)
    
    
    # Prior for x_0
    part_1 <- dnorm(x = x_lat[1], mean = mu_prop, sd = sqrt(sig2_x[i-1]/(1 - phi_prop^2)), log = T) -
      dnorm(x = x_lat[1], mean = mu[i-1], sd = sqrt(sig2_x[i-1]/(1 - phi[i-1]^2)), log = T)
    
    # Prior for phi
    if(Prior$Prior_type == "beta"){
      
      part_2 <- dbeta(x = phi_prop, shape1 = a_0, shape2 = b_0, log = T) - 
        dbeta(x = phi[i-1], shape1 = a_0, shape2 = b_0, log = T) +
        dnorm(x = gamma_prop, mean = b_mu*(1 - phi_prop), sd = sqrt(B_mu*(1 - phi_prop)^2), log = T) -
        dnorm(x = gamma[i-1], mean = b_mu*(1 - phi[i-1]), sd = sqrt(B_mu*(1 - phi[i-1])^2), log = T)
      
    }else if(Prior$Prior_type == "reference"){
      
      part_2 <- log(f_prior_Berger(phi = phi_prop)) - 
        log(f_prior_Berger(phi = phi[i-1])) +
        dnorm(x = gamma_prop, mean = b_mu*(1 - phi_prop), sd = sqrt(B_mu*(1 - phi_prop)^2), log = T) -
        dnorm(x = gamma[i-1], mean = b_mu*(1 - phi[i-1]), sd = sqrt(B_mu*(1 - phi[i-1])^2), log = T)
      
    }else if(Prior$Prior_type == "jeffreys_3"){
      part_2 <- log(f_prior_Jef_3_parameters(phi = phi_prop, T_len = N)) -
        log(f_prior_Jef_3_parameters(phi = phi[i-1], T_len = N))
      
    }else{
      print("You forgot to specify the prior!!")
    }
    
    # Auxiliary priors part
    part_3 <- ( (gamma_phi_prop - t(b_0_mu_phi)) %*% B_0_inv %*% ( t(gamma_phi_prop) - b_0_mu_phi) - 
                  (cbind(gamma[i-1], phi[i-1]) - t(b_0_mu_phi)) %*% B_0_inv %*% ( rbind(gamma[i-1], phi[i-1]) - b_0_mu_phi))/ (2*sig2_x_tilde[i])
    
    # Correction for \beta
    Sigma_0_beta_prop <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi_prop, 
                                          sig2_v_in   = sig2_v[i-1], 
                                          sig2_x_in   = sig2_x[i-1],
                                          psi_in      = psi[i-1])
    
    part_4         <- dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_prop), log = T) -
      dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta[i-1]), log = T)
    
    # Combine all parts
    r_phi_gamma <-  part_1 + part_2 + part_3 + part_4
    
    
    # Save if we have NA
    if(is.na(r_phi_gamma)){
      warnings_message <- rbind(warnings_message, c(i, phi_prop, phi[i-1], gamma_prop, gamma[i-1]))
      r_phi_gamma <- -Inf
    }
    
    # Draw an indicator
    u_phi_gamma <- log(runif(n = 1, min = 0, max = 1))
    
    # Accept/reject
    if((r_phi_gamma > u_phi_gamma)){
      acceptance_phi_gamma = acceptance_phi_gamma + 1
      phi[i]   <- phi_prop
      gamma[i] <- gamma_prop
    }else{
      phi[i]   <- phi[i-1]
      gamma[i] <- gamma[i-1] 
    }
    
    # Find the mu
    mu[i] <- gamma[i]/(1 - phi[i])
    

    #############################
    ##           PSI           ##
    #############################
    
    # Find the errors
    u_1 <- y_lat - gamma[i] - x_lat * phi[i]
    
    # Find the posterior mean and variance
    Sigma_t_psi <- 1 / (1 / Sigma_0_psi + sum(u_1^2) / sig2_v[i-1])  
    mu_t_psi    <- Sigma_t_psi * ( mu_0_psi / Sigma_0_psi + sum(u_1 * u_2) / sig2_v[i-1] )
    
    # Sample proposed psi  
    psi_prop    <- rnorm(n    = 1, 
                         mean = mu_t_psi, 
                         sd   = sqrt(Sigma_t_psi))
    
    
    # MH step
    Sigma_0_beta_prop <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi[i], 
                                          sig2_v_in   = sig2_v[i-1], 
                                          sig2_x_in   = sig2_x[i-1],
                                          psi_in      = psi_prop)
    Sigma_0_beta_old  <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi[i], 
                                          sig2_v_in   = sig2_v[i-1], 
                                          sig2_x_in   = sig2_x[i-1],
                                          psi_in      = psi[i-1])
    
    r_psi          <- dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_prop), log = T) -
      dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_old), log = T)
    
    # Draw an indicator
    u_psi <- log(runif(n = 1, min = 0, max = 1))
    
    # Accept/reject
    if((r_psi > u_psi)){
      acceptance_psi = acceptance_psi + 1
      psi[i]   <- psi_prop
    }else{
      psi[i]   <- psi[i-1]
    }

    #############################
    ##          Sig2_x         ##
    #############################
    
    sig2_x_prop <- 1/rgamma(n     = 1, 
                            shape = shape_x_n_con, 
                            rate  = Sigma_0_x + 0.5*sum(u_1^2)) 
    
    # MH step
    Sigma_0_beta_prop <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi[i], 
                                          sig2_v_in   = sig2_v[i-1], 
                                          sig2_x_in   = sig2_x_prop,
                                          psi_in      = psi[i])
    Sigma_0_beta_old  <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi[i], 
                                          sig2_v_in   = sig2_v[i-1], 
                                          sig2_x_in   = sig2_x[i-1],
                                          psi_in      = psi[i])
    
    r_sig2_x       <- dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_prop), log = T) -
      dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_old), log = T)
    
    # Draw an indicator
    u_sig2_x       <- log(runif(n = 1, min = 0, max = 1))
    
    # Accept/reject
    if((r_sig2_x > u_sig2_x)){
      acceptance_sig2_x = acceptance_sig2_x + 1
      sig2_x[i]   <- sig2_x_prop
    }else{
      sig2_x[i]   <- sig2_x[i-1]
    }

    #############################
    ##          Sig2_v         ##
    #############################
    
    
    # Find an updated v_2
    v_2 <- y_obs - x_obs*beta[i] - u_1*psi[i] - alpha_2[i]
    
    # Sample sig2_v
    sig2_v_prop <- 1/rgamma(n     = 1, 
                            shape = shape_v_n_con, 
                            rate  = Sigma_0_v + 0.5*sum(v_2^2))
    
    # MH step
    Sigma_0_beta_prop <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi[i], 
                                          sig2_v_in   = sig2_v_prop, 
                                          sig2_x_in   = sig2_x[i],
                                          psi_in      = psi[i])
    
    Sigma_0_beta_old  <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i-1], 
                                          phi_in      = phi[i], 
                                          sig2_v_in   = sig2_v[i-1], 
                                          sig2_x_in   = sig2_x[i],
                                          psi_in      = psi[i])
    
    r_sig2_v       <- dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_prop), log = T) -
      dnorm(x = beta[i], mean = 0, sd = sqrt(Sigma_0_beta_old), log = T)
    
    # Draw an indicator
    u_sig2_v       <- log(runif(n = 1, min = 0, max = 1))
    
    # Accept/reject
    if((r_sig2_v > u_sig2_v)){
      acceptance_sig2_v = acceptance_sig2_v + 1
      sig2_v[i]   <- sig2_v_prop
    }else{
      sig2_v[i]   <- sig2_v[i-1]
    }
    

    #############################
    ##     Sigma_eta^2         ##
    #############################
    
    # Sample Z_eta
    Z_eta[i] <- rgamma(n     = 1, 
                       shape = a_R[i-1] + b_R[i-1], 
                       rate  = 1 + 1/sig2_eta[i-1])
    
    # Compute needed term
    S_eta_1 <- beta[i]^2 / (2 * (1 - phi[i]^2))
    S_eta_2 <- sig2_v[i] / sig2_x[i] + psi[i]^2
    S_eta   <- S_eta_1 / S_eta_2
    
    # Sample sig2_eta
    const_c <- Z_eta[i] + S_eta
    sig2_eta[i] <- const_c / rgamma(n     = 1, 
                                    shape = b_R[i-1] + 0.5, 
                                    rate  = 1)
    
    # Finally impute Sigma_0_beta
    Sigma_0_beta[i]  <- fun.Sigma_0_beta(sig2_eta_in = sig2_eta[i], 
                                         phi_in      = phi[i], 
                                         sig2_v_in   = sig2_v[i], 
                                         sig2_x_in   = sig2_x[i],
                                         psi_in      = psi[i])
    
    
    
    #############################
    ##           a_R           ##
    #############################
    
    # If aR and bR both random
    if(Prior_a_R & Prior_b_R){
      
      # Propose a new value
      a_R_b_R_prop <- fun.proposal_a_R_b_R(a_R_old  = a_R[i-1],
                                           a_R_in   = a_R_in, 
                                           b_R_in   = b_R_in)
      
      r_a_R_b_R       <- dbetapr(x = sig2_eta[i], shape1 = a_R_b_R_prop[1], shape2 = a_R_b_R_prop[2], log = T) -
        dbetapr(x = sig2_eta[i], shape1 = a_R[i-1], shape2 = b_R[i-1], log = T)
      
      
      if(is.na(r_a_R_b_R)){
        r_a_R_b_R <- -Inf
      }
      
      # Draw an indicator
      u_a_R_b_R       <- log(runif(n = 1, min = 0, max = 1))
      
      # Accept/reject
      if((r_a_R_b_R > u_a_R_b_R)){
        a_R[i]   <- a_R_b_R_prop[1]
        b_R[i]   <- a_R_b_R_prop[2]
      }else{
        a_R[i]   <- a_R[i-1]
        b_R[i]   <- b_R[i-1]
      } 
      
    }else if(Prior_a_R){
      
      # Propose a new value
      a_R_prop <- fun.proposal_a_R_v2(a_R_old = a_R[i-1], 
                                      a_R_1   = a_R_1, 
                                      a_R_2   = a_R_2)
      
      a_R_efficient_sampl <- TRUE
      if(a_R_efficient_sampl){
        
        r_a_R       <- dbetapr(x = sig2_eta[i], shape1 = a_R_prop, shape2 = b_R[i], log = T) -
          dbetapr(x = sig2_eta[i], shape1 = a_R[i-1], shape2 = b_R[i-1], log = T) +
          fun.prior_a_R(a_R = a_R_prop, a_R_1 = a_R_1) - 
          fun.prior_a_R(a_R = a_R[i-1], a_R_1 = a_R_1)  
        
      }else{
        
        r_a_R       <- dgamma(x = Z_eta[i], shape = a_R_prop, rate = 1, log = T) -
          dgamma(x = Z_eta[i], shape = a_R[i-1], rate = 1, log = T) + 
          fun.prior_a_R(a_R = a_R_prop, a_R_1 = a_R_1) - 
          fun.prior_a_R(a_R = a_R[i-1], a_R_1 = a_R_1)
      }
      
      
      if(is.na(r_a_R)){
        r_a_R <- -Inf
      }
      
      # Draw an indicator
      u_a_R       <- log(runif(n = 1, min = 0, max = 1))
      
      # Accept/reject
      if((r_a_R > u_a_R)){
        a_R[i]   <- a_R_prop
      }else{
        a_R[i]   <- a_R[i-1]
      } 
      
    }
    
    
  }
  
  # Save warning message
  warnings_message <- as_tibble(warnings_message)
  colnames(warnings_message) <- c("Iteration", "phi_prop", "phi_old", "gamma_prop", "gamma_old")
  
  # Type some information
  cat(green("Borys, I'm done at the" ,Sys.time() ,"Have a nice day","\n"))
  
  # Thining
  indicator_seq <- seq(from = burnin, to = mc, by = thin)
  
  # Delete warnings
  if(isFALSE(save_warnings)){
    warnings_message <- NULL 
  }
  
  # Save the choice of prior for a_R
  Prior$Prior_a_R <- Prior_a_R
  Prior$Prior_b_R <- Prior_b_R
  
  # Return
  return(list(alpha_2              = alpha_2[indicator_seq],
              beta                 = beta[indicator_seq],
              psi                  = psi[indicator_seq],
              
              sig2_x               = sig2_x[indicator_seq],
              sig2_v               = sig2_v[indicator_seq],
              sig2_x_tilde         = sig2_x_tilde[indicator_seq],
              Sigma_0_beta         = Sigma_0_beta[indicator_seq],
              sig2_eta             = sig2_eta[indicator_seq],
              Z_eta                = Z_eta[indicator_seq],
              a_R                  = a_R[indicator_seq],
              b_R                  = b_R[indicator_seq],
              
              gamma                = gamma[indicator_seq],
              phi                  = phi[indicator_seq],
              mu                   = mu[indicator_seq],
              
              mu_2_beta_n          = mu_2_beta_n[indicator_seq],    
              Sigma_2_beta_n       = Sigma_2_beta_n[indicator_seq],
              
              ols                  = ols,
              ols_corrected        = ols_corrected,
              res_bias_corrected   = res_bias_corrected,
              ols_true             = ols_true,
              
              Data                 = Data,
              true_data            = true_data,
              Prior                = Prior,
              mc                   = mc,
              burnin               = burnin,
              thin                 = thin,
              
              acceptance_phi_gamma = acceptance_phi_gamma,
              warnings_message     = warnings_message
  ))
  
}

#### Results ####

# Create the results
fun.create_result <- function(data_in, level, prior_r2_type){
  
  # Safe the full length
  N_full    <- length(data_in)
  N_eff_ols <- N_full
  
  index_ols <- 1e9
  
  # Estimate the ESS
  index_bay     <- f_index(data_in, level)
  
  if(index_bay[1] == 1e9){
    N_eff <- N_full 
  }else{
    N_eff <- N_full - length(index_bay)
  }
  
  # Save the true data
  if(is.null(data_in[[1]]$true_data)){
    true_data <- data_4
  }else{
    true_data <- data_in[[1]]$true_data  
  }
  
  
  
  
  # Create a table for phi
  table_results_phi             <- as.data.frame(matrix(NA_real_, nrow = 3, ncol = 11))
  table_results_phi[,1]         <- c("OLS", "RBE", "BAY")
  colnames(table_results_phi)   <- c("Method", "Prior_phi", "Prior_R2" , "N_EFF" , "MAE", 
                                     "MSE", "Viol", "Exc_1" , "Len", "Bias", "Var")
  
  # Save Effective size and prior type
  table_results_phi[1 ,"N_EFF"]    <- N_eff_ols
  table_results_phi[2 ,"N_EFF"]    <- N_eff_ols
  table_results_phi[3 ,"N_EFF"]    <- N_eff
  
  # Prior for R2
  table_results_phi$Prior_R2        <- prior_r2_type
  
  # Prior for phi
  table_results_phi$Prior_phi       <- data_in[[1]]$Prior$Prior_type
  
  # If beta = 0, then we have FP, If beta = 0.1, then we have FN
  if(true_data$beta == 0){
    
    # Create a table for beta
    table_results_beta            <- as.data.frame(matrix(NA_real_, nrow = 3, ncol = 13))
    table_results_beta[,1]        <- c("OLS", "RBE", "BAY")
    colnames(table_results_beta)  <- c("Method", "Prior_phi", "Prior_R2", "N_EFF", "MAE", "MSE",
                                       "Len", "Bias", "Var","Viol", "SAV_med","SAV_avg", "SAV_sam")
  }else{
    
    # Create a table for beta
    table_results_beta            <- as.data.frame(matrix(NA_real_, nrow = 3, ncol = 14))
    table_results_beta[,1]        <- c("OLS", "RBE", "BAY")
    colnames(table_results_beta)  <- c("Method", "Prior_phi", "Prior_R2", "N_EFF", "MAE", "MSE", 
                                       "Len", "Bias", "Var", "Viol", "FN", "SAV_med", "SAV_avg", "SAV_sam")
  }
  
  # Save Effective size and prior type
  table_results_beta[1 ,"N_EFF"]    <- N_eff_ols
  table_results_beta[2 ,"N_EFF"]    <- N_eff_ols
  table_results_beta[3 ,"N_EFF"]    <- N_eff
  
  # Prior for R2
  table_results_beta$Prior_R2        <- prior_r2_type
  
  # Prior for phi
  table_results_beta$Prior_phi       <- data_in[[1]]$Prior$Prior_type
  
  
  # Create tables for the Bayes Factor
  table.BF              <- as.data.frame(matrix(NA_real_, nrow = N_full, ncol = 5))
  colnames(table.BF)    <- c("BF_med", "BF_avg", "BF_sam", "BF_ols", "BF_rbe")
  
  # Create tables for the Beta 
  table.beta            <- as.data.frame(matrix(NA_real_, nrow = N_full, ncol = 3))
  colnames(table.beta)  <- c("BAY", "OLS", "RBE")
  
  # Create tables for the aR
  table.ar            <- as.data.frame(matrix(NA_real_, nrow = N_full, ncol = 1))
  colnames(table.ar)  <- c("BAY")
  
  ### PHI ###
  
  # MAE error
  table_results_phi[1 ,"MAE"] <- f_mae(res = data_in, object = "x$res_bias_corrected$res$coefficients[2,1]", index = index_ols, true_value = true_data$phi)$mean_error
  table_results_phi[2 ,"MAE"] <- f_mae(res = data_in, object = "x$res_bias_corrected$phi_corrected", index = index_ols, true_value = true_data$phi)$mean_error
  table_results_phi[3 ,"MAE"] <- f_mae(res = data_in, object = "x$phi", index = index_bay, true_value = true_data$phi)$mean_error
  
  # MSE error
  table_results_phi[1 ,"MSE"] <- f_mse(res = data_in, object = "x$res_bias_corrected$res$coefficients[2,1]", index = index_ols, true_value = true_data$phi)$mse 
  table_results_phi[2 ,"MSE"] <- f_mse(res = data_in, object = "x$res_bias_corrected$phi_corrected", index = index_ols, true_value = true_data$phi)$mse
  table_results_phi[3 ,"MSE"] <- f_mse(res = data_in, object = "x$phi", index = index_bay, true_value = true_data$phi)$mse
  
  # Bias
  table_results_phi[1 ,"Bias"] <- f_mse(res = data_in, object = "x$res_bias_corrected$res$coefficients[2,1]", index = index_ols, true_value = true_data$phi)$bias
  table_results_phi[2 ,"Bias"] <- f_mse(res = data_in, object = "x$res_bias_corrected$phi_corrected", index = index_ols, true_value = true_data$phi)$bias
  table_results_phi[3 ,"Bias"] <- f_mse(res = data_in, object = "x$phi", index = index_bay, true_value = true_data$phi)$bias
  
  # Variance
  table_results_phi[1 ,"Var"] <- f_mse(res = data_in, object = "x$res_bias_corrected$res$coefficients[2,1]", index = index_ols, true_value = true_data$phi)$var_est
  table_results_phi[2 ,"Var"] <- f_mse(res = data_in, object = "x$res_bias_corrected$phi_corrected", index = index_ols, true_value = true_data$phi)$var_est
  table_results_phi[3 ,"Var"] <- f_mse(res = data_in, object = "x$phi", index = index_bay, true_value = true_data$phi)$var_est
  
  
  # Adjust results
  table_results_phi[,"MAE"]      <- table_results_phi[,"MAE"] * 100
  table_results_phi[,"MSE"]      <- table_results_phi[,"MSE"] * 1000
  table_results_phi[,"Bias"]     <- table_results_phi[,"Bias"] * 100
  table_results_phi[,"Var"]      <- table_results_phi[,"Var"] * 1000
  
  # Credible/Confidence intervals
  # Save the results
  hpd_bay <- f_hpd(res = data_in, object = "x$phi", index = index_bay, true_value = true_data$phi)
  
  conf_ols  <- f_conf_interval_freq(res        = data_in, 
                                    true_value = true_data$phi, 
                                    index      = index_ols, 
                                    FUN        = f_conf_interval_ols_single_phi, 
                                    type       = "ols")
  
  conf_bias <- f_conf_interval_freq(res        = data_in, 
                                    true_value = true_data$phi, 
                                    index      = index_ols, 
                                    FUN        = f_conf_interval_bias_cor_single_phi, 
                                    type       = "ols")
  
  # Save the number of violations
  table_results_phi[1 , "Viol"] <- sum(conf_ols$k_above)  +  sum(conf_ols$k_below)
  table_results_phi[2 , "Viol"] <- sum(conf_bias$k_above) +  sum(conf_bias$k_below)
  table_results_phi[3 , "Viol"] <- sum(hpd_bay$k_above)   +  sum(hpd_bay$k_below)
  
  table_results_phi[, "Viol"] <- 100 * table_results_phi[, "Viol"] / table_results_phi[, "N_EFF"]
  
  # Save the length of the interval
  table_results_phi[1 , "Len"] <- conf_ols$conf_length
  table_results_phi[2 , "Len"] <- conf_bias$conf_length
  table_results_phi[3 , "Len"] <- hpd_bay$hpd_length
  
  # Save the number of credible intervals that exceed 1
  table_results_phi[1 , "Exc_1"] <- sum(conf_ols$conf_inter[,"upper"] > 1)
  table_results_phi[2 , "Exc_1"] <- sum(conf_bias$conf_inter[,"upper"] > 1)
  table_results_phi[3 , "Exc_1"] <- 0
  
  
  
  ### BETA ###
  
  # MAE error
  table_results_beta[1 , "MAE"] <- f_mae(res = data_in, object = "x$ols$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$mean_error
  table_results_beta[2 , "MAE"] <- f_mae(res = data_in, object = "x$ols_corrected$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$mean_error
  table_results_beta[3 , "MAE"] <- f_mae(res = data_in, object = "x$beta", index = index_bay, true_value = true_data$beta)$mean_error
  
  
  # MSE error
  table_results_beta[1 , "MSE"] <- f_mse(res = data_in, object = "x$ols$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$mse
  table_results_beta[2 , "MSE"] <- f_mse(res = data_in, object = "x$ols_corrected$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$mse
  table_results_beta[3 , "MSE"] <- f_mse(res = data_in, object = "x$beta", index = index_bay, true_value = true_data$beta)$mse
  
  
  # Bias
  table_results_beta[1 , "Bias"] <- f_mse(res = data_in, object = "x$ols$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$bias
  table_results_beta[2 , "Bias"] <- f_mse(res = data_in, object = "x$ols_corrected$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$bias
  table_results_beta[3 , "Bias"] <- f_mse(res = data_in, object = "x$beta", index = index_bay, true_value = true_data$beta)$bias
  
  
  # MSE error
  table_results_beta[1 , "Var"] <- f_mse(res = data_in, object = "x$ols$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$var_est
  table_results_beta[2 , "Var"] <- f_mse(res = data_in, object = "x$ols_corrected$coefficients[2,1]", index = index_ols, true_value = true_data$beta)$var_est
  table_results_beta[3 , "Var"] <- f_mse(res = data_in, object = "x$beta", index = index_bay, true_value = true_data$beta)$var_est
  
  # Adjust results
  table_results_beta[,"MAE"]     <- table_results_beta[,"MAE"]  * 100
  table_results_beta[,"MSE"]     <- table_results_beta[,"MSE"]  * 1000
  table_results_beta[,"Bias"]    <- table_results_beta[,"Bias"] * 100
  table_results_beta[,"Var"]     <- table_results_beta[,"Var"]  * 1000
  
  # Credible/Confidence intervals
  # Save the results
  hpd_bay <- f_hpd(res = data_in, object = "x$beta", index = index_bay, true_value = true_data$beta)
  
  conf_ols  <- f_conf_interval_freq(res        = data_in, 
                                    true_value = true_data$beta, 
                                    index      = index_ols, 
                                    FUN        = f_conf_interval_ols_single_beta,
                                    type       = "ols")
  
  conf_bias <- f_conf_interval_freq(res        = data_in, 
                                    true_value = true_data$beta, 
                                    index      = index_ols, 
                                    FUN        = f_conf_interval_bias_cor_single_beta,
                                    type       = "ols")
  
  # Save the number of violations
  table_results_beta[1 , "Viol"] <- sum(conf_ols$k_above)  +  sum(conf_ols$k_below)
  table_results_beta[2 , "Viol"] <- sum(conf_bias$k_above) +  sum(conf_bias$k_below)
  table_results_beta[3 , "Viol"] <- sum(hpd_bay$k_above)   +  sum(hpd_bay$k_below)
  
  table_results_beta[, "Viol"] <- 100 * table_results_beta[, "Viol"] / table_results_beta[, "N_EFF"]
  
  # Save the length of the interval
  table_results_beta[1 , "Len"] <- conf_ols$conf_length
  table_results_beta[2 , "Len"] <- conf_bias$conf_length
  table_results_beta[3 , "Len"] <- hpd_bay$hpd_length
  
  # If beta != 0 then we need to fill a column FN
  if(true_data$beta != 0){
    
    table_results_beta[1 , "FN"] <- sum(conf_ols$conf_inter[,"lower"] < 0)
    table_results_beta[2 , "FN"] <- sum(conf_bias$conf_inter[,"lower"] < 0)
    table_results_beta[3 , "FN"] <- sum(hpd_bay$hpd[,"lower"] < 0)
    
    
    table_results_beta[, "FN"] <- 100 * table_results_beta[, "FN"] / table_results_beta[, "N_EFF"]
    
  }
  
  if(prior_r2_type == "sig2" || prior_r2_type == "sig2_x"){
    
    prior_beta <- dnorm(x = 0, mean = 0, sd = sqrt(data_in[[1]]$Prior$Sigma_0_beta))
    
    # The same for all since here we have no difference
    prior_beta_med <- prior_beta
    
    prior_beta_avg <- prior_beta
    
    prior_beta_sam <- prior_beta
    
    posterior_beta <- fun.sim_posterior_beta(res = data_in, index = index_bay)
    
    print("code is not adjusted in posterior beta")
    
  }else if(prior_r2_type == "r2"){
    
    # In the case of ols we have for every dataset different prior
    if(data_in[[1]]$Prior$b_R == "ols" || data_in[[1]]$Prior$b_R == "fixed"){
      
      # We apply function to every element of list
      prior_beta         <- sapply(data_in, function(x) fun.sim_prior_beta(T_len = 1e6, res = x))
      
      # Delete non-needed elements for the median results
      prior_beta_med     <- unlist(prior_beta["beta_est_0_med",])[-index_bay]
      
      # Delete non-needed elements for the mean results
      prior_beta_avg     <- unlist(prior_beta["beta_est_0_avg",])[-index_bay]
      
      # Delete non-needed elements for the sample results
      prior_beta_sam     <- unlist(prior_beta["beta_est_0_sam",])[-index_bay]
      
      
    }else{
      
      prior_beta         <- fun.sim_prior_beta(T_len = 1e6, res = data_in[[1]])
      
      # Median results
      prior_beta_med     <- prior_beta$beta_est_0_med
      
      # Average results
      prior_beta_avg     <- prior_beta$beta_est_0_avg
      
      # Sample results
      prior_beta_sam     <- prior_beta$beta_est_0_sam
      
    }
    
    posterior_beta <- fun.sim_posterior_beta_given_b_t(res = data_in)
    
    
  }else{
    print("Please specify the correct prior type")
  }
  
  # Ratio for median
  K_med <- posterior_beta$post_dens_med / prior_beta_med  
  
  # Ratio for mean
  K_avg <- posterior_beta$post_dens_avg / prior_beta_avg  
  
  # Ratio for sampling
  K_sam <- posterior_beta$post_dens_sam / prior_beta_sam 
  
  # If beta == 0
  if(true_data$beta == 0){
    
    
    # We need the Savage ratio to be LARGER than 1 to call it False Positive.
    # In this case our Prior is larger than our Posterior and hence we impose less mass to 0 after observing the data
    table_results_beta[1 , "SAV_med"] <- table_results_beta[1, "Viol"]
    table_results_beta[2 , "SAV_med"] <- table_results_beta[2, "Viol"]
    table_results_beta[3 , "SAV_med"] <- 100 * mean(K_med[-index_bay] < 1)
    
    table_results_beta[1 , "SAV_avg"] <- table_results_beta[1, "Viol"]
    table_results_beta[2 , "SAV_avg"] <- table_results_beta[2, "Viol"]
    table_results_beta[3 , "SAV_avg"] <- 100 * mean(K_avg[-index_bay] < 1)
    
    table_results_beta[1 , "SAV_sam"] <- table_results_beta[1, "Viol"]
    table_results_beta[2 , "SAV_sam"] <- table_results_beta[2, "Viol"]
    table_results_beta[3 , "SAV_sam"] <- 100 * mean(K_sam[-index_bay] < 1)
    
  }else{
    
    # We need the Savage ratio to be SMALLER than 1 to call it False Negative
    # In this case our Prior is smaller than our Posterior and hence we impose more mass to 0 after observing the data
    table_results_beta[1 , "SAV_med"] <- table_results_beta[1, "FN"]
    table_results_beta[2 , "SAV_med"] <- table_results_beta[2, "FN"]
    table_results_beta[3 , "SAV_med"] <- 100 * mean(K_med[-index_bay] > 1)
    
    table_results_beta[1 , "SAV_avg"] <- table_results_beta[1, "FN"]
    table_results_beta[2 , "SAV_avg"] <- table_results_beta[2, "FN"]
    table_results_beta[3 , "SAV_avg"] <- 100 * mean(K_avg[-index_bay] > 1)
    
    table_results_beta[1 , "SAV_sam"] <- table_results_beta[1, "FN"]
    table_results_beta[2 , "SAV_sam"] <- table_results_beta[2, "FN"]
    table_results_beta[3 , "SAV_sam"] <- 100 * mean(K_sam[-index_bay] > 1)
    
  }
  
  # Store BF
  table.BF$BF_med  <- K_med
  table.BF$BF_avg  <- K_avg
  table.BF$BF_sam  <- K_sam
  table.BF$BF_ols  <- sapply(data_in, function(x) x$ols$coefficients[2,"t value"])
  table.BF$BF_rbe  <- sapply(data_in, function(x) fun.extract_p_value_rbe(x)$rbe.t)
  
  # Store Betas
  table.beta$BAY     <- sapply(data_in, function(x) mean(x$beta))
  table.beta$OLS     <- sapply(data_in, function(x) x$ols$coefficients[2,1])
  table.beta$RBE     <- sapply(data_in, function(x) x$ols_corrected$coefficients[2,1])
  
  # Store aRs
  if(length(data_in[[1]]$Prior$a_R) > 1){
    table.ar$BAY       <- sapply(data_in, function(x) mean(x$a_R))  
  }
  
  
  
  return(res = list(res_phi       = as_tibble(table_results_phi), 
                    res_beta      = as_tibble(table_results_beta), 
                    table.beta    = as_tibble(table.beta), 
                    table.ar      = as_tibble(table.ar), 
                    table.BF      = as_tibble(table.BF),
                    index_bay     = index_bay,
                    level         = level,
                    Prior         = data_in[[1]]$Prior))
  
}


# OOS analysis
fun.oos_estimation <- function(data_in, type, i, window){
  
  # Take a subset of the data
  if(type == "exp"){
    data_in_sub <- data_in[1:i,]
    
  }else if(type == "rol"){
    data_in_sub <- data_in[(i - window_size):i,]
    
  }else{
    return("Specify the correct type")
  }
  
  # Prepare the data
  data_in_sub <- list(x_latent   = data_in_sub$x_in[-nrow(data_in_sub)],
                      y_latent   = data_in_sub$x_in[-1],
                      x_observe  = data_in_sub$x_in[-nrow(data_in_sub)],
                      y_observe  = data_in_sub$ret.crsp.l[-1],
                      dg_observe = NULL)
  
  
  # Estimate the model
  res.bayes <- f_Bayesian_control_function_Cholesky_AR_prior_R2(Data             = data_in_sub, 
                                                                Prior            = Prior_cont_ar_in,
                                                                mc               = MC_draws,
                                                                burnin           = burnin,
                                                                start_true_ind   = F,
                                                                thin             = thining,
                                                                true_data        = NULL)
  
  # Calculate the Returns
  ret_bay      <- res.bayes$alpha_2 + res.bayes$beta * last(data_in_sub$y_latent)
  ret_bay_mean <- mean(ret_bay)
  
  ret_rbe <- res.bayes$ols_corrected$coefficients[1,1] + 
    res.bayes$ols_corrected$coefficients[2,1]* last(data_in_sub$y_latent)
  
  ret_ols <- res.bayes$ols$coefficients[1,1] + 
    res.bayes$ols$coefficients[2,1]* last(data_in_sub$y_latent)
  
  # We need data_in since it contains future returns as well.
  ret_fut <- data_in$ret.crsp.l[i+1]
  
  # Calculate the Errors
  errors    <- {}
  errors[1] <- i
  errors[2] <- ret_fut
  errors[3] <- ret_fut - ret_bay_mean
  errors[4] <- ret_fut - ret_rbe
  errors[5] <- ret_fut - ret_ols
  errors[6] <- ret_fut - mean(data_in_sub$y_observe)
  
  return(list(errors = errors, res.bayes = res.bayes))
  
}

# MAE
f_mae <- function(res, object, index, true_value){
  
  object_est <- sapply(res, function(x) mean(eval(parse(text=paste0(object)))))
  
  res <- abs(object_est - true_value)[-index]
  
  return(list(errors = res, mean_error = mean(res)))
}

# MSE
f_mse <- function(res, object, index, true_value){
  
  # Save elements
  object_est <- sapply(res, function(x) mean(eval(parse(text=paste0(object)))))
  
  # Delete non-convergent MCMC
  object_est <- object_est[-index]
  
  # Calculate squared error
  res        <- (object_est - true_value)^2
  
  # Calculate bias squared
  bias       <- (mean(object_est) - true_value)
  
  # Calculate the variance of the estimator
  var_est    <- var(object_est)
  
  return(list(errors = res, mse = mean(res), bias = bias, var_est = var_est))
}

# We specify the rule when we exclude the "bad samples".
f_index_rule <- function(res){
  
  index <- sapply(res, function(x) effectiveSize(x$beta))
  
  return(index)
  
}

# Using the f_index_rule we find all number of datasets which have to excluded
f_index <- function(data_in, level){
  
  # List the data
  results <- list(data_in)
  
  res_index <- sapply(results, function(x) f_index_rule(x))
  
  index_mix <- unlist(apply(res_index, 2, function(x) which(x < level)))
  
  index_mix <- unique(index_mix)
  
  # Estimate the Convergence
  Z_score       <- sapply(data_in, function(x) geweke.diag(x = x$beta)$z)
  index_conv    <- which(Z_score < qnorm(0.025) | Z_score > qnorm(0.975))
  
  index <- sort(union(index_mix, index_conv))
  
  if(sum(index) == 0) index <- 1e9
  
  return(index)
}

# k = TRUE, means if violated
f_hpd <- function(res, object, index = NULL, true_value){
  hdi_mcmc <- t(sapply(res, function(x) hdi(eval(parse(text=paste0(object))))))
  mean     <- sapply(res, function(x) mean(eval(parse(text=paste0(object)))))
  N        <- length(mean)
  
  if(!is.null(index)){
    hdi_mcmc   <- hdi_mcmc[-index,,drop = FALSE]
    mean       <- mean[-index]
    N          <- length(mean)
  }
  
  hpd_length <- round(mean(hdi_mcmc[,2] - hdi_mcmc[,1]),3)
  
  k_above  = true_value > hdi_mcmc[,2]
  k_below  = true_value < hdi_mcmc[,1]
  
  return(list(hpd = hdi_mcmc, mean = mean, k_above = k_above, k_below = k_below , N = N, hpd_length = hpd_length))
}

# Given OLS objects compute the confidence intervals for phi
f_conf_interval_ols_single_phi <- function(object, type){
  
  object <- object$res_bias_corrected$res
  
  mean_in <- object$coef[2,1]
  
  #Compute the total number of observations
  N       <- object$df[1] + object$df[2]
  
  if(type == "ols"){
    conf_inter    <- c("lower" = mean_in - qnorm(0.975) * object$coef[2, 2],
                       "upper" = mean_in + qnorm(0.975) * object$coef[2, 2])   
  }else if(type == "t"){
    conf_inter    <- c("lower" = mean_in - qt(0.975, df = object$df[2]) * object$coef[2, 2],
                       "upper" = mean_in + qt(0.975, df = object$df[2]) * object$coef[2, 2])
    
  }else if(type == "asy"){
    if(mean_in < 1){
      conf_inter   <- c("lower" = mean_in - qnorm(0.975) * sqrt((1 - mean_in^2)/N),
                        "upper" = mean_in + qnorm(0.975) * sqrt((1 - mean_in^2)/N))   
    }else{
      # We take the ols result and the conditional on phi_est < 1 we change to asymptotic
      conf_inter             <- c("lower" = mean_in - qnorm(0.975) * object$coef[2, 2],
                                  "upper" = mean_in + qnorm(0.975) * object$coef[2, 2])   
    }
    
  }else{
    print("Wrong name")
  }
  
  return(list(mean = mean_in, conf = conf_inter))
}

# Given bias corrected objects compute the confidence intervals for phi
f_conf_interval_bias_cor_single_phi <- function(object, type){  
  # Save the object
  object  <- object$res_bias_corrected
  
  # Take the mean
  mean_in <- object$phi_corrected
  
  # Compute the total number of observations
  N       <- object$res$df[1] + object$res$df[2]
  
  # Estimate the se
  se <- object$res$coef[2, 2]*(1 + 3/N + 9/(N^2) )
  
  if(type == "ols"){
    conf_inter <- c("lower" = mean_in - qnorm(0.975) * se,
                    "upper" = mean_in + qnorm(0.975) * se)
    
  }else if(type == "t"){
    conf_inter <- c("lower" = mean_in - qt(0.975, df = object$res$df[2]) * se,
                    "upper" = mean_in + qt(0.975, df = object$res$df[2]) * se)
    
  }else if(type == "asy"){
    if(mean_in < 1){
      conf_inter    <- c("lower" = mean_in - qnorm(0.975) * sqrt((1 - mean_in^2)/N)*(1 + 3/N + 9/N^2),
                         "upper" = mean_in + qnorm(0.975) * sqrt((1 - mean_in^2)/N)*(1 + 3/N + 9/N^2)) 
    }else{
      # We take the ols result and the conditional on phi_est < 1 we change to asymptotic
      conf_inter <- c("lower" = mean_in - qnorm(0.975) * se,
                      "upper" = mean_in + qnorm(0.975) * se)  
    }
    
  }else{
    print("Wrong name")
  }
  
  
  
  return(list(mean = mean_in, conf = conf_inter, se = se))
}

# Given OLS objects compute the confidence intervals for beta
f_conf_interval_ols_single_beta <- function(object, type, alpha_in = 0.05){
  
  object <- object$ols
  
  mean_in <- object$coef[2,1]
  conf_inter    <- c("lower" = mean_in - qnorm(1 - alpha_in/2) * object$coef[2, 2],
                     "upper" = mean_in + qnorm(1 - alpha_in/2) * object$coef[2, 2])
  
  return(list(mean = mean_in, conf = conf_inter))
}

# Given bias corrected objects compute the confidence intervals for beta
f_conf_interval_bias_cor_single_beta <- function(object, type, alpha_in = 0.05){
  
  # Take the number of observations
  N       <- object$ols_corrected$df[1] + object$ols_corrected$df[2]
  
  # Compute the mean
  mean_in <- object$ols_corrected$coefficients[2,1]
  
  # Compute the standard errors
  hat_var_rho <- ( object$res_bias_corrected$res$coefficients[2,2]*( 1 + 3/N + 9/(N^2) ) )^2
  hat_psi     <- object$ols_corrected$coefficients[3, 1] 
  se          <- sqrt(hat_var_rho * hat_psi^2 + object$ols_corrected$coefficients[2, 2]^2)
  
  conf_inter <- c("lower" = mean_in - qnorm(1 - alpha_in/2) * se,
                  "upper" = mean_in + qnorm(1 - alpha_in/2) * se) 
  
  return(list(mean = mean_in, conf = conf_inter, se = se))
}

# Computes the confidence intervals for the MCMC datasets
f_conf_interval_freq <- function(res, true_value, index = NULL, FUN, type){
  
  conf_inter <- t(sapply(res, function(x) FUN(x, type = type)$conf))
  mean       <- sapply(res, function(x) FUN(x, type = type)$mean)
  N          <- length(res) 
  
  if(!is.null(index)){
    conf_inter <- conf_inter[-index,]
    mean       <- mean[-index]
    N          <- length(mean)
  }
  
  conf_length <- round(mean(conf_inter[,2] - conf_inter[,1]),3)
  
  k_above  = true_value > conf_inter[,2]
  k_below  = true_value < conf_inter[,1]
  
  return(list(conf_inter = conf_inter, mean = mean, k_above = k_above, k_below = k_below, N = N, conf_length = conf_length))
}

f_conf_interval_ols_plot <- function(res, true_value, plot_name, index = NULL, FUN, type){
  
  conf_inter <- t(sapply(res, function(x) FUN(x, type = type)$conf))
  mean       <- sapply(res, function(x) FUN(x, type = type)$mean)
  N          <- length(res)  
  
  if(!is.null(index)){
    conf_inter <- conf_inter[-index,]
    mean       <- mean[-index]
    N          <- length(mean)
  }
  
  conf_length <- round(mean(conf_inter[,2] - conf_inter[,1]),3)
  
  k = (true_value < conf_inter[,1])|(true_value > conf_inter[,2])
  print(paste0(sum(k), " Empirical number of interval which don't contain the true value "))
  print(paste0(N*0.05, " Theoretical number of interval which don't contain the true value "))
  
  # Plot means n_data times, we choose suitable range according the conf interval
  plot(mean, 1:N, xlim = range(conf_inter), axes = FALSE,
       xlab = "", ylab = "", pch = 19, cex = .7,
       col = c("blue", "red")[1+k], 
       main = paste("Var:", plot_name, ", Viol:", sum(k), ", N points:", N, ", Avg length:", conf_length))
  axis(1)
  segments(conf_inter[,1],1:N,conf_inter[,2],1:N,col=c("blue","red")[1+k])
  abline(v=true_value)
}

f_conf_interval_bias_cor_plot <- function(res, true_value, plot_name, index = NULL, FUN, type){
  
  conf_inter <- t(sapply(res, function(x) FUN(x, type = type)$conf))
  mean       <- sapply(res, function(x) FUN(x, type = type)$mean)
  N          <- length(res) 
  
  if(!is.null(index)){
    conf_inter <- conf_inter[-index,]
    mean       <- mean[-index]
    N          <- length(mean)
  }
  
  conf_length <- round(mean(conf_inter[,2] - conf_inter[,1]),3)
  
  k = (true_value < conf_inter[,1])|(true_value > conf_inter[,2])
  print(paste0(sum(k), " Empirical number of interval which don't contain the true value "))
  print(paste0(N*0.05, " Theoretical number of interval which don't contain the true value "))
  
  # Plot means n_data times, we choose suitable range according the conf interval
  plot(mean, 1:N, xlim = range(conf_inter), axes = FALSE,
       xlab = "", ylab = "", pch = 19, cex = .7,
       col = c("blue", "red")[1+k], 
       main = paste("Var:", plot_name, ", Viol:", sum(k), ", N points:", N, ", Avg length:", conf_length))
  axis(1)
  segments(conf_inter[,1],1:N,conf_inter[,2],1:N,col=c("blue","red")[1+k])
  abline(v=true_value)
}

# Extract R2 for RBE
fun.extract_r2_rbe <- function(res){
  
  y_obs          <- res$Data$y_observe
  x_obs          <- res$Data$x_observe
  res_bias_error <- res$res_bias_corrected$error_corrected
  ols_full       <- res$ols_corrected
  
  ols_reduced    <- summary(lm(y_obs ~ res_bias_error))
  
  r2_rbe <- 1 - sum(ols_full$residuals^2) / sum(ols_reduced$residuals^2)
  
  return(r2_rbe)
  
}

# Extract p-value for RBE
fun.extract_p_value_rbe <- function(res){
  
  rbe.beta    <- res$ols_corrected$coefficients[2,1]
  rbe.se      <- f_conf_interval_bias_cor_single_beta(res)$se
  rbe.t       <- rbe.beta / rbe.se
  T_len       <- length(res$Data$y_observe)
  
  rbe.p_value <- pt(q          = abs(rbe.t), 
                    df         = T_len, 
                    lower.tail = F) +
    pt(q          = -abs(rbe.t), 
       df         = T_len)
  
  return(list(rbe.t       = rbe.t,
              rbe.p_value = rbe.p_value))
  
}

# ESS summary for the individual variable
fun.effective_summary <- function(var, res, names_in){
  
  # Length
  N_row <- length(names_in)
  
  # Save the variables
  var <- paste0("x$", var)
  
  # Create Table
  table_result <- tibble(Sample = names_in, 
                         Eff    = rep(NA_real_, N_row),
                         Length = rep(NA_real_, N_row))
  
  
  # Calculate the results
  table_result[,2] <- apply(sapply(res, function(x) eval(parse(text=paste0(var)))), 2, effectiveSize)
  table_result[,3] <- apply(sapply(res, function(x) eval(parse(text=paste0(var)))), 2, length)
  
  return(table_result)
  
}

# ESS summary for all variables
fun.effective_summary_table <- function(res_in, names_in){
  
  N_row <- length(names_in)
  
  # Create Table
  table_result   <- as.data.frame(matrix(NA_real_, nrow = N_row*2, ncol = 5))
  
  # Names
  table_result[,1] <- rep(names_in,2)
  
  # Number of observations
  table_result[1:N_row,5] <- fun.effective_summary(res      = res_in, 
                                                   var      = "phi",
                                                   names_in = names_in)[,3]
  
  # ESS
  table_result[1:N_row,2] <- fun.effective_summary(res      = res_in, 
                                                   var      = "psi",
                                                   names_in = names_in)[,2]
  table_result[1:N_row,3] <- fun.effective_summary(res      = res_in, 
                                                   var      = "phi",
                                                   names_in = names_in)[,2]
  table_result[1:N_row,4] <- fun.effective_summary(res      = res_in, 
                                                   var      = "beta",
                                                   names_in = names_in)[,2]
  
  # Convergence
  
  Z_score_psi   <- sapply(res_in, function(x) geweke.diag(x = x$psi)$z)
  Z_score_phi   <- sapply(res_in, function(x) geweke.diag(x = x$phi)$z)
  Z_score_beta  <- sapply(res_in, function(x) geweke.diag(x = x$beta)$z)
  
  table_result[(N_row+1):(N_row*2),2] <- Z_score_psi
  table_result[(N_row+1):(N_row*2),3] <- Z_score_phi
  table_result[(N_row+1):(N_row*2),4] <- Z_score_beta
  
  colnames(table_result) <- c("Sample", "Psi", "Phi", "Beta", "Length")
  
  return(as_tibble(table_result))
  
}

# Summary for the posterior density results for beta
fun.posterior_results_beta <- function(data_in, prior_type){
  
  if(prior_type == "sig2" || prior_type == "sig2_x"){
    
    prior_beta <- dnorm(x    = 0, 
                        mean = 0, 
                        sd   = sqrt(data_in$Prior$Sigma_0_beta))
    
    posterior_beta <- fun.sim_posterior_beta_given_b_t(res   = list(data_in), 
                                                       index = NULL)
    
    return(list(K_med      = posterior_beta$post_dens_med / prior_beta,
                K_avg      = posterior_beta$post_dens_avg / prior_beta,
                K_sam      = posterior_beta$post_dens_sam / prior_beta))  
    
  }else if(prior_type == "r2"){
    
    prior_beta     <- fun.sim_prior_beta(T_len = 1e6, 
                                         res   = data_in)
    
    posterior_beta <- fun.sim_posterior_beta_given_b_t(res   = list(data_in))
    
    return(list(K_med      = posterior_beta$post_dens_med / prior_beta$beta_est_0_med,
                K_avg      = posterior_beta$post_dens_avg / prior_beta$beta_est_0_avg,
                K_sam      = posterior_beta$post_dens_sam / prior_beta$beta_est_0_sam,
                prior_beta = prior_beta$beta_sample))  
    
  }  
  
  
}


# Combine the results
fun.combine_results <- function(folder){
  
  files_dgp_1 <- list.files(paste0(path,"data/1.sim_data/", folder,"/dgp_1"))
  files_dgp_2 <- list.files(paste0(path,"data/1.sim_data/", folder,"/dgp_2"))
  
  res_all_phi_dgp_1  <- {}
  res_all_phi_dgp_2  <- {}
  res_all_beta_dgp_1 <- {}
  res_all_beta_dgp_2 <- {}
  
  split.dgp_1 <- lapply(files_dgp_1, function(x) unlist(strsplit(x, split = "[.]")))
  split.dgp_2 <- lapply(files_dgp_2, function(x) unlist(strsplit(x, split = "[.]")))
  
  names.col_dgp_1 <- sapply(split.dgp_1, function(x) paste(x[2], x[3], sep = "_"))
  names.col_dgp_2 <- sapply(split.dgp_2, function(x) paste(x[2], x[3], sep = "_"))
  
  max_n           <- max(sapply(split.dgp_1, function(x) as.double(x[3])), na.rm = T)
  
  data.BF_1 <- as_tibble(matrix(NA_real_, nrow = 1e3, ncol = length(files_dgp_1) + 2))
  colnames(data.BF_1)   <- c("OLS", "RBE", paste0("BAY_", names.col_dgp_1))
  
  data.BF_2 <- as_tibble(matrix(NA_real_, nrow = 1e3, ncol = length(files_dgp_2) + 2))
  colnames(data.BF_2)   <- c("OLS", "RBE", paste0("BAY_", names.col_dgp_2))
  
  # Beta
  data.beta_1 <- as_tibble(matrix(NA_real_, nrow = 1e3, ncol = length(files_dgp_1) + 2))
  colnames(data.beta_1)   <- c("OLS", "RBE", paste0("BAY_", names.col_dgp_1))
  data.beta_2 <- as_tibble(matrix(NA_real_, nrow = 1e3, ncol = length(files_dgp_2) + 2))
  colnames(data.beta_2)   <- c("OLS", "RBE", paste0("BAY_", names.col_dgp_2))
  
  # a_R
  data.a_R_1 <- as_tibble(matrix(NA_real_, nrow = 1e3, ncol = length(files_dgp_1) + 2))
  colnames(data.a_R_1)   <- c("OLS", "RBE", paste0("BAY_", names.col_dgp_1))
  data.a_R_2 <- as_tibble(matrix(NA_real_, nrow = 1e3, ncol = length(files_dgp_2) + 2))
  colnames(data.a_R_2)   <- c("OLS", "RBE", paste0("BAY_", names.col_dgp_2))
  
  
  for(i in 1:length(files_dgp_1)){
    
    # Load the data
    data_1     <- readRDS(file = paste0(path,"data/1.sim_data/", folder,"/dgp_1/", files_dgp_1[i]))
    
    # a_R could be fixed or random
    if(isTRUE(data_1$Prior$Prior_a_R)){
      
      a_R_in <- paste(data_1$Prior$a_R[1], data_1$Prior$a_R[2])
      
    }else{
      
      a_R_in <- data_1$Prior$a_R
      
    }
    
    res.phi_1  <- cbind(paste0("BAY_", names.col_dgp_1)[i], a_R_in, data_1$Prior$b_R, data_1$res_phi)
    
    res.beta_1 <- cbind(paste0("BAY_", names.col_dgp_1)[i], a_R_in, data_1$Prior$b_R, data_1$res_beta)
    
    res_all_phi_dgp_1  <- rbind(res_all_phi_dgp_1, res.phi_1)
    
    res_all_beta_dgp_1 <- rbind(res_all_beta_dgp_1, res.beta_1)
    
    data.BF_1[,i+2]    <- data_1$table.BF$BF_med
    
    
    if(!is.null(data_1$table.beta)){
      data.beta_1[,i+2]       <- data_1$table.beta$BAY
      data.a_R_1[,i+2]        <- data_1$table.ar$BAY
    }
    
  }
  
  for(i in 1:length(files_dgp_2)){
    
    data_2     <- readRDS(file = paste0(path,"data/1.sim_data/", folder,"/dgp_2/", files_dgp_2[i]))
    
    # a_R could be fixed or random
    if(isTRUE(data_2$Prior$Prior_a_R)){
      
      a_R_in <- paste(data_2$Prior$a_R[1], data_2$Prior$a_R[2])
      
    }else{
      
      a_R_in <- data_2$Prior$a_R
      
    }
    
    res.phi_2  <- cbind(paste0("BAY_", names.col_dgp_2)[i], a_R_in, data_2$Prior$b_R, data_2$res_phi)
    
    res.beta_2 <- cbind(paste0("BAY_", names.col_dgp_2)[i], a_R_in, data_2$Prior$b_R, data_2$res_beta)
    
    res_all_phi_dgp_2  <- rbind(res_all_phi_dgp_2, res.phi_2)
    
    res_all_beta_dgp_2 <- rbind(res_all_beta_dgp_2, res.beta_2)
    
    data.BF_2[,i+2]    <- data_2$table.BF$BF_med
    
    
    if(!is.null(data_2$table.beta)){
      data.beta_2[,i+2]       <- data_2$table.beta$BAY
      data.a_R_2[,i+2]        <- data_2$table.ar$BAY
    }
    
    
  }
  
  # Add OLS
  data.BF_1[,"OLS"]     <- data_1$table.BF$BF_ols
  data.BF_2[,"OLS"]     <- data_2$table.BF$BF_ols
  data.beta_1[,"OLS"]   <- data_1$table.beta$OLS
  data.beta_2[,"OLS"]   <- data_2$table.beta$OLS
  
  # Add RBE
  data.BF_1[,"RBE"] <- data_1$table.BF$BF_rbe
  data.BF_2[,"RBE"] <- data_2$table.BF$BF_rbe
  data.beta_1[,"RBE"]   <- data_1$table.beta$RBE
  data.beta_2[,"RBE"]   <- data_2$table.beta$RBE
  
  # Order in columns
  if(isTRUE(data_1$Prior$Prior_a_R)){
    
    data.BF_1 <- data.BF_1 %>%
      relocate(any_of(paste0("BAY_5_", 10:max_n)), .after = last_col())
    data.BF_2 <- data.BF_2 %>%
      relocate(any_of(paste0("BAY_5_", 10:max_n)), .after = last_col())
    
  }else{
    
    data.BF_1 <- data.BF_1 %>%
      relocate(any_of(paste0("BAY_4_", 10:max_n)), .after = last_col())
    data.BF_2 <- data.BF_2 %>%
      relocate(any_of(paste0("BAY_4_", 10:max_n)), .after = last_col())
    
  }
  
  # Change names
  colnames(res_all_phi_dgp_1)[1]  <- "Prior_type"
  colnames(res_all_phi_dgp_2)[1]  <- "Prior_type"
  colnames(res_all_beta_dgp_1)[1] <- "Prior_type"
  colnames(res_all_beta_dgp_2)[1] <- "Prior_type"
  
  # Delete/Change for a_R and b_R 
  colnames(res_all_phi_dgp_1)[2]  <- "a_R"
  colnames(res_all_phi_dgp_2)[2]  <- "a_R"
  colnames(res_all_beta_dgp_1)[2] <- "a_R"
  colnames(res_all_beta_dgp_2)[2] <- "a_R"
  
  colnames(res_all_phi_dgp_1)[3]  <- "b_R"
  colnames(res_all_phi_dgp_2)[3]  <- "b_R"
  colnames(res_all_beta_dgp_1)[3] <- "b_R"
  colnames(res_all_beta_dgp_2)[3] <- "b_R"  
  
  #### Beta #### 
  beta_all <- res_all_beta_dgp_1 %>% 
    inner_join(select(res_all_beta_dgp_2, -Prior_R2, -a_R, -b_R, -Prior_phi), by = c("Prior_type", "Method"), suffix = c(".dgp_1", ".dgp_2"))
  
  beta_all_methods <- beta_all %>% 
    filter(!(Method %in% c("OLS", "RBE")))
  
  beta_all_OLS_RBE_avg <- beta_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi" , "Prior_R2"), ~first(.x)),
              across(N_EFF.dgp_1:SAV_sam.dgp_2, ~ mean(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  beta_all_OLS_RBE_max <- beta_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi" , "Prior_R2"), ~first(.x)),
              across(N_EFF.dgp_1:SAV_sam.dgp_2, ~ max(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  beta_all_methods <- rbind(beta_all_OLS_RBE_avg,  beta_all_methods)
  
  #### PHI #### 
  
  phi_all <- res_all_phi_dgp_1 %>% 
    inner_join(select(res_all_phi_dgp_2, -Prior_R2, -a_R, -b_R, -Prior_phi), by = c("Prior_type", "Method"), suffix = c(".dgp_1", ".dgp_2"))
  
  phi_all_methods <- phi_all %>% 
    filter(!(Method %in% c("OLS", "RBE")))
  
  phi_all_OLS_RBE_avg <- phi_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi",  "Prior_R2"), ~first(.x)),
              across(N_EFF.dgp_1:Var.dgp_2, ~ mean(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  phi_all_OLS_RBE_max <- phi_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi",  "Prior_R2"), ~first(.x)),
              across(N_EFF.dgp_1:Var.dgp_2, ~ max(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  phi_all_methods <- rbind(phi_all_OLS_RBE_avg,  phi_all_methods)
  
  
  return(list(phi_all_methods  = phi_all_methods,
              beta_all_methods = beta_all_methods,
              data.BF_1        = data.BF_1,
              data.BF_2        = data.BF_2,
              data.beta_1      = data.beta_1,
              data.beta_2      = data.beta_2,
              data.a_R_1       = data.a_R_1,
              data.a_R_2       = data.a_R_2
  ))
}

# Combine the results V2
fun.combine_results_v2 <- function(folder){
  
  # Create the path
  path_in    <- paste0(path,"data/1.sim_data/", folder)
  
  files_dgp <- list.files(path_in)
  
  res_all_phi_dgp  <- {}
  res_all_beta_dgp <- {}
  
  split.dgp <- lapply(files_dgp, function(x) unlist(strsplit(x, split = "[.]")))
  
  names.col_dgp <- sapply(split.dgp, function(x) paste(x[2], x[3], sep = "_"))
  
  max_n           <- max(sapply(split.dgp, function(x) as.double(x[3])), na.rm = T)
  
  # Need to know how many Datasets
  data_try     <- readRDS(file = paste0(path_in, "/" ,files_dgp[1]))
  N_data       <- nrow(data_try$table.BF)
  
  data.BF <- as_tibble(matrix(NA_real_, nrow = N_data, ncol = length(files_dgp) + 2))
  colnames(data.BF)   <- c("OLS", "RBE", paste0("BAY_",names.col_dgp))
  
  # Beta
  data.beta <- as_tibble(matrix(NA_real_, nrow = N_data, ncol = length(files_dgp) + 2))
  colnames(data.beta)   <- c("OLS", "RBE", paste0("BAY_",names.col_dgp))
  
  # a_R
  data.a_R <- as_tibble(matrix(NA_real_, nrow = N_data, ncol = length(files_dgp) + 2))
  colnames(data.a_R)   <- c("OLS", "RBE", paste0("BAY_",names.col_dgp))
  
  # effective Beta
  data.beta_eff <- as_tibble(matrix(NA_real_, nrow = N_data, ncol = length(files_dgp) + 2))
  colnames(data.beta_eff)   <- c("OLS", "RBE", paste0("BAY_",names.col_dgp))
  
  for(i in 1:length(files_dgp)){
    
    # Load the data
    data     <- readRDS(file = paste0(path_in, "/" ,files_dgp[i]))
    
    # a_R could be fixed or random
    if(isTRUE(data$Prior$Prior_a_R) & isTRUE(data$Prior$Prior_b_R)){
      
      a_R_in <- paste(data$Prior$a_R[1], data$Prior$a_R[2])
      b_R_in <- paste(data$Prior$b_R[1], data$Prior$b_R[2])
      
    }else if(isTRUE(data$Prior$Prior_a_R)){
      
      a_R_in <- paste(data$Prior$a_R[1], data$Prior$a_R[2])
      b_R_in <- data$Prior$b_R
      
    }else{
      
      a_R_in <- data$Prior$a_R
      b_R_in <- data$Prior$b_R
      
    }
    
    res.phi  <- cbind(paste0("BAY_",names.col_dgp)[i], a_R_in, b_R_in, data$res_phi)
    
    res.beta <- cbind(paste0("BAY_",names.col_dgp)[i], a_R_in, b_R_in, data$res_beta)
    
    res_all_phi_dgp  <- rbind(res_all_phi_dgp, res.phi)
    
    res_all_beta_dgp <- rbind(res_all_beta_dgp, res.beta)
    
    data.BF[,i+2]    <- data$table.BF$BF_med
    
    if(!is.null(data$table.beta)){
      data.beta[,i+2]       <- data$table.beta$BAY
      data.a_R[,i+2]        <- data$table.ar$BAY
      #data.beta_eff[,i+2]   <- data$mean_eff_beta  
    }
    
  }
  
  
  # Add OLS
  data.BF[,"OLS"]     <- data$table.BF$BF_ols
  data.beta[,"OLS"]   <- data$table.beta$OLS
  
  # Add RBE
  data.BF[,"RBE"]     <- data$table.BF$BF_rbe
  data.beta[,"RBE"]   <- data$table.beta$RBE
  
  
  # If we have few results
  if(max_n > 10){
    
    # Order in columns
    if(isTRUE(data$Prior$Prior_a_R)){
      
      data.BF <- data.BF %>%
        relocate(any_of(paste0("5_", 10:max_n)), .after = last_col())
      
    }else{
      
      data.BF <- data.BF %>%
        relocate(any_of(paste0("4_", 10:max_n)), .after = last_col())
      
    }
    
  }
  
  # Change names
  colnames(res_all_phi_dgp)[1]  <- "Prior_type"
  colnames(res_all_beta_dgp)[1] <- "Prior_type"
  
  # Delete/Change for a_R and b_R 
  colnames(res_all_phi_dgp)[2]  <- "a_R"
  colnames(res_all_beta_dgp)[2] <- "a_R"
  
  colnames(res_all_phi_dgp)[3]  <- "b_R"
  colnames(res_all_beta_dgp)[3] <- "b_R"
  
  #### Beta #### 
  beta_all <- res_all_beta_dgp
  
  beta_all_methods <- beta_all %>% 
    filter(!(Method %in% c("OLS", "RBE")))
  
  beta_all_OLS_RBE_avg <- beta_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi" , "Prior_R2"), ~first(.x)),
              across(N_EFF:SAV_sam, ~ mean(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  beta_all_OLS_RBE_max <- beta_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi" , "Prior_R2"), ~first(.x)),
              across(N_EFF:SAV_sam, ~ max(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  beta_all_methods <- rbind(beta_all_OLS_RBE_avg,  beta_all_methods)
  
  #### PHI #### 
  
  phi_all <- res_all_phi_dgp
  
  phi_all_methods <- phi_all %>% 
    filter(!(Method %in% c("OLS", "RBE")))
  
  phi_all_OLS_RBE_avg <- phi_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi",  "Prior_R2"), ~first(.x)),
              across(N_EFF:Var, ~ mean(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  phi_all_OLS_RBE_max <- phi_all %>% 
    group_by(Method) %>% 
    summarise(across(c("Prior_type", "a_R", "b_R", "Prior_phi",  "Prior_R2"), ~first(.x)),
              across(N_EFF:Var, ~ max(.x, na.rm = TRUE))) %>% 
    filter(Method %in% c("OLS", "RBE")) %>% 
    relocate(Method, .after = Prior_type)
  
  phi_all_methods <- rbind(phi_all_OLS_RBE_avg,  phi_all_methods)
  
  
  
  return(list(phi_all_methods  = phi_all_methods,
              beta_all_methods = beta_all_methods,
              data.BF          = data.BF,
              data.beta        = data.beta,
              data.a_R         = data.a_R,
              data.beta_eff    = data.beta_eff))
}

# Plot Generalized histogram for beta
fun.plot.sim_data.beta <- function(date_in, DGP, index_bay){
  
  date_in$BAY_5_8[index_bay] <- NA_real_
  
  tbl.gen.hist <- tibble(BAY = date_in$BAY_5_8,
                         RBE = date_in$RBE,
                         OLS = date_in$OLS)
  
  q_low  <- min(apply(tbl.gen.hist, 2, function(x) quantile(x, 0.01, na.rm = T)))
  q_high <- max(apply(tbl.gen.hist, 2, function(x) quantile(x, 0.99, na.rm = T)))
  
  tbl.gen.hist <- tbl.gen.hist %>% 
    pivot_longer(cols = BAY:OLS, names_to = "Variable", values_to = "Value") %>% 
    mutate(Variable = as.factor(Variable))
  
  res <- ggplot(tbl.gen.hist, aes(x = Value, fill = Variable)) + 
    xlim(c(q_low,
           q_high)) +
    geom_density(alpha = 0.3) + 
    xlab(expression(widehat(beta))) +
    # xlab("") +
    ylab("Density") + 
    ggtitle(paste0("DGP ", DGP)) +
    theme(text = element_text(size = 20)) +
    theme(legend.key.size = unit(1, 'cm')) +
    labs(fill = "Estimation Method")
  
  return(res)
  
}

# Used to plot the aggregated densities.
fun.plot_aggregated_density_beta <- function(var_in, 
                                             left_in, 
                                             right_in, 
                                             y_lim_in, 
                                             title_in){
  
  dens_1   <- density(x    = var_in[[1]], 
                      from = left_in, 
                      to   = right_in, 
                      n    = 2e3)
  
  dens_all <- sapply(1:length(var_in), function(x) density(x    = var_in[[x]], 
                                                           from = left_in, 
                                                           to   = right_in, 
                                                           n    = 2e3)$y)

  Prior_cont_ar_in$a_R <- c(0.1, 0.5)
  Prior_cont_ar_in$b_R <- "1"
  Prior_cont_ar_in$Prior_variance <- "r2"
  Prior_cont_ar_in$Prior_a_R <- TRUE
  
  prior.beta <- fun.sim_prior_beta(T_len = 1e6, res = list(Prior = Prior_cont_ar_in))
  
  #dens_all.avg <- rowMeans(dens_all)
  dens_all.med <- apply(dens_all, 1, function(x) quantile(x, 0.5))
  dens_all.low <- apply(dens_all, 1, function(x) quantile(x, 0.05))
  dens_all.hig <- apply(dens_all, 1, function(x) quantile(x, 0.95))
  
  data_in <- tibble(x          = dens_1$x, 
                    dens.med   = dens_all.med,
                    dens.low   = dens_all.low,
                    dens.hig   = dens_all.hig,
                    dens.prior = density(x    = prior.beta$beta_sample, 
                                         from = left_in, 
                                         to   = right_in, 
                                         n    = 2e3)$y)
  
  data_in.longer <- data_in %>% 
    pivot_longer(cols = dens.med:dens.prior, names_to = "Variable", values_to = "Value")
  
  
  lines_name <- c("dens.low"   = "dotted",
                  "dens.med"   = "dashed",
                  "dens.hig"   = "dotted",
                  "dens.prior" = "solid")
  
  labels <- list("5% quantile",
                 "50% quantile",
                 "95% quantile",
                 "Prior")
  
  p <- ggplot(data = data_in.longer, aes(x = x, y = Value)) +
    geom_line(aes(linetype = Variable, color = Variable)) + 
    scale_linetype_manual(values = lines_name,
                          labels = labels) +
    scale_color_manual(values = c("black", "black", "black", "red"),
                       labels  = labels) +
    coord_cartesian(xlim = c(left_in, right_in),
                    ylim = c(0, y_lim_in)) +
    xlab(expression(beta)) +
    ylab("Density") +
    ggtitle(title_in) +
    theme(text = element_text(size = 20)) +
    theme(legend.key.size = unit(1, 'cm'))
  
  return(p)
  
  
}
# Plot BF for DGP 1 and DGP 2 
fun.sim_data_plot_BF <- function(data_in, 
                                 Prior_type,
                                 k,
                                 index_bay_1 = NULL,
                                 index_bay_2 = NULL,
                                 a_R_in      = NULL,
                                 b_R_in      = NULL){
  
  
  # Plot name
  prior_number <- as.double(strsplit(Prior_type, "_")[[1]][3])
  
  if(Prior_type == "1_6"){
    
    plot_name <- paste0("BF for Prior: ", Prior_type)
    
    data_0  <- data_in$data.BF_1[,Prior_type, drop = T]
    data_1  <- data_in$data.BF_2[,Prior_type, drop = T]
    
  } else if(Prior_type %in% c("OLS", "RBE")){
    
    plot_name <- paste0("t-values: ", Prior_type)  
    
    xlab_in <- "t-value"
    
    data_0  <- data_in$data.BF_1[,Prior_type, drop = T]
    data_1  <- data_in$data.BF_2[,Prior_type, drop = T]
    
  } else{
    
    plot_name <- paste0("Bayes Factor")  
    
    xlab_in <- "Bayes Factor"
    
    data_0  <- data_in$data.BF_1[,Prior_type, drop = T][-index_bay_1]
    data_1  <- data_in$data.BF_2[,Prior_type, drop = T][-index_bay_2]
  }
  
  shape_in <- c("OLS" = 3, "RBE" = 4, "Prior" = 2, "Posterior" = 17)
  
  
  
  q_left  <-  min(quantile(data_0, 0.01, na.rm = T), quantile(data_1, 0.01, na.rm = T))
  q_right <-  max(quantile(data_0, 0.99, na.rm = T), quantile(data_1, 0.99, na.rm = T))
  
  # For the OLS and RBE we have slightly different look  
  if(Prior_type %in% c("OLS", "RBE")){
    
    dens_0  <- density(data_0, na.rm = T)
    dens_1  <- density(data_1, na.rm = T)
    
    # We also differently compute errors since we have t-stats
    data_0  <- tibble(x = dens_0$x, y = dens_0$y, Variable = "DGP 0") %>% 
      mutate(variable = case_when(
        (x >= k) ~ "Error_0",
        (x <= -k) ~ "Error_0_small",
        TRUE ~ "Okay"))
    
    data_1  <- tibble(x = dens_1$x, y = dens_1$y, Variable = "DGP 1") %>% 
      mutate(variable = case_when(
        (x <= k & x >= -k) ~ "Error_1",
        TRUE ~ "Okay"))
    
  }else{
    
    dens_0  <- density(data_0, from = 0, na.rm = T)
    dens_1  <- density(data_1, from = 0, na.rm = T)
    
    data_0  <- tibble(x = dens_0$x, y = dens_0$y, Variable = "DGP 0") %>% 
      mutate(variable = case_when(
        (x <= k) ~ "Error_0",
        TRUE ~ "Okay"))
    
    data_1  <- tibble(x = dens_1$x, y = dens_1$y, Variable = "DGP 1") %>% 
      mutate(variable = case_when(
        (x >= k) ~ "Error_1",
        TRUE ~ "Okay"))
    
    q_left  <- 0
    q_right <- 15
    
  }
  
  
  data_all <- bind_rows(data_0, data_1)
  
  # Plot
  res <- ggplot(data_all, aes(x = x, y = y, col = Variable)) + geom_line() +
    geom_area(data = filter(data_all, variable == 'Error_0'), fill = 'red', alpha = 0.1, show.legend = F) +
    geom_area(data = filter(data_all, variable == 'Error_0_small'), fill = 'red', alpha = 0.1, show.legend = F) +
    geom_area(data = filter(data_all, variable == 'Error_1'), fill = 'blue', alpha = 0.1, show.legend = F) +
    geom_vline(xintercept = c(-k,k)) +
    xlim(q_left, q_right) +
    ggtitle(plot_name) + 
    ylab("Density") +
    xlab(xlab_in) +
    theme(text = element_text(size = 20)) +
    theme(legend.key.size = unit(1, 'cm'))
  
  return(res)
  
  
}

# Plot CSSE for the real data
fun.real_data_csse <- function(data_in, name_y){
  
  data.plot <- data_in %>% 
    mutate(err_mean = cumsum(err_mean^2),
           err_bay  = cumsum(err_bay^2),
           err_ols  = cumsum(err_ols^2),
           err_rbe  = cumsum(err_rbe^2)) %>% 
    mutate(err_bay = err_bay - err_mean,
           err_ols = err_ols - err_mean,
           err_rbe = err_rbe - err_mean) %>% 
    rename(BAY = err_bay,
           OLS = err_ols,
           RBE = err_rbe) %>% 
    dplyr::select(BAY, OLS, RBE, dates) %>% 
    pivot_longer(BAY:RBE) %>% 
    rename(CSSE     = value,
           Variable = name)
  
  # Plot
  res <- ggplot(data = data.plot, aes(x = dates, y = CSSE, color = Variable)) + 
    geom_line() +
    ggtitle(name_y) +
    geom_hline(yintercept = 0) +
    xlab("") + 
    ylab(expression(paste(Delta ,'CSSE'))) +
    theme(text = element_text(size = 20)) +
    theme(legend.key.size = unit(1, 'cm')) +
    scale_x_date(date_breaks = "10 year", date_labels =  "%Y") 
  
  
  return(res)
  
  
}


# Plot b prior and posterior
fun.real_data_plot_beta <- function(data_in, 
                                    plot_name, 
                                    name_x   = NULL, 
                                    name_y   = "Density",
                                    shape_in = NULL,
                                    x_lim    = 0.25,
                                    y_lim    = NULL){
  
  res.summary <- fun.summary_beta(list(data_in))
  
  res.prior   <- fun.sim_prior_beta(T_len = 1e6, 
                                    res   = data_in)
  # Create a data 
  data.plot <- as_tibble(rbind(
    cbind("Posterior", data_in$beta),
    cbind("Prior", res.prior$beta_sample)
  )) %>% 
    rename(Variable    = V1,
           Value       = V2) %>% 
    mutate(Value = as.double(Value)) %>% 
    filter(between(Value, left = -1, right = 1))
  
  data.plot.prior <- data.plot %>% 
    filter(Variable == "Prior")
  
  data.plot.posterior <- data.plot %>% 
    filter(Variable == "Posterior")
  
  data.mean <- as_tibble(cbind(c("Prior", "Posterior" , "RBE", "OLS"),
                               c(0, res.summary$Mean)
  )) %>% 
    rename(Variable    = V1,
           Value       = V2) %>% 
    mutate(Value       = as.double(Value),
           Variable    = as_factor(Variable))
  
  if(is.null(shape_in)){
    shape_in <- c("OLS" = 0, "RBE" = 1, "Prior" = 2, "Posterior" = 17)  
  }
  
  
  # Plot
  res <- ggplot() + 
    stat_density(data = data.plot.posterior, 
                 aes(x = Value),
                 color = "black",
                 alpha = 0.1,
                 linetype = "dashed",
                 bw = c(0.009938),
                 n  = c(3001)) +
    stat_density(data = data.plot.prior, 
                 aes(x = Value),
                 color = "black", 
                 trim = T,
                 alpha = 0,
                 bw = c(0.015),
                 n  = c(1e4)) +
    geom_point(data      = data.mean, 
               aes(x     = Value, 
                   y     = 0, 
                   shape = Variable),
               alpha = 1,
               size  = 5) +
    scale_shape_manual(values = shape_in) +
    coord_cartesian(xlim = c(-x_lim, x_lim),
                    ylim = c(0, y_lim)) +
    ggtitle(plot_name) +
    theme(text = element_text(size = 20)) +
    theme(legend.key.size = unit(1, 'cm'))
  
  
  if(!is.null(name_x)){
    res <- res +
      xlab(name_x) +
      ylab(name_y)
  }
  
  
  return(res)
  
  
}
# Traceplots
fun.real_data_traceplot <- function(value_in, name_y){
  
  data.plot <- tibble(Value = value_in, 
                      Draws = 1:length(value_in))
  
  
  res <- ggplot(data  = data.plot, 
                aes(x = Draws,
                    y = Value)) +
    geom_line() + 
    ylab(name_y) +  
    theme(text = element_text(size = 20)) +
    theme(legend.key.size = unit(1, 'cm'))
  
  
  return(res)
  
  
}

# Combine all results from the function fun.posterior_results_beta
fun.posterior_results_beta_table <- function(data_in, prior_in, alpha_in = 0.05){
  
  
  res     <- list(fun.posterior_results_beta(data_in = data_in, prior_type = "r2"))
  N_prior <- 1
  
  # Create Table
  table.result             <- as.data.frame(matrix(NA_real_, nrow = N_prior + 2, ncol = 6))
  colnames(table.result)   <- c("Method", "K_med", "K_avg", "K_sam", "HDI_l", "HDI_h")
  table.result$Method      <- c(paste0("Prior_", 1:N_prior), "OLS", "RBE")
  table.result$K_med[1:N_prior]  <- sapply(res, function(x) x$K_med)
  table.result$K_avg[1:N_prior]  <- sapply(res, function(x) x$K_avg)
  table.result$K_sam[1:N_prior]  <- sapply(res, function(x) x$K_sam)
  table.result[1:N_prior, c("HDI_l", "HDI_h")]   <- hdi(data_in$beta, credMass = 1 - alpha_in)
  
  
  table.result[N_prior + 1, c("HDI_l", "HDI_h")] <- f_conf_interval_ols_single_beta(data_in, alpha_in = alpha_in)$conf
  table.result[N_prior + 2, c("HDI_l", "HDI_h")] <- f_conf_interval_bias_cor_single_beta(data_in, alpha_in = alpha_in)$conf
  
  table.result[N_prior + 1, 2:4] <- data_in$ols$coefficients[2,4]
  table.result[N_prior + 2, 2:4] <- fun.extract_p_value_rbe(data_in)$rbe.p_value
  
  
  return(table.result)
  
}



#### Supplementary ####

# Bias correction model
f_bias_correction <- function(data, type){
  
  # Fit first stage regression
  lm   <- lm(y ~ x, data = data)
  res  <- summary(lm)
  
  # Needed values
  n             <- length(data$x)
  phi_estimated <- res$coefficients[2,1]
  
  # Calculate corrected rho
  if(type == "estimated"){
    phi_corrected <- phi_estimated  
  } else if(type == "corrected"){
    phi_corrected <- phi_estimated + (1  + 3*phi_estimated)/n + 3*(1 + 3*phi_estimated)/n^2
  } else if(type == "true"){
    phi_corrected <- 0.95
    print(paste0("I choose the default option:", phi_corrected))
  }
  
  theta_corrected <- (1 - phi_corrected)*mean(c(data$x[1],data$y))
  
  error_corrected <- data$y - data$x*phi_corrected - theta_corrected
  
  return(list(res             = res,
              phi_corrected   = phi_corrected,
              theta_corrected = theta_corrected,
              error_corrected = error_corrected))
}

# Reference prior, where const - standardization
f_prior_Berger <- function(phi, const = 0.25){
  ifelse(abs(phi) < 1,  (1/(2*pi*sqrt(1 - phi^2))) / const,  (1/(2*pi*abs(phi)*sqrt(phi^2 - 1))) / const )
}

# Function to calculate Sigma_0_beta
fun.Sigma_0_beta <- function(sig2_eta_in, phi_in, sig2_v_in, sig2_x_in, psi_in){
  return(sig2_eta_in * (1 - phi_in^2) * (sig2_v_in/sig2_x_in + psi_in^2))
}

# Prior for a_R
fun.prior_a_R <- function(a_R, a_R_1){
  ifelse(a_R == a_R_1, log(0.5), log(0.5))
}

# Propose a_R
fun.proposal_a_R <- function(a_R_old, a_R_1, a_R_2){
  
  if(a_R_old == a_R_1) sample(x = c(a_R_1, a_R_2), size = 1, prob = c(0.7, 0.3))
  else sample(x = c(a_R_1, a_R_2), size = 1, prob = c(0.3, 0.7))
  
}

# Proposal density that always proposes another value
fun.proposal_a_R_v2 <- function(a_R_old, a_R_1, a_R_2){
  
  if(a_R_old == a_R_1) return(a_R_2)
  else return(a_R_1)
}

# Proposal density that always proposes another value
fun.proposal_a_R_b_R <- function(a_R_old, a_R_in, b_R_in){
  
  if(a_R_old == a_R_in[1]) return(c(a_R_in[2], b_R_in[2]))
  else return(c(a_R_in[1], b_R_in[1]))
  
}

# Evaluate density of proposal
fun.proposal_density_a_R <- function(a_R_prop, a_R_old, a_R_1, a_R_2){
  
  if(a_R_prop == a_R_1 & a_R_old == a_R_1){
    return(log(0.7))
  }else if(a_R_prop == a_R_2 & a_R_old == a_R_1){
    return(log(0.3))
  }else if(a_R_prop == a_R_1 & a_R_old == a_R_2){
    return(log(0.3))
  }else if(a_R_prop == a_R_2 & a_R_old == a_R_2){
    return(log(0.7))
  }
  
}

# Prior beta
fun.sim_prior_beta               <- function(T_len, res){
  
  # Save prior
  Prior <- res$Prior
  
  ## R2 ##
  if((length(Prior$b_R) == 1) & (length(Prior$a_R) == 1)){
    
    R2                     <- rbeta(n      = T_len, 
                                    shape1 = Prior$a_R, 
                                    shape2 = as.double(Prior$b_R)) 
    
  }else if((length(Prior$b_R) == 1) & (length(Prior$a_R) > 1)){
    
    a_R <- sample(x = Prior$a_R, size = T_len, replace = TRUE)
    
    R2                     <- rbeta(n      = T_len, 
                                    shape1 = a_R, 
                                    shape2 = as.double(Prior$b_R)) 
    
  }else if((length(Prior$b_R) > 1) & (length(Prior$a_R) > 1)){
    
    # Important that they a_R and b_R respective.
    ind.beta <- sample(x = 1:length(Prior$a_R), size = T_len, replace = TRUE)
    
    R2                     <- rbeta(n      = T_len, 
                                    shape1 = Prior$a_R[ind.beta], 
                                    shape2 = Prior$b_R[ind.beta])
    
  }else{
    
    print("No Specification found")
    
  }
  
  
  
  sig2_x_sample          <- 1/rgamma(n      = T_len, 
                                     shape  = Prior$v_0_x, 
                                     rate   = Prior$Sigma_0_x)
  
  
  sig2_y_tilde_sample    <- 1/rgamma(n      = T_len, 
                                     shape  = Prior$v_0_v, 
                                     rate   = Prior$Sigma_0_v)
  
  
  
  
  Sigma_beta <- R2/(1 - R2)
  
  psi_sample <- rnorm(n = T_len, mean = 0, sd = sqrt(Prior$Sigma_0_psi))
  
  phi_sample_unif              <- runif(n = T_len)
  phi_sample_refe              <- sin(pi * phi_sample_unif / 2)
  
  sig2_sample      <- sig2_x_sample/(1 - phi_sample_refe^2)
  sig2_y_sample    <- sig2_y_tilde_sample + sig2_x_sample * psi_sample^2
  
  Sigma_0_beta <- Sigma_beta * sig2_y_sample * (1 - phi_sample_refe^2) / sig2_x_sample
  
  res <- 1/sqrt(Sigma_0_beta)
  
  beta_est_0_med <- median(res) / sqrt(2*pi)
  beta_est_0_avg <- mean(res) / sqrt(2*pi)
  
  beta_sample <- rnorm(n    = T_len, 
                       mean = 0, 
                       sd   = sqrt(Sigma_0_beta))
  
  
  beta_est_0_sam <- fun.density.beta(beta_sample, x_left = -1, x_right = 1)

  return(list(beta_sample    = beta_sample, 
              res            = res,
              beta_est_0_med = beta_est_0_med,
              beta_est_0_avg = beta_est_0_avg,
              beta_est_0_sam = beta_est_0_sam))
  
}

# Estimates the density for beta using the density function
fun.density.beta                 <- function(object, x_left = -1, x_right = 1){
  
  density_est <- density(object, from = x_left, to = x_right)
  
  return(approx(density_est$x, density_est$y, xout = 0, n = 1000)$y)
}

# Posterior beta
fun.sim_posterior_beta             <- function(res, index = NULL){
  
  # Priors
  Sigma_0_alpha_2 <- res[[1]]$Prior$Sigma_0_alpha_2
  mu_0_alpha_2    <- res[[1]]$Prior$mu_0_alpha_2
  mu_0_beta       <- res[[1]]$Prior$mu_0_beta
  mu_2_0          <- matrix(c(mu_0_alpha_2, mu_0_beta), ncol = 1)
  
  post_dens_avg <- {}
  post_dens_med <- {}
  post_dens_sam <- {}
  
  #k=1
  for(k in 1:length(res)){
    
    # Create the regressor matrix
    X_2_con  <- cbind(1, res[[k]]$Data$x_observe) 
    y_obs    <- res[[k]]$Data$y_observe
    cross_XX <- crossprod(X_2_con, X_2_con)
    cross_XY <- crossprod(X_2_con, y_obs)
    
    
    # Extract parameters
    # dim(Sigma_0_beta)
    Sigma_0_beta    <- res[[k]]$Prior$Sigma_0_beta
    sig2_v          <- res[[k]]$sig2_v
    sig2_x          <- res[[k]]$sig2_x
    psi             <- res[[k]]$psi
    
    # Find the sig2_y (okay, matrix summation)
    sig2_y <- sig2_v + sig2_x * psi^2
    
    #+ Create Prior for the variance
    Sigma_2_0_inv <- lapply(Sigma_0_beta, function(x) diag(c(1/Sigma_0_alpha_2, 1/x)))
    
    # Find the parameters
    term_1    <- lapply(sig2_y, function(x) cross_XX/x)
    term_2    <- Map("+", term_1, Sigma_2_0_inv)
    Sigma_2_n <- lapply(term_2, function(x) solve(x))
    
    term_3    <- lapply(sig2_y, function(x) cross_XY/x)
    term_4    <- lapply(Sigma_2_0_inv, function(x) x %*% mu_2_0)
    term_5    <- Map("+", term_3, term_4)
    mu_2_n    <- Map("%*%", Sigma_2_n, term_5)
    
    b_N <- unlist(map(mu_2_n, 2))
    B_N <- unlist(map(Sigma_2_n, 4))
    
    
    post_dens_avg[k] <- mean(exp(-b_N^2 / (2*B_N)) / sqrt(B_N)) / sqrt(2*pi)
    post_dens_med[k] <- median(exp(-b_N^2 / (2*B_N)) / sqrt(B_N)) / sqrt(2*pi)
    
    post_dens_sam[k]   <- fun.density.beta(res[[k]]$beta)
    
  }
  
  # Delete bad samples
  if(!is.null(index)){
    post_dens_med <- post_dens_med[-index]
    post_dens_avg <- post_dens_avg[-index]
    post_dens_sam <- post_dens_sam[-index]
  }
  
  
  return(list(post_dens_med = post_dens_med, 
              post_dens_avg = post_dens_avg, 
              post_dens_sam = post_dens_sam))
  
} 
fun.sim_posterior_beta_given_b_t   <- function(res){
  
  # Priors
  Sigma_0_alpha_2 <- res[[1]]$Prior$Sigma_0_alpha_2
  mu_0_alpha_2    <- res[[1]]$Prior$mu_0_alpha_2
  mu_0_beta       <- res[[1]]$Prior$mu_0_beta
  mu_2_0          <- matrix(c(mu_0_alpha_2, mu_0_beta), ncol = 1)
  
  post_dens_avg <- {}
  post_dens_med <- {}
  post_dens_sam <- {}
  
  #k=1
  for(k in 1:length(res)){
    
    b_N <- res[[k]]$mu_2_beta_n   
    B_N <- res[[k]]$Sigma_2_beta_n
    
    
    post_dens_avg[k] <- mean(exp(-b_N^2 / (2*B_N)) / sqrt(B_N)) / sqrt(2*pi)
    post_dens_med[k] <- median(exp(-b_N^2 / (2*B_N)) / sqrt(B_N)) / sqrt(2*pi)
    
    beta_in            <- res[[k]]$beta
    post_dens_sam[k]   <- fun.density.beta(beta_in, x_left = -1, x_right = 1)
    
  }
  
  
  return(list(post_dens_med = post_dens_med, 
              post_dens_avg = post_dens_avg, 
              post_dens_sam = post_dens_sam))
  
} 


# Summary for the beta for the real data
fun.summary_beta <- function(data_in){
  
  res_beta <- cbind(sapply(data_in, function(x) mean(x$psi)),
                    sapply(data_in, function(x) mean(x$phi)),
                    sapply(data_in, function(x) mean(x$beta)),
                    sapply(data_in, function(x) sd(x$beta)),
                    sapply(data_in, function(x) estimate_mode(x$beta)),
                    sapply(data_in, function(x) skewness(x$beta)),
                    sapply(data_in, function(x) kurtosis(x$beta)))
  res_beta <- rbind(res_beta,
                    c(data_in[[1]]$ols_corrected$coefficients[3,1],
                      data_in[[1]]$res_bias_corrected$phi_corrected,
                      data_in[[1]]$ols_corrected$coefficients[2,1],
                      f_conf_interval_bias_cor_single_beta(object = data_in[[1]])$se,
                      data_in[[1]]$ols_corrected$coefficients[2,1], 
                      0, 
                      3),
                    c(0,
                      data_in[[1]]$res_bias_corrected$res$coefficients[2,1],
                      data_in[[1]]$ols$coefficients[2,1], 
                      data_in[[1]]$ols$coefficients[2,2],
                      data_in[[1]]$ols$coefficients[2,1],
                      0, 
                      3))
  
  
  res_beta <- cbind(c(paste0("Prior_", 1:length(data_in)), "RBE", "OLS"),
                    res_beta)
  
  res_beta <- as_tibble(res_beta) %>% 
    set_names(c("Method", "Psi", "Phi","Mean", "SD", "Mode", "Skewness", "Kurtosis")) %>% 
    mutate(across(Psi:Kurtosis, as.double))
  
  return(res_beta)
  
}






