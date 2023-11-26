# Load the prepared data
data_load <- readRDS(paste0(path, "data/2.real_data/tbl.goyal_annual_v3.rds"))

# Choose which varaible to analyze
indicator_var_oos <- 1

# Settings
MC_draws <- 1e5
burnin   <- 1e4
thining  <- 30

Prior_cont_ar_in$a_R            <- c(0.1, 0.5) 
Prior_cont_ar_in$b_R            <- as.character(1)
Prior_cont_ar_in$Prior_variance <- "r2"

# Data
tbl.data <- data_load %>% 
  filter(between(date, left = as.Date("1926-01-01"), right = as.Date("2022-12-01"))) 


vars_in <- c("dp.crsp", "b.m", "ep.sp", "dg.sp", "dfy", "tms")


# Choose 1 variable
var_in <- sym(vars_in[indicator_var_oos])

tbl.data <- tbl.data %>% 
  mutate(x_in = !!var_in)

# FOR loop
T_len               <- nrow(tbl.data)
t_start             <- 41
window_size         <- 40
j                   <- 1
tbl.errors.rol      <- matrix(NA_real_, nrow = length(t_start:(T_len-1)), ncol = 6)
tbl.errors.exp      <- matrix(NA_real_, nrow = length(t_start:(T_len-1)), ncol = 6)
res.rol.all         <- list()
res.exp.all         <- list()

for(i in t_start:(T_len-1)){
  
  
  res.rol <- fun.oos_estimation(data_in = tbl.data, 
                                type    = "rol", 
                                i       = i, 
                                window  = window)
  
  res.exp <- fun.oos_estimation(data_in = tbl.data, 
                                type    = "exp", 
                                i       = i, 
                                window  = window)
  
  
  
  # Save results
  res.rol.all[[j]]        <- res.rol$res.bayes
  tbl.errors.rol[j, ]     <- res.rol$errors
  
  res.exp.all[[j]]        <- res.exp$res.bayes
  tbl.errors.exp[j, ]     <- res.exp$errors
  
  # Progress
  print(j)
  
  # Save indicator
  j <- j + 1
  
}

