#### Table 1 ####

# Save the results
results <- fun.combine_results(folder = "v4")

results.phi <- results$phi_all_methods %>% 
  mutate(Prior_type  = case_when(Method %in% c("OLS", "RBE") ~ "None",
                                 TRUE ~ Prior_type)) %>% 
  mutate(SD.dgp_1    = 100 * sqrt(Var.dgp_1 / 1000),
         RMSE.dgp_1  = 100 * sqrt(MSE.dgp_1 / 1000),
         SD.dgp_2    = 100 * sqrt(Var.dgp_2 / 1000),
         RMSE.dgp_2  = 100 * sqrt(MSE.dgp_2 / 1000))


results.phi.1 <- results.phi %>% 
  filter(Prior_type %in% c("None", "BAY_5_8")) %>% 
  select(Method, 
         Bias.dgp_1,
         SD.dgp_1,
         MAE.dgp_1,
         RMSE.dgp_1) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))

results.phi.2 <- results.phi %>% 
  filter(Prior_type %in% c("None", "BAY_5_8")) %>% 
  select(Method, 
         Bias.dgp_2,
         SD.dgp_2,
         MAE.dgp_2,
         RMSE.dgp_2) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))

results.phi.1
results.phi.2

rm(list = ls())

#### Table 2 ####

# Save the results
results             <- fun.combine_results(folder = "v4")
results.fixed_dgp_0 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_1")
results.fixed_dgp_1 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_2")


results.beta <- results$beta_all_methods %>% 
  mutate(Prior_type  = case_when(Method %in% c("OLS", "RBE") ~ "None",
                                 TRUE ~ Prior_type)) %>% 
  select(-Prior_phi) %>% 
  mutate(SD.dgp_1    = 100 * sqrt(Var.dgp_1 / 1000),
         RMSE.dgp_1  = 100 * sqrt(MSE.dgp_1 / 1000),
         SD.dgp_2    = 100 * sqrt(Var.dgp_2 / 1000),
         RMSE.dgp_2  = 100 * sqrt(MSE.dgp_2 / 1000))

results.beta.1.prior <- results.beta %>% 
  filter(Prior_type %in% c("None", "BAY_5_8")) %>% 
  select(Method, 
         Bias.dgp_1,
         SD.dgp_1,
         MAE.dgp_1,
         RMSE.dgp_1,
         SAV_med.dgp_1) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))

results.beta.1.fixed <- results.fixed_dgp_0$beta_all_methods %>% 
  filter(Method == "BAY") %>% 
  mutate(SD.dgp_1    = 100 * sqrt(Var / 1000),
         RMSE.dgp_1  = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD.dgp_1,
         MAE,
         RMSE.dgp_1,
         SAV_med)

results.beta.2.prior <- results.beta %>% 
  filter(Prior_type %in% c("None", "BAY_5_8")) %>% 
  select(Method, 
         Bias.dgp_2,
         SD.dgp_2,
         MAE.dgp_2,
         RMSE.dgp_2,
         SAV_med.dgp_2) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))

results.beta.2.fixed <- results.fixed_dgp_1$beta_all_methods %>% 
  filter(Method == "BAY") %>% 
  mutate(SD.dgp_2    = 100 * sqrt(Var / 1000),
         RMSE.dgp_2  = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD.dgp_2,
         MAE,
         RMSE.dgp_2,
         SAV_med)

results.beta.1.prior
results.beta.1.fixed

results.beta.2.prior
results.beta.2.fixed

rm(list = ls())

#### Table 3(4) ####

# Take necessary names
names_files <- list.files(paste0(path, "data/2.real_data/res/"))

# Object to save
res_all <- list()

# Save
for(i in 1:length(names_files)){
  
  res_all[[i]] <- read_rds(paste0(path, "data/2.real_data/res/", names_files[i]))
  
  print(i)
  
}

res.beta <- lapply(res_all,
                   function(x) as_tibble(fun.posterior_results_beta_table(data_in  = x,  
                                                                          prior_in = rep("r2", 2), alpha_in = 0.05)))

# Summary
summary.beta <- lapply(1:length(names_files), function(x) fun.summary_beta(list(res_all[[x]])))

# Prepare results
for(i in 1:length(names_files)){
  
  if(i == 1){
    res.beta.table.all <- summary.beta[[i]] %>% 
      inner_join(res.beta[[i]], by = "Method") %>% 
      select(Method, Psi ,Phi, Mean, SD, K_med) %>% 
      mutate(across(Psi:K_med, ~round(., 3))) %>% 
      mutate(Method = gsub("_", " ", Method, fixed = TRUE)) %>% 
      mutate(Var = str_split(names_files[i], "[.]")[[1]][3])
  }else{
    
    res.beta.table <- summary.beta[[i]] %>% 
      inner_join(res.beta[[i]], by = "Method") %>% 
      select(Method, Psi ,Phi, Mean, SD, K_med) %>% 
      mutate(across(Psi:K_med, ~round(., 3))) %>% 
      mutate(Method = gsub("_", " ", Method, fixed = TRUE)) %>% 
      mutate(Var = str_split(names_files[i], "[.]")[[1]][3])
    
    res.beta.table.all <- res.beta.table.all %>% 
      bind_rows(res.beta.table)
    
    
  }
  
}

res.beta.table.all

rm(list = ls())

#### Table 5 ####


# Names
names_in <- c("dp.crsp", "ep.sp", "b.m", "dg.sp", "dfy","tms")

# Empty tables
tbl.rol <- matrix(NA_real_, nrow = length(names_in),  ncol = 3)
tbl.exp <- matrix(NA_real_, nrow = length(names_in),  ncol = 3)

# Save the results
for(i in 1:length(names_in)){
  
  tbl.errors.rol <- read_rds(paste0(path, "/data/2.real_data/oos/",names_in[i],"/tbl.errors.rol.rds"))
  
  tbl.errors.exp <- read_rds(paste0(path, "/data/2.real_data/oos/",names_in[i],"/tbl.errors.exp.rds"))
  
  # OOS R2
  tbl.rol[i,] <- 1 - colMeans(tbl.errors.rol[,3:5]^2) / mean(tbl.errors.rol[,6]^2)
  tbl.exp[i,] <- 1 - colMeans(tbl.errors.exp[,3:5]^2) / mean(tbl.errors.exp[,6]^2)
  
}

# Prepare for the output
res.exp <- as_tibble(tbl.exp) %>% 
  mutate(across(V1:V3, ~round(100 * ., 3))) %>% 
  mutate(name = c("dp", "ep", "bm", "dg", "dfy", "tms"), .before = V1) %>% 
  set_names(c("Variable", "BAY", "OLS", "RBE"))

res.rol <- as_tibble(tbl.rol) %>% 
  mutate(across(V1:V3, ~round(100 * ., 3))) %>% 
  mutate(name = c("dp", "ep", "bm", "dg", "dfy", "tms"), .before = V1) %>% 
  set_names(c("Variable", "BAY", "OLS", "RBE"))

# Res
res.exp
res.rol

rm(list = ls())

#### Table 6 ####

# We use another notation for multiple DGPs

# beta == 0      - DGP1
# beta == 0.1    - DGP2
# beta == 0.05   - DGP3
# beta == 0.2    - DGP4
# beta == 0.025  - DGP5
# beta == 0.075  - DGP6

results.prior.1 <- fun.combine_results_v2(folder = "v4/dgp_1")
results.prior.2 <- fun.combine_results_v2(folder = "v4/dgp_2")
results.prior.3 <- fun.combine_results_v2(folder = "v8/prior/bR_1/dgp_3")
results.prior.4 <- fun.combine_results_v2(folder = "v8/prior/bR_1/dgp_4")
results.prior.5 <- fun.combine_results_v2(folder = "v8/prior/bR_1/dgp_5")
results.prior.6 <- fun.combine_results_v2(folder = "v8/prior/bR_1/dgp_6")

results.fixed.1 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_1")
results.fixed.2 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_2")
results.fixed.3 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_3")
results.fixed.4 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_4")
results.fixed.5 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_5")
results.fixed.6 <- fun.combine_results_v2(folder = "v8/fixed/bR_1/dgp_6")

results.fixed.1$beta_all_methods$DGP <- "1"
results.fixed.2$beta_all_methods$DGP <- "2"
results.fixed.3$beta_all_methods$DGP <- "3"
results.fixed.4$beta_all_methods$DGP <- "4"
results.fixed.5$beta_all_methods$DGP <- "5"
results.fixed.6$beta_all_methods$DGP <- "6"

results.prior.1$beta_all_methods$DGP <- "1"
results.prior.2$beta_all_methods$DGP <- "2"
results.prior.3$beta_all_methods$DGP <- "3"
results.prior.4$beta_all_methods$DGP <- "4"
results.prior.5$beta_all_methods$DGP <- "5"
results.prior.6$beta_all_methods$DGP <- "6"

results.fixed.1 <- results.fixed.1$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.2 <- results.fixed.2$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.3 <- results.fixed.3$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.4 <- results.fixed.4$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.5 <- results.fixed.5$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.6 <- results.fixed.6$beta_all_methods %>% 
  filter(Method == "BAY")


results.beta.1 <- results.prior.1$beta_all_methods %>%
  filter(Prior_type %in% c("BAY_5_8") | Method %in% c("OLS", "RBE")) %>% 
  rbind(results.fixed.1) %>% 
  relocate(DGP, .before = everything()) %>% 
  select(-Prior_phi) %>% 
  mutate(SD   = 100 * sqrt(Var / 1000),
         RMSE = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD,
         MAE,
         RMSE,
         SAV_med) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))


results.beta.rest <- results.prior.5$beta_all_methods %>%
  rbind(results.fixed.5) %>% 
  
  rbind(results.prior.3$beta_all_methods) %>%
  rbind(results.fixed.3) %>% 
  
  rbind(results.prior.6$beta_all_methods) %>%
  rbind(results.fixed.6) %>% 
  
  rbind(results.prior.2$beta_all_methods %>% 
          filter(Prior_type %in% c("BAY_5_8") | Method %in% c("OLS", "RBE"))
        ) %>%
  rbind(results.fixed.2) %>% 

  rbind(results.prior.4$beta_all_methods) %>%
  rbind(results.fixed.4) %>% 
  
  relocate(DGP, .before = everything()) %>% 
  select(-Prior_phi) %>% 
  mutate(SD   = 100 * sqrt(Var / 1000),
         RMSE = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD,
         MAE,
         RMSE,
         SAV_med) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))


results.beta.1
results.beta.rest


rm(list = ls())


#### Table 7 ####

# We do not report OLS and RBE since they are the same as in Table 6

# We use another notation for multiple DGPs

# beta == 0      - DGP1
# beta == 0.1    - DGP2
# beta == 0.05   - DGP3
# beta == 0.2    - DGP4
# beta == 0.025  - DGP5
# beta == 0.075  - DGP6

results.prior.1 <- fun.combine_results_v2(folder = "v4/dgp_1")
results.prior.2 <- fun.combine_results_v2(folder = "v4/dgp_2")
results.prior.3 <- fun.combine_results_v2(folder = "v8/prior/bR_0.5/dgp_3")
results.prior.4 <- fun.combine_results_v2(folder = "v8/prior/bR_0.5/dgp_4")
results.prior.5 <- fun.combine_results_v2(folder = "v8/prior/bR_0.5/dgp_5")
results.prior.6 <- fun.combine_results_v2(folder = "v8/prior/bR_0.5/dgp_6")

results.fixed.1 <- fun.combine_results_v2(folder = "v8/fixed/bR_0.5/dgp_1")
results.fixed.2 <- fun.combine_results_v2(folder = "v8/fixed/bR_0.5/dgp_2")
results.fixed.3 <- fun.combine_results_v2(folder = "v8/fixed/bR_0.5/dgp_3")
results.fixed.4 <- fun.combine_results_v2(folder = "v8/fixed/bR_0.5/dgp_4")
results.fixed.5 <- fun.combine_results_v2(folder = "v8/fixed/bR_0.5/dgp_5")
results.fixed.6 <- fun.combine_results_v2(folder = "v8/fixed/bR_0.5/dgp_6")

results.fixed.1$beta_all_methods$DGP <- "1"
results.fixed.2$beta_all_methods$DGP <- "2"
results.fixed.3$beta_all_methods$DGP <- "3"
results.fixed.4$beta_all_methods$DGP <- "4"
results.fixed.5$beta_all_methods$DGP <- "5"
results.fixed.6$beta_all_methods$DGP <- "6"

results.prior.1$beta_all_methods$DGP <- "1"
results.prior.2$beta_all_methods$DGP <- "2"
results.prior.3$beta_all_methods$DGP <- "3"
results.prior.4$beta_all_methods$DGP <- "4"
results.prior.5$beta_all_methods$DGP <- "5"
results.prior.6$beta_all_methods$DGP <- "6"

results.fixed.1 <- results.fixed.1$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.2 <- results.fixed.2$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.3 <- results.fixed.3$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.4 <- results.fixed.4$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.5 <- results.fixed.5$beta_all_methods %>% 
  filter(Method == "BAY")
results.fixed.6 <- results.fixed.6$beta_all_methods %>% 
  filter(Method == "BAY")


results.beta.1 <- results.prior.1$beta_all_methods %>%
  filter(Prior_type %in% c("BAY_5_2")) %>% 
  rbind(results.fixed.1) %>% 
  relocate(DGP, .before = everything()) %>% 
  select(-Prior_phi) %>% 
  mutate(SD   = 100 * sqrt(Var / 1000),
         RMSE = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD,
         MAE,
         RMSE,
         SAV_med) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))


results.beta.rest <- results.prior.5$beta_all_methods %>%
  rbind(results.fixed.5) %>% 
  
  rbind(results.prior.3$beta_all_methods) %>%
  rbind(results.fixed.3) %>% 
  
  rbind(results.prior.6$beta_all_methods) %>%
  rbind(results.fixed.6) %>% 
  
  rbind(results.prior.2$beta_all_methods %>% 
          filter(Prior_type %in% c("BAY_5_2"))
  ) %>%
  rbind(results.fixed.2) %>% 
  
  rbind(results.prior.4$beta_all_methods) %>%
  rbind(results.fixed.4) %>% 
  
  filter(Method == "BAY") %>% 
  relocate(DGP, .before = everything()) %>% 
  select(-Prior_phi) %>% 
  mutate(SD   = 100 * sqrt(Var / 1000),
         RMSE = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD,
         MAE,
         RMSE,
         SAV_med) %>% 
  mutate(across(is.numeric, ~round(.x, 2)))


results.beta.1
results.beta.rest


rm(list = ls())

#### Table 8 ####

# We do not report OLS and RBE since they are the same as in Table 6
# We do not report results for non-fixed R^2 since the results are the same as in Table 6

results_1 <- fun.combine_results_v2(folder = "v11/dgp_1")
results_2 <- fun.combine_results_v2(folder = "v11/dgp_2")
results_3 <- fun.combine_results_v2(folder = "v11/dgp_3")
results_4 <- fun.combine_results_v2(folder = "v11/dgp_4")
results_5 <- fun.combine_results_v2(folder = "v11/dgp_5")
results_6 <- fun.combine_results_v2(folder = "v11/dgp_6")


results <- results_1$beta_all_methods %>% 
  bind_rows(results_2$beta_all_methods) %>% 
  bind_rows(results_3$beta_all_methods) %>% 
  bind_rows(results_4$beta_all_methods) %>% 
  bind_rows(results_5$beta_all_methods) %>% 
  bind_rows(results_6$beta_all_methods) %>% 
  mutate(SD   = 100 * sqrt(Var / 1000),
         RMSE = 100 * sqrt(MSE / 1000)) %>% 
  select(Method, 
         Bias,
         SD,
         MAE,
         RMSE,
         SAV_med) %>% 
  filter(Method == "BAY") %>% 
  mutate(across(is.numeric, ~round(.x, 2)))

results

rm(list = ls())

#### Table 9 ####

# Load data
data_load <- read_rds(paste0(path, "data/2.real_data/tbl.goyal_annual_v3.rds"))

# Data
tbl.data <- data_load %>% 
  filter(between(date, left  = as.Date("1926-01-01"), 
                 right = as.Date("2022-12-01"))) 


vars_in <- c("ep.sp", "b.m", "dg.sp", "dfy", "tms")

data_save <- tbl.data %>% 
  dplyr::select(date, vars_in)

# Please refer to the PDF in the folder EViews that explains how to perform stationarity tests.

rm(list = ls())

#### Table 10 ####

# Take necessary names
names_files <- list.files(paste0(path, "data/2.real_data/res/"))

# Save names
names_in <- sapply(names_files, function(x) strsplit(x, split = "[.]")[[1]][3])

# Object to save
res_all <- list()

# Save
for(i in 1:length(names_files)){
  
  res_all[[i]] <- read_rds(paste0(path, "data/2.real_data/res/", names_files[i]))
  
  print(i)
  
}

# Produces the output
fun.effective_summary_table(res_all, names_in = names_in)

rm(list = ls())






#### Figure 1 ####

results       <- fun.combine_results(folder = "v4")
results_dgp_1 <- read_rds(paste0(path, "data/1.sim_data/v4/dgp_1/res.5.8.1-3.1.rds"))
results_dgp_2 <- read_rds(paste0(path, "data/1.sim_data/v4/dgp_2/res.5.8.1-3.2.rds"))

plot.beta_1 <- fun.plot.sim_data.beta(date_in   = results$data.beta_1, 
                                      DGP       = 0,
                                      index_bay = results_dgp_1$index_bay)
plot.beta_2 <- fun.plot.sim_data.beta(date_in   = results$data.beta_2, 
                                      DGP       = 0,
                                      index_bay = results_dgp_2$index_bay)

plot.beta <- ggarrange(plot.beta_1,
                       plot.beta_2,
                       ncol = 2, 
                       nrow = 1, 
                       common.legend = TRUE, 
                       legend="bottom")
plot.beta

rm(list = ls())

#### Figure 2 ####

results       <- fun.combine_results(folder = "v4")
results_dgp_1 <- read_rds(paste0(path, "data/1.sim_data/v4/dgp_1/res.5.8.1-3.1.rds"))
results_dgp_2 <- read_rds(paste0(path, "data/1.sim_data/v4/dgp_2/res.5.8.1-3.2.rds"))

a_R_1 <- results$data.a_R_1$BAY_5_8
a_R_1[results_dgp_1$index_bay] <- NA_real_

a_R_2 <- results$data.a_R_2$BAY_5_8
a_R_2[results_dgp_2$index_bay] <- NA_real_

tbl.gen.hist <- tibble('DGP 0' = a_R_1,
                       'DGP 1' = a_R_2)

tbl.gen.hist <- tbl.gen.hist %>% 
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value") %>% 
  mutate(Variable = as.factor(Variable))

plot.a_R <- ggplot(tbl.gen.hist, aes(x = Value, fill = Variable)) + 
  xlim(0.1, 0.5) +
  geom_density(alpha = 0.3) + 
  xlab(expression(widehat(a[0]^R))) +
  ylab("Density") + 
  geom_point(x = 0.1, y = 0, show.legend = F, shape = 21, colour = "black", fill = "white", size = 3, stroke = 3) + 
  geom_point(x = 0.5, y = 0, show.legend = F, shape = 21, colour = "black", fill = "white", size = 3, stroke = 3) +
  geom_point(x = mean(a_R_1, na.rm = T), y = 0, shape = 21, colour = "black", fill = "red", size = 3, stroke = 3) + 
  geom_point(x = mean(a_R_2, na.rm = T), y = 0, shape = 21, colour = "black", fill = "blue", size = 3, stroke = 3) +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

plot.a_R

rm(list = ls())

#### Figure 3 ####

load(paste0(path, "data/1.sim_data/v9/res_main.RData"))

p1 <- fun.plot_aggregated_density_beta(var_in   = beta_1[-index_1], 
                                       left_in  = -0.1, 
                                       right_in = 0.1, 
                                       y_lim_in = 15, 
                                       title_in = "DGP 0")


p2 <- fun.plot_aggregated_density_beta(var_in   = beta_2[-index_2], 
                                       left_in  = -0.1, 
                                       right_in = 0.3, 
                                       y_lim_in = 15, 
                                       title_in = "DGP 1")


p.all <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

p.all

rm(list = ls())

#### Figure 4(5) ####

names_in   <- c("dp.crsp", "ep.sp", "b.m", "dg.sp", "dfy","tms")
names_plot <- c("DP", "EP", "BM", "DG", "DFY", "TMS")


plot.exp <- list()
plot.rol <- list()

for(i in 1:length(names_in)){
  
  tbl.errors.rol <- read_rds(paste0(path, "data/2.real_data/oos/",names_in[i],"/tbl.errors.rol.rds"))
  
  tbl.errors.exp <- read_rds(paste0(path, "data/2.real_data/oos/",names_in[i],"/tbl.errors.exp.rds"))
  
  tbl.dates <- seq(from = as.Date("1926-01-01"),
                   to   = as.Date("2021-01-01"),
                   by   = "year")[42:96]
  
  # Rolling
  colnames(tbl.errors.rol) <- c("i", "ret_fut", "err_bay", "err_rbe", "err_ols", "err_mean")
  
  tbl.errors <- as_tibble(tbl.errors.rol)
  
  tbl.errors.rol <- tbl.errors.rol %>% 
    bind_cols(dates = tbl.dates)
  
  # Expanding
  colnames(tbl.errors.exp) <- c("i", "ret_fut", "err_bay", "err_rbe", "err_ols", "err_mean")
  
  tbl.errors <- as_tibble(tbl.errors.exp)
  
  tbl.errors.exp <- tbl.errors.exp %>% 
    bind_cols(dates = tbl.dates)
  
  
  plot.exp[[i]] <- fun.real_data_csse(data_in = tbl.errors.exp, 
                                      name_y  = names_plot[i])
  
  plot.rol[[i]] <- fun.real_data_csse(data_in = tbl.errors.rol, 
                                      name_y  = names_plot[i])
  
  
}

plot.all.exp <- ggarrange(plotlist = plot.exp, ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom")
plot.all.rol <- ggarrange(plotlist = plot.rol, ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom")

plot.all.exp
plot.all.rol

rm(list = ls())

#### Figure 6 ####

T_len <- 1e4

Prior_R2_all <- list(
  list(a_R = 0.5, b_R = 1),
  list(a_R = 0.1, b_R = 1),
  list(a_R = sample(x = c(0.1, 0.5), size = T_len, replace = T),
       b_R = 1)
)

# Simulate priors
prior.r2     <- sapply(Prior_R2_all, function(x) density(rbeta(n      = T_len, 
                                                               shape1 = x$a_R,
                                                               shape2 = x$b_R), 
                                                         from = 0, 
                                                         to   = 0.2, 
                                                         n    = 1e4)$y)
prior.r2 <- tibble(x = seq(from = 0, to = 0.2, length = 1e4),
                   unif = 1,
                   as_tibble(prior.r2))

colnames(prior.r2) <- c("x", paste0("Prior_", 1:4))


prior.r2 <- pivot_longer(data      = prior.r2, 
                         cols      = Prior_1:Prior_4, 
                         names_to  = "Variable", 
                         values_to = "Value")

labels_in <- list(expression(paste(a[0]^R == 1)),
                  expression(paste(a[0]^R == 0.5)),
                  expression(paste(a[0]^R == 0.1)),
                  expression(paste(a[0]^R, " random")))


p1 <- ggplot(data = prior.r2, aes(x = x, y = Value)) +
  geom_line(aes(color = Variable)) + 
  scale_color_hue(labels = labels_in) +
  coord_cartesian(xlim = c(0, 0.2),
                  ylim = c(0, 15)) +
  xlab(expression(R^2)) +
  ylab("Density") +
  theme(text = element_text(size = 30)) +
  theme(legend.key.size = unit(1, 'cm'))

p1


rm(list = ls())

#### Figure 7 ####

T_len <- 1e5

dgp_0 <- rnorm(n = T_len, mean = 2, sd = 0.5)
dgp_1 <- rnorm(n = T_len, mean = 0, sd = 0.5)


dens_0  <- density(dgp_0, from = 0, na.rm = T)
dens_1  <- density(dgp_1, from = 0, na.rm = T)

data_0  <- tibble(x = dens_0$x, y = dens_0$y, Variable = "DGP 0") %>% 
  mutate(variable = case_when(
    (x <= 1) ~ "Error_0",
    TRUE ~ "Okay"))

data_1  <- tibble(x = dens_1$x, y = dens_1$y, Variable = "DGP 1") %>% 
  mutate(variable = case_when(
    (x >= 1) ~ "Error_1",
    TRUE ~ "Okay"))


data_all <- bind_rows(data_0, data_1)

# Plot
res <- ggplot(data_all, aes(x = x, y = y, col = Variable)) + geom_line() +
  geom_area(data = filter(data_all, variable == 'Error_0'), fill = 'red', alpha = 0.1, show.legend = F) +
  geom_area(data = filter(data_all, variable == 'Error_0_small'), fill = 'red', alpha = 0.1, show.legend = F) +
  geom_area(data = filter(data_all, variable == 'Error_1'), fill = 'blue', alpha = 0.1, show.legend = F) +
  geom_vline(xintercept = 1) +
  xlim(0, 2) + 
  xlab("Bayes Factor") + 
  ylab("Density") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

res

rm(list = ls())

#### Figure 8 ####

# Save the results
results       <- fun.combine_results(folder = "v4")
results_dgp_1 <- read_rds(paste0(path, "data/1.sim_data/v4/dgp_1/res.5.8.1-3.1.rds"))
results_dgp_2 <- read_rds(paste0(path, "data/1.sim_data/v4/dgp_2/res.5.8.1-3.2.rds"))

res.bay <- fun.sim_data_plot_BF(data_in      = results, 
                                Prior_type   = "BAY_5_8",
                                index_bay_1  = results_dgp_1$index_bay,
                                index_bay_2  = results_dgp_2$index_bay,
                                k            = 1)

res.ols <- fun.sim_data_plot_BF(data_in    = results, 
                                Prior_type = "OLS", 
                                k          = 1.96)

res.rbe <- fun.sim_data_plot_BF(data_in    = results, 
                                Prior_type = "RBE", 
                                k          = 1.96)


fig.BF <- ggarrange(res.bay,
                    ggarrange(res.ols, res.rbe, ncol = 2, nrow = 1, common.legend = FALSE, legend=FALSE), ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom")

fig.BF

rm(list = ls())


#### Figure 9 ####

load(paste0(path, "data/1.sim_data/v9/res.a_R_prior_spec.beta_0.1.T_len_1000.RData"))

beta_1    <- beta
index_1   <- index_bay

load(paste0(path, "data/1.sim_data/v9/res.a_R_prior_spec.beta_0.15.T_len_500.RData"))

beta_2    <- beta
index_2   <- index_bay

load(paste0(path, "data/1.sim_data/v9/res.a_R_prior_spec.beta_0.2.T_len_100.RData"))

beta_3    <- beta
index_3   <- index_bay

p1 <- fun.plot_aggregated_density_beta(var_in   = beta_1[-index_1], 
                                       left_in  = -0.1, 
                                       right_in = 0.3, 
                                       y_lim_in = 30, 
                                       title_in = expression(paste(beta == 0.1, " and ", T == 1000)))

p2 <- fun.plot_aggregated_density_beta(var_in   = beta_2[-index_2], 
                                       left_in  = -0.1, 
                                       right_in = 0.3, 
                                       y_lim_in = 20, 
                                       title_in = expression(paste(beta == 0.15, " and ", T == 500)))

p3 <- fun.plot_aggregated_density_beta(var_in   = beta_3[-index_3], 
                                       left_in  = -0.1, 
                                       right_in = 0.3, 
                                       y_lim_in = 10, 
                                       title_in = expression(paste(beta == 0.2, " and ", T == 100)))

p.all <- ggarrange(p1, p2, p3, ncol = 1, nrow = 3, common.legend = TRUE, legend="bottom")

p.all

rm(list = ls())

#### Figure 10(11) ####

data_load <- read_rds(paste0(path, "data/2.real_data/tbl.goyal_annual_v3.rds"))

# Data
tbl.data.1 <- data_load %>% 
  filter(between(date, left  = as.Date("1926-01-01"), 
                       right = as.Date("2004-01-01"))) 

tbl.data.2 <- data_load %>% 
  filter(between(date, left  = as.Date("1953-01-01"), 
                       right = as.Date("2022-12-01"))) 


p1 <- ggAcf(tbl.data.1$ret.crsp.l, lag.max = 30, main = "Sample 1") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p2 <- ggAcf(tbl.data.2$ret.crsp.l, lag.max = 30, main = "Sample 2") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p3 <- ggPacf(tbl.data.1$ret.crsp.l, lag.max = 30, main = "Sample 1") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p4 <- ggPacf(tbl.data.2$ret.crsp.l, lag.max = 30, main = "Sample 2") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

p.all <- ggarrange(p1, p3, p2 , p4, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

p.all


p1 <- ggAcf(tbl.data.1$dp.crsp, lag.max = 30, main = "Sample 1") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p2 <- ggAcf(tbl.data.2$dp.crsp, lag.max = 30, main = "Sample 2") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p3 <- ggPacf(tbl.data.1$dp.crsp, lag.max = 30, main = "Sample 1") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p4 <- ggPacf(tbl.data.2$dp.crsp, lag.max = 30, main = "Sample 2") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

p.all <- ggarrange(p1, p3, p2 , p4, ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

p.all


rm(list = ls())

#### Figure 12(13) ####

data_load <- read_rds(paste0(path, "data/2.real_data/tbl.goyal_annual_v3.rds"))

# Data
tbl.data <- data_load %>% 
  filter(between(date, left  = as.Date("1926-01-01"), 
                       right = as.Date("2022-12-01"))) 


vars_in <- c("ep.sp", "b.m", "dg.sp", "dfy", "tms")

data_save <- tbl.data %>% 
  dplyr::select(date, vars_in)


## EP ##
p.acf.ep <- ggAcf(tbl.data$ep.sp, lag.max = 30, main = "EP") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p.pacf.ep <- ggPacf(tbl.data$ep.sp, lag.max = 30, main = "EP") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

## BM ##
p.acf.bm <- ggAcf(tbl.data$b.m, lag.max = 30, main = "BM") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p.pacf.bm <- ggPacf(tbl.data$b.m, lag.max = 30, main = "BM") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

## DG ##
p.acf.dg <- ggAcf(tbl.data$dg.sp, lag.max = 30, main = "DG") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p.pacf.dg <- ggPacf(tbl.data$dg.sp, lag.max = 30, main = "DG") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

## DFY ##
p.acf.dfy <- ggAcf(tbl.data$dfy, lag.max = 30, main = "DFY") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p.pacf.dfy <- ggPacf(tbl.data$dfy, lag.max = 30, main = "DFY") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

## TMS ##
p.acf.tms <- ggAcf(tbl.data$tms, lag.max = 30, main = "TMS") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))
p.pacf.tms <- ggPacf(tbl.data$tms, lag.max = 30, main = "TMS") +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, 'cm'))

p.all.1 <- ggarrange(p.acf.ep, p.pacf.ep,
                     p.acf.bm, p.pacf.bm,
                     p.acf.dfy, p.pacf.dfy,
                     ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom")

p.all.1

p.all.2 <- ggarrange(p.acf.dg, p.pacf.dg,
                     p.acf.tms, p.pacf.tms,
                     ncol = 2, nrow = 2, common.legend = TRUE, legend="bottom")

p.all.2


rm(list = ls())


#### Figure 14 ####


res_1 <- read_rds(paste0(path, "data/2.real_data/res/res.real_data.dp_crsp_1.rds")) 
res_2 <- read_rds(paste0(path, "data/2.real_data/res/res.real_data.dp_crsp_2.rds")) 


p1 <- fun.real_data_plot_beta(data_in   = res_1, 
                              plot_name = "Sample 1", 
                              name_x    = expression(beta), 
                              y_lim     = 13)

p2 <- fun.real_data_plot_beta(data_in   = res_2, 
                              plot_name = "Sample 2", 
                              name_x    = expression(beta), 
                              y_lim     = 13)

p.all <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")
p.all


rm(list = ls())



#### Figure 15 ####

res_1 <- read_rds(paste0(path, "data/2.real_data/res/res.real_data.dp_crsp_1.rds")) 
res_2 <- read_rds(paste0(path, "data/2.real_data/res/res.real_data.dp_crsp_2.rds")) 


traceplot.all <- ggarrange(
  fun.real_data_traceplot(value_in = res_1$beta, name_y = expression(beta)),
  fun.real_data_traceplot(value_in = res_2$beta , name_y = expression(beta)),
  
  fun.real_data_traceplot(value_in = res_1$phi , name_y = expression(phi)),
  fun.real_data_traceplot(value_in = res_2$phi , name_y = expression(phi)),
  
  fun.real_data_traceplot(value_in = res_1$psi , name_y = expression(psi)),
  fun.real_data_traceplot(value_in = res_2$psi , name_y = expression(psi)),
  ncol = 2, nrow = 3, common.legend = TRUE, legend="bottom")

traceplot.all
