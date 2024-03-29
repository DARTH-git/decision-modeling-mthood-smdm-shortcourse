#install darthtools in advance

rm(list = ls())      # clear memory (removes all the variables from the workspace)

# 01 Load packages
 
library(devtools)
install_github("DARTH-git/darthtools", force = TRUE)
p_load_gh("DARTH-git/darthtools")
library(darthtools)
library(dampack)
library(reshape2)
library(diagram)

# 02 Load functions

# all functions are in the darthtools package

# 03 Model input

n_time_horizon_yr <- 60                               # time horizon (in years)
cycle_length    <- 1                                  # cycle length in years (use 1/12 for monthly)
n_cycles        <- n_time_horizon_yr / cycle_length   # number of cycles
v_names_cycles  <- paste("cycle", 0:n_cycles)         # cycle names
v_names_states  <- c("Healthy", "Sick", "Dead")       # state names
n_states        <- length(v_names_states)             # number of health states 

### Discounting factors 
d_c <- 0.03 # annual discount rate for costs 
d_e <- 0.03 # annual discount rate for QALYs

### Strategies 
v_names_str     <- c("Standard of Care", # store the strategy names
                     "Treatment A", 
                     "Treatment B")  

n_str <- length(v_names_str)           # number of strategies


### Within-cycle correction (WCC) using Simpson's 1/3 rule 
v_wcc <- gen_wcc(n_cycles = n_cycles,  method = "Simpson1/3")

### Transition probabilities
p_HD_yr <- 0.01  # annual probability of dying when healthy
p_SD_yr <- 0.10  # annual probability of dying when sick
# Annual probability of becoming sick from Healthy, conditional on surviving cycle
p_HS_yr_SoC   <- 0.05   # under standard of care
p_HS_yr_trtA  <- 0.04  # under treatment A
p_HS_yr_trtB  <- 0.02  # under treatment B


### State rewards

#### Costs 
c_H_yr     <- 400   # cost of one year in healthy state
c_S_yr     <- 1000  # cost of one year in sick state
c_D_yr     <- 0     # cost of one year in dead state
c_trtA_yr  <- 800   # cost of treatment A (per year) in healthy state
c_trtB_yr  <- 1500  # cost of treatment B (per year) in healthy state
#### Utilities
u_H        <- 1     # utility when healthy 
u_S        <- 0.5   # utility when sick
u_D        <- 0     # utility when dead

## 03.2 Calculate internal model parameters

### Cycle-specific discount weight for costs and effects 
v_dwc   <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
v_dwe   <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))

### Calculate cycle-specific transition probabilities
p_HD <- 1 - exp(-(-log(1 - p_HD_yr) * cycle_length)) # probability of dying from Healthy
p_SD <- 1 - exp(-(-log(1 - p_SD_yr) * cycle_length)) # probability of dying from Sick
# probability of becoming sick from Healthy
p_HS_SoC  <- 1 - exp(-(-log(1 - p_HS_yr_SoC) * cycle_length))  # Standard of Care
p_HS_trtA <- 1 - exp(-(-log(1 - p_HS_yr_trtA) * cycle_length)) # Treatment A
p_HS_trtB <- 1 - exp(-(-log(1 - p_HS_yr_trtB) * cycle_length)) # Treatment B

### Calculate cycle-specific state rewards
#### Costs
v_c_SoC  <- c(c_H_yr, c_S_yr, c_D_yr) * cycle_length             # Standard of Care
v_c_trtA <- c(c_H_yr + c_trtA_yr, c_S_yr, c_D_yr) * cycle_length # Treatment A
v_c_trtB <- c(c_H_yr + c_trtB_yr, c_S_yr, c_D_yr) * cycle_length # Treatment B
#### QALYs
v_q_SoC  <- c(u_H, u_S, u_D) * cycle_length # Standord of Care
v_q_trtA <- v_q_trtB <- v_q_SoC             # Treatments A and B have same utilities as SoC

# 04 Construct state-transition models

m_P_diag <- matrix(0, nrow = n_states, ncol = n_states, 
                   dimnames = list(v_names_states, v_names_states))
m_P_diag["Healthy", "Sick" ]     = "" 
m_P_diag["Healthy", "Dead" ]     = ""
m_P_diag["Healthy", "Healthy" ]  = ""
m_P_diag["Sick"   , "Dead" ]     = ""
m_P_diag["Sick"   , "Sick" ]     = ""
m_P_diag["Dead"   , "Dead" ]     = ""
layout.fig <- c(2, 1)
plotmat(t(m_P_diag), t(layout.fig), self.cex = 0.5, curve = 0, arr.pos = 0.8,  
        latex = T, arr.type = "curved", relsize = 0.85, box.prop = 0.8, 
        cex = 0.8, box.cex = 0.7, lwd = 1)


## 04.1 Initial state vector

# All starting healthy
v_m_init <- c("Healthy" = 1, "Sick" = 0, "Dead" = 0)  
v_m_init

## 04.2 Initialize cohort traces

### Initialize cohort trace for SoC 
m_M_SoC <- matrix(0, 
                  nrow = (n_cycles + 1), ncol = n_states, 
                  dimnames = list(v_names_cycles, v_names_states))
# Store the initial state vector in the first row of the cohort trace
m_M_SoC[1, ] <- v_m_init

## Initialize cohort traces for treatments A and B
# Structure and initial states are the same as for SoC
m_M_trtA <- m_M_trtB <- m_M_SoC

## 04.3 Create transition probability arrays

## Create transition probability matrices for strategy SoC 
### Initialize transition probability matrix 
# All transitions to a non-death state are assumed to be conditional on survival 
m_P_SoC  <- matrix(0,
                   nrow = n_states, ncol = n_states,
                   dimnames = list(v_names_states, v_names_states)) 
### Fill in matrix 
# from Healthy
m_P_SoC["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_SoC)
m_P_SoC["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_SoC
m_P_SoC["Healthy", "Dead"]    <-      p_HD
# from Sick
m_P_SoC["Sick", "Sick"] <- 1 - p_SD
m_P_SoC["Sick", "Dead"] <-     p_SD
# from Dead
m_P_SoC["Dead", "Dead"] <- 1

## Treatment A
# Start with same matrix as SoC, but replace parameters that differ for trtA
m_P_trtA <- m_P_SoC 
m_P_trtA["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trtA)
m_P_trtA["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_trtA

## Treatment B
# Start with same matrix as SoC, but replace parameters that differ for trtB
m_P_trtB <- m_P_SoC
m_P_trtB["Healthy", "Healthy"] <- (1 - p_HD) * (1 - p_HS_trtB)
m_P_trtB["Healthy", "Sick"]    <- (1 - p_HD) *      p_HS_trtB

## Check if transition probability matrices are valid
### Check that transition probabilities are in [0, 1]
check_transition_probability(m_P_SoC,  verbose = TRUE)
check_transition_probability(m_P_trtA, verbose = TRUE)
check_transition_probability(m_P_trtB, verbose = TRUE)
### Check that all rows sum to 1
check_sum_of_transition_array(m_P_SoC,  n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(m_P_trtA, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)
check_sum_of_transition_array(m_P_trtB, n_states = n_states, n_cycles = n_cycles, verbose = TRUE)

# 05 Run cohort state transition model

## Loop over time
# Calculating cohort state based on previous state and transition matrix
for (t in 1:n_cycles) {  
  # For SoC
  m_M_SoC[t + 1, ] <- m_M_SoC[t, ] %*% m_P_SoC  
  # For treatment A
  m_M_trtA[t + 1, ] <- m_M_trtA[t, ] %*% m_P_trtA  
  # For treatment B
  m_M_trtB[t + 1, ] <- m_M_trtB[t, ] %*% m_P_trtB 
}



# 06 Plot Outputs

## 06.1 Plot the cohort trace for strategies SoC

plot_trace(m_M_SoC) 

## 06.2 Overall Survival (OS)

#Print the overall survival for the Standard of Care

v_os_SoC <- 1 - m_M_SoC[, "Dead"]    # calculate the overall survival (OS) probability
v_os_SoC <- rowSums(m_M_SoC[, 1:2])  # alternative way of calculating the OS probability   

plot(v_os_SoC, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")  # create a simple plot showing the OS

# add grid 
grid(nx = n_cycles, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), 
     equilogs = TRUE) 

## 06.2.1 Life Expectancy (LE)

le_SoC <- sum(v_os_SoC)  # summing probability of OS over time  (i.e. life expectancy)
le_SoC

## 06.2.2 Disease prevalence

v_prev <- m_M_SoC[, "Sick"]/v_os_SoC
plot(v_prev,
     ylim = c(0, 1),
     ylab = "Prevalence",
     xlab = "Cycle",
     main = "Disease prevalence")

# 07 State Rewards 

## Scale by the cycle length 
# Vector of state utilities under strategy SoC
v_u_SoC    <- c(H  = u_H, 
                S  = u_S,
                D  = u_D) * cycle_length
# Vector of state costs under strategy SoC
v_c_SoC    <- c(H  = c_H_yr, 
                S  = c_S_yr,
                D  = c_D_yr) * cycle_length
# Vector of state utilities under treatment A
v_u_trtA   <- c(H  = u_H, 
                S  = u_S, 
                D  = u_D) * cycle_length
# Vector of state costs under treatment A
v_c_trtA   <- c(H  = c_H_yr + c_trtA_yr, 
                S  = c_S_yr, 
                D  = c_D_yr) * cycle_length
# Vector of state utilities under treatment B
v_u_trtB   <- c(H  = u_H, 
                S  = u_S, 
                D  = u_D) * cycle_length
# Vector of state costs under treatment B
v_c_trtB   <- c(H  = c_H_yr + c_trtB_yr, 
                S  = c_S_yr, 
                D  = c_D_yr) * cycle_length

## Store state rewards 
# Store the vectors of state utilities for each strategy in a list 
l_u   <- list(SQ = v_u_SoC,
              A  = v_u_trtA,
              B  = v_u_trtB)
# Store the vectors of state cost for each strategy in a list 
l_c   <- list(SQ = v_c_SoC,
              A  = v_c_trtA,
              B  = v_c_trtB)

# assign strategy names to matching items in the lists
names(l_u) <- names(l_c) <- v_names_str

# 08 Compute expected outcomes 


# Create empty vectors to store total costs and QALYs
v_tot_cost <- v_tot_qaly <- vector(mode = "numeric", length = n_str)
names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str


### Expected discounted cost for each strategy
v_tot_cost["Standard of Care"] <- t(m_M_SoC %*% v_c_SoC) %*% (v_dwe * v_wcc)
v_tot_cost["Treatment A"] <- t(m_M_trtA %*% v_c_trtA) %*% (v_dwe * v_wcc)
v_tot_cost["Treatment B"] <- t(m_M_trtB %*% v_c_trtB) %*% (v_dwe * v_wcc)

### Expected discounted QALYs for each strategy
v_tot_qaly["Standard of Care"] <- t(m_M_SoC %*% v_q_SoC) %*% (v_dwe * v_wcc)
v_tot_qaly["Treatment A"] <- t(m_M_trtA %*% v_q_trtA) %*% (v_dwe * v_wcc)
v_tot_qaly["Treatment B"] <- t(m_M_trtB %*% v_q_trtB) %*% (v_dwe * v_wcc)

# 09 Cost-effectiveness analysis (CEA) 

## Incremental cost-effectiveness ratios (ICERs) 
df_cea <- dampack::calculate_icers(cost       = v_tot_cost, 
                                   effect     = v_tot_qaly,
                                   strategies = v_names_str)
df_cea

## CEA table in proper format 
table_cea <- format_table_cea(df_cea) 
table_cea

## CEA frontier using dampack
plot(df_cea, label = "all", txtsize = 14) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.3))

# 10 Deterministic Sensitivity Analysis (DSA)

## Load model, CEA and PSA functions 
source('Functions_cSTM_3state.R')

## 10.1 Model input for SA

l_params_all <- list(
    # Transition probabilities
    # probability of dying
  p_HD_yr       = 0.01,  # annual probability of dying when healthy
  p_SD_yr       = 0.10,  # annual probability of dying when sick
  # Annual probability of becoming sick from Healthy, conditional on surviving cycle
  p_HS_yr_SoC   = 0.05,  # under standard of care
  p_HS_yr_trtA  = 0.04,  # under treatment A
  p_HS_yr_trtB  = 0.02,  # under treatment B
  
  ### State rewards
  
  #### Costs 
  c_H_yr     = 400,   # cost of one year in healthy state
  c_S_yr     = 1000,  # cost of one year in sick state
  c_D_yr     = 0,     # cost of one year in dead state
  c_trtA_yr  = 800,   # cost of treatment A (per year) in healthy state
  c_trtB_yr  = 1500,  # cost of treatment B (per year) in healthy state
  #### Utilities
  u_H        = 1,     # utility when healthy 
  u_S        = 0.5,   # utility when sick
  u_D        = 0,     # utility when dead
  # Discount rates
  d_e       = 0.03,  # discount rate per cycle equal discount of costs and QALYs by 3%
  d_c       = 0.03,  # discount rate per cycle equal discount of costs and QALYs by 3%
  # Time horizon
  n_time_horizon_yr = 60,                 # time horizon (in years)
  cycle_length      = 1                   # cycle length in years (use 1/12 for monthly)

)

#**Test model functions**
#
#A function is defined in the `Functions_cSTM_3state.R` file. 
#The first is the `calculate_ce_out()` function which runs the decision model 
#and uses the resulting Markov trace to compute the total costs, QALYs, 
#and net monetary benefit (NMB) for each strategy, returning a data frame of 
#costs and QALYs.
#

# Try the calculate_ce_out() function
df_ce <- calculate_ce_out(l_params_all)
df_ce

# Get strategies names (will be used to label plots)
v_names_str <- df_ce$Strategy
n_str <- length(v_names_str)

## 10.2 One-way sensitivity analysis (OWSA)

options(scipen = 999) # disabling scientific notation in R
# dataframe containing all parameters, their base case values, and the min and 
# max values of the parameters of interest 
df_params_owsa <- data.frame(pars = c("c_trtA_yr", "c_trtB_yr", "c_S_yr"),
                             min  = c(300 , 500,   500), # min parameter values
                             max  = c(1200, 2000, 2000)  # max parameter values
                             )
owsa_nmb  <- run_owsa_det(params_range     = df_params_owsa,   # data.frame with parameters for OWSA
                          params_basecase  = l_params_all,     # list with all parameters
                          nsamp            = 100,              # number of parameter values
                          FUN              = calculate_ce_out, # function to compute outputs
                          outcomes         = c("NMB"),         # output to do the OWSA on
                          strategies       = v_names_str,      # names of the strategies
                          n_wtp            = 5000)             # extra argument to pass to FUN

plot(owsa_nmb, txtsize = 10, n_x_ticks = 4, 
     facet_scales = "free") +
     theme(legend.position = "bottom")

### 10.2.1 Optimal strategy with OWSA

owsa_opt_strat(owsa = owsa_nmb, txtsize = 10)

### 10.2.2 Tornado plot

tornado <- owsa_tornado(owsa = owsa_nmb)
tornado + labs(title = "Tornado plot",
               subtitle = "c_trtA_yr: min = €300, max = €1,200\nc_trtB_yr: min = €500, max = €2,000\nc_S_yr:    min = €500, max = €2,000",
               x = "Outcome (Net Monetary Benefit in €)") +
  theme(plot.subtitle = element_text(size = 12, 
                                     face = "plain", 
                                     color = "black"))

## 10.3 Two-way sensitivity analysis (TWSA)

# data.frame containing all parameters, their base-case values, and the min and 
# max values of the parameters of interest
df_params_twsa <- data.frame(pars = c("c_trtA_yr", "c_trtB_yr"),
                             min  = c(300, 500),  # min parameter values
                             max  = c(1200, 2000) # max parameter values
                             )

twsa_nmb <- run_twsa_det(params_range    = df_params_twsa,    # dataframe with parameters for TWSA
                         params_basecase = l_params_all,      # list with all parameters
                         nsamp           = 40,                # number of parameter values
                         FUN             = calculate_ce_out,  # function to compute outputs
                         outcomes        = "NMB",             # output to do the TWSA on
                         strategies      = v_names_str,       # names of the strategies
                         n_wtp           = 5000)              # extra argument to pass to FUN

### 10.3.1 Plot TWSA

plot(twsa_nmb)

# 11 Probabilistic Sensitivity Analysis (PSA) 

## 11.1 Model input

# Store the parameter names into a vector
v_names_params <- names(l_params_all)

## Test functions to generate CE outcomes and PSA dataset 
# Test function to compute CE outcomes
calculate_ce_out(l_params_all) 

# Test function to generate PSA input dataset
generate_psa_params(10) 

## Generate PSA dataset 
# Number of simulations
n_sim <- 1000

# Generate PSA input dataset
df_psa_input <- generate_psa_params(n_sim = n_sim)
# First six observations
head(df_psa_input)

### Histogram of parameters 
ggplot(melt(df_psa_input, variable.name = "Parameter"), aes(x = value)) +
  facet_wrap(~Parameter, scales = "free") +
  geom_histogram(aes(y = after_stat(density))) +
  ylab("") +
  theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()) 

## 11.2 Run PSA

# Initialize data.frames with PSA output 
# data.frame of costs
df_c <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_c) <- v_names_str
# data.frame of effectiveness
df_e <- as.data.frame(matrix(0, 
                             nrow = n_sim,
                             ncol = n_str))
colnames(df_e) <- v_names_str

# Conduct probabilistic sensitivity analysis
# Run Markov model on each parameter set of PSA input dataset
n_time_init_psa_series <- Sys.time()
for (i in 1:n_sim) { # i <- 1
  l_psa_input <- update_param_list(l_params_all, df_psa_input[i,])
  # Outcomes
  l_out_ce_temp  <- calculate_ce_out(l_psa_input)
  df_c[i, ]  <- l_out_ce_temp$Cost  
  df_e[i, ]  <- l_out_ce_temp$Effect
  # Display simulation progress
  if (i/(n_sim/100) == round(i/(n_sim/100), 0)) { # display progress every 5%
    cat('\r', paste(i/n_sim * 100, "% done", sep = " "))
  }
}
n_time_end_psa_series <- Sys.time()
n_time_total_psa_series <- n_time_end_psa_series - n_time_init_psa_series
print(paste0("PSA with ", scales::comma(n_sim), " simulations run in series in ", 
             round(n_time_total_psa_series, 2), " ", 
             units(n_time_total_psa_series)))

# 11.3 Visualize PSA results for CEA 

### Create PSA object 
l_psa <- dampack::make_psa_obj(cost          = df_c, 
                               effectiveness = df_e, 
                               parameters    = df_psa_input, 
                               strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost) <- v_names_str

# Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 30000, by = 1000)

## 11.3.1 Cost-Effectiveness Scatter plot 

### Cost-Effectiveness Scatter plot 
txtsize <- 13
gg_scattter <- plot_psa(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter

## 11.3.2 Incremental cost-effectiveness ratios (ICERs) with probabilistic output

### Incremental cost-effectiveness ratios (ICERs) with probabilistic output 
# Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)
df_cea_psa <- dampack::calculate_icers(cost       = df_out_ce_psa$meanCost, 
                                       effect     = df_out_ce_psa$meanEffect,
                                       strategies = df_out_ce_psa$Strategy)
df_cea_psa

## 11.3.3 Plot cost-effectiveness frontier with probabilistic output

### Plot cost-effectiveness frontier with probabilistic output 
plot_icers(df_cea_psa, label = "all", txtsize = txtsize) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.3))

## 11.3.4 Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF)

### Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) 
ceac_obj <- dampack::ceac(wtp = v_wtp, psa = l_psa)
# Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
# CEAC & CEAF plot
gg_ceac <- plot_ceac(ceac_obj, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))
gg_ceac

## 11.3.5 Expected Loss Curves (ELCs)

### Expected Loss Curves (ELCs) 
elc_obj <- dampack::calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

# ELC plot
gg_elc <- plot_exp_loss(elc_obj, log_y = FALSE, 
                        txtsize = txtsize, 
                        xlim = c(0, NA), 
                        n_x_ticks = 14,
                        col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7),)
gg_elc

## 11.3.6 Expected value of perfect information (EVPI) 

### Expected value of perfect information (EVPI) 
evpi <- dampack::calc_evpi(wtp = v_wtp, psa = l_psa)
# EVPI plot
gg_evpi <- plot_evpi(evpi, effect_units = "QALY", 
                     txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     breaks = number_ticks(10),
                     labels = function(x) x/1000)
gg_evpi