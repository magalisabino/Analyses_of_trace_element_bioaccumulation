##--------------------------------------------------------------------------------------------------------
## SCRIPT : Process of model selection for mathematical correction of d13C values in
##          non-lipid-free swordfish samples.
##
## As part of :
##        Magali SABINO PhD "Bioaccumulation of trace elements in Seychelles marine food webs"
##        supervised by Paco BUSTAMANTE and Nathalie BODIN
##
## Author : Magali Sabino
## First created : 2022-01-13
## Last update : 2022-01-13
##
##
## For more information on the tested models and on the method, see:
##    - Fry (2002) Stable isotopic indicators of habitat use by Mississippi River fish. J North Am Benthol Soc 21:676–685. doi: 10.2307/1468438
##    - Logan et al. (2008) Lipid corrections in carbon and nitrogen stable isotope analyses: comparison of chemical extraction and modelling methods. J Anim Ecol 77:838–846. doi: 10.1111/j.1365-2656.2008.01394.x
##    - McConnaughey and McRoy (1979) Food-Web structure and the fractionation of Carbon isotopes in the Bering Sea. Mar Biol 53:257–262. doi: 10.1007/BF00952434
##    - Sardenne et al. (2015) Methods of lipid-normalization for multi-tissue stable isotope analyses in tropical tuna. Rapid Commun Mass Spectrom 29:1253-1267. doi: 10.1002/rcm.7215
##    - Sabino (2021) Appendix 2.3. in PhD manuscript "Bioaccumulation of trace elements in Seychelles marine food webs" pp. 230-232
##
####
##
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## Copyright (C) 2021 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------


### 0 // Packages ##########################################################################################

## Open libraries
lapply(c("tidyverse", "openxlsx", "nlstools", "Metrics"),
       library, character.only = TRUE)

## Creating function to calculate standard deviation
std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))



### I // Creation of database with SI values for lipid-free and non-lipid-free samples ##############################################################################################

## 1 / Dataset importation

SWO_SI_norm <- read.xlsx("C:/Users/msabin02/Nextcloud/PhD/R_projects/Final_scripts/SWO_SI_data.xlsx")

## 2 / Dataframe manipulation prior to data  analyses

SI_liprem <- SWO_SI_norm %>% 
  filter(lipid_removed == "yes") %>% 
  rename(d13C_delip = d13C,
         C_delip = C,
         d15N_delip = d15N,
         N_delip = N,
         C_N_delip = C_N) %>% 
  select(-lipid_removed)

SI_lipnotrem <- SWO_SI_norm %>% 
  filter(lipid_removed == "no",
         sample_identifier %in% SI_liprem$sample_identifier) %>% 
  rename(d13C_bulk = d13C,
         C_bulk = C,
         d15N_bulk = d15N,
         N_bulk = N,
         C_N_bulk = C_N) %>% 
  select(-lipid_removed)

SI_norm_data <- SI_liprem %>% 
  left_join(SI_lipnotrem, by = "sample_identifier") %>% 
  mutate(delta_d13C= d13C_delip-d13C_bulk)

rm(SI_liprem,SI_lipnotrem)



### II // Testing models to determine which is best to correct d13C values ##############################################################################################

# Five different models were tested in this method :
#   - The McConnaughey and McRoy (1979) model (Eq 1)
#   - The Logan et al. (2008) model, derived from the McConnaughey and McRoy (1979) (Eq 2)
#   - The Fry (2002) model (Eq 3)
#   - A log-linear model (Eq 4)
#   - A classic linear model (Eq 5)
#
# Cross-validation was used to assess the predictive performance of each model, with the following steps :
#   - Step 1: model was trained with a random subset of the data
#   - Step 2: optimal parameters were calculated for this subset using the nls function of the ntools package.
#     For Eq 1 to 4, starting values were found in Logan et al. (2008); for the linear model (Eq 5), a linear
#     regression was first modelled on all the dataset.
#   - Step 3: d13C values were calculated from the validation dataset and compared to measured d13C lipid-free
#     using paired wilcoxon tests.
#   - Step 4: mean squared error (MSE) and mean absolute error (MAE) were calculated using the Metrics package
# These steps were repeated 500 times for each model, and then the percentage of good prediction was calculated.
#
# The selected best model was the one with the lowest MSE, the lowest MAE and the highest percentage of good
# predictions.


## 1 / Equation 1 - McConnaughey and McRoy (1979) model

formula_eq1 <- as.formula(delta_d13C ~ D * (theta + (3.90/(1 + 287 / (93/(1 + (0.246 * C_N_bulk - 0.775)^-1))))))
CV_results_eq1 <- NULL

for (i in 1:500) {
  #i = 1
  
  # STEP 1: Database creation
  CV_results <- data.frame(iter = i,
                           D = NA,
                           theta = NA,
                           mean_d13C_valid = NA,
                           mean_d13C_predict = NA,
                           p_wilcox = NA,
                           mean_delta_d13C = NA,
                           MSE = NA,
                           MAE = NA)
  
  # STEP 2: Data selection
  SI_train <- SI_norm_data %>% 
    sample_n(22)
  SI_valid <- SI_norm_data %>% 
    filter(!sample_identifier %in% SI_train$sample_identifier)
  
  # STEP 3: Calculation of optimal parameter on train dataset
  opt.par <- nls(formula_eq1, data = SI_train, start = list(D = 7.49, theta = 0.015))
  opt.par <- summary(opt.par)
  D <- opt.par$parameters[1,1]
  theta <- opt.par$parameters[2,1]
  
  CV_results$D <- D
  CV_results$theta <- theta
  
  # STEP 4: Calculation of d13C corrected
  SI_valid <- SI_valid %>% 
    mutate(d13C_corr = d13C_bulk + (D * (theta + (3.90/(1 + 287 / (93/(1 + (0.246 * C_N_bulk - 0.775)^-1)))))))
  
  CV_results$mean_d13C_valid <- mean(SI_valid$d13C_delip, na.rm = TRUE)
  CV_results$mean_d13C_predict <- mean(SI_valid$d13C_corr, na.rm = TRUE)
  
  # STEP 5: Test difference between lipid-free and corrected
  p_wilcox <- wilcox.test(SI_valid$d13C_delip,SI_valid$d13C_corr, paired = TRUE)
  CV_results$p_wilcox <- p_wilcox$p.value
  
  delta_d13C <- SI_valid$d13C_corr-SI_valid$d13C_delip
  CV_results$mean_delta_d13C <- mean(delta_d13C, na.rm = TRUE)
  
  CV_results$MSE <- Metrics::rmse(SI_valid$d13C_delip,SI_valid$d13C_corr)
  CV_results$MAE <- Metrics::mae(SI_valid$d13C_delip,SI_valid$d13C_corr)
  
  # STEP 6: Merge results
  CV_results_eq1 <- rbind(CV_results_eq1, CV_results)
  
  # rm variables
  rm(CV_results)
  rm(SI_train,SI_valid)
  rm(opt.par,D,theta)
  rm(p_wilcox)
  rm(delta_d13C)
}
rm(i)

CV_eq1 <- CV_results_eq1 %>% 
  summarise(mean_D = mean(D, na.rm = TRUE),
            SE_D = std(D),
            mean_theta = mean(theta, na.rm = TRUE),
            SE_theta = std(theta),
            mean_mdelta_d13C = mean(mean_delta_d13C, na.rm = TRUE),
            SE_delta_d13C = std(mean_delta_d13C),
            mean_p_wilcox = mean(p_wilcox, na.rm = TRUE),
            SE_p_wilcox = std(p_wilcox),
            mean_MSE = mean(MSE, na.rm = TRUE),
            SE_MSE = std(MSE),
            mean_MAE = mean(MAE, na.rm = TRUE),
            SE_MAE = std(MAE)) %>% 
  mutate(mean_D = round(mean_D, 2),
         SE_D = round(SE_D, 2),
         mean_theta = round(mean_theta, 4),
         SE_theta = round(SE_theta, 4),
         mean_mdelta_d13C = round(mean_mdelta_d13C, 4),
         SE_delta_d13C = round(SE_delta_d13C, 4),
         mean_p_wilcox = round(mean_p_wilcox, 2),
         SE_p_wilcox = round(SE_p_wilcox, 2),
         mean_MSE = round(mean_MSE, 3),
         SE_MSE = round(SE_MSE, 3),
         mean_MAE = round(mean_MAE, 3),
         SE_MAE = round(SE_MAE, 3)) %>% 
  mutate(model = "eq1",
         D = paste0(mean_D," ± ",SE_D),
         theta = paste0(mean_theta," ± ",SE_theta),
         delta_d13C = paste0(mean_mdelta_d13C," ± ",SE_delta_d13C),
         p_wilcox = paste0(mean_p_wilcox," ± ",SE_p_wilcox),
         MSE = paste0(mean_MSE," ± ",SE_MSE),
         MAE = paste0(mean_MAE," ± ",SE_MAE)) %>% 
  select(model,D,theta,delta_d13C,p_wilcox,MSE,MAE)

# Calculation % good prediction
numb_good_pred <- length(which(CV_results_eq1$p_wilcox >= 0.05))
CV_eq1$percent_good_pred <- numb_good_pred*100/500

rm(formula_eq1, CV_results_eq1, numb_good_pred)


## 2 / Equation 2 - Logan et al. (2008) model

formula_eq2 <- as.formula(delta_d13C ~ (D * C_N_bulk + a)/(C_N_bulk + b))
CV_results_eq2 <- NULL

for (i in 1:500) {
  #i = 1
  
  # STEP 1: Database creation
  CV_results <- data.frame(iter = i,
                           D = NA,
                           a = NA,
                           b = NA,
                           mean_d13C_valid = NA,
                           mean_d13C_predict = NA,
                           p_wilcox = NA,
                           mean_delta_d13C = NA,
                           MSE = NA,
                           MAE = NA)
  
  # STEP 2: Data selection
  SI_train <- SI_norm_data %>% 
    sample_n(22)
  SI_valid <- SI_norm_data %>% 
    filter(!sample_identifier %in% SI_train$sample_identifier)
  
  # STEP 3: Calculation of optimal parameter on train dataset
  opt.par <- nls(formula_eq2, data = SI_norm_data, start = list(D = 7.415, a = -22.732, b = 0.746))
  opt.par <- summary(opt.par)
  D <- opt.par$parameters[1,1]
  a <- opt.par$parameters[2,1]
  b <- opt.par$parameters[3,1]
  
  CV_results$D <- D
  CV_results$a <- a
  CV_results$b <- b
  
  # STEP 4: Calculation of d13C corrected
  SI_valid <- SI_valid %>% 
    mutate(d13C_corr = d13C_bulk + ((D * C_N_bulk + a)/(C_N_bulk + b)))
  
  CV_results$mean_d13C_valid <- mean(SI_valid$d13C_delip, na.rm = TRUE)
  CV_results$mean_d13C_predict <- mean(SI_valid$d13C_corr, na.rm = TRUE)
  
  # STEP 5: Test difference between lipid-free and corrected
  p_wilcox <- wilcox.test(SI_valid$d13C_delip,SI_valid$d13C_corr, paired = TRUE)
  CV_results$p_wilcox <- p_wilcox$p.value
  
  delta_d13C <- SI_valid$d13C_corr-SI_valid$d13C_delip
  CV_results$mean_delta_d13C <- mean(delta_d13C, na.rm = TRUE)
  
  CV_results$MSE <- Metrics::rmse(SI_valid$d13C_delip,SI_valid$d13C_corr)
  CV_results$MAE <- Metrics::mae(SI_valid$d13C_delip,SI_valid$d13C_corr)
  
  # STEP 6: Merge results
  CV_results_eq2 <- rbind(CV_results_eq2, CV_results)
  
  # rm variables
  rm(CV_results)
  rm(SI_train,SI_valid)
  rm(opt.par,D,a,b)
  rm(p_wilcox)
  rm(delta_d13C)
}
rm(i)

CV_eq2 <- CV_results_eq2 %>% 
  summarise(mean_D = mean(D, na.rm = TRUE),
            SE_D = std(D),
            mean_a = mean(a, na.rm = TRUE),
            SE_a = std(a),
            mean_b = mean(b, na.rm = TRUE),
            SE_b = std(b),
            mean_mdelta_d13C = mean(mean_delta_d13C, na.rm = TRUE),
            SE_delta_d13C = std(mean_delta_d13C),
            mean_p_wilcox = mean(p_wilcox, na.rm = TRUE),
            SE_p_wilcox = std(p_wilcox),
            mean_MSE = mean(MSE, na.rm = TRUE),
            SE_MSE = std(MSE),
            mean_MAE = mean(MAE, na.rm = TRUE),
            SE_MAE = std(MAE)) %>% 
  mutate(mean_D = round(mean_D, 2),
         SE_D = round(SE_D, 2),
         mean_a = round(mean_a, 2),
         SE_a = round(SE_a, 2),
         mean_b = round(mean_b, 2),
         SE_b = round(SE_b, 2),
         mean_mdelta_d13C = round(mean_mdelta_d13C, 4),
         SE_delta_d13C = round(SE_delta_d13C, 4),
         mean_p_wilcox = round(mean_p_wilcox, 2),
         SE_p_wilcox = round(SE_p_wilcox, 2),
         mean_MSE = round(mean_MSE, 3),
         SE_MSE = round(SE_MSE, 3),
         mean_MAE = round(mean_MAE, 3),
         SE_MAE = round(SE_MAE, 3)) %>% 
  mutate(model = "eq2",
         D = paste0(mean_D," ± ",SE_D),
         a = paste0(mean_a," ± ",SE_a),
         b = paste0(mean_b," ± ",SE_b),
         delta_d13C = paste0(mean_mdelta_d13C," ± ",SE_delta_d13C),
         p_wilcox = paste0(mean_p_wilcox," ± ",SE_p_wilcox),
         MSE = paste0(mean_MSE," ± ",SE_MSE),
         MAE = paste0(mean_MAE," ± ",SE_MAE)) %>% 
  select(model,D,a,b,delta_d13C,p_wilcox,MSE,MAE)

# Calculation % good prediction
numb_good_pred <- length(which(CV_results_eq2$p_wilcox >= 0.05))
CV_eq2$percent_good_pred <- numb_good_pred*100/500

rm(formula_eq2, CV_results_eq2, numb_good_pred)


## 3 / Equation 3 - Fry (2002) model

formula_eq3 <- as.formula(delta_d13C ~ P - ((P * eF)/C_N_bulk))
CV_results_eq3 <- NULL

for (i in 1:500) {
  #i = 1
  
  # STEP 1: Database creation
  CV_results <- data.frame(iter = i,
                           P = NA,
                           eF = NA,
                           mean_d13C_valid = NA,
                           mean_d13C_predict = NA,
                           p_wilcox = NA,
                           mean_delta_d13C = NA,
                           MSE = NA,
                           MAE = NA)
  
  # STEP 2: Data selection
  SI_train <- SI_norm_data %>% 
    sample_n(22)
  SI_valid <- SI_norm_data %>% 
    filter(!sample_identifier %in% SI_train$sample_identifier)
  
  # STEP 3: Calculation of optimal parameter on train dataset
  opt.par <- nls(formula_eq3, data = SI_norm_data, start = list(P = 6.699, eF = 3.098))
  opt.par <- summary(opt.par)
  P <- opt.par$parameters[1,1]
  eF <- opt.par$parameters[2,1]
  
  CV_results$P <- P
  CV_results$eF <- eF
  
  # STEP 4: Calculation of d13C corrected
  SI_valid <- SI_valid %>% 
    mutate(d13C_corr = d13C_bulk + (P - ((P * eF)/C_N_bulk)))
  
  CV_results$mean_d13C_valid <- mean(SI_valid$d13C_delip, na.rm = TRUE)
  CV_results$mean_d13C_predict <- mean(SI_valid$d13C_corr, na.rm = TRUE)
  
  # STEP 5: Test difference between lipid-free and corrected
  p_wilcox <- wilcox.test(SI_valid$d13C_delip,SI_valid$d13C_corr, paired = TRUE)
  CV_results$p_wilcox <- p_wilcox$p.value
  
  delta_d13C <- SI_valid$d13C_corr-SI_valid$d13C_delip
  CV_results$mean_delta_d13C <- mean(delta_d13C, na.rm = TRUE)
  
  CV_results$MSE <- Metrics::rmse(SI_valid$d13C_delip,SI_valid$d13C_corr)
  CV_results$MAE <- Metrics::mae(SI_valid$d13C_delip,SI_valid$d13C_corr)
  
  # STEP 6: Merge results
  CV_results_eq3 <- rbind(CV_results_eq3, CV_results)
  
  # rm variables
  rm(CV_results)
  rm(SI_train,SI_valid)
  rm(opt.par,P,eF)
  rm(p_wilcox)
  rm(delta_d13C)
}
rm(i)

CV_eq3 <- CV_results_eq3 %>% 
  summarise(mean_P = mean(P, na.rm = TRUE),
            SE_P = std(P),
            mean_eF = mean(eF, na.rm = TRUE),
            SE_eF = std(eF),
            mean_mdelta_d13C = mean(mean_delta_d13C, na.rm = TRUE),
            SE_delta_d13C = std(mean_delta_d13C),
            mean_p_wilcox = mean(p_wilcox, na.rm = TRUE),
            SE_p_wilcox = std(p_wilcox),
            mean_MSE = mean(MSE, na.rm = TRUE),
            SE_MSE = std(MSE),
            mean_MAE = mean(MAE, na.rm = TRUE),
            SE_MAE = std(MAE)) %>% 
  mutate(mean_P = round(mean_P, 2),
         SE_P = round(SE_P, 2),
         mean_eF = round(mean_eF, 2),
         SE_eF = round(SE_eF, 2),
         mean_mdelta_d13C = round(mean_mdelta_d13C, 4),
         SE_delta_d13C = round(SE_delta_d13C, 4),
         mean_p_wilcox = round(mean_p_wilcox, 2),
         SE_p_wilcox = round(SE_p_wilcox, 2),
         mean_MSE = round(mean_MSE, 3),
         SE_MSE = round(SE_MSE, 3),
         mean_MAE = round(mean_MAE, 3),
         SE_MAE = round(SE_MAE, 3)) %>% 
  mutate(model = "eq3",
         P = paste0(mean_P," ± ",SE_P),
         eF = paste0(mean_eF," ± ",SE_eF),
         delta_d13C = paste0(mean_mdelta_d13C," ± ",SE_delta_d13C),
         p_wilcox = paste0(mean_p_wilcox," ± ",SE_p_wilcox),
         MSE = paste0(mean_MSE," ± ",SE_MSE),
         MAE = paste0(mean_MAE," ± ",SE_MAE)) %>% 
  select(model,P,eF,delta_d13C,p_wilcox,MSE,MAE)

# Calculation % good prediction
numb_good_pred <- length(which(CV_results_eq3$p_wilcox >= 0.05))
CV_eq3$percent_good_pred <- numb_good_pred*100/500

rm(formula_eq3, CV_results_eq3, numb_good_pred)


## 4 / Equation 4 - log-linear model

formula_eq4 <- as.formula(delta_d13C ~ a + b * log(C_N_bulk))
CV_results_eq4 <- NULL

for (i in 1:500) {
  #i = 1
  
  # STEP 1: Database creation
  CV_results <- data.frame(iter = i,
                           a = NA,
                           b = NA,
                           mean_d13C_valid = NA,
                           mean_d13C_predict = NA,
                           p_wilcox = NA,
                           mean_delta_d13C = NA,
                           MSE = NA,
                           MAE = NA)
  
  # STEP 2: Data selection
  SI_train <- SI_norm_data %>% 
    sample_n(22)
  SI_valid <- SI_norm_data %>% 
    filter(!sample_identifier %in% SI_train$sample_identifier)
  
  # STEP 3: Calculation of optimal parameter on train dataset
  opt.par <- nls(formula_eq4, data = SI_norm_data, start = list(a = -4.763, b = 4.401))
  opt.par <- summary(opt.par)
  a <- opt.par$parameters[1,1]
  b <- opt.par$parameters[2,1]
  
  CV_results$a <- a
  CV_results$b <- b
  
  # STEP 4: Calculation of d13C corrected
  SI_valid <- SI_valid %>% 
    mutate(d13C_corr = d13C_bulk + (a + b * log(C_N_bulk)))
  
  CV_results$mean_d13C_valid <- mean(SI_valid$d13C_delip, na.rm = TRUE)
  CV_results$mean_d13C_predict <- mean(SI_valid$d13C_corr, na.rm = TRUE)
  
  # STEP 5: Test difference between lipid-free and corrected
  p_wilcox <- wilcox.test(SI_valid$d13C_delip,SI_valid$d13C_corr, paired = TRUE)
  CV_results$p_wilcox <- p_wilcox$p.value
  
  delta_d13C <- SI_valid$d13C_corr-SI_valid$d13C_delip
  CV_results$mean_delta_d13C <- mean(delta_d13C, na.rm = TRUE)
  
  CV_results$MSE <- Metrics::rmse(SI_valid$d13C_delip,SI_valid$d13C_corr)
  CV_results$MAE <- Metrics::mae(SI_valid$d13C_delip,SI_valid$d13C_corr)
  
  # STEP 6: Merge results
  CV_results_eq4 <- rbind(CV_results_eq4, CV_results)
  
  # rm variables
  rm(CV_results)
  rm(SI_train,SI_valid)
  rm(opt.par,a,b)
  rm(p_wilcox)
  rm(delta_d13C)
}
rm(i)

CV_eq4 <- CV_results_eq4 %>% 
  summarise(mean_a = mean(a, na.rm = TRUE),
            SE_a = std(a),
            mean_b = mean(b, na.rm = TRUE),
            SE_b = std(b),
            mean_mdelta_d13C = mean(mean_delta_d13C, na.rm = TRUE),
            SE_delta_d13C = std(mean_delta_d13C),
            mean_p_wilcox = mean(p_wilcox, na.rm = TRUE),
            SE_p_wilcox = std(p_wilcox),
            mean_MSE = mean(MSE, na.rm = TRUE),
            SE_MSE = std(MSE),
            mean_MAE = mean(MAE, na.rm = TRUE),
            SE_MAE = std(MAE)) %>% 
  mutate(mean_a = round(mean_a, 2),
         SE_a = round(SE_a, 2),
         mean_b = round(mean_b, 2),
         SE_b = round(SE_b, 2),
         mean_mdelta_d13C = round(mean_mdelta_d13C, 4),
         SE_delta_d13C = round(SE_delta_d13C, 4),
         mean_p_wilcox = round(mean_p_wilcox, 2),
         SE_p_wilcox = round(SE_p_wilcox, 2),
         mean_MSE = round(mean_MSE, 3),
         SE_MSE = round(SE_MSE, 3),
         mean_MAE = round(mean_MAE, 3),
         SE_MAE = round(SE_MAE, 3)) %>% 
  mutate(model = "eq4",
         a = paste0(mean_a," ± ",SE_a),
         b = paste0(mean_b," ± ",SE_b),
         delta_d13C = paste0(mean_mdelta_d13C," ± ",SE_delta_d13C),
         p_wilcox = paste0(mean_p_wilcox," ± ",SE_p_wilcox),
         MSE = paste0(mean_MSE," ± ",SE_MSE),
         MAE = paste0(mean_MAE," ± ",SE_MAE)) %>% 
  select(model,a,b,delta_d13C,p_wilcox,MSE,MAE)

# Calculation % good prediction
numb_good_pred <- length(which(CV_results_eq4$p_wilcox >= 0.05))
CV_eq4$percent_good_pred <- numb_good_pred*100/500

rm(formula_eq4, CV_results_eq4, numb_good_pred)


## 5 / Equation 5 - Classic linear model

formula_eq5 <- as.formula(delta_d13C ~ a + b * C_N_bulk)
CV_results_eq5 <- NULL

# Linear regression to determine the starting values to give to the nls function
lin.reg <- lm(SI_norm_data$delta_d13C ~ SI_norm_data$C_N_bulk)
summary(lin.reg)
# a = -0.59729 et b = 0.65 ; adjusted r-squared = 0.7765
rm(lin.reg)

for (i in 1:500) {
  #i = 1
  
  # STEP 1: Database creation
  CV_results <- data.frame(iter = i,
                           a = NA,
                           b = NA,
                           mean_d13C_valid = NA,
                           mean_d13C_predict = NA,
                           p_wilcox = NA,
                           mean_delta_d13C = NA,
                           MSE = NA,
                           MAE = NA)
  
  # STEP 2: Data selection
  SI_train <- SI_norm_data %>% 
    sample_n(22)
  SI_valid <- SI_norm_data %>% 
    filter(!sample_identifier %in% SI_train$sample_identifier)
  
  # STEP 3: Calculation of optimal parameter on train dataset
  opt.par <- nls(formula_eq5, data = SI_norm_data, start = list(a = -0.59729, b = 0.65))
  opt.par <- summary(opt.par)
  a <- opt.par$parameters[1,1]
  b <- opt.par$parameters[2,1]
  
  CV_results$a <- a
  CV_results$b <- b
  
  # STEP 4: Calculation of d13C corrected
  SI_valid <- SI_valid %>% 
    mutate(d13C_corr = d13C_bulk + (a + b * C_N_bulk))
  
  CV_results$mean_d13C_valid <- mean(SI_valid$d13C_delip, na.rm = TRUE)
  CV_results$mean_d13C_predict <- mean(SI_valid$d13C_corr, na.rm = TRUE)
  
  # STEP 5: Test difference between lipid-free and corrected
  p_wilcox <- wilcox.test(SI_valid$d13C_delip,SI_valid$d13C_corr, paired = TRUE)
  CV_results$p_wilcox <- p_wilcox$p.value
  
  delta_d13C <- SI_valid$d13C_corr-SI_valid$d13C_delip
  CV_results$mean_delta_d13C <- mean(delta_d13C, na.rm = TRUE)
  
  CV_results$MSE <- Metrics::rmse(SI_valid$d13C_delip,SI_valid$d13C_corr)
  CV_results$MAE <- Metrics::mae(SI_valid$d13C_delip,SI_valid$d13C_corr)
  
  # STEP 6: Merge results
  CV_results_eq5 <- rbind(CV_results_eq5, CV_results)
  
  # rm variables
  rm(CV_results)
  rm(SI_train,SI_valid)
  rm(opt.par,a,b)
  rm(p_wilcox)
  rm(delta_d13C)
}
rm(i)

CV_eq5 <- CV_results_eq5 %>% 
  summarise(mean_a = mean(a, na.rm = TRUE),
            SE_a = std(a),
            mean_b = mean(b, na.rm = TRUE),
            SE_b = std(b),
            mean_mdelta_d13C = mean(mean_delta_d13C, na.rm = TRUE),
            SE_delta_d13C = std(mean_delta_d13C),
            mean_p_wilcox = mean(p_wilcox, na.rm = TRUE),
            SE_p_wilcox = std(p_wilcox),
            mean_MSE = mean(MSE, na.rm = TRUE),
            SE_MSE = std(MSE),
            mean_MAE = mean(MAE, na.rm = TRUE),
            SE_MAE = std(MAE)) %>% 
  mutate(mean_a = round(mean_a, 2),
         SE_a = round(SE_a, 2),
         mean_b = round(mean_b, 2),
         SE_b = round(SE_b, 2),
         mean_mdelta_d13C = round(mean_mdelta_d13C, 4),
         SE_delta_d13C = round(SE_delta_d13C, 4),
         mean_p_wilcox = round(mean_p_wilcox, 2),
         SE_p_wilcox = round(SE_p_wilcox, 2),
         mean_MSE = round(mean_MSE, 3),
         SE_MSE = round(SE_MSE, 3),
         mean_MAE = round(mean_MAE, 3),
         SE_MAE = round(SE_MAE, 3)) %>% 
  mutate(model = "eq5",
         a = paste0(mean_a," ± ",SE_a),
         b = paste0(mean_b," ± ",SE_b),
         delta_d13C = paste0(mean_mdelta_d13C," ± ",SE_delta_d13C),
         p_wilcox = paste0(mean_p_wilcox," ± ",SE_p_wilcox),
         MSE = paste0(mean_MSE," ± ",SE_MSE),
         MAE = paste0(mean_MAE," ± ",SE_MAE)) %>% 
  select(model,a,b,delta_d13C,p_wilcox,MSE,MAE)

# Calculation % good prediction
numb_good_pred <- length(which(CV_results_eq5$p_wilcox >= 0.05))
CV_eq5$percent_good_pred <- numb_good_pred*100/500

rm(formula_eq5, CV_results_eq5, numb_good_pred)


## 6 / Plot all curves

SI_norm_data %>% 
  ggplot()+
  geom_point(aes(x = C_N_bulk, y = delta_d13C))+
  stat_function(aes(color = "green"), fun = function(C_N_bulk) 7.69 * (0.007 + (3.90/(1 + 287 / (93/(1 + (0.246 * C_N_bulk - 0.775)^-1))))))+ #eq1
  stat_function(aes(color = "red"),fun = function(C_N_bulk) (7.05 * C_N_bulk + -22.4)/(C_N_bulk + -0.44))+ #eq2
  stat_function(aes(color = "blue"),fun = function(C_N_bulk) 7.46 - ((7.46 * 3.12)/C_N_bulk))+ #eq3
  stat_function(aes(color = "pink"), fun = function(C_N_bulk) -4.11 + 4.23 * log(C_N_bulk))+ #eq4
  stat_function(aes(color = "orange"), fun = function(C_N_bulk) -0.6 + 0.65 * C_N_bulk)+ #eq5
  labs(x = "Bulk C:N", y = "d13C lipid-free - d13C bulk")+
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"))


rm(SWO_SI_norm, SI_norm_data, CV_eq1, CV_eq2, CV_eq3, CV_eq4, CV_eq5)
