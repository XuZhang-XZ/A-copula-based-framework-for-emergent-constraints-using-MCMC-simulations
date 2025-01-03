

# Obs Data -------------------------------------------------------------

## Data
Final_EC = read.csv("Output Data/CMIP6/tas/Final.series_EC.csv",row.names = 1)
Final_Copula = read.csv("Output Data/CMIP6/tas/Final_Copula.csv",row.names = 1)
Obs_IVAR_Summary = read.csv("Output Data/CMIP6_EC/PMIP.files_IVAR/Obs_IVAR_Summary.csv",row.names = 1)
Univariate.best = read.csv("Output Data/CMIP6/tas/Univariate.best.csv",row.names = 1)
set.seed(1)

## Correlation
cor(Final_EC$x_1_p,Final_EC$y_1_p,method = "kendall")

## All series and Covariance
All.covariance = array(
  c(
    Obs_IVAR_Summary$Var[1],
    Obs_IVAR_Summary$Cov[1],
    Obs_IVAR_Summary$Cov[2],
    Obs_IVAR_Summary$Var[2]
  ),
  c(2, 2)
)

# 3-d Copula ------------------------------------------------------

## Copula
Copula.family = expand.grid(
  "Copula.Num_1" = c(1,4:6),
  "Copula.Num_2" = c(1,4:6),
  "Copula.Num_3" = c(1,4:6)
)
Copula.paras = Copula.family %>% 
  dplyr::mutate(Var_1 = 1, Var_2 = 3, Var_3 = 2) %>% 
  dplyr::mutate(AIC = NA)
Initial.par = c(0.5,0.5,1.5,1.5,1.5,1.5)

## Loop
i = 1

## Results
for(i in 1:NROW(Copula.paras)) {
  
  ## Probabilities
  u <- as.matrix(Final_EC[,c("x_1_p","x_2_p","y_1_p")])
  
  ## Vine Copula
  d = 3
  order = c(Copula.paras$Var_1[i],Copula.paras$Var_2[i],Copula.paras$Var_3[i])
  family = c(Copula.paras$Copula.Num_1[i],Copula.paras$Copula.Num_2[i],Copula.paras$Copula.Num_3[i])
  RVM <- D2RVine(order = order, family = family, par = Initial.par[family], par2 = rep(NA, d*(d-1)/2))
  
  ## MLE
  temp = RVineMLE(u, RVM, grad = TRUE)
  RVM_Fitted = temp$RVM
  
  ## Copula
  Copula.paras[i,"AIC"] = RVM_Fitted$AIC
  
  ## Print
  print(i)
}

## Select the best Copula with minimum AIC
i = which.min(Copula.paras$AIC)

## Probabilities
u <- as.matrix(Final_EC[,c("x_1_p","x_2_p","y_1_p")])

## Vine Copula
d = 3
order = c(Copula.paras$Var_1[i],Copula.paras$Var_2[i],Copula.paras$Var_3[i])
family = c(Copula.paras$Copula.Num_1[i],Copula.paras$Copula.Num_2[i],Copula.paras$Copula.Num_3[i])
RVM <- D2RVine(order = order, family = family, par = Initial.par[family], par2 = rep(NA, d*(d-1)/2))

## MLE
# plot(RVM_Fitted)
temp = RVineMLE(u, RVM, grad = TRUE)
RVM_Fitted = temp$RVM
RVM_Fitted

## Structure of Vine Copula
# C-vine copula with the following pair-copulas:
#   Tree 1:
#   3,2  Frank (par = 4.19, tau = 0.4) 
#   1,3  Gaussian (par = 0.58, tau = 0.4) 
# 
#   Tree 2:
#   1,2;3  Frank (par = 1.42, tau = 0.15) 


# MCMC ------------------------------------------------------------

## MCMC
MCMC_CEC_Vine = data.frame("y_1" = NA, "CEC" = NA)

## EC for TAS.var.His
beta = c(0.3,0.1,4)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log = function(beta) {
  
  ## Marginal Distribution
  x_1_pre = beta[1]
  x_2_pre = beta[2]
  y_1_pre = beta[3]
  x_1_obs = Obs[1]
  x_2_obs = Obs[2]
  
  ## Fitted Density
  x_1_d = Marignal.distribution(x_1_pre,Type = "d", Num = "x_1")
  x_2_d = Marignal.distribution(x_2_pre,Type = "d", Num = "x_2")
  y_1_d = Marignal.distribution(y_1_pre,Type = "d", Num = "y_1")
  
  ## Fitted Probability
  x_1_p = Marignal.distribution(x_1_pre,Type = "p", Num = "x_1")
  x_2_p = Marignal.distribution(x_2_pre,Type = "p", Num = "x_2")
  y_1_p = Marignal.distribution(y_1_pre,Type = "p", Num = "y_1")
  
  ## Copula Function
  u <- c(x_1_p,x_2_p,y_1_p)
  Copula_PDF = RVinePDF(u, RVM_Fitted, verbose = FALSE)
  prior = Copula_PDF * x_1_d * x_2_d * y_1_d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(x_1_obs,x_2_obs),
                    log = FALSE,
                    mean = c(x_1_pre,x_2_pre),
                    sigma = All.covariance)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])
MCMC_Vine_1 = data.frame("y_1" = Postrior_Re,
                         "CEC" = "CEC")
MCMC_CEC_Vine[1:30000,] = MCMC_Vine_1

## Write
write.csv(MCMC_CEC_Vine,"Output Data/CMIP6/tas/MCMC_CEC_Vine.csv")

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Summary 
mean(samp.1$samples[-c(1:8000), 3])
var(samp.1$samples[-c(1:8000), 3])
cor(samp.1$samples[-c(1:8000), 1],samp.1$samples[-c(1:8000), 2])

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)
gel_2 = gelman.plot(sample.list)

## Convergence
ts = gel_2$shrink
Shrink_factor_Vine = data.frame(
  last.iter = gel_2$last.iter,
  Factor_1 = ts[, 1, 1],
  Factor_1_Upper = ts[, 1, 2],
  Factor_2 = ts[, 2, 1],
  Factor_2_Upper = ts[, 2, 2],
  Factor_3 = ts[, 3, 1],
  Factor_3_Upper = ts[, 3, 2]
)
write.csv(Shrink_factor_Vine,"Output Data/CMIP6_EC/EC.summary/Shrink_factor_Vine.csv")

## Test
plot(density(samp.1$samples[5000:10000,1]))
lines(density(samp.2$samples[5000:10000,1]))
lines(density(samp.3$samples[5000:10000,1]))

## Test
plot(density(samp.1$samples[5000:10000,2]))
lines(density(samp.2$samples[5000:10000,2]))
lines(density(samp.3$samples[5000:10000,2]))

## Test
plot(density(samp.1$samples[5000:10000,3]),xlim = c(0,10))
lines(density(samp.2$samples[5000:10000,3]))
lines(density(samp.3$samples[5000:10000,3]))

## SD and Mean
Sel.rows = c(8001:18000)
mean(samp.1$samples[Sel.rows,3])
var(samp.1$samples[Sel.rows,3])
mean(samp.2$samples[Sel.rows,3])
var(samp.2$samples[Sel.rows,3])
mean(samp.3$samples[Sel.rows,3])
var(samp.3$samples[Sel.rows,3])

# 100 chains for uncertainty ----------------------------------------------

## Initial Values
mean(samp.1$samples[-c(1:8000), 1])
mean(samp.1$samples[-c(1:8000), 2])
mean(samp.1$samples[-c(1:8000), 3])
MCMC.series = data.frame("mean" = NA, "var" = NA, "sd" = NA)
i = 1

## MCMC
for(i in 1:100) {
  samp.1 <- MCMC(EC_Con_Log, n=18000, init= c(0.29,0.43,5.4), adapt=TRUE, acc.rate=0.2)
  MCMC.series[i,"mean"] = mean(samp.1$samples[-c(1:8000), 3])
  MCMC.series[i,"var"] = var(samp.1$samples[-c(1:8000), 3])
  MCMC.series[i,"sd"] = sd(samp.1$samples[-c(1:8000), 3])
  print(i)
}

## SD
apply(MCMC.series, 2, sd)

## Write
write.csv(MCMC.series,"Output Data/CMIP6_EC/EC.summary/MCMC.series_100_Chains.csv")


# 1. MCMC for Cor = 0.90 ------------------------------------------------------------

## All MCMC
MCMC_Exp = data.frame("y_1" = NA, "CEC" = NA)

## Parameters
RVM_Fitted
RVM_Fitted$par
RVM_Fitted$Matrix
RVM_Fitted$family

## Replace Parameters
ts_1 = BiCopTau2Par(family = 5, tau = 0.90, check.taus = TRUE)

## Copula
RVM_Fitted_Ind = D2RVine(order = c(1,3,2), family = c(1,5,5), par = c(0.58,4.19,35))
RVM_Fitted
RVM_Fitted_Ind

## EC for TAS.var.His
beta = c(0.3,0.1,4)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  x_1_pre = beta[1]
  x_2_pre = beta[2]
  y_1_pre = beta[3]
  x_1_obs = Obs[1]
  x_2_obs = Obs[2]
  
  ## Fitted Density
  x_1_d = Marignal.distribution(x_1_pre,Type = "d", Num = "x_1")
  x_2_d = Marignal.distribution(x_2_pre,Type = "d", Num = "x_2")
  y_1_d = Marignal.distribution(y_1_pre,Type = "d", Num = "y_1")
  
  ## Fitted Probability
  x_1_p = Marignal.distribution(x_1_pre,Type = "p", Num = "x_1")
  x_2_p = Marignal.distribution(x_2_pre,Type = "p", Num = "x_2")
  y_1_p = Marignal.distribution(y_1_pre,Type = "p", Num = "y_1")
  
  ## Copula Function
  u <- c(x_1_p,x_2_p,y_1_p)
  Copula_PDF = RVinePDF(u, RVM_Fitted_Ind, verbose = FALSE)
  prior = Copula_PDF * x_1_d * x_2_d * y_1_d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(x_1_obs,x_2_obs),
                    log = FALSE,
                    mean = c(x_1_pre,x_2_pre),
                    sigma = All.covariance)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs)

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])
MCMC_Vine_1 = data.frame("y_1" = Postrior_Re,
                         "CEC" = "Exp A")
MCMC_Exp[1:30000 + 0 * 30000,] = MCMC_Vine_1

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Summary
mean(samp.1$samples[-c(1:8000), 3])
var(samp.1$samples[-c(1:8000), 3])
sd(samp.1$samples[-c(1:8000), 3])
cor(samp.1$samples[-c(1:8000), 1],samp.1$samples[-c(1:8000), 2],method = "kendall")
cor(samp.1$samples[-c(1:8000), 1],samp.1$samples[-c(1:8000), 3],method = "kendall")
cor(samp.1$samples[-c(1:8000), 2],samp.1$samples[-c(1:8000), 3],method = "kendall")

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)

## SD and Mean
Sel.rows = c(8001:18000)
plot(density((samp.1$samples[Sel.rows,3])))
mean(samp.1$samples[Sel.rows,3])
var(samp.1$samples[Sel.rows,3])
lines(density((samp.2$samples[Sel.rows,3])))
mean(samp.2$samples[Sel.rows,3])
var(samp.2$samples[Sel.rows,3])
lines(density((samp.3$samples[Sel.rows,3])))
mean(samp.3$samples[Sel.rows,3])
var(samp.3$samples[Sel.rows,3])


# 2. MCMC for EC Strength = 0.90 ------------------------------------------------------------

## Parameters
ts_1 = BiCopTau2Par(family = 1, tau = 0.90, check.taus = TRUE)
ts_2 = BiCopTau2Par(family = 5, tau = 0.90, check.taus = TRUE)

## Copula
RVM_Fitted_Ind = D2RVine(order = c(1,3,2), family = c(1,5,5), par = c(0.98,35,1.42))
RVM_Fitted
RVM_Fitted_Ind

## EC for TAS.var.His
beta = c(0.3,0.1,4)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  x_1_pre = beta[1]
  x_2_pre = beta[2]
  y_1_pre = beta[3]
  x_1_obs = Obs[1]
  x_2_obs = Obs[2]
  
  ## Fitted Density
  x_1_d = Marignal.distribution(x_1_pre,Type = "d", Num = "x_1")
  x_2_d = Marignal.distribution(x_2_pre,Type = "d", Num = "x_2")
  y_1_d = Marignal.distribution(y_1_pre,Type = "d", Num = "y_1")
  
  ## Fitted Probability
  x_1_p = Marignal.distribution(x_1_pre,Type = "p", Num = "x_1")
  x_2_p = Marignal.distribution(x_2_pre,Type = "p", Num = "x_2")
  y_1_p = Marignal.distribution(y_1_pre,Type = "p", Num = "y_1")
  
  ## Copula Function
  u <- c(x_1_p,x_2_p,y_1_p)
  Copula_PDF = RVinePDF(u, RVM_Fitted_Ind, verbose = FALSE)
  prior = Copula_PDF * x_1_d * x_2_d * y_1_d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(x_1_obs,x_2_obs),
                    log = FALSE,
                    mean = c(x_1_pre,x_2_pre),
                    sigma = All.covariance)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs)

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])
MCMC_Vine_1 = data.frame("y_1" = Postrior_Re,
                         "CEC" = "Exp B")
MCMC_Exp[1:30000 + 1 * 30000,] = MCMC_Vine_1

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Summary
mean(samp.1$samples[-c(1:8000), 3])
var(samp.1$samples[-c(1:8000), 3])
sd(samp.1$samples[-c(1:8000), 3])
cor(samp.1$samples[-c(1:8000), 1],samp.1$samples[-c(1:8000), 2])

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)

# 3. MCMC for IVA ------------------------------------------------------------

## EC for TAS.var.His
beta = c(0.3,0.1,4)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  x_1_pre = beta[1]
  x_2_pre = beta[2]
  y_1_pre = beta[3]
  x_1_obs = Obs[1]
  x_2_obs = Obs[2]
  
  ## Fitted Density
  x_1_d = Marignal.distribution(x_1_pre,Type = "d", Num = "x_1")
  x_2_d = Marignal.distribution(x_2_pre,Type = "d", Num = "x_2")
  y_1_d = Marignal.distribution(y_1_pre,Type = "d", Num = "y_1")
  
  ## Fitted Probability
  x_1_p = Marignal.distribution(x_1_pre,Type = "p", Num = "x_1")
  x_2_p = Marignal.distribution(x_2_pre,Type = "p", Num = "x_2")
  y_1_p = Marignal.distribution(y_1_pre,Type = "p", Num = "y_1")
  
  ## Copula Function
  u <- c(x_1_p,x_2_p,y_1_p)
  Copula_PDF = RVinePDF(u, RVM_Fitted, verbose = FALSE)
  prior = Copula_PDF * x_1_d * x_2_d * y_1_d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(x_1_obs,x_2_obs),
                    log = FALSE,
                    mean = c(x_1_pre,x_2_pre),
                    sigma = All.covariance / 5)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs)

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])
MCMC_Vine_1 = data.frame("y_1" = Postrior_Re,
                         "CEC" = "Exp C")
MCMC_Exp[1:30000 + 2 * 30000,] = MCMC_Vine_1

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)

## Summary
mean(samp.1$samples[-c(1:8000), 3])
var(samp.1$samples[-c(1:8000), 3])
sd(samp.1$samples[-c(1:8000), 3])
cor(samp.1$samples[-c(1:8000), 1],samp.1$samples[-c(1:8000), 2])


# 4. MCMC for Different EC ------------------------------------------------------------

## EC for TAS.var.His
beta = c(0.3,0.1,4)
Obs = Obs_IVAR_Summary$Obs

## Change
Obs = c(0.2,0.3)

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  x_1_pre = beta[1]
  x_2_pre = beta[2]
  y_1_pre = beta[3]
  x_1_obs = Obs[1]
  x_2_obs = Obs[2]
  
  ## Fitted Density
  x_1_d = Marignal.distribution(x_1_pre,Type = "d", Num = "x_1")
  x_2_d = Marignal.distribution(x_2_pre,Type = "d", Num = "x_2")
  y_1_d = Marignal.distribution(y_1_pre,Type = "d", Num = "y_1")
  
  ## Fitted Probability
  x_1_p = Marignal.distribution(x_1_pre,Type = "p", Num = "x_1")
  x_2_p = Marignal.distribution(x_2_pre,Type = "p", Num = "x_2")
  y_1_p = Marignal.distribution(y_1_pre,Type = "p", Num = "y_1")
  
  ## Copula Function
  u <- c(x_1_p,x_2_p,y_1_p)
  Copula_PDF = RVinePDF(u, RVM_Fitted, verbose = FALSE)
  prior = Copula_PDF * x_1_d * x_2_d * y_1_d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(x_1_obs,x_2_obs),
                    log = FALSE,
                    mean = c(x_1_pre,x_2_pre),
                    sigma = All.covariance)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs)

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])
MCMC_Vine_1 = data.frame("y_1" = Postrior_Re,
                         "CEC" = "Exp D")
MCMC_Exp[1:30000 + 3 * 30000,] = MCMC_Vine_1

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)

## Plot hist
hist(samp.1$samples[-c(1:8000), 3])
hist(samp.2$samples[-c(1:8000), 3])
hist(samp.3$samples[-c(1:8000), 3])
mean(samp.1$samples[-c(1:8000), 3])
mean(samp.2$samples[-c(1:8000), 3])
mean(samp.3$samples[-c(1:8000), 3])
var(samp.1$samples[-c(1:8000), 3])
var(samp.2$samples[-c(1:8000), 3])
var(samp.3$samples[-c(1:8000), 3])

## Obs
var(Final_EC$y_1)

# Summary -----------------------------------------------------------------

## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp A")
mean(MCMC_1$y_1)
var(MCMC_1$y_1)

## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp B")
mean(MCMC_1$y_1)
var(MCMC_1$y_1)

## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp C")
mean(MCMC_1$y_1)
var(MCMC_1$y_1)

## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp D")
mean(MCMC_1$y_1)
var(MCMC_1$y_1)

## Write
write.csv(MCMC_Exp,"Output Data/CMIP6/tas/MCMC_Exp.csv")


# Test for difference ECs --------------------------------------------------------------------

## Copula
RVM_Fitted_Ind = RVM_Fitted
RVM_Fitted
RVM_Fitted$Matrix
RVM_Fitted$par

## Select Tau
Tau.copula = RVM_Fitted$par
Tau.copula
Paras.copula = RVM_Fitted$par
Paras.copula

## Replace
ts = BiCopTau2Par(family = 5, tau = 0.95, check.taus = TRUE)
Tau.copula[2,1] = 0.95
Paras.copula[2,1] = ts
RVM_Fitted_Ind$tau = Tau.copula
RVM_Fitted_Ind$par = Paras.copula

## EC for TAS.var.His
beta = c(0.1,0.04,0.6)
Obs = Obs_IVAR_Summary$Obs

## EC for TAS.var.His
beta = c(0.1,0.04,0.6)
Obs = Obs_IVAR_Summary$Obs

## Change
Obs = c(0.05,0.1)

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  His_1 = beta[1]
  His_2 = beta[2]
  Future = beta[3]
  Obs_1 = Obs[1]
  Obs_2 = Obs[2]
  
  ## Fitted Density
  His_1.d = Marignal.distribution(His_1,Type = "d", Num = 1)
  His_2.d = Marignal.distribution(His_2,Type = "d", Num = 2)
  Future.d = Marignal.distribution(Future,Type = "d", Num = 3)
  
  ## Fitted Probability
  His_1.pdf = Marignal.distribution(His_1,Type = "p", Num = 1)
  His_2.pdf = Marignal.distribution(His_2,Type = "p", Num = 2)
  Future.pdf = Marignal.distribution(Future,Type = "p", Num = 3)
  
  ## Copula Function
  u <- c(His_1.pdf, His_2.pdf, Future.pdf)
  Copula_PDF = RVinePDF(u, RVM_Fitted, verbose = FALSE)
  prior = Copula_PDF*His_1.d* His_2.d* Future.d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(Obs_1,Obs_2),
                    log = FALSE,
                    mean = c(His_1,His_2),
                    sigma = All.covariance)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs)

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)

## Test
plot(density(samp.1$samples[5000:10000,1]))
lines(density(samp.2$samples[5000:10000,1]))
lines(density(samp.3$samples[5000:10000,1]))

## Test
plot(density(samp.1$samples[5000:10000,2]))
lines(density(samp.2$samples[5000:10000,2]))
lines(density(samp.3$samples[5000:10000,2]))

## Test
plot(density(samp.1$samples[5000:10000,3]))
lines(density(samp.2$samples[5000:10000,3]))
lines(density(samp.3$samples[5000:10000,3]))

## SD and Mean
mean(samp.1$samples[10000:13000,3])
var(samp.1$samples[10000:13000,3])

mean(samp.2$samples[10000:13000,3])
var(samp.2$samples[10000:13000,3])

mean(samp.3$samples[10000:13000,3])
var(samp.3$samples[10000:13000,3])



# Test for greater Magnitudes ------------------------------------------------------------

## EC for TAS.var.His
beta = c(0.1,0.04,0.5) * 100
Obs = Obs_IVAR_Summary$Obs * 100

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  x_1_pre = beta[1]
  x_2_pre = beta[2]
  y_1_pre = beta[3]
  x_1_obs = Obs[1]
  x_2_obs = Obs[2]
  
  ## Fitted Density
  x_1_d = Marignal.distribution(x_1_pre/100,Type = "d", Num = 1)/100
  x_2_d = Marignal.distribution(x_2_pre/100,Type = "d", Num = 2)/100
  y_1_d = Marignal.distribution(y_1_pre/100,Type = "d", Num = 3)/100
  
  ## Fitted Probability
  x_1_p = Marignal.distribution(x_1_pre/100,Type = "p", Num = 1)
  x_2_p = Marignal.distribution(x_2_pre/100,Type = "p", Num = 2)
  y_1_p = Marignal.distribution(y_1_pre/100,Type = "p", Num = 3)
  
  ## Copula Function
  u <- c(x_1_p,x_2_p,y_1_p)
  Copula_PDF = RVinePDF(u, RVM_Fitted, verbose = FALSE)
  prior = Copula_PDF * x_1_d * x_2_d * y_1_d
  
  ## Conditional Distribution
  Con_cdf = dmvnorm(x = c(x_1_obs,x_2_obs),
                    log = FALSE,
                    mean = c(x_1_pre,x_2_pre),
                    sigma = All.covariance * 100 * 100)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs)

## Saving
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 3], 
                samp.2$samples[-c(1:8000), 3], 
                samp.3$samples[-c(1:8000), 3])
MCMC_Vine_1 = data.frame("y_1" = Postrior_Re,
                         "CEC" = "CEC")

## Summary
mean(MCMC_Vine_1$y_1)
var(MCMC_Vine_1$y_1)
sd(MCMC_Vine_1$y_1)

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)

## Test
plot(density(samp.1$samples[5000:10000,1]))
lines(density(samp.2$samples[5000:10000,1]))
lines(density(samp.3$samples[5000:10000,1]))

## Test
plot(density(samp.1$samples[5000:10000,2]))
lines(density(samp.2$samples[5000:10000,2]))
lines(density(samp.3$samples[5000:10000,2]))

## Test
plot(density(samp.1$samples[5000:10000,3]))
lines(density(samp.2$samples[5000:10000,3]))
lines(density(samp.3$samples[5000:10000,3]))

## SD and Mean
Sel.rows = c(8001:18000)
mean(samp.1$samples[Sel.rows,3])
sd(samp.1$samples[Sel.rows,3])
mean(samp.2$samples[Sel.rows,3])
sd(samp.2$samples[Sel.rows,3])
mean(samp.3$samples[Sel.rows,3])
sd(samp.3$samples[Sel.rows,3])


# End ---------------------------------------------------------------------
