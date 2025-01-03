

# Obs Data -------------------------------------------------------------

## Data
Final_EC = read.csv("Output Data/CMIP6/tas/Final.series_EC.csv",row.names = 1)
Final_Copula = read.csv("Output Data/CMIP6/tas/Final_Copula.csv",row.names = 1)
Obs_IVAR_Summary = read.csv("Output Data/CMIP6_EC/PMIP.files_IVAR/Obs_IVAR_Summary.csv",row.names = 1)
Univariate.best = read.csv("Output Data/CMIP6/tas/Univariate.best.csv",row.names = 1)
set.seed(1)

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

# Hierarchical EC for x_1 ----------------------------------------------------

## EC for x_1
EC_1 = EC.fun(His = Final_EC$x_1, 
              Future  = Final_EC$y_1,
              Obs.mean = Obs_IVAR_Summary$Obs[1],
              Int.Var = All.covariance[1,1],
              Type = "x_1")

## EC for x_2
EC_2 = EC.fun(His = Final_EC$x_2, 
              Future  = Final_EC$y_1,
              Obs.mean = Obs_IVAR_Summary$Obs[2],
              Int.Var = All.covariance[2,2],
              Type = "x_2")

## EC
HEC_Summary = rbind(EC_1,EC_2)
HEC_Summary

## Write
write.csv(HEC_Summary,"Output Data/CMIP6_EC/Summary/HEC_Summary.csv")


# MCMC for HEC for x_1 ------------------------------------------------------------

## MCMC
HEC_MCMC = data.frame("y_1" = NA, "CEC" = NA)

## EC for TAS.var.His
beta = c(0.5,5)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  His = beta[1]
  Future = beta[2]
  
  ## Prior
  temp = dmvnorm(
    x = c(His, Future),
    mean = c(mean(Final_EC$x_1),  mean(Final_EC$y_1)),
    sigma = cov(Final_EC[, c("x_1", "y_1")])
  )
  prior = temp
  
  ## Conditional Distribution
  Con_cdf = dnorm(Obs, mean = His, sd = All.covariance[1,1]^0.5)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs[1])

## MCMC
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 2], 
                samp.2$samples[-c(1:8000), 2], 
                samp.3$samples[-c(1:8000), 2])
CEC_Two_1 = data.frame("y_1" = Postrior_Re, 
                       "CEC" = "x_1")
HEC_MCMC[1:30000 + 0 * 30000,] = CEC_Two_1

## Mean and Var
mean(CEC_Two_1$y_1)
var(CEC_Two_1$y_1)

# MCMC for HEC for x_2 ------------------------------------------------------------

## EC for TAS.var.His
beta = c(0.15,5)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log_pattern = function(Obs) function(beta) {
  
  ## Marginal Distribution
  His = beta[1]
  Future = beta[2]
  
  ## Prior
  temp = dmvnorm(
    x = c(His, Future),
    mean = c(mean(Final_EC$x_2),  mean(Final_EC$y_1)),
    sigma = cov(Final_EC[, c("x_2", "y_1")])
  )
  prior = temp
  
  ## Conditional Distribution
  Con_cdf = dnorm(Obs, mean = His, sd = All.covariance[2,2]^0.5)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

## Function
EC_Con_Log = EC_Con_Log_pattern(Obs = Obs[2])

## MCMC
samp.1 <- MCMC(EC_Con_Log, n=18000, init= beta * 1.2, adapt=TRUE, acc.rate=0.2)
samp.2 <- MCMC(EC_Con_Log, n=18000, init= beta * 1, adapt=TRUE, acc.rate=0.2)
samp.3 <- MCMC(EC_Con_Log, n=18000, init= beta * 0.8, adapt=TRUE, acc.rate=0.2)

## SD and Mean
Postrior_Re = c(samp.1$samples[-c(1:8000), 2], 
                samp.2$samples[-c(1:8000), 2], 
                samp.3$samples[-c(1:8000), 2])
CEC_Two_1 = data.frame("y_1" = Postrior_Re, 
                       "CEC" = "x_2")
HEC_MCMC[1:30000 + 1 * 30000,] = CEC_Two_1

## Write
write.csv(HEC_MCMC,"Output Data/CMIP6_EC/EC.summary/HEC_MCMC.csv")

# Plot HEC ------------------------------------------------------------

## Read
HEC_MCMC = read.csv("Output Data/CMIP6_EC/EC.summary/HEC_MCMC.csv",row.names = 1)
Sel_Colors = c("red",brewer.pal(8,"Dark2")[c(1,2)])

## HEC_Future
Sim_Density = data.frame("X_Value" = seq(0,10,length.out = 1000))
Sim_Density$D1.HEC = dnorm(Sim_Density$X_Value,mean = HEC_Summary$Future_mean_EC[1],sd = HEC_Summary$Future_Var_EC[1]^0.5)

## Mean and SD
CEC_Two_Sel = subset(HEC_MCMC, CEC == "x_1")
Vel_mean = c(HEC_Summary$Future_mean_EC[1],mean(CEC_Two_Sel$y_1))
Vel_Var = c(HEC_Summary$Future_Var_EC[1], var(CEC_Two_Sel$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p1 = ggplot() +
  geom_histogram(
    data = CEC_Two_Sel,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.5,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[2],
    linewidth = 1,
    aes(x = X_Value, y = D1.HEC)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 1)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[2],linewidth = 1)+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[2])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[2])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,0.6,0.1)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.6),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )


## HEC_Future
Sim_Density = data.frame("X_Value" = seq(0,10,length.out = 1000))
Sim_Density$D1.HEC = dnorm(Sim_Density$X_Value,mean = HEC_Summary$Future_mean_EC[2],sd = HEC_Summary$Future_Var_EC[2]^0.5)

## Mean and SD
CEC_Two_Sel = subset(HEC_MCMC, CEC == "x_2")
Vel_mean = c(HEC_Summary$Future_mean_EC[2],mean(CEC_Two_Sel$y_1))
Vel_Var = c(HEC_Summary$Future_Var_EC[2], var(CEC_Two_Sel$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p2 = ggplot() +
  geom_histogram(
    data = CEC_Two_Sel,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.5,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[2],
    linewidth = 1,
    aes(x = X_Value, y = D1.HEC)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 1)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[2],linewidth = 1)+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[2])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[2])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,0.6,0.1)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.6),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )

## Legend
p.legend.1 = ggplot() +
  geom_path(data = data.frame(
    "x" = 1,
    "y" = 1,
    type = c("Analytical equations in HEC", "MCMC simulations for HEC")
  ), 
  aes(x = x, y = y , color = type)) +
  scale_color_manual(
    breaks = c("Analytical equations in HEC", "MCMC simulations for HEC"),
    values = Sel_Colors[c(2, 3)]
  ) +
  theme_figure1+
  theme(
    legend.position = "top",
    legend.justification.top = c(0.5,0.5),
    legend.direction = "vertical"
  )
legend1 = cowplot::get_plot_component(p.legend.1, 'guide-box-top', return_all = TRUE)

## Combine
p2 = p2+annotation_custom(legend1, xmin = 6, ymin = 0.45)

## Total Plot
p.total = plot_grid(p1,p2,
                    labels = "auto",nrow = 1)

## Saving
ggsave(filename = "Figures/Final Figures/Fig_Temp_6_HEC.tiff",
       plot = p.total,
       width = 9, 
       height = 4.5)
ggsave(filename = "Figures/Final Figures/Fig_Temp_6_HEC.pdf",
       plot = p.total,
       device = cairo_pdf,
       width = 9, 
       height = 4)

# End ---------------------------------------------------------------------


