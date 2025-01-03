

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


# Copula-based EC for x_1 --------------------------------------------------

## MCMC
CEC_Two = data.frame("y_1" = NA, "CEC" = NA)

## EC for TAS.var.His
beta = c(0.4,5)
Obs = Obs_IVAR_Summary$Obs

# ## Test
# beta = c(0.4,5) ## 0.595 
# beta = c(0.4,8) ## Clayton -2.84  Gaussian -3.65

## Functions
EC_Con_Log = function(beta) {
  
  ## Marginal Distribution
  His = beta[1]
  Future = beta[2]
  
  ## Prior CDF
  His.pdf = Marignal.distribution(His,Type = "p", Num = "x_1")
  Future.pdf = Marignal.distribution(Future,Type = "p", Num = "y_1")
  
  ## Prior Density
  His.d = Marignal.distribution(His,Type = "d", Num = "x_1")
  Future.d = Marignal.distribution(Future,Type = "d", Num = "y_1")

  # ## Change
  # if(His.pdf < 0.001) {His.pdf = 0.001}
  # if(His.pdf > 0.999) {His.pdf = 0.999}
  # if(Future.pdf < 0.001) {Future.pdf = 0.001}
  # if(Future.pdf > 0.999) {Future.pdf = 0.999}
  
  ## Copula Function
  cop <- BiCop(family = Final_Copula$Family[1], par = Final_Copula$Para[1])
  u <- c(His.pdf, Future.pdf)
  prior = BiCopPDF(u1 = u[1], u2 = u[2], cop)*His.d*Future.d

  ## Conditional Distribution
  Con_cdf = dnorm(Obs[1], mean = His, sd = All.covariance[1,1]^0.5)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

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
CEC_Two[1:30000 + 0 * 30000,] = CEC_Two_1

## Summary
mean(CEC_Two_1$y_1)
var(CEC_Two_1$y_1)
sd(CEC_Two_1$y_1)

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)
gel_1 = gelman.plot(sample.list)

## Convergence
ts = gel_1$shrink
Shrink_factor_1 = data.frame(
  last.iter = gel_1$last.iter,
  Factor_1 = ts[, 1, 1],
  Factor_1_Upper = ts[, 1, 2],
  Factor_2 = ts[, 2, 1],
  Factor_2_Upper = ts[, 2, 2]
)

## SD and Mean
Sel.rows = c(1:8000)
plot(density(samp.1$samples[Sel.rows,2]))
mean(samp.1$samples[Sel.rows,2])
var(samp.1$samples[Sel.rows,2])
mean(samp.2$samples[Sel.rows,2])
var(samp.2$samples[Sel.rows,2])
mean(samp.3$samples[Sel.rows,2])
var(samp.3$samples[Sel.rows,2])

## Test
Sel.rows = c(1:8000)
plot(density(samp.1$samples[Sel.rows,1]))
plot(density(samp.1$samples[Sel.rows,2]))


# Copula-based EC for x_2 --------------------------------------------------

## Parameters
beta = c(0.15,5)
Obs = Obs_IVAR_Summary$Obs

## Functions
EC_Con_Log = function(beta) {
  
  ## Marginal Distribution
  His = beta[1]
  Future = beta[2]
  
  ## Prior CDF
  His.pdf = Marignal.distribution(His,Type = "p", Num = "x_2")
  Future.pdf = Marignal.distribution(Future,Type = "p", Num = "y_1")
  
  ## Prior Density
  His.d = Marignal.distribution(His,Type = "d", Num = "x_2")
  Future.d = Marignal.distribution(Future,Type = "d", Num = "y_1")

  ## Copula Function
  cop <- BiCop(family = Final_Copula$Family[2], par = Final_Copula$Para[2])
  u <- c(His.pdf, Future.pdf)
  prior = BiCopPDF(u1 = u[1], u2 = u[2], cop)*His.d*Future.d
  
  ## Conditional Distribution
  Con_cdf = dnorm(Obs[2], mean = His, sd = All.covariance[2,2]^0.5)
  
  ## Return
  log_density = log(prior) + log(Con_cdf)
  return(log_density)
}

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
CEC_Two[1:30000 + 1 * 30000,] = CEC_Two_1

## Summary
mean(CEC_Two_1$y_1)
var(CEC_Two_1$y_1)
sd(CEC_Two_1$y_1)

## Convergence
samp.coda1 <- convert.to.coda(samp.1)
samp.coda2 <- convert.to.coda(samp.2)
samp.coda3 <- convert.to.coda(samp.3)
sample.list = list(samp.coda1,samp.coda2,samp.coda3)
gelman.diag(sample.list)
gel_2 = gelman.plot(sample.list)

## Convergence
ts = gel_2$shrink
Shrink_factor_2 = data.frame(
  last.iter = gel_2$last.iter,
  Factor_1 = ts[, 1, 1],
  Factor_1_Upper = ts[, 1, 2],
  Factor_2 = ts[, 2, 1],
  Factor_2_Upper = ts[, 2, 2]
)

## SD and Mean
Sel.rows = c(1:8000)
plot(density(samp.1$samples[Sel.rows,2]))
mean(samp.1$samples[Sel.rows,2])
var(samp.1$samples[Sel.rows,2])
mean(samp.2$samples[Sel.rows,2])
var(samp.2$samples[Sel.rows,2])
mean(samp.3$samples[Sel.rows,2])
var(samp.3$samples[Sel.rows,2])

## Write
write.csv(CEC_Two,"Output Data/CMIP6_EC/EC.summary/CEC_Two.csv")
write.csv(Shrink_factor_1,"Output Data/CMIP6_EC/EC.summary/Shrink_factor_1.csv")
write.csv(Shrink_factor_2,"Output Data/CMIP6_EC/EC.summary/Shrink_factor_2.csv")

# Plot CEC results -----------------------------------------------------

## Read
Sel_Colors = brewer.pal(9,"Set1")
CEC_Two = read.csv("Output Data/CMIP6_EC/EC.summary/CEC_Two.csv")

## Raw Density
Sim_Density = data.frame("X_Value" = seq(0,10,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_Value, Type = "d", Num = "y_1")

## Mean and SD
CEC_Two_Sel = subset(CEC_Two, CEC == "x_1")
Vel_mean = c(mean(Final_EC$y_1),mean(CEC_Two_Sel$y_1))
Vel_Var = c(var(Final_EC$y_1),var(CEC_Two_Sel$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))
(Vel_Var[1] - Vel_Var[2]) / Vel_Var[1]

## Plot TAS.var.His
p1 = ggplot() +
  geom_histogram(data = CEC_Two_Sel,
                 breaks = seq(0,10,0.25),
                 fill=Sel_Colors[3], color="#e9ecef", alpha= 0.5,
                 aes(x = y_1, y = after_stat(density)))+
  geom_path(data = Sim_Density,color = Sel_Colors[1],linewidth = 0.8, ## Unconstrained
            aes(x = X_Value, y = D1))+
  geom_point(
    data = Final_EC,
    shape = 8,
    color = Sel_Colors[1],
    aes(x = y_1, y = 0.02)
  ) +
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[1])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[1])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,0.6,0.1)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.6),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )

## Raw Density
Sim_Density = data.frame("X_Value" = seq(0,10,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_Value, Type = "d", Num = "y_1")

## Mean and SD
CEC_Two_Sel = subset(CEC_Two, CEC == "x_2")
Vel_mean = c(mean(Final_EC$y_1),mean(CEC_Two_Sel$y_1))
Vel_Var = c(var(Final_EC$y_1),var(CEC_Two_Sel$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))
(Vel_Var[1] - Vel_Var[2]) / Vel_Var[1]

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
    color = Sel_Colors[1],
    linewidth = 0.8,
    aes(x = X_Value, y = D1)
  ) +
  geom_point(
    data = Final_EC,
    shape = 8,
    color = Sel_Colors[1],
    aes(x = y_1, y = 0.02)
  ) +
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[1])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[1])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,0.6,0.1)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.6),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8)
  )

## Legend
p.legend.1 = ggplot() +
  geom_path(data = data.frame("x" = 1,"y" = 1, type = c("Unconstrained","CEC-constrained")),linewidth = 0.8,
            aes(x = x, y = y ,color = type)) +
  scale_color_manual(breaks = c("Unconstrained", "CEC-constrained"),
                     values = Sel_Colors[c(1, 3)]) +
  theme_figure1+
  theme(
    legend.position = "top",
    legend.justification.top = c(0.9,0.88),
    legend.direction = "vertical"
  )
legend1 = cowplot::get_plot_component(p.legend.1, 'guide-box-top', return_all = TRUE)

## Combine
p2 = p2+annotation_custom(grob = legend1, xmin = 7.5, xmax = 8, ymin = 0.45,ymax = 0.55)

## Total Plot
p.total = plot_grid(p1,p2,
                    labels = "auto",nrow = 1)

## Saving
ggsave(filename = "Figures/Final Figures/Fig_Temp_5_CEC.tiff",
       plot = p.total,
       width = 8, 
       height = 4)
ggsave(filename = "Figures/Final Figures/Fig_Temp_5_CEC.pdf",
       plot = p.total,
       device = cairo_pdf,
       width = 10, 
       height = 4.5)



# Plot Convergence --------------------------------------------------------

## Data
Shrink_factor_1 = read.csv("Output Data/CMIP6_EC/EC.summary/Shrink_factor_1.csv",row.names = 1)

## Plot
p1 = ggplot() +
  geom_path(data = Shrink_factor_1,
            aes(x = last.iter, y = Factor_2)) +
  geom_path(data = Shrink_factor_1,
            color = "red",linetype = "dashed",
            aes(x = last.iter, y = Factor_2_Upper)) +
  scale_x_continuous(name = "Number of iterations",breaks = c(0,4000,8000,12000,18000)) +
  scale_y_continuous(name = "Gelman and Rubin's diagnostics") +
  theme_figure1+
  coord_cartesian(xlim = c(0,18000), ylim = c(0,5),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8)
  )

## Data
Shrink_factor_2 = read.csv("Output Data/CMIP6_EC/EC.summary/Shrink_factor_1.csv",row.names = 1)

## Plot
p2 = ggplot() +
  geom_path(data = Shrink_factor_2,
            aes(x = last.iter, y = Factor_2)) +
  geom_path(data = Shrink_factor_2,
            color = "red",linetype = "dashed",
            aes(x = last.iter, y = Factor_2_Upper)) +
  scale_x_continuous(name = "Number of iterations",breaks = c(0,4000,8000,12000,18000)) +
  scale_y_continuous(name = "Gelman and Rubin's diagnostics") +
  theme_figure1+
  coord_cartesian(xlim = c(0,18000), ylim = c(0,5),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8)
  )

## Summary
p.total = plot_grid(p1, p2, labels = "auto", nrow = 1)

## Saving
ggsave(filename = "Figures/Supplementary Figures/G_R_two.tiff",
       plot = p.total,
       width = 8, 
       height = 4)

# End -------------------------------------------------------------
 
