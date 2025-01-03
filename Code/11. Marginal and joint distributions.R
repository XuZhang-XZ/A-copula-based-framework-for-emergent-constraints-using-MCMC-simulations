

# Marginal Distributions -------------------------------------------

## Data
Final_EC = read.csv("Output Data/CMIP6_EC/Summary/Final_EC_Sel.csv",row.names = 1)

## Functions
Univariate.summary = expand.grid(
  "Variable" = c("x_1", "x_2", "y_1"),
  "Functions" = c("normal", "logistic", "cauchy","lognormal","weibull","gamma"),
  stringsAsFactors = FALSE
)
i = 1

## Best functions
for(i in 1:NROW(Univariate.summary)) {
  
  ## AIC
  temp1 = Univariate.summary$Variable[i]
  temp2 = Univariate.summary$Functions[i]
  ts = fitdistr(Final_EC[,temp1], temp2)
  Univariate.summary[i, "AIC"] = 2 * 2 - 2 * ts$loglik
  
  ## Paras
  Paras = ts$estimate
  Univariate.summary[i,"para1"] = Paras[1]
  Univariate.summary[i,"para2"] = Paras[2]
}

## The Best One
Univariate.best = Univariate.summary %>%
  dplyr::group_by(Variable) %>%
  dplyr::summarise(across(everything(), ~.x[which.min(AIC)]))

## Write
write.csv(Univariate.summary,"Output Data/CMIP6/tas/Univariate.summary.csv")
write.csv(Univariate.best,"Output Data/CMIP6/tas/Univariate.best.csv")

## Write
Univariate.summary$AIC_1 = sprintf("%.2f", round(Univariate.summary$AIC,2))
Univ_Longer = pivot_wider(data = Univariate.summary,id_cols = Functions,names_from = Variable,values_from = AIC_1)
write_xlsx(Univ_Longer,"Output Data/CMIP6/tas/Univ_Longer.xlsx")

## Density
for(i in 1:NROW(Univariate.best)) {
  temp1 = Univariate.best$Variable[i]
  temp2 = Univariate.best$Functions[i]
  Final_EC[, paste0(temp1, "_d")] = Marignal.distribution(x = Final_EC[, temp1], Type = "d", Num = temp1)
  Final_EC[, paste0(temp1, "_p")] = Marignal.distribution(x = Final_EC[, temp1], Type = "p", Num = temp1)
}

## Cor
cor.test(Final_EC$x_1_p,Final_EC$y_1_p, method = "kendall")
cor.test(Final_EC$x_2_p,Final_EC$y_1_p, method = "kendall")
cor.test(Final_EC$x_1_p,Final_EC$x_2_p, method = "kendall")


# 2-d Joint Distributions ---------------------------------------------

## Copula
Copula.paras = data.frame(
  Copula.family = c(1, 2, 4, 5, 6), 
  Copula = c("Gaussian","Student t", "Gumbel", "Frank", "Joe")
)
i = 1

## Results
for(i in 1:NROW(Copula.paras)) {
  
  ## x_1_p & y_1_p
  u <- as.matrix(Final_EC[,c("x_1_p","y_1_p")])
  fit.copula = BiCopEst(u1 = u[,1], u2 = u[,2], family = Copula.paras$Copula.family[i], method = "mle")
  Copula.paras[i,"para_1"] = fit.copula$par
  Copula.paras[i,"AIC_1"] = fit.copula$AIC
  
  ## x_2_p & y_1_p
  u <- as.matrix(Final_EC[,c("x_2_p","y_1_p")])
  fit.copula = BiCopEst(u1 = u[,1], u2 = u[,2], family = Copula.paras$Copula.family[i], method = "mle")
  Copula.paras[i,"para_2"] = fit.copula$par
  Copula.paras[i,"AIC_2"] = fit.copula$AIC
}

## Change  Akaike information criterion AIC
AIC_1_Shift = sum(log(c(Final_EC$x_1_d,Final_EC$y_1_d)))
AIC_2_Shift = sum(log(c(Final_EC$x_2_d,Final_EC$y_1_d)))
Copula.paras$AIC_1 = Copula.paras$AIC_1 - 2*AIC_1_Shift + 2*4
Copula.paras$AIC_2 = Copula.paras$AIC_2 - 2*AIC_2_Shift + 2*4

## MVNORM
i = NROW(Copula.paras) + 1
Copula.paras[i,"Copula"] = "mvnorm"

## x_1 ~ y_1
temp = dmvnorm(x = Final_EC[,c("x_1","y_1")],log = FALSE,
               mean = c(mean(Final_EC$x_1),mean(Final_EC$y_1)),
               sigma = cov(Final_EC[,c("x_1","y_1")]))
Copula.paras[i,"AIC_1"] = 2*5 - 2*sum(log(temp))

## x_2 ~ y_1
temp = dmvnorm(x = Final_EC[,c("x_2","y_1")],
               mean = c(mean(Final_EC$x_2),mean(Final_EC$y_1)),
               sigma = cov(Final_EC[,c("x_2","y_1")]))
Copula.paras[i,"AIC_2"] = 2*5 - 2*sum(log(temp))

## Empirical Probability
Final_EC[,"Empri_Prob_x1"] = empirical_bivariate(x = Final_EC$x_1, y = Final_EC$y_1)
Final_EC[,"Empri_Prob_x2"] = empirical_bivariate(x = Final_EC$x_2, y = Final_EC$y_1)

# Empirical Copula --------------------------------------------------------

## Final Copula
C_Sel_1 = which.min(Copula.paras$AIC_1)
C_Sel_2 = which.min(Copula.paras$AIC_2)
Final_Copula = data.frame("Family" = Copula.paras$Copula.family[c(C_Sel_1,C_Sel_2)],
                          "Para" = c(Copula.paras$para_1[C_Sel_1],Copula.paras$para_2[C_Sel_2]),
                          "Copula.name" = Copula.paras$Copula[c(C_Sel_1,C_Sel_2)],
                          "AIC" = c(Copula.paras$AIC_1[C_Sel_1],Copula.paras$AIC_2[C_Sel_2]))

## Print
Copula.paras
Final_Copula

## Copula x_1
u <- as.matrix(Final_EC[,c("x_1_p","y_1_p")])
cop <- BiCop(family = Final_Copula$Family[1], par = Final_Copula$Para[1])
Final_EC[,"Copula_Prob_x1"] = BiCopCDF(u1 = u[,1], u2 = u[,2], cop)

## Copula x_2
u <- as.matrix(Final_EC[,c("x_2_p","y_1_p")])
cop <- BiCop(family = Final_Copula$Family[2], par = Final_Copula$Para[2])
Final_EC[,"Copula_Prob_x2"] = BiCopCDF(u1 = u[,1], u2 = u[,2], cop)

# Write -------------------------------------------------------------------

## Write
write.csv(Final_EC,"Output Data/CMIP6/tas/Final.series_EC.csv")

## Write
write.csv(Final_Copula,"Output Data/CMIP6/tas/Final_Copula.csv")

## Write
Copula.paras$AIC_1 = sprintf("%.2f", round(Copula.paras$AIC_1,2))
Copula.paras$AIC_2 = sprintf("%.2f", round(Copula.paras$AIC_2,2))
write.csv(Copula.paras,"Output Data/CMIP6/tas/Copula.paras.csv")
write_xlsx(Copula.paras,"Output Data/CMIP6/tas/Copula.paras.xlsx")

# Plot marginal distributions -------------------------------------

## Colors
Sel_Colors = brewer.pal(9,"Set1")

## Data
Final_EC = read.csv("Output Data/CMIP6/tas/Final.series_EC.csv",row.names = 1)

## x_1
Sim_Density = data.frame("X_value" = seq(0,1,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_value, Type = "d", Num = "x_1")

## Plot
p1 = ggplot() +
  geom_histogram(
    data = Final_EC,
    binwidth = 0.07,
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.6,
    aes(x = x_1, y = after_stat(density))
  ) +
  geom_path(data = Sim_Density,color = Sel_Colors[2],
            aes(x = X_value, y = D1))+
  scale_x_continuous(name =  expression("1984-2014 " * italic(T["Trend,SH"]) * " (°C per decade)")) +
  scale_y_continuous(name = "Probability density") +
  theme_figure1+
  coord_cartesian(xlim = c(0,0.8), ylim = c(0,6),expand = FALSE)+
  theme(
    plot.margin = margin(7, 10, 7, 10, "pt"),
    legend.position = "inside",
    legend.justification.inside = c(0.98,0.92)
  )

## x_2
Sim_Density = data.frame("X_value" = seq(0,1,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_value, Type = "d", Num = "x_2")

## Plot
p2 = ggplot() +
  geom_histogram(
    data = Final_EC,
    binwidth = 0.04,
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.6,
    aes(x = x_2, y = after_stat(density))
  ) +
  geom_path(data = Sim_Density, color = Sel_Colors[2], aes(x = X_value, y = D1)) +
  scale_x_continuous(name =  expression("1984-2014 " * italic(T["Trend,SH"]) * " (°C per decade)")) +
  scale_y_continuous(name = "Probability density") +
  theme_figure1+
  coord_cartesian(xlim = c(0,0.5), ylim = c(0,10),expand = FALSE)+
  theme(
    plot.margin = margin(7, 10, 7, 10, "pt"),
    legend.position = "inside",
    legend.justification.inside = c(0.98,0.92)
  )

## y_1
Sim_Density = data.frame("X_value" = seq(0,10,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_value, Type = "d", Num = "y_1")

## Plot
p3 = ggplot() +
  geom_histogram(
    data = Final_EC,
    binwidth = 0.6,
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.6,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(data = Sim_Density,color = Sel_Colors[2],
            aes(x = X_value, y = D1))+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,0.6,0.1)) +
  theme_figure1+
  coord_cartesian(xlim = c(2,10), ylim = c(0,0.6),expand = FALSE)+
  theme(
    plot.margin = margin(7, 10, 7, 10, "pt"),
    legend.position = "inside",
    legend.justification.inside = c(0.98,0.92)
  )


## Total Plot
p.total = plot_grid(p1, p2, p3, labels = "auto", nrow = 1)

## Saving
ggsave(
  filename = "Figures/Final Figures/Fig_Temp_3_Marginal Distribution.pdf",
  plot = p.total,
  device = cairo_pdf,
  width = 11,
  height = 4
)


# Plot Joint Distribution ------------------------------------------------------

## Colors
Sel_Colors = brewer.pal(9,"Set1")

## Data
Final_EC = read.csv("Output Data/CMIP6/tas/Final.series_EC.csv",row.names = 1)

## Plot
p1 = ggplot() +
  geom_abline(intercept = 0, slope = 1, color = Sel_Colors[2], linewidth = 0.5, linetype = "dashed") +
  geom_point(data = Final_EC,
             aes(x = Copula_Prob_x1, y = Empri_Prob_x1))+
  scale_x_continuous(name = "Theoretical probability") +
  scale_y_continuous(name = "Empirical probability") +
  theme_figure1+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1),expand = FALSE)+
  theme(
    plot.margin = margin(10, 20, 10, 10, "pt"),
    legend.position = "inside",
    legend.justification.inside = c(0.12,0.85)
  )

## Plot
p2 = ggplot() +
  geom_abline(intercept = 0, slope = 1, color = Sel_Colors[2], linewidth = 0.5, linetype = "dashed") +
  geom_point(data = Final_EC,
             aes(x = Copula_Prob_x2, y = Empri_Prob_x2))+
  scale_x_continuous(name = "Theoretical probability") +
  scale_y_continuous(name = "Empirical probability") +
  theme_figure1+
  coord_cartesian(xlim = c(0,1), ylim = c(0,1),expand = FALSE)+
  theme(
    plot.margin = margin(10, 20, 10, 10, "pt"),
    legend.position = "inside",
    legend.justification.inside = c(0.12,0.85)
  )

## Total Plot
p.total = plot_grid(p1, p2, labels = "auto", nrow = 1)

## Saving
ggsave(
  filename = "Figures/Final Figures/Fig_Temp_4_Joint Distributions.pdf",
  plot = p.total,
  device = cairo_pdf,
  width = 8,
  height = 4
)

# End -------------------------------------------------------------




