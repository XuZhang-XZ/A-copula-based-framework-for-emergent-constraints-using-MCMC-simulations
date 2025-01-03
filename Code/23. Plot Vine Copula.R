


# Obs Data -------------------------------------------------------------

## Data
Final_EC = read.csv("Output Data/CMIP6/tas/Final.series_EC.csv",row.names = 1)
Final_Copula = read.csv("Output Data/CMIP6/tas/Final_Copula.csv",row.names = 1)
Obs_IVAR_Summary = read.csv("Output Data/CMIP6_EC/PMIP.files_IVAR/Obs_IVAR_Summary.csv",row.names = 1)
Univariate.best = read.csv("Output Data/CMIP6/tas/Univariate.best.csv",row.names = 1)

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

## Read
Sel_Colors = brewer.pal(9,"Set1")

## Cor
cor.test(Final_EC$x_1_p,Final_EC$y_1_p, method = "kendall")
cor.test(Final_EC$x_2_p,Final_EC$y_1_p, method = "kendall")
cor.test(Final_EC$x_1_p,Final_EC$x_2_p, method = "kendall")

# Plot MCMC ------------------------------------------------------------

## Data
MCMC_CEC_Vine = read.csv("Output Data/CMIP6/tas/MCMC_CEC_Vine.csv",row.names = 1)

## Raw Density
Sim_Density = data.frame("X_Value" = seq(0,10,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_Value, Type = "d", Num = "y_1")

## SD and Mean
MCMC_1 = subset(MCMC_CEC_Vine,CEC == "CEC")

## Mean and SD
Vel_mean = c(mean(Final_EC$y_1),mean(MCMC_1$y_1))
Vel_Var = c(var(Final_EC$y_1),var(MCMC_1$y_1))
(Vel_Var[1] - Vel_Var[2]) / (Vel_Var[1])

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p1 = ggplot() +
  geom_histogram(
    data = MCMC_1,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.4,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[1],
    linewidth = 0.8,
    aes(x = X_Value, y = D1)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
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

## Legend
p.legend.1 = ggplot() +
  geom_path(data = data.frame("x" = 1,"y" = 1, type = c("Unconstrained","CEC-constrained")),linewidth = 0.8,
            aes(x = x, y = y ,color = type)) +
  scale_color_manual(breaks = c("Unconstrained","CEC-constrained"), values = Sel_Colors[c(1,3)])+
  theme_figure1+
  theme(
    legend.position = "top",
    legend.justification.top = c(0.9,0.88),
    legend.direction = "vertical"
  )
legend1 = cowplot::get_plot_component(p.legend.1, 'guide-box-top', return_all = TRUE)

## Combine
p2 = p1 +
  annotation_custom(grob = legend1,
                    xmin = 8.0,
                    xmax = 9.0,
                    ymin = 0.5,
                    ymax = 0.6)

## Void
p.left = ggplot() +
  theme_void()

## Total Plot
p.total = plot_grid(p.left, p1, rel_widths = c(2.5, 5.5), labels = "auto")

## Saving
ggsave(filename = "Figures/Final Figures/Fig_Temp_7_Joint_CEC.tiff",
       plot = p2,
       width = 5.5, 
       height = 4.5)
ggsave(filename = "Figures/Final Figures/Fig_Temp_7_Joint_CEC.svg",
       plot = p2,
       width = 5.5, 
       height = 4.5)


# Plot GR metric ---------------------------------------------------------

## Data
Shrink_factor_Vine = read.csv("Output Data/CMIP6_EC/EC.summary/Shrink_factor_Vine.csv",row.names = 1)

## Plot
p1 = ggplot() +
  geom_path(data = Shrink_factor_Vine,
            aes(x = last.iter, y = Factor_3)) +
  geom_path(data = Shrink_factor_Vine,
            color = "red",linetype = "dashed",
            aes(x = last.iter, y = Factor_3_Upper)) +
  scale_x_continuous(name = "Number of iterations",breaks = c(0,4000,8000,12000,18000)) +
  scale_y_continuous(name = "Gelman and Rubin's diagnostics") +
  theme_figure1+
  coord_cartesian(xlim = c(0,18000), ylim = c(0,5),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8)
  )

## Saving
ggsave(filename = "Figures/Supplementary Figures/G_R_Vine.tiff",
       plot = p1,
       width = 7, 
       height = 5)


# Plot for 4 experiments ------------------------------------------------------------

## Data
MCMC_CEC_Vine = read.csv("Output Data/CMIP6/tas/MCMC_Exp.csv",row.names = 1)

## Raw Density
Sim_Density = data.frame("X_Value" = seq(0,10,length.out = 1000))
Sim_Density$D1 = Marignal.distribution(x = Sim_Density$X_Value, Type = "d", Num = "y_1")

## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp A")

## Mean and SD
Vel_mean = c(mean(Final_EC$y_1),mean(MCMC_1$y_1))
Vel_Var = c(var(Final_EC$y_1),var(MCMC_1$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p1 = ggplot() +
  geom_histogram(
    data = MCMC_1,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.4,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[1],
    linewidth = 0.8,
    aes(x = X_Value, y = D1)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.2, label = "Exp A", color = "red")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[1])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[1])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,1,0.2)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.8),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )




## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp B")

## Mean and SD
Vel_mean = c(mean(Final_EC$y_1),mean(MCMC_1$y_1))
Vel_Var = c(var(Final_EC$y_1),var(MCMC_1$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p2 = ggplot() +
  geom_histogram(
    data = MCMC_1,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.4,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[1],
    linewidth = 0.8,
    aes(x = X_Value, y = D1)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.2, label = "Exp B", color = "red")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[1])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[1])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,1,0.2)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,1),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )




## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp C")

## Mean and SD
Vel_mean = c(mean(Final_EC$y_1),mean(MCMC_1$y_1))
Vel_Var = c(var(Final_EC$y_1),var(MCMC_1$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p3 = ggplot() +
  geom_histogram(
    data = MCMC_1,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.4,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[1],
    linewidth = 0.8,
    aes(x = X_Value, y = D1)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.2, label = "Exp C", color = "red")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[1])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[1])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,1,0.2)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.8),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )




## SD and Mean
MCMC_1 = subset(MCMC_Exp,CEC == "Exp D")

## Mean and SD
Vel_mean = c(mean(Final_EC$y_1),mean(MCMC_1$y_1))
Vel_Var = c(var(Final_EC$y_1),var(MCMC_1$y_1))

## Print
Vel_mean1 = sprintf("%.3f", round(Vel_mean,3))
Vel_Var1 = sprintf("%.3f", round(Vel_Var,3))

## Plot TAS.var.His
p4 = ggplot() +
  geom_histogram(
    data = MCMC_1,
    breaks = seq(0,10,0.25),
    fill = Sel_Colors[3],
    color = "#e9ecef",
    alpha = 0.4,
    aes(x = y_1, y = after_stat(density))
  ) +
  geom_path(
    data = Sim_Density,
    color = Sel_Colors[1],
    linewidth = 0.8,
    aes(x = X_Value, y = D1)
  ) +
  geom_vline(xintercept = Vel_mean[1], color = Sel_Colors[1],linewidth = 0.8)+
  geom_vline(xintercept = Vel_mean[2], color = Sel_Colors[3],linewidth = 0.8)+
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.2, label = "Exp D", color = "red")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -10.5, label = "Mean", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -9.0, label = Vel_mean1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -7.5, label = Vel_mean1[1], color = Sel_Colors[1])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -4.5, label = "Var ", color = "black")+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -3.0, label = Vel_Var1[2], color = Sel_Colors[3])+
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.5, label = Vel_Var1[1], color = Sel_Colors[1])+
  scale_x_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)")) +
  scale_y_continuous(name = "Probability density",breaks = seq(0,0.6,0.2)) +
  theme_figure1+
  coord_cartesian(xlim = c(1,10), ylim = c(0,0.6),expand = FALSE)+
  theme(
    plot.margin = margin(10, 15, 5, 10, "pt"),
    legend.position = "none"
  )


## Legend
p.legend.1 = ggplot() +
  geom_path(data = data.frame("x" = 1,"y" = 1, type = c("Unconstrained","CEC-constrained")),linewidth = 0.8,
            aes(x = x, y = y ,color = type)) +
  scale_color_manual(breaks = c("Unconstrained","CEC-constrained"), values = Sel_Colors[c(1,3)])+
  theme_figure1+
  theme(
    legend.position = "top",
    legend.justification.top = c(0.9,0.88),
    legend.direction = "vertical"
  )
legend1 = cowplot::get_plot_component(p.legend.1, 'guide-box-top', return_all = TRUE)

## Combine
p4 = p4 +
  annotation_custom(grob = legend1,
                    xmin = 7.5,
                    xmax = 8.0,
                    ymin = 0.5,
                    ymax = 0.6)

## Summary
p.total = plot_grid(p1, p2, p3, p4, nrow = 2, labels = "auto")

## Saving
ggsave(filename = "Figures/Final Figures/Fig_Temp_8.tiff",
       plot = p.total,
       width = 9, 
       height = 8)


# End -------------------------------------------------------------
