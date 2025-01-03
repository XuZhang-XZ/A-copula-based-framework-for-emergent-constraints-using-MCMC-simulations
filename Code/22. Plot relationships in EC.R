

# Plot Relationships in EC ----------------------------------------------

## Colors
Sel_Colors = brewer.pal(9,"Set1")

## Data
Final_EC = read.csv("Output Data/CMIP6_EC/Summary/Final_EC_Sel.csv",row.names = 1)
Obs_IVAR_Summary = read.csv("Output Data/CMIP6_EC/PMIP.files_IVAR/Obs_IVAR_Summary.csv",row.names = 1)

## Cor
cor.test(Final_EC$x_1,Final_EC$y_1)
cor.test(Final_EC$x_2,Final_EC$y_1)

## Partial Correlation
ppcor::pcor.test(Final_EC$x_1,Final_EC$y_1,Final_EC$x_2)
ppcor::pcor.test(Final_EC$x_2,Final_EC$y_1,Final_EC$x_1)

## Plot
p1 = ggplot() +
  # geom_rect(
  #   fill = Sel_Colors[2],
  #   alpha = 0.2,
  #   aes(
  #     xmin = Obs_IVAR_Summary$L1[1],
  #     xmax = Obs_IVAR_Summary$U1[1],
  #     ymin = 0,
  #     ymax = 50
  #   )
  # ) +
  geom_vline(xintercept = Obs_IVAR_Summary$Obs[1], color = Sel_Colors[2], linewidth = 0.6) +
  geom_point(data = Final_EC, color = Sel_Colors[3], aes(x = x_1, y = y_1)) +
  geom_smooth(
    data = Final_EC,
    color = Sel_Colors[3],
    method = "lm",
    se = FALSE,
    linewidth = 0.6,
    linetype = "dashed",
    fullrange = TRUE,
    aes(x = x_1, y = y_1)
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -0.5, label = expression(italic(r)* "= 0.557"))+
  scale_x_continuous(name = expression("1984-2014 "*italic(T["Trend,NH"])*" (°C per decade)"),breaks = seq(0,0.7,0.1),
                     limits = c(0.1,0.7)) +
  scale_y_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)"),breaks = seq(0, 10, 1)) +
  coord_cartesian(xlim = c(0.1,0.7), ylim = c(2,8),expand = FALSE)+
  theme_figure1+
  theme(
    plot.margin = margin(10, 20, 10, 10, "pt")
  )

## Plot
p2 = ggplot() +
  # geom_rect(
  #   fill = Sel_Colors[2],
  #   alpha = 0.2,
  #   aes(
  #     xmin = Obs_IVAR_Summary$L1[2],
  #     xmax = Obs_IVAR_Summary$U1[2],
  #     ymin = 0,
  #     ymax = 50
  #   )
  # ) +
  geom_vline(xintercept = Obs_IVAR_Summary$Obs[2], color = Sel_Colors[2], linewidth = 0.6) +
  geom_point(data = Final_EC, color = Sel_Colors[3], 
             aes(x = x_2, y = y_1)) +
  geom_smooth(
    data = Final_EC,
    color = Sel_Colors[3],
    method = "lm",
    se = FALSE,
    linewidth = 0.6,
    linetype = "dashed",
    fullrange = TRUE,
    aes(x = x_2, y = y_1)
  ) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.2, vjust = -0.5, label = expression(italic(r)* "= 0.567"))+
  scale_x_continuous(name = expression("1984-2014 "*italic(T["Trend,SH"])*" (°C per decade)"),breaks = seq(0,0.7,0.1),
                     limits = c(0.0,0.4)) +
  scale_y_continuous(name = expression("2081-2100 "*Delta*italic(T["GL"])*" (°C relative to 1850-1900)"),breaks = seq(0, 10, 1)) +
  coord_cartesian(xlim = c(0.0,0.4), ylim = c(2,8),expand = FALSE)+
  theme_figure1+
  theme(
    plot.margin = margin(10, 20, 10, 10, "pt")
  )

## Total Plot
p.total = plot_grid(p1,p2,align = "h",labels = "auto")

## Saving
ggsave(filename = "Figures/Final Figures/Fig_Temp_2.tiff",
       plot = p.total,
       width = 9, 
       height = 4.5)
ggsave(filename = "Figures/Final Figures/Fig_Temp_2.pdf",
       plot = p.total,
       device = cairo_pdf,
       width = 9, 
       height = 4.5)

# End -------------------------------------------------------------


