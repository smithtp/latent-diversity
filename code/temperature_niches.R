####################################################################
## Figures for comparing temperature profiles of isolated strains ##
##
## Reproduce figure 2 and related analyses from manuscript
##
## Tom Smith 2021
####################################################################

rm(list = ls())
library("ggplot2")
library(grid)
library(gridExtra)
dataset <- read.csv("data/growth_summary.csv")


main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title=element_text(size = 18),
        strip.text.x = element_text(size = 16))


# best_fits <- dataset[dataset$r_sq_sch > 0.5,]
best_fits <- dataset[!is.na(dataset$strain),]

assign("k", 8.617 * 10^-5, envir = .GlobalEnv) # incase we need boltzmann


#############################################
### does isolation temperature affect Tpk? ##
#############################################

# first add in the streptomyces OD fits
strep_fits <- read.csv("data/summary_OD_R.csv")

best_fits$iso_temp <- as.character(best_fits$iso_temp)
strep_fits$iso_temp <- as.character(strep_fits$iso_temp)

best_fits <- bind_rows(best_fits, strep_fits)

# add a column seperating the iso temp and standard temp experiments
best_fits$iso_incu <- NA

best_fits[best_fits$iso_temp == "RT",]$iso_incu <- "RT"
best_fits[best_fits$iso_temp != "RT",]$iso_incu <- "Incu"

isotemp_data <- best_fits[best_fits$iso_incu == "Incu",]
isotemp_data$iso_temp <- as.numeric(as.character(isotemp_data$iso_temp))
isotemp_data <- isotemp_data %>%
  filter(temps_after_peak > 0)

isotemp_plot <- ggplot(isotemp_data, aes(x = iso_temp, y = T_pk_est_sch-273.15)) +
  geom_point(size = 3) +
  geom_abline(slope = 1) +
  xlim(0, 60) +
  ylim(0, 60) +
  xlab("Isolation Temperature (Celsius)") +
  ylab("Peak Growth Temperature (Celsius)") +
  main_theme +
  geom_smooth(data = isotemp_data[isotemp_data$iso_temp > 20,], method = lm, fill = "blue") +
  geom_smooth(method = lm, color = "red", fill = "red") +
  annotate("text", x = 1, y = 58, label = "A", size = 10)
isotemp_plot

isotemp_plot_poly <- ggplot(isotemp_data, aes(x = iso_temp, y = T_pk_est_sch-273.15)) +
  geom_point(size = 4, aes(fill = Phylum), shape = 21, alpha = 0.6) +
  scale_fill_manual(values = c("red", "blue", "orange")) +
  geom_abline(slope = 1, linetype = "dotted") +
  xlim(0, 60) +
  ylim(0, 60) +
  #xlab("Isolation Temperature (Celsius)") +
  #ylab("Peak Growth Temperature (Celsius)") +
  main_theme +
  geom_smooth(data = isotemp_data, method = lm, formula = y ~ poly(x, 2, raw=TRUE), colour = "black", fill = "grey") +
  #geom_smooth(method = lm, color = "black", fill = "black") +
  #annotate("text", x = 1, y = 58, label = "A", size = 10) +
  theme(legend.position = c(0.75, 0.2),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
isotemp_plot_poly

ggsave("../Results/Figures/IsoTemp_Tpk_poly.svg", isotemp_plot_poly, width = 5, height = 5)

lm_isotemp_all <- lm(isotemp_data$T_pk_est_sch ~ isotemp_data$iso_temp)
summary(lm_isotemp_all)

lm_isotemp_high <- lm(isotemp_data[isotemp_data$iso_temp > 20,]$T_pk_est_sch ~ isotemp_data[isotemp_data$iso_temp > 20,]$iso_temp)
summary(lm_isotemp_high)

poly_isotemp <- lm(isotemp_data$T_pk_est_sch ~ poly(isotemp_data$iso_temp, 2))
summary(poly_isotemp)
plot(poly_isotemp)

anova(lm_isotemp_all, poly_isotemp)


############################################################
## did pre-incubation change Tpk for RT isolated strains? ##
############################################################

incutemp_data <- best_fits[best_fits$iso_incu == "RT",]

incutemp_plot <- ggplot(incutemp_data, aes(x = incu_temp, y = T_pk_est_sch-273.15)) +
  geom_point(size = 4, aes(fill = Phylum), alpha = 0.6, shape = 21) +
  scale_fill_manual(values = c("red", "blue", "orange")) +
  xlim(0, 60) +
  ylim(0, 60) +
  # xlab("Incubation Temperature (Celsius)") +
  # ylab("Peak Growth Temperature (Celsius)") +
  geom_smooth(method = lm, col = "black") +
  main_theme +
  #annotate("text", x = 1, y = 59.5, label = "B", size = 10) +
  theme(legend.position = c(0.75, 0.85),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
incutemp_plot

ggsave("../Results/Figures/IncuTemp_Tpk.svg", incutemp_plot, width = 5, height = 5)

lm_incutemp <- lm(incutemp_data$T_pk_est_sch ~ incutemp_data$incu_temp)
summary(lm_incutemp)

g2 <- arrangeGrob(isotemp_plot_poly, incutemp_plot, ncol=2, nrow=1, widths=c(3, 3), heights=4, clip = TRUE)
ggsave(file="../Results/Figures/incu_vs_iso_plots.png", g2, width = 12, height = 6) #saves g
