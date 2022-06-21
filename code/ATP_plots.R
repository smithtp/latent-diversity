#########################################################
### Compare ATP with Growth and Respirtion Rates      ###
##
## Reproduce fig 4B
##
## Tom Smith 2022
#########################################################

rm(list = ls())
library("ggplot2")
library("gridExtra")
library(dplyr)
library(tidyr)
#source("multiplot.R")


main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.text.x = element_text(size = 16))



# read the data
ATP_data <- read.csv("data/ATP-trait-data.csv")

## do the linear models

lm_normal_firmi <- lm(ATP_per_biomass_all ~ Resp_rate_per_biomass, data = ATP_data %>% filter(Phylum == "Firmicutes"))
summary(lm_normal_firmi)

lm_log_firmi <- lm(log(ATP_per_biomass_all) ~ log(Resp_rate_per_biomass), data =  ATP_data %>% filter(Phylum == "Firmicutes"))
summary(lm_log_firmi)

lm_normal_proteo <- lm(ATP_per_biomass_all ~ Resp_rate_per_biomass, data =  ATP_data %>% filter(Phylum == "Proteobacteria"))
summary(lm_normal_proteo)

lm_log_proteo <- lm(log(ATP_per_biomass_all) ~ log(Resp_rate_per_biomass), data =  ATP_data %>% filter(Phylum == "Proteobacteria"))
summary(lm_log_proteo)

lm_log_actino <- lm(log(ATP_per_biomass_all) ~ log(Resp_rate_per_biomass), data =  ATP_data %>% filter(Phylum == "Actinobacteria"))
summary(lm_log_actino)

confint(lm_log_proteo)
confint(lm_log_firmi)
confint(lm_log_actino)

# power-law relationships
# lets plot proteobacteria and firmicutes seperately
# to see the difference in magnitude

ATP_vs_resp <- ggplot(ATP_data %>% filter(Phylum != "Actinobacteria"), 
       aes(x = log(Resp_rate_per_biomass), y = log(ATP_per_biomass_all), col = Phylum)) +
  geom_point(size = 3, alpha = 0.5) +
  #ylim(0, 40) +
  scale_color_manual(values=c("blue", "orange")) +
  geom_smooth(method = lm) +
  #geom_abline(slope = 1) +
  labs(x = expression(paste("log(Respiration Rate (", h^-1, "))", sep = "")),
       y = expression(paste("log(nmols ATP/", mu, "g C biomass)", sep = ""))) +
  main_theme +
  theme(legend.title = element_blank(),
        legend.position = "none",
        aspect.ratio = 1)
ATP_vs_resp

ggsave("results/ATP_vs_respiration_phyla_combined.svg", ATP_vs_resp,  height = 6, width = 7, dpi = 300)

