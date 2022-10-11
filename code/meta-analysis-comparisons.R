###########################################################
### Compare Firmicutes and Proteobacteria in other      ###
## meta-analysis datasets - Delong 2010 and Smith 2019
##
## Reproduce fig 5
##
## Tom Smith 2022
#########################################################

# load packages

library(dplyr)
library(ggpubr)
library(ggplot2)

# load data
delong_smith_data <- read.csv("data/delong-and-smith-data.csv")

# get Topt data directly from git repository
Topt_data <- read.csv(url("https://raw.githubusercontent.com/smithtp/hotterbetterprokaryotes/master/Data/summaries/summary.csv")) %>%
  filter(ConPhylum %in% c("Firmicutes", "Proteobacteria"), Points_After_Peak > 0)

########################
## Do the stats tests
########################

# Wilcoxon rank sum tests because underlying data isn't normally distributed

compare_means(Rate ~ Phylum, data = delong_smith_data %>% filter(RateType == "Active Metabolic Rate", Dataset == "DeLong2010"))
compare_means(Rate ~ Phylum, data = delong_smith_data %>% filter(RateType == "Passive Metabolic Rate", Dataset == "DeLong2010"))
compare_means(Rate ~ Phylum, data = delong_smith_data %>% filter(RateType == "20C Growth Rate", Dataset == "DeLong2010"))
compare_means(Rate ~ Phylum, data = delong_smith_data %>% filter(RateType == "20C Growth Rate", Dataset == "Smith2019"))


#####################
## Make the plots 
#####################

# main plotting theme
main_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        plot.title = element_text(size=18, vjust=1),
        legend.text=element_text(size=18),
        legend.title = element_text(size = 18),
        strip.text.x = element_text(size = 16))


# Fig 5 A
delong_metabolic <- ggplot(delong_smith_data %>% filter(RateType %in% c("Active Metabolic Rate", "Passive Metabolic Rate"), 
                                                                            Dataset == "DeLong2010"), 
                               aes(x = Phylum, y = log(Rate), fill = Phylum)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + 
  geom_point(size = 5, shape = 21, position = position_jitterdodge(), alpha = 0.6) +
  main_theme +
  scale_fill_manual(values = c("blue", "orange")) +
  ylim(-10, 2) +
  ylab("log(mass specific metabolic rate)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank()) +
  facet_wrap(~RateType)
delong_metabolic

# Fig 5 B
delong_rmax <- ggplot(delong_smith_data %>% filter(RateType == "20C Growth Rate", Dataset == "DeLong2010"), 
                      aes(x = Phylum, y = log(Rate), fill = Phylum)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + 
  geom_point(size = 5, shape = 21, position = position_jitterdodge(), alpha = 0.6) +
  main_theme +
  ylim(c(0, 5)) +
  scale_fill_manual(values = c("blue", "orange")) +
  ylab("log(Specific Growth Rate)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank())
delong_rmax

# Fig 5 C
smith_rmax <- ggplot(delong_smith_data %>% filter(RateType == "20C Growth Rate", Dataset == "Smith2019"), 
                     aes(x = Phylum, y = log(Rate), fill = Phylum)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) + 
  geom_point(size = 5, shape = 21, position = position_jitterdodge(), alpha = 0.6) +
  main_theme +
  ylim(-5, 5) +
  scale_fill_manual(values = c("blue", "orange")) +
  ylab("log(Specific Growth Rate)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90),
        axis.title.x = element_blank())
smith_rmax

# Fig 5 D
smith_Topt <- ggplot(Topt_data, aes(x = T_pk-273.15, fill = ConPhylum)) + 
  geom_histogram(position = "identity", alpha = 0.6) +
  main_theme +
  scale_fill_manual(values = c("blue", "orange")) +
  xlab("Optimum Growth Temperature") +
  geom_vline(xintercept = 40.5, linetype = "dotted") +
  theme(legend.position = c(0.7, 0.8),
        legend.title = element_blank())
smith_Topt 

# save plots

ggsave("Results/DeLong_metabolic.svg", delong_metabolic, width = 6, height = 6)
ggsave("Results/DeLong_growth.svg", delong_rmax, width = 3, height = 6)
ggsave("Results/Smith_growth.svg", smith_rmax, width = 3, height = 6)
ggsave("Results/Smith_Tpk.svg", smith_Topt, width = 6, height = 6)