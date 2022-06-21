#######################################################################
# Do a PCA on all of the traits to see if the phyla cluster together
#
#
# Reproduces Fig 4A
#
# Tom Smith 2021
#############################################

rm(list = ls())

library("ggplot2")
library("plyr")
library("ggbiplot")
library("tidyr")
library("dplyr")

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


growth_summary <- read.csv("data/growth_summary.csv")
# rename some columns
names(growth_summary)[4] <- "E"
names(growth_summary)[7] <- "Topt"
names(growth_summary)[9] <- "Growth"

growth_summary <- growth_summary[growth_summary$temps_after_peak > 0,]

# add the respiration data
resp_summary <- read.csv("data/summary_resp_R_biomass.csv")
resp_summary <- resp_summary[resp_summary$temps_after_peak > 0,]

names(resp_summary)[9] <- "Respiration"

# join
combined_data <- merge(growth_summary, resp_summary[,c("strain", "Respiration")], by = "strain")

# Add ATP?
# load all data
ATP_data <- read.csv("data/ATP-trait-data.csv")

# get average ATP for each strain
meanATP <- data.frame(ATP_data[,c("Strain", "ATP_per_biomass_all")] %>%
  group_by(Strain) %>%
    dplyr::summarize(mean_ATP = mean(ATP_per_biomass_all, na.rm = TRUE)))

# join with the original data
combined_data <- merge(combined_data, meanATP, by.x = "strain", by.y = "Strain")

combined_data$logATP <- log(combined_data$mean_ATP)

# add CUE
combined_data$CUE <- combined_data$Growth/(combined_data$Growth+combined_data$Respiration)

# do the PCA
PCA_data <- combined_data[combined_data$Phylum != "Actinobacteria",c("Topt", "Growth", "Respiration", "E", "niche_width_pk", "operational_range", "cc_OD", "logATP", "Class", "Phylum")]

test_data <- PCA_data[,c("Topt", "Growth", "Respiration", "niche_width_pk", "cc_OD", "logATP")]
names(test_data) <- c("Topt", "Growth", "Respiration", "Niche Width", "Carrying Capacity", "log(ATP)")

# do the PCA
pca_model <- prcomp(test_data, scale. = TRUE, center = TRUE)

pca_plot <- ggbiplot(pca_model, groups = PCA_data$Phylum, varname.size = 4, varname.adjust = 1, ellipse = TRUE, choices = c(1,2)) + 
  geom_point(aes(colour = PCA_data$Phylum), alpha = 0.6, size = 3) +
  scale_colour_manual(values=c("blue", "orange")) +
  theme_bw() + 
  #xlim(c(-2, 3)) + 
  #ylim(c(-3, 3)) + 
  main_theme +
  theme(aspect.ratio=1,
        legend.position = c(0.2, 0.85),
        legend.title = element_blank())
pca_plot

pca_model
summary(pca_model)

ggsave("results/pca_plot.svg", pca_plot, height = 6, width = 7, dpi = 300)
