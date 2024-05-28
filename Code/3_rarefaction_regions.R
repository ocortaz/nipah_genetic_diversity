##################################################################
##                                                              ##
##              RAREFACTION ANALYSIS PER REGION                 ##
##                                                              ##
## Author: Oscar Cortes Azuero                                  ## 
## Last update : 15/11/2022                                     ##
##                                                              ##  
##################################################################

## Load necessary packages
library(dplyr)
library(ggplot2)
library(iNEXT)

set.seed(199)

## Importing data, labelling regions more clearly and keeping only regions with >10 sequences
metadata <- read.csv('../Data/seq_metadata.csv', stringsAsFactors = F)
metadata$short_province[which(metadata$short_country=='THA')] <- metadata$short_region[which(metadata$short_country=='THA')]
metadata$short_region <- paste(metadata$short_country, metadata$short_province, sep="_")
metadata <- metadata[metadata$short_region %in% c("KHM_KA", "BGD_DH", "BGD_RJ", "BGD_RN", "THA_CE", "THA_EA"),]

###### Conducting Rarefaction Analysis for each region #############################
ab_region <- as.data.frame.matrix(table(metadata$cluster, metadata$short_region))
reg_names <- c("BGD - Dhaka", "BGD - Rajshahi", "BGD - Rangpur",  "KHM - Kandal", "THA - Central", "THA - Eastern")
colnames(ab_region) <- reg_names

reg_raw_est <- iNEXT(ab_region, q=0, datatype="abundance", se=T, endpoint=152)
reg_est <- reg_raw_est$iNextEst$size_based
asy_est <- reg_raw_est$AsyEst

reg_est$Method <- factor(reg_est$Method,
                         levels = c("Rarefaction", "Observed", "Extrapolation"))

reg_est$Assemblage <- factor(reg_est$Assemblage, 
                               levels = c( "BGD - Rajshahi",  "BGD - Dhaka", "KHM - Kandal", "BGD - Rangpur", "THA - Central", "THA - Eastern"))

## Plotting diversity estimates, lines only, for all regions in one plot 
reg_labels <- data.frame(Region=unique(reg_est$Assemblage))
reg_labels$x <- c(45, 45, 35, 60, 30, 45)
reg_labels$y <- c(5, 3, 1.5, 3, 4, 7)
col_pal_1 <- c("#2271B2",  "#3DB7E9",  "#F748A5",  "#359B73", "#d55e00", "#e69f00")
g <- ggplot(reg_est) +
  geom_line(aes(x=m, y=qD, color=Assemblage, linetype=Method), data=reg_est[which(reg_est$Method!="Observed"),], linewidth=0.8) +
  geom_point(aes(x=m, y=qD, color=Assemblage), data=reg_est[which(reg_est$Method == "Observed"),], shape=18, size=3) +
  theme_bw(base_size=14) + labs(x="Number of sequences", y="Number of unique genetic clusters observed", linetype="Type", color="") +
  geom_label(data=reg_labels, aes(x=x, y=y, label=Region, color=Region), size=3) +
  # guides(fill=guide_legend(nrow=3)) +
  scale_colour_manual(values=col_pal_1) + #met.brewer(name="Archambault", n=6, type="discrete")
  scale_x_continuous(limits=c(0,100), breaks=c(0, 20, 40, 60, 80, 100)) +
  scale_y_continuous(limits=c(1,8), breaks=seq(2,8,2)) +
  theme(legend.position = "none", axis.text.x=element_text(colour="black"), 
        axis.text.y=element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
g

## Diversity estimates with ribbons in a grid
g <- ggplot(reg_est) + facet_wrap(. ~ Assemblage, ncol=3) +
  geom_ribbon(aes(x=m, ymin=qD.LCL, ymax=qD.UCL, fill=Assemblage), alpha=0.4) +
  geom_line(aes(x=m, y=qD, color=Assemblage, linetype=Method), data=reg_est[which(reg_est$Method!="Observed"),], linewidth=0.8) +
  geom_point(aes(x=m, y=qD, color=Assemblage), data=reg_est[which(reg_est$Method == "Observed"),], shape=18, size=3) +
  theme_bw(base_size=14) + labs(x="Number of sequences", y="Number of unique genetic clusters observed", linetype="Type", color="") +
  scale_colour_manual(values=col_pal_1, guide='none') + #met.brewer(name="Archambault", n=6, type="discrete")
  scale_fill_manual(values=col_pal_1, guide = 'none') +
  scale_x_continuous(limits=c(0,150), breaks=seq(0,150,50)) +
  scale_y_continuous(limits=c(1,12), breaks=seq(3,12,3)) +
  theme(legend.position = "bottom", axis.text.x=element_text(colour="black"), 
        axis.text.y=element_text(colour="black"), axis.ticks = element_line(colour = 'black'),
        strip.background = element_rect(colour="black", fill="white", linewidth=1.5, linetype="solid"),
        strip.text = element_text(size=15, colour="black"))
g

## Linear regressions of richness in function of pop density, veg index and mean distance ######################################################################
library(dplyr)
library(stringr)
library(MASS)

asy_est <- readRDS("../Data/asymptotic_rarefaction_estimates_per_region_with_environmental_data.rda")

## Model on Mean Distance ####
model <- glm(Estimator ~ mean_dist, data=asy_est, family = poisson(link="log"))

sim_data <- data.frame(mean_dist=seq(0, 65000, by=5e3))
sim_data$fit <- exp(coef(model)[1] + coef(model)[2]*sim_data$mean_dist)

coeffs_conf_int <- as.data.frame(t(confint(model, method="boot")))
sim_data$lwr <- exp(coeffs_conf_int$`(Intercept)`[1] + coeffs_conf_int$mean_dist[1]*sim_data$mean_dist)
sim_data$upr <- exp(coeffs_conf_int$`(Intercept)`[2] + coeffs_conf_int$mean_dist[2]*sim_data$mean_dist)
sim_data$mean_dist <- sim_data$mean_dist/1e3
asy_est$Assemblage <- factor(asy_est$Assemblage, 
                             levels = c( "BGD - Rajshahi",  "BGD - Dhaka", "KHM - Kandal", "BGD - Rangpur", "THA - Central", "THA - Eastern"))

asy_est$mean_dist <- asy_est$mean_dist/1e3

g <- ggplot(asy_est, aes(x=mean_dist)) +
  geom_ribbon(data=sim_data, aes(x=mean_dist, ymin=lwr, ymax=upr), alpha=0.1, fill="black") +
  geom_line(data=sim_data, aes(x=mean_dist, y=fit), size=0.8) +
  geom_point(aes(y=Estimator, color=Assemblage),   shape=18, size=3) +
  geom_linerange(aes(y=Estimator, ymin=LCL, ymax=UCL, color=Assemblage), linewidth=0.8) +
  theme_bw(base_size=14) + scale_y_continuous(breaks=seq(0,20,5), limits = c(0, 20)) +
  theme(legend.position = "none", axis.text = element_text(colour="black"), axis.ticks = element_line(colour="black")) +
  labs(x="Mean Pairwise Distance (km)",  y="Number of genetic clusters") +
  scale_color_manual(values=col_pal_1)
  # scale_y_continuous(breaks = c(0,5,10)) + scale_x_continuous(breaks = c(0.45, 0.55, 0.65))
g

## Model on Human Population Density ####
model <- glm(Estimator ~ pop_dens , data=asy_est, family = poisson(link="log"))

sim_data <- data.frame(pop_dens =seq(0, 2200, by=100))
sim_data$fit <- exp(coef(model)[1] + coef(model)[2]*sim_data$pop_dens)
coeffs_conf_int <- as.data.frame(t(confint(model, method="boot")))
sim_data$lwr <- exp(coeffs_conf_int$`(Intercept)`[1] + coeffs_conf_int$pop_dens[1]*sim_data$pop_dens)
sim_data$upr <- exp(coeffs_conf_int$`(Intercept)`[2] + coeffs_conf_int$pop_dens[2]*sim_data$pop_dens)

g <- ggplot(asy_est, aes(x=pop_dens)) +
  geom_ribbon(data=sim_data, aes(x=pop_dens, ymin=lwr, ymax=upr), alpha=0.1, fill="black") +
  geom_line(data=sim_data, aes(x=pop_dens, y=fit), size=0.8) +
  geom_point(aes(y=Estimator, color=Assemblage),   shape=18, size=3) +
  geom_linerange(aes(y=Estimator, ymin=LCL, ymax=UCL, color=Assemblage), linewidth=0.8) +
  theme_bw(base_size=14) + scale_y_continuous(breaks=seq(0,20,5), limits=c(0, 20)) +
  scale_x_continuous(breaks=seq(0, 2e3, 1e3)) +
  theme(legend.position = "none", axis.text = element_text(colour="black"), axis.ticks = element_line(colour="black")) +
  labs(x="Human Population Density",  y="Number of genetic clusters") +
  scale_color_manual(values=col_pal_1)
# scale_y_continuous(breaks = c(0,5,10)) + scale_x_continuous(breaks = c(0.45, 0.55, 0.65))
g


## Model on Tree Cover ####
model <- glm(Estimator ~ tree_cover , data=asy_est, family = poisson(link="log"))

sim_data <- data.frame(tree_cover =seq(0, 0.35, by=0.01))
sim_data$fit <- exp(coef(model)[1] + coef(model)[2]*sim_data$tree_cover)

coeffs_conf_int <- as.data.frame(t(confint(model, method="boot")))
sim_data$lwr <- exp(coeffs_conf_int$`(Intercept)`[1] + coeffs_conf_int$tree_cover[1]*sim_data$tree_cover)
sim_data$upr <- exp(coeffs_conf_int$`(Intercept)`[2] + coeffs_conf_int$tree_cover[2]*sim_data$tree_cover)

sim_data$tree_cover <- sim_data$tree_cover*100
asy_est$tree_cover <- asy_est$tree_cover*100
sim_data$upr[sim_data$upr >= 20] <- 20
sim_data$fit[sim_data$fit >= 20] <- 20

g <- ggplot(asy_est, aes(x=tree_cover)) +
  geom_ribbon(data=sim_data, aes(x=tree_cover, ymin=lwr, ymax=upr), alpha=0.1, fill="black") +
  geom_line(data=sim_data, aes(x=tree_cover, y=fit), size=0.8) +
  geom_point(aes(y=Estimator, color=Assemblage),   shape=18, size=3) +
  geom_linerange(aes(y=Estimator, ymin=LCL, ymax=UCL, color=Assemblage), linewidth=0.8) +
  theme_bw(base_size=14) +
  theme(legend.position = "bottom", axis.text = element_text(colour="black"), axis.ticks = element_line(colour="black")) +
  labs(x="% Tree Cover",  y="Number of genetic clusters") +
  scale_color_manual(values=col_pal_1)
# scale_y_continuous(breaks = c(0,5,10)) + scale_x_continuous(breaks = c(0.45, 0.55, 0.65))
g



