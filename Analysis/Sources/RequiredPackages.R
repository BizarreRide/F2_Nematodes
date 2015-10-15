###########################
# F2 Nematodes
# Required Packages
# Quentin Schorpp
# 06.08.2015
###########################

# Packages:
library(vegan)
#library(faraway)
library(ggplot2)
library(grid)
#library(ggbiplot)
library(multcomp)
library(nlme)
library(lme4)
library(car)
library(extrafont)
library(BiodiversityR)
library(reshape2)
library(dplyr)
library(tidyr)

# Functions:
source("Analysis/sources/cleanplot.pca.R")
source("Analysis/sources/evplot.R")
source("Analysis/Sources/panelutils.R")


veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}


#### ggplot2 - Theme ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mytheme = 
  theme_bw() + 
  theme(strip.background = element_rect(color = "grey", fill="black", size=0.1),
        strip.text.x = element_text(size=8,  colour="white", face="italic"),
        strip.text.y = element_text(size=8,  colour="white", face="italic"),
        axis.text.x = element_text(size=7),
        axis.title.x = element_text(size=8, family="Times New Roman"),
        axis.text.y = element_text(size=7,  family="Times New Roman"),
        axis.title.y = element_text(size=8, family="Times New Roman"),
        axis.line = element_line(size=0.25),
        axis.ticks = element_line(size=0.25),
        plot.title = element_text(size=11,face="bold", family="Times New Roman"),
        panel.margin = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", size=0.2, fill=NA),
        legend.key=element_blank(),
        legend.background=element_blank(),
        legend.text=element_text(size=8,face="italic", family="Times New Roman"),
        legend.title=element_text(size=8))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Additional functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Indices ####
# Calculate and add biodiversity Indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

biodiv.fun <- function(data) data.frame(
  SR = rowSums(data >0),                              # families Richness
  rarefy = vegan::rarefy(data,90,se=F, MARGIN=1),     # rarefaction
  H = vegan::diversity(data, index="shannon"),        # Shannon entropy
  D = vegan::diversity(data, index="simpson"),        # simpson dominance
  J = vegan::diversity(data, index="shannon")/log(rowSums(data >0)),     # Pielou Evenness
  H1 = exp(vegan::diversity(data, index="shannon")),  # Hill's N1
  N = rowSums(data))                     

#************************************************************

# Nematode Channel Ratio (NCR) ####
NCR.fun <- function(data,bacterivore,fungivore) data$fungivore/(data$bacterivore+data$fungivore)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

