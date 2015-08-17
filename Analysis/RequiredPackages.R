###########################
# F2 Nematodes
# Required Packages
# Quentin Schorpp
# 06.08.2015
###########################

# Packages:
library(vegan)
library(faraway)
library(ggplot2)
library(grid)
#library(ggbiplot)
library(multcomp)
library(nlme)
library(car)

# Functions:
source("analysis/cleanplot.pca.R")
source("analysis/evplot.R")

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}



