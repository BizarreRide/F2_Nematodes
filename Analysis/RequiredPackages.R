###########################
# F2 Nematodes
# Required Packages
# Quentin Schorpp
# 06.08.2015
###########################

library(vegan)
library(faraway)
library(ggplot2)
library(grid)
library(ggbiplot)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}



