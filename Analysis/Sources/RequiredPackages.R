###########################
# F2 Nematodes
# Required Packages
# Quentin Schorpp
# 06.08.2015
###########################

# Packages:
#library(vegan)
library(BiodiversityR)
#library(faraway)
library(reshape2)
library(ggplot2)
library(extrafont)
#library(ggbiplot)
library(multcomp)
library(car)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(grid)
library(influence.ME)
#library(MASS)

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
        strip.text.x = element_text(size=9,  colour="white", face="italic"),
        strip.text.y = element_text(size=9,  colour="white", face="italic"),
        axis.text.x = element_text(size=9),
        axis.title.x = element_text(size=9, family="Times New Roman"),
        axis.text.y = element_text(size=9,  family="Times New Roman"),
        axis.title.y = element_text(size=9, family="Times New Roman"),
        axis.line = element_line(size=0.25),
        axis.ticks = element_line(size=0.25),
        plot.title = element_text(size=11,face="bold", family="Times New Roman"),
        panel.margin = unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black", size=0.2, fill=NA),
        legend.key=element_blank(),
        legend.background=element_blank(),
        legend.text=element_text(size=9,face="italic", family="Times New Roman"),
        legend.title=element_text(size=9))
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

panel.cor <- function(x,y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr =c(0,1,0,1))
  r <- abs(cor(x,y, use="complete.obs"))
  txt <- format(c(r,0,123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  test <- cor.test(x,y)
  # borrowed from prinCoefmat
  Signif <- symnum(test$p.value, corr=FALSE, na=FALSE,
                   cutpoints=c(0,0.05,0.1,1),
                   symbols=c("*", "."," "))
  text(0.5, 0.5, txt, cex=cex *r)
  text(.8, .8, Signif, cex=cex, col=2)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asinTransform <- function(p) { asin(sqrt(p)) }

dispersion_glmer<- function(modelglmer)
{   
  # computing  estimated scale  ( binomial model)
  #following  D. Bates :
  #That quantity is the square root of the penalized residual sum of
  #squares divided by n, the number of observations, evaluated as:
  
  n <- length(modelglmer@resid)
  
  return(  sqrt( sum(c(modelglmer@resid, modelglmer@u) ^2) / n ) )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

