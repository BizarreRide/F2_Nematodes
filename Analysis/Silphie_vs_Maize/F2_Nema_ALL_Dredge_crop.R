#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Nematodes
# GLMM with Dredge
# Quentin Schorpp
# 12.02.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") # Load original datasets nema, env, counts and define factors, etc.
rm(counts.env,counts.org,counts3SC,env.org,nema) # remove unnecessary datasets

source("Data/DataProcessing/EnvDataProcessing.R") # mainly slice and categorize, extract orthogonal variables in subset env.fin
rm(climate1, climate30, env.fin, groups, mngmnt, soil, spa) # remove unnecessary datasets

source("Data/DataProcessing/FamDatProcessing.R") # upscaled (counts * relative abundance in subset of 100 Ind.), presene absence data and bioiv indices on upscales data 
source("Data/DataProcessing/AverageData.R")
biodiv.av <- biodiv.fun(round(fam.av.usc,0))
biodiv.av$N2 <- rowSums(fam.av.org)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Data procsessing ####

# Explanatory Variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# environmental subset
env <- data.frame(field.ID = env.av$field.ID,
                  age_class=env.av$age_class,
                  age = env.av$age,
                  location=env.av$location,
                  cnratio = env.av$cn,
                  mc = env.av$mc,
                  pH = env.av$pH,
                  ats1 = env.av$ats1)

df.exp <- env
row.names(df.exp) <- 1:nrow(df.exp)

# Vector of all variables of interest to define "q"
explanatory <- c("age_class") #,"pH", "mc", "cnratio", "ats1") # include "intercept" when using Anova type III
q <- length(explanatory)#+1
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Response Variables ####
# 1. Nematode Indices ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.av.org
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FaunalProfileIndices.R") 

df.indices <- cbind(FaPro[,-c(1:5)], MaturityIndices, biodiv.av[,-8])
df.rsp1 <- df.indices
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.av.org
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores  <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


df.rsp2 <- fety[,-c(1:4)]
df.rsp3 <- df.rsp2/biodiv.av$N2     # Percentage data
df.rsp4 <- round(df.rsp3*counts.av[, "counts"],0)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. Selected Taxa ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Upscaled Abundance
df.taxa <- round(fam.av.usc[,c("Tylenchidae", "Aphelenchidae", "Hoplolaimidae", "Cephalobidae", "Plectidae", "Telotylenchidae", "Rhabditidae", "Aporcelaimidae", "Aphelenchoididae", "Panagrolaimidae")],0)
df.rsp5 <- df.taxa
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Put together #### 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.rspX <- cbind(df.rsp1, df.rsp4, df.rsp5)
row.names(df.rspX) <- 1:nrow(df.rspX)
p <- ncol(df.rspX)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#****************************************************************************************

# Analysis ####

# Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Detect Outliers 
windows(record=TRUE)
for(i in 1:p) {
  par(mfrow = c(1,1),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  df.exp$y <- df.rspX[,i]
  car::Boxplot(df.exp$y ~ df.exp$age_class)
  title(names(df.rspX)[i],outer=TRUE)
}
# ? nicht unbedingt!
# normal: aus qq plots zu normal distribution (all 1:18)
# lognormal: aus qq plots zu log-normal distribution (all 1:18)
# poisson: aus qq plots zu poisson distribution (all 1:18)
# bocplot: wenn real daten boxplots zeigen outlier
# NA: keine outlier
# NA2: keine outlier, da massive abweicheung

# Vector list to drop Outliers 
outlier <- list(nema.1.BI <- 1:18, #lognormal: NA, normal: -18?
                nema.2.SI <- 1:18, #lognormal: NA,  normal: NA
                nema.3.EI <- -18, # lognormal: -18, normal: -18
                nema.4.CI <- -8, # lognormal: -8?, normal: -8?
                nema.5.MI <- 1:18, # lognormal: NA, normal: NA
                nema.PPI <- -6, # lognormal: -6?, normal: -6, 
                nema.MI25 <- 1:18, # lognormal: NA, normal: NA
                nema.sigmaMI <- 1:18, #lognormal: NA, normal: -18?
                nema.sigmaMI25 <- 1:18, # lognormal: , normal: NA
                nema.10.PPI1 <- -6, # lognormal: -6?, normal: -6, 
                nema.SR <- 1:18, # lognormal: NA, normal: -16?
                nema.rarefy <- 1:18, # lognormal: NA, normal: NA
                nema.H <- 1:18, # lognormal: NA, normal: NA
                nema.D <- -18, # lognormal: -18, normal: -18
                nema.15.J <- 1:18, #lognormal: NA, normal: -12?
                nema.H1 <- 1:18, # lognormal: NA, normal: NA
                nema.N <- -2, # poisson: -2, normal: -2
                fety.fungi <- -15, # poisson: -15 # normal: -15, boxplot: -15
                fety.bacti <- -8, # poisson: -8?,normal: -8, 
                fety.20.carni <- 1:18, # poisson: NA(2), normal: NA2, 
                fety.omni <- 1:18, # poisson: NA , normal: NA, # boxplot: -15
                fety.Tyli <- 1:18, # poisson: NA,normal: -5?, #boxplot: -15
                fety.herbi <- -15, # poisson: -15? ,normal:-2, #boxplot: -15
                fety.herbi2 <- -15, # poisson: -15?,normal:-2, #boxplot: -15
                fety.25.fungi2 <- -15, # poisson: NA,normal: NA(2), #boxplot: -15
                spec.Tyli <- 1:18, #poisson: NA,normal: -5? boxplot: -15
                spec.Aph <- 1:18, # poisson: NA,normal:-9,
                spec.Hop <- -15, # poisson: NA,normal: -2,#boxplot: -15
                spec.Cph <- 1:18, # poisson: NA,normal: NA
                spec.30.Plec <- -17,# poisson: -17,normal: NA2, boxplot: -17? 
                spec.Telo <- -15, # poisson: -c(6,15) # normal:-c(15,6), # boxplot: -15
                spec.Rha <- -8, # poisson: -8?, normal:-8,
                spec.Apc <- 1:18, # poisson: NA,normal: NA
                spec.Aphdd <- -15, # poisson: -15,normal:-15, boxplot: -15
                spec.35.Pan <- -9) # poisson: NA,normal: -9)

# No Outlier to check with Model validation
# for (i in 1:35){
#   outlier[[i]] <- 1:18
# }
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# change factor properties (optional) ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.exp$age_class2 <- df.exp$age_class
df.exp$age_class <- as.ordered(df.exp$age_class)
df.exp$age_class <- df.exp$age_class2

df.exp$samcam2 <- df.exp$samcam
df.exp$samcam <- df.exp$nsamcam
df.exp$samcam <- df.exp$samcam2

str(df.exp)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Global Models *********************************************####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Ease code writing 
# v1 <- c()
# for (i in 1:ncol(df.rspX)) {
#   v1[i] <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], " <- ", sep="."))
# }
# v1

# Decision Frame for model family
df.families <- read.delim("Data/F2_Nema_GLMMFamilies.txt", header=TRUE)  

# Change family assumptions manually #####
#? nicht unbedingt, Wackelkandidat, gucken ob alternative besser, wenn nicht dann normal
# Alle Änderungen ausgehend von Diagnostics of normal without outliers
df.families$V6 <- rep("normal", nrow(df.families))
df.families$V6[1] <- "normal" 
df.families$V6[2] <- "normal" 
df.families$V6[3] <- "normal"
df.families$V6[4] <- "normal" #difficult, residuals more sky at night with normal; but, qqplot better with lognormal
df.families$V6[5] <- "normal" # No Difference!!!!
df.families$V6[6] <- "normal" # only very small differences
df.families$V6[7] <- "normal" 
df.families$V6[8] <- "normal" 
df.families$V6[9] <- "normal" 
df.families$V6[10] <- "lognormal" # only small differences
df.families$V6[11] <- "normal" 
df.families$V6[12] <- "normal" 
df.families$V6[13] <- "normal" 
df.families$V6[14] <- "normal" # Neither nor; no adequate model
df.families$V6[15] <- "normal" 
df.families$V6[16] <- "lognormal" 
df.families$V6[17] <- "normal" 
df.families$V6[18] <- "negbin" 
df.families$V6[19] <- "negbin" 
df.families$V6[20] <- "negbin" 
df.families$V6[21] <- "negbin"  
df.families$V6[22] <- "negbin"  #?
df.families$V6[23] <- "negbin"  
df.families$V6[24] <- "negbin"  
df.families$V6[25] <- "negbin" 
df.families$V6[26] <- "negbin"  
df.families$V6[27] <- "negbin" 
df.families$V6[28] <- "negbin"  
df.families$V6[29] <- "negbin" 
df.families$V6[30] <- "negbin"  
df.families$V6[31] <- "negbin" # 6 still outlier?
df.families$V6[32] <- "negbin" 
df.families$V6[33] <- "negbin" #?
df.families$V6[34] <- "negbin"  
df.families$V6[35] <- "negbin"  
df.families$V6[36] <- "binomial"

df.families[,7] <- c(rep("lognormal",16), rep("poisson",20))

# GLobal Model formulation ####

# Formulas
fml.glb <- as.formula(y ~ age_class) # + pH + cnratio + mc + ats1 
fml.glb.log <- as.formula(log1p(y) ~ age_class) #+ pH + cnratio + mc + ats1 

con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = TRUE, check.conv.grad="ignore")

ls.glbModels <- list()
j=6
for ( i in 1:ncol(df.rspX)) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  if(df.families[i,j] == "normal")
    G1 <- lm(fml.glb, df.exp2)
  if(df.families[i,j] == "lognormal")
    G1 <- lm(fml.glb.log, df.exp2)
  if(df.families[i,j] == "poisson") #{
    G1 <- glm(fml.glb, family = "poisson", df.exp2)
  if(df.families[i,j] == "quasipoisson") #{
    G1 <- glm(fml.glb, family = "quasipoisson", df.exp2)
  if(df.families[i,j] == "negbin") #{
    G1 <- glm.nb(fml.glb, link="log", df.exp2)
#   print(1-pchisq(G1$deviance,G1$df.residual))
#   print(G1$deviance/G1$df.residual)}
  name <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
  ls.glbModels[[i]] <- assign(name, G1)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Global Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### Residual Plots ####
windows(record=TRUE)

for(k in 1:ncol(df.rspX)){ 
  df.exp$y <- df.rspX[,k]
  df.exp2 <- df.exp[outlier[[k]],]
  
  E1 <- resid(ls.glbModels[[k]], type="pearson")
  if(df.families[i,j] %in% c("normal","lognormal","poisson","negbin"))
    E2 <- resid(ls.glbModels[[k]], type="response")
  else     {
    mu <- predict(ls.glbModels[[k]], type="response")
    E <- df.exp2$y - mu
    E2 <- E /sqrt(ls.glbModels[[i]]$deviance/ls.glbModels[[33]]$df.residual*mu)
    }
  E3 <- residuals(ls.glbModels[[k]], type="deviance")
  
  F1 <- fitted(ls.glbModels[[k]], type="response")
  P1 <- predict(ls.glbModels[[k]], type="response")
  
#     mypath <- file.path("Analysis","Silphie_vs_Maize","Model_Validation","LogNormalPoisson_withoutOutliers",
#                         paste("ModVal", k,abbreviate(colnames(df.rspX),3)[k], ".jpeg", sep = ""))
#     
#     jpeg(file=mypath)
  
  par(mfrow=c(3,3),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
  
  # Plot fitted vs. residuals
  p1 <- scatter.smooth(F1, E1, cex.lab = 1.5, xlab="Fitted values", ylab="Pearson Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.exp2), pos = 4)
  scatter.smooth(P1, E2, cex.lab = 1.5, xlab="Fitted values", ylab="Response Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.exp2), pos = 2)
  scatter.smooth(F1, E3, cex.lab = 1.5, xlab="Fitted values", ylab="Deviance Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(log(predict(ls.glbModels[[k]])), E1, cex.lab = 1.5, xlab="log Predicted values", ylab="Pearson Residuals")
  abline(h = 0, v=0, lty=2); text(log(predict(ls.glbModels[[k]])), E1, labels = row.names(df.exp2), pos = 4)
  scatter.smooth(log(predict(ls.glbModels[[k]])), E2, cex.lab = 1.5, xlab="log Predicted values", ylab="Response Residuals")
  abline(h = 0, v=0, lty=2)
  scatter.smooth(log(predict(ls.glbModels[[k]])), E3, cex.lab = 1.5, xlab="log Predicted values", ylab="Deviance Residuals")
  abline(h = 0, v=0, lty=2)
  
  # Normal QQ Plots
  qq <- qqnorm(E1)
  qqline(E1)
  text(qq$x, qq$y, row.names(df.exp2), pos=4)
  
  qq <- qqnorm(E2)
  qqline(E2)
  text(qq$x, qq$y, row.names(df.exp2), pos=2)
  
  qq <- qqnorm(E3)
  qqline(E3)
  text(qq$x, qq$y, row.names(df.exp2), pos=2)
  
  title(paste(k,colnames(df.rspX)[k]), outer=TRUE)
  #dev.off()
}
### Plots including observed values ####

windows(record=TRUE)

for(k in 1:ncol(df.rspX)){ 
  E1 <- resid(ls.glbModels[[k]], type="pearson")
  E2 <- resid(ls.glbModels[[k]], type="response")
  E3 <- residuals(ls.glbModels[[k]], type="deviance")
  
  F1 <- fitted(ls.glbModels[[k]], type="response")
  P1 <- predict(ls.glbModels[[k]], type="response")
  
  par(mfrow=c(3,3),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  ) 
  
  # plot fitted vs. predicted
  plot(F1, P1, cex.lab = 1.5, xlab="Fitted values", ylab="Predicted")
  abline(h = 0, v=0, lty=2)
  
  # Histogram of Residuals
  hist(E1, prob=TRUE, main = "", breaks = 20, cex.lab = 1.5, xlab = "Response Residuals", col="PapayaWhip")
  lines(density(E1), col="light blue", lwd=3)
  lines(density(E1, adjust=2), lty="dotted", col="darkgreen", lwd=2) 
  
  df.exp$y <- df.rspX[,k]
  scatter.smooth(F1,df.exp$y[outlier[[k]]], cex.lab = 1.5, xlab="Fitted values", ylab="Original values")
  #text(F1,df.exp$y,label=interaction(df.exp$field.ID,df.exp$samcam),col='red')
  abline(h = 0, v=0, lty=2)
  
  # plot age_class vs. residuals
  boxplot(E1 ~ df.exp$age_class[outlier[[k]]], cex.lab = 1.5, xlab="age_class", ylab="Residuals")
  
  title(colnames(df.rspX)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Multimodel Inference # demo(dredge.subset)*****************####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Generate Cluster ####
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterEvalQ(clust, c(library("lme4")))
#clusterExport(clust, varlist=c("dt.exp", "con"))

# Subsets of models excluded from dredge: ####
opo <- env[,c("mc","pH","cnratio","ats1")]
opo <- as.data.frame(opo)

# Dredge ####
library(MuMIn)
options(na.action = na.fail)
p <- ncol(df.rspX)
ls.dredge <- list()

for(i in 1:p) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  #   dt.exp$y <- dt.rsp.abn[,i, with=F]
  clusterExport(clust, varlist=c("df.exp2", "con"))
  GM.dredge <- pdredge(ls.glbModels[[i]], cluster=clust)
  #fixed=c("age_class", "samcam"),
  name <- c(paste("Dredge",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
  assign(name, GM.dredge)
  ls.dredge[[i]] <- assign(name, GM.dredge)
}


# Get the best models **************************************####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ls.bestmodels <- list()

# Extract Best Models with lowes AICc
for ( i in 1:p) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  #dt.exp$y <- dt.rsp[,i, with=F]
  M.best <- get.models(ls.dredge[[i]], 1)[[1]]
  name <- c(paste("BM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
  assign(name, M.best)
  ls.bestmodels[[i]] <- assign(name, M.best)
  names(ls.bestmodels)[[i]] <- name
}

# See Model Formulae of best Models
for ( i in 1:p) {
  print(formula(ls.bestmodels[[i]]))
}

# See Sets of Candidate Models
df.compM2 <- vector()
for(i in 1:p){
  delta4 <- subset(ls.dredge[[i]], delta < 4)
  df.compM <- data.frame(delta4)
  df.compM[,ncol(df.compM)+1] <- rep(colnames(df.rspX)[i], length(rownames(df.compM)))
  df.compM2 <- rbind(df.compM2, df.compM)
}

# Other Methods alternative to Multimodel Inference ####

# Leave One Out ####
for ( i in 1:p) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  print(drop1(ls.glbModels[[i]]))
}

# Stepwise Selection ####
step(ls.glbModels[[1]])


# Use no Model selection Procedure ####
ls.bestModels <- ls.glbModels

# Get P-Values for Variables from Likelihood ratio Tests ####

df.ANOVA <- matrix(NA,1+q,2+(2*p))
colnames(df.ANOVA) <- c("Env", "DF", rep(colnames(df.rspX)[1:p], each=2))
df.ANOVA[1,] <- c("X", "X", rep(c("LR-CHI2/F", "p-value"),p))

for(i in 1:p){
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  if(df.families[i,j] %in% c("normal","lognormal")) {
    aov1 <- car::Anova(ls.bestModels[[i]], type="II")
    df.ANOVA[2,2+((i*2)-1)] <- round(aov1[[3]][1],2)
    df.ANOVA[2,2+((i*2)-0)] <- round(aov1[[4]][1],3)
  }
  if(df.families[i,j] %in% c("poisson","negbin")) {
    aov1 <- car::Anova(ls.bestModels[[i]], type="II", test.statistic="LR");aov1
    df.ANOVA[2,2+((i*2)-1)] <- round(aov1[[1]],2)
    df.ANOVA[2,2+((i*2)-0)] <- round(aov1[[3]],3)
  }
  if(df.families[i,j]=="quasipoisson") {
    aov1 <- car::Anova(ls.bestModels[[i]], type="II", test.statistic="F");aov1
    df.ANOVA[2,2+((i*2)-1)] <- round(aov1[[3]][1],2)
    df.ANOVA[2,2+((i*2)-0)] <- round(aov1[[4]][1],3)
  }}

df.ANOVA[2:q,1] <- row.names(aov1)[1]
df.ANOVA[2:q,2] <- aov1[[2]][1]


# R2  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Post Hoc ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detach("package:piecewiseSEM", unload=TRUE)
detach("package:lmerTest", unload=TRUE)
detach("package:afex", unload=TRUE)
detach("package:lsmeans", unload=TRUE)
library(lsmeans)

df.posthoc <- matrix(NA,1,1+2*p)
df.posthoc2 <- matrix(NA,4,1+2*p)

for(i in 1:p) {
  lsm <- lsmeans::lsmeans(ls.bestModels[[i]], ~ age_class);lsm
  summary(lsm, type="response")
  summary(pairs(lsm), type="response")
  summary(regrid(pairs(lsm)), type="reponse")
  lsm2 <- contrast(regrid(lsm), "trt.vs.ctrl",ref=1,adjust="fdr");lsm2
  summary(lsm2,adjust="fdr")
  cld(lsm)
  lsm2 <- contrast(regrid(lsm), "trt.vs.ctrl", ref=c(2:5));lsm2
  lsm2 <- test(contrast(regrid(lsm), "trt.vs.ctrl", ref=c(2:5)),side = "nonsup", delta = .25);lsm2
  summary(lsm2, type="response")
  summary(pairs(lsm2), type="response")
  summary(regrid(lsm2), type="response")
  print(lsm2)
  df.posthoc[,(2*i)] <- round(summary(lsm2)$estimate,3)
  df.posthoc[,(2*i)+1] <- round(summary(lsm2)$p.value,3)
  lsm2 <- contrast(regrid(lsm), "trt.vs.ctrl", ref=c(2:5)) # Why do the p-values differ if i use regrid() ???
  df.posthoc[,(2*i)] <- round(summary(lsm2)$estimate,3)
  lsm3 <- contrast(lsm, "trt.vs.ctrl", ref=c(1))
  print(lsm3)
  df.posthoc2[,(2*i)] <- round(summary(lsm3)$estimate,3)
  df.posthoc2[,(2*i)+1] <- round(summary(lsm3)$p.value,3)
  lsm3 <- contrast(regrid(lsm), "trt.vs.ctrl", ref=c(1))
  df.posthoc2[,(2*i)] <- round(summary(lsm3)$estimate,3)
}

df.posthoc[,1] <- paste(summary(lsm2)$"contrast")
df.posthoc2[,1] <- paste(summary(lsm3)$"contrast")

colnames(df.posthoc) <- c("Contrast", rep(c("estimate", "p-value"), p))
colnames(df.posthoc2) <- c("Contrast", rep(c("estimate", "p-value"), p))

# write.csv(df.posthoc, file="Results/ANOVATables/PostHocC_All_LMM_crop.csv")
# write.csv(df.posthoc2, file="Results/ANOVATables/PostHocAC_All_LMM_crop.csv")
# write.table(df.posthoc, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
# write.table(df.posthoc2, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



