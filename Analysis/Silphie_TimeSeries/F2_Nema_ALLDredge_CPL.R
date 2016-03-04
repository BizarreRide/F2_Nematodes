#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Nematodes
# GLMM with Dredge
# Quentin Schorpp
# 12.02.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") # Load original datasets nema, env, counts and define factors, etc.
rm(counts.env,counts.org,counts3SC,env.org,nema) # remove unnecessary datasets

source("Data/DataProcessing/EnvDataProcessing.R") # mainly slice and categorize, extract orthogonal variables in subset env.fin
env1$ACSC <- interaction(env1$age_class,env1$samcam)

source("Data/DataProcessing/FamDatProcessing.R") # upscaled (counts * relative abundance in subset of 100 Ind.), presene absence data and bioiv indices on upscales data 
biodiv <- biodiv.fun(round(fam.usc,0))
biodiv$N2 <- rowSums(fam.org)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.1 Data procsessing ####

# 1.1.1 Explanatory Variables ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# environmental subset
env <- data.frame(ID = 1:nrow(env1),
                  age_class=env1$age_class,
                  age = env1$age,
                  cpl = env1$cpl,
                  samcam = env1$samcam,
                  nsamcam = as.numeric(factor(env1$samcam)),
                  ACSC <- env1$ACSC,
                  location=env1$location,
                  field.ID = env1$field.ID,
                  cnratio = env1$cn,
                  mc = env1$mc,
                  pH = env1$pH,
                  ats1 = env1$ats1)

df.exp <- droplevels(env[!env$age_class %in% "A_Cm",])
row.names(df.exp) <- 1:nrow(df.exp)

# Vector of all variables of interest to define "q"
explanatory <- c("age_class", "samcam", "age_class:samcam")#, "pH", "mc", "cnratio", "ats1") # include "intercept" when using Anova type III
q <- length(explanatory)#+1

# change factor properties (optional) 
df.exp$age_class2 <- df.exp$age_class
df.exp$age_class <- as.ordered(df.exp$age_class)
df.exp$age_class <- df.exp$age_class2

df.exp$samcam2 <- df.exp$samcam
df.exp$samcam <- df.exp$nsamcam
df.exp$samcam <- df.exp$samcam2

str(df.exp)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1.1.2 Response Variables ####
# A. Nematode Indices ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.org
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FaunalProfileIndices.R") 

df.indices <- cbind(FaPro[,-c(1:5)], MaturityIndices, biodiv[,-8])
df.rsp1 <- df.indices[!env$age_class %in% "A_Cm",]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# B. FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=fam.org
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores  <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


df.rsp2 <- fety[!env$age_class %in% "A_Cm",-c(1:4)]
df.rsp3 <- df.rsp2/biodiv$N2[!env$age_class %in% "A_Cm"]     # Percentage data
df.rsp4 <- round(df.rsp3*counts[!env$age_class %in% "A_Cm", "counts"],0)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# C. Selected Taxa ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Upscaled Abundance
df.taxa <- round(fam.usc[,c("Tylenchidae", "Aphelenchidae", "Hoplolaimidae", "Cephalobidae", "Plectidae", "Telotylenchidae", "Rhabditidae", "Aporcelaimidae", "Aphelenchoididae", "Panagrolaimidae")],0)
df.rsp5 <- df.taxa[!env$age_class%in% "A_Cm",]
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# D. Put together #### 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
df.rspX <- cbind(df.rsp1, df.rsp4, df.rsp5)
row.names(df.rspX) <- 1:nrow(df.rspX)
p <- ncol(df.rspX)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#****************************************************************************************

# 2. Analysis ####

# 2.1 Outliers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Detect Outliers 
windows(record=TRUE)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  df.exp$y <- df.rspX[,i]
  car::Boxplot(df.exp$y ~ df.exp$age_class)
  car::Boxplot(df.exp$y ~ df.exp$samcam)
  car::Boxplot(df.exp$y ~ df.exp$ACSC)
  title(names(df.rspX)[i],outer=TRUE)
}

# Vector list to drop Outliers 
outlier <- list(nema.BI <- -c(24,9), #lognormal:-c(24,9), normal: -c(9,24),boxplot: -c(24,9,15,3)
                nema.SI <- 1:24,
                nema.EI <- 1:24,
                nema.CI <- -20, #lognormal:-20,boxplot:-24
                nema.MI.5 <- 1:24,#normal: -24?,boxplot:-17
                nema.PPI <- 1:24, #normal: c(6,14,15),boxplot:-c(6,19,13,14)
                nema.MI25 <- 1:24,
                nema.sigmaMI <- 1:24,#boxplot: -c(17,3)
                nema.sigmaMI25 <- -c(5,17),#lognormal:-c(6,5,10,17),normal: -c(6,5,10,17),boxplot:-23
                nema.PPI1.10 <- -15, #lognormal:-15, normal: -c(6,14,15),boxplot:-c(13,14,19,6)
                nema.SR <- 1:24, 
                nema.rarefy <- 1:24, #lognormal:-9?, normal: ,-9?
                nema.H <- 1:24,
                nema.D <- 1:24,
                nema.J.15 <- -3, #lognormal:-3,normal: , -3
                nema.H1 <- 1:24, # normal: -9,
                nema.N <- 1:24, # negbin: -14?, poisson: -c(1,13)?,boxplot:-c(2,14,20)
                fety.fungi <- 1:24, #poisson: -20, normal: -13?,boxplot: -9
                fety.bacti <- -20, #negbin:-20, poisson: -7, normal: -c(19,20),boxplot: -c(14,20)
                fety.carni.20 <- 1:24,#normal: -c(14,20),boxplot: -c(6,9,11,14,20,17)
                fety.omni <- -10, #negbin:-14?, poisson: -c(10,17), normal: -10,boxplot: -c(2,20)
                fety.Tyli <- 1:24, #negbin:-22, boxplot:-5
                fety.herbi <- -c(2,6), #negbin: -c(6,14), poisson: -17, normal: -c(2,6),boxplot: -c(2,6,14)
                fety.herbi2 <- -c(7,21), #negbin:-2,14, normal: -2,boxplot: -c(2,6,13,14)
                fety.fungi2.25 <- 1:24, #boxplot: -5
                spec.Tyli <- 1:24, #negbin:- 22, boxplot: -5
                spec.Aph <- 1:24, #poisson: -c(4,7,20), normal: -c(7,15,20),boxplot: -9
                spec.Hop <- -4, #negbin:-2,14, normal: -2,boxplot: -c(2,6,14)
                spec.Cph <- -20, #negbin:-20, normal: -c(19,20),boxplot:-c(2,5,20,14)
                spec.Plec.30 <- -21, #negbin:-14, poisson: -21?, normal: -c(14,6,20),boxplot: -c(6,14,20)
                spec.Telo <- -6, #poisson: -10,17, normal: -4,5,6,boxplot: -6
                spec.Rha <- 1:24, #poisson: -3, normal: -14,20,boxplot: -c(5,14,20)
                spec.Apc <- -16,#poisson: -c(5,18,9,17,6,12), normal: -c(6,17,22,20),boxplot: -c(6,20,22)
                spec.Aphdd <- -15, #normal: -10,15,boxplot: -10
                spec.Pan.35 <- -9 )# negbin:-c(10,8), poisson: -c(19,20,4), normal: NA2 )

# No Outlier to check with Model validation
for (i in 1:35){
  outlier[[i]] <- 1:24
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.2 Global Models *********************************************####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Family decision strategy:
#   1. all normal
#   2. lognormal (no-count), poisson(counts)
#   3. lognormal (no-count), negbin (count)
#   4. delete outlier
#   5. check again

# 2.2.1 Model Family ####
# Decision Frame for model family
df.families <- read.delim("Data/F2_Nema_GLMMFamilies.txt", header=TRUE)  
df.families[,6] <- c(rep("lognormal",16), rep("poisson",20))
df.families[,7] <- c(rep("lognormal",16), rep("negbin",20))


r <- ncol(df.families)

# Change family assumptions manually 
df.families[,r+1] <- rep("normal", nrow(df.families))
df.families[,r+1][1] <- "normal" 
df.families[,r+1][2] <- "normal" 
df.families[,r+1][3] <- "normal"
df.families[,r+1][4] <- "lognormal" 
df.families[,r+1][5] <- "normal"
df.families[,r+1][6] <- "lognormal" 
df.families[,r+1][7] <- "normal" 
df.families[,r+1][8] <- "normal" 
df.families[,r+1][9] <- "normal" 
df.families[,r+1][10] <- "lognormal" 
df.families[,r+1][11] <- "normal" 
df.families[,r+1][12] <- "normal" 
df.families[,r+1][13] <- "normal" 
df.families[,r+1][14] <- "normal" 
df.families[,r+1][15] <- "normal" 
df.families[,r+1][16] <- "lognormal" 
df.families[,r+1][17] <- "normal" 
df.families[,r+1][18] <- "normal" 
df.families[,r+1][19] <- "quasipoisson1" 
df.families[,r+1][20] <- "quasipoisson1" 
df.families[,r+1][21] <- "quasipoisson1" 
df.families[,r+1][22] <- "quasipoisson1"  
df.families[,r+1][23] <- "normal" # -14?
df.families[,r+1][24] <- "quasipoisson1" 
df.families[,r+1][25] <- "negbin" #negbin
df.families[,r+1][26] <- "normal" #-3?
df.families[,r+1][27] <- "negbin" # negbin, -11,23
df.families[,r+1][28] <- "quasipoisson1" #-14?
df.families[,r+1][29] <- "negbin" #negbin
df.families[,r+1][30] <- "quasipoisson1" #-20?
df.families[,r+1][31] <- "negbin" #negbin -22?
df.families[,r+1][32] <- "negbin" #negbin -4?
df.families[,r+1][33] <- "quasipoisson1" 
df.families[,r+1][34] <- "normal" 
df.families[,r+1][35] <- "normal" 
df.families[,r+1][36] <- "binomial"

# 2.2.2 GLobal Model formulation ####

# Formulas
fml.glb          <- as.formula(y ~ cpl + (1|field.ID)) # + pH + cnratio + mc + ats1 
fml.glb.log      <- update(fml.glb, log1p(y) ~ .)
fml.glb.qpoisson <- update(fml.glb, . ~ . + (1|ID))
fml.glb.pql.fxd  <-update(fml.glb, . ~ . -(1|field.ID))
fml.glb.pql.rnd  <- as.formula(~ 1|field.ID) #+ pH + cnratio + mc + ats1 
# including (1|ID) is changing the residuals, mainly a decrease with fitted values occurs... 

con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = TRUE, check.conv.grad="ignore")
con2 = glmerControl(optimizer=c("bobyqa", "Nelder_Mead"),tolPwrss=1e-3)

ls.glbModels <- list()

# Which families should be used?
j=r+1 #!!!!!!!!!!!!!!!!!!!!!!!!!
j=6
for ( i in 1:ncol(df.rspX)) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  df.exp2$ID  <- 1:nrow(df.exp2)
  
  if(df.families[i,j] == "normal")
    G1 <- lmer(fml.glb, df.exp2)
  if(df.families[i,j] == "lognormal")
    G1 <- lmer(fml.glb.log, df.exp2)
  if(df.families[i,j] == "poisson") {
    G1 <- glmer(fml.glb, family = "poisson", df.exp2, control = con)
    print(overdisp_fun(G1))
  }
  if(df.families[i,j] == "quasipoisson1"){
    G1 <- glmer(fml.glb.qpoisson, family = "poisson", df.exp2, control = con)
    print(overdisp_fun(G1))
  }
  if(df.families[i,j] == "quasipoisson2")
    G1 <- glmmPQL(fml.glb.pql.fxd, fml.glb.pql.rnd, family = quasipoisson, df.exp2, control = con2)
  if(df.families[i,j] == "negbin")
    #G1 <- glmer.nb(fml.glb, df.exp2, con2)
    G1 <- glmmADMB::glmmadmb(fml.glb, family="nbinom2", df.exp2)
  
  name <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
  ls.glbModels[[i]] <- assign(name, G1)
}

# 2.2.3 [Special Case, no continuous covariates] Likelihood ratio tests by Hand ####
# for (i in 1:p) {
#   df.exp$y <- df.rspX[,i]
#   df.exp2 <- df.exp[outlier[[i]],]
#   df.exp2$ID  <- 1:nrow(df.exp2)
#   G2 <- update(ls.glbModels[[i]], .~.-age_class:samcam,df.exp2)
#   G3 <- update(G2, .~.-age_class)
#   G4 <- update(G2, .~.-samcam)
#   G5 <- update(G2, .~.-age_class - samcam)
#   aov1 <- anova(ls.glbModels[[i]],G2,G3,G4,G5)
#   aov1 <- anova(ls.glbModels[[i]],G2);aov1
#   aov1 <- anova(ls.glbModels[[i]],G5);aov1
#   aov1 <- Anova(ls.glbModels[[i]], method="Wald");aov1
#   aov1 <- anova(G2,G3);aov1
#   aov1 <- anova(G2,G4);aov1
#   aov1 <- anova(G5,G3);aov1
#   aov1 <- anova(G5,G4);aov1
#   print(aov1)
# }
# Tasks:
# FInd description of using anova with more than one model!!!

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.3 Global Model Validation ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2.3.1 Residual Plots ####
windows(record=TRUE)

for(k in 1:ncol(df.rspX)){ 
  df.exp$y <- df.rspX[,k]
  df.exp2 <- df.exp[outlier[[k]],]
  
  E1 <- resid(ls.glbModels[[k]], type="pearson")
  E2 <- resid(ls.glbModels[[k]], type="response")
  if(df.families[k,j] == "quasipoisson")
    E3 <- residuals(ls.glbModels[[k]], type="normalized")
  if(df.families[k,j] == "negbin")
    E3 <- residuals(ls.glbModels[[k]], type="pearson")
  else
    E3 <- residuals(ls.glbModels[[k]], type="deviance")
  
  F1 <- fitted(ls.glbModels[[k]], type="response")
  P1 <- predict(ls.glbModels[[k]], type="response")
  
  #   mypath <- file.path("Analysis","Silphie_TimeSeries","Model_Validation","XXX",
  #                       paste("ModVal", k,abbreviate(colnames(df.rspX),3)[k], ".jpeg", sep = ""))
  #   
  #   jpeg(file=mypath)
  
  par(mfrow=c(3,3),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0)
  )
  
  # Plot fitted vs. residuals
  p1 <- scatter.smooth(F1, E1, cex.lab = 1.5, xlab="Fitted values", ylab="Pearson Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.exp2), pos = 4)
  scatter.smooth(F1, E2, cex.lab = 1.5, xlab="Fitted values", ylab="Response Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.exp2), pos = 2)
  scatter.smooth(F1, E3, cex.lab = 1.5, xlab="Fitted values", ylab="Deviance Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.exp2), pos = 2)
  
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
# 2.3.2 Plots including observed values ####

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
  
  # plot samcam vs. residuals
  boxplot(E1 ~ df.exp$samcam[outlier[[k]]],cex.lab = 1.5, xlab="samcam", ylab="Residuals")
  
  # plot samcam vs. residuals
  boxplot(E1 ~ df.exp$ACSC[outlier[[k]]],cex.lab = 1.5, xlab="agsam", ylab="Residuals")
  
  plot(df.exp2$cpl, df.exp2$y)
  points(df.exp2$cpl,P1, col="red")
  
  plot(df.exp2$cpl, df.exp2$y)
  abline(ls.glbModels[[i]], col="red")
  
  title(colnames(df.rspX)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 3. Multimodel Inference # demo(dredge.subset)*****************####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 3.1 Generate Cluster ####
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterEvalQ(clust, c(library("lme4")))
#clusterExport(clust, varlist=c("dt.exp", "con"))

# 3.2 Subsets of models excluded from dredge (optional): ####
opo <- env[,c("mc","pH","cnratio","ats1")]
opo <- as.data.frame(opo)

# 3.3 Dredge ####
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

# 3.4 Get the best models ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ls.bestmodels <- list()

# Extract Best Models with lowest AICc
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

# 3.5 Sets of Candidate Models
df.compM2 <- vector()
for(i in 1:p){
  delta4 <- subset(ls.dredge[[i]], delta < 4)
  df.compM <- data.frame(delta4)
  df.compM[,ncol(df.compM)+1] <- rep(colnames(df.rspX)[i], length(rownames(df.compM)))
  df.compM2 <- rbind(df.compM2, df.compM)
}

# 4. Other Model Selection Technics  ####

# 4.1 Leave One Out ####
for ( i in 1:p) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  print(drop1(ls.glbModels[[i]]))
}

# 4.2 Stepwise Selection ####
step(ls.glbModels[[1]])


# 5. Use no Model selection Procedure ####
ls.bestModels <- ls.glbModels

# 6. Get P-Values for Variables from Likelihood ratio Tests ####

df.LRT <- matrix(NA,1+q,2+(2*p))
colnames(df.LRT) <- c("Env", "DF", rep(colnames(df.rspX)[1:p], each=2))
df.LRT[1,] <- c("X", "X", rep(c("CHI2", "p-value"),p))


for(i in 1:p){
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  df.exp2$ID <- 1:nrow(df.exp2)
  if(df.families[i,j] %in% c("poisson", "quasipoisson1")){
    obj.afex <- afex::mixed(ls.bestModels[[i]],  df.exp2, family="poisson", method="LRT")
    df.LRT[2:4,2+(i*2)-1] <- round(obj.afex[[1]][1:3,2],2) # Chisq
    df.LRT[2:4,2+((i*2))] <- round(obj.afex[[1]][1:3,4],3)} #p-value
  if(df.families[i,j] %in% c("normal", "lognormal")) {
    obj.afex <- afex::mixed(ls.bestModels[[i]],  df.exp2,  method="LRT")
    df.LRT[2:4,2+(i*2)-1] <- round(obj.afex[[1]][1:3,2],2) # Chisq
    df.LRT[2:4,2+((i*2))] <- round(obj.afex[[1]][1:3,4],3)} #p-value
  if(df.families[i,j] == "negbin"){
    wald1 <- car::Anova(ls.bestModels[[i]],  type="II", method="LR")
    #obj.afex <- afex::mixed(ls.bestModels[[i]], df.exp2, family=negative.binomial(theta=getME(ls.bestModels[[i]], "glmer.nb.theta")), control=con,method="LRT") # getME(ls.bestModels[[i]]); lme4:::getNBdisp(ls.bestModels[[i]]))
    df.LRT[2:4,2+((i*2)-1)] <- round(wald1$"Chisq",2)[1:3]
    df.LRT[2:4,2+(i*2)] <- round(wald1$"Pr(>Chisq)",3)[1:3]}
  #   if(df.families[i,j] == "quasipoisson1"){
  #     wald1 <- car::Anova(ls.bestModels[[i]],  type="II", method="LR");wald1
  #     df.LRT[2:4,2+((i*2)-1)] <- round(wald1$"Chisq",2)[1:3]
  #     df.LRT[2:4,2+(i*2)] <- round(wald1$"Pr(>Chisq)",3)[1:3]}
}
df.LRT[2:4,1] <- row.names(obj.afex[[1]][1])
df.LRT[2:4,2] <- obj.afex[[1]]$"Chi Df"

# 7. Post Hoc data inspection with lsmeans package *************************************####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detach("package:piecewiseSEM", unload=TRUE)
detach("package:lmerTest", unload=TRUE)
detach("package:afex", unload=TRUE)
detach("package:lsmeans", unload=TRUE)
library(lsmeans)
library(multcompView)

# 7.1 Selected contrasts ####
ls.lsmAC <- list()
df.posthocAC <- matrix(NA,12,2+(2*p))
colnames(df.posthocAC) <- c("contrast", "samcam", rep(colnames(df.rspX),each=2))

for (i in 1:p) {
  if (df.LRT[4,2+((i*2))] < 0.056) 
    lsm <- lsmeans(ls.bestModels[[i]],  ~ age_class|samcam, at=list(samcam=c("2","4")))
  else 
    lsm <- lsmeans(ls.bestModels[[i]],  ~ age_class)
  
  if(df.families[i,j] %in% c("poisson", "quasipoisson1", "negbin"))
    x1 <- contrast(regrid(lsm), "pairwise" )
  else
    x1 <- contrast(lsm,"pairwise")
  xx1 <- summary(x1, type="response")
  df.posthocAC[,2+((2*i)-1)] <- round(xx1$"estimate",2)
  df.posthocAC[,2+((2*i)-0)] <- round(xx1$"p.value",3)
  name <- paste("lsm",i,names(df.rspX)[i], sep = ".")
  ls.lsmAC[[i]] <- assign(name, lsm)
} 

df.posthocAC[,1] <- paste(xx1$"contrast")
df.posthocAC[,2] <- paste(xx1$"samcam")


#####

ls.lsmSC <- list()
df.posthocSC <- matrix(NA,4,2+(2*p))
colnames(df.posthocSC) <- c("contrast", "age_class", rep(colnames(df.rspX),each=2))

for (i in 1:p) {
  if (df.LRT[4,2+((i*2))] < 0.056) 
    lsm <- lsmeans(ls.bestModels[[i]],  ~ samcam|age_class, at=list(samcam=c("2","4")))
  else 
    lsm <- lsmeans(ls.bestModels[[i]],  ~ samcam)
  
  if(df.families[i,j] %in% c("poisson", "quasipoisson1", "negbin"))
    x1 <- contrast(regrid(lsm), "pairwise" )
  else
    x1 <- contrast(lsm,"pairwise")
  xx1 <- summary(x1, type="response");xx1
  df.posthocSC[,2+((2*i)-1)] <- round(xx1$"estimate",2)
  df.posthocSC[,2+((2*i)-0)] <- round(xx1$"p.value",3)
  name <- paste("lsm",i,names(df.rspX)[i], sep = ".")
  ls.lsmAC[[i]] <- assign(name, lsm)
} 


df.posthocSC[,1] <- paste(xx1$"contrast")
df.posthocSC[,2] <- paste(xx1$"age_class")

fam <- c("XX", "XX",rep(df.families[1:35,r+1], each=2))
result <- rbind(df.LRT,df.posthocAC,df.posthocSC,fam)

# 7.2 Multicomparisons ####



#save(list=c("ls.bestModels","ls.lsm", "ls.lsmAC", "ls.lsmSC",  "df.FpvalueR2", "df.posthoc", "df.posthocAC", "df.posthocSC"), file="Results/ANOVATables/All_LMM.rda")
# write.csv(df.posthoc, file="Results/ANOVATables/PostHoc_All_LMM.csv")
# write.csv(df.posthocAC, file="Results/ANOVATables/PostHocAC_All_LMM.csv")
# write.csv(df.posthocSC, file="Results/ANOVATables/PostHocSC_All_LMM.csv")
#  write.table(df.posthoc, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
#  write.table(df.posthocAC, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
#  write.table(df.posthocSC, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Depot ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Ease manual code writing for Models
# v1 <- c()
# for (i in 1:ncol(df.rspX)) {
#   v1[i] <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], " <- ", sep="."))
# }
# v1

mypath <- file.path("Analysis","Silphie_TimeSeries","Model_Validation","PostHoc")
save(list=c("df.posthocAC.P", "df.posthocAC.QP", "df.posthocAC.all_N", "df.posthocAC.all_P_logN", "df.posthocAC.all_QP_logN", "df.posthocAC.all_NB_logN"), file="Analysis/Silphie_TimeSeries/Model_Validation/PostHoc/posthoctables.rda")
load(file="Analysis/Silphie_TimeSeries/Model_Validation/PostHoc/posthoctables.rda")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








