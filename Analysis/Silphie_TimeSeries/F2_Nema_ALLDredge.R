#%%%%%%%%%%%%%%%%%%%%%%%%%%%
# F2_Nematodes
# GLMM with Dredge
# Quentin Schorpp
# 12.02.2015
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Load Data ______________________________________________________________________####
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 1.1 Data procsessing ####

# 1.1.1 Explanatory Variables =======================================================####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# environmental subset
env <- data.frame(ID = 1:nrow(env1),
                  age_class=env1$age_class,
                  age = env1$age,
                  cpl = env1$cpl,
                  samcam = env1$samcam,
                  nsamcam = as.numeric(factor(env1$samcam)),
                  ACSC = env1$ACSC,
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

# 1.1.2 Response Variables ==========================================================####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. Nematode Indices ####
#****************************************************************************************
data=fam.org
source("Data/DataProcessing/MaturityIndices.R")
source("Data/DataProcessing/FaunalProfileIndices.R") 

df.indices <- cbind(FaPro[,-c(1:5)], MaturityIndices, biodiv[,-8])
df.rsp1 <- df.indices[!env$age_class %in% "A_Cm",]

# B. FeedingTypes ####
#****************************************************************************************
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

# C. Selected Taxa ####
#****************************************************************************************

# Upscaled Abundance
df.taxa <- round(fam.usc[,c("Tylenchidae", "Aphelenchidae", "Hoplolaimidae", "Cephalobidae", "Plectidae", "Telotylenchidae", "Rhabditidae", "Aporcelaimidae", "Aphelenchoididae", "Panagrolaimidae")],0)
df.rsp5 <- df.taxa[!env$age_class%in% "A_Cm",]

# D. Put together #### 
#****************************************************************************************
df.rspX <- cbind(df.rsp1, df.rsp4, df.rsp5)
row.names(df.rspX) <- 1:nrow(df.rspX)
p <- ncol(df.rspX)
n <- nrow(df.rspX)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Plot the data ####
#////////////////////////////////////////////////////////////////////////////////////////
# Plots for random effects
# library(lattice)
# i=19
# df.exp$y <- df.rspX[,i]
# xyplot(y ~ samcam | field.ID, df.exp)#, type="1", strip=FALSE)
#////////////////////////////////////////////////////////////////////////////////////////


# 2. Analysis _______________________________________________________________________####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# 2.1 Outliers ======================================================================####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# A. Detect Outliers ####
#****************************************************************************************

windows(record=TRUE)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  df.exp$y <- df.rspX[,i]
  #df.exp2 <- df.exp[outlier[[i]],]
  car::Boxplot(df.exp$y ~ df.exp$age_class)
  car::Boxplot(df.exp$y ~ df.exp$samcam)
  car::Boxplot(df.exp$y ~ df.exp$ACSC)
  title(names(df.rspX)[i],outer=TRUE)
}

# B. Drop Outliers by Vector lists #####
#****************************************************************************************

all=c(1:n)
# Outliers based on Residual plots ####
# outlier.resid <- list(nema.BI <- -c(24,9), #lognormal:-c(24,9), normal: -c(9,24),boxplot: -c(24,9,15,3)
#                 nema.SI <- all,
#                 nema.EI <- all,
#                 nema.CI <- -20, #lognormal:-20,boxplot:-24
#                 nema.MI.5 <- all,#normal: -24?,boxplot:-17
#                 nema.PPI <- all, #normal: c(6,14,15),boxplot:-c(6,19,13,14)
#                 nema.MI25 <- all,
#                 nema.sigmaMI <- all,#boxplot: -c(17,3)
#                 nema.sigmaMI25 <- -c(5,17),#lognormal:-c(6,5,10,17),normal: -c(6,5,10,17),boxplot:-23
#                 nema.PPI1.10 <- -15, #lognormal:-15, normal: -c(6,14,15),boxplot:-c(13,14,19,6)
#                 nema.SR <- all, 
#                 nema.rarefy <- all, #lognormal:-9?, normal: ,-9?
#                 nema.H <- all,
#                 nema.D <- all,
#                 nema.J.15 <- -3, #lognormal:-3,normal: , -3
#                 nema.H1 <- all, # normal: -9,
#                 nema.N <- all, # negbin: -14?, poisson: -c(1,13)?,boxplot:-c(2,14,20)
#                 fety.fungi <- all, #poisson: -20, normal: -13?,boxplot: -9
#                 fety.bacti <- -20, #negbin:-20, poisson: -7, normal: -c(19,20),boxplot: -c(14,20)
#                 fety.carni.20 <- all,#normal: -c(14,20),boxplot: -c(6,9,11,14,20,17)
#                 fety.omni <- -10, #negbin:-14?, poisson: -c(10,17), normal: -10,boxplot: -c(2,20)
#                 fety.Tyli <- all, #negbin:-22, boxplot:-5
#                 fety.herbi <- -c(2,6), #negbin: -c(6,14), poisson: -17, normal: -c(2,6),boxplot: -c(2,6,14)
#                 fety.herbi2 <- -c(7,21), #negbin:-2,14, normal: -2,boxplot: -c(2,6,13,14)
#                 fety.fungi2.25 <- all, #boxplot: -5
#                 spec.Tyli <- all, #negbin:- 22, boxplot: -5
#                 spec.Aph <- all, #poisson: -c(4,7,20), normal: -c(7,15,20),boxplot: -9
#                 spec.Hop <- 1:24, #negbin:-2,14, normal: -2,boxplot: -c(2,6,14)
#                 spec.Cph <- -20, #negbin:-20, normal: -c(19,20),boxplot:-c(2,5,20,14)
#                 spec.Plec.30 <- -21, #negbin:-14, poisson: -21?, normal: -c(14,6,20),boxplot: -c(6,14,20)
#                 spec.Telo <- -4, #poisson: -10,17, normal: -4,5,6,boxplot: -6
#                 spec.Rha <- all, #poisson: -3, normal: -14,20,boxplot: -c(5,14,20)
#                 spec.Apc <- -16,#poisson: -c(5,18,9,17,6,12), normal: -c(6,17,22,20),boxplot: -c(6,20,22)
#                 spec.Aphdd <- -15, #normal: -10,15,boxplot: -10
#                 spec.Pan.35 <- -9 )# negbin:-c(10,8), poisson: -c(19,20,4), normal: NA2, boxplot: -c(18,20)

# Outlier based on boxplots and Cooks Distances #### 
# resids = Residuals
# Cook = Cooks Distance > 4/(24-3-1)=0.2
outlier.bxplCook <- list(nema.BI <- -24, # resids: good; Cook: 24
                nema.SI <- all,# resids: good; Cook: 0
                nema.EI <- -3,# resids: good; Cook:3
                nema.CI <- -24, # resids: -20; Cook: 24
                nema.MI.5 <- all, # resids: -3,24,10; Cook: 0
                nema.PPI <- -c(6,14), # resids: -15; Cook:14,15
                nema.MI25 <- -24,# resids: -24; Cook:24
                nema.sigmaMI <- -c(17,3),# resids: good; Cook:5,17
                nema.sigmaMI25 <- -c(17,5), # resids: good; Cook:17,5
                nema.PPI1.10 <- -c(6,14), # resids: good; Cook:14,15
                nema.SR <- all, # resids: good; Cook:9
                nema.rarefy <- -9, # resids: good; Cook:9
                nema.H <- all,# resids: good; Cook:0
                nema.D <- all,# resids: good; Cook:11
                nema.J.15 <- -3, # resids: -3; Cook:3
                nema.H1 <- all, # resids: good; Cook:8,9
                nema.N <- all, # resids: good; Cook:0
                fety.fungi <- -c(8,20), # resids: -23; Cook: 8,20
                fety.bacti <- -20, # resids: -14; Cook: 7,20,19
                fety.carni.20 <- all,# resids: good; Cook: 24
                fety.omni <- -22, # resids: good; Cook:22
                fety.Tyli <- -22, # resids: -23,24; Cook:6,22
                fety.herbi <- -c(2,17), # resids: -14; Cook:17,5
                fety.herbi2 <- all, ## resids: good; Cook:7,19
                fety.fungi2.25 <- -22, # resids: good; Cook:23,22
                spec.Tyli <- -6, # resids: good; Cook:6,22
                spec.Aph <- -20, # resids: good; Cook:8,4,20
                spec.Hop <- -c(6,18), # resids: -14; Cook:18,6
                spec.Cph <- -20, # resids: -14; Cook:2,19,20
                spec.Plec.30 <- all, ## resids: -20; Cook: 6,20
                spec.Telo <- -c(24,22), # resids: good; Cook: 20,22,24
                spec.Rha <- all, # resids: good; Cook:0
                spec.Apc <- all,# resids: -16; Cook:16
                spec.Aphdd <- -10, # resids: good; Cook:11,10
                spec.Pan.35 <- -19 )# resids: -9; Cook:19

# No Outlier to check with Model validation ####
no.outlier=list()
for (i in 1:35){
  no.outlier[[i]] <- 1:24
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 2.2 Model Family =================================================================#####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Family decision strategy:
#   1. all normal
#   2. lognormal (no-count), poisson(counts)
#   3. lognormal (no-count), negbin (count)
#   4. delete outlier
#   5. check again


# A. Decision Frame for model family ####
df.families <- read.delim("Data/F2_Nema_GLMMFamilies.txt", header=TRUE)  
df.families[,6] <- c(rep("lognormal",16), rep("poisson",20))
df.families[,7] <- c(rep("normal",16), rep("quasipoisson1",20))
df.families[,8] <- c(rep("lognormal",16), rep("negbin",20))
r <- ncol(df.families)

# B. Change family assumptions manually ####
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
df.families[,r+1][17] <- "quasipoisson1" 
df.families[,r+1][18] <- "normal" 
df.families[,r+1][19] <- "quasipoisson1" 
df.families[,r+1][20] <- "quasipoisson1" 
df.families[,r+1][21] <- "quasipoisson1" 
df.families[,r+1][22] <- "quasipoisson1"  # Not Choose!
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 2.3 GLobal Models =================================================================####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.3.1 Formulas #####
#****************************************************************************************
fml.glb          <- as.formula(y ~ age_class + samcam + age_class:samcam + (1|field.ID)) # + pH + cnratio + mc + ats1 
fml.glb.log      <- update(fml.glb, log1p(y) ~ .)
fml.glb.qpoisson <- update(fml.glb, . ~ . + (1|ID))
fml.glb.pql.fxd  <-update(fml.glb, . ~ . -(1|field.ID))
fml.glb.pql.rnd  <- as.formula(~ 1|field.ID) #+ pH + cnratio + mc + ats1 
# including (1|ID) is changing the residuals, mainly a decrease with fitted values occurs... 

# 2.3.2 Optimization Functions #####
#****************************************************************************************

con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = TRUE, check.conv.grad="ignore")
con2 = glmerControl(optimizer=c("bobyqa", "Nelder_Mead"),tolPwrss=1e-3)

# 2,3,3 Calculate Models #####
#****************************************************************************************

# Which families should be used?
j=r+1 #!!!!!!!!!!!!!!!!!!!!!!!!!
outlier <- outlier.bxplCook
j=3
ls.glbModels <- list()

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
  #alt.est <- influence.ME::influence(G1, obs=TRUE)
  #print(tail(cooks.distance.estex(alt.est, sort=TRUE),3))
  name <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
  ls.glbModels[[i]] <- assign(name, G1)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.4 Global Model Validation =======================================================####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 2.3.1 Residual Plots ####
#****************************************************************************************

windows(record=TRUE)

for(i in 1:ncol(df.rspX)){ 
  
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  
  E1 <- resid(ls.glbModels[[i]], type="pearson")
  E2 <- resid(ls.glbModels[[i]], type="response")
  if(df.families[i,j] == "quasipoisson")
    E3 <- residuals(ls.glbModels[[i]], type="normalized")
  if(df.families[i,j] == "negbin")
    E3 <- residuals(ls.glbModels[[i]], type="pearson")
    else
    E3 <- residuals(ls.glbModels[[i]], type="deviance")
  
  F1 <- fitted(ls.glbModels[[i]], type="response")
  P1 <- predict(ls.glbModels[[i]], type="response")

# Save Plots extern  
#   mypath <- file.path("Analysis","Silphie_TimeSeries","Model_Validation","AllLog_Negbin_NoOutlier",
#                       paste("ModVal", i,abbreviate(colnames(df.rspX),3)[i], ".jpeg", sep = ""))
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
  scatter.smooth(log(predict(ls.glbModels[[i]])), E1, cex.lab = 1.5, xlab="log Predicted values", ylab="Pearson Residuals")
  abline(h = 0, v=0, lty=2); text(log(predict(ls.glbModels[[i]])), E1, labels = row.names(df.exp2), pos = 4)
  scatter.smooth(log(predict(ls.glbModels[[i]])), E2, cex.lab = 1.5, xlab="log Predicted values", ylab="Response Residuals")
  abline(h = 0, v=0, lty=2)
  scatter.smooth(log(predict(ls.glbModels[[i]])), E3, cex.lab = 1.5, xlab="log Predicted values", ylab="Deviance Residuals")
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
  
  title(paste(i,colnames(df.rspX)[i]), outer=TRUE)
  #dev.off()
}

# 2.3.2 Plots including observed values ####
#****************************************************************************************

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
  
  car::Boxplot(F1 ~ df.exp2$age_class)
  car::Boxplot(F1 ~ df.exp2$samcam)
  car::Boxplot(F1 ~ df.exp2$ACSC)
  
  title(colnames(df.rspX)[k], outer=TRUE)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# # Methods from influence.ME #####
#////////////////////////////////////////////////////////////////////////////////////////
# for(i in 1:p){
#   df.exp$y <- df.rspX[,i]
#   df.exp2 <- df.exp[outlier[[i]],]
#   alt.est <- influence(ls.glbModels[[i]], obs=TRUE)
#   print(tail(cooks.distance.estex(alt.est, sort=TRUE),3))
#   #dfbetas(alt.est)
# }
#////////////////////////////////////////////////////////////////////////////////////////


# 3. Get P-Values for Variables from Likelihood ratio Tests _________________________####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ls.bestModels1 <- ls.glbModels
ls.bestModels <- ls.bestModels1

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
df.LRT1 <- df.LRT
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 4. Model selection for significant p-values ______________________________________#####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 4.1 Remove non significant interaction #####
#****************************************************************************************

for(i in 1:p){
    df.exp$y <- df.rspX[,i]
    df.exp2 <- df.exp[outlier[[i]],]
    df.exp2$ID  <- 1:nrow(df.exp2)
    if (df.LRT[4,2+((i*2))] > 0.056){
      G2 <- update(ls.glbModels[[i]], .~.-age_class:samcam)
    #alt.est <- influence.ME::influence(G2, obs=TRUE)
   #print(colnames(df.rspX)[i])
    #print(tail(cooks.distance.estex(alt.est, sort=TRUE),3))
    name <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
    ls.glbModels[[i]] <- assign(name, G2)}
    }

ls.bestModels2 <- ls.glbModels
ls.bestModels <- ls.bestModels2

# 4.2 Repeat Likelihood Ratio tests #####
#****************************************************************************************

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
df.LRT2 <- df.LRT

# 4.3 Remove non-significant main effects (age_class and samcam only) ####
#****************************************************************************************

for(i in 1:p){
    df.exp$y <- df.rspX[,i]
    df.exp2 <- df.exp[outlier[[i]],]
    df.exp2$ID  <- 1:nrow(df.exp2)
    if (df.LRT[3,2+((i*2))] > 0.056 & df.LRT[2,2+((i*2))] < 0.056 &  is.na(df.LRT[4,2+((i*2))])==T){
    G3 <- update(ls.glbModels[[i]], . ~ . -samcam)
    name <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
    ls.glbModels[[i]] <- assign(name, G3)}
    if (df.LRT[2,2+((i*2))] > 0.056 & df.LRT[3,2+((i*2))] < 0.056 &  is.na(df.LRT[4,2+((i*2))])==T){
    G3 <- update(ls.glbModels[[i]], . ~ . -age_class)
    name <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
    ls.glbModels[[i]] <- assign(name, G3)}
    }

# 4.4 Repeat Lieklihood ratio Tests ####
#****************************************************************************************

ls.bestModels3 <- ls.glbModels
ls.bestModels <- ls.bestModels3

df.LRT3 <- matrix(NA,1+q,2+(2*p))
colnames(df.LRT3) <- c("Env", "DF", rep(colnames(df.rspX)[1:p], each=2))
df.LRT3[1,] <- c("X", "X", rep(c("CHI2", "p-value"),p))

for(i in 1:p){
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  df.exp2$ID <- 1:nrow(df.exp2)
  if(df.families[i,j] %in% c("poisson", "quasipoisson1")){
    obj.afex <- afex::mixed(ls.bestModels[[i]],  df.exp2, family="poisson", method="LRT")
    if(rownames(obj.afex[[1]])[1]=="samcam"){
      df.LRT3[3,2+(i*2)-1] <- round(obj.afex[[1]][1,2],2) # Chisq
      df.LRT3[3,2+((i*2))] <- round(obj.afex[[1]][1,4],3)} #p-value
    if(rownames(obj.afex[[1]])[1]=="age_class"){
    df.LRT3[2:4,2+(i*2)-1] <- round(obj.afex[[1]][1:3,2],2) # Chisq
    df.LRT3[2:4,2+((i*2))] <- round(obj.afex[[1]][1:3,4],3)}} #p-value
  if(df.families[i,j] %in% c("normal", "lognormal")) {
    obj.afex <- afex::mixed(ls.bestModels[[i]],  df.exp2,  method="LRT")
    if(rownames(obj.afex[[1]])[1]=="samcam"){
    df.LRT3[3,2+(i*2)-1] <- round(obj.afex[[1]][1,2],2) # Chisq
    df.LRT3[3,2+((i*2))] <- round(obj.afex[[1]][1,4],3)} #p-value
    if(rownames(obj.afex[[1]])[1]=="age_class"){
    df.LRT3[2:4,2+(i*2)-1] <- round(obj.afex[[1]][1:3,2],2) # Chisq
    df.LRT3[2:4,2+((i*2))] <- round(obj.afex[[1]][1:3,4],3)}} #p-value
  if(df.families[i,j] == "negbin"){
    wald1 <- car::Anova(ls.bestModels[[i]],  type="II", method="LR")
    #obj.afex <- afex::mixed(ls.bestModels[[i]], df.exp2, family=negative.binomial(theta=getME(ls.bestModels[[i]], "glmer.nb.theta")), control=con,method="LRT") # getME(ls.bestModels[[i]]); lme4:::getNBdisp(ls.bestModels[[i]]))
    if(rownames(obj.afex[[1]])[1]=="samcam"){
      df.LRT3[3,2+(i*2)-1] <- round(obj.afex[[1]][1,2],2) # Chisq
      df.LRT3[3,2+((i*2))] <- round(obj.afex[[1]][1,4],3)} #p-value
    if(rownames(obj.afex[[1]])[1]=="age_class"){
    df.LRT3[2:4,2+((i*2)-1)] <- round(wald1$"Chisq",2)[1:3]
    df.LRT3[2:4,2+(i*2)] <- round(wald1$"Pr(>Chisq)",3)[1:3]}}
  #   if(df.families[i,j] == "quasipoisson1"){
  #     wald1 <- car::Anova(ls.bestModels[[i]],  type="II", method="LR");wald1
  #     df.LRT3[2:4,2+((i*2)-1)] <- round(wald1$"Chisq",2)[1:3]
  #     df.LRT3[2:4,2+(i*2)] <- round(wald1$"Pr(>Chisq)",3)[1:3]}
}
df.LRT3[2:4,1] <- row.names(obj.afex[[1]][1])
df.LRT3[2:4,2] <- obj.afex[[1]]$"Chi Df"
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# 5. Post Hoc data inspection with lsmeans package __________________________________####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(lsmeans)
library(multcompView)

ls.bestModels <- ls.bestModels3
df.LRT <- df.LRT3

# 5.1 Multicomparisons ####
#****************************************************************************************

df.posthoc <- matrix(NA,8,2+(2*p))
colnames(df.posthoc) <- c("AC*SC", "AC+SC", rep(colnames(df.rspX),each=2))

for (i in 1:p) {
  if (df.LRT[3,2+((i*2))] < 0.056 & is.na(df.LRT[3,2+((i*2))])==F){
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ samcam)
    x <- cld(lsm, type = "response", sort=FALSE, Letters=c("abcdefg"), adjust="None")
    df.posthoc[5:8,2+((2*i)-1)] <- round(x[,2],2)
    df.posthoc[5:8,2+((2*i)-0)] <- x$".group"}
  if (df.LRT[2,2+((i*2))] < 0.056 & is.na(df.LRT[2,2+((i*2))])==F){
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ age_class)
    x <- cld(lsm, type = "response", sort=FALSE, Letters=c("abcdefg"), adjust="None")
    df.posthoc[1:4,2+((2*i)-1)] <- round(x[,2],2)
    df.posthoc[1:4,2+((2*i)-0)] <- x$".group"}
  if (df.LRT[4,2+((i*2))] < 0.056 & is.na(df.LRT[4,2+((i*2))])==F) {
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ age_class*samcam)
    x <- cld(lsm, type = "response", sort=FALSE, Letters=c("abcdefg"), adjust="None")
    df.posthoc[,2+((2*i)-1)] <- round(x[,3],2)
    df.posthoc[,2+((2*i)-0)] <- x$".group"}
  # to see the results graphically
  #p1 <- plot(lsm, by = "samcam", intervals = TRUE, type = "response")
  #print(p1)
  #title(names(ncr.biglmer)[i], outer=TRUE)
}

df.posthoc[,1] <- paste(paste(x[,1]),rep(c(2,4),each=4))
df.posthoc[,2] <- c(paste(x[1:4,1]),2,4,2,4)

fam <- c("XX", "XX",rep(df.families[1:35,j], each=2))
result <- rbind(df.LRT,df.posthoc,fam)

for (i in 1:p){ 
  if(Result_adjFam_deletedOlier_NoIntpsel[13,2+((2*i)-1)]=="lognormal")
    for(k in 5:12){
      Result_adjFam_deletedOlier_NoIntpsel[k,2+((2*i)-1)] = round(expm1(as.numeric(try1[k,2+((2*i)-1)])),2)#, na.omit==T)
    }}

View(result)

# write.csv(result,file="Analysis/Silphie_TimeSeries/Results/Result_adjFam_deletedOlier_NoIntpsel.csv")
# save(result,file="Analysis/Silphie_TimeSeries/Results/Result_adjFam_deletedOlier_NoIntpsel.rda")
# Result_adjFam_deletedOlier_NoIntpsel <- result

# 5.2 Selected contrasts ####
#****************************************************************************************

# A. Age Class ####
ls.lsmAC <- list()
df.posthocAC <- matrix(NA,12,2+(2*p))
colnames(df.posthocAC) <- c("contrast", "samcam", rep(colnames(df.rspX),each=2))

for (i in 1:p) {
  if (df.LRT[4,2+((i*2))] < 0.056 & is.na(df.LRT[2,2+((i*2))])==F) 
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ age_class|samcam, at=list(samcam=c("2","4")))
  if (df.LRT[2,2+((i*2))] < 0.056 & is.na(df.LRT[2,2+((i*2))])==F) 
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ age_class)
  
  x1 <- contrast(lsm,"pairwise")
  xx1 <- summary(x1, type="response",adjust = "none");xx1 # trusted p-values? rate.ratio for poisson
  df.posthocAC[,2+((2*i)-0)] <- round(xx1$"p.value",3)
  
  if(df.families[i,j] %in% c("poisson", "quasipoisson1", "negbin")) {
    x2 <- summary(pairs(regrid(lsm)));x2 # highest p-values; differences in counts
    df.posthocAC[,2+((2*i)-1)] <- round(x2$estimate,2)
    #xx2 <- summary(regrid(x1), type="response");xx2 # lowest p-values
    #df.posthocAC[,2+((2*i)-0)] <- round(x2$"p.value",3)
  }
  else
    df.posthocAC[,2+((2*i)-1)] <- round(xx1$"estimate",2)
  name <- paste("lsm",i,names(df.rspX)[i], sep = ".")
  ls.lsmAC[[i]] <- assign(name, lsm)
} 

df.posthocAC[,1] <- paste(xx1$"contrast")
df.posthocAC[,2] <- paste(xx1$"samcam")

# B. Sampling Campaign ####
ls.lsmSC <- list()
df.posthocSC <- matrix(NA,4,2+(2*p))
colnames(df.posthocSC) <- c("contrast", "age_class", rep(colnames(df.rspX),each=2))

for (i in 1:p) {
  if (df.LRT[4,2+((i*2))] < 0.056 & is.na(df.LRT[4,2+((i*2))])==F) 
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ samcam|age_class, at=list(samcam=c("2","4")))
  if (df.LRT[3,2+((i*2))] < 0.056 & is.na(df.LRT[3,2+((i*2))])==F) 
    lsm <- lsmeans::lsmeans(ls.bestModels[[i]],  ~ samcam)
  
  x1 <- contrast(lsm,"pairwise")
  xx1 <- summary(x1, type="response",adjust = "none");xx1 # trusted p-values? rate.ratio for poisson
  df.posthocSC[,2+((2*i)-0)] <- round(xx1$"p.value",3)
  
  if(df.families[i,j] %in% c("poisson", "quasipoisson1", "negbin")) {
    x2 <- summary(pairs(regrid(lsm)));x2 # highest p-values; differences in counts
    df.posthocSC[,2+((2*i)-1)] <- round(x2$estimate,2)
    #xx2 <- summary(regrid(x1), type="response");xx2 # lowest p-values
    #df.posthocSC[,2+((2*i)-0)] <- round(x2$"p.value",3)
    }
  else
   df.posthocSC[,2+((2*i)-1)] <- round(xx1$"estimate",2)
  name <- paste("lsm",i,names(df.rspX)[i], sep = ".")
  ls.lsmSC[[i]] <- assign(name, lsm)
} 


df.posthocSC[,1] <- paste(xx1$"contrast")
df.posthocSC[,2] <- paste(xx1$"age_class")

fam <- c("XX", "XX",rep(df.families[1:35,j], each=2))
result <- rbind(df.LRT,df.posthocAC,df.posthocSC,fam)
View(result)
# write.csv(result,file="Analysis/Silphie_TimeSeries/Results/Result_AllLog_Negbin_NoOutlier.csv")
# save(result,file="Analysis/Silphie_TimeSeries/Results/Result_AllLog_Negbin_NoOutlier.rda")
# Result_AllLog_Negbin_NoOutlier <- result
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
windows(record=TRUE)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  car::Boxplot(df.exp2$y ~ df.exp2$age_class)
  car::Boxplot(df.exp2$y ~ df.exp2$samcam)
  car::Boxplot(df.exp2$y ~ df.exp2$ACSC)
  title(names(df.rspX)[i],outer=TRUE)
}
windows(record=TRUE)
for(i in 1:p) {
  par(mfrow = c(2,2),
      mar = c(3,3,0,1),
      mgp = c(1.5,0.5,0),
      tck = -0.03,
      oma = c(0,0,2,0))
  df.exp2 <- df.exp[outlier[[i]],]
  df.exp2$y <- predict(ls.glbModels[[i]], type="response")
  df.exp2$y <- fitted(ls.glbModels[[i]], type="response")
  car::Boxplot(df.exp2$y ~ df.exp2$age_class)
  car::Boxplot(df.exp2$y ~ df.exp2$samcam)
  car::Boxplot(df.exp2$y ~ df.exp2$ACSC)
  title(names(df.rspX)[i],outer=TRUE)
}

#df.posthoc.allInt <- df.posthoc
#save(list=c("ls.bestModels","ls.lsm", "ls.lsmAC", "ls.lsmSC",  "df.FpvalueR2", "df.posthoc", "df.posthocAC", "df.posthocSC"), file="Results/ANOVATables/All_LMM.rda")
# write.csv(df.posthoc, file="Results/ANOVATables/PostHoc_All_LMM.csv")
# write.csv(df.posthocAC, file="Results/ANOVATables/PostHocAC_All_LMM.csv")
# write.csv(df.posthocSC, file="Results/ANOVATables/PostHocSC_All_LMM.csv")
#  write.table(df.posthoc, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
#  write.table(df.posthocAC, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
#  write.table(df.posthocSC, file="Results/ANOVATables/FpR2_All_LMM.csv", append=TRUE, sep=",", dec=".", qmethod="double", col.names=NA)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 6. Boxplots of original and  fitted values according to the best Models____________####

windows(record=T)

for (i in 1:p) {
  # Save Plots extern  
#   mypath <- file.path("Analysis","Silphie_TimeSeries","results","Boxplots_adjFam_delOlier_NoIntpsel",
#                       paste("ModVal", i,abbreviate(colnames(df.rspX),3)[i], ".pdf", sep = ""))
#   
#   pdf(file=mypath)
  par(mfrow=c(2,3),
      oma=c(0,0,3,0))
  
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  F1 <- fitted(ls.bestModels[[i]], type="response")
  P1 <- predict(ls.bestModels[[i]], type="response")
  car::Boxplot(F1 ~ df.exp2$age_class)
  car::Boxplot(F1 ~ df.exp2$samcam)
  car::Boxplot(F1 ~ df.exp2$ACSC)
  car::Boxplot(df.exp2$y ~ df.exp2$age_class)
  car::Boxplot(df.exp2$y ~ df.exp2$samcam)
  car::Boxplot(df.exp2$y ~ df.exp2$ACSC)
  title(paste(i,colnames(df.rspX)[i]), outer=TRUE)
  #dev.off()
}


# Depot _____________________________________________________________________________####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
windows(record=T)
for(i in 1:p){
  par(mfrow=c(2,2))
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  df.exp2$ID <- 1:nrow(df.exp2)
  a <- paste(formula(ls.bestModels[[i]])[3])
  if("age_class:samcam" %in% unlist(strsplit(a, " "))){
    M1 <- update(ls.bestModels[[i]], . ~ . + ACSC - age_class - samcam - age_class:samcam)
    posthoc <- glht(M1, linfct=mcp(ACSC="Tukey"))
    letters <- cld(posthoc)
    print(plot(letters, main=colnames(df.rspX)[i]))}
  else {
    if("age_class" %in% unlist(strsplit(a, " "))){
      M1 <- ls.bestModels[[i]]
      posthoc <- glht(M1, linfct=mcp(age_class="Tukey"))
      letters <- cld(posthoc)
      print(plot(letters, main=colnames(df.rspX)[i]))
    }
    if("samcam" %in% unlist(strsplit(a, " "))){
      M1 <- ls.bestModels[[i]]
      posthoc <- glht(M1, linfct=mcp(samcam="Tukey"))
      letters <- cld(posthoc)
      print(plot(letters, main=colnames(df.rspX)[i]))
  }
}}

windows(record=T)

for (i in 1:p) {
 # Save Plots extern  
    mypath <- file.path("Analysis","Silphie_TimeSeries","results","Boxplots_adjFam_delOlier_NoIntpsel",
                        paste("ModVal", i,abbreviate(colnames(df.rspX),3)[i], ".pdf", sep = ""))
    
    pdf(file=mypath)
    par(mfrow=c(2,3),
        oma=c(0,0,3,0))
    
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  F1 <- fitted(ls.bestModels[[i]], type="response")
  P1 <- predict(ls.bestModels[[i]], type="response")
  car::Boxplot(F1 ~ df.exp2$age_class)
  car::Boxplot(F1 ~ df.exp2$samcam)
  car::Boxplot(F1 ~ df.exp2$ACSC)
  car::Boxplot(df.exp2$y ~ df.exp2$age_class)
  car::Boxplot(df.exp2$y ~ df.exp2$samcam)
  car::Boxplot(df.exp2$y ~ df.exp2$ACSC)
  title(paste(i,colnames(df.rspX)[i]), outer=TRUE)
  dev.off()
}

# Ease manual code writing for Models
# v1 <- c()
# for (i in 1:ncol(df.rspX)) {
#   v1[i] <- c(paste("glbM",i,abbreviate(colnames(df.rspX),3)[i], " <- ", sep="."))
# }
# v1

# 2.2.3 [Special Case, no continuous covariates] Likelihood ratio tests by Hand 
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


mypath <- file.path("Analysis","Silphie_TimeSeries","Model_Validation","PostHoc")
save(list=c("df.posthocAC.P", "df.posthocAC.QP", "df.posthocAC.all_N", "df.posthocAC.all_P_logN", "df.posthocAC.all_QP_logN", "df.posthocAC.all_NB_logN"), file="Analysis/Silphie_TimeSeries/Model_Validation/PostHoc/posthoctables.rda")
load(file="Analysis/Silphie_TimeSeries/Model_Validation/PostHoc/posthoctables.rda")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 6. Multimodel Inference # demo(dredge.subset)  ____________________________________####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 6.1 Generate Cluster ####
#****************************************************************************************
library(parallel)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))
clusterEvalQ(clust, c(library("lme4")))
#clusterExport(clust, varlist=c("dt.exp", "con"))

# 6.2 Subsets of models excluded from dredge (optional): ####
#****************************************************************************************
opo <- env[,c("mc","pH","cnratio","ats1")]
opo <- as.data.frame(opo)

# 6.3 Dredge ####
#****************************************************************************************
library(MuMIn)
options(na.action = na.fail)
p <- ncol(df.rspX)
ls.dredge <- list()

for(i in 1:p) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  #   dt.exp$y <- dt.rsp.abn[,i, with=F]
  clusterExport(clust, varlist=c("df.exp2", "con"))
  GM.dredge <- pdredge(ls.glbModels[[i]], m.min = 1, cluster=clust)
  #fixed=c("age_class", "samcam"),
  name <- c(paste("Dredge",i,abbreviate(colnames(df.rspX),3)[i], sep="."))
  assign(name, GM.dredge)
  ls.dredge[[i]] <- assign(name, GM.dredge)
}

# 6.4 Get the best models ####
#****************************************************************************************
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
  print(i)
  print(colnames(df.rspX)[i])
  print(formula(ls.bestModels[[i]]))
}

# 6.5 Sets of Candidate Models
#****************************************************************************************
df.compM2 <- vector()
for(i in 1:p){
  delta4 <- subset(ls.dredge[[i]], delta < 4)
  df.compM <- data.frame(delta4)
  df.compM[,ncol(df.compM)+1] <- rep(colnames(df.rspX)[i], length(rownames(df.compM)))
  df.compM2 <- rbind(df.compM2, df.compM)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 7. Other Model Selection Technics  ________________________________________________####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# 7.1 Leave One Out ####
#****************************************************************************************
for ( i in 1:p) {
  df.exp$y <- df.rspX[,i]
  df.exp2 <- df.exp[outlier[[i]],]
  print(drop1(ls.glbModels[[i]]))
}

# 7.2 Stepwise Selection ####
#****************************************************************************************
step(ls.glbModels[[1]])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







