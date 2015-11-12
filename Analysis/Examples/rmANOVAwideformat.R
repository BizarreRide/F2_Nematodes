#§§§§§§§§§§§§§§§§§§§§§§§§§§
# F2 Nematodes
# Univaraite multiple repeated measurements Tests
# Quentin Schorpp
# 11.08.2015
#§§§§§§§§§§§§§§§§§§§§§§§§§§



# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# as needed
# source("Data/DataProcessing/MaturityIndices.R") 
# source("Data/DataProcessing/FaunalProfileIndices.R") 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



# Data processing ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create interaction factor for interaction between age_class and samcam, 
# to see if age_class effects stay constant

env.fin$agsam <- with(env.fin, interaction(age_class,samcam))
#env.fin$c <-  env1$c # include carbon

#fam.rel <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)

indices <- droplevels(cbind(groups,env.fin, location=env1$location, N=rowSums(fam.org)))
indices$ID <- 1:nrow(indices)
indices.backup <- indices
indices <- indices.backup
indices$nsamcam <- as.numeric(factor(indices$samcam))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# 1. Analysis of FeedingTypes ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Data/DataProcessing/FeedingTypes.R") 
colnames(fety)[2:10] <- c("herbivoresb","herbivoresc","herbivoresd","herbivorese","fungivores","bacterivores", "carnivores", "omnivores", "Tylenchidae")
fety <- as.data.frame(fety[,-1])

fety$herbivores <- rowSums(fety[,c(1,2,3,4)])
fety$herbivores2 <- fety$herbivores + fety$Tylenchidae
fety$fungivores2 <- fety$fungivores + fety$Tylenchidae


fety.backup <- fety/indices$N
fety.backup <- cbind(fety.backup, group=indices$age_class, time=indices$samcam, subject=indices$field.ID, crop = indices$crop)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fety3 <- fety.backup[!fety.backup$group %in% "A_Cm",-c(1:12,16)]
fety3$response <- fety.backup$herbivores2[!fety.backup$group %in% "A_Cm"]


dfTemp2   <- reshape(fety3, v.names="response", timevar="time",
                    idvar=c("subject", "group"), direction="wide")
dfRBFpqW2 <- reshape(dfTemp2, v.names=c("response.2", "response.4"),
                    timevar="group", idvar="subject", direction="wide")
#dfRBFpqW3 <- reshape(dfTemp2, v.names=c("response.2", "response.4"),
                    timevar="subject", idvar="group", direction="wide")

fitRBFpq2   <- lm(cbind(response.2.E_Sp_old, response.4.E_Sp_old, response.2.D_Sp_int2, response.4.D_Sp_int2, response.2.C_Sp_int1, response.4.C_Sp_int1, response.2.B_Sp_young, response.4.B_Sp_young) ~ 1, data=dfRBFpqW2)

inRBFpq2    <- fety3[,c("group", "time")]

AnovaRBFpq2 <- Anova(fitRBFpq2, idata=inRBFpq, idesign=~group*time)
summary(AnovaRBFpq, multivariate=FALSE, univariate=TRUE)

colnames(dfRBFpqW2)
cind(response.2.E_Sp_old, response.4.E_Sp_old, response.2.D_Sp_int2, response.4.D_Sp_int2, response.2.C_Sp_int1, response.4.C_Sp_int1, response.2.B_Sp_young, response.4.B_Sp_young)

set.seed(123)
N    <- 10
P    <- 2
Q    <- 3
muJK <- c(rep(c(1, -2), N), rep(c(2, 0), N), rep(c(3, 3), N))
dfRBFpqL <- data.frame(id =factor(rep(1:N, times=P*Q)),
                       IV1=factor(rep(rep(1:P, each=N), times=Q)),
                       IV2=factor(rep(rep(1:Q, each=N*P))),
                       DV =rnorm(N*P*Q, muJK, 2))

dfTemp   <- reshape(dfRBFpqL, v.names="DV", timevar="IV1",
                    idvar=c("id", "IV2"), direction="wide")
dfRBFpqW <- reshape(dfTemp, v.names=c("DV.1", "DV.2"),
                    timevar="IV2", idvar="id", direction="wide")

fitRBFpq   <- lm(cbind(DV.1.1, DV.2.1, DV.1.2, DV.2.2, DV.1.3, DV.2.3) ~ 1,
                 data=dfRBFpqW)

inRBFpq    <- expand.grid(IV1=gl(P, 1), IV2=gl(Q, 1))

AnovaRBFpq <- Anova(fitRBFpq, idata=inRBFpq, idesign=~IV1*IV2)
summary(AnovaRBFpq, multivariate=FALSE, univariate=TRUE)


set.seed(123)
N      <- 10
P      <- 4
muJ    <- rep(c(-1, 0, 1, 2), each=N)
dfRBpL <- data.frame(id=factor(rep(1:N, times=P)),
                     IV=factor(rep(1:P,  each=N)),
                     DV=rnorm(N*P, muJ, 3))



