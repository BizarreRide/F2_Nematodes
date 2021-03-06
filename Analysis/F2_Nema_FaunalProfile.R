###########################
# F2 Nematodes
# Faunal Profile
# Quentin Schorpp
# 07.10.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("Analysis/Sources/RequiredPackages.R")
source("Data/DataProcessing/DataProcessing.R") 
env1 <- droplevels(env.org[16:45,])
source("Data/DataProcessing/EnvDataProcessing.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Calculate faunal profile Indices ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data <- fam.org
source("Data/DataProcessing/FamDatProcessing.R") 

# Faunal Profile on subsample abundances
data <- round(fam.rel*100,0) # choose basis data for faunal profile (.org, .rel, .usc)
#data1 <- aggregate(data,list(env1$age_class,env1$samcam), mean)
#data <- data1[,-c(1,2)]
#str(data)
source("Data/DataProcessing/FaunalProfileIndices.R") 
FaPro <- data.frame(FaPro)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Publication plots

# Altogether
p1 <- ggplot(FaPro, aes(x=SI, y=EI, shape=env1$age_class)) + 
  geom_point(size=2, fill="black") + #aes(size=env1$n)
  scale_shape_manual(values=c(1,22,21,24,2), labels=c("Cm", "Sp_Y", "Sp_I1","Sp_I2","Sp_O"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Faunal Profile") +
  mytheme + theme(legend.position = "bottom") 
p1

#ggsave("FaunalProfileTotal.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)


# Only Maize
p2 <- ggplot(FaPro[env1$crop=="Maize",], aes(x=SI, y=EI, fill=env1$samcam[env1$crop=="Maize"], shape=env1$age_class[env1$crop=="Maize"])) + 
  geom_point(size=2) + #aes(size=env1$n)
  scale_shape_manual(values=c(21), labels=c("Cm"), guide=guide_legend(title=NULL)) + 
  scale_fill_manual(values=c("black","white"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Maize 2012/13") +
  mytheme + theme(legend.position = "none") 
p2

#ggsave("FaunalProfileTotal.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)

FaPro2 <- droplevels(FaPro[!env1$crop=="Maize",])
env2 <- env1[!env1$crop=="Maize",]

# Autumn 2012
p3 <- ggplot(FaPro2[env2$samcam==2,], aes(x=SI, y=EI, shape=env2$age_class[env2$samcam==2])) + 
  geom_point(size=2, fill="black") + #aes(size=env1$n)
  scale_shape_manual(values=c(22,23,25,24), labels=c("Sp_Y", "Sp_I1","Sp_I2","Sp_O"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Cup plant 2012") +
  mytheme + theme(legend.position = "none") 
p3

#ggsave("Results/FaunalProfileSamCam2.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)


# Autumn 2013
p4 <- ggplot(FaPro2[env2$samcam==4,], aes(x=SI, y=EI, shape=env2$age_class[env2$samcam==4])) + 
  geom_point(size=2) + #aes(size=env1$n)
  scale_shape_manual(values=c(22,23,25,24), labels=c("Sp_Y", "Sp_I1","Sp_I2","Sp_O"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Cup plant 2013") +
  mytheme + theme(legend.position = "none") 
p4

#ggsave("Results/FaunalProfileSamCam4.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)


gridExtra::grid.arrange(p3,p4,p2, ncol=3)
dev.copy2pdf(file="Results/FaunalProfile3parts.pdf", width=19/2.54, height=7/2.54, useDingbats=F, out.type = "pdf")# 


FaPro.av <- aggregate(FaPro,list(age_class=env1$age_class,samcam=env1$samcam), mean)

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}
FaPro.se <- aggregate(FaPro,list(age_class=env1$age_class,samcam=env1$samcam), st.err)

# Average and standard error of the mean
ggplot(FaPro.av, aes(x=SI, y=EI, shape=interaction(age_class,samcam))) + 
  geom_point(size=2, fill="black") + #aes(size=env1$n)
  geom_errorbar(aes(ymin=EI-FaPro.se$EI,ymax=EI+FaPro.se$EI), height=1.2,size=0.2) +
  geom_errorbarh(aes(xmin=SI-FaPro.se$SI,xmax=SI+FaPro.se$SI),height=1.2, size=0.2) +
  scale_shape_manual(values=c(1,0,5,6,2,21,22,23,25,24), labels=c("Cm-2012", "Sp_Y-2012", "Sp_I1-2012","Sp_I2-2012","Sp_O-2012","Cm-2013", "Sp_Y-2013", "Sp_I1-2013","Sp_I2-2013","Sp_O-2013"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Faunal Profile") +
  mytheme + theme(legend.position = "bottom") 

ggsave("FaunalProfileTotalAvSe.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)



# Miscallaneous Plot
ggplot(FaPro, aes(x=SI, y=EI, shape=env1$age_class, fill=env1$age, group=env1$field.ID)) + 
  geom_point(size=4) + #aes(size=env1$n)
  geom_path(aes(group=env1$field.ID),arrow = arrow(angle = 15, type = "closed")) +
  scale_shape_manual(values=c(1,22,21,24,23)) + 
  scale_fill_gradient(low="blue", high="red") +
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50) + geom_hline(y=50) + 
  ggtitle("Faunal Profile") +
  #facet_grid(env1$samcam ~ .) +
  theme_bw()
ggsave("Results/FaunalProfileArrows.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)

ggplot(FaPro[env1$samcam==2,], aes(x=SI, y=EI, shape=env1$age_class[env1$samcam==2], fill=env1$age[env1$samcam==2])) + 
  geom_point(size=4) + #aes(size=env1$n)
  scale_shape_manual(values=c(1,22,21,24,23)) + 
  scale_fill_gradient(low="blue", high="red") +
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50) + geom_hline(y=50) + 
  ggtitle("Faunal Profile") +
  #facet_grid(env1$samcam ~ .) +
  theme_bw()

ggplot(FaPro[env1$samcam==4,], aes(x=SI, y=EI, shape=env1$age_class[env1$samcam==4], fill=env1$age[env1$samcam==4])) + 
  geom_point(size=4) + #aes(size=env1$n)
  scale_shape_manual(values=c(1,22,21,24,23)) + 
  scale_fill_gradient(low="blue", high="red") +
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50) + geom_hline(y=50) + 
  ggtitle("Faunal Profile") +
  #facet_grid(env1$samcam ~ .) +
  theme_bw()

# Quadrat A: poorly structured, n-enrichment
# Quadrat B: maturing, structured food-web and N-enrichment
# Quadrat C: Undisturbed strucured food webs, relatively low primary production
# Quadrat D: basal, poorly developed food web condition


# Which other covariates influence EI and SI ? 
par(mar=c(1,1,1,1),
    mfrow=c(5,3))
for(i in 1:15) {
  scatter.smooth(FaPro$EI ~ env.fin[,i])
}

par(mar=c(1,1,1,1),
    mfrow=c(5,3))
for(i in 1:15) {
  scatter.smooth(FaPro$SI ~ env.fin[,i])
}



