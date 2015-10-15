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
data <- fam.rel # choose basis data for faunal profile (.org, .rel, .usc)
source("Data/DataProcessing/FaunalProfileIndices.R") 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Publication plot
ggplot(FaPro, aes(x=SI, y=EI, shape=env1$age_class)) + 
  geom_point(size=2, fill="black") + #aes(size=env1$n)
  scale_shape_manual(values=c(1,22,21,24,2), labels=c("Cm", "Sp_Y", "Sp_I1","Sp_I2","Sp_O"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Faunal Profile") +
  mytheme + theme(legend.position = "bottom") 

ggsave("FaunalProfileTotal.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)


ggplot(FaPro[env1$samcam==2,], aes(x=SI, y=EI, shape=env1$age_class[env1$samcam==2])) + 
  geom_point(size=2, fill="black") + #aes(size=env1$n)
  scale_shape_manual(values=c(1,22,21,24,2), labels=c("Cm", "Sp_Y", "Sp_I1","Sp_I2","Sp_O"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Faunal Profile") +
  mytheme + theme(legend.position = "bottom") 

ggsave("FaunalProfileSamCam2.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)


ggplot(FaPro[env1$samcam==4,], aes(x=SI, y=EI, shape=env1$age_class[env1$samcam==4])) + 
  geom_point(size=2, fill="black") + #aes(size=env1$n)
  scale_shape_manual(values=c(1,22,21,24,2), labels=c("Cm", "Sp_Y", "Sp_I1","Sp_I2","Sp_O"), guide=guide_legend(title=NULL)) + 
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50, size=0.2) + geom_hline(y=50, size=0.2) + 
  ggtitle("Faunal Profile") +
  mytheme + theme(legend.position = "bottom") 

ggsave("FaunalProfileSamCam4.pdf", width = 7, height= 8.6, units="cm", useDingbats=FALSE)



# Miscallaneous Plot
ggplot(FaPro, aes(x=SI, y=EI, shape=env1$age_class, fill=env1$age, group=env1$field.ID)) + 
  geom_point(size=4) + #aes(size=env1$n)
  geom_line(arrow = arrow(angle = 15, ends = "first", type = "closed")) +
  scale_shape_manual(values=c(1,22,21,24,23)) + 
  scale_fill_gradient(low="blue", high="red") +
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  geom_vline(x=50) + geom_hline(y=50) + 
  ggtitle("Faunal Profile") +
  #facet_grid(env1$samcam ~ .) +
  theme_bw()

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



