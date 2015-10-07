###########################
# F2 Nematodes
# Faunal Profile
# Quentin Schorpp
# 07.10.2015
###########################

# Load Data ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
source("data/RMAkeLikeFile.R")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




ggplot(indices, aes(x=SI, y=EI, shape=env1$age_class, fill=env1$age, group=env1$field.ID)) + 
  geom_point(aes(size=env1$n)) + 
  geom_line(arrow = arrow(angle = 15, ends = "first", type = "closed")) +
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

ggplot(indices, aes(x=SI, y=EI)) + 
  geom_point() + 
  #geom_line(arrow = arrow(angle = 15, ends = "first", type = "closed")) +
  #scale_shape_manual(values=c(1,22,21,24,23)) + 
  #scale_fill_gradient(low="blue", high="red") +
  scale_x_continuous(limits=c(0,100), name="Structure Index") +
  scale_y_continuous(limits=c(0,100), name="Enrichment Index") +
  #geom_vline(x=50) + geom_hline(y=50) + 
  ggtitle("Faunal Profile") +
  facet_grid(env1$samcam ~ .) +
  theme_bw()