###########################
# F2 Nematodes
# GGplots for capscale
# Quentin Schorpp
# 19.08.2015
###########################



sc1 <- scores(dbRDA, display=c("sp", "lc","cn", "bp"), scaling=1)
sc2 <- scores(dbRDA, display=c("sp", "lc","cn", "bp"), scaling=2)

# age_class centroids
ac.cid1 <- data.frame(sc1$centroids)
ac.cid1$txt <- c("C.mays", "S.p. - young", "S.p. - int1", "S.p. - int2", "S.p. - old")
ac.cid2 <- data.frame(sc2$centroids)
ac.cid2$txt <- c("C.mays", "S.p. - young", "S.p. - int1", "S.p. - int2", "S.p. - old")


# species
spec.cid1 <- data.frame(sc1$species)
spec.cid1$txt <- rownames(spec.cid1)
spec.cid2 <- data.frame(sc2$species)
spec.cid2$txt <- rownames(spec.cid2) 
# you can also put label = rownames() etc in the ggplot code. 
# But sometimes you may want it in a column for later

# arrows for continous environmental
env.ar1 <- data.frame(sc1$biplot)
env.ar1$txt <- rownames(env.ar1)
env.ar2 <- data.frame(sc2$biplot)
env.ar2$txt <- rownames(env.ar2)

mult <- attributes(scores(dbRDA))$const

#sites
site_scores1 <- data.frame(scores(dbRDA, scaling=1)$sites)
site_scores1$age_class <- dbenv$age_class

site_scores2 <- data.frame(scores(dbRDA, scaling=2)$sites)
site_scores2$age_class <- dbenv$age_class

# site constraints 
site_constraints1 <- data.frame(sc1$constraints)
site_constraints1$age_class <- dbenv$age_class
site_constraints2 <- data.frame(sc2$constraints)
site_constraints2$age_class <- dbenv$age_class


#sc.mean=aggregate(site_scores[,1:2],list(group=site_scores$age_class),mean)

df_ell1.1 <- data.frame()
for(g in levels(site_scores1$age_class)){
  df_ell1.1 <- rbind(df_ell1.1, cbind(as.data.frame(with(site_scores1[site_scores1$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}

df_ell1.2 <- data.frame()
for(g in levels(site_constraints1$age_class)){
  df_ell1.2 <- rbind(df_ell1.2, cbind(as.data.frame(with(site_constraints1[site_constraints1$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}


df_ell2.1 <- data.frame()
for(g in levels(site_scores2$age_class)){
  df_ell2.1 <- rbind(df_ell2.1, cbind(as.data.frame(with(site_scores2[site_scores2$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}

df_ell2.2 <- data.frame()
for(g in levels(site_constraints2$age_class)){
  df_ell2.2 <- rbind(df_ell2.2, cbind(as.data.frame(with(site_constraints2[site_constraints2$age_class==g,],
                                                         veganCovEllipse(cov.wt(cbind(CAP1,CAP2),wt=rep(1/length(CAP1),length(CAP1)))$cov,center=c(mean(CAP1),mean(CAP2)))))
                                      ,age_class=g))
}


# Scaling 1
Scaling1 <- ggplot(ac.cid1, aes(x = CAP1, y = CAP2))+
  geom_text(aes(label = txt),fontface="italic", cex=3)+ # label tells geom_text which text you want to plot
  geom_point(pch=17, size=5, col="red")+ # label tells geom_text which text you want to plot
  geom_text(data = spec.cid1, aes(label = txt), colour = "orange") +
  coord_cartesian(y=c(1.2*min(spec.cid1$CAP2),1.2*max(spec.cid1$CAP2)), x=c(1.2*min(spec.cid1$CAP1),1.3*max(spec.cid1$CAP1))) + 
  #NB note that this is a convenience wrapper and may cut data out of your plot
  #important if you are calculating stats in the plot - these will be excluded
  geom_segment(data = env.ar1[-c(1:4),] , aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "light blue") + #grid is required for arrow to work.
  geom_text(data = env.ar1[-c(1:4),] , aes(x= (mult + mult/20) * CAP1, y = (mult + mult/20) * CAP2,  label = env.ar1[-c(1:4),]$txt), size = 5, hjust = 0.5) +
  # we add 10% to the text to push it slightly out from arrows
  # otherwise you could use hjust and vjust. I prefer this option
  #geom_point(data=site_scores1, aes(x = CAP1, y = CAP2, shape=dbenv$age_class)) +
  geom_point(data=site_constraints1, aes(x = CAP1, y = CAP2, shape=dbenv$age_class)) +
  #geom_path(data=df_ell1.1, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  geom_path(data=df_ell1.2, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  #geom_text(data = site_scores, aes(x = CAP1, y = CAP2, label = rownames(site_scores)), size = 3, colour = "grey50")+    
  theme_bw()

# use env.ar1[-c(1:4),]  to ommit factor levels



#Scaling 2
Scaling2 <- ggplot(ac.cid2, aes(x = CAP1, y = CAP2))+
  geom_text(aes(label = txt),fontface="italic", cex=3)+ # label tells geom_text which text you want to plot
  geom_point(pch=17, size=5, col="red")+ # label tells geom_text which text you want to plot
  geom_text(data = spec.cid2, aes(label = txt), colour = "orange") +
  coord_cartesian(y=c(1.2*mult*min(env.ar2$CAP2),1.2*max(site_scores2$CAP2)), x=c(1.2*mult*min(env.ar2$CAP1),1.3*max(site_scores2$CAP1))) + 
  #NB note that this is a convenience wrapper and may cut data out of your plot
  #important if you are calculating stats in the plot - these will be excluded
  geom_segment(data = env.ar2[-c(1:4),], aes(x = 0, xend = mult * CAP1, y = 0, yend = mult * CAP2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "light blue") + #grid is required for arrow to work.
  geom_text(data = env.ar2[-c(1:4),], aes(x= (mult + mult/20) * CAP1, y = (mult + mult/20) * CAP2,  label = env.ar2[-c(1:4),]$txt), size = 5, hjust = 0.5) +
  # we add 10% to the text to push it slightly out from arrows
  # otherwise you could use hjust and vjust. I prefer this option
  #geom_point(data=site_scores2, aes(x = CAP1, y = CAP2, shape=dbenv$age_class)) +
  geom_point(data=site_constraints2, aes(x = CAP1, y = CAP2, shape=dbenv$age_class)) +
  #geom_path(data=df_ell2.1, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  geom_path(data=df_ell2.2, aes(x=CAP1, y=CAP2, colour=age_class), size=0.5, linetype=2) +
  #geom_text(data = site_scores, aes(x = CAP1, y = CAP2, label = rownames(site_scores)), size = 3, colour = "grey50")+    
  theme_bw()