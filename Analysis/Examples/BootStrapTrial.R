




my.fit = fety.biglmer[[1]]


my.bootstrap.predictions.f <- function(data, indices){
  return(mean(predict(my.fit, newdata = data[indices, ], type = "response", allow.new.levels=TRUE), na.rm=TRUE))
}

## predict for year 1996 to year 2014
new.df <- indices[sample(nrow(indices), replace = TRUE), ] # totally chaotic data frame

agsam <- unique(indices$agsam)

my.results <- matrix(nrow=length(agsam), ncol = 4) # 4 column for index, mean proportion, upper ci and lower ci
                                                         # as many rows as predictions should be made

for(x in 1:length(agsam)){
  my.results[x, 1] <- agsam[x]
  
  new.df$agsam <- agsam[x]
  
  ## bootstrap using a realistic number of samples per year, say 20000
  my.boot.obj <- boot(data = new.df[sample(nrow(new.df), 20000, replace = TRUE), ], # repeated use of sample
                      statistic = my.bootstrap.predictions.f, 
                      R = 100)
  my.results[x, 2] <- my.boot.obj[[1]]
  my.results[x, 3:4] <- quantile(my.boot.obj[[2]], c(0.025, 0.975))
}
colnames(my.results) <- c("Year", "mean proportion", "lower.ci", "upper.ci")





DF.interaction <- Anova(fety.biglmer[[1]])$Df[3]

par.agsam <- fixef(fety.biglmer[[1]])[6:8]

vc.agsam <- vcov(fety.biglmer[[1]])[6:8, 6:8]

library(multcomp)
csimint(par.agsam, df=DF.interaction, covm=vc.agsam)


indices$hb.suc <- fety[!indices.backup$age_class %in% "A_Cm",11] 

indices$hb.fail <- indices$N-indices$herbivores.suc
str(indices)

model <- glmer(cbind(hb.suc, hb.fail) ~ age_class*samcam  + (1|ID) + (1|field.ID), family=binomial(link="logit"), indices, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

#confint(model)

testdata = expand.grid(age_class=unique(indices$age_class),
                       samcam = unique(indices$samcam))


X <- model.matrix(~ age_class*samcam, data = testdata)
testdata$fit <- X %*% fixef(model)
testdata$SE <- sqrt(  diag(X %*%vcov(fety.biglmer[[1]]) %*% t(X))  )
testdata$upr=testdata$fit+1.96*testdata$SE
testdata$lwr=testdata$fit-1.96*testdata$SE
endad.pred <- testdata


predfig.endad1 <- ggplot(endad.pred, aes(x = age_class, y = exp(fit), ymin = exp(lwr), ymax = exp(upr))) + 
  geom_bar(stat="identity",position = position_dodge(1), col="454545", size=0.15, fill="grey") +
  geom_errorbar(position = position_dodge(1),col="black",width=0.15, size=0.15) + 
  facet_grid(.~samcam) +
  geom_hline(xintercept = 1, size=0.15) +
  ylab("Anecic Abundance Ind./m²") +
  xlab("Age Class") +
  scale_x_discrete(labels=c("Cm", "Sp_Y", "Sp_I1", "Sp_I2", "Sp_O")) +
  mytheme +
  theme(axis.text.x =element_text(angle=30, hjust=1, vjust=1))
predfig.endad1

fety2[,"N"]




