




#pvalues <- rbind(p.fety.aov,p.fety.adis0,p.fety.lm, p.fety.adis,p.fety.adissc,p.fety.biglm,p.fety.biglmer,p.fety.lmer,p.fety.nparLD)
pvalues <- pvalues[-10,]
pvalues <- data.frame(pvalues, stringsAsFactors = FALSE)
test <- rep(c("anova",
              "pAnova0",
              "arcsinelm",
              "pAnova.field",
              "pAnova.samcam",
              "binomialGLM",
              "binomialGLMM",
              "ArcsineLMM",
              "nparLD"
),each=3)
pvalues$test <- test

# Compare normal and arcsine aov, repmes aov and binomial models and non parametric, without CM class
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("pValuesNorm.rda")
load("pValuesBin.rda")

pvalues <- rbind(p.fety.aov,p.fety.aovas,p.fety.rmaov, p.fety.rmaovas,p.fety.nparLD,p.fety.biglm,p.fety.biglmer)
pvalues <- data.frame(pvalues, stringsAsFactors = FALSE)
for(i in 2:10) {
  pvalues[,i] <-  as.numeric(pvalues[,i])
}
test <- rep(c("ANOVA",
              "ANOVAarcsine",
              "rmANOVA",
              "rmANOVAarcsine",
              "nparLD",
              "binomialGLM",
              "binomialGLMM"
),each=3)
pvalues$test <- test
pvalues$Env <- plyr::revalue(pvalues$Env, replace=c("group"="age_class", "time"="samcam", "group:time"="age_class:samcam"))
str(pvalues)


pvalues2 <- melt(pvalues,id.vars = c(1,2,11))
pvalues2$value <- as.numeric(pvalues2$value)

ggplot(pvalues2, aes(Env, value)) + 
  geom_point() +
  geom_point(aes(col=value<0.05)) +
  facet_grid(variable ~ test)


