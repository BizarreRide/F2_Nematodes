

df.ex1 <- data.frame(age_class = df.exp$age_class,
                     samcam = df.exp$samcam,
                      field.ID = df.exp$field.ID,
                      ID = df.exp$ID)
                      
df.rsp1 <- data.frame(CI = df.rspX$CI,
                      Apc = df.rspX$Aporcelaimidae,
                      Hoplo = df.rspX$Hoplolaimidae)

outlier <- list(nema.CI <- 1:24, #-24?
                spec.Apc <- 1:24,
                spec.Hoplo <- -2)

# Formulas for "global models"
fml.glb <- as.formula(y ~ age_class + samcam + age_class:samcam + (1|field.ID)) # + pH + cnratio + mc + ats1 
fml.glb.log <- as.formula(log1p(y) ~ age_class + samcam + age_class:samcam + (1|field.ID)) #+ pH + cnratio + mc + ats1 

# Set optimizer function gor glmm
con = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs = TRUE, check.conv.grad="ignore")

# Case1
df.ex1$y <- df.rsp1$CI
df.ex2 <- df.ex1[outlier[[1]],]
M1.A <- lmer(fml.glb, df.ex2)
M1.B <- lmer(fml.glb.log, df.ex2)

#vCase2
df.ex1$y <- df.rsp1$Apc
df.ex2 <- df.ex1[outlier[[2]],]
M2.A <- lmer(fml.glb, df.ex2)
M2.B <- glmer(fml.glb, family = "poisson", df.ex2, control = con)
M2.C <- lmer(fml.glb.log, df.ex2)

# Case3
df.ex1$y <- df.rsp1$Hoplo
df.ex2 <- df.ex1[outlier[[3]],]
M3.A <- lmer(fml.glb, df.ex2)
M3.B <- glmer(fml.glb, family = "poisson", df.ex2, control = con)
M3.C <- lmer(fml.glb.log, df.ex2)

# Store Models in a list
ls.glbModels <- list(M1.A, M1.B, M2.A, M2.B, M2.C, M3.A, M3.B, M3.C)
names(ls.glbModels) <- c("M1.A_CI", "M1.B_CI", "M2.A_Apc", "M2.B_Apc", "M2.C_Apc", "M3.A_Hoplo", "M3.B_Hoplo", "M3.C_Hoplo")



windows(record=TRUE)

for(k in 1:8){ 

  E1 <- resid(ls.glbModels[[k]], type="pearson")
  E2 <- resid(ls.glbModels[[k]], type="response")
  E3 <- residuals(ls.glbModels[[k]], type="deviance")
  F1 <- fitted(ls.glbModels[[k]], type="response")
  P1 <- predict(ls.glbModels[[k]], type="response")
  
  par(mfrow=c(3,3),
      mar=c(4,4.5,1,2),
      oma=c(0,0,2,0))
 
  # Plot fitted vs. residuals
  p1 <- scatter.smooth(F1, E1, cex.lab = 1.5, xlab="Fitted values", ylab="Pearson Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.ex1), pos = 4)
  scatter.smooth(F1, E2, cex.lab = 1.5, xlab="Fitted values", ylab="Response Residuals")
  abline(h = 0, v=0, lty=2); text(F1, E1, labels = row.names(df.ex1), pos = 2)
  scatter.smooth(F1, E3, cex.lab = 1.5, xlab="Fitted values", ylab="Deviance Residuals")
  abline(h = 0, v=0, lty=2)
  
  # plot predicted vs. residuals
  scatter.smooth(log(predict(ls.glbModels[[k]])), E1, cex.lab = 1.5, xlab="log Predicted values", ylab="Pearson Residuals")
  abline(h = 0, v=0, lty=2); text(log(predict(ls.glbModels[[k]])), E1, labels = row.names(df.ex1), pos = 4)
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
  
  qqnorm(E3)
  qqline(E3)
  
  title(paste(k,names(ls.glbModels)[k]), outer=TRUE)
}


