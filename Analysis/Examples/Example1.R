

# Factors that Subject is nested within, like Training group, are called between-subjects factors. 
# Factors that Subject is crossed with, like Time, are called within-subjects factors.


set.seed(123)
P     <- 2               # Xb1 b= between subjects, subject is nested in this factor
Q     <- 2               # Xb2
R     <- 3               # Xw1 w= within subject, subject is crossed with this facor
S     <- 3               # Xw2
Njklm <- 20              # obs per cell; Cell is what?
Njk   <- Njklm*P*Q       # number of subjects
N     <- Njklm*P*Q*R*S   # number of observations
id    <- gl(Njk,         R*S, N, labels=c(paste("s", 1:Njk, sep="")))
Xb1   <- gl(P,   Njklm*Q*R*S, N, labels=c("CG", "T")) # Control, Treatment
Xb2   <- gl(Q,   Njklm  *R*S, N, labels=c("f", "m"))  # male, female
Xw1   <- gl(R,             S, N, labels=c("A", "B", "C")) # every subject should have A, B, and C
Xw2   <- gl(S,   1,           N, labels=c("-", "o", "+")) # every subject should have....

mu      <- 100
eB1     <- c(-5, 5)
eB2     <- c(-5, 5)
eW1     <- c(-5, 0, 5)
eW2     <- c(-5, 0, 5)
eB1B2   <- c(-5, 5, 5, -5)
eB1W1   <- c(-5, 5, 2, -2, 3, -3)
eB1W2   <- c(-5, 5, 2, -2, 3, -3)
eB2W1   <- c(-5, 5, 2, -2, 3, -3)
eB2W2   <- c(-5, 5, 2, -2, 3, -3)
eW1W2   <- c(-5, 2, 3, 2, 3, -5, 2, -5, 3)
eB1B2W1 <- c(-5, 5, 5, -5, 2, -2, -2, 2, 3, -3, -3, 3)
eB1B2W2 <- c(-5, 5, 5, -5, 2, -2, -2, 2, 3, -3, -3, 3)
eB1W1W2 <- c(-5, 5, 2, -2, 3, -3, 3, -3, -5, 5, 2, -2, 2, -2, 3, -3, -5, 5)
eB2W1W2 <- c(-5, 5, 2, -2, 3, -3, 3, -3, -5, 5, 2, -2, 2, -2, 3, -3, -5, 5)

names(eB1)     <- levels(Xb1)
names(eB2)     <- levels(Xb2)
names(eW1)     <- levels(Xw1)
names(eW2)     <- levels(Xw2)
names(eB1B2)   <- levels(interaction(Xb1, Xb2))
names(eB1W1)   <- levels(interaction(Xb1, Xw1))
names(eB1W2)   <- levels(interaction(Xb1, Xw2))
names(eB2W1)   <- levels(interaction(Xb2, Xw1))
names(eB2W2)   <- levels(interaction(Xb2, Xw2))
names(eW1W2)   <- levels(interaction(Xw1, Xw2))
names(eB1B2W1) <- levels(interaction(Xb1, Xb2, Xw1))
names(eB1B2W2) <- levels(interaction(Xb1, Xb2, Xw2))
names(eB1W1W2) <- levels(interaction(Xb1, Xw1, Xw2))
names(eB2W1W2) <- levels(interaction(Xb2, Xw1, Xw2))

# hier wird quasi das lineare Modell geschrieben
muJKLM <- mu +
  eB1[Xb1] + eB2[Xb2] + eW1[Xw1] + eW2[Xw2] +
  eB1B2[interaction(Xb1, Xb2)] +
  eB1W1[interaction(Xb1, Xw1)] +
  eB1W2[interaction(Xb1, Xw2)] +
  eB2W1[interaction(Xb2, Xw1)] +
  eB2W2[interaction(Xb2, Xw2)] +
  eW1W2[interaction(Xw1, Xw2)] +
  eB1B2W1[interaction(Xb1, Xb2, Xw1)] +
  eB1B2W2[interaction(Xb1, Xb2, Xw2)] +
  eB1W1W2[interaction(Xb1, Xw1, Xw2)] +
  eB2W1W2[interaction(Xb2, Xw1, Xw2)]
muId  <- rep(rnorm(Njk, 0, 3), each=R*S)
mus   <- muJKLM + muId # lineares Model plus Zufallszahl pro achtzig subjects, 80 = R*S, 
sigma <- 50

Y  <- round(rnorm(N, mus, sigma), 1)
d2 <- data.frame(id, Xb1, Xb2, Xw1, Xw2, Y)

str(d2)
head(d2)

d1 <- aggregate(Y ~ id + Xw1 + Xb1 + Xb2, data=d2, FUN=mean)

# conventional Analysis using aov()
summary(aov(Y ~ Xw1 + Error(id/Xw1), data=d1))
