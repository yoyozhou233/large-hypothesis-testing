---
  title: "project 1"
output: html_document
date: "2022-11-16"
---
library(ggplot2)

n <- 100
R <- 10000
I <- 30
alpha <- 0.05
q <- 0.1


N <- c(50, 100, 200, 300, 500, 700, 1000)
Is <- c(20, 50, 100, 150, 200)
FWER_Bon <- c()
TRR_Bon <- c()
FWER_Holm <- c()
TRR_Holm <- c()
TRR_BH <- c()
FDR_BH <- c()

FWER_BonI <- c()
TRR_BonI <- c()
FWER_HolmI <- c()
TRR_HolmI <- c()
TRR_BHI <- c()
FDR_BHI <- c()
P<- 400
mus <- c(0.5,1)

FWER_Bon_mu <- c()
TRR_Bon_mu <- c()
FWER_Holm_mu <- c()
TRR_Holm_mu <- c()
TRR_BH_mu <- c()
FDR_BH_mu <- c()

for (p in N) 
{mu <- c(rep(0.5, I), rep(0, p-I))
H0False <- c(rep(1,I), rep(0, p-I))


X <- matrix(rnorm(p*R, 0 , 1/sqrt(n)), p, R) + mu
pValue <- 2*(1-pnorm(abs(X), 0,  1/sqrt(n)))


#Bonferroni

TestBon <- (pValue <= alpha/p)
FRBon <- colSums(TestBon[(I+1):p,])
FWEBon <- colSums(TestBon[(I+1):p,]) >= 1
FWERBon <- mean(FWEBon)
TRBon <- colSums(TestBon[1:I,])
TRRBon <- mean(TRBon)

FWER_Bon <- c(FWER_Bon, FWERBon)
TRR_Bon <- c(TRR_Bon, TRRBon)

pValue_ordered <- apply(pValue, 2, sort)
pValue_order <- apply(pValue, 2, order)

#Holm's

m <-  rep(0, p)
for (i in 1:p) m[i] <- alpha/(p-i+1)
TestHolm <- (pValue_ordered <= m)
for (i in 1:R) {for (j in 1:p) if (TestHolm[j, i] == F) {TestHolm[j,] == F}}

FRHolm <- rep(0, R)
for (i in 1:R)
  FRHolm[i] <- sum(TestHolm[which(pValue_order[,i] > I),i])
FWEHolm<- FRHolm >= 1
FWERHolm <- mean(FWEHolm)

TRHolm <- rep(0, R)
for (i in 1:R)
  TRHolm[i] <- sum(TestHolm[which(pValue_order[,i] <= I),i])
TRRHolm <- mean(TRHolm)

FWER_Holm <- c(FWER_Holm, FWERHolm)
TRR_Holm <- c(TRR_Holm, TRRHolm)

#B-H

c <- rep(0, p)
for (i in 1:p) c[i] <- i*q/p
TestBH <- (pValue_ordered <= c)
for (i in 1:R) {for (j in 1:p) if (TestBH[j, i] == F) {TestBH[j,] == F}}

TotalR_BH <- colSums(TestBH)
FRBH <- rep(0, R)
for (i in 1:R)
  FRBH[i] <- sum(TestBH[which(pValue_order[,i] > I),i])
FDR <- mean(FRBH/TotalR_BH)

TRBH <- rep(0, R)
for (i in 1:R)
  TRBH[i] <- sum(TestBH[which(pValue_order[,i] <= I),i])
TRRBH <- mean(TRBH)

FDR_BH <- c(FDR_BH, FDR)
TRR_BH <- c(TRR_BH, TRRBH)
}
for (i in Is) 
{mu <- c(rep(0.5, i), rep(0, P-i))
H0False <- c(rep(1,i), rep(0, P-i))


X <- matrix(rnorm(P*R, 0 , 1/sqrt(n)), P, R) + mu
pValue <- 2*(1-pnorm(abs(X), 0,  1/sqrt(n)))


#Bonferroni

TestBon <- (pValue <= alpha/P)
FRBon <- colSums(TestBon[(i+1):P,])
FWEBon <- colSums(TestBon[(i+1):P,]) >= 1
FWERBon <- mean(FWEBon)
TRBon <- colSums(TestBon[1:i,])
TRRBon <- mean(TRBon)

FWER_BonI <- c(FWER_BonI, FWERBon)
TRR_BonI <- c(TRR_BonI, TRRBon)

pValue_ordered <- apply(pValue, 2, sort)
pValue_order <- apply(pValue, 2, order)

#Holm's

m <-  rep(0, P)
for (j in 1:P) m[j] <- alpha/(P-j+1)
TestHolm <- (pValue_ordered <= m)
for (j in 1:R) {for (k in 1:P) if (TestHolm[k, j] == F) {TestHolm[k,] == F}}

FRHolm <- rep(0, R)
for (j in 1:R)
  FRHolm[j] <- sum(TestHolm[which(pValue_order[,j] > i),j])
FWEHolm<- FRHolm >= 1
FWERHolm <- mean(FWEHolm)

TRHolm <- rep(0, R)
for (j in 1:R)
  TRHolm[j] <- sum(TestHolm[which(pValue_order[,j] <= i),j])
TRRHolm <- mean(TRHolm)

FWER_HolmI <- c(FWER_HolmI, FWERHolm)
TRR_HolmI <- c(TRR_HolmI, TRRHolm)

#B-H

c <- rep(0, P)
for (j in 1:P) c[j] <- j*q/P
TestBH <- (pValue_ordered <= c)
for (j in 1:R) {for (k in 1:P) if (TestBH[k, j] == F) {TestBH[k,] == F}}

TotalR_BH <- colSums(TestBH)
FRBH <- rep(0, R)
for (j in 1:R)
  FRBH[j] <- sum(TestBH[which(pValue_order[,j] > i),j])
FDR <- mean(FRBH/TotalR_BH)

TRBH <- rep(0, R)
for (j in 1:R)
  TRBH[j] <- sum(TestBH[which(pValue_order[,j] <= i),j])
TRRBH <- mean(TRBH)

FDR_BHI <- c(FDR_BHI, FDR)
TRR_BHI <- c(TRR_BHI, TRRBH)
}

for (Mu in mus) 
{mu <- c(rep(Mu, I), rep(0, P-I))
H0False <- c(rep(1,I), rep(0, P-I))

X <- matrix(rnorm(P*R, 0 , 1/sqrt(n)), P, R) + mu
pValue <- 2*(1-pnorm(abs(X), 0,  1/sqrt(n)))


#Bonferroni

TestBon <- (pValue <= alpha/P)
FRBon <- colSums(TestBon[(I+1):P,])
FWEBon <- colSums(TestBon[(I+1):P,]) >= 1
FWERBon <- mean(FWEBon)
TRBon <- colSums(TestBon[1:I,])
TRRBon <- mean(TRBon)

FWER_Bon_mu <- c(FWER_Bon_mu, FWERBon)
TRR_Bon_mu <- c(TRR_Bon_mu, TRRBon)

pValue_ordered <- apply(pValue, 2, sort)
pValue_order <- apply(pValue, 2, order)

#Holm's

m <-  rep(0, P)
for (j in 1:P) m[j] <- alpha/(P-j+1)
TestHolm <- (pValue_ordered <= m)
for (j in 1:R) {for (k in 1:P) if (TestHolm[k, j] == F) {TestHolm[k,] == F}}

FRHolm <- rep(0, R)
for (j in 1:R)
  FRHolm[j] <- sum(TestHolm[which(pValue_order[,j] > I),j])
FWEHolm<- FRHolm >= 1
FWERHolm <- mean(FWEHolm)

TRHolm <- rep(0, R)
for (j in 1:R)
  TRHolm[j] <- sum(TestHolm[which(pValue_order[,j] <= I),j])
TRRHolm <- mean(TRHolm)

FWER_Holm_mu <- c(FWER_Holm_mu, FWERHolm)
TRR_Holm_mu <- c(TRR_Holm_mu, TRRHolm)

#B-H

c <- rep(0, P)
for (j in 1:P) c[j] <- j*q/P
TestBH <- (pValue_ordered <= c)
for (j in 1:R) {for (k in 1:P) if (TestBH[k, j] == F) {TestBH[k,] == F}}

TotalR_BH <- colSums(TestBH)
FRBH <- rep(0, R)
for (j in 1:R)
  FRBH[j] <- sum(TestBH[which(pValue_order[,j] > I),j])
FDR <- mean(FRBH/TotalR_BH)

TRBH <- rep(0, R)
for (j in 1:R)
  TRBH[j] <- sum(TestBH[which(pValue_order[,j] <= I),j])
TRRBH <- mean(TRBH)

FDR_BH_mu <- c(FDR_BH_mu, FDR)
TRR_BH_mu <- c(TRR_BH_mu, TRRBH)
}

#Plot for FWER control

#N plot
plot(N, FWER_Bon, ylim = c(0, 0.06), ylab = "FWER", col = "blue", pch = 4)
points(N, FWER_Holm, ylim = c(0, 0.06), col = "red", pch = 4)
legend(x = "bottomright", legend = c("Bonferroni", "Holms"), col = c("blue", "red"), pch = 4)
abline(h = alpha, lty = 2)
#I plot
plot(Is, FWER_IBon, ylim = c(0, 0.06), xlab = "I", ylab = "FWER", col = "blue", pch = 4)
points(Is, FWER_IHolm, ylim = c(0, 0.06), col = "red", pch = 4)
legend(x = "bottomright", legend = c("Bonferroni", "Holms"), col = c("blue", "red"), pch = 4)
abline(h = alpha, lty = 2)
#Plots for FDR control
#N changing plot
plot(N, FDR_BH, ylim = c(0, 0.1), col = "blue", pch = 4)
abline(h = q, lty = 2, col = 'red')
#I changing plot
plot(Is, FDR_IBH, ylim = c(0, 0.1), xlab = "I", ylab = "FDR",col = "blue", pch = 4)
abline(h = q, lty = 2, col = 'red')

Bonferroni_N <- TRR_Bon
Holm_N <- TRR_Holm
b <- rbind("N" = signif(N, 3), "Bonferroni" = signif(Bonferroni_N, 3), "Holm's" = signif(Holm_N, 3))
True_Rejection <- data.frame(b)

BH_N <- TRR_BH
c <- rbind("N" = signif(N, 3), "Bonferroni" = signif(Bonferroni_N, 3), "Holm's" = signif(Holm_N, 3), "Benjamini-Hochberg" = signif(BH_N, 3))
TRBH <- data.frame(c)

FDR_N_Est <- data.frame(N, "FDR" = FDR_BH, Theory = ((N-I)*q)/N)
TRR_I_Change <- data.frame("I" = Is, "Bonferroni" = signif(TRR_IBon, 3), "Holms" = signif(TRR_IHolm,3), "Benjamini-Hochberg" = signif(TRR_IBH,3))
FDR_I_Est <- data.frame("I" = Is, "FDR" = FDR_IBH, Theory = ((400-Is)*q)/400)
mu_change <- data.frame("mu" = mus, "Bonferroni" = FWER_Bon_mu, "Holm's" = FWER_Holm_mu)
TRR_mu_change <- data.frame(mus, "Bonferroni" = signif(TRR_Bon_mu, 3), "Holm's" = signif(TRR_Holm_mu, 3))
FDRTRR_mu_change <- data.frame(mus, "Bonferroni" = signif(TRR_Bon_mu, 3), "Holm's" = signif(TRR_Holm_mu, 3), "Benjamini-Hochberg" = signif(TRR_BH_mu, 3))
FDR_mu_change <- data.frame("mu" = mus, "Bonferroni" = FWER_Bon_mu, "Holm's" = FWER_Holm_mu, "Benjamini-Hochberg" = FDR_BH_mu)

TRR_I_Change <- data.frame(Is, TRR_BonI, "Percentage.Bonferroni" = 100*TRR_BonI/Is, TRR_HolmI, "Percentage.Holm's" = 100*TRR_HolmI/Is) 


