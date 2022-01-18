# Title     : Local-Level Model with explanatory variable
# Objective : Corona Inzidenz und tageshöchst Temperatur
# Created by: Timo Kuehne, Jonathan Laib
# Created on: 01.04.2021
#Paket einbinden
library(dlm)
#Daten einlesen und Transformieren
data <- read.csv("CoronaTempWoche.CSV", header=TRUE, sep = ";")
corona <- ts(log(data$Neuinfektionen))
temperatur <- ts(data$Temperatur)
#Modell erstellen
fn <- function(parameter){
  mod <- dlmModPoly(order = 1) + dlmModReg(temperatur, addInt = FALSE)
  V(mod) <- exp(parameter[1])
  diag(W(mod))[1:2] <- exp(parameter[2:3])
  return(mod)
}
#Maximum Likelihood Schätzung
fit <- dlmMLE(corona, rep(0,3), fn)
mod <- fn(fit$par)
#AIC
AIC(lm(temperatur ~ corona))
loglik <- dlmLL(corona, dlmModPoly(1) + dlmModReg(temperatur, addInt = FALSE))
n.coef <- 2 + 2
r.aic <- (2 * (loglik)) + 2 * (sum(n.coef))
#Varianzen auslesen
obs.error.var <- V(mod)
state.error.var <- diag(W(mod))
#Kalman Filter und Smoother
filtered <- dlmFilter(corona, mod = mod)
smoothed <- dlmSmooth(filtered)
resids <- residuals(filtered, sd = FALSE)
#Restlichen Werte auslesen
mu <- dropFirst(smoothed$s[, 1])
betas <- dropFirst(smoothed$s[, 2])
mu.1 <- mu[1]
mu.end <- mu[length(mu)]
beta.1 <- betas[1]
beta.end <- betas[length(mu)]
#Diagramme und Graphen
par(mfrow=c(1,1))
val <- mu + betas*temperatur
plot.ts(ts(cbind(corona,val)), plot.type="single" , col =c("darkgrey","blue"),
        lty=c(1,2), xlab="Wochen nach 30.03.2020", ylab = "log Corona Neuinfektionen")
legend("bottomright",
       leg = c("log Corona Neuinfektionen","stochastic level + beta*Temperatur"),
       cex = 1, lty = c(1, 2), col = c("darkgrey","blue"))
par(mfrow=c(1,1))
plot.ts(resids, col = "darkgrey", ylab="", xlab="Wochen nach 30.03.2020", cex.main = 0.8)
abline(h=0, col = "sienna")
legend("bottomright",
       leg = "irregular ",
       cex = 1, lty = 1, col = "darkgrey")
par(mfrow=c(1,1))
plot.ts(corona - (mu + betas*temperatur), ylim=c(-0.000005,0.000005),
        col = "darkgrey", ylab="", xlab="Wochen nach 30.03.2020", cex.main = 0.8)
abline(h=0, col = "sienna")
legend("bottomright",
       leg = "irregular ",
       cex = 1, lty = 1, col = "darkgrey")







