library(fda)
library(fda.usc)
library(ggplot2)

data("CanadianWeather")
str(CanadianWeather)

#Let's focus on log-transformed average daily precipitation of Vancouver (one observed curve).

y.lprec <- CanadianWeather$dailyAv[,,3]
l <- which(CanadianWeather$place == "Vancouver")
y <- y.lprec[,l]
day <- 1:365
plot(day, y, 
     type='o', pch = 16, cex = 0.5, col='royalblue2',
     xlab="day", ylab="log-precipitation",
     main="Log of daily average precipitation of Vancouver")
plot(day, y,  type='o', col="blue",
     xlab="day", ylab="precipitation",
     main="Daily precipitation of Vancouver")
# Fourier basis
T=range(day)
basis.fourier3 <- create.fourier.basis(rangeval = T, nbasis = 3)
bvals3 <- eval.basis(day, basis.fourier3)
lm.fit = lm(y ~ 0 + bvals3)
y.fit = lm.fit$fitted.values; coef= lm.fit$coefficient

par(mfrow=c(1,2))
plot(basis.fourier3)
plot(day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="precipitation")
points(day, y, col="blue", lwd=0.5)
lines(day, lm.fit$fitted.values, lwd=2, col="red")
mtext("Fourier Basis K = 3", outer = TRUE, cex = 1.2, line = -2)


basis.fourier7 <- create.fourier.basis(rangeval = T, nbasis = 9)
bvals7 <- eval.basis(day, basis.fourier7)
lm.fit = lm(y ~ 0 + bvals7)
y.fit = lm.fit$fitted.values; coef= lm.fit$coefficient


par(mfrow=c(1,2))
plot(basis.fourier7)
plot(day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="precipitation")
points(day, y, col="blue", lwd=0.5)
lines(day, lm.fit$fitted.values, lwd=2, col="red")
mtext("Fourier Basis K = 9", outer = TRUE, cex = 1.2, line = -2)


# B-spline basis
T=range(day)
basis.spline3 <- create.bspline.basis(rangeval = T, nbasis = 9, norder = 4)
bs3 <- eval.basis(day, basis.spline3)
lm.fit = lm(y ~ 0 + bs3)
y.fit = lm.fit$fitted.values; coef= lm.fit$coefficient

par(mfrow=c(1,2))
plot(basis.spline3)
plot(day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="precipitation")
points(day, y, col="blue", lwd=0.5)
lines(day, lm.fit$fitted.values, lwd=2, col="red")
mtext("B-spline Basis K = 9, m = 4", outer = TRUE, cex = 1.2, line = -2)


basis.spline7 <- create.bspline.basis(rangeval = T, nbasis = 19, norder = 4)
bs7 <- eval.basis(day, basis.spline7)
lm.fit = lm(y ~ 0 + bs7)
y.fit = lm.fit$fitted.values; coef= lm.fit$coefficient

par(mfrow=c(1,2))
plot(basis.spline7)
plot(day, y, type="n",lwd=4, col="black",
     xlab="day", ylab="precipitation")
points(day, y, col="blue", lwd=0.5)
lines(day, lm.fit$fitted.values, lwd=2, col="red")
mtext("B-spline Basis K = 19, m = 4 ", outer = TRUE, cex = 1.2, line = -2)


