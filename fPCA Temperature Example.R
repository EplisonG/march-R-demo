library(fda)
library(ggplot2)
library(reshape2)


daybasis365 <- create.fourier.basis(c(0,365), 365)
temp <- daily$tempav
temp.long <- melt(temp, id.vars=c("days", "city","temp"))


ggplot(data = temp.long, aes(x=Var1, y=value,group = Var2, color = Var2)) +geom_line(show.legend = FALSE)+
  xlab("Days") +ylab("Daily Temperature") +theme_bw()

#  carry out a PCA of temperature
daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
#  penalize harmonic acceleration, use varimax rotation

harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), c(0, 365))
harmfdPar <- fdPar(daybasis65, harmaccelLfd, lambda=1e5)
daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],
                          daybasis65, fdnames=list("Day", "Station", "Deg C"))$fd

plot(daytempfd, xlab = 'day', ylab = "Daily Temperature", cex.lab =1.5)

# variance matrix
temp.var <- var.fd(daytempfd)
tvvals=eval.bifd(1:365, 1:365, temp.var)
contour(1:365, 1:365, tvvals, xlab='day', ylab = 'day', cex.lab =1.5)
library(fields)
image.plot(1:365, 1:365, tvvals, xlab ="Day", ylab = "Day", cex.lab =1.5)

# Do functional principal components analysis
daytemppcaobj <- pca.fd(daytempfd, nharm=4, harmfdPar)
names(daytemppcaobj)
plot(daytemppcaobj$values[1:8], xlab = 'component', ylab = 'variance')

#show the mean curve 
plot(daytemppcaobj$meanfd)
#functional principal components
harmfd <- daytemppcaobj$harmonics
harmvals <- eval.fd(1:365, harmfd)

# plot the top 4 FPCs
par(mfrow=c(2,2))
plot(1:365, harmvals[,1], xlab = 'day', ylim = c(-0.1, 0.1), ylab = "Principal Components", col="blue",  lwd =2, type = "l")
abline(h=0, col="red", lty=2)
plot(1:365, harmvals[,2], xlab = 'day', ylim = c(-0.1, 0.1), ylab = "Principal Components", col="blue",  lwd =2, type = "l")
abline(h=0, col="red", lty=2)
plot(1:365, harmvals[,3], xlab = 'day', ylim = c(-0.1, 0.1), ylab = "Principal Components", col="blue",  lwd =2, type = "l")
abline(h=0, col="red", lty=2)
plot(1:365, harmvals[,4], xlab = 'day', ylim = c(-0.1, 0.1), ylab = "Principal Components", col="blue",  lwd =2, type = "l")
abline(h=0, col="red", lty=2)


# plot the 1 vs 2 scores
par(mfrow=c(1,1))
plot(daytemppcaobj$scores[, 1:2], xlab= 'PC Score 1', ylab = 'PC Score 2', col="blue", cex =1)
text(daytemppcaobj$scores[,1], daytemppcaobj$scores[,2], labels = daily$place, cex = 1)

#  plot harmonics
op <- par(mfrow=c(2,2))
plot.pca.fd(daytemppcaobj, cex.main=0.9)
#rotation of PC
daytemppcaVarmx <- varmx.pca.fd(daytemppcaobj)
plot.pca.fd(daytemppcaVarmx, cex.main=0.9)
par(op)

plot(daytemppcaobj$harmonics)
plot(daytemppcaVarmx$harmonics)
