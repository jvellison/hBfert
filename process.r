##Process Model hB forecasts
#Plotting set-up
library(sp)
bpy_org     <- bpy.colors(15)[12]
bpy_blue    <- bpy.colors(15)[3]
windowsFonts(Calibri = windowsFont("Calibri"))

#General set-up
ages        <- 15:44
hBforecast  <- dget(paste(this.dir, "/", forecastFile, "/", "forecast.dput", sep = "")) #load forecast file
name.red    <- names(hBforecast)[which(unlist(lapply(hBforecast, length))  == 28)] #remove countries without forecasts from "name"
namef.red   <- unlist(lapply(hBforecast, function(x) x$namef))[name.red] #remove countries without forecasts from "namef"
n.countries <- length(name.red) #number of countries
comm.ind    <- unlist(lapply(hBforecast, function(x) x$comm.ind))[name.red] #communist indicator

## 1. Generate CRPS and LogS and save plots of the average CRPS and LogS by country
#First, check how many countries have at least 5 scores (equivalently at least 5 observed future CFR values)
#so that the average is reliable
n.cfr          <- unlist(lapply(hBforecast, function(x) x$n.cfr)) #Extract number of observed future CFR values for each country
name.score     <- name.red[n.cfr > 4]   #extract country names with at least 5 scores
namef.score    <- namef.red[name.score] #extract full country names with at least 5 scores
comm.ind.score <- comm.ind[name.score]  #extract comm.ind for countries with at least 5 scores
n.score        <- length(name.score)    #extract the number of countries with at least 5 scores

#Only generate score plots if there is at least one country with at least 5 scores (n.score > 0)
if (n.score > 0) {
#Extract Model hB CRPS scores for each country from hBforecast
  CRPShB <- matrix(NA, n.score, 40)
  rownames(CRPShB) <- name.score
  colnames(CRPShB) <- con
  for (k in name.score)
    CRPShB[k, hBforecast[[k]]$coh.cfr] <- hBforecast[[k]]$CRPS[hBforecast[[k]]$coh.cfr]
  CRPShBavg <- apply(CRPShB, 1, mean, na.rm = T)

#Compute freeze rates CFR forecasts
freezerates <- list()
for (k in name.score) {
  freezerates[[k]] <- hBforecast[[k]]$presrates
  for (i in 2:30) {
    freezerates[[k]][i, (length(which(!is.na(freezerates[[k]][i,])))+1):40] <- freezerates[[k]][i, length(which(!is.na(freezerates[[k]][i,])))] 
  }
  if (hBforecast[[k]]$n.complete==10)
    freezerates[[k]][1, 40] <- freezerates[[k]][1, 39]
}

#Compute freeze rates CRPS scores for each country (absolute error as deterministic)
#We also compute the absolute percentage error here for convenience
AEfreeze <- matrix(NA, n.score, 40)
rownames(AEfreeze) <- name.score
colnames(AEfreeze) <- con
APEfreeze <- AEfreeze
for (k in name.score) {
  AEfreeze[k, hBforecast[[k]]$coh.cfr]  <- abs(hBforecast[[k]]$obs.cfr - apply(freezerates[[k]], 2, sum)[hBforecast[[k]]$coh.cfr])
  APEfreeze[k, hBforecast[[k]]$coh.cfr] <- abs(hBforecast[[k]]$obs.cfr - apply(freezerates[[k]], 2, sum)[hBforecast[[k]]$coh.cfr])/hBforecast[[k]]$obs.cfr
}
AEfreezeavg <- apply(AEfreeze, 1, mean, na.rm = T)

#Compute Model hB LogS scores for each country using a Gaussian approximation
LogShB <- matrix(NA, n.score, 40)
rownames(LogShB) <- name.score
colnames(LogShB) <- con
for (k in name.score)
  LogShB[k, hBforecast[[k]]$coh.cfr] <- -log(dnorm(hBforecast[[k]]$obs.cfr, hBforecast[[k]]$cfrmean[hBforecast[[k]]$coh.cfr], hBforecast[[k]]$cfrsd[hBforecast[[k]]$coh.cfr]))
LogShBavg <- apply(LogShB, 1, mean, na.rm = T)

#Save plot of average CRPS scores for Model hB and freeze rates for all countries (see Figure 5)
png(file=paste(forecastFile, "/CRPS_plot.png", sep = ""), width = 1000, height = 500, family = "Calibri", res = 100)
ind <- order(CRPShBavg, decreasing = T)
par(mfrow = c(1, 1), mar = c(8, 4, 1, 1), family = "Calibri")
plot(1:n.score, rep(0, n.score), type = "n", las = 1, xaxt = "n", xlab = "", ylab = "Average CRPS",
     ylim = range(c(CRPShBavg, AEfreezeavg)), main = " ")
abline(v=seq(0.5, n.score + 0.5, 1), lty = 2, col = "grey")
points(1:n.score, AEfreezeavg[ind], pch =  8, type="o")
points(1:n.score, CRPShBavg[ind],   pch = 16, type="o", col = bpy_blue)
axis(side = 1, at = which(comm.ind.score[ind] == 0), label = (namef.score[ind])[which(comm.ind.score[ind] == 0)], las=2, col = 1)
axis(side = 1, at = which(comm.ind.score[ind] == 1), label = (namef.score[ind])[which(comm.ind.score[ind] == 1)], las=2, col.axis = 1, font = 2)
dev.off()

#Save plot of average LogS scores for Model hB for all countries (see Figure 6)
png(file=paste(forecastFile, "/LogS_plot.png", sep = ""), width = 1000, height = 500, family = "Calibri", res = 100)
ind <- order(LogShBavg, decreasing = T)
par(mfrow = c(1, 1), mar = c(8, 4, 1, 1), family = "Calibri")
plot(1:n.score, rep(0, n.score), type = "n", las = 1, xaxt = "n", xlab = "", ylab = "Average LogS",
     ylim = range(LogShBavg), main = " ")
abline(v=seq(0.5, n.score + 0.5, 1), lty = 2, col = "grey")
points(1:n.score, LogShBavg[ind], pch = 16, type="o", col = bpy_blue)
axis(side = 1, at = which(comm.ind.score[ind] == 0), label = (namef.score[ind])[which(comm.ind.score[ind] == 0)], las=2, col = 1)
axis(side = 1, at = which(comm.ind.score[ind] == 1), label = (namef.score[ind])[which(comm.ind.score[ind] == 1)], las=2, col.axis = 1, font = 2)
dev.off()
}

## 2. Print summary statistics
#MAE, MAPE, RMSE, 90% coverage (COV90hB) and 50% coverage (COV50hB) 
#Only compute summary statistics if there is at least one country with at least 5 observed future CFR values (n.score > 0)
if (n.score > 0) {
  #Model hB
  AEhB <- matrix(NA, n.score, 40)
  rownames(AEhB) <- name.score
  colnames(AEhB) <- con
  APEhB <- AEhB
  for (k in name.score) {
    AEhB[k, hBforecast[[k]]$coh.cfr] <-  abs(hBforecast[[k]]$obs.cfr - hBforecast[[k]]$cfrq["50%",hBforecast[[k]]$coh.cfr])
    APEhB[k, hBforecast[[k]]$coh.cfr] <- abs(hBforecast[[k]]$obs.cfr - hBforecast[[k]]$cfrq["50%",hBforecast[[k]]$coh.cfr])/hBforecast[[k]]$obs.cfr
  }
  MAEhB  <- mean(AEhB, na.rm = T)
  MAPEhB <- mean(APEhB,  na.rm = T)*100
  RMSEhB <- sqrt(mean(AEhB^2, na.rm = T))
  
  allobscfr <- numeric()
  for (k in name.score)
    allobscfr <- c(allobscfr, hBforecast[[k]]$obs.cfr)
  allhBcfrq <- list()
  for (k in name.score) {
    allhBcfrq$q05 <- c(allhBcfrq$q05, hBforecast[[k]]$cfrq["5%",  hBforecast[[k]]$coh.cfr])
    allhBcfrq$q95 <- c(allhBcfrq$q95, hBforecast[[k]]$cfrq["95%", hBforecast[[k]]$coh.cfr])
    allhBcfrq$q25 <- c(allhBcfrq$q25, hBforecast[[k]]$cfrq["25%", hBforecast[[k]]$coh.cfr])
    allhBcfrq$q75 <- c(allhBcfrq$q75, hBforecast[[k]]$cfrq["75%", hBforecast[[k]]$coh.cfr])
  }
  COV90hB <- 100*mean((allobscfr >=  allhBcfrq$q05) & (allobscfr <= allhBcfrq$q95))
  COV50hB <- 100*mean((allobscfr >=  allhBcfrq$q25) & (allobscfr <= allhBcfrq$q75))
  
  #Freeze rates
  MAEfreeze  <- mean(AEfreeze, na.rm = T)
  MAPEfreeze <- mean(APEfreeze,  na.rm = T)*100
  RMSEfreeze <- sqrt(mean(AEfreeze^2, na.rm = T))
  
  #Print summary statistics
  sumstats <- round(data.frame(hB = c(MAEhB, MAPEhB, RMSEhB, COV90hB, COV50hB),
                   freeze = c(MAEfreeze, MAPEfreeze, RMSEfreeze, NA, NA)), 3)
  rownames(sumstats) <- c("MAE", "MAPE", "RMSE", "COV90", "COV50")
  cat("Summary statistics", "\n")
  print(sumstats)
  
  #Save summary statistics
  write.table(sumstats, paste(forecastFile, "/sumstats.txt", sep = ""))
}

## 3. Save age-specific forecast plots
plot.age <- function(age, k) {
  a <- which(ages == paste(age))
  plot(con, rep(0, 40), type = "n", ylim = c(0, 0.22), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2,
       xlab = "Cohort Year of Birth", ylab = expression(paste(theta)), main = bquote(paste("Age ", .(age))))
  polygon(c(con, rev(con)), c(hBforecast[[k]]$theta05[a,], rev(hBforecast[[k]]$theta95[a,])), col = rgb(121, 121, 255, alpha = 200, maxColorValue = 255), border = rgb(121, 121, 255, alpha = 200, maxColorValue = 255))
  polygon(c(con, rev(con)), c(hBforecast[[k]]$theta25[a,], rev(hBforecast[[k]]$theta75[a,])), col = rgb( 79,  79, 255, alpha = 200, maxColorValue = 255), border = rgb( 79,  79, 255, alpha = 200, maxColorValue = 255))
  lines(con, hBforecast[[k]]$theta05[a,], col = bpy_blue,lty=1)
  lines(con, hBforecast[[k]]$theta95[a,], col = bpy_blue,lty=1)
  lines(con, hBforecast[[k]]$theta25[a,], col = bpy_blue,lty=2)
  lines(con, hBforecast[[k]]$theta75[a,], col = bpy_blue,lty=2)
  points(con, hBforecast[[k]]$allrates[a,], pch = 16, cex = 1.2)
  abline(v = min(con) - 0.5 + sum(!is.na(hBforecast[[k]]$presrates[a,])), lty = 2)
  abline(h = exp(hBforecast[[k]]$Xbetaq)["50%", a])
}

for (k in name.red) {
  png(file=paste(forecastFile, "/", namef.red[k], "_age_plot.png", sep = ""), width = 1200, height = 400, family = "Calibri", res = 100)
  par(mfrow = c(1, 5), mar = c(4, 4, 2, 1))
  sapply(c(20, 25, 30, 35, 40), plot.age, k = k)
  dev.off()
}

## 4. Save CFR forecast plots (5 countries per plot)
plot.cohort <- function(k) {
  plot(con, rep(0, 40), type = "n", ylim = c(0.6, 3), cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2, las = 1,
       xlab = "Cohort Year of Birth", ylab = "CFR", main = bquote(paste(.(namef.red[k]))))
  polygon(c(con, rev(con)), c(hBforecast[[k]]$cfrq["5%",],  rev(hBforecast[[k]]$cfrq["95%",])), col = rgb(121, 121, 255, alpha = 200, maxColorValue = 255),border = rgb(121, 121, 255, alpha = 200, maxColorValue = 255))
  polygon(c(con, rev(con)), c(hBforecast[[k]]$cfrq["25%",], rev(hBforecast[[k]]$cfrq["75%",])), col = rgb( 79,  79, 255, alpha = 200, maxColorValue = 255),border = rgb( 79,  79, 255, alpha = 200, maxColorValue = 255))
  lines(con, hBforecast[[k]]$cfrq["5%",],  col = bpy_blue, lty = 1)
  lines(con, hBforecast[[k]]$cfrq["95%",], col = bpy_blue, lty = 1)
  lines(con, hBforecast[[k]]$cfrq["25%",], col = bpy_blue, lty = 2)
  lines(con, hBforecast[[k]]$cfrq["75%",], col = bpy_blue, lty = 2)
  points(con, apply(hBforecast[[k]]$allrates, 2, sum), pch = 16, cex = 1.2)
  abline(v = min(con) - 0.5 + hBforecast[[k]]$n.complete, lty = 2)
  abline(h = 2)
  if (k %in% name.score) { #add average scores to plot if country has at least 5 scores
    LogShBavgleg <- ifelse(LogShBavg[k] < 0, stringi::stri_pad_right(round(LogShBavg[k], 3), 6, "0"),
                           paste("+", stringi::stri_pad_right(round(LogShBavg[k], 3), 5, "0"), sep = ""))
    legend.text1 <- c(stringi::stri_pad_right(round(CRPShBavg[k], 3), 5, "0"), LogShBavgleg)
    legend("bottomleft",  paste(legend.text1[1]), pch = 20, col = "white", bg = "white", text.col = bpy_blue, horiz = T, title = "CRPS", title.col = 1, cex = 1.2)
    legend("bottomright", paste(legend.text1[2]), pch = 20, col = "white", bg = "white", text.col = bpy_blue, horiz = T, title = "LogS", title.col = 1, cex = 1.2)
  }
}

for (k in 1:ceiling(n.countries/5)) {
  png(file=paste(forecastFile, "/", "CFR_plot", k, ".png", sep = ""), width = 1200, height = 400, family = "Calibri", res = 100)
  par(mfrow = c(1, 5), mar = c(4, 4, 2, 1))
  sapply(name.red[(5*(k-1)+1):min(5*k, n.countries)], plot.cohort)
  dev.off()
}

## 5. Save rho plots
png(file=paste(forecastFile, "/", "rho_plot.png", sep = ""), width = 2000, height = 400, family = "Calibri", res = 100)
col.comm <- comm.ind + 1
par(mfrow=c(1, 4),mar=c(6, 4, 2, 1))
for (l in 1:4) {
  if (l == 1)
    plot(0, 0, xlim = c(1, n.countries), ylim = c(0.9, 1.1), type = "n", xaxt = "n", xlab = "", ylab = bquote(paste(rho[.(l)])))
  if (l %in% c(2:3))
    plot(0, 0, xlim = c(1, n.countries), ylim = c(0, 1),   type = "n", xaxt = "n", xlab = "", ylab = bquote(paste(rho[.(l)])))
  if (l == 4)
    plot(0, 0, xlim = c(1, n.countries), ylim = c(0.9, 1), type = "n", xaxt = "n", xlab = "", ylab = bquote(paste(rho[2]+rho[3])))
  axis(side = 1, at = which(comm.ind == 0), label = name.red[which(comm.ind == 0)], las = 2, col = 1)
  axis(side = 1, at = which(comm.ind == 1), label = name.red[which(comm.ind == 1)], las = 2, col.axis = 2)
  j <- 1
  for (k in name.red) {
    if (l %in% c(1:3)) {
      points(j, hBforecast[[k]]$rhoq["50%", l], pch = 19, col = col.comm[k])
      segments(j - 0.5, hBforecast[[k]]$rhoq["5%", l],  j + 0.5, hBforecast[[k]]$rhoq["5%", l],  col = col.comm[k])
      segments(j - 0.5, hBforecast[[k]]$rhoq["95%", l], j + 0.5, hBforecast[[k]]$rhoq["95%", l], col = col.comm[k])
      segments(j, hBforecast[[k]]$rhoq["5%", l], j, hBforecast[[k]]$rhoq["95%", l], col = col.comm[k])   
    }
    if (l==4) {
      points(j, hBforecast[[k]]$rhosumq["50%"], pch = 19, col = col.comm[k])
      segments(j - 0.5, hBforecast[[k]]$rhosumq["5%"], j + 0.5, hBforecast[[k]]$rhosumq["5%"],  col = col.comm[k])
      segments(j - 0.5, hBforecast[[k]]$rhosumq["95%"],j + 0.5, hBforecast[[k]]$rhosumq["95%"], col = col.comm[k])
      segments(j, hBforecast[[k]]$rhosumq["5%"], j, hBforecast[[k]]$rhosumq["95%"], col = col.comm[k])
    }
    j <- j + 1
  }
  if (l==4)
    abline(h=1,lty=2)
}
dev.off()

## 6. Print summary of computation times (in minutes)
cat("Summary of computation times (in minutes)", "\n")
print(summary(unlist(lapply(hBforecast[name.red], function(x) x$time[[3]]))/60))