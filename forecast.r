##Generate Model hB forecasts
start_time <- Sys.time()

## 1. Prepare data (adapted from the code provided in Schmertmann et al. (2014b))
#Set-up
name <- names(ASFR) #all countries
ages <- 15:44       #reproductive age range
if (diff(range(his)) < 49)  stop("Cohort range for historical data needs to be at least 50 years")
if (diff(range(con)) != 39) stop("Cohort range for contemporary data needs to be 40 years")
if ((min(con) - max(his)) < -4) warning("Overlap of more than 5 years between historical and contemporary cohort ranges")

#Create lists to contain important information for all countries
rates    <- list() #age x cohort ASFR matrices for each country
expos    <- list() #age x cohort Exposure matrices for each country
birth    <- list() #age x cohort implied births matrices for each country
cohort   <- list() #numeric vectors containing the available cohorts for each country
complete <- list() #logical vectors indicating which cohort schedules are complete for each country

#Fill these lists for each country
for (k in name) {
  sel           <- (ASFR[[k]]$ARDY %in% ages) # selected subset of ASFR (and Exposure) data
  rates[[k]]    <- tapply(ASFR[[k]]$ASFR[sel], list(ASFR[[k]]$ARDY[sel], ASFR[[k]]$Cohort[sel]), sum)
  cohort[[k]]   <- as.numeric(colnames(rates[[k]]))
  complete[[k]] <- apply(rates[[k]], 2, function(x) !any(is.na(x)))
  expos[[k]]    <- tapply(Exposure[[k]]$Exposure[sel], list(Exposure[[k]]$ARDY[sel], Exposure[[k]]$Cohort[sel]), sum)
  birth[[k]]    <- round(rates[[k]]*expos[[k]], 0)
}

#Create the vector "namef" which contains the full names of all the countries (this is correct as of 27/11/19, but you
#may need to add new countries following further updates to the HFD)
namef <- c(
  AUT         = "Austria",
  BLR         = "Belarus",
  BGR         = "Bulgaria",
  CAN         = "Canada",
  CHL         = "Chile",
  HRV         = "Croatia",
  CZE         = "Czechia",
  DNK         = "Denmark",
  EST         = "Estonia",
  FIN         = "Finland",
  FRATNP      = "France",
  DEUTNP      = "Germany",
  DEUTW       = "Western Germany",
  DEUTE       = "Eastern Germany",
  HUN         = "Hungary",
  ISL         = "Iceland",
  ITA         = "Italy",
  JPN         = "Japan",
  LTU         = "Lithuania",
  NLD         = "Netherlands",
  NOR         = "Norway",
  POL         = "Poland",
  PRT         = "Portugal",
  RUS         = "Russia",
  SVK         = "Slovakia",        
  SVN         = "Slovenia",
  ESP         = "Spain",
  SWE         = "Sweden",
  CHE         = "Switzerland",
  TWN         = "Taiwan",
  UKR         = "Ukraine",
  GBR_NP      = "Great Britain",
  GBRTENW     = "England and Wales",
  GBR_SCO     = "Scotland",
  GBR_NIR     = "Northern Ireland",
  USA         = "USA",
  KOR         = "Republic of Korea"
)[name]

#Create the vector "comm.ind" which contains a communist indicator for each country 
comm.ind <- c(
  AUT         = 0,
  BLR         = 1,
  BGR         = 1,
  CAN         = 0,
  CHL         = 0,
  HRV         = 1,
  CZE         = 1,
  DNK         = 0,
  EST         = 1,
  FIN         = 0,
  FRATNP      = 0,
  DEUTNP      = 0,
  DEUTW       = 0,
  DEUTE       = 1,
  HUN         = 1,
  ISL         = 0,
  ITA         = 0,
  JPN         = 0,
  LTU         = 1,
  NLD         = 0,
  NOR         = 0,
  POL         = 1,
  PRT         = 0,
  RUS         = 1,
  SVK         = 1,        
  SVN         = 1,
  ESP         = 0,
  SWE         = 0,
  CHE         = 0,
  TWN         = 0,
  UKR         = 1,
  GBR_NP      = 0,
  GBRTENW     = 0,
  GBR_SCO     = 0,
  GBR_NIR     = 0,
  USA         = 0,
  KOR         = 0
)[name]

## 2. Compute X, the matrix of the first three principal components obtained from the SVD of the matrix of
#     log historical cohort schedules (adapted from the code provided in Schmertmann et al. (2014b))
#Create PHI, the matrix of complete historical cohort schedules from all countries
PHI <- NULL
for (k in name) {
  keep <- complete[[k]] & (cohort[[k]] %in% his)
  tmp  <- rates[[k]][, keep]
  PHI  <- cbind(PHI, tmp)
}

#Create PI, the matrix of complete log historical cohort schedules from all countries
PI <- log(PHI)

#Create X
X <- svd(PI)$u[,1:3]

#Create Px to implement the identifiability constraint (see Section 2.2 of paper)
Px <- X %*% solve(t(X) %*% X) %*% t(X)

#Create "cov", which multiplies eta to give the vector of jump-off epsilons (see Section 2.2 of paper)
cov <- diag(30) - Px

## 3. Generate Model hB forecasts for all countries
#Load rstan package, set options and set fitting parameters
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
warmup <- 1000 #number of warmup iterations
iter   <- 9000 #total number of iterations including warmup
thin   <- 2    #factor to thin samples by
chains <- 1    #number of chains to run
cores  <- 1    #number of cores to use
if((((iter - warmup)/thin)/100)%%1 != 0) stop("Number of retained iterations needs to be divisible by 100")

#Create a folder called "forecastFile" within the current working directory where the output will be saved
forecastFile <- paste(format(Sys.time(), "%Y_%m_%d_%H.%M"), "_Model_hB_", pres, "_Forecast", sep = "")
dir.create(forecastFile)

#Create a list "hBforecast" that the forecasts will be saved in
hBforecast <- list()

#Obtain forecasts for each country
for (k in name) {
  #Save full name of country in "kf" and add to hBforecast
  kf <- (namef[k])[[1]]
  hBforecast[[k]]$namef <- kf
  
  #Print "kf" to know which country forecasts are currently being obtained for
  cat(kf, "\n")
  
  #Check the first 39 contemporary cohorts are present in rates, and add sufficiency comment (insufficient
  #cohorts if less than 39 are present) to hBforecast; skip this country if insufficient cohorts
  con.ind <- as.numeric(paste(con) %in% colnames(rates[[k]]))
  hBforecast[[k]]$comment$cohorts <- ifelse(sum(con.ind[1:39]) < 39, "Insufficient", "Sufficient")
  if(hBforecast[[k]]$comment$cohorts == "Insufficient") next
  
  #Subset data to contemporary cohorts in "allrates", "allexpos" and "allbirth" (add an extra column
  #of NAs if only the first 39 contemporary cohorts are present)
  n.con <- sum(con.ind)
  allrates  <- rates[[k]][, paste(con)[1:n.con]]
  allexpos  <- expos[[k]][, paste(con)[1:n.con]]
  allbirth  <- birth[[k]][, paste(con)[1:n.con]]
  if (n.con == 39) {
    allrates <- cbind(allrates, rep(NA, 30))
    allexpos <- cbind(allexpos, rep(NA, 30))
    allbirth <- cbind(allbirth, rep(NA, 30))
    colnames(allrates)[40] <- paste(con)[40]
    colnames(allexpos)[40] <- paste(con)[40]
    colnames(allbirth)[40] <- paste(con)[40]
  }
  
  #Make values observed after the present ("pres") NA in "presrates", "presexpos" and "presbirth"
  presrates <- allrates
  presexpos <- allexpos
  presbirth <- allbirth
  for (i in ages) {
    for (j in con) {
      if ((i+j) > pres) {
        presrates[(i - (min(ages) - 1)), (j - (min(con) - 1))] <- NA
        presbirth[(i - (min(ages) - 1)), (j - (min(con) - 1))] <- NA
        presexpos[(i - (min(ages) - 1)), (j - (min(con) - 1))] <- NA
      }
    }
  }
  
  #Backcast rates under 20, if any are missing for the first 5 cohorts (as in Schmertmann et al. (2014b))
  for (i in 1:5) {
    s <- sum(is.na(presrates[i, 1:5]))
    if (s > 0) {
      presrates[i, 1:s] <- presrates[i, (s + 1)]
      presexpos[i, 1:s] <- presexpos[i, (s + 1)]
      presbirth[i, 1:s] <- presbirth[i, (s + 1)]
    }}
  
  #Identify complete cohorts, compute the number of complete cohorts and add to hBforecast
  complete.ind <- as.numeric(apply(presrates, 2, function(x) !any(is.na(x))))
  n.complete <- sum(complete.ind)
  hBforecast[[k]]$n.complete <- n.complete

  #Add "allx" and "presx" matrices, and comm.ind to hBforecast
  hBforecast[[k]]$allrates  <- allrates
  hBforecast[[k]]$allexpos  <- allexpos
  hBforecast[[k]]$allbirth  <- allbirth
  hBforecast[[k]]$presrates <- presrates
  hBforecast[[k]]$presexpos <- presexpos
  hBforecast[[k]]$presbirth <- presbirth
  hBforecast[[k]]$comm.ind  <- as.numeric(comm.ind[k])
  
  #Add sufficiency comment (insufficient data if there are not 10 or 11 complete cohorts) to hBforecast;
  #skip this country if insufficient data
  hBforecast[[k]]$comment$data <- ifelse((n.complete != 10 & n.complete != 11), "Insufficient", "Sufficient")
  if (hBforecast[[k]]$comment$data == "Insufficient") next
  
  #Prepare data 
  N <- presbirth
  N[is.na(N)] <- 0
  W <- presexpos
  for (i in 2:30)
    W[i, (length(which(!is.na(W[i,])))+1):40] <- W[i, length(which(!is.na(W[i,])))]
  if (n.complete == 10)
    W[1, 40] <- W[1, 39]

  data <- list(X   = array(as.vector(X),dim = c(30, 3)),
               cov = array(as.vector(cov), dim = c(30, 30)),
               N   = array(as.vector(N),   dim = c(30, 40)),
               W   = array(as.vector(W),   dim = c(30, 40)))
  
  #Fit model
  if (n.complete == 10)
    fit <- stan(file="c10.stan", data = data, chains = chains, warmup = warmup, iter = iter, thin = thin, cores = cores)
  if (n.complete == 11)
    fit <- stan(file="c11.stan", data = data, chains = chains, warmup = warmup, iter = iter, thin = thin, cores = cores)
  
  #Add elapsed time to hBforecast
  hBforecast[[k]]$time <- c(get_elapsed_time(fit)[, 1:2], total = sum(get_elapsed_time(fit)))
  
  #Save trace plots with and without warmup
  windowsFonts(Calibri = windowsFont("Calibri"))
  png(file=paste(forecastFile, "/", kf, "_trace_plots_nowarmup.png", sep = ""),width=1500,height=1000,family="Calibri",res=100)
  par(mfrow=c(3, 4),mar=c(4, 4, 2, 1))
  plot(extract(fit, "beta[1]", permuted = FALSE), type = "l", ylab =  expression(paste(beta[1])), main =  expression(paste(beta[1])), cex.main = 1.2)
  plot(extract(fit, "beta[2]", permuted = FALSE), type = "l", ylab =  expression(paste(beta[2])), main =  expression(paste(beta[2])), cex.main = 1.2)
  plot(extract(fit, "beta[3]", permuted = FALSE), type = "l", ylab =  expression(paste(beta[3])), main =  expression(paste(beta[3])), cex.main = 1.2)
  par(mfg=c(2, 1))
  plot(extract(fit,"sigma[1]", permuted = FALSE), type = "l", ylab = expression(paste(sigma[1])), main = expression(paste(sigma[1])), cex.main = 1.2)
  plot(extract(fit,"sigma[2]", permuted = FALSE), type = "l", ylab = expression(paste(sigma[2])), main = expression(paste(sigma[2])), cex.main = 1.2)
  plot(extract(fit,"sigma[3]", permuted = FALSE), type = "l", ylab = expression(paste(sigma[3])), main = expression(paste(sigma[3])), cex.main = 1.2)
  plot(extract(fit,"sigma[4]", permuted = FALSE), type = "l", ylab = expression(paste(sigma[4])), main = expression(paste(sigma[4])), cex.main = 1.2)
  plot(extract(fit,  "rho[1]", permuted = FALSE), type = "l", ylab =   expression(paste(rho[1])), main =   expression(paste(rho[1])), cex.main = 1.2)
  plot(extract(fit,  "rho[2]", permuted = FALSE), type = "l", ylab =   expression(paste(rho[2])), main =   expression(paste(rho[2])), cex.main = 1.2)
  plot(extract(fit,  "rho[3]", permuted = FALSE), type = "l", ylab =   expression(paste(rho[3])), main =   expression(paste(rho[3])), cex.main = 1.2)
  dev.off()
  png(file=paste(forecastFile, "/", kf, "_trace_plots_warmup.png", sep = ""),width=1500,height=1000,family="Calibri",res=100)
  par(mfrow=c(3, 4),mar=c(4, 4, 2, 1))
  plot(extract(fit, "beta[1]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab =  expression(paste(beta[1])), main =  expression(paste(beta[1])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit, "beta[2]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab =  expression(paste(beta[2])), main =  expression(paste(beta[2])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit, "beta[3]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab =  expression(paste(beta[3])), main =  expression(paste(beta[3])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  par(mfg=c(2, 1))
  plot(extract(fit,"sigma[1]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab = expression(paste(sigma[1])), main = expression(paste(sigma[1])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit,"sigma[2]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab = expression(paste(sigma[2])), main = expression(paste(sigma[2])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit,"sigma[3]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab = expression(paste(sigma[3])), main = expression(paste(sigma[3])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit,"sigma[4]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab = expression(paste(sigma[4])), main = expression(paste(sigma[4])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit,  "rho[1]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab =   expression(paste(rho[1])), main =   expression(paste(rho[1])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit,  "rho[2]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab =   expression(paste(rho[2])), main =   expression(paste(rho[2])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  plot(extract(fit,  "rho[3]", permuted = FALSE, inc_warmup = TRUE), type = "l", ylab =   expression(paste(rho[3])), main =   expression(paste(rho[3])), cex.main = 1.2); abline(v = warmup/thin, lty = 2)
  dev.off()
  
  #Extract samples
  mu      <- extract(fit,      "mu")$mu
  beta    <- extract(fit,    "beta")$beta
  rho     <- extract(fit,     "rho")$rho
  epsilon <- extract(fit, "epsilon")$epsilon
  sigma   <- extract(fit,   "sigma")$sigma
  S       <- dim(mu)[1] #number of samples    
  
  #Add additional Poisson variation (see Section 2.3)
  set.seed(1)
  mu_r <- rpois(S*30*40, mu)

  #Create matrices to store theta quantiles
  theta_05 <- matrix(0, 30, 40)
  theta_95 <- matrix(0, 30, 40)
  theta_25 <- matrix(0, 30, 40)
  theta_75 <- matrix(0, 30, 40)
  theta_50 <- matrix(0, 30, 40)
  
  #Extract theta quantiles
  for (i in 1:1200) {
    a <- sort(mu_r[(((i-1)*S) + 1):(i*S)])
    theta_05[i] <- a[S*0.05]/W[i]
    theta_95[i] <- a[S*0.95]/W[i]
    theta_25[i] <- a[S*0.25]/W[i]
    theta_75[i] <- a[S*0.75]/W[i]
    theta_50[i] <- a[S*0.50]/W[i]
  }
  
  #Add parameter quantiles to hBforecast
  hBforecast[[k]]$theta05 <- theta_05
  hBforecast[[k]]$theta95 <- theta_95
  hBforecast[[k]]$theta25 <- theta_25
  hBforecast[[k]]$theta75 <- theta_75
  hBforecast[[k]]$theta50 <- theta_50
  
  hBforecast[[k]]$betaq   <- apply(beta,  2, quantile, p = c(0.05, 0.95, 0.25, 0.75, 0.5), type = 1)
  hBforecast[[k]]$sigmaq  <- apply(sigma, 2, quantile, p = c(0.05, 0.95, 0.25, 0.75, 0.5), type = 1)
  hBforecast[[k]]$rhoq    <- apply(rho,   2, quantile, p = c(0.05, 0.95, 0.25, 0.75, 0.5), type = 1)
  hBforecast[[k]]$rhosumq <- quantile(apply(rho[,2:3], 1, sum), p = c(0.05, 0.95, 0.25, 0.75, 0.5), type = 1)
 
  #Create "surface", a list containing the S sampled complete contemporary Lexis surfaces
  surface <- list()
  for (i in 1:S)
    surface[[i]] <- matrix(mu_r[seq(from = i, to = (S*((30*40) - 1) + i), length.out = 30*40)], nrow = 30, ncol = 40, byrow = F)/W
  
  #Extract and sort the S CFR samples for each of the 40 cohorts
  cfr <- matrix(0, S, 40)
  colnames(cfr) <- con
  for (j in 1:S) {
    cfr[j,] <- apply(surface[[j]], 2, sum)
  }
  for (j in 1:40) {
    cfr[,j] <- sort(cfr[,j])
  }
  
  #Add CFR quantiles, mean and standard deviation to hBforecast
  hBforecast[[k]]$cfrq    <- apply(cfr, 2, quantile, p = c(0.05, 0.95, 0.25, 0.75, 0.5), type = 1)
  hBforecast[[k]]$cfrmean <- apply(cfr, 2, mean)
  hBforecast[[k]]$cfrsd   <- apply(cfr, 2, sd)
  
  #Extract Xbeta samples
  Xbeta_mat <- matrix(0,S,30)
  for (j in 1:S) {
    Xbeta_mat[j,] <- X%*%beta[j,]
  }
  
  #Add Xbeta quantiles to hBforecast
  hBforecast[[k]]$Xbetaq <- apply(Xbeta_mat, 2, quantile, p = c(0.05, 0.95, 0.25, 0.75, 0.5), type = 1)
  
  #Calculate CRPS using the empirical CDF-based approximation given in equation (3) of Jordan et al. (2017)
  obs.cfr     <- apply(allrates[, paste((min(con) + n.complete):max(con))], 2, sum) #observed future CFR values (after pres)
  n.cfr       <- sum(!is.na(obs.cfr)) #number of observed future CFR values
  CRPS        <- rep(NA, 40) #vector to store CRPS values
  names(CRPS) <- con
  if (n.cfr > 0) { #compute CFR values for cohorts with observed future CFR values
    coh.cfr     <- names(obs.cfr[1:n.cfr]) #cohorts with observed future CFR values
    for (h in coh.cfr)
      CRPS[h]   <- (2/(S^2))*sum((cfr[,h] - obs.cfr[h])*(S*as.numeric(rep(obs.cfr[h], S) < cfr[,h]) - (1:S) + 0.5))
  }

  #Add obs.cfr, n.cfr, coh.cfr and CRPS values to hBforecast
  hBforecast[[k]]$obs.cfr <- ifelse(n.cfr > 0, obs.cfr[1:n.cfr], NA)
  hBforecast[[k]]$n.cfr   <- n.cfr
  hBforecast[[k]]$coh.cfr <- ifelse(n.cfr > 0, coh.cfr, NA)
  hBforecast[[k]]$CRPS    <- CRPS
  
  #Save plot of normal approximations to CFR posterior distribution to check suitability for calculating the LogS
  png(file=paste(forecastFile, "/", kf, "_normal_approx.png", sep = ""), width = 1500, height = 1000, family = "Calibri", res = 100)
  par(mfrow=c(6, 7), mar = c(4, 4, 2, 1))
  for (i in 1:40) {
    hist(cfr[,i], breaks = 30, freq = F, main = paste(con[i], "Cohort"), xlab = "CFR")
    curve(dnorm(x, hBforecast[[k]]$cfrmean[i], hBforecast[[k]]$cfrsd[i]), add = T)
  }
  dev.off()
  
  #Save hBforecast
  dput(hBforecast, file = paste(forecastFile, "/", "forecast.dput", sep = ""))
}

end_time <- Sys.time()
print(end_time - start_time) #duration