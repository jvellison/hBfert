##Source file
## 1. Set-up
#Set working directory ("this.dir") to the location of these files on your computer
rm(list=ls()) # remove any objects currently in your workspace
this.dir <- c(".")
setwd(this.dir)

#Install required packages (only need to do this when running R code for first time)
source("install_packages.r")

#Download data from HFD (only need to do this when running R code for the first time)
replicate <- TRUE
#Change this to FALSE if you want to download the most recent HFD data rather than
#the data that was used to generate the results in the paper (downloaded on 15/2/19)
source("download_data.r")

#Load data
if (replicate)  load("data_replicate.RData")
if (!replicate) load("data_current.RData")

## 2. Generate Model hB forecasts
replicate2004 <- TRUE  #to replicate the 2004 forecasts
replicate2014 <- FALSE #change this to TRUE (and replicate2004 to FALSE) to replicate
                       #the 2014 forecasts

#Specify historical and contemporary cohorts, e.g., to replicate the 2004 forecasts
if (replicate2004) {
  his  <- 1904:1953     #historical cohorts (range must be at least 50 years)
  con  <- 1950:1989     #contemporary cohorts (range must be 40 years)
  pres <- max(con) + 15 #forecast year (when the youngest cohort is 15)
}

#And to replicate the 2014 forecasts
if (replicate2014) {
  his <- 1904:1953
  con <- 1960:1999
  pres <- max(con) + 15
}

#Source forecast file
source("forecast.r")
# This script creates a folder called "YYYY_MM_DD_HH.MM_Model_hB_pres_Forecast" in 
# "this.dir" where for each country it saves trace plots with and without warmup, and
# a CFR normal approximation plot. It also saves a file called "forecast.dput" which
# contains the forecast information for all countries stored in the list "hBforecast".
# NB: if a forecast has not been generated for a particular country k, look at
# hBforecast[[k]]$comment which contains an element called "cohorts" and possibly an
# element called "data". If either of these elements are "Insufficient (indicating that
# the first 39 cohorts of "con" are not present or there are not 10 or 11 complete cohorts
# respectively), then the country would have been skipped. If the "cohorts" element is
# "Sufficient" then the "presrates" element can be inspected to confirm the exclusion.

## 3. Process Model hB forecasts
#Specify the name of the folder ("forecastFile" in "forecast.r") containing the set of
#forecasts that you want to process
if (replicate2004) forecastFile <- "YYYY_MM_DD_HH.MM_Model_hB_2004_Forecast"
if (replicate2014) forecastFile <- "YYYY_MM_DD_HH.MM_Model_hB_2014_Forecast"

#Source process file
source("process.r")
# This script performs 6 tasks to process the chosen set of forecasts:
# 1. If there are at least 5 observed future CFR values for at least one country, the
#    average CRPS and LogS are computed and plotted against the respective countries,
#    with the average CRPS plotted for the freeze rates method as well (see Figures 5
#    and 6). These plots are saved as "CRPS_plot.png" and "LogS_plot.png" in "forecastFile".
# 2. If the condition for task 1 is satisfied, the summary statistics presented in Table 1
#    of the paper will be computed for Model hB and freeze rates, with the table printed
#    in the R session and also saved as "sumstats.txt" in "forecastFile".
# 3. Age-specific forecast plots like those in Figure 9 are saved for each country as
#    "Country_age_plot.png" in "forecastFile".
# 4. CFR forecast plots like those in Figures 7 and 10 are saved for each country (in
#    groups of 5) as "CFR_plotx.png" in "forecastFile". If the condition for task 1 is
#    satisfied then the average CRPS and LogS are added to the relevant plots.
# 5. Plots of the posterior distributions of rho1, rho2, rho3 and rho2+rho3 are saved
#    as "rho_plot.png" in "forecastFile" (see Figures 8 and 11).
# 6. A summary of the computation times (in minutes) is printed.
