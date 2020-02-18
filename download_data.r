##Download HFD data
library(HMDHFDplus)

myusername <- readline(prompt = "Enter HFD username:") # provide your HFD username
mypassword <- readline(prompt = "Enter HFD password:") # provide your HFD password

if (replicate==TRUE) {
  #Create the vector "name" which contains the country codes
  name <- c("AUT", "BLR", "BGR", "CAN", "CHL", "HRV", "CZE", "DNK", "EST", "FIN",
            "FRATNP", "DEUTNP", "DEUTW", "DEUTE", "HUN", "ISL", "ITA", "JPN", "LTU",
            "NLD", "NOR", "POL", "PRT", "RUS", "SVK", "SVN", "ESP", "SWE", "CHE",
            "TWN", "UKR", "GBR_NP", "GBRTENW", "GBR_SCO", "GBR_NIR", "USA")
  
  #Create the vector "Date" which contains the most recent update as of 15/2/2019
  Date <- c("20181012", "20180525", "20111101", "20160311", "20140519", "20181130",
            "20180416", "20180906", "20181121", "20170220", "20180924", "20160620",
            "20160620", "20160620", "20181121", "20181130", "20170918", "20181130",
            "20181121", "20180815", "20160301", "20171006", "20170327", "20160608",
            "20180525", "20160208", "20181121", "20180131", "20180815", "20160511",
            "20160118", "20180906", "20180906", "20180906", "20180906", "20180815")
  names(Date) <- name
  
  #Create the list "ASFR" and fill it with the HFD fertility rates for each country
  ASFR <- list()
  for (k in name)
    ASFR[[k]] <- readHFDweb(CNTRY = k, item="asfrVV",
                             username=myusername, password=mypassword, Update=Date[k])
  
  #Create the list "Exposure" and fill it with the HFD exposures for each country
  Exposure <- list()
  for (k in name)
    Exposure[[k]] <- readHFDweb(CNTRY = k, item="exposVV",
                             username=myusername, password=mypassword, Update=Date[k])
  
  #Save ASFR and Exposure in the file "data_replicate.RData"
  save(ASFR, Exposure, file="data_replicate.RData")
  rm(myusername, mypassword, name, Date, ASFR, Exposure)
}

if (replicate==FALSE) {
  #Create the vector "name" which contains the currently available HFD country codes
  name <- getHFDcountries()
  
  #Create the list "ASFR" and fill it with the HFD fertility rates for each country
  ASFR <- list()
  for (k in name)
    ASFR[[k]] <- readHFDweb(CNTRY = k, item="asfrVV",
                             username=myusername, password=mypassword)
  
  #Create the list "Exposure" and fill it with the HFD exposures for each country
  Exposure <- list()
  for (k in name)
    Exposure[[k]] <- readHFDweb(CNTRY = k, item="exposVV",
                             username=myusername, password=mypassword)
  
  #Save ASFR and Exposure in the file "data_current.RData"
  save(ASFR, Exposure, file="data_current.RData")
  rm(myusername, mypassword, name, ASFR, Exposure)
}