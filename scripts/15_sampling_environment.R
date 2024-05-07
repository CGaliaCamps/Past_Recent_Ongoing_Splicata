
library(ncdf4)
library(RANN)
library(reshape)
library(dplyr)
library(tidyverse)

####This is an example for mean salinity

wd <- "G:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/ambientals/"
  
setwd(wd)

nc <- list.files(wd, pattern = ".nc")

    ncin  <- nc_open(paste(nc[1], sep = ""))
  

    lt <- ncvar_get(ncin, "latitude")
    ln <- ncvar_get(ncin, "longitude")
    time <- ncvar_get(ncin, "time")
    depth <- ncvar_get(ncin, "depth")
   
    ##Aprox. projected 10km ->
    latdist <- 30*(1/111.325)
    londist <- 30*(1/85.27990)
    ddist <- (latdist+londist)/2
    
    #Read sampling points 
    ptt <- read.table("points.txt", sep = "\t", header = T)
    
    #Extract salinity from the file nc
  
    sal <- ncvar_get(ncin, "to")
     
    print(dim(sal))
    print(dim(depth))
    print(dim(time))
    
    #for each locality, extract the mean values 
    sal <- apply(sal, c(1, 2), max, na.rm = TRUE)
    
    rownames(sal) <- ln
    colnames(sal) <- lt
    
    print(dim(sal))
    print(dim(depth))
    print(dim(time))
  
    #Transform to table tabla and name columns
    sal <- melt(sal)
    colnames(sal)<-c("lon", "lat", "to")
    
    #make sure both files are equally sorted (lon/lat) 
    head(sal)
    head(ptt)
    
    sal <- as.data.frame(sal)
    ptt <- as.data.frame(ptt)
    #Look for the moints below "ddist" (see above) distance from sampling points (in degrees)
    dt <- data.frame()
    for (i in 1:length(ptt$lon)){
      dmed1i <- RANN::nn2(ptt[i,1:2], na.omit(sal[,1:2]),  k = 1, searchtype = "radius", radius = ddist)
      dmed2i <- cbind(sal, dmed1i)
      dmed3i <- dmed2i[dmed2i$nn.idx!=0,]
      dmed4i <- cbind(ptt[i,], max(dmed3i$to, na.rm = T))
      
      dt <- rbind(dt,dmed4i)
    }
  
    colnames(dt) <- c("lon", "lat", "pop", "maxto")

    dt
#Export table
write.table(dt, "maxto.txt", sep = "\t", row.names = F)
