# CSAA-IA experiment: Wall et al. 
# using M. capitata
# PAR Odyssey

rm(list=ls())
ls()

library(reshape2)
library(plyr)
library(ggplot2)
library(lubridate)
library(zoo)
library(plotrix)
library(devtools)
library(tools)

main<-setwd(getwd())
setwd("~/Desktop/Research and Teaching/UH MANOA/Research/CSAA-IA")


list.files(path="Nov tanks_light", pattern = "CSV$", full.names = T)

#####################
# LiCor
# 4 pi sensor
#####################
data<-read.csv("Wall_indoor_licor4pi.csv", skip=13)
df<-data[, c(2:3)] # remove trash columns
colnames(df)<-c("timestamp", "Licor.RAW")
df$timestamp<-mdy_hm(as.character(df$timestamp)) #  timestamp as character to POSIXct date
df$timestamp<-strptime(df$timestamp, format="%Y-%m-%d %H:%M:%S")

#remove all periods where loggers were not logging with the LiCor
df<-df[!(df$timestamp < "2017-12-15 06:15:00"),] # start at this time
df<-df[!(df$timestamp > "2017-12-15 17:45:00"),]

df$Licor.RAW<-as.numeric(as.character(df$Licor.RAW)) 
df$PAR<-(df$Licor.RAW*(10^6)/(15*60))

str(df)

Licor.data<-df[,c(1,3)] # cleaned Licor data with night removed (night is not present due to logging program)

#### FIGURE
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(y=df$PAR, x=df$timestamp, type="o",
     cex=0.5, col="black",
     ylab=expression(paste(mu*mol~photons~m^-2~s^1), sep=""), ylim=c(0,500), 
     xlab = "Date")

dev.copy(pdf, "Tank6_PAR.pdf", encod="MacRoman", height=4, width=7)
dev.off()
####

# umol/s.. /1,000,000 * 12h * 60m * 60s = 0.0432...
# gives the average of mol photons/d

mean((df$PAR)*0.0432) ## ~12 mol photons/d
