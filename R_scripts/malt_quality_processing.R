#file to read in data for malting quality
#need to be aware of positioning on the excel file, as that is how we have info on TMC checks
library(openxlsx)
library(dplyr)

t2<-as.data.frame(read.xlsx("data/Malting_quality/SWE/SWE 01-28-22 21CYGGS1216-endspringtp2_WMB6011-6051-as.xlsx",
                           sheet = "Data",startRow = 9,detectDates = TRUE))%>%rename(Set=1,Extract_number=2,Entry=3,Treatment=4,PLOT=5,Date=6,DP=7,AA=8)%>%select(1:8)


SWE<-function(t){
  
total<-as.data.frame(table(t$Date))[1,2]

t[t$PLOT=="Tradition Malt Check"&!is.na(t$PLOT),c("Entry")]<-"TMC"
total

if (total==94) {
  t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t[49:72,]$Set<-"Set 3";t[73:96,]$Set<-"Set 4";t$Date<-t$Date[3]
  t<-t%>%slice(1:96)
  return(t)
} 

else {
t[1:24,]$Set<-"Set 1"; t[25:48,]$Set<-"Set 2";t[49:72,]$Set<-"Set 3";t$Date<-t$Date[3];t<-t%>%slice(1:72)
return(t)
}

t$Check<-"E";t[t$Entry%in%c("R","TMC"),]$Check<-"C";t$Check<-as.factor(t$Check)

#specific to split spring and winter
return(t)

}
SWE(t)
t[1:56,]$Treatment<-"SpringTP2-4";t[57:72,]$Treatment<-"WinterTP1-1"

t<-t%>%filter(Treatment=="WinterTP1-1"|Check=="C")%>%arrange(Extract_number)

