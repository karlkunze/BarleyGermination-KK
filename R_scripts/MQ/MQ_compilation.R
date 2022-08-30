load("data/MQ/CM/WMB21_CM.Rdata")
load("data/MQ/SWE/WMB21_SWE.Rdata")
load("data/MQ/WMB21_Moistures.Rdata")
library(dplyr)
CM$Date
colnames(CM)[!colnames(CM)%in%colnames(SWE)]
SWE$PLOT[!SWE$PLOT%in%CM$PLOT]

Experiment_SWE<-SWE%>%filter(!Entry=="TMC")

Experiment_CM<-CM%>%filter(!Entry=="TMC")
Exp_Moisture<-WMB21_Moisture%>%rename(Treatment=Treatment_ID)%>%filter(!PLOT=="TMC")%>%arrange(Date,Order)%>%group_by(PLOT)%>%arrange(Order)%>%mutate(numbering=1,Date_Moisture=Date,Order_Moisture=Order)

Check_SWE<-SWE%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())
Check_SWE
Check_CM<-CM%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())

Check_Moisture<-WMB21_Moisture%>%rename(Treatment=Treatment_ID)%>%filter(PLOT=="TMC")%>%arrange(Date,Order)%>%group_by(Treatment,PLOT)%>%arrange(Order)%>%
  mutate(numbering=row_number())%>%rename(Entry=PLOT,Date_Moisture=Date,Order_Moisture=Order)


check<-full_join(Check_CM,Check_SWE,by=c("Entry","Treatment","numbering"))%>%
  full_join(Check_Moisture,by=c("Entry","Treatment","numbering"))%>%
  mutate(PLOT=paste0(PLOT.x,"-",PLOT.y))%>%select(!c(PLOT.x,PLOT.y))%>%filter(!Treatment=="SpringTP2-4")

Exper<-full_join(Experiment_CM,Experiment_SWE,by=c("Treatment","PLOT","Entry"))%>%full_join(Exp_Moisture,by=c("Treatment","PLOT"))%>%
  mutate(numbering=1)#%>%filter(!Treatment=="Rahr")
colnames(check)
colnames(Exper)
colnames(check)[!colnames(check)%in%colnames(Exper)]
colnames(Exper)[!colnames(Exper)%in%colnames(Exper)]
t<-Exper%>%plyr::rbind.fill(check)%>%arrange(Treatment,Date.x,Extract_number.x)
CTMC<-t%>%filter(Entry=="TMC")
CTMC$PLOT
t$E

