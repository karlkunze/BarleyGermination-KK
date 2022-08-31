load("data/MQ/CM/WMB21_CM.Rdata")
load("data/MQ/SWE/WMB21_SWE.Rdata")
load("data/MQ/WMB21_Moistures.Rdata")
load("data/MQ/Leco/Leco.Rdata")
library(dplyr)

#some last minute formatting, this change was to accomodate an out of order TMC check in TP1-1 of the the Congress mash
CM[CM$Date=="2022-01-31"&CM$Extract_number==52&CM$Entry=="TMC",c("Treatment","Extract_number")]<-c("WinterTP1-1",0.5)




Experiment_SWE<-SWE%>%filter(!Entry=="TMC")
Experiment_CM<-CM%>%filter(!Entry=="TMC")
Exp_Moisture<-WMB21_Moisture%>%rename(Treatment=Treatment_ID)%>%filter(!PLOT=="TMC")%>%arrange(Date,Order)%>%group_by(PLOT)%>%arrange(Order)%>%mutate(Date_Moisture=Date,Order_Moisture=Order)
Exp_Leco<-Leco%>%filter(!PLOT=="TMC")%>%arrange(Treatment,PLOT)



Check_SWE<-SWE%>%filter(!Treatment=="SpringTP2-4")%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())
Check_CM<-CM%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())%>%ungroup()
Check_Moisture<-WMB21_Moisture%>%rename(Treatment=Treatment_ID)%>%filter(PLOT=="TMC")%>%arrange(Date,Order)%>%group_by(Treatment,PLOT)%>%arrange(Order)%>%
  mutate(numbering=row_number())%>%rename(Entry=PLOT,Date_Moisture=Date,Order_Moisture=Order)
Check_Leco<-Leco%>%filter(PLOT=="TMC")%>%arrange(Treatment,TB)%>%group_by(Treatment)%>%arrange(Treatment,TB)%>%mutate(numbering=row_number())%>%ungroup()%>%select(!PLOT)
Check_Leco
#View(Check_Leco)
check<-full_join(Check_CM,Check_SWE,by=c("Entry","Treatment","numbering"))%>%
  full_join(Check_Moisture,by=c("Entry","Treatment","numbering"))%>%full_join(Check_Leco,by=c("Entry","Treatment","numbering"))%>%
  mutate(PLOT=paste0(PLOT.x,"-",PLOT.y))%>%select(!c(PLOT.x,PLOT.y))%>%filter(!Treatment=="SpringTP2-4")
#View(check)



Exper<-full_join(Experiment_CM,Experiment_SWE,by=c("Treatment","PLOT","Entry"))%>%full_join(Exp_Moisture,by=c("Treatment","PLOT"))%>%full_join(Exp_Leco,by=c("Treatment","PLOT","Entry"))%>%
  mutate(numbering=1)#%>%filter(!Treatment=="Rahr")

colnames(Exper)
colnames(check)[!colnames(check)%in%colnames(Exper)]
colnames(Exper)[!colnames(Exper)%in%colnames(Exper)]

MQ_WMB21<-Exper%>%plyr::rbind.fill(check)%>%arrange(Treatment,TB)


#
save(MQ_WMB21,file="data/MQ/WMB21_Master_Cornell.Rdata")#complete form

#Now were going to make a reduced form with calculated percentages


