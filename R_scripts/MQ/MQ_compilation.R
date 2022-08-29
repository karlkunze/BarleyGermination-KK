load("data/phenotypes/MQ_2021/WMB21_CM.Rdata")
load("data/phenotypes/MQ_2021/WMB21_Moistures.Rdata")
load("data/phenotypes/MQ_2021/WMB21_SWE.Rdata")
table(CM$Treatment)
colnames(CM)[!colnames(CM)%in%colnames(SWE)]
table(SWE$Check)
Experiment_SWE<-SWE%>%filter(!Entry=="TMC")

Experiment_CM<-CM%>%filter(!Entry=="TMC")

Check_SWE<-SWE%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())
Check_SWE
Check_CM<-CM%>%filter(Entry=="TMC")%>%arrange(Entry,Date,Extract_number)%>%group_by(Treatment)%>%arrange(Date)%>%mutate(numbering = row_number())
Check_SWE
check<-full_join(Check_CM,Check_SWE,by=c("Entry","Treatment","numbering"))%>%mutate(PLOT=paste0(PLOT.x,"-",PLOT.y))%>%select(!c(PLOT.x,PLOT.y))
ncol(check)
t<-full_join(Experiment_CM,Experiment_SWE,by=c("Treatment","PLOT","Entry"))%>%
  mutate(numbering=1)%>%rbind(check)%>%arrange(Treatment,PLOT)

View(t)