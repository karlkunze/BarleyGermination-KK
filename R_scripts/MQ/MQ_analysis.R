#analysis of MQ data

load("data/MQ/WMB21_Moistures.Rdata")


library(readxl)
library(dplyr)
library(stringr)
WMB21_leco<- read_excel("data/MQ/21CY Cornell Winter TB Malting Leco Final Report.xlsx",skip=2)%>%as.data.frame()%>%
  select(-Treatment)%>%rename(TB=1,Treatment="Entry ID",N_percent="Nitrogen %",Date="Analysis Date",Year="Crop Year")%>%
  mutate(Date=as.Date(Date))


WMB21_leco[WMB21_leco$Treatment%in%c("TMC"),]$Entry<-"TMC"
WMB21_leco[WMB21_leco$Treatment%in%c("TMC"),]$PLOT<-"TMC"
WMB21_leco[1:64,]$Treatment<-"WinterTP1-1"
WMB21_leco[65:128,]$Treatment<-"WinterTP1-2"
WMB21_leco[129:192,]$Treatment<-"WinterTP1-3"
WMB21_leco[193:256,]$Treatment<-"WinterTP1-4"
WMB21_leco[257:320,]$Treatment<-"WinterTP2-1"
WMB21_leco[321:384,]$Treatment<-"WinterTP2-2"
WMB21_leco[385:448,]$Treatment<-"WinterTP2-3"
WMB21_leco[449:512,]$Treatment<-"WinterTP2-4"
WMB21_leco<-WMB21_leco%>%mutate(Order=1:nrow(WMB21_leco))#
TMC_sets<-WMB21_leco%>%filter(Entry%in%c("TMC"))%>%group_by(Treatment)%>%mutate(n=1:n())%>%ungroup()%>%
  select(Order,Entry,Treatment,Date,n)
WMB21_leco<-WMB21_leco%>%left_join(TMC_sets,by=c("Order","Treatment","Date","Entry"))

WMB21_leco$TP<-substr(WMB21_leco$Treatment, start = 9, stop = 11)
WMB21_leco$PLOT_TP<-paste0(WMB21_leco$PLOT,"-",WMB21_leco$TP)
WMB21_leco[WMB21_leco$Entry%in%c("TMC"),]$PLOT_TP<-paste0(WMB21_leco[WMB21_leco$Entry%in%c("TMC"),]$PLOT_TP,"-",WMB21_leco[WMB21_leco$Entry%in%c("TMC"),]$n)

WMB21_M<-WMB21_Moisture%>%rename(Treatment=Treatment_ID,Date_Moisture=Date)%>%filter(!is.na(PLOT))%>%group_by(Treatment)

TMC_sets<-WMB21_M%>%filter(PLOT%in%c("TMC"))%>%group_by(Treatment)%>%mutate(n=1:n())%>%ungroup()%>%
  select(Order,PLOT,Treatment,n)
WMB21_M<-WMB21_M%>%left_join(TMC_sets,by=c("Order","Treatment","PLOT"))

WMB21_M$TP<-substr(WMB21_M$Treatment, start = 9, stop = 11)
WMB21_M$PLOT_TP<-paste0(WMB21_M$PLOT,"-",WMB21_M$TP)
WMB21_M[WMB21_M$PLOT%in%c("TMC"),]$PLOT_TP<-paste0(WMB21_M[WMB21_M$PLOT%in%c("TMC"),]$PLOT_TP,"-",WMB21_M[WMB21_M$PLOT%in%c("TMC"),]$n)

WMB21_Leco_M<-full_join(WMB21_leco,WMB21_M,by=c("PLOT_TP","Treatment","TP","PLOT"))%>%select(PLOT,Entry,Treatment,PLOT_TP,N_percent,Moisture_percent)
###

load("data/MQ/SWE/WMB21_SWE.Rdata")
SWE<-SWE%>%filter(!Treatment%in%c("SpringTP2-4"))
SWE[SWE$Entry%in%c("TMC"),]$PLOT<-SWE[SWE$Entry%in%c("TMC"),]$Entry

TMC_sets<-SWE%>%filter(Entry%in%c("TMC"))%>%group_by(Treatment)%>%mutate(n=1:n())%>%ungroup()%>%
  select(order_column,Entry,Treatment,Date,n)
table(TMC_sets$n,TMC_sets$Treatment)

SWE<-SWE%>%left_join(TMC_sets,by=c("order_column","Treatment","Date","Entry"))

SWE$TP<-substr(SWE$Treatment, start = 9, stop = 11)
SWE[SWE$Entry%in%c("R"),]$TP<-"R"

SWE$PLOT_TP<-paste0(SWE$PLOT,"-",SWE$TP)
SWE[SWE$PLOT%in%c("TMC"),]$PLOT_TP<-paste0(SWE[SWE$PLOT%in%c("TMC"),]$PLOT_TP,"-",SWE[SWE$PLOT%in%c("TMC"),]$n)
SWE$PLOT_TP
table(SWE$Treatment)


load("data/MQ/CM/WMB21_CM.Rdata")
CM<-CM%>%arrange(Treatment,Date,Extract_number)
View(CM)
CM[CM$Entry%in%c("TMC"),]$PLOT<-CM[CM$Entry%in%c("TMC"),]$Entry

MQ_list<-as.data.frame(read.xlsx("data/Phenotype_Data/2021/MQ_samples_schedule.xlsx",
                                 sheet = "WinterSamples"))%>%filter(rep==1)
#View(MQ_list)
CM$order_column<-1:nrow(CM)
CM[CM$Entry%in%c("TMC"),]$PLOT<-CM[CM$Entry%in%c("TMC"),]$Entry
#View(CM)
#View(CM)
TMC_sets<-CM%>%filter(PLOT%in%c("TMC"))%>%group_by(Treatment)%>%mutate(n=1:n())%>%ungroup()%>%
  select(order_column,Entry,Treatment,Date,n)
table(TMC_sets$n,TMC_sets$Treatment)


CM<-CM%>%left_join(TMC_sets,by=c("order_column","Treatment","Date","Entry"))

CM$TP<-substr(CM$Treatment, start = 9, stop = 11)

CM[CM$Entry%in%c("NaCl"),]$TP<-"N"

CM$PLOT_TP<-paste0(CM$PLOT,"-",CM$TP)
CM[CM$PLOT%in%c("TMC"),]$PLOT_TP<-paste0(CM[CM$PLOT%in%c("TMC"),]$PLOT_TP,"-",CM[CM$PLOT%in%c("TMC"),]$n)
CM$PLOT_TP
t<-CM%>%full_join(SWE,by=c("PLOT_TP","Entry","PLOT"))
View(t)
