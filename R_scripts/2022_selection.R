#script for some preliminary 2022 selections
load("data/Malting_quality/SWE/WMB21_CM.Rdata")
load("data/Malting_quality/SWE/WMB21_SWE.Rdata")
SWE
library(readxl)
WMB22_Reg_seneca <- as.data.frame(read_excel("C:/Users/kars8/Dropbox/Cornell Small Grains Breeding/2022/Winter Grains22/WMB Repl22/WMB22- Reg_kk.xlsx", 
                        sheet = "Seneca Co."))%>%select(Entry,PLOT,ws_percent)%>%arrange(PLOT)

WMB22_Reg_Snyder <- as.data.frame(read_excel("C:/Users/kars8/Dropbox/Cornell Small Grains Breeding/2022/Winter Grains22/WMB Repl22/WMB22- Reg_kk.xlsx", 
                                             sheet = "Snyder 2"))%>%select(Entry,PLOT,ws_percent)%>%arrange(PLOT)

WMB22_Reg_Helfer <- as.data.frame(read_excel("C:/Users/kars8/Dropbox/Cornell Small Grains Breeding/2022/Winter Grains22/WMB Repl22/WMB22- Reg_kk.xlsx", 
                                             sheet = "Helfer 1"))%>%select(Entry,PLOT,ws_percent)%>%arrange(PLOT)
WMB22_Reg_Avon<- as.data.frame(read_excel("C:/Users/kars8/Dropbox/Cornell Small Grains Breeding/2022/Winter Grains22/WMB Repl22/WMB22- Reg_kk.xlsx", 
                                             sheet = "Monroe Co."))%>%select(Entry,PLOT,ws_percent)%>%arrange(PLOT)
WMB22_Reg_seneca$location<-"Seneca"
WMB22_Reg_Snyder$location<-"Snyder"
WMB22_Reg_Helfer$location<-"Helfer"
WMB22_Reg_Avon$location<-"Livingston"
#add spatial corretion
WMB22_Reg_Snyder$Row<-rep(c(1,2,3),each=36)
WMB22_Reg_Snyder$Column<-c(1:36,36:1,1:36)
WMB22_Reg_Helfer$Row<-rep(c(1,2,3),each=36)
WMB22_Reg_Helfer$Column<-c(1:36,36:1,1:36)

WMB22_Reg_seneca$Row<-rep(c(1:6),each=18)
WMB22_Reg_seneca$Column<-rep(c(1:18,18:1),times=3)

WMB22_Reg_Avon$Row<-rep(c(1:6),each=18)
WMB22_Reg_Avon$Column<-rep(c(1:18,18:1),times=3)

total<-plyr::rbind.fill(WMB22_Reg_seneca,WMB22_Reg_Snyder,WMB22_Reg_Helfer,WMB22_Reg_Avon)
total<-total%>%arrange(location,Row,Column)
total$Row<-as.factor(total$Row)
total$Column<-as.factor(total$Column)
total$location<-as.factor(total$location)
total$Entry<-as.factor(total$Entry)
library(asreml)

ws<-asreml(fixed = ws_percent~Entry + location ,random=~Entry:location,residual =~dsum(~ar1(Row):ar1(Column)|location),data = total)
summary(ws)$varcomp
model<-as.data.frame(predict(ws,classify = "Entry"))
model
ws
ws_scores<-model%>%select(pvals.Entry,pvals.predicted.value)%>%dplyr::rename(Entry=1,ws=2)%>%arrange(desc(Entry))
View(ws_scores)
DH<-ws_scores %>%filter(grepl('BS', Entry))%>%filter(ws>87)
DH
CM$ME<-round(as.numeric(CM$ME),2)  
CM$ME

congress_mash<-CM%>%right_join(alat_select,by="Entry")#%>%filter(Entry%in%DH$Entry)%>%arrange(Entry)%>%select(Entry,PLOT,Treatment,ME,BG,FAN,sp_mean)
View(congress_mash)
alat_select
View(congress_mash)
SWE_lines<-SWE%>%filter(Entry%in%DH$Entry)%>%arrange(Entry)%>%select(Entry,PLOT,Treatment,AA,DP)
SWE_lines



load("data/Genotype_data/AlaT_markers.Rdata")
alat_select<-AlaT_markers%>%select(Entry,AlaT_Allele)
alat_select
total_s<-ws_scores%>%full_join(alat_select,by="Entry")

View(total)

t<-total%>%full_join(alat_select,by="Entry")
View(t)
