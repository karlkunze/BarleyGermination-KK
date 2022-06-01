#File for moisture 
library(openxlsx)
library(readxl)



M_0130<-as.data.frame(read_xls("data/MQ/Moisture/21CY GGTB 1200-1264 21CY WBTP1 6011-6038.xls",skip = 7))%>%select(1:12)%>%
select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  filter(!Treatment_ID%in%c("SpringTP2-4"))%>%select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-01-30"))


M_0131<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP1 6039-6321.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-01-31"))

M_0202_01<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP1 6327-7123.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-02-02"))

M_0202_02<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP1 7129-7419.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-02-02"))

M_0207<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP1 WBTP2 7419-6184.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-02-07"))

M_0208<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP2 6188-7019.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-02-08"))
M_0301<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP2 7021-7301.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-03-01"))%>%mutate(Moisture_percent=as.numeric(Moisture_percent))
M_0302<-as.data.frame(read_xls("data/MQ/Moisture/21CY WBTP2 7305-7477.xls",skip = 7))%>%select(1:12)%>%
  select(Order,Dish_number,Sample_number,Dish_weight,Initial_sample_weight,Final_sample_weight,Moisture_percent,As_is_weight_170g_DB,Treatment_ID,Actual_sample_weight)%>%
  select(Order,Dish_number,Sample_number,Treatment_ID,Moisture_percent)%>%
  mutate(Date=as.Date("2022-03-02"))%>%mutate(Moisture_percent=as.numeric(Moisture_percent))

WMB21_Moisture<-plyr::rbind.fill(M_0130,M_0131,M_0202_01,M_0202_02,M_0207,M_0208,M_0301,M_0302)#
WMB21_Moisture<-WMB21_Moisture%>%mutate(Treatment_ID=gsub(" ", "", WMB21_Moisture$Treatment_ID, fixed = TRUE))%>%
  arrange(Treatment_ID,Date,Order,Sample_number)
WMB21_Moisture$Order<-1:nrow(WMB21_Moisture)
save(WMB21_Moisture,file="data/MQ/WMB21_Moistures.Rdata")
