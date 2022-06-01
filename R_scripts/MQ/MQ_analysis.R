#analysis of MQ data

load("data/MQ/WMB21_Moistures.Rdata")
load("data/MQ/CM/WMB21_CM.Rdata")
load("data/MQ/SWE/WMB21_SWE.Rdata")
SWE
library(readxl)
library(dplyr)
WMB21_leco<- read_excel("data/MQ/21CY Cornell Winter TB Malting Leco Final Report.xlsx",skip=2)%>%as.data.frame()%>%rename()

names(WMB21_leco)
