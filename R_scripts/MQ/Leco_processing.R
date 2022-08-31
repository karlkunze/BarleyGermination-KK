#Leco compilation
library(readxl)
Leco <- read_excel("data/MQ/Leco/21CY Cornell Winter TB Malting Leco Final Report.xlsx", 
                                                                skip = 2)%>%rename(TB=1,NitrogenPercent=5,Date_leco=6)
str(Leco)
save(Leco,file="data/MQ/Leco/Leco.Rdata")
