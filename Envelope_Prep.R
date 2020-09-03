library(readxl)
setwd('C:/Users/Karl/Google Drive/Grad School_/Field 2020/Germination_assays')

packetsList = data.frame(read_excel(path = 'GGS_EntriesandPlots.xlsx'))

set.seed(1)
packetsList$Sort = runif(length(packetsList$PLOT),0,1)

packetsList = packetsList[order(packetsList$Sort),]

packetsList$SampleGroup = c(rep('a',1122/2),rep('b',1122/2))
 
packetsList$GroupOrder = paste(packetsList$SampleGroup,packetsList$PLOT)
packetsList = packetsList[order(packetsList$GroupOrder),]
packetsList$number = c(1:561,1:561)
packetsList$SampleID = paste(packetsList$SampleGroup,packetsList$number)

packetsListFull = rbind(packetsList,packetsList,packetsList,packetsList)
packetsListFull=packetsListFull[order(packetsListFull$GroupOrder),]
packetsListFull$Purpose = '30kernel'

write.csv(packetsListFull,file = 'PacketsToPrint30kernelAssays.csv')

DHList_Snyder= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='Snyder'))
DHList_Ketola= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='Ketola'))
DHList_Ketola_PHS= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='PHS_list_Ketola'))
DHList_Ketola<-merge(DHList_Ketola,DHList_Ketola_PHS,by=c('PLOT','Entry'),all=TRUE)
DHList_Ketola
DHList<-merge(DHList_Snyder,DHList_Ketola,by=c('Entry','Check','Source','Pedigree','Facultative_Spring_planted'),all = TRUE)
library(dplyr)
DHList$Thresh_sample.y

DHList_Check<- DHList[which(DHList$Check=="C" & DHList$),]
nrow(DHList_Check)
DHList_Sprout<- DHList[which(DHList$Sprout=="Y"),]

nrow(DHList_Check_Sprout)
View(DHList_Check_Sprout)
DHList_Other<-DHList[-which(DHList$Check=="C" | DHList$Sprout=="Y"),]  
nrow(DHList_Other)


typeof(DHList$Thresh_sample.x)

View(DHList[,c('Entry','Thresh_sample.x','Thresh_sample.y')])
by_cyl <- DHList %>% group_by(Pedigree)
c<-sample_frac(by_cyl,0.5)
table(c$Pedigree)
table(DHList$Pedigree)
# write.csv(packetsListMQ, file = 'Data/GerminationTests_2020/PrepForAssays/MaltingQualityPackets.csv')


######################
BagLocations = data.frame(read_excel('Data/GerminationTests_2020/PrepForAssays/PHSsubsample_PMdates_MassCollected.xlsx', sheet = 'MassCollected'))
bagLocationsWithSampleNumer = merge(x = BagLocations, y = packetsList[,c('PLOT','SampleID','GroupOrder')], by = 'PLOT')
write.csv(bagLocationsWithSampleNumer, file = 'Data/GerminationTests_2020/PrepForAssays/Sorting.csv')