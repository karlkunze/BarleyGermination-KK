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
nrow(DHList_Snyder)
nrow(DHList_Ketola)
DHList_Ketola_Check<-DHList_Ketola[which(DHList_Ketola$Check=="C"),]
nrow(DHList_Ketola_Check)


DHList_Ketola<-DHList_Ketola[which(DHList_Ketola$Check=="B"),]
nrow(DHList_Ketola)
DHList_Snyder_Check<-DHList_Snyder[which(DHList_Snyder$Check=="C"),]
nrow(DHList_Snyder_Check)

DHList_Snyder<-DHList_Snyder[which(DHList_Snyder$Check=="B"),]
head(DHList_Snyder)

DHList_Snyder_Check$Check_ID<-substr(DHList_Snyder_Check$PLOT,3,5)
DHList_Ketola_Check$Check_ID<-substr(DHList_Ketola_Check$PLOT,3,5)

DHList_Check<-merge(DHList_Snyder_Check,DHList_Ketola_Check,by=c('Entry','Check','Source','Pedigree','Facultative_Spring_planted','Check_ID'),all = TRUE)

nrow(DHList_Snyder)
nrow(DHList_Ketola)
DHList<-merge(DHList_Snyder,DHList_Ketola,by=c('Entry','Check','Source','Pedigree','Facultative_Spring_planted'),all = TRUE)
write.csv(DHList, file="Dhlist_0.csv")


DHList_PHS<-DHList[which(DHList$Sprout=="Y"),]
nrow(DHList_PHS)
DHList<-DHList[-which(DHList$Sprout=="Y"),]

DHlist_common<-DHList[which(DHList$Thresh_sample.x=="Y"),]

nrow(DHList_Check)/(0.5*(nrow(DHList) +nrow(DHList_PHS)))
nrow(DHList_PHS)

0.5*nrow(DHList)
nrow(DHlist_common)
nrow(DHlist_common)+nrow(DHList_PHS)
250-nrow(DHList_PHS)-nrow(DHList_Check)

table(DHList$Pedigree,DHList$Facultative_Spring_planted)
table(DHlist_common$Pedigree,DHlist_common$Facultative_Spring_planted)
#120 from the common list need to be selected

c1<-nrow(DHList_PHS)
c2<-nrow(DHList_Check)



nrow(DHList)
nrow(DHlist_common)
value<-120/nrow(DHlist_common)
set.seed(1)
library(dplyr)
by_cyl <- DHlist_common %>% group_by(Pedigree)
Common_sample<-sample_frac(by_cyl,value)
nrow(Common_sample)

nrow(DHList_Check_Sprout)
View(DHList_Check_Sprout)
DHList_Other<-DHList[-which(DHList$Check=="C" | DHList$Sprout=="Y"),]  
nrow(DHList_Other)


typeof(DHList$Thresh_sample.x)

View(DHList[,c('Entry','Thresh_sample.x','Thresh_sample.y')])

c<-sample_frac(by_cyl,0.5)
table(c$Pedigree)
table(DHList$Pedigree)
# write.csv(packetsListMQ, file = 'Data/GerminationTests_2020/PrepForAssays/MaltingQualityPackets.csv')


######################
BagLocations = data.frame(read_excel('Data/GerminationTests_2020/PrepForAssays/PHSsubsample_PMdates_MassCollected.xlsx', sheet = 'MassCollected'))
bagLocationsWithSampleNumer = merge(x = BagLocations, y = packetsList[,c('PLOT','SampleID','GroupOrder')], by = 'PLOT')
write.csv(bagLocationsWithSampleNumer, file = 'Data/GerminationTests_2020/PrepForAssays/Sorting.csv')