library(readxl)








DHList_Snyder= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='Snyder'))
DHList_Ketola= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='Ketola'))

DHList_Ketola_PHS= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='PHS_list_Ketola'))
DHList_Ketola_PHS<-na.omit(DHList_Ketola_PHS)
DHList_Ketola_PHS

DHList_Ketola<-merge(DHList_Ketola,DHList_Ketola_PHS,by=c('PLOT','Entry'),all=TRUE)
View(DHList_Ketola)
nrow(DHList_Snyder)
nrow(DHList_Ketola)
which(is.na(DHList_Ketola$Location))

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
nrow(DHList_Check)
nrow(DHList_Snyder)
nrow(DHList_Ketola)
DHList<-merge(DHList_Snyder,DHList_Ketola,by=c('Entry','Check','Source','Pedigree','Facultative_Spring_planted'),all = TRUE)
write.csv(DHList, file="Dhlist_0.csv")


DHList_PHS<-DHList[which(DHList$Sprout=="Y"),]
nrow(DHList_PHS)
DHList<-DHList[-which(DHList$Sprout=="Y"),]
nrow(DHList)
DHlist_common<-DHList[which(DHList$Thresh_sample.x=="Y" & DHList$Thresh_sample.y=="Y"),]
nrow(DHlist_common)
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
value<-120.5/nrow(DHlist_common)
value
set.seed(1)
library(dplyr)
by_cyl <- DHlist_common %>% group_by(Pedigree)
Common_sample<-as.data.frame(sample_frac(by_cyl,value))
nrow(Common_sample)

nrow(Common_sample)+nrow(DHList_PHS)+nrow(DHList_Check)
str(Common_sample)
str(DHList_Check)
str(DHList_PHS)

Common_sample$Check_ID<-substring(Common_sample$PLOT.y,3,5)
DHList_PHS$Check_ID<-substring(DHList_PHS$PLOT.y,3,5)

Full_samples<-rbind(Common_sample,DHList_PHS,DHList_Check)
Full_samples$Full_sample<-"Full"

DHList_Snyder2= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='Snyder'))
DHList_Ketola2= data.frame(read_excel(path = 'DH_2020_germination_assay_list.xlsx',sheet='Ketola'))
DHList_Ketola2<-merge(DHList_Ketola2,DHList_Ketola_PHS,by=c('PLOT','Entry'),all=TRUE)

DHList_Snyder2$Check_ID<-substr(DHList_Snyder2$PLOT,3,5)
DHList_Ketola2$Check_ID<-substr(DHList_Ketola2$PLOT,3,5)
str(DHList2)
DHList2<-merge(DHList_Snyder2,DHList_Ketola2,by=c('Entry','Check','Source','Pedigree','Facultative_Spring_planted','Check_ID'),all = TRUE)
str(Common_sample)

Full_samples
Full_SamplesSub<-Full_samples[,c('Entry','Check','Pedigree','Check_ID','Full_sample')]

DHList2_all<-merge(DHList2,Full_SamplesSub,by=c('Entry','Check','Pedigree','Check_ID'),all=TRUE)
DHList2_all$Full_sample
DHList2_all$Full_sample[is.na(DHList2_all$Full_sample)] <- "Half"
write.csv(DHList2_all, file="DHlist_2.csv")
str(DHList2_all)
DHList_for_print<-DHList2_all[,c('PLOT.x','PLOT.y','Entry','Pedigree','Check','Facultative_Spring_planted','Full_sample','Thresh_sample.x','Thresh_sample.y')]

DHList_for_print11000<-DHList_for_print[!is.na(DHList_for_print$PLOT.x),]
dropY12 <- c("PLOT.y","Thresh_sample.y")


DHList_for_print11000<-DHList_for_print11000[,!c(names(DHList_for_print11000) %in% dropY12)]
colnames(DHList_for_print11000)[c(1,7)]<-c("PLOT","Thresh_sample")

dropX11 <- c("PLOT.x","Thresh_sample.x")
DHList_for_print12000<-DHList_for_print12000[,!c(names(DHList_for_print12000) %in% dropX11)]
colnames(DHList_for_print12000)[c(1,7)]<-c("PLOT","Thresh_sample")

DH_ListT<-rbind(DHList_for_print11000,DHList_for_print12000)
nrow(DH_ListT)
DH_ListT=DH_ListT[order(DH_ListT$PLOT),]
DH_ListT$Sample_ID<-c(1:nrow(DH_ListT))
DH_Full_amount = rbind(DH_ListT,DH_ListT,DH_ListT,DH_ListT)

DH_Full_amount=DH_Full_amount[order(DH_Full_amount$PLOT),]
nrow(DH_Full_amount)
DH_Full_amount
DH_Full_amount$Purpose = '30kernel'


write.csv(DH_Full_amount,"DH_Full_amount.csv")
View(DH_Full_amount)


### make envelopes

setwd('C:/Users/Karl/Google Drive/Grad School_/Field 2020/Germination_assays')

packetsList = data.frame(read_excel(path = 'GGS_EntriesandPlots.xlsx'))

set.seed(1)
packetsList$Sort = runif(length(packetsList$PLOT),0,1)

#packetsList = packetsList[order(packetsList$Sort),]

#packetsList$SampleGroup = c(rep('a',1122/2),rep('b',1122/2))

packetsList$GroupOrder = paste(packetsList$SampleGroup,packetsList$PLOT)
packetsList = packetsList[order(packetsList$GroupOrder),]
packetsList$number = c(1:561,1:561)
packetsList$SampleID = paste(packetsList$SampleGroup,packetsList$number)

packetsListFull = rbind(packetsList,packetsList,packetsList,packetsList)
packetsListFull=packetsListFull[order(packetsListFull$GroupOrder),]
packetsListFull$Purpose = '30kernel'

write.csv(packetsListFull,file = 'PacketsToPrint30kernelAssays.csv')


# write.csv(packetsListMQ, file = 'Data/GerminationTests_2020/PrepForAssays/MaltingQualityPackets.csv')


######################
BagLocations = data.frame(read_excel('Data/GerminationTests_2020/PrepForAssays/PHSsubsample_PMdates_MassCollected.xlsx', sheet = 'MassCollected'))
bagLocationsWithSampleNumer = merge(x = BagLocations, y = packetsList[,c('PLOT','SampleID','GroupOrder')], by = 'PLOT')
write.csv(bagLocationsWithSampleNumer, file = 'Data/GerminationTests_2020/PrepForAssays/Sorting.csv')