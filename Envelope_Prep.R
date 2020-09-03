library(readxl)
setwd('C:/Users/Karl/Google Drive/Grad School_/Field 2020/Germination_assays')
masses = data.frame(read_excel(path = 'C:/Users/Karl/Google Drive/Grad School_/Field 2020/Germination_assays/PHSsubsample_PMdates_MassCollected.xlsx', sheet = 'MassCollected'))

library(ggplot2)
ggplot(data = masses, aes(x = MassGrams))+geom_histogram()
masses$numkernels = masses$MassGrams/0.05

masses$assays = trunc(masses$numkernels/30)
ggplot(data = masses, aes(x = assays))+geom_histogram()
sum(ifelse(masses$assays<6,masses$PLOT,0))
ifelse(masses$assays<6,masses$PLOT,NA)

Helfer = subset(masses,masses$PLOT>1999& masses$PLOT<2449)
ggplot(Helfer, aes(x = MassGrams))+geom_histogram()
Helfer$massAfterMalt = Helfer$MassGrams - 5
Helfer$assaysAfterMalt = trunc(Helfer$massAfterMalt/0.05/30)
Helfer$assaysAfterMalt

ggplot(Helfer, aes(x = assaysAfterMalt))+geom_histogram()

CaldwellGGS = subset(masses,masses$PLOT>1000& masses$PLOT<1449)
ggplot(CaldwellGGS, aes(x = MassGrams))+geom_histogram()
CaldwellGGS$massAfterMalt = CaldwellGGS$MassGrams - 5
CaldwellGGS$assaysAfterMalt = trunc(CaldwellGGS$massAfterMalt/0.05/30)
CaldwellGGS$assaysAfterMalt

ggplot(CaldwellGGS, aes(x = assaysAfterMalt))+geom_histogram()

library(gridExtra)
grid.arrange(ggplot(Helfer, aes(x = assaysAfterMalt))+geom_histogram()+
               xlim(-3,35)+ggtitle('Helfer GGS'),
             ggplot(CaldwellGGS, aes(x = assaysAfterMalt))+geom_histogram()+
               xlim(-3,35)+ggtitle('caldwell'))

# Take the MQ sample for micro malting from HELFER
# need 6 packets of everything, and 2 extra for malting quality for 
# the GGS population of lines. 
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

packetsListMQ = subset(packetsList,packetsList$PLOT %in% c(2001:2448))
packetsListMQ = rbind(packetsListMQ,packetsListMQ)
packetsListMQ = packetsListMQ[order(packetsListMQ$GroupOrder),]
packetsListMQ$Purpose = 'MQ 2.5g'

# write.csv(packetsListMQ, file = 'Data/GerminationTests_2020/PrepForAssays/MaltingQualityPackets.csv')


######################
BagLocations = data.frame(read_excel('Data/GerminationTests_2020/PrepForAssays/PHSsubsample_PMdates_MassCollected.xlsx', sheet = 'MassCollected'))
bagLocationsWithSampleNumer = merge(x = BagLocations, y = packetsList[,c('PLOT','SampleID','GroupOrder')], by = 'PLOT')
write.csv(bagLocationsWithSampleNumer, file = 'Data/GerminationTests_2020/PrepForAssays/Sorting.csv')