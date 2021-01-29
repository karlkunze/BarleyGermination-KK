#aligning genome
wdh_raw<-read.delim("C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_raw.hmp.txt",header = T,sep = "\t",fill = TRUE)
original_config<-as.data.frame(str(wdh_raw))
wdh_raw0<-wdh_raw
original_config
colnames(wdh_raw)=gsub("\\.","-",colnames(wdh_raw))
str(wdh_raw)
orginal_colnames<-colnames(wdh_raw)[1:11]
orginal_colnames
#read.delim
#write as txt file
wdh_raw[1:5,1:5]

######
library(readxl)
library(openxlsx)
colnames(wdh_raw)[ncol(wdh_raw)]
Barley50K<-read_xlsx("data/Genotype_data/50K_Physical_positions.xlsx", sheet = 'EXCAP markers')
#there is also a 9K sheet if needed, although not all positions from the 9K are validated

Barley50K=Barley50K[,c(1:3)]
str(Barley50K)
library(stringr)
Barley50K$`EXCAP marker name`<-str_split_fixed(Barley50K$`EXCAP marker name`,"-0_" , 2)[,1]
colnames(Barley50K)<-c("Name","Chr","bp")
Barley50K$Type<-"50K"

#9K
Barley9K<-as.data.frame(read_xlsx("data/Genotype_data/50K_Physical_positions.xlsx", sheet = '9k markers'))
str(Barley9K)
Barley9K<-Barley9K[,c(1,2)]
str(Barley9K)
Barley9K$chr<-sapply(strsplit(Barley9K$`final chromosome:position`,":"), "[", 1)
Barley9K$chr<-sapply(strsplit(Barley9K$chr,"r"), "[", 2)
Barley9K$bp<-sapply(strsplit(Barley9K$`final chromosome:position`,":"), "[", 2)
Barley9K$bp<-sapply(strsplit(Barley9K$bp,"\\("), "[", 1)
Barley9K<-Barley9K[,-2]
str(Barley9K)
colnames(Barley9K)<-c("Name","Chr","bp")
Barley9K$Type<-"9K"
Barley50_9K<-rbind(Barley50K,Barley9K)
table(Barley50_9K$Chr)
Barley50_9K<-Barley50_9K[with(Barley50_9K, order(Chr, bp)),]
Barley50_9K$Name
##

str(wdh_raw)
Barley50_9K$Name
colnames(wdh_raw)
wdh_raw$`rs-`
Alignment_list<-Barley50_9K[Barley50_9K$Name%in%wdh_raw$`rs-`,]
colnames(Alignment_list)[1]<-"Marker"
colnames(wdh_raw)[1]<-"Marker"
str(wdh_raw)
colnames(wdh_raw)
wdh_raw<-merge(wdh_raw,Alignment_list,by="Marker",all.x = TRUE)
wdh_raw$chrom
#need bp to be numeric for correct ordering
wdh_raw$bp<-as.numeric(wdh_raw$bp)
#order based on Chr and base pair position
wdh_raw<-wdh_raw[with(wdh_raw, order(Chr, bp)),]
#find corresponding names
colnames(wdh_raw)[c(3,4,ncol(wdh_raw)-2,ncol(wdh_raw)-1)]
#replcace original Chromosome/position data with 50K/9K data
wdh_raw[,c(3,4)]<-wdh_raw[,c(ncol(wdh_raw)-2,ncol(wdh_raw)-1)]
#remove attached 50K/9K data at end of dataframe
wdh_raw<-wdh_raw[,-c(ncol(wdh_raw)-2,ncol(wdh_raw)-1,ncol(wdh_raw))]
#rename
orginal_colnames
colnames(wdh_raw)[1:11]<-orginal_colnames
colnames(wdh_raw)[1:11]
table(wdh_raw$chrom)
wdh_raw[1:10,1:12]
str(wdh_raw)
str(wdh_raw0)
wdh_raw0[1:5,1:5]
wdh_raw0
colnames(wdh_raw)=gsub("-","\\.",colnames(wdh_raw))

table(wdh_raw$alleles)
colnames(wdh_raw)[c(1,6,10)]
colnames(wdh_raw)[1:12]
colnames(wdh_raw0)[1:12]
colnames(Tutorial)[1:12]
colnames(Misc)[1:12]
#colnames(wdh_raw)[c(1,6,10)]<-c("rs#","assembly#","panel")
wdh_raw$chrom<-as.integer(gsub("H","",wdh_raw$chrom))
wdh_raw[1:5,12:14]
colnames(wdh_raw)[1:14]
library(janitor)
library(diffdf)
wdh_raw$chrom<-as.numeric(wdh_raw$chrom)
wdh_raw$chrom<-as.integer(wdh_raw$chrom)
wdh_raw$pos<-as.integer(wdh_raw$pos)
#wdh_raw$`assembly-`<-as.character(wdh_raw$`assembly-`)
wdh_raw$center<-as.character(wdh_raw$center)
#wdh_raw$panel<-as.character(wdh_raw$panel)



#wdh_raw<-as.data.frame(wdh_raw)
str(wdh_raw)
rownames(Tutorial)[1:10]
rownames(wdh_raw)[1:10]
wdh_raw["3138",]
wdh_raw[12,]
Tutorial[12,]
rownames(wdh_raw)<-1:nrow(wdh_raw)
wdh_raw$`rs-`<-gsub("-","\\.",wdh_raw$`rs-`)
#colnames(wdh_raw)[12:ncol(wdh_raw)]<-paste("X",colnames(wdh_raw)[12:ncol(wdh_raw)],sep = "")
colnames(wdh_raw)
wdh_raw[-1,]
wdh_raw[2,]
table(Tutorial[2,])
colnames(wdh_raw)[ncol(wdh_raw)-3]
wdh_raw<-wdh_raw[with(wdh_raw, order("chrom", "pos")),]
wdh_raw[is.na(c(12:ncol(wdh_raw)))]



write.table(wdh_raw,"C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_rwas_positions.hmp.txt",sep = "\t",quote=FALSE)
checkFrame[1,1:14]
wdh_raw0[1,1:14]
checkFrame<-read.delim("C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_rwas_positions.hmp.txt",header=T,sep = "\t")
Tutorial<-read.delim("C:/Users/Karl/50k_data/TASSELTutorialData5/TASSELTutorialData5/mdp_genotype.hmp.txt",header = T,sep = "\t")

Misc<-read.delim("data/Genotype_data/CU_WMB_SMB_50k.hmp.txt",header = T,sep ="\t")
str(Misc)
write.table(Misc,"data/Genotype_data/CU_WMB_SMB_50k_test.hmp.txt",sep ="\t")
checkFrame[1:2,1:15]
Tutorial[1:2,1:15]
Misc[1:2,1:15]
length(checkFrame[1,])
length(Tutorial[1,])
str(checkFrame)
str(Tutorial)
diffdf(checkFrame,Tutorial)
diffdf(wdh_raw,wdh_raw0)


str(Tutorial)
str(wdh_raw)
#wdh_raw0[1:5,1:5]
#checkFrame[1:5,1:5]
