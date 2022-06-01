#aligning genome
wdh_raw<-as.data.frame(read.delim("C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_raw.hmp.txt",header = T,sep = "\t"))
#WDH_50K_raw.hmp <- read.delim("C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_raw.hmp.txt", header=FALSE, row.names=NULL, quote="")
library(janitor)
WMB_SMB<-as.data.frame(read.delim("C:/Users/Karl/git/BarleyGermination-KK/data/Genotype_data/CU_WMB_SMB_50k.hmp.txt",header = T,sep = "\t"))
WMB_SMB$chrom<-as.integer(gsub("H","",WMB_SMB$chrom))
colnames(wdh_raw)
colnames(WMB_SMB)=gsub("\\.","-",colnames(WMB_SMB))

#Tutorial<-read.delim("C:/Users/Karl/50k_data/TASSELTutorialData5/TASSELTutorialData5/mdp_genotype.hmp.txt",header = T,sep = "\t")
#View(wdh_raw)
#write.table(cu_geno, file="GGS_all_50k_nuc.hmp.txt", sep="\t",row.names=F, quote=F)
#write.table(wdh_raw,"C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_raw4.hmp.txt",sep = "\t",quote =F,row.names = F)
#Tutorial_check<-read.delim("C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_raw4.hmp.txt",sep = "\t",header = TRUE)
#library(diffdf)
#diffdf(Tutorial,Tutorial_check)
#ncol(wdh_raw[1,])

#nrow(wdh_raw)
#original_config<-as.data.frame(str(wdh_raw))
wdh_raw0<-wdh_raw
#original_config
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
Morex2019<-read.delim("data/Genotype_data/barley50kMarkerPositions_Morex2019Assembly.gff",skip=1,header = FALSE)
head(Morex2019)
#there is also a 9K sheet if needed, although not all positions from the 9K are validated
Morex2019=Morex2019[,c(1,2,4,7,9)]
head(Morex2019)


library(stringr)
Morex2019$V1<-str_split_fixed(Morex2019$V1,"r" , 2)[,2]
Morex2019$V9<-str_split_fixed(Morex2019$V9,"=" , 2)[,2]
head(Morex2019)
colnames(Morex2019)
colnames(Morex2019)<-c("Chrom","Type","pos","strand","Marker")


#9K
#Barley9K<-as.data.frame(read_xlsx("data/Genotype_data/50K_Physical_positions.xlsx", sheet = '9k markers'))
#str(Barley9K)
#Barley9K<-Barley9K[,c(1,2)]
#str(Barley9K)
#Barley9K$chr<-sapply(strsplit(Barley9K$`final chromosome:position`,":"), "[", 1)
#Barley9K$chr<-sapply(strsplit(Barley9K$chr,"r"), "[", 2)
#Barley9K$bp<-sapply(strsplit(Barley9K$`final chromosome:position`,":"), "[", 2)
#Barley9K$bp<-sapply(strsplit(Barley9K$bp,"\\("), "[", 1)
#Barley9K<-Barley9K[,-2]
#str(Barley9K)
#colnames(Barley9K)<-c("Name","Chr","bp")
#Barley9K$Type<-"9K"
#Barley50_9K<-rbind(Barley50K,Barley9K)
str(Morex2019)
Morex2019[1:10,]

table(Morex2019$Chrom,useNA = "ifany")
Morex2019<-Morex2019[with(Morex2019, order(Chrom, pos)),]
head(Morex2019)
##


Alignment_list<-Morex2019[Morex2019$Marker%in%wdh_raw$`rs-`,]
Alignment_list_WMB<-Morex2019[Morex2019$Marker%in%WMB_SMB$`rs-`,]

str(Alignment_list_WMB)
table(Barley50_9K$Chr,useNA = "ifany")
table(Alignment_list$Chr,useNA = "ifany")
table(Alignment_list_WMB$Chr,useNA = "ifany")

colnames(wdh_raw)[1]<-"Marker"
colnames(WMB_SMB)[1]<-"Marker"
str(wdh_raw)
colnames(Alignment_list)
colnames(wdh_raw)[1:15]
#
table(Alignment_list$Chr,useNA = "ifany")
sum(table(Alignment_list$Chr,useNA = "ifany"))
#
table(wdh_raw$chrom,useNA = "ifany")
sum(table(wdh_raw$chrom,useNA = "ifany"))
hist(Alignment_list$pos)
wdh_raw<-merge(wdh_raw,Alignment_list,by="Marker",all.x = TRUE)
WMB_SMB<-merge(WMB_SMB,Alignment_list_WMB,by="Marker",all.x = TRUE)
table(wdh_raw$Chrom,useNA = "ifany")
table(wdh_raw$Chr,useNA = "ifany")
wdh_raw$chrom
hist(wdh_raw$pos.y)
#need bp to be numeric for correct ordering
wdh_raw$pos.y<-as.numeric(wdh_raw$pos.y)
WMB_SMB$pos.y<-as.numeric(WMB_SMB$pos.y)
colnames(wdh_raw)
wdh_raw[,c(493)]
#remove type column
wdh_raw<-wdh_raw[,-c(493)]
#find corresponding names
ncol(wdh_raw)
colnames(wdh_raw)[c(3:5)]
#replcace original Chromosome/position data with Morex2019 data
wdh_raw[,c(3,4,5)]
wdh_raw[,c(492)]

wdh_raw[,c(3,4,5)]<-wdh_raw[,c(492:494)]

#remove data at end of data frame
wdh_raw<-wdh_raw[,-c(492,493,494)]
colnames(wdh_raw)[1:6]
str(wdh_raw)

length(table(wdh_raw$pos.x,useNA = "ifany"))-length(table(wdh_raw$pos.x))

ncol(wdh_raw)

#rename
orginal_colnames
colnames(wdh_raw)[1:12]
colnames(wdh_raw)[1:11]<-orginal_colnames
wdh_raw<-wdh_raw[with(wdh_raw, order(chrom, pos)),]
#WMB_SMB<-WMB_SMB[with(WMB_SMB, order(Chr, bp)),]
table(wdh_raw$chrom,useNA = "ifany")
wdh_raw[1:10,1:13]
str(wdh_raw)
str(wdh_raw0)
wdh_raw[1:5,1:5]
wdh_raw0
colnames(wdh_raw)=gsub("-","\\.",colnames(wdh_raw))

table(wdh_raw$alleles)
colnames(wdh_raw)[c(1,6,10)]
colnames(wdh_raw)[1:12]
colnames(wdh_raw0)[1:12]
colnames(Tutorial)[1:12]
colnames(Misc)[1:12]
#colnames(wdh_raw)[c(1,6,10)]<-c("rs#","assembly#","panel")
#wdh_raw$chrom<-as.integer(gsub("H","",wdh_raw$chrom))
wdh_raw[1:5,12:14]
colnames(wdh_raw)[1:14]
library(janitor)
library(diffdf)
#wdh_raw$chrom<-as.numeric(wdh_raw$chrom)
#wdh_raw$chrom<-as.integer(wdh_raw$chrom)
wdh_raw$pos<-as.integer(wdh_raw$pos)
#wdh_raw$`assembly-`<-as.character(wdh_raw$`assembly-`)
wdh_raw$center<-as.character(wdh_raw$center)
#wdh_raw$panel<-as.character(wdh_raw$panel)
wdh_raw[1:5,1:14]

#############################################
#wdh_raw<-as.data.frame(wdh_raw)
str(wdh_raw)
rownames(Tutorial)[1:10]
rownames(wdh_raw)[1:10]
wdh_raw["3138",]
wdh_raw[12,]
#Tutorial[12,]

#rownames(wdh_raw)<-1:nrow(wdh_raw)
wdh_raw
wdh_raw$rs.<-gsub("-","\\.",wdh_raw$rs.)
wdh_raw[1:5,1:15]
#colnames(wdh_raw)[12:ncol(wdh_raw)]<-paste("X",colnames(wdh_raw)[12:ncol(wdh_raw)],sep = "")
colnames(wdh_raw)
typeof(wdh_raw$chrom)

#PROBLEM HEREE
str(wdh_raw)
wdh_raw[with(wdh_raw, order(chrom,pos)),]
wdh_raw$chrom
wdh_raw<-as.data.frame(wdh_raw[with(wdh_raw, order(chrom,pos)),])
wdh_raw[is.na(wdh_raw$chrom),]
wdh_raw[!is.na(wdh_raw$alleles),]
wdh_raw<-as.data.frame(wdh_raw[!is.na(wdh_raw$alleles),])
wdh_raw<-as.data.frame(wdh_raw[!is.na(wdh_raw$pos),])

#wdh_raw worked, I removed NA pos and Chr positions and allele
dim(wdh_raw)
#wdh_raw[is.na(wdh_raw[,4])]
wdh_raw[1:4,1:14]
wdh_raw[1:10,1:10]
table(wdh_raw$chrom,useNA = "ifany")
colnames(wdh_raw)[1:15]
str(wdh_raw)
#WMB_SMB<-as.data.frame(read.delim("C:/Users/Karl/git/BarleyGermination-KK/data/Genotype_data/CU_WMB_SMB_50k.hmp.txt",header = T,sep = "\t"))
str(WMB_SMB$rs.)
#WMB_SMB$rs.<-gsub()
str(wdh_raw$rs.)
WMB_SMB$chrom<-as.integer(gsub("H","",WMB_SMB$chrom))
str(WMB_SMB)
table(WMB_SMB$chrom)
wdh_raw[1:4,1:15]
Alignment_list_WMB
WMB_SMB$rs.<-gsub("-","\\.",WMB_SMB$rs.)
wdh_raw[wdh_raw$rs.%in%WMB_SMB$rs.,]

View(wdh_raw)
wdh_raw$rs.<-gsub("\\.","-",wdh_raw$rs.)

write.table(wdh_raw,"C:/Users/Karl/50k_data/Winter DH population/Raw_tassel/WDH_50K_v2_aligned_positions_test3.hmp.txt",sep = "\t",quote=F,row.names = F)
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
