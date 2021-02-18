setwd("C:/Users/Karl/git/BarleyGermination-KK")
library(here)
load(here("data/Genotype_data/GWAS.Field.hits.Rdata"))
load(here("data/Genotype_data/GWAS.TP.hits.Rdata"))
#chromosome 0 not 1
#end position not inclusive
#format example  chr1  11873  14409  uc001aaa.3
getwd()
str(Field_sign)
library(bedr)
str(Ger)
Field_sign$start=Field_sign$Position-1
Field_sign$end=Field_sign$Position+1
Field_sign$Chromosome=Field_sign$Chromosome-1
Field=Field_sign[,c('Chromosome','start','end','SNP')]
Field$Chromosome=paste0("chr",Field$Chromosome)
colnames(Field)<-c("chr","start","end","name")
Field=Field[with(Field,order(chrom,start)),]
Field=distinct(Field)
library(dplyr)

Field_bed3.bed=convert2bed(Field[,1:3],
  set.type = TRUE,
  check.zero.based = TRUE,
  check.chr = TRUE,
  check.valid = TRUE,
  check.sort = TRUE,
  check.merge = TRUE,
  verbose = TRUE
)


write.table(Field_bed3.bed,file = "data/Genotype_data/Field_sig_SNPS.bed",sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

str(GWAS.sign.TP_all)
GWAS.sign.TP_all$start=GWAS.sign.TP_all$Position-1
GWAS.sign.TP_all$end=GWAS.sign.TP_all$Position+1
GWAS.sign.TP_all$Chromosome=GWAS.sign.TP_all$Chromosome-1
TP_all=GWAS.sign.TP_all[,c('Chromosome','start','end','SNP')]
TP_all$Chromosome=paste0("chr",TP_all$Chromosome)
colnames(TP_all)<-c("chr","start","end","name")
TP_all=TP_all[with(TP_all,order(chr,start)),]

TP_all=distinct(TP_all)
library(dplyr)
TP_all_bed3.bed=convert2bed(TP_all[,1:3])
write.table(TP_all,file = "data/Genotype_data/GermTP_SNPS.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
