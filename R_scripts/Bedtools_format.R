setwd("C:/Users/Karl/git/BarleyGermination-KK")
library(here)
load(here("data/Genotype_data/GWAS.Field.hits.Rdata"))
load(here("data/Genotype_data/GWAS.TP.hits.Rdata"))
#chromosome 0 not 1
#end position not inclusive
#format example  chr1  11873  14409  uc001aaa.3
getwd()
str(Field_sign)
Field_sign$start=Field_sign$Position-1
Field_sign$end=Field_sign$Position+1
Field_sign$Chromosome=Field_sign$Chromosome-1
Field=Field_sign[,c('Chromosome','start','end','SNP')]
Field$Chromosome=paste0("chr",Field$Chromosome)
colnames(Field)<-c("chrom","start","end","name")
write.table(Field,file = "data/Genotype_data/Field_sig_SNPS.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

str(GWAS.sign.TP_all)
GWAS.sign.TP_all$start=GWAS.sign.TP_all$Position-1
GWAS.sign.TP_all$end=GWAS.sign.TP_all$Position+1
GWAS.sign.TP_all$Chromosome=GWAS.sign.TP_all$Chromosome-1
TP_all=GWAS.sign.TP_all[,c('Chromosome','start','end','SNP')]
TP_all$Chromosome=paste0("chr",TP_all$Chromosome)
colnames(TP_all)<-c("chrom","start","end","name")
TP_all
write.table(Field,file = "data/Genotype_data/GermTP_SNPS.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
