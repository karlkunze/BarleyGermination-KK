rTASSEL::startLogger(fullPath = NULL, fileName = NULL)
#options(java.parameters = c("-Xmx6g"))
library(rTASSEL)
path_all<-"data/Genotype_data/DH_RIL_parents.hmp.txt"

# getwd()
# path_DH<-"/home/karl/Documents/Grad_school_local/WMB_genotype_data/DH_2020_sort.vcf"
# #need to sort the vcf file before uploadinf to Rtassel, you can currently do this with the Tassel GUI
# path_RIL<-"/home/karl/Documents/Grad_school_local/WMB_genotype_data/RIL_sort.vcf"
# path_WMB_parents2019<-"data/genotypes/2019/Winter_2019.vcf"

# DH_VCF <- rTASSEL::readGenotypeTableFromPath(
#   path = path_DH
# )
# RIL_VCF <- rTASSEL::readGenotypeTableFromPath(
#   path = path_RIL
# )
# parents_VCF<-rTASSEL::readGenotypeTableFromPath(
#   path = path_WMB_parents2019
# )
all_HMP<-rTASSEL::readGenotypeTableFromPath(
  path = path_all
)


all_HMP
679*0.2
?filterGenotypeTableSites()
# DH_filter<-filterGenotypeTableSites(tasObj = DH_VCF,siteMinAlleleFreq = 0.05, siteMaxAlleleFreq = 1.0,siteRangeFilterType = "none",siteMinCount = 38,maxHeterozygous = 0.01)
# RIL_filter<-filterGenotypeTableSites(tasObj = RIL_VCF,siteMinAlleleFreq = 0.05, siteMaxAlleleFreq = 1.0,siteRangeFilterType = "none",siteMinCount = 96
#                          ,maxHeterozygous = 0.05)
# parents_filter<-filterGenotypeTableSites(tasObj = parents_VCF,siteMinAlleleFreq = 0.05, siteMaxAlleleFreq = 1.0,siteMinCount = 2,siteRangeFilterType = "none"
#                                      ,maxHeterozygous = 0.05)
all_HMP
679*0.2
library(rTASSEL)
all_filter<-filterGenotypeTableSites(tasObj = all_HMP,siteMinAlleleFreq = 0.05, siteMaxAlleleFreq = 1.0,siteMinCount = 136,siteRangeFilterType = "none"
                                     ,maxHeterozygous = 0.05)
all_filter
all_filter<-filterGenotypeTableTaxa(all_filter,minNotMissing = 0.5)
all_filter
# we will have to filter on hets/DHs/parent lines later for heterozygosity

# exportGenotypeTable(DH_filter,file = "data/genotypes/DH/DH_tassel_filter",format = "hapmap",)
# exportGenotypeTable(RIL_filter,file = "data/genotypes/RIL/RIL_tassel_filter",format = "hapmap")
# exportGenotypeTable(parents_filter,file = "data/genotypes/parents_2019_filter",format = "hapmap")
exportGenotypeTable(all_filter,file = "data/Genotype_data/all_filter",format = "hapmap")
?exportGenotypeTable()
str(RIL_filter)
DH_filter
