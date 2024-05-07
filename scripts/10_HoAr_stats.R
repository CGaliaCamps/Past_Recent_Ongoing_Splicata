if (!require("vcfR")) {
    install.packages("vcfR", dependencies = TRUE)
    library(vcfR)}
  
  if (!require("adegenet")) {
      install.packages("adegenet", dependencies = TRUE)
      library(adegenet)}
    
  
  if (!require("hierfstat")) {
      install.packages("hierfstat", dependencies = TRUE)
      library(hierfstat)}
    
    
########### A)PREAPRE YOUR FILE (NOTE THAT EACH DATASET HAS TO BE SPECIFIED)############

setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/vcfs/Chr_noInv")

#read your vcf and transform it to genind
vcf <- read.vcfR("Chr_noInv_2mac_5-172DP_70_2bRAD_nou.recode.vcf")
sdata<- vcfR2genind(vcf)

##Open your metadata
metadata <- read.table("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/sampleinfo.txt", sep = "\t", header = TRUE)

#####################locality-level analyses###################
rdata <- sdata
hdata <- sdata

colnames(metadata) <- c("Indv", "Hap", "Pop", "Region")
pop(sdata) <- metadata$Pop
pop(rdata) <- metadata$Region
pop(hdata) <- metadata$Hap

sdatahf<-genind2hierfstat(sdata,pop=NULL)
rdatahf<-genind2hierfstat(rdata,pop=NULL)
hdatahf<-genind2hierfstat(hdata,pop=NULL)
    
#####BASIC STATS
    
basics<-basic.stats(sdatahf)
alleleRich<-allelic.richness(sdatahf)

    
write.csv(basics$Ho, file = "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat/Ho/Ho_Chr_noInv.csv")
write.csv(alleleRich, file = "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat/Allele richness/Allele_Chr_noInv.csv")

Ho<- read.csv("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat/Ho/Ho_Chr_noInv.csv", sep = ",", header = TRUE)
ALE<- read.csv("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat/Allele richness/Allele_Chr_noInv.csv", sep = ",", header = TRUE)

meanHo <- colMeans(Ho[sapply(Ho, is.numeric)], na.rm=TRUE)
meanAle <- colMeans(ALE[sapply(ALE, is.numeric)], na.rm=TRUE)

write.csv(meanHo, file = "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat/Ho/meanHo_Chr_noInv.csv")
write.csv(meanAle, file = "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat/Allele richness/meanAllele_Chr_noInv.csv")

boot.ppfis(dat=sdatahf,nboot = 100,quant = c(0.25,0.75),diploid = TRUE)    

    
