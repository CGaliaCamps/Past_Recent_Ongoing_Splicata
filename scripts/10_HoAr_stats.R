if (!require("vcfR")) {
    install.packages("vcfR", dependencies = TRUE)
    library(vcfR)}
  
  if (!require("adegenet")) {
      install.packages("adegenet", dependencies = TRUE)
      library(adegenet)}
    
  
  if (!require("hierfstat")) {
      install.packages("hierfstat", dependencies = TRUE)
      library(hierfstat)}

library(ggplot2)
library(ggnewscale)
    
    
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

### Data:
```{r, message=FALSE, warning=FALSE}
setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Diversitat")

ArA<-read.table("./Ar/meanAllele_Chr_noInv.csv", header = TRUE, sep = ",")
Ar2<-read.table("./Ar/meanAllele_Chr2.csv", header = TRUE, sep = ",", dec = ".")
Ar4<-read.table("./Ar/meanAllele_Chr4.csv", header = TRUE, sep = ",", dec = ".")
Ar11<-read.table("./Ar/meanAllele_Chr11.csv", header = TRUE, sep = ",", dec = ".")
Ar16<-read.table("./Ar/meanAllele_Chr16.csv", header = TRUE, sep = ",", dec = ".")

HoA<-read.table("./Ho/meanHo_Chr_noInv.csv", header = TRUE, sep = ",", dec = ".")
Ho2<-read.table("./Ho/meanHo_Chr2.csv", header = TRUE, sep = ",", dec = ".")
Ho4<-read.table("./Ho/meanHo_Chr4.csv", header = TRUE, sep = ",", dec = ".")
Ho11<-read.table("./Ho/meanHo_Chr11.csv", header = TRUE, sep = ",", dec = ".")
Ho16<-read.table("./Ho/meanHo_Chr16.csv", header = TRUE, sep = ",", dec = ".")

```

### Add metadata to each file:

```{r, message=FALSE, warning=FALSE}

ArA<-cbind(ArA, "DATASET"= "1NoInv") 
Ar2<-cbind(Ar2, "DATASET"= "2Chr2") 
Ar4<-cbind(Ar4, "DATASET"= "3Chr4") 
Ar11<-cbind(Ar11, "DATASET"= "Chr11") 
Ar16<-cbind(Ar16, "DATASET"= "Chr16")

HoA<-cbind(HoA, "DATASET"= "1NoInv") 
Ho2<-cbind(Ho2, "DATASET"= "2Chr2") 
Ho4<-cbind(Ho4, "DATASET"= "3Chr4") 
Ho11<-cbind(Ho11, "DATASET"= "Chr11") 
Ho16<-cbind(Ho16, "DATASET"= "Chr16")

ArA<-cbind(ArA, "ANALYSIS"= "1Ar") 
Ar2<-cbind(Ar2, "ANALYSIS"= "1Ar") 
Ar4<-cbind(Ar4, "ANALYSIS"= "1Ar") 
Ar11<-cbind(Ar11, "ANALYSIS"= "1Ar") 
Ar16<-cbind(Ar16, "ANALYSIS"= "1Ar") 

HoA<-cbind(HoA, "ANALYSIS"= "3Ho") 
Ho2<-cbind(Ho2, "ANALYSIS"= "3Ho") 
Ho4<-cbind(Ho4, "ANALYSIS"= "3Ho")
Ho11<-cbind(Ho11, "ANALYSIS"= "3Ho") 
Ho16<-cbind(Ho16, "ANALYSIS"= "3Ho") 

ArA$X<-gsub("Ar.", "", ArA$X)
Ar2$X<-gsub("Ar.", "", Ar2$X)
Ar4$X<-gsub("Ar.", "", Ar4$X)
Ar11$X<-gsub("Ar.", "", Ar11$X)
Ar16$X<-gsub("Ar.", "", Ar16$X)

ArA<-ArA[-1,]
Ar2<-Ar2[-1,]
Ar4<-Ar4[-1,]
Ar11<-Ar11[-1,]
Ar16<-Ar16[-1,]

sum_stats<-rbind(ArA,Ar2,Ar4,Ar11,Ar16,HoA,Ho2,Ho4,Ho11,Ho16)
colnames(sum_stats) <- c("population", "value", "Dataset", "Analysis")

sum_stats$population<-gsub("agadir", "17AGA", sum_stats$population)
sum_stats$population<-gsub("barcelona", "14BAR", sum_stats$population)
sum_stats$population<-gsub("australia", "28AU", sum_stats$population)
sum_stats$population<-gsub("blanes", "13BLA", sum_stats$population)
sum_stats$population<-gsub("brasil", "19BRA", sum_stats$population)
sum_stats$population<-gsub("california", "22CA", sum_stats$population)
sum_stats$population<-gsub("ferrol", "16FE", sum_stats$population)
sum_stats$population<-gsub("knysna", "20KN", sum_stats$population)
sum_stats$population<-gsub("misaki", "24MIS", sum_stats$population)
sum_stats$population<-gsub("northcarolina", "11NC", sum_stats$population)
sum_stats$population<-gsub("southcarolina", "12SC", sum_stats$population)
sum_stats$population<-gsub("okinashima", "26OKI", sum_stats$population)
sum_stats$population<-gsub("porthelisabeth", "21PE", sum_stats$population)
sum_stats$population<-gsub("sakushima", "23SAK", sum_stats$population)
sum_stats$population<-gsub("tenerife", "18TE", sum_stats$population)
sum_stats$population<-gsub("vilanova", "15VIL", sum_stats$population)
sum_stats$population<-gsub("wakayama", "25WAK", sum_stats$population)
sum_stats$population<-gsub("hongkong", "27HK", sum_stats$population)


pdf(file="bubbleplot_genomic_stats.pdf", width=4.5, height=9)
ggplot(data=sum_stats, aes(x=Analysis, y= population))+

  new_scale("fill") + new_scale("size") +
  
  geom_point(data=sum_stats[sum_stats$Analysis == "1Ar",], aes(x=Analysis, y=population, fill=value, size=value), alpha=1, shape=21, color="black", stroke=1) +
  scale_size(range = c(2,5), name="1Ar") +
scale_fill_gradient("1Ar",  limits = c(1.1,1.3), low="#ffc40c", high="#ffc40c")+
 
   new_scale("fill") + new_scale("size") +

    geom_point(data=sum_stats[sum_stats$Analysis == "3Ho",], aes(x=Analysis, y=population, fill=value, size=value), alpha=1, shape=21, color="black", stroke=1) +
  scale_size(range = c(1,8), name="3Ho", breaks = c("0.10", "0.15", "0.20", "0.25")) +
  scale_fill_gradient("2He", low= "#d0a016",  high="#d0a016")+

   theme_classic()+
  scale_y_discrete(limits=rev)+
   theme(text = element_text(size = 12))+
   theme(legend.position = "none", axis.text.x = element_text(angle = 90, size = 8, face="bold"))+
  
  facet_grid(~Dataset)
dev.off()




pdf(file="bubbleplot_genomic_stats_legends.pdf", width=9, height=16)
ggplot(data=sum_stats, aes(x=Analysis, y=population))+

  new_scale("fill") + new_scale("size") +
  
  geom_point(data=sum_stats[sum_stats$Analysis == "1Ar",], aes(x=Analysis, y=population, fill=value, size=value), alpha=1, shape=21, color="black", stroke=1) +
  scale_size(range = c(2,5), name="1Ar") +
scale_fill_gradient("1Ar",  limits = c(1.1,1.3), low= "#ffc40c", high="#ffc40c")+
 
new_scale("fill") + new_scale("size") +
  
     geom_point(data=sum_stats[sum_stats$Analysis == "3Ho",], aes(x=Analysis, y=population, fill=value, size=value), alpha=1, shape=21, color="black", stroke=1) +
  scale_size(range = c(1,8), name="3Ho") +
  scale_fill_gradient("2He", low= "#d0a016",  high="#d0a016")+
   
  theme_classic()+
 scale_y_discrete(limits=rev)+
  theme(text = element_text(size = 15))+
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, size = 10, face="bold"))+
  
  facet_grid(~Dataset)
dev.off()
  
  

##download for manual editing
write.table(sum_stats, "taula_stats_summary.txt", sep="\t")


    
