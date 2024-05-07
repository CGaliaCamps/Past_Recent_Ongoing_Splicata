
##Primer aneu a Session, set working directory i seleccioneu el directori on es trobaran tots els fitxers d'entrada i sortida

##A continuaci? carreguem els paquets que farem servir

if (!require("vcfR")) {
  install.packages("vcfR", dependencies = TRUE)
library(vcfR)
}
if (!require("adegenet")) {
    install.packages("adegenet", dependencies = TRUE)
library(adegenet)
}
if (!require("hierfstat")) {
  install.packages("hierfstat", dependencies = TRUE)
  library(hierfstat)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
library(RColorBrewer)
}
  if (!require("gplots")) {
    install.packages("gplots", dependencies = TRUE)
library(gplots)
}
    if (!require("ggplot2")) {
      install.packages("ggplot2", dependencies = TRUE)
library(ggplot2)
    }
if (!require("mice")) {
  install.packages("mice", dependencies = TRUE)
  library(mice)
}


  
##Si alg?n paquet no s'instala b?, aneu a Packages i busqueu-los a Install

########### A) PREPARACI? DEL FITXER GENIND I ESTRUCTURA GENERAL############

setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/vcfs/Chr_noInv")

#Ara llegirem el fitxer de genotips vcf i el passarem a genind
vcf <- read.vcfR("Chr_noInv_2mac_5-172DP_70_2bRAD_nou_NoNCSC.recode.vcf")
sdata<- vcfR2genind(vcf)

##Ara llegirem el fitxer de metadata
metadata <- read.table("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/sampleinfo.txt", sep = "\t", header = TRUE)

#####################an?lisis per poblacionals###################

colnames(metadata) <- c("Indv", "Hap", "Pop", "Region")
metadata <- metadata[!grepl(c("NC"), metadata$Indv),]
metadata <- metadata[!grepl(c("SC"), metadata$Indv),]
pop(sdata) <- metadata$Pop

######## Comandes ?tils #############

# Number of individuals
nInd(sdata)

# Number of loci
nLoc(sdata)

# Number of Populations
nPop(sdata)

#Sample size
table(pop(sdata))
barplot(table(pop(sdata)), ylab = "N? Individus")

### DAPC  ######

##xvalDapc per a determinar el nombre de PCoAs a retindre
#crossvalidation with xvalDa
x = tab(sdata, NA.method = "mean") ##Aix? ?s per imputar NAs

xval_Results <- xvalDapc(x, metadata$Pop, n.pca.max = 87, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 1000, xval.plot = TRUE)
pdf("valResults.pdf")
xval_Results
dev.off()

xval_Results[2:6] # Number of PCs Achieving Lowest MSE = 80 ens quedem aquesta xifra. Triga una mica. Est? paral?lelitzat per windows (perque vagi m?s depressa), en linux canvia snow per multicore
xval_Results
#The number of PCs associated with the lowest Mean Squared Error is then retained in the DAPC
pdf("valResults_boxplot__.pdf")
boxplot(xval_Results$`Cross-Validation Results`$success ~ xval_Results$`Cross-Validation Results`$n.pca , xlab="Number of PCA components",
        ylab="Classification succes", main="DAPC - cross-validation")
dev.off()

npca <- as.numeric(xval_Results$`Number of PCs Achieving Lowest MSE`)

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/AIC_POP_Chr_noInv_NoNCSC.pdf")
grpA <- find.clusters(sdata, max.n.clust=18, stat="AIC",n.pca=npca)#amb el gr?fic decidim el nombre de cl?sters i ho posem
#grpA
dev.off()

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/BIC_POP_Chr_noInv_NoNCSC.pdf")
grpB3 <- find.clusters(sdata, max.n.clust=18, stat="BIC",n.pca = npca)#amb el gr?fic decidim el nombre de cl?sters i ho posem
grpB4 <- find.clusters(sdata, max.n.clust=18, stat="BIC",n.pca = npca)#amb el gr?fic decidim el nombre de cl?sters i ho posem
#grpB
dev.off()

table(pop(sdata), grpA$grp)#quants individus per poblaci? hi ha per cluster
table(pop(sdata), grpB3$grp)#quants individus per poblaci? hi ha per cluster
assignacions <- table(pop(sdata), grpB3$grp)#quants individus per poblaci? hi ha per cluster

write.table(assignacions, "assignacions_dapc_noSCNC.txt") 

#Representaci? gr?fica
d <- dapc(sdata, n.pca=15, n.da = 15) ##El nombre de DAPCs retinguts es N/3=136

mycol <- c("#9FCAE6","#9C0824","#559F52","#B3E0A6","#5F90BB", #AGA,AU,BAR,BLA,BRA
           "#EAC0BD","#B9DDF1","#AD0B1E","#3F6E9A","#EC6756", #CA,FE,HK,KN,MIS
           #"#000000", #NC,
           "#CC1619","#2A5783","#F8968C", #OKI,PE,SAK,
           #"grey50",  #SC
           "#7CADD2","#24693D","#DA3E2F") #TE,VIL,WAK

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/dapc_tots_Chr_noInv_NoNCSC_1vs2.pdf", height=4, width=4)
scatter(d,1,2, scree.da=T, scree.pca=T, bg="white", 
        pch=20, cell=0, cstar=0, solid=.8,
        cex=4,clab=0, leg=F, col = mycol, posi.da="topleft", posi.pca="topright") 
dev.off()

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/dapc_tots_Chr_noInv_NoNCSC_2vs3.pdf", height=4, width=4)
scatter(d,2,3, scree.da=T, scree.pca=T, bg="white", 
        pch=20, cell=0, cstar=0, solid=.8,
        cex=4,clab=0, leg=F, col = mycol, posi.da="topleft", posi.pca="topright") 
dev.off()

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/dapc_tots_Chr_noInv_NoNCSC_3vs4.pdf", height=4, width=4)
scatter(d,3,4, scree.da=T, scree.pca=T, bg="white", 
        pch=20, cell=0, cstar=0, solid=.8,
        cex=4,clab=0, leg=F, col = mycol, posi.da="topleft", posi.pca="topright") 
dev.off()

#Plot on how each individual assigns to each population
pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/assi_DAPC_Chr_noInv_NoNCSC.pdf")
assignplot(d, cex.lab=0.4,pch=15)
dev.off()

#STRUCTURE like DAPC plot
lab <- as.vector(metadata$Indv)

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/comp_DAPC_AKA_Chr_noInv_NoNCSC.pdf")
compoplot(d, posi=list(x=20,y=1.17),lab=lab, show.lab=T, cleg=1, col.pal=mycol,
                  txt.leg=popNames(sdata))
dev.off()

#Marcadors que més cobntribueixen a l’estructura

boxplot<-d$var.contr
boxplot1 <- boxplot[,1]
summary(boxplot1)
thres <- quantile(boxplot1, 0.99)

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/contri_DAPC_Chr_noInv_NoNCSC_eix1.pdf", height=3.5, width=5)
loadingplot(d$var.contr, axis=1, thres=thres, lab.jitter=0)
dev.off()

loading1 <- loadingplot(d$var.contr, axis=1, thres=thres, lab.jitter=0)
write.table(loading1, "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/axis1_noNCSN.txt")

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/contri_DAPC_Chr_noInv_NoNCSC_eix2.pdf", height=3.5, width=5)
loadingplot(d$var.contr, axis=2, thres=thres, lab.jitter=0)
dev.off()

loading2 <- loadingplot(d$var.contr, axis=2, thres=thres, lab.jitter=0)
write.table(loading2, "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/axis2_noNCSN.txt")

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/contri_DAPC_Chr_noInv_NoNCSC_eix3.pdf", height=3.5, width=5)
loadingplot(d$var.contr, axis=3, thres=thres, lab.jitter=0, lab="")
dev.off()

loading3 <- loadingplot(d$var.contr, axis=3, thres=thres, lab.jitter=0)
write.table(loading3, "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/axis3_noNCSN.txt")

pdf("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/contri_DAPC_Chr_noInv_NoNCSC_eix4.pdf", height=3.5, width=5)
loadingplot(d$var.contr, axis=4, thres=thres, lab.jitter=0, lab="")
dev.off()

loading4 <- loadingplot(d$var.contr, axis=4, thres=thres, lab.jitter=0)
write.table(loading4, "E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC/axis4_noNCSN.txt")



Chromosome_3_26369514 <- tab(genind2genpop(sdata[loc=c("Chromosome_3_26369514")]),freq=TRUE)

par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,.1),las=3)
matplot(Chromosome_3_26369514, pch=c("a","c"), type="b",
        xlab="locality",ylab="allele frequency", xaxt="n",
        cex=1.5, main="26369514")
axis(side=1, at=1:16, lab=c("AGA","AU","BAR","BLA","BRA", "CA","FE","HK","KN","MIS","OKI","PE","SAK","TE","VIL","WAK"))



dev.off()








