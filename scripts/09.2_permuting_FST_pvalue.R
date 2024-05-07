library("adegenet")
library("hierfstat")
library("vcfR")

setwd("/piec2/cgalia/Styela_mundial/2bRAD_world/pvalue_FST/")

metadata <- read.table("sampleinfo.txt", sep = "\t", header = TRUE)
colnames(metadata) <- c("Indv", "Hap", "Pop", "Region")

vcf2 <- read.vcfR("Chr2_2mac_5-172DP_70_2bRAD_nou.recode.vcf")
vcf4 <- read.vcfR("Chr4_2mac_5-172DP_70_2bRAD_nou.recode.vcf")
vcf11 <- read.vcfR("Chr11_2mac_5-172DP_70_2bRAD_nou.recode.vcf")
vcf16 <- read.vcfR("Chr16_2mac_5-172DP_70_2bRAD_nou.recode.vcf")
vcfNI <- read.vcfR("Chr_noInv_2mac_5-172DP_70_2bRAD_nou.recode.vcf")

     
sdata2<- vcfR2genind(vcf2)
sdata4<- vcfR2genind(vcf4)
sdata11<- vcfR2genind(vcf11)
sdata16<- vcfR2genind(vcf16)
sdataNI<- vcfR2genind(vcfNI)

rdata2 <- sdata2
rdata4 <- sdata4
rdata11 <- sdata11
rdata16 <- sdata16
rdataNI <- sdataNI

hdata2 <- sdata2
hdata4 <- sdata4
hdata11 <- sdata11
hdata16 <- sdata16
hdataNI <- sdataNI

pop(sdata2) <- metadata$Pop
pop(sdata4) <- metadata$Pop
pop(sdata11) <- metadata$Pop
pop(sdata16) <- metadata$Pop
pop(sdataNI) <- metadata$Pop

pop(rdata2) <- metadata$Region
pop(rdata4) <- metadata$Region
pop(rdata11) <- metadata$Region
pop(rdata16) <- metadata$Region
pop(rdataNI) <- metadata$Region

pop(hdata2) <- metadata$Hap
pop(hdata4) <- metadata$Hap
pop(hdata11) <- metadata$Hap
pop(hdata16) <- metadata$Hap
pop(hdataNI) <- metadata$Hap

sdatahf2<-genind2hierfstat(sdata2,pop=NULL)
sdatahf4<-genind2hierfstat(sdata4,pop=NULL)
sdatahf11<-genind2hierfstat(sdata11,pop=NULL)
sdatahf16<-genind2hierfstat(sdata16,pop=NULL)
sdatahfNI<-genind2hierfstat(sdataNI,pop=NULL)

rdatahf2<-genind2hierfstat(rdata2,pop=NULL)
rdatahf4<-genind2hierfstat(rdata4,pop=NULL)
rdatahf11<-genind2hierfstat(rdata11,pop=NULL)
rdatahf16<-genind2hierfstat(rdata16,pop=NULL)
rdatahfNI<-genind2hierfstat(rdataNI,pop=NULL)

hdatahf2<-genind2hierfstat(hdata2,pop=NULL)
hdatahf4<-genind2hierfstat(hdata4,pop=NULL)
hdatahf11<-genind2hierfstat(hdata11,pop=NULL)
hdatahf16<-genind2hierfstat(hdata16,pop=NULL)
hdatahfNI<-genind2hierfstat(hdataNI,pop=NULL)

popFst2<-pairwise.neifst(sdatahf2,diploid=TRUE)
popFst4<-pairwise.neifst(sdatahf4,diploid=TRUE)
popFst11<-pairwise.neifst(sdatahf11,diploid=TRUE)
popFst16<-pairwise.neifst(sdatahf16,diploid=TRUE)
popFstNI<-pairwise.neifst(sdatahfNI,diploid=TRUE)

ratFst2<-pairwise.neifst(rdatahf2,diploid=TRUE)
ratFst4<-pairwise.neifst(rdatahf4,diploid=TRUE)
ratFst11<-pairwise.neifst(rdatahf11,diploid=TRUE)
ratFst16<-pairwise.neifst(rdatahf16,diploid=TRUE)
ratFstNI<-pairwise.neifst(rdatahfNI,diploid=TRUE)

hatFst2<-pairwise.neifst(hdatahf2,diploid=TRUE)
hatFst4<-pairwise.neifst(hdatahf4,diploid=TRUE)
hatFst11<-pairwise.neifst(hdatahf11,diploid=TRUE)
hatFst16<-pairwise.neifst(hdatahf16,diploid=TRUE)
hatFstNI<-pairwise.neifst(hdatahfNI,diploid=TRUE)

####p-valor FST####
library(parallel)
mat.obs <- as.matrix(popFst2)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf2[,1]),(sdatahf2[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=18,ncol=18)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr2.txt")


#### ####
library(parallel)
mat.obs <- as.matrix(popFst4)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf4[,1]),(sdatahf4[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=18,ncol=18)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr4.txt")

#####
library(parallel)
mat.obs <- as.matrix(popFst11)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf11[,1]),(sdatahf11[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=18,ncol=18)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr11.txt")

#####
library(parallel)
mat.obs <- as.matrix(popFst16)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf16[,1]),(sdatahf16[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=18,ncol=18)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr16.txt")

#####
library(parallel)
mat.obs <- as.matrix(popFstNI)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahfNI[,1]),(sdatahfNI[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=18,ncol=18)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr_noInv.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst2)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf2[,1]),(rdatahf2[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=18,ncol=18)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr2.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst4)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf4[,1]),(rdatahf4[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr4.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst11)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf11[,1]),(rdatahf11[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr11.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst16)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf16[,1]),(rdatahf16[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr16.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFstNI)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahfNI[,1]),(rdatahfNI[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr_noInv.txt")


#####
library(parallel)
mat.obs <- as.matrix(hatFst2)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst2[,1]),(hatFst2[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr2.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFst4)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst4[,1]),(hatFst4[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr4.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFst11)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst11[,1]),(hatFst11[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hCrh11.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFst16)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst16[,1]),(hatFst16[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr16.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFstNI)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFstNI[,1]),(hatFstNI[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr_noInv.txt")



library("adegenet")
library("hierfstat")
library("vcfR")

setwd("/piec2/cgalia/Styela_mundial/2bRAD_world/pvalue_FST_NCSC")

metadata <- read.table("sampleinfo.txt", sep = "\t", header = TRUE)
colnames(metadata) <- c("Indv", "Hap", "Pop", "Region")

vcf2 <- read.vcfR("2_2mac_5-172DP_70_2bRAD_nou_NoNCSC.recode.vcf")
vcf4 <- read.vcfR("4_2mac_5-172DP_70_2bRAD_nou_NoNCSC.recode.vcf")
vcf11 <- read.vcfR("11_2mac_5-172DP_70_2bRAD_nou_NoNCSC.recode.vcf")
vcf16 <- read.vcfR("16_2mac_5-172DP_70_2bRAD_nou_NoNCSC.recode.vcf")
vcfNI <- read.vcfR("noInv_2mac_5-172DP_70_2bRAD_nou_NoNCSC.recode.vcf")

     
sdata2<- vcfR2genind(vcf2)
sdata4<- vcfR2genind(vcf4)
sdata11<- vcfR2genind(vcf11)
sdata16<- vcfR2genind(vcf16)
sdataNI<- vcfR2genind(vcfNI)

rdata2 <- sdata2
rdata4 <- sdata4
rdata11 <- sdata11
rdata16 <- sdata16
rdataNI <- sdataNI

hdata2 <- sdata2
hdata4 <- sdata4
hdata11 <- sdata11
hdata16 <- sdata16
hdataNI <- sdataNI

pop(sdata2) <- metadata$Pop
pop(sdata4) <- metadata$Pop
pop(sdata11) <- metadata$Pop
pop(sdata16) <- metadata$Pop
pop(sdataNI) <- metadata$Pop

pop(rdata2) <- metadata$Region
pop(rdata4) <- metadata$Region
pop(rdata11) <- metadata$Region
pop(rdata16) <- metadata$Region
pop(rdataNI) <- metadata$Region

pop(hdata2) <- metadata$Hap
pop(hdata4) <- metadata$Hap
pop(hdata11) <- metadata$Hap
pop(hdata16) <- metadata$Hap
pop(hdataNI) <- metadata$Hap

sdatahf2<-genind2hierfstat(sdata2,pop=NULL)
sdatahf4<-genind2hierfstat(sdata4,pop=NULL)
sdatahf11<-genind2hierfstat(sdata11,pop=NULL)
sdatahf16<-genind2hierfstat(sdata16,pop=NULL)
sdatahfNI<-genind2hierfstat(sdataNI,pop=NULL)

rdatahf2<-genind2hierfstat(rdata2,pop=NULL)
rdatahf4<-genind2hierfstat(rdata4,pop=NULL)
rdatahf11<-genind2hierfstat(rdata11,pop=NULL)
rdatahf16<-genind2hierfstat(rdata16,pop=NULL)
rdatahfNI<-genind2hierfstat(rdataNI,pop=NULL)

hdatahf2<-genind2hierfstat(hdata2,pop=NULL)
hdatahf4<-genind2hierfstat(hdata4,pop=NULL)
hdatahf11<-genind2hierfstat(hdata11,pop=NULL)
hdatahf16<-genind2hierfstat(hdata16,pop=NULL)
hdatahfNI<-genind2hierfstat(hdataNI,pop=NULL)

popFst2<-pairwise.neifst(sdatahf2,diploid=TRUE)
popFst4<-pairwise.neifst(sdatahf4,diploid=TRUE)
popFst11<-pairwise.neifst(sdatahf11,diploid=TRUE)
popFst16<-pairwise.neifst(sdatahf16,diploid=TRUE)
popFstNI<-pairwise.neifst(sdatahfNI,diploid=TRUE)

ratFst2<-pairwise.neifst(rdatahf2,diploid=TRUE)
ratFst4<-pairwise.neifst(rdatahf4,diploid=TRUE)
ratFst11<-pairwise.neifst(rdatahf11,diploid=TRUE)
ratFst16<-pairwise.neifst(rdatahf16,diploid=TRUE)
ratFstNI<-pairwise.neifst(rdatahfNI,diploid=TRUE)

hatFst2<-pairwise.neifst(hdatahf2,diploid=TRUE)
hatFst4<-pairwise.neifst(hdatahf4,diploid=TRUE)
hatFst11<-pairwise.neifst(hdatahf11,diploid=TRUE)
hatFst16<-pairwise.neifst(hdatahf16,diploid=TRUE)
hatFstNI<-pairwise.neifst(hdatahfNI,diploid=TRUE)

####p-valor FST####
library(parallel)
mat.obs <- as.matrix(popFst2)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf2[,1]),(sdatahf2[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr2.txt")


#### ####
library(parallel)
mat.obs <- as.matrix(popFst4)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf4[,1]),(sdatahf4[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr4.txt")

#####
library(parallel)
mat.obs <- as.matrix(popFst11)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf11[,1]),(sdatahf11[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr11.txt")

#####
library(parallel)
mat.obs <- as.matrix(popFst16)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahf16[,1]),(sdatahf16[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr16.txt")

#####
library(parallel)
mat.obs <- as.matrix(popFstNI)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(sdatahfNI[,1]),(sdatahfNI[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_pChr_noInv.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst2)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf2[,1]),(rdatahf2[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr2.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst4)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf4[,1]),(rdatahf4[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr4.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst11)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf11[,1]),(rdatahf11[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr11.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFst16)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahf16[,1]),(rdatahf16[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr16.txt")


#####
library(parallel)
mat.obs <- as.matrix(ratFstNI)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(rdatahfNI[,1]),(rdatahfNI[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=16,ncol=16)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_rChr_noInv.txt")


#####
library(parallel)
mat.obs <- as.matrix(hatFst2)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst2[,1]),(hatFst2[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr2.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFst4)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst4[,1]),(hatFst4[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr4.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFst11)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst11[,1]),(hatFst11[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hCrh11.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFst16)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFst16[,1]),(hatFst16[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr16.txt")



#####
library(parallel)
mat.obs <- as.matrix(hatFstNI)
NBPERM <- 1000 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- mclapply(1:NBPERM, function(i) genet.dist(cbind(sample(hatFstNI[,1]),(hatFstNI[,-1])),method = "WC84"))

library(ade4)
allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) as.matrix(mat.perm[[k]])[i,j])), mat.obs[i,j], alter="greater")
  }
}
pvals <- matrix(data=NA,nrow=11,ncol=11)
x=0
for (i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    x=x+1
    #     pvals[i,j] <- allTests[[x]][[6]]
    pvals[i,j] <- allTests[[x]]$pvalue
  }
}
pvals[lower.tri(pvals)] <- NA
fst_pval <- mat.obs
fst_pval[upper.tri(fst_pval)] <- pvals[upper.tri(pvals)]
diag(fst_pval) <- NA
write(x=fst_pval,file="Fst-pvalues_hChr_noInv.txt")





