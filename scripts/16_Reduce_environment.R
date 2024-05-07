
##Load packages

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
if (!require("vegan")) {
  install.packages("vegan", dependencies = TRUE)
  library(vegan)
} 
if (!require("poppr")) {
  install.packages("poppr", dependencies = TRUE)
  library(poppr)
}
if (!require("pegas")) {
  install.packages("pegas", dependencies = TRUE)
  library(pegas)
}
if (!require("psych")) {
  install.packages("psych", dependencies = TRUE)
  library(psych)
}


setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Adaptacio/Dades ambientals")

env <- read.table("all_env_data.txt", sep = "\t", header = TRUE)
metadata <- read.csv("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/sampleinfo.txt", sep = "\t", header = TRUE)
row.names(env)<-env$pop

#check correlation values among variables
pdf(file="correlacions_environment.pdf", width=10, height=10)
pairs.panels(env[4:22], scale=F, density=T, stars=F, ellipse=F, hist.col="grey80", 
             lm=F, smooth=T, digits=1, rug=T, cex=1.2, cex.cor=2.5, smoother=T, method="pearson") 
dev.off()

####Select only non-correlated variables####

env3 <- subset(env, select=c(pop,meanto,meanso,rangeto,rangeso,chl,fe,no3,si,o2,po4,spco2
                          ))

pdf(file="correlacions_environment_3.pdf", width=10, height=10)
pairs.panels(env3[2:12], scale=F, density=T, stars=F, ellipse=F, hist.col="grey80", 
             lm=F, smooth=T, digits=1, rug=T, cex=1.2, cex.cor=2.5, smoother=T, method="spearman") 
dev.off()


#get the loading of each of your variables to the first PCA axes.
env_PCA <- princomp(env3[,2:12], cor = T)
env_PCA_plot <- as.data.frame(env_PCA$scores)
PCA_loadings <- env_PCA$loadings[1:11,1:11]
PCA_loadings <- as.data.frame(PCA_loadings)
library(factoextra)
eig.val <- get_eigenvalue(env_PCA)


#####Discard non-important varibales (do not appear on the first 3 axes)####

env4 <- subset(env3, select=-c(chl,fe,no3,po4
))

colnames(env4) <- c("pop", "Mean T","Mean S","Range T","Range S","Si","O2","CO2")

write.table(env4,file="environmental_final.txt")  


pdf(file="correlacions_environment_4.pdf", width=10, height=10)
pairs.panels(env4[2:8], scale=F, density=T, stars=F, ellipse=F, hist.col="grey80", 
             lm=F, smooth=T, digits=1, rug=T, cex=1.2, cex.cor=1.5, smoother=T) 
dev.off()



env_PCA2 <- princomp(env4[,2:8], cor = F)

env_PCA2_plot <- as.data.frame(env_PCA2$scores)
PCA2_loadings <- env_PCA2$loadings[1:7,1:7]
PCA2_loadings <- as.data.frame(PCA2_loadings)
library(factoextra)
eig.val <- get_eigenvalue(env_PCA2)


###Plot

colnames(env_PCA_plot) <- c("Comp.1", "Comp.2", "Comp3", "Comp4", "Comp5", "Comp6", "Comp7", "Comp8", "Comp9","Comp10","Comp11")
colnames(PCA_loadings) <- c("x.Comp.1", "x.Comp.2", "Comp3", "Comp4", "Comp5", "Comp6", "Comp7", "Comp8", "Comp9","Comp10","Comp11")

mycol <- c("#9FCAE6","#9C0824","#559F52","#B3E0A6","#5F90BB", #AGA,AU,BAR,BLA,BRA
           "#EAC0BD","#B9DDF1","#AD0B1E","#3F6E9A","#EC6756", #CA,FE,HK,KN,MIS
           "#000000","#CC1619","#2A5783","#F8968C","grey50",  #NC,OKI,PE,SAK,SC
           "#7CADD2","#24693D","#DA3E2F")                     #TE,VIL,WAK

pdf("PCA_environmental_finalreduction1_2.pdf", width = 6, height = 4)
ggplot(data = env_PCA2_plot, aes(x = Comp.1, y = Comp.2))+
    geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) +   
  geom_vline(xintercept = 0, lty = "dotted", color="grey", cex=1) +
  geom_point(data = env_PCA2_plot, size = 4, aes(color=rownames(env_PCA2_plot), fill=rownames(env_PCA2_plot)), stroke = 1.5)+
  geom_segment(data = PCA2_loadings, aes(x = 0, xend = 13*Comp.1, y = 0, yend = 13*Comp.2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1) +
  geom_text(data = PCA2_loadings, aes(x= 15*Comp.1, y = 15*Comp.2, label = rownames(PCA2_loadings)), 
            size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_manual(values=mycol)+
  scale_shape_manual(values=c(21)) +
  scale_color_manual(values=mycol)+
  theme_classic()+
  labs(x = paste("PCA1 (", round(eig.val[1,2], 2), "%)", sep = ""), y = paste("PCA2 (", round(eig.val[2,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

pdf("PCA_environmental_finalreduction3_4.pdf", width = 6, height = 4)
ggplot(data = env_PCA2_plot, aes(x = Comp.3, y = Comp.4))+
  geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) +   
  geom_vline(xintercept = 0, lty = "dotted", color="grey", cex=1) +
  geom_point(data = env_PCA2_plot, size = 4, aes(color=rownames(env_PCA2_plot), fill=rownames(env_PCA2_plot)), stroke = 1.5)+
  geom_segment(data = PCA2_loadings, aes(x = 0, xend = 13*Comp.3, y = 0, yend = 13*Comp.4),
               arrow = arrow(length = unit(0.3, "cm")), colour = "grey20", cex=1) +
  geom_text(data = PCA2_loadings, aes(x= 15*Comp.3, y = 15*Comp.4, label = rownames(PCA2_loadings)), 
            size = 5, hjust = 0.5)+
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_manual(values=mycol)+
  scale_shape_manual(values=c(21)) +
  scale_color_manual(values=mycol)+
  theme_classic()+
  labs(x = paste("PCA3 (", round(eig.val[3,2], 2), "%)", sep = ""), y = paste("PCA4 (", round(eig.val[4,2], 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))
dev.off()

