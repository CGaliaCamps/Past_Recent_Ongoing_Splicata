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

setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Adaptacio/Fora inveriso/RDA")

###Read files####
gnl <- read.PLINK("Chr_noInv_2mac_5-172DP_70_2bRAD_nou.plink.ped.raw", parallel = F)
gn <- as.data.frame(gnl)
env4 <- read.table("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Adaptacio/Dades ambientals/reduccio_dimensions/environmental_final_sample.txt", sep = "\t",header=T)
env4 <- env4[,c(3:9)]
metadata <- read.csv("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/sampleinfo.txt", sep = "\t", header = TRUE)

####impute NA####
sum(is.na(gn))
gen.imp <- apply(gn, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp))

sty.rda <- rda(gen.imp~.,data=env4,scale=T)#vector=marcadors+variables

summary(sty.rda)## Check the accumulated constrained eigenvalues


summary(eigenvals(sty.rda, model = "constrained")) #can see the % explained by each axis
pdf("Chr_noInv_eigenvalues_RDA.pdf")
screeplot(sty.rda) #plot
dev.off()

####Test if the RDA test is significative overall and for each axis. This is slow####
 signif.full <- anova.cca(sty.rda, parallel=getOption("mc.cores")) # default is permutation=999 PER veure si el model es significatiu
 write.table(signif.full,file="Chr_noInv_model_anova_fullaxes_RDA.txt", sep="\t")
# signif.axis <- anova.cca(sty.rda, by="axis", parallel=getOption("mc.cores")) # PER veure si cada eix del model es significatiu
# write.table(signif.axis,file="Chr_noInv_model_anova_eachaxis_RDA.txt", sep="\t")
# #vif.cca(sty.rda) #


####Plot RDA ####
rdaplot<-summary(sty.rda)
arrows <- as.data.frame(rdaplot$biplot) #Get arrow information
write.table(arrows,file="Chr_noInv_envs_asso_RDA.txt", sep="\t")
rdaplot <- as.data.frame(rdaplot$sites) # Get individuals information
rdaplot<- cbind(rdaplot, metadata$pop) # Assign pop to individual

head(rdaplot)

colnames(rdaplot) <- c("RDA1","RDA2","RDA3","RDA4","RDA5","RDA6","Pop") # Change RDA columns name
row.names(arrows) <- c("Mean T","Mean S","Range T","Range S","Si","O2","CO2") #Change row names to enviromental variables

##set your colors
mycol <- c("#9FCAE6","#9C0824","#559F52","#B3E0A6","#5F90BB", #AGA,AU,BAR,BLA,BRA
                    "#EAC0BD","#B9DDF1","#AD0B1E","#3F6E9A","#EC6756", #CA,FE,HK,KN,MIS
                    "#000000","#CC1619","#2A5783","#F8968C","grey50",  #NC,OKI,PE,SAK,SC
                    "#7CADD2","#24693D","#DA3E2F")                     #TE,VIL,WAK


#pcorrelative plots for first 4 axes
pdf("Chr_noInv_RDA_1vs2.pdf", width = 7, height = 7)
ggplot(data = rdaplot, aes(x = RDA1, y = RDA2))+
  geom_point(data = rdaplot, pch = 21, size = 8, aes(fill = Pop, shape=Pop, alpha=1))+
  geom_segment(data = arrows, aes(x = 0, xend = 4*RDA1, y = 0, yend = 4*RDA2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "steelblue", cex=1) +
  geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) + geom_vline(xintercept = 0, 
                                                                               lty = "dotted", color="grey", cex=1) +
  geom_text(data = arrows, aes(x= 4.1*RDA1, y = 4.1*RDA2, label = rownames(arrows)), 
            size = 7, hjust = 0.5)+
  #geom_text(data = prdav, aes(x = RDA1, y = RDA2, label = rownames(prdav)), 
  #         size = 2.5, col = "black", hjust = 1.2)+
  scale_fill_manual(values=mycol, labels=c("AGA","AU","BAR","BLA","BRA","CA","FE","HK","KN","MIS","NC","OKI","PE","SAK","SC","TE","VIL","WAK"))+
  theme_classic()+
  labs(x = paste("RDA1 (", round(summary(sty.rda)$cont$importance[2,1]*100, 2), "%)", sep = ""), y = paste("RDA2 (", round(summary(sty.rda)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
dev.off()

# pdf("Chr_noInv_RDA_2vs3.pdf", width = 4, height = 4)
# ggplot(data = rdaplot, aes(x = RDA2, y = RDA3))+
#   geom_point(data = rdaplot, pch = 21, size = 8, aes(fill = Pop, shape=Pop, alpha=1))+
#   geom_segment(data = arrows, aes(x = 0, xend = 4*RDA2, y = 0, yend = 4*RDA3),
#                arrow = arrow(length = unit(0.3, "cm")), colour = "steelblue", cex=1) +
#   geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) + geom_vline(xintercept = 0, 
#                                                                                lty = "dotted", color="grey", cex=1) +
#   geom_text(data = arrows, aes(x= 4.1*RDA2, y = 4.1*RDA3, label = rownames(arrows)), 
#             size = 5, hjust = 0.5)+
#   #geom_text(data = prdav, aes(x = RDA1, y = RDA2, label = rownames(prdav)), 
#   #         size = 2.5, col = "black", hjust = 1.2)+
#   scale_fill_manual(values=mycol, labels=c("AGA","AU","BAR","BLA","BRA","CA","FE","HK","KN","MIS","NC","OKI","PE","SAK","SC","TE","VIL","WAK"))+
#   theme_classic()+
#   labs(x = paste("RDA2 (", round(summary(sty.rda)$cont$importance[2,2]*100, 2), "%)", sep = ""), 
#        y = paste("RDA3 (", round(summary(sty.rda)$cont$importance[2,3]*100, 2), "%)", sep = ""))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
# dev.off()
# 
# pdf("Chr_noInv_RDA_3vs4.pdf", width = 6, height = 3)
# ggplot(data = rdaplot, aes(x = RDA3, y = RDA4))+
#   geom_point(data = rdaplot, pch = 21, size = 8, aes(fill = Pop, shape=Pop, alpha=1))+
#   geom_segment(data = arrows, aes(x = 0, xend = 4*RDA3, y = 0, yend = 4*RDA4),
#                arrow = arrow(length = unit(0.3, "cm")), colour = "steelblue", cex=1) +
#   geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) + geom_vline(xintercept = 0, 
#                                                                                lty = "dotted", color="grey", cex=1) +
#   geom_text(data = arrows, aes(x= 4.1*RDA3, y = 4.1*RDA4, label = rownames(arrows)), 
#             size = 5, hjust = 0.5)+
#   #geom_text(data = prdav, aes(x = RDA1, y = RDA2, label = rownames(prdav)), 
#   #         size = 2.5, col = "black", hjust = 1.2)+
#   scale_fill_manual(values=mycol, labels=c("AGA","AU","BAR","BLA","BRA","CA","FE","HK","KN","MIS","NC","OKI","PE","SAK","SC","TE","VIL","WAK"))+
#   theme_classic()+
#   labs(x = paste("RDA3 (", round(summary(sty.rda)$cont$importance[2,3]*100, 2), "%)", sep = ""), 
#        y = paste("RDA4 (", round(summary(sty.rda)$cont$importance[2,4]*100, 2), "%)", sep = ""))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
# dev.off()

# ####Identify candidate SNPs####
# load.rda <- scores(sty.rda, choices=c(1:4), display="species")
# 
# hist(load.rda[,1], main="Loadings on RDA1") #cand de l'eix 1
# hist(load.rda[,2], main="Loadings on RDA2") #cand de l'eix 2
# hist(load.rda[,3], main="Loadings on RDA3") #cand de l'eix 3
# hist(load.rda[,4], main="Loadings on RDA4") #cand de l'eix 4
# 
# outliers <- function(x,z){
#   lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
#   x[x < lims[1] | x > lims[2]]               # locus names in these tails
# }
# 
# cand1 <- outliers(load.rda[,1],3) # number of outliers for each axis
# cand2 <- outliers(load.rda[,2],3) # 
# cand3 <- outliers(load.rda[,3],3) # 
# cand4 <- outliers(load.rda[,4],3) # 
# 
# ncand <- length(cand1) + length(cand2) + length(cand3) + length(cand4)
# ncand
# 
# cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1)) #Define SNPs of each axis
# cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
# cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
# cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
# 
# colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4)<- c("axis","snp","loading")
# 
# cand <- rbind(cand1, cand2, cand3, cand4) #cconcatenate candidate SNPs
# cand$snp <- as.character(cand$snp)
# 
# foo <- matrix(nrow=(ncand), ncol=7)  # 8 columns for 8 predictors
# colnames(foo) <- c("meanto","meanso","rangeto","rangeso","si","o2","spco2")
# 
# for (i in 1:length(cand$snp)) {
#   nam <- cand[i,2]
#   snp.gen <- gen.imp[,nam]
#   foo[i,] <- apply(env4,2,function(x) cor(x,snp.gen))
# }
# 
# cand <- cbind.data.frame(cand,foo)  
# head(cand)
# cand <- cand[!duplicated(cand$snp),] #get rid of duplicate candidates 
# 
# #### It gives to you a table with each candidate SNP, the environmental variable correlated with, and dthe correlation value####
# for (i in 1:length(cand$snp)) { 
#   bar <- cand[i,]
#   cand[i,11] <- names(which.max(abs(bar[4:10]))) # gives the variable
#   cand[i,12] <- max(abs(bar[4:10]))              # gives the correlation
# }
# 
# colnames(cand)[11] <- "predictor"
# colnames(cand)[12] <- "correlation"
# 
# table(cand$predictor) 
# write.table(cand,file="Chr_noInv_Candidate_4eixos.txt", sep="\t")
# 
# ####set the file to plot####
# sel <- cand$snp
# env <- cand$predictor
# #env[env=="meanto"] <- 'orangered4'
#  # env[env=="meanso"] <- '#1f78b4'
#   #  env[env=="rangeto"] <- 'orangered1'
#    #   env[env=="rangeso"] <- '#a6cee3'
#     #    env[env=="si"] <- '#33a02c'
#      #     env[env=="o2"] <- 'grey25'
#       #      env[env=="spco2"] <- 'grey50'
#           
#           # color by predictor:
#           col.pred <- rownames(sty.rda$CCA$v) # pull the SNP names
#           
#           for (i in 1:length(sel)) {           # color code candidate SNPs
#             foo <- match(sel[i],col.pred)
#             col.pred[foo] <- env[i]
#           }
#           
#           col.pred[grep("Chr",col.pred)] <- 'non-candidate' # non-candidate SNPs
#             empty <- col.pred
#             empty[grep("white",empty)] <- rgb(0,1,0, alpha=0) # transparent
#             empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","grey")
#             bg <- c('blue','red',"grey95",'black','#1f78b4','orange','purple2','grey75') #c("Mean_sal","Mean_Temp","non-candidate","CO2","Range_Sal","Range_Temp","Si","O2"))+
#             bg_col <- c('black','black',"NA",'black','black','black','black','black')
#             
#           colors <- as.data.frame(col.pred)
#           col.pred <- as.vector(col.pred)
# candplot <- summary(sty.rda) 
# candplot <- as.data.frame(candplot$species)  
# candplot<- cbind(candplot, colors)
# head(candplot)
# colnames(candplot) <- c("RDA1","RDA2","RDA3","RDA4","RDA5","RDA6","cand")
# 
# order <- c("meanto", "rangeto", "meanso", "rangeso", "spco2", "o2","si")
# 
# top <- candplot[grep("non-candidate", candplot$cand),]
# bottom <- candplot[-grep("non-candidate", candplot$cand),]
# bottom <- bottom[order(match(bottom$cand,order)),]
# candplot <- rbind(top, bottom) 
# 
# tail(candplot)
# 
# # 
# pdf("Chr_noInv_RDA_1vs2_candidate.pdf", width = 7, height = 7)     
# ggplot(data = candplot, aes(x = RDA1, y = RDA2))+
#   geom_point(data = candplot, pch = 21, size = 4, aes(color=candplot$cand, fill=candplot$cand, alpha=1))+
#   geom_segment(data = arrows, aes(x = 0, xend = 0.4*RDA1, y = 0, yend = 0.4*RDA2),
#                arrow = arrow(length = unit(0.3, "cm")), colour = "steelblue", cex=1) +
#   geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) + geom_vline(xintercept = 0, 
#                                                                                lty = "dotted", color="grey", cex=1) +
#   geom_text(data = arrows, aes(x= 0.35*RDA1, y = 0.35*RDA2, label = rownames(arrows)), 
#             size = 5, hjust = 0.5)+
#   #geom_text(data = candplot, aes(x = RDA1, y = RDA2, label = rownames(candplot)), 
#   #        size = 2.5, col = "black", hjust = 1.2)+
#   scale_fill_manual(values=bg,labels=c("Mean_sal","Mean_Temp","non-candidate","CO2","Range_Sal","Mean_Temp","Si","O2"))+
#   scale_color_manual(values=bg_col)+
#   theme_classic()+
#   labs(x = paste("RDA1 (", round(summary(sty.rda)$cont$importance[2,1]*100, 2), "%)", sep = ""), 
#        y = paste("RDA2 (", round(summary(sty.rda)$cont$importance[2,2]*100, 2), "%)", sep = ""))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
#   theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
# 
# dev.off()
#             
#             
# pdf("Chr_noInv_RDA_2vs3_candidate.pdf", width = 7, height = 7)
#             ggplot(data = candplot, aes(x = RDA2, y = RDA3))+
#               geom_point(data = candplot, pch = 21, size = 4, aes(color=candplot$cand, fill=candplot$cand, alpha=1))+
#               geom_segment(data = arrows, aes(x = 0, xend = 0.4*RDA2, y = 0, yend = 0.4*RDA3),
#                            arrow = arrow(length = unit(0.3, "cm")), colour = "steelblue", cex=1) +
#               geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) + geom_vline(xintercept = 0, 
#                                                                                            lty = "dotted", color="grey", cex=1) +
#               geom_text(data = arrows, aes(x= 0.35*RDA2, y = 0.35*RDA3, label = rownames(arrows)), 
#                         size = 5, hjust = 0.5)+
#               #geom_text(data = candplot, aes(x = RDA1, y = RDA2, label = rownames(candplot)), 
#               #        size = 2.5, col = "black", hjust = 1.2)+
#               scale_fill_manual(values=bg,labels=c("Mean_sal","Mean_Temp","non-candidate","CO2","Range_Sal","Mean_Temp","Si","O2"))+
#               scale_color_manual(values=bg_col)+
#               theme_classic()+
#               labs(x = paste("RDA2 (", round(summary(sty.rda)$cont$importance[2,2]*100, 2), "%)", sep = ""), 
#                    y = paste("RDA3 (", round(summary(sty.rda)$cont$importance[2,3]*100, 2), "%)", sep = ""))+
#               theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#               theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
#               theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
#             dev.off()
#             
# pdf("Chr_noInv_RDA_3vs4_candidate.pdf", width = 7, height = 7)
#       ggplot(data = candplot, aes(x = RDA3, y = RDA4))+
#         geom_point(data = candplot, pch = 21, size = 4, aes(color=candplot$cand, fill=candplot$cand, alpha=1))+
#         geom_segment(data = arrows, aes(x = 0, xend = 0.4*RDA3, y = 0, yend = 0.4*RDA4),
#                      arrow = arrow(length = unit(0.3, "cm")), colour = "steelblue", cex=1) +
#         geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) + geom_vline(xintercept = 0, 
#                      lty = "dotted", color="grey", cex=1) +
#         geom_text(data = arrows, aes(x= 0.35*RDA3, y = 0.35*RDA4, label = rownames(arrows)), 
#                      size = 5, hjust = 0.5)+
#        #geom_text(data = candplot, aes(x = RDA1, y = RDA2, label = rownames(candplot)), 
#             #        size = 2.5, col = "black", hjust = 1.2)+
#         scale_fill_manual(values=bg,labels=c("Mean_sal","Mean_Temp","non-candidate","CO2","Range_Sal","Range_Temp","Si","O2"))+
#         scale_color_manual(values=bg_col)+
#          theme_classic()+
#                      labs(x = paste("RDA3 (", round(summary(sty.rda)$cont$importance[2,3]*100, 2), "%)", sep = ""), 
#                      y = paste("RDA4 (", round(summary(sty.rda)$cont$importance[2,4]*100, 2), "%)", sep = ""))+
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#         theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
#         theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
# dev.off()           
#            
# 
# candplot$cand
              
  
