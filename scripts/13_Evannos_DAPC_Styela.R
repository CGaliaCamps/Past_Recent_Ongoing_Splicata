library(ggplot2)

setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Pop structure/DAPC")

si <- read.table("evannoPC_siNCSC.txt", header=T)
no <- read.table("evannoPC_noNCSC.txt", header=T)


si$Delta <-round(si$Delta, digits=1)
no$Delta <-round(no$Delta, digits=1)

pdf("Evanno_PC_si.pdf", width=6, height=3)
ggplot(data=si, aes(x=PCs, y=Delta^2)) +
  geom_line(lwd=1.5, color="navy") + 
  geom_point(size=2) +  
  #geom_text(label=si$Delta, vjust = -0.75, size=3) +
  ylab("Delta BIC") +
  xlab("Number of clusters") + 
  theme_classic() +
  theme(axis.text.x=element_text(size=10, face = "bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        axis.title.x =element_text(size=15, face = "bold"),
        axis.title.y =element_text(size=15, face = "bold"))
dev.off()

pdf("Evanno_PC_no.pdf", width=6, height=3)
ggplot(data=no, aes(x=PCs, y=Delta^2)) +
  geom_line(lwd=1.5, color="navy") + 
  geom_point(size=2) +  
  #geom_text(label=no$Delta, vjust = -0.75, size=3) +
  ylab("Delta BIC") +
  xlab("Number of clusters") + 
  theme_classic() +
  theme(axis.text.x=element_text(size=10, face = "bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        axis.title.x =element_text(size=15, face = "bold"),
        axis.title.y =element_text(size=15, face = "bold"))
dev.off()

pdf("BIC_PC_si.pdf", width=6, height=3)
ggplot(data=si, aes(x=PCs, y=BIC)) +
  geom_line(lwd=1, color="grey30") +
  geom_point(size=5,shape=21,stroke=2,fill="white") +  
  geom_errorbar(aes(ymin=BIC-SD, ymax=BIC+SD), width=0.5,
                position=position_dodge(1), color="red3", lwd=1)+
  ylab("BIC") +
  xlab("Number of clusters") + 
  theme_classic() +
  theme(axis.text.x=element_text(size=10, face = "bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        axis.title.x =element_text(size=15, face = "bold"),
        axis.title.y =element_text(size=15, face = "bold"))
dev.off()

pdf("BIC_PC_no.pdf", width=6, height=3)
ggplot(data=no, aes(x=PCs, y=BIC)) +
  geom_line(lwd=1, color="grey30") +
  geom_point(size=5,shape=21,stroke=2,fill="white") +  
  geom_errorbar(aes(ymin=BIC-SD, ymax=BIC+SD), width=0.5,
                position=position_dodge(1), color="red3", lwd=1)+
  ylab("BIC") +
  xlab("Number of clusters") + 
  theme_classic() +
  theme(axis.text.x=element_text(size=10, face = "bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        axis.title.x =element_text(size=15, face = "bold"),
        axis.title.y =element_text(size=15, face = "bold"))
dev.off()






