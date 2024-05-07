# PREVIOUSLY, REMEMBER TO "UNPHASE" SNPs THAT HAVE BEEN PHASED BY ACCIDENT!!! IF NOT, IT WILL CAUSE PROBLEMS AFTERWARDS

# Important steps to modify: 	a) To optine pure/hybrid individuals --> steps 4 & 9
#				b) To optine the level of confidence --> steps 6,7,8 & 11. If you want to keep more SNPs, for ex, you can do: diff > 0.85 (85% or more difference between allele frequencies of pure populations).
#				c) To rename or modify plots --> steps 16 & 19

#--------------------------------------------------------------------------------------------------------------

# Load all required packages
#install.packages("introgress")
#install.packages("Rtools")
#install.packages("introgress", dependencies=TRUE, repos="https://cran.r-project.org/")
#install.packages("introgress", type="source")
#install.packages("remotes")
#library(remote)
#install.packages("devtools")
#remotes::install_github("cran/introgress")

library(devtools)
library(ggplot2)
library(nnet)
library(genetics)
library(combinat)
library(gdata)
library(gtools)
library(MASS)
library(mvtnorm)
library(RColorBrewer)
library(vcfR)
library(adegenet)
library(ggsci)
library(scales)
library(ggpubr)
library(introgress)
library(ggrepel)

# Read in vcf as vcfR and population info

setwd("E:/styela/Molecular/Genomica poblacions/DNA mundial/2b_mundial/Coemncem a mirar pel paper/Colonització/introgress")
vcf <- read.vcfR("Chr_noInv_2mac_5-172DP_70_2bRAD_nou.recode.vcf")

pops <- as.data.frame(read.table("sampleinfo_intro.txt", stringsAsFactors = T, sep="\t", header=T))

# 1. optine colour palette for latter plots (?)

mypal <- c("steelblue", "black", "indianred", "grey")
show_col(mypal)

#---------------------------------------------------------------------------------------------------------------

# Start introgress analysis using the whole dataset 
# 2. Create SNP matrices
mat <- extract.gt(vcf)

conv.mat <- mat
conv.mat[conv.mat == "0/0"]<-0
conv.mat[conv.mat == "0/1"]<-1
conv.mat[conv.mat == "1/1"]<-2
conv.mat<-as.data.frame(conv.mat)

# 3. Convert to numeric
for (i in 1:ncol(conv.mat)){
  conv.mat[,i]<-as.numeric(as.character(conv.mat[,i]))
}

# 4. Calc AF for the samples you will use to call fixed differences, as well as the potential hybrids - refer to individuals by columns in both the numerator and denominator

nc.af <- (rowSums(conv.mat[,c(47:51)], na.rm=T)/(rowSums(is.na(conv.mat[,c(47:51)]) == FALSE)))/2
sc.af <- (rowSums(conv.mat[,c(67:71)], na.rm=T)/(rowSums(is.na(conv.mat[,c(67:71)]) == FALSE)))/2
atl.af <- (rowSums(conv.mat[,c(1,7:22,28:30,36:40,72:82)], na.rm=T)/(rowSums(is.na(conv.mat[,c(1,7:22,28:30,36:40,72:82)]) == FALSE)))/2
pac.af <- (rowSums(conv.mat[,c(2:6,23:27,31:35,41:46,52:56,62:66,83:87)], na.rm=T)/(rowSums(is.na(conv.mat[,c(2:6,23:27,31:35,41:46,52:56,62:66,83:87)]) == FALSE)))/2
atpa.af <- (rowSums(conv.mat[,c(1:46,52:66,72:87)], na.rm=T)/(rowSums(is.na(conv.mat[,c(1:46,52:66,72:87)]) == FALSE)))/2

# 5. Find fixed SNPs
#NCAT <- abs(nc.af - atl.af)
NCAT <- abs(nc.af - pac.af)
#NCAT <- abs(nc.af - atpa.af)


# 6. How many SNPs are fixed
table(is.na(NCAT) == FALSE & NCAT >= 0.9)
vcf@fix[,1][is.na(NCAT) == FALSE & NCAT >= 0.9]

# 7. Subsample original matrix to only fixed diff SNPs
gen.mat <- mat[is.na(NCAT) == FALSE & NCAT >= 0.9,]
dim(gen.mat)

# 8. Subsample matrix converted for AF calcs to only fixed SNPS
conv.mat<-conv.mat[is.na(NCAT) == FALSE & NCAT >= 0.9,]
dim(conv.mat)

# 9. Write a logical test to convert alleles so that a single number represents ONE PARENTAL ANCESTRY (one of the pure set of indvs.)
for (i in 1:nrow(gen.mat)){
  #if 1 is the NC allele
  if((sum(conv.mat[i,c(47:51)], na.rm=TRUE)/(sum(is.na(conv.mat[i,c(47:51)]) == FALSE)))/2 == 0){
    #swap all '0/0' cells with '2/2'
    gen.mat[i,][gen.mat[i,] == "0/0"] <- "2/2"
    #swap all '1/1' cells with '0/0'
    gen.mat[i,][gen.mat[i,] == "1/1"] <- "0/0"
    #finally convert all '2/2' cells (originally 0/0) into '1/1'
    gen.mat[i,][gen.mat[i,] == "2/2"] <- "1/1"
    #no need to touch hets
  }
}

# 10. Convert R class NAs to the string "NA/NA"
gen.mat[is.na(gen.mat) == TRUE] <- "NA/NA"

# 11. Make locus info df

locus.info <- data.frame(locus=rownames(gen.mat),
                       type=rep("C", times=nrow(gen.mat)),
                       lg=vcf@fix[,1][is.na(NCAT) == FALSE & NCAT >= 0.9],
                       marker.pos=vcf@fix[,2][is.na(NCAT) == FALSE & NCAT >= 0.9])

# 12 Make linkage group numeric

locus.info$lg <- gsub("Chromosome_", "", locus.info$lg)
locus.info$lg <- as.numeric(locus.info$lg)
locus.info$marker.pos <- as.numeric(as.character(locus.info$marker.pos))

# 13. Make bpcum
nCHR <- length(unique(locus.info$lg))
locus.info$BPcum <- NA
s <- 0
nbp <- c()
for (i in sort(unique(locus.info$lg))){
  nbp[i] <- max(locus.info[locus.info$lg == i,]$marker.pos)
  locus.info[locus.info$lg == i,"BPcum"] <- locus.info[locus.info$lg == i,"marker.pos"] + s
  s <- s + nbp[i]
}

# --------------------------------------------------------------------------------------------------------------------

# We now have a gt matrix in proper format for introgress

# 14. Convert genotype data into a matrix of allele counts
count.matrix <- prepare.data(admix.gen=gen.mat, loci.data=locus.info,
                           parental1="1",parental2="0", pop.id=F,
                           ind.id=F, fixed=T)

# 15. Estimate hybrid index values
hi.index.sim <- est.h(introgress.data=count.matrix,loci.data=locus.info,
                    fixed=T, p1.allele="1", p2.allele="0")

locus.info$locus<-rep("", times=nrow(locus.info))

# 16. Create "local ancestry" plot

mk.image(introgress.data=count.matrix, loci.data=locus.info,
         marker.order=order(locus.info$BPcum),hi.index=hi.index.sim, ylab.image="Individuals",
         xlab.h="population 2 ancestry", pdf=T, out.file="anc.NCvsPac.pdf",
         col.image=c(rgb(1,0,0,alpha=.6),rgb(0,0,0,alpha=.9),rgb(0,0,1,alpha=.5)))


# 17. Calculate mean heterozygosity across fixed markers for each sample
het <- calc.intersp.het(introgress.data=count.matrix)

# 18. Make triangle plot
triangle.plot(hi.index=hi.index.sim, int.het=het, pdf = T)

# 19. Create a triangle plot and merge dataframes

d_tr_plot <- merge(pops$region, hi.index.sim, by=0)
rownames(d_tr_plot) <- d_tr_plot$Row.names
d_tr_plot <- d_tr_plot[,-1] 
d_tr_plot <- d_tr_plot[order(as.numeric(row.names(d_tr_plot))), ]
d_tr_plot$intersp_het <- het
names<-pops$sample
d_tr_plot <- d_tr_plot[which(d_tr_plot$x!="atl"),]

p_triangle_all <- ggplot(d_tr_plot, aes(x=h, y=intersp_het, colour=x)) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  theme(axis.text = element_text(color="black", size=14),
        axis.text.x=element_text(angle = 90, vjust=0.5),
        axis.title = element_text(size=14)) +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = +2)
  ) +
  theme(legend.position="none",
        legend.title = element_text(size = 18, face = "bold"), legend.text=element_text(size=14),legend.key.height=unit(1.0,'cm')) +
  geom_segment(x = 0, y = 0, xend = 0.5, yend = 1, colour="black", size=0.1) +
  geom_segment(x = 0.5, y = 1, xend = 1, yend = 0, colour="black", size=0.1) +
  geom_segment(x = 0, y = 0, xend = 1, yend = 0, colour="black", size=0.1) +
  geom_point(size=3) +
  scale_color_manual(values=c( "black", "indianred","grey")) +
  #coord_fixed(ratio = 1) +
  xlim(0,1) + ylim(0,1) +
  labs(x="Introgression index", y="Heterozygosity") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) 
  #scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))
  #geom_text_repel(aes(x=h, y=intersp_het, label=names), colour = "black", na.rm = TRUE, size = 3, max.overlaps = getOption("ggrepel.max.overlaps", default = 6), min.segment.length = 0.25, box.padding = 0.75)  

ggsave("triangle.NCPac.pdf",plot=last_plot(),device="pdf", path=getwd(), height = 3 , width = 4)
  
png("triangle.lep_mci_ben_r70p3_opt_LmB.png", width=2000, height=750, res=300)
p_triangle_all
dev.off()
