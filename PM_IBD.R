###########
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)

#Load in vcf file for %85 of data 
wd <- "."
setwd(wd)

#load in vcf file 
pm_75 <- read.vcfR("data/PM_final_all.recode.vcf")
head(pm_75)               #check the vcf object
pm_75@fix[1:10,1:5]       #check


#load in pop information 
pm_pop <- read.csv("data/PM_Seq_locality.csv", header=T)

### convert to genlight
PM_75.genlight <- vcfR2genlight(pm_75, n.cores=2)
PM_75.genlight

#add pop names
pop(PM_75.genlight)<- pm_pop[,2]

#Calculation and visualization of Nei’s distances
### Calculate Nei's distances between individuals/pops
PM.75.ind <- stamppNeisD(PM_75.genlight, pop = FALSE) # Nei's 1972 distance between indivs
PM.75.pop <- stamppNeisD(PM_75.genlight, pop = TRUE) # Nei's 1972 distance between pops

### Calculate pairwise Fst among populations
PM_75.genlight@ploidy <- as.integer(ploidy(PM_75.genlight))

PM.fst <- stamppFst(PM_75.genlight, nboots = 1, percent =95, nclusters=3)
PM.fst

### Isolation by distance
coords <- read.csv ("data/PM_pop_grps.csv") # tab-separated file for all pops
xy.coords.only<- subset(coords, select=c("Latitude","Longitude"))
Dgeo <- dist(xy.coords.only)

## Calculate distance in km 
DistMat <- rdist.earth(x1 = coords[2:3], miles = FALSE)
DistMat <- as.dist(DistMat)

# create the dist objects used in analyses below
colnames(PM.75.ind) <- rownames(PM.75.ind)
PM.75.ind.dist<-as.dist(PM.75.ind, diag=T)
attr(PM.75.ind.dist, "Labels")<-rownames(PM.75.ind) # name the rows of a matrix


#test IBD 
IBD_PM <- mantel.randtest(DistMat, PM.75.ind.dist)
IBD_PM
plot(DistMat, PM.75.ind.dist, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM.75.ind.dist~DistMat))
title("Correlation of Genetic and Geographical distances (p-value=0.001)")

#plot (come back tothis later)
pdf("Figures/PM_IBD.pdf")
par(mfrow = c(2, 2))

#all
plot(DistMat, PM.75.ind.dist, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM.75.ind.dist~DistMat))
title("Correlation of Genetic and Geographical distances (p-value=0.001)", cex.main=0.7)

#sor
plot(DistMat_sor, PM_75.dist_sor, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM_75.dist_sor~DistMat_sor))
title("Correlation of Genetic and Geographical distances (p-value=0.001)", cex.main=0.7)

#gam
plot(DistMat_gam, PM.75.dist_gam, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM.75.dist_gam~DistMat_gam))
title("Correlation of Genetic and Geographical distances (p-value=0.001)", cex.main=0.7)

#keeni
plot(DistMat_k, PM.75.dist_k, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM.75.dist_k~DistMat_k))
title("Correlation of Genetic and Geographical distances (p-value=0.001)", cex.main=0.7)
dev.off()
      

#############################
##IBD for PM sor pop 
###plot AFS per one pop
PM_75.genlight.sep <- seppop(PM_75.genlight, drop=TRUE)

#check if it worked 
PM_75.genlight.sep$`1`  
### Calculate Nei's distances between individuals/pops for 85
PM_75_sor <- stamppNeisD(PM_75.genlight.sep$`1`, pop = FALSE) # Nei's 1972 distance between indivs
colnames(PM_75_sor) <- rownames(PM_75_sor)
PM_75.dist_sor<-as.dist(PM_75_sor, diag=T)
attr(PM_75.dist_sor, "Labels")<-rownames(PM_75_sor) # name the rows of a matrix

Sor <- read.csv("data/PMS_gen.csv")[2:3]
Dgeo_sor <- dist(Sor)

## Calculate distance in km 
DistMat_sor <- rdist.earth(x1 = Sor, miles = FALSE)
DistMat_sor <- as.dist(DistMat_sor)

#test IBD
IBD_75.sor <- mantel.randtest(DistMat_sor, PM_75.dist_sor)
IBD_75.sor
plot(DistMat_sor, PM_75.dist_sor, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM_75.dist_sor~DistMat_sor))
title("Correlation of Genetic and Geographical distances (p-value=0.001)")

#####IBD for PMG
PM_75.genlight.sep$`2`  
### Calculate Nei's distances between individuals/pops
PM.75_gam <- stamppNeisD(PM_75.genlight.sep$`2`, pop = FALSE) # Nei's 1972 distance between indivs
colnames(PM.75_gam) <- rownames(PM.75_gam)
PM.75.dist_gam<-as.dist(PM.75_gam, diag=T)
attr(PM.75.dist_gam, "Labels")<-rownames(PM.75.dist_gam) # name the rows of a matrix

gam <- read.csv("data/PMG_gen.csv")[2:3]
Dgeo_gam <- dist(gam)

## Calculate distance in km 
DistMat_gam <- rdist.earth(x1 = gam, miles = FALSE)
DistMat_gam <- as.dist(DistMat_gam)


#test IBD
IBD_75.gam <- mantel.randtest(DistMat_gam, PM.75.dist_gam)
IBD_75.gam
plot(DistMat_gam, PM.75.dist_gam, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM.75.dist_gam~DistMat_gam))
title("Correlation of Genetic and Geographical distances (p-value=0.001)")


#####IBD for PK
PM_75.genlight.sep$`3`  
### Calculate Nei's distances between individuals/pops
PM.75_k <- stamppNeisD(PM_75.genlight.sep$`3`, pop = FALSE) # Nei's 1972 distance between indivs
colnames(PM.75_k) <- rownames(PM.75_k)
PM.75.dist_k<-as.dist(PM.75_k, diag=T)
attr(PM.75.dist_k, "Labels")<-rownames(PM.75.dist_k) # name the rows of a matrix

k <- read.csv("data/Pk_gen.csv")[2:3]
Dgeo_k <- dist(k)

## Calculate distance in km 
DistMat_k <- rdist.earth(x1 = k, miles = FALSE)
DistMat_k <- as.dist(DistMat_k)

#test IBD
IBD_75.k <- mantel.randtest(DistMat_k, PM.75.dist_k)
IBD_75.k
plot(DistMat_k, PM.75.dist_k, pch=20,cex=.5, xlab="Geographic Distance (km)", 
     ylab="Genetic Distance")
abline(lm(PM.75.dist_k~DistMat_k))
title("Correlation of Genetic and Geographical distances (p-value=0.001)")

















