###########Prelim analyses of PM
library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)

#Load in vcf file for %90 of data 
wd <- "."
setwd(wd)
pm_all <- read.vcfR("data/PM_final_all.recode.vcf")
head(pm_all)               #check the vcf object
pm_all@fix[1:10,1:5]       #check


#quick check read depth distribution per individual
pdf("read_depth_per_ind.pdf")
dp <- extract.gt(pm_all, element='DP', as.numeric=TRUE)
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.4, cex.axis=0.5)
dev.off()

#zoom to smaller values
pdf("PM_DP_RAD_data.pdf", width = 10, height=3) # boxplot
par(mar=c(8,4,1,1))
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.4, cex.axis=0.5, ylim=c(0,50))
abline(h=8, col="red")
dev.off()

#load in pop information 
pm_pop <- read.csv("data/PM_Seq_locality.csv", header=T)
### convert to genlight
PM_all.genlight <- vcfR2genlight(pm_all, n.cores=4)

#add pop names
pop(PM_all.genlight)<- pm_pop[,5]

# check the genlight
PM_all.genlight                        # check the basic info on the genlightobject
indNames(PM_all.genlight)              # check individual names
as.matrix(PM_all.genlight)[1:16,1:10]  # see tiny bit of the data
pop(PM_all.genlight)                   # population assignment
# look at the total data matrix (0,1,2; white = missing data)
glPlot (PM_75.genlight)  # takes some time
# N missing SNPs per sample
x_75 <- summary(t(as.matrix(PM_75.genlight)))

###plot total AFS of the dataset
mySum_75 <- glSum(PM_75.genlight, alleleAsUnit = TRUE)
barplot(table(mySum_75), col="blue", space=0, xlab="Allele counts",
        main="Distribution of ALT allele counts in total dataset")

##PCA
toRemove_75 <- is.na(glMean(PM_all.genlight, alleleAsUnit = FALSE)) # TRUE where NA
which(toRemove_75) # position of entirely non-typed loci
b_75 <- PM_all.genlight[, !toRemove_75]
PM_pca_all <- glPca(b_75, nf=300, n.cores=4) # this should work

# proportion of explained variance by first three axes
PM_pca_all$eig[1]/sum(PM_pca_all$eig) # proportion of variation explained by 1st axis
PM_pca_all$eig[2]/sum(PM_pca_all$eig) # proportion of variation explained by 2nd axis
PM_pca_all$eig[3]/sum(PM_pca_all$eig) # proportion of variation explained by 3rd axis

# fast plot
scatter(PM_pca_all, posi="bottomright")

# save fig
pdf ("Figures/PM_all_PCA.pdf", width=14, height=7)
col <- c("yellow", "red", "blue")
g1 <- s.class(PM_pca_all$scores, pop(PM_all.genlight), xax=1, yax=2,
              col=transp(col,.6),
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T,
              pgrid.draw =F, plot = FALSE)
g2 <- s.label (PM_pca_all$scores, xax=1, yax=2, ppoints.col = "red", plabels =
                 list(box = list(draw = FALSE),
                      optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#DAPC
grp_PM <- find.clusters(PM_all.genlight, max.n.clust=10, glPca = PM_pca_all, method= "ward", 
                        perc.pca = 100, n.iter=1e7, n.start=1000)
dapc_PM <- dapc(PM_all.genlight, grp_PM$grp, glPca =PM_pca_all)

groups <- data.frame(grp_PM$grp)
col <- c("yellow", "blue", "red")

pdf ("Figures/PM_all_DAPC_plot.pdf", width=14, height=7)
scatter(dapc_PM, col = col)
dev.off()

pdf ("../Figures/PM_all_DAPC.pdf", width=14, height=7)
#par(mfrow = c(1, 2))
scatter(dapc_PM, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=col, solid=.4,
        cex=3,clab=0, leg=TRUE)
dev.off()

