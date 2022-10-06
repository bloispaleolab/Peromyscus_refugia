####################################################################################
####################### Range Expansion Peromyscus populations #####################
####################################################################################

#installing Range expansion package 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("snpStats")
library("snpStats")


devtools::install_github("BenjaminPeter/rangeExpansion", ref="package")
library(rangeExpansion)
library(sp)


wd <- "."
setwd(wd)

############ running it on the P. keeni

ploidy <- 2 #diploid individuals
region <- list("REGION_1")

#load in data 
snp.file <- "data/PK.snapp"
coord.file <- "data/Pk_locs.csv" 

raw_data <- load.data.snapp(snp.file, coord.file, sep = ",", ploidy=ploidy) 

#we calculate the population-level data from individual data and calculate all pairwise statistics:
pop <- make.pop(raw_data, ploidy)
psi <- get.all.psi(pop)

#Find origin 
Pk_results <- run.regions(region=region, pop=pop, psi=psi, xlen=10, ylen=20)

#Not wroking on newer versions of packages 
#summary(PM_results)
#plot(PM_results)

#Work around to get results 
summary.origin.results(Pk_results$tbl[[1]])

#create data table of range expansion locs 
north <- summary.origin.results(Pk_results$tbl[[1]])
north


############ running it on the P. m gambelii

ploidy <- 2 #diploid individuals
region <- list("REGION_1")

#load in data 
snp.file <- "PMG.snapp"
coord.file <- "PMG_locs.csv" 

raw_data <- load.data.snapp(snp.file, coord.file, sep = ",", ploidy=ploidy) 

#we calculate the population-level data from individual data and calculate all pairwise statistics:
pop <- make.pop(raw_data, ploidy)
psi <- get.all.psi(pop)

#Find origin 
PMG_results <- run.regions(region=region, pop=pop, psi=psi, xlen=10, ylen=20)

#Not wroking on newer versions of packages 
#summary(PM_results)
#plot(PM_results)

#Work around to get results 
summary.origin.results(PMG_results$tbl[[1]])

#create data table of range expansion locs 
scali <- summary.origin.results(PMG_results$tbl[[1]])
scali


############ running it on the P. m sonoriensis

ploidy <- 2 #diploid individuals
region <- list("REGION_1", "REGION_2", "REGION_3")

#load in data 
snp.file <- "PMS.snapp"
coord.file <- "PMS_locs.csv" 

raw_data <- load.data.snapp(snp.file, coord.file, sep = ",", ploidy=ploidy) 

#we calculate the population-level data from individual data and calculate all pairwise statistics:
pop <- make.pop(raw_data, ploidy)
psi <- get.all.psi(pop)

#Find origin 
PMS_results <- run.regions(region=region, pop=pop, psi=psi, xlen=10, ylen=20)

#Not wroking on newer versions of packages 
#summary(PM_results)
#plot(PM_results)

#Work around to get results 
summary.origin.results(PMS_results$tbl[[1]])
summary.origin.results(PMS_results$tbl[[2]])
summary.origin.results(PMS_results$tbl[[3]])
#create data table of range expansion locs 
PMS_1 <- summary.origin.results(PMS_results$tbl[[1]])
PMS_2 <- summary.origin.results(PMS_results$tbl[[2]])
PMS_3 <- summary.origin.results(PMS_results$tbl[[3]])


range_expan <- rbind (PMS_1, PMS_2, PMS_3, scali, north)
range_expan

write.csv(range_expan, "PM_range_expansion.csv")
