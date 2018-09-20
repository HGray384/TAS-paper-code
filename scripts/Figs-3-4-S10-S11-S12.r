#############################################################################
#
# Title: script data-based simulation
#
# Authors: Harry Gray, Gwenael Leday and Catalina Vallejos
#
# Date: Sep 20, 2018
#
#############################################################################




#############################################################################
#############################################################################
### PRELIMINARIES

# clean environment
rm(list=ls());gc()

# Set the wd to the base github folder
base.dir <-
setwd(base.dir)

# model-based simulations
library("TAS")# library(devtools);install_github("HGray384/TAS")
library("corpcor")
library("ShrinkCovMat")
library("abind")
library("data.table")
library("ggplot2")
source('./scripts/scm_un.R')
# install.packages("cgdsr")
library("cgdsr")
#source("https://bioconductor.org/biocLite.R");biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")


#############################################################################
#############################################################################
## Obtain gene symbols for all genes within the p53 and apoptosis pathways

# Obtain data bases linking gene symbols to entrez ID and KEGG pathways
maploc <- toTable(org.Hs.egPATH)
maploc2 <- toTable(org.Hs.egSYMBOL)
head(maploc)
head(maploc2)

# Obtain gene symbols for p53 and apoptosis pathways
p53entrez <- maploc$gene_id[maploc$path_id=="04115"] # Entrez ids of genes in p53 signaling pathway 
p53symbols <- maploc2$symbol[match(p53entrez, maploc2$gene_id)]
apoptentrez <- maploc$gene_id[maploc$path_id=="04210"] # Entrez ids of genes in apoptosis pathway 
apoptsymbols <- maploc2$symbol[match(apoptentrez, maploc2$gene_id)]

p53 <- data.frame("entrezID" = p53entrez, "symbol" = p53symbols, stringsAsFactors=FALSE)
apoptosis <- data.frame("entrezID" = apoptentrez, "symbol" = apoptsymbols, stringsAsFactors=FALSE)

head(p53)
head(apoptosis)

dim(p53)
dim(apoptosis)

#############################################################################
#############################################################################
## Download TCGA gene expression data

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/")

# Identify cancer types for which microarray data are available
f <- function(x){substr(x, nchar(x)-4, nchar(x))}
ind <- sapply(getCancerStudies(mycgds)[,1], f) == "_tcga"
allcancers <- sort(getCancerStudies(mycgds)[ind,1])

# Download for each cancer type gene expression of genes in p53 and apoptosis pathways
TCGAp53 <- rep(list(NULL), length(allcancers))
TCGAapopt <- rep(list(NULL), length(allcancers))
names(TCGAp53) <- names(TCGAapopt) <- unlist(lapply(strsplit(allcancers, "_"), function(x){x[1]}))
f2 <- function(x){substr(x, nchar(x)-8, nchar(x))}
for(i in 1:length(allcancers)){
  tp <- getCaseLists(mycgds, allcancers[i])[,c(1:2)]
  ind <- sapply(tp[,1], f2) == "tcga_mrna"
  if(sum(ind)>0){
    print(tp[ind,])
    
    # Get expression of genes in p53 pathway
    TCGAp53[[i]] <- data.matrix(getProfileData(mycgds, genes=p53$symbol, geneticProfiles=tp[ind,1], caseList=tp[ind,1]))
    
    # Get expression of genes in apoptosis pathway
    TCGAapopt[[i]] <- data.matrix(getProfileData(mycgds, genes=apoptosis$symbol, geneticProfiles=tp[ind,1], caseList=tp[ind,1]))
  }
}

# Remove cancer types for which there was no microarray data
TCGAp53 <- TCGAp53[!unlist(lapply(TCGAp53, is.null))]
TCGAapopt <- TCGAapopt[!unlist(lapply(TCGAapopt, is.null))]

# No missing values in p53 data, all good
lapply(TCGAp53, function(x){any(is.na(x))})
lapply(TCGAp53, dim)

# Missing values present for IL3RA gene in apoptosis data
lapply(TCGAapopt, function(x){any(is.na(x))})
lapply(TCGAapopt, function(x){sum(is.na(x))})
TCGAapopt <- lapply(TCGAapopt, function(x){x[,!colnames(x)%in%"IL3RA"]})
lapply(TCGAapopt, dim)

# Check same genes in present in data matrices
all(table(unlist(lapply(TCGAp53, colnames)))==length(TCGAp53))
all(table(unlist(lapply(TCGAapopt, colnames)))==length(TCGAapopt))


#############################################################################
#############################################################################
### SIMULATIONS BREAST CANCER DATA

# Center the data
mat <- t(scale(TCGAp53$brca, scale=FALSE, center = TRUE)) 
dim(mat)
extmat <- lapply(TCGAp53[!names(TCGAp53)%in%"brca"], scale, center = TRUE, scale=FALSE)
length(extmat) # number of additional cancers
# Transpose to have the format required for TAS
extmat <- lapply(extmat, t)

# look at sample covariance matrix
sam <- tcrossprod(mat)/ncol(mat)
tcond <- kappa(sam, exact=T)

# get the value of p
p <- nrow(mat)

# set the high-dimensional values of n
n0 <- round(c(3*p/4, p/2, p/4))

# number of reps
reps <- 1000

# number of columns: useful for sampling
N <- 1:ncol(mat)

# set shrinkage intensity grid
alpha <- seq(0.01, 0.99, 0.01)

# initialise variables to store condition numbers and frobenius losses
frobsBRCA <- array(0, c(8, reps, length(n0))) # 3 TAS + 3 scm + cpc + MLE
weightsBRCA <- array(0, c(9, length(alpha), reps, length(n0)))
weightsBRCA2 <- array(0, c(10, length(alpha), reps, length(n0)))
weightsBRCA3 <- array(0, c(18, length(alpha), reps, length(n0)))

# Function to return regularized covariance estimate
varcorcodes <- which(matrix(TRUE, 3, 3), arr.ind=TRUE)
regS <- function(y, varcorcodes){
	thelistTargets <- sapply(1:nrow(varcorcodes), function(x){getTarget(y, var=varcorcodes[x, 1], cor=varcorcodes[x, 2])}, simplify=FALSE)
	thearrayTargets <- simplify2array(thelistTargets)
	taShrink(y, targets = thearrayTargets, without = 0, alpha = alpha, plots = FALSE)$sigmahat
}

# Get data-based target on pooled data
datPooled <- Reduce("cbind", extmat)
dim(datPooled)
dataBasedTargetPooled <- regS(datPooled, varcorcodes)

# Obtained target from each cancer type (regularized sample covariance)
dataBasedTargetsList <- lapply(extmat, regS, varcorcodes = varcorcodes)
dataBasedTargetsArray <- simplify2array(dataBasedTargetsList)

# loop
for (i in 1:length(n0)){
  # sample columns in a reproducible way (different seed for each n)
  set.seed(reps+n0[i])
  n0cols <- replicate(reps, sample(N, n0[i])) # each column contains indices of columns to sample from mat
  for (r in 1:reps){
    cat("n0=", n0[i], "\t", "rep=", r, " of ", reps, "\n", sep="")
    
    # get the partitioned matrices
    hd <- mat[,n0cols[,r]]
    truemat <- mat[,-c(n0cols[,r])]
    
    # get the sample covariances
    hdsam <- tcrossprod(hd)/ncol(hd)
    truesam <- tcrossprod(truemat)/ncol(truemat)

    ## TAS with 9 targets
    listTargets <- sapply(1:nrow(varcorcodes), function(x){getTarget(hd, var=varcorcodes[x, 1], cor=varcorcodes[x, 2])}, simplify=FALSE)
    arrayTargets <- simplify2array(listTargets) # include 9 default targets
    tas <- taShrink(hd, targets = arrayTargets, without = 0, alpha = alpha, plots = FALSE)
    weightsBRCA[,,r, i] <- tas$weights

    ## TAS with 9 default targets + 1 data-based target (regularized pooled sample covariance of external data)
    arrayTargets2 <- abind(arrayTargets, dataBasedTargetPooled)
    tas2 <- taShrink(hd, targets = arrayTargets2, without = 0, alpha = alpha, plots = FALSE)
    weightsBRCA2[,,r, i] <- tas2$weights

    ## TAS with 9 default targets + 9 data-based targets (regularized sample covariances in each cancer)
    arrayTargets3 <- abind(arrayTargets, dataBasedTargetsArray)
    tas3 <- taShrink(hd, targets = arrayTargets3, without = 0, alpha = alpha, plots = FALSE)
    weightsBRCA3[,,r, i] <- tas3$weights

    ## shrinkcovmat identity
    at1 <- shrinkcovmat.identity(hd, centered = TRUE)

    ## shrinkcovmat equal
    at2 <- shrinkcovmat.equal(hd, centered = TRUE)

    ## shrinkcovmat unequal
    at3 <- scm.un(hd, centered = TRUE)
    
    ## corpcor
    cpc <- cov.shrink(t(hd), verbose = FALSE)
    
    ## Frobenius losses
    frobsBRCA[1, r, i] <- norm(truesam-tas$sigmahat, type = "F")
    frobsBRCA[2, r, i] <- norm(truesam-tas2$sigmahat, type = "F")
    frobsBRCA[3, r, i] <- norm(truesam-tas3$sigmahat, type = "F")
    frobsBRCA[4, r, i] <- norm(truesam-at1$Sigmahat, type = "F")
    frobsBRCA[5, r, i] <- norm(truesam-at2$Sigmahat, type = "F")
    frobsBRCA[6, r, i] <- norm(truesam-at3$Sigmahat, type = "F")
    frobsBRCA[7, r, i] <- norm(truesam-cpc, type = "F")
    frobsBRCA[8, r, i] <- norm(truesam-hdsam, type = "F")
  }
}


#############################################################################
#############################################################################
### SIMULATIONS OVARIAN CANCER DATA

# Center the data
mat <- t(scale(TCGAapopt$ov, scale=FALSE, center = TRUE)) 
dim(mat)
extmat <- lapply(TCGAapopt[!names(TCGAapopt)%in%"ov"], scale, center = TRUE, scale=FALSE)
length(extmat) # number of additional cancers
# Transpose to have the format required for TAS
extmat <- lapply(extmat, t)

# look at sample covariance matrix
sam <- tcrossprod(mat)/ncol(mat)
tcond <- kappa(sam, exact=T)

# get the value of p
p <- nrow(mat)

# set the high-dimensional values of n
n0 <- round(c(3*p/4, p/2, p/4))

# number of reps
reps <- 1000

# number of columns: useful for sampling
N <- 1:ncol(mat)

# set shrinkage intensity grid
alpha <- seq(0.01, 0.99, 0.01)

# initialise variables to store condition numbers and frobenius losses
frobsOV <- array(0, c(8, reps, length(n0))) # TAS + 3 scm + cpc + MLE
weightsOV <- array(0, c(9, length(alpha), reps, length(n0)))
weightsOV2 <- array(0, c(10, length(alpha), reps, length(n0)))
weightsOV3 <- array(0, c(18, length(alpha), reps, length(n0)))

# Get data-based target on pooled data
datPooled <- Reduce("cbind", extmat)
dim(datPooled)
varcorcodes <- which(matrix(TRUE, 3, 3), arr.ind=TRUE)
dataBasedTargetPooled <- regS(datPooled, varcorcodes)

# Obtained target from each cancer type (regularized sample covariance)
dataBasedTargetsList <- lapply(extmat, regS, varcorcodes)
dataBasedTargetsArray <- simplify2array(dataBasedTargetsList)

# loop
for (i in 1:length(n0)){
  # sample columns in a reproducible way (different seed for each n)
  set.seed(reps+n0[i])
  n0cols <- replicate(reps, sample(N, n0[i])) # each column contains indices of columns to sample from mat
  for (r in 1:reps){
    cat("n0=", n0[i], "\t", "rep=", r, "of", reps, "\n", sep = " ")
    
    # get the partitioned matrices
    hd <- mat[,n0cols[,r]]
    truemat <- mat[,-c(n0cols[,r])]
    
    # get the sample covariances
    hdsam <- tcrossprod(hd)/ncol(hd)
    truesam <- tcrossprod(truemat)/ncol(truemat)

    ## TAS with 9 targets
    listTargets <- sapply(1:nrow(varcorcodes), function(x){getTarget(hd, var=varcorcodes[x, 1], cor=varcorcodes[x, 2])}, simplify=FALSE)
    arrayTargets <- simplify2array(listTargets) # include 9 default targets
    tas <- taShrink(hd, targets = arrayTargets, without = 0, alpha = alpha, plots = FALSE)
    weightsOV[,,r, i] <- tas$weights

    ## TAS with 9 default targets + 1 data-based target (regularized pooled sample covariance of external data)
    arrayTargets2 <- abind(arrayTargets, dataBasedTargetPooled)
    tas2 <- taShrink(hd, targets = arrayTargets2, without = 0, alpha = alpha, plots = FALSE)
    weightsOV2[,,r, i] <- tas2$weights

    ## TAS with 9 default targets + 9 data-based targets (regularized sample covariances in each cancer)
    arrayTargets3 <- abind(arrayTargets, dataBasedTargetsArray)
    tas3 <- taShrink(hd, targets = arrayTargets3, without = 0, alpha = alpha, plots = FALSE)
    weightsOV3[,,r, i] <- tas3$weights

    ## shrinkcovmat identity
    at1 <- shrinkcovmat.identity(hd, centered = TRUE)

    ## shrinkcovmat equal
    at2 <- shrinkcovmat.equal(hd, centered = TRUE)

    ## shrinkcovmat unequal
    at3 <- scm.un(hd, centered = TRUE)
    
    ## corpcor
    cpc <- cov.shrink(t(hd), verbose = FALSE)
    
    ## Frobenius losses
    frobsOV[1, r, i] <- norm(truesam-tas$sigmahat, type = "F")
    frobsOV[2, r, i] <- norm(truesam-tas2$sigmahat, type = "F")
    frobsOV[3, r, i] <- norm(truesam-tas3$sigmahat, type = "F")
    frobsOV[4, r, i] <- norm(truesam-at1$Sigmahat, type = "F")
    frobsOV[5, r, i] <- norm(truesam-at2$Sigmahat, type = "F")
    frobsOV[6, r, i] <- norm(truesam-at3$Sigmahat, type = "F")
    frobsOV[7, r, i] <- norm(truesam-cpc, type = "F")
    frobsOV[8, r, i] <- norm(truesam-hdsam, type = "F")
    
  }
}


#############################################################################
#############################################################################
### PLOT RESULTS

# PRIAL function
prial <- function(frobs, baseline){
  (sum(frobs[baseline, ])-rowSums(frobs[-baseline,]))/sum(frobs[baseline, ])*100
}

# Function used to personalise barplot 
mybarplot <- function(labs, vals){
  tp <- data.frame("labels"=labs, "scenario"=vals)
  tp$labels <- factor(tp$labels, levels = labs)
  myggplot <- ggplot(tp, aes(x=tp$labels, y=tp$scenario)) +
    geom_bar(stat="identity") +
    theme_bw() +
    ylim(min(c(0,vals)), 100) +
    xlab("") + ylab("PRIAL") + 
    theme(axis.title.y = element_text(size = rel(1.8)),
          axis.text.y = element_text(size = rel(1.6)),
          axis.text.x = element_text(size = rel(1.4)))
  myggplot
}

# Function used to personalise barplot (TAS)
myboxplot <- function(vals){
  tp <- data.table(t(vals))
  names(tp) <- c(paste0("T", seq_len(9)), "S")
  tp$Iteration <- seq_len(nrow(tp))
  tp2 <- data.table::melt(tp, id.vars = c("Iteration"),
                          measure.vars = c(paste0("T", seq_len(9)), "S"))
  
  myggplot <- ggplot(tp2, aes(x = variable, y = value)) + geom_boxplot() + 
    theme_classic() + ylim(0, 1) + xlab("Target") + ylab("Posterior weight") +
    theme(axis.title.y = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.text.y = element_text(size = rel(1.6)),
          axis.text.x = element_text(size = rel(1.4)))
  myggplot
}

# Function used to personalise barplot (TAS-info)
myboxplot2 <- function(vals){
  tp <- data.table(t(vals))
  names(tp) <- c(paste0("T", seq_len(9)), "ext", "S")
  tp$Iteration <- seq_len(nrow(tp))
  tp2 <- data.table::melt(tp, id.vars = c("Iteration"),
                          measure.vars = c(paste0("T", seq_len(9)), "ext", "S"))
  
  myggplot <- ggplot(tp2, aes(x = variable, y = value)) + geom_boxplot() + 
    theme_classic() + ylim(0, 1) + xlab("Target") + ylab("Posterior weight") +
    theme(axis.title.y = element_text(size = rel(1.8)),
          axis.title.x = element_text(size = rel(1.8)),
          axis.text.y = element_text(size = rel(1.6)),
          axis.text.x = element_text(size = rel(1.4)))
  myggplot
}

##############################
#### Plots for BRCA data
##############################

# get p and n0
mat <- t(scale(TCGAp53$brca, scale=FALSE, center = TRUE)) 
p <- nrow(mat)
n0 <- round(c(3*p/4, p/2, p/4))

##-------------##
## Plot PRIAL
##-------------##

# Compute PRIAL for BRCA
frobsBRCAsq <- frobsBRCA^2 # squared loss
prial1 <- prial(frobsBRCAsq[,,1], dim(frobsBRCAsq)[1])[-3]
prial2 <- prial(frobsBRCAsq[,,2], dim(frobsBRCAsq)[1])[-3]
prial3 <- prial(frobsBRCAsq[,,3], dim(frobsBRCAsq)[1])[-3]

# Plot PRIAL for BRCA
nm <- c("TAS", "TAS-info", "AT1", "AT2", "AT3", "CPC")

par(mar=c(2,4,1,2)+0.1)
barplot(rbind(prial1, prial2, prial3), beside=TRUE, ylim=c(0, 100), names.arg=nm, ylab="PRIAL", xlab="", main="")
legend("topright", legend=c("n=3p/4", "n=p/2", "n=p/4"), fill=gray.colors(3), cex=1.2)
box()


##----------------------------------------##
## Plot TAS target weights
##----------------------------------------##

# plots of target-specific shrinkage intensities
tweights <- apply(weightsBRCA, c(1, 3, 4), function(x){x*alpha})
tweights <- apply(tweights, c(2, 3, 4), sum)
sweights <- 1-apply(tweights, c(2, 3), sum)

# Plot and save results - weights
for(ii in 1:3){
  assign(paste0("boxplot_target_weights_BRCA_n",n0[ii]),
         myboxplot(vals = rbind(tweights[,,ii], sweights[,ii]) ))
  plot(get(paste0("boxplot_target_weights_BRCA_n",n0[ii])))
}

##--------------------------------------------##
## Plot TAS-info target weights
##--------------------------------------------##

# plots of target-specific shrinkage intensities
tweights2 <- apply(weightsBRCA2, c(1, 3, 4), function(x){x*alpha})
tweights2 <- apply(tweights2, c(2, 3, 4), sum)
sweights2 <- 1-apply(tweights2, c(2, 3), sum)

# Plot and save results - weights
for(ii in 1:3){
  assign(paste0("boxplot_target_weights2_BRCA_n",n0[ii]),
         myboxplot2(vals = rbind(tweights2[,,ii], sweights2[,ii]) ))
  plot(get(paste0("boxplot_target_weights2_BRCA_n",n0[ii])))
}

##############################
#### Plots for OV data
##############################

# get p and n0
mat <- t(scale(TCGAapopt$ov, scale=FALSE, center = TRUE)) 
p <- nrow(mat)
n0 <- round(c(3*p/4, p/2, p/4))

##-------------##
## Plot PRIAL
##-------------##

# Compute PRIAL for OV
frobsOVsq <- frobsOV^2 # squared loss
prial1 <- prial(frobsOVsq[,,1], dim(frobsOVsq)[1])[-3]
prial2 <- prial(frobsOVsq[,,2], dim(frobsOVsq)[1])[-3]
prial3 <- prial(frobsOVsq[,,3], dim(frobsOVsq)[1])[-3]

# Plot PRIAL for OV
nm <- c("TAS", "TAS-info", "AT1", "AT2", "AT3", "CPC")
par(mar=c(2,4,1,2)+0.1)
barplot(rbind(prial1, prial2, prial3), beside=TRUE, ylim=c(0, 100), names.arg=nm, ylab="PRIAL", xlab="", main="")
legend("topright", legend=c("n=3p/4", "n=p/2", "n=p/4"), fill=gray.colors(3), cex=1.2)
box()


##----------------------------------------##
## Plot TAS target weights
##----------------------------------------##

# plots of target-specific shrinkage intensities
tweights <- apply(weightsOV, c(1, 3, 4), function(x){x*alpha})
tweights <- apply(tweights, c(2, 3, 4), sum)
sweights <- 1-apply(tweights, c(2, 3), sum)

# Plot and save results - weights
for(ii in 1:3){
  assign(paste0("boxplot_target_weights_OV_n",n0[ii]),
         myboxplot(vals = rbind(tweights[,,ii], sweights[,ii]) ))
  plot(get(paste0("boxplot_target_weights_OV_n",n0[ii])))
}

##--------------------------------------------##
## Plot TAS-info target weights
##--------------------------------------------##

# plots of target-specific shrinkage intensities
tweights2 <- apply(weightsOV2, c(1, 3, 4), function(x){x*alpha})
tweights2 <- apply(tweights2, c(2, 3, 4), sum)
sweights2 <- 1-apply(tweights2, c(2, 3), sum)

# Plot and save results - weights
for(ii in 1:3){
  assign(paste0("boxplot_target_weights2_OV_n",n0[ii]),
         myboxplot2(vals = rbind(tweights2[,,ii], sweights2[,ii]) ))
  plot(get(paste0("boxplot_target_weights2_OV_n",n0[ii])))
}

