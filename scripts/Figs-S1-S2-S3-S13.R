#############################################################################
#
# Title: supplementary material figures
#
# Authors: Harry Gray, Gwenael Leday and Catalina Vallejos
#
# Date: Sep 18, 2018
#
#############################################################################




#############################################################################
#############################################################################
### PRELIMINARIES

# clean environment
rm(list=ls());gc()

########################### MLE plot
####################################
library("abind")
library("MCMCpack")

# data dimensions
p <- seq(200, 1000, 200)

# generate data
mle_sim <- function()
{
  X <- lapply(p, function(x){lapply(c(x/0.1,x/0.5,x/1,x/2,x/10), 
                                    function(y){mvrnorm(n=y, mu=rep(0, x), Sigma=diag(x))})})
  S <- lapply(X, function(x){lapply(x, function(y){crossprod(scale(y, scale=FALSE, center = TRUE))/nrow(y)})})
  frobs <- lapply(S, function(x){lapply(x, function(y){norm(diag(ncol(y))-y, type = "F")^2})})
  condslog10 <- lapply(S, function(x){lapply(x, function(y){log10(kappa(y, exact=TRUE))})})
  rm(X)
  rm(S)
  return(list(frobs=simplify2array(frobs), log10conds=simplify2array(condslog10)))
}

# robustness
reps <- 100
# reproducibility 
set.seed(45)
result <- replicate(n = reps, mle_sim()) # takes about an hour

## create plots
# format frobenius losses and condition numbers
fl <- array(simplify2array(unlist(result[1,])), c(5, 5, 100))
cn <- array(simplify2array(unlist(result[2,])), c(5, 5, 100))

lab1 <- seq(1, 25, 6)
nlab <- sapply(lab1, function(x){return(c(x, x+1, x+2, x+3, x+4))})
nlab <- c(t(nlab))


rbc <- rep(rainbow(5), 5)
rbc <- sort(rbc)

# Frobenius loss
par(bty='n')
boxplot(t(log10(abind(fl[,1,], fl[,2,], fl[,3,], fl[,4,], fl[,5,], along=1))), at=nlab, col=t(rbc), 
        border=t(rbc),
        xaxt='n', xlab="n", ylab="Squared Frobenius loss (log-10)", main="Frobenius loss of S for identity covariance", 
        cex.lab=1.4, cex.main=1.2, outcex=0.3, pch=16)
axis(1,at=seq(3, 29, 6),labels=c("10p", "2p", "p", "p/2", "p/10"),cex.axis=1)
legend("topleft", inset=.05, title="p",
       legend=seq(200, 1000, 200), fill=sort(rainbow(5)), bty='n', cex=0.8)
mtext("a)", side=3, adj=-0.23, cex=1.5, line=1.5)

# log10 condition number
boxplot(t(abind(cn[,1,], cn[,2,], cn[,3,], cn[,4,], cn[,5,], along=1)), at=nlab, col=t(rbc), 
        border=t(rbc),
        xaxt='n', xlab="n", ylab="Condition number (log-10)", main="Condition number of S", 
        cex.lab=1.4, cex.main=1.2, outcex=0.3, pch=16)
abline(h=log10(1/.Machine$double.eps), col='grey', lwd=3)
axis(1,at=seq(3, 29, 6),labels=c("10p", "2p", "p", "p/2", "p/10"),cex.axis=1)
legend("topleft", inset=.05, title="p",
       legend=seq(200, 1000, 200), fill=sort(rainbow(5)), bty='n', cex=0.8)
mtext("b)", side=3, adj=-0.23, cex=1.5, line=1.5)


############################ uncertain marginal likelihood
##########################################################

# model-based simulations
library("TAS")
library("MCMCpack")
library("MASS")
library("Matrix")

# dimension values
p <- 200
n <- 20

# true covariance matrix with 0 correlations and unit variances
sigma <- 2*diag(p)
set.seed(10)
X <- mvrnorm(n=n, mu=rep(0, p), Sigma=sigma)
X <- t(scale(X, scale=FALSE, center = TRUE))

# target and alpha grid
target <- diag(p)
alpha <- seq(0.01, 0.99, 0.01)

# compute log-marginal likelihood
lml <- logML(X, target = target, alpha = alpha)
BF <- exp(lml-max(lml))
inds <- which(1/BF<40)

# create bayes factor plot
plot(1/BF[inds], x=c(inds)/100, ylab="Bayes Factor (BF)", xlab = expression(alpha),
     bty='n', type='l', lwd = 3, ylim=c(1, 40))
abline(h=3)
abline(h=20, lty=2)
legend("topleft", legend=c("BF=3", "BF=20"), lty=c(1, 2), cex=0.9, bty='n')

# get range of values with different evidence classifications
which(1/BF<3)
which(1/BF<20)
which(1/BF>3)[which(1/BF<20)]

######################################### grid of A versus PRIAL
################################################################

library("TAS")
library("MCMCpack")
library("MASS")
library("Matrix")

# grids to compare
alpha1 <- c(0.2, 0.4, 0.6, 0.8)
alpha2 <- seq(0.1, 0.9, 0.1)
alpha3 <- seq(0.05, 0.95, 0.05)
alpha4 <- seq(0.01, 0.99, 0.01)
alpha5 <- seq(0.005, 0.995, 0.005)
alpha6 <- seq(0.001, 0.999, 0.001)

# dimensions of data
p <- 100
n <- 25

# true covariance matrix with 0 correlations and equal variances
corrmat <- diag(p)
sdmat <- 2*diag(p)
sigma <- sdmat %*% corrmat %*% sdmat

# generate data
simu <- function()
{
  X <- lapply(n, function(x){mvrnorm(n=x, mu=rep(0, p), Sigma=sigma)})
  X <- lapply(X, function(x){t(scale(x, scale=FALSE, center = TRUE))})
  
  # perform shrinkage
  frobs <- lapply(X, 
                  function(x){
                    ta1 <- taShrink(x, alpha = alpha1, plots = FALSE)$sigmahat;
                    ta2 <- taShrink(x, alpha = alpha2, plots = FALSE)$sigmahat;
                    ta3 <- taShrink(x, alpha = alpha3, plots = FALSE)$sigmahat;
                    ta4 <- taShrink(x, alpha = alpha4, plots = FALSE)$sigmahat;
                    ta5 <- taShrink(x, alpha = alpha5, plots = FALSE)$sigmahat;
                    ta6 <- taShrink(x, alpha = alpha6, plots = FALSE)$sigmahat;
                    S <- tcrossprod(x)/ncol(x);
                    ests <- list(ta1, ta2, ta3, ta4, ta5, ta6, S)
                    return(lapply(ests, function(y){Matrix::norm(sigma-y, type = "F")^2}))
                  }
  )
  return(simplify2array(frobs))
}

set.seed(25)
losses <- replicate(100, simu())
losses <- matrix(unlist(losses), nrow=7, ncol=100)

# PRIAL function
prial <- function(frobs, baseline){
  (sum(frobs[baseline, ])-rowSums(frobs[-baseline,]))/sum(frobs[baseline, ])*100
}

# make plots
nm <- c("4", "9", "19", "99", "199", "999")
plot(prial(losses, 7), ylim=c(94, 100), xaxt='n', main="Equidistant sets for A", ylab="PRIAL", xlab="# of values in A", bty='n')
axis(1, at=1:6, labels = nm, las=1)

####################### multivariate normality of TCGA datasets
###############################################################

# Libraries
# install.packages("cgdsr")
library(cgdsr)
#source("https://bioconductor.org/biocLite.R");biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
library('MVN')

###### run the same processing script as in TCGA data files
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

# Test the CGDS endpoint URL using a few simple API tests
#test(mycgds)

# List of cancer studies
#getCancerStudies(mycgds)[,1:2]

# Data available for Breast cancer
#getCaseLists(mycgds, "brca_tcga")[,c(1:2)]

# Get available genetic profiles
#mygeneticprofile = getGeneticProfiles(mycgds, "brca_tcga_mrna")[1,1]

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

####### end processing

ov <- TCGAapopt$ov
ov <- scale(ov, scale = FALSE, center = TRUE)
brca <- TCGAp53$brca
brca <- scale(brca, scale = FALSE, center = TRUE)

# 3 different tests for ovarian - none reject the null
ovmardia <- mvn(ov, mvnTest = "mardia")
ovhz <- mvn(ov, mvnTest = "hz")
ovroyston <- mvn(ov, mvnTest = "royston")

# 3 different tests for breast - none reject the null
brcamardia <- mvn(brca, mvnTest = "mardia")
brcahz <- mvn(brca, mvnTest = "hz")
brcaroyston <- mvn(brca, mvnTest = "royston")

# qq plots
qqnorm(ov[,"AKT3"], main="AKT3 gene ovarian cancer")
qqline(ov[,"AKT3"], main="AKT3 gene ovarian cancer")
mtext("a)", side=3, adj=-0.23, cex=1.5, line=1.5)

qqnorm(ov[,"IL1A"], main="IL1A gene ovarian cancer")
qqline(ov[,"IL1A"], main="IL1A gene ovarian cancer")
mtext("b)", side=3, adj=-0.23, cex=1.5, line=1.5)


qqnorm(brca[,"CCNE1"], main="CCNE1 gene breast cancer")
qqline(brca[,"CCNE1"], main="CCNE1 gene breast cancer")
mtext("c)", side=3, adj=-0.23, cex=1.5, line=1.5)

qqnorm(brca[,"CDK4"], main="CDK4 gene breast cancer")
qqline(brca[,"CDK4"], main="CDK4 gene breast cancer")
mtext("d)", side=3, adj=-0.23, cex=1.5, line=1.5)
