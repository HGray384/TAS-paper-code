#############################################################################
#
# Title: script TCPA application
#
# Authors: Harry Gray, Gwenael Leday and Catalina Vallejos
#
# Date: Oct 31, 2018
#
#############################################################################





#############################################################################
#############################################################################
### PRELIMINARIES

# clean environment
rm(list=ls());gc()

# model-based simulations
if(!require("TAS")){
  if(!require("devtools")){
    install.packages("devtools")
  }
  library(devtools)
  install_github("HGray384/TAS")
}
library(TAS)
if(!require("abind")){
  install.packages("abind")
}
library(abind)
if(!require("corrplot")){
  install.packages("corrplot")
}
library(corrplot)

# Import pancancer proteomic data
# the data can be found at 'http://tcpaportal.org/tcpa/download.html'
# using the /Downloads/4.2/TCGA/Pan-Can 32/Level 4/TCGA-PANCAN32-L4.zip
# folder of the drop down menu
# The data used here was downloaded on 21/Aug/2018

# set the directory for the saved TCPA data on your local machine
datadir <- 
PAN <- read.csv(datadir, stringsAsFactors = FALSE)

# check dimension - last version:
dim(PAN)
# [1] 7694  261

# Remove annotations
colnames(PAN)[1:4]
pandat <- data.matrix(PAN[,-c(1:4)])

# add sample names
rownames(pandat) <- PAN$Sample_ID

# Remove protein containing missing values in pancancer data set
pandat <- pandat[, !apply(pandat, 2, function(x){any(is.na(x))})]

# Function to return regularized covariance estimate
varcorcodes <- which(matrix(TRUE, 3, 3), arr.ind=TRUE) # all combinations
# of c(1, 2, 3) with c(1, 2, 3) - used for variance and correlation
# parameter codes in getTarget()
alpha <- seq(0.01, 0.99, 0.01) # shrinkage grid
regS <- function(y, varcorcodes){
  # returns the regularised sample covariance matrix of matrix input y
  # by running TAS with default target set
  # the target set is indexed by varcorcodes according to variance 
  # and correlation structure (?taShrink)
	thelistTargets <- sapply(1:nrow(varcorcodes), 
	                         function(x) {getTarget(t(y), var=varcorcodes[x, 1], cor=varcorcodes[x, 2])}, 
	                         simplify=FALSE)
	thearrayTargets <- simplify2array(thelistTargets)
	mytas <- taShrink(t(y), targets = thearrayTargets, 
	                  without = 0, alpha = alpha, plots = FALSE)
	mytas$sigmahat
}

# Get single data-based target from pooled data
pandat <- scale(pandat, scale=FALSE, center = TRUE) # scaling necessary
# since getTarget is called directly (and does not automatically center)

# Create list of pan cancer data sets
cancerTypes <- sort(unique(PAN$Cancer_Type))
pandatList <- sapply(1:length(cancerTypes), 
                     function(x) { pandat[PAN$Cancer_Type==cancerTypes[x],] }, 
                     simplify=FALSE)
names(pandatList) <- cancerTypes
# Sample size per cancer
simplify2array(lapply(pandatList, nrow))

# Get 31 data-based targets from each cancer type
dataBasedTargetsList <- lapply(pandatList, regS, varcorcodes=varcorcodes)
dataBasedTargetsArray <- simplify2array(dataBasedTargetsList)
# Check dimensions are ok
simplify2array(lapply(dataBasedTargetsList, dim))

#############################################################################
#############################################################################
### Apply TAS to three selected cancers  (Fig 5)
### 9 default targets + cancers derived from other (single) cancers + pancan

FigFive <- function(pandatList, Cancer, varcorcodes, 
                    dataBasedTargetsList, dataBasedTargetsArray) {
  ## generates Figure 5
  ## see above variables for the argument descriptions
  
  # Get 9 'default' targets
  dat <- scale(pandatList[[Cancer]], center = TRUE, scale = FALSE)
  listTargets <- sapply(1:nrow(varcorcodes), 
                        function(x){getTarget(t(dat), var=varcorcodes[x, 1], 
                                              cor=varcorcodes[x, 2])}, simplify=FALSE)  
  defaultTargets <- simplify2array(listTargets)
  
  # Get regularized pooled sample estimate (other cancer types)
  tpdat <- scale(pandat[PAN$Cancer_Type != Cancer,], scale=FALSE, center = TRUE)
  dataBasedTargetPooled <- regS(tpdat, varcorcodes = varcorcodes)
  
  # Combine targets
  arrayTargets <- abind(defaultTargets,
                        dataBasedTargetsArray[,,-which(names(dataBasedTargetsList)== Cancer)],
                        dataBasedTargetPooled)
  
  # Apply TAS
  tasCancer <- taShrink(t(dat), targets = arrayTargets, 
                        without = 0, alpha = alpha, plots = FALSE)
  
  # Obtain posterior weights
  weightsCancer <- targetWeights(tasCancer)
  labelsCancer <- format(c(paste0("T", 1:9), cancerTypes[cancerTypes != Cancer], 
                         "PANCAN", "S"), justify="right")
  names(weightsCancer) <- labelsCancer
  
  # Plot weights
  par(mar=c(4.5, 4, 0.5, 0.5) + 0.1)
  xx <- barplot(weightsCancer, ylim=c(0, 1.1), ylab="Posterior weight", 
                xlab="", main="", col=rgb(48, 48, 48, maxColorValue = 255), 
                axisnames=FALSE)
  axis(1, at=xx[,1], labels = rep("", length(xx[,1])))
  text(xx[,1], par("usr")[3], srt = 60, adj= c(1.3,1.3), xpd = TRUE,
       labels = names(weightsCancer), cex=1)
  box()

}

for(Cancer in c("CHOL", "LIHC", "READ")) {
  cat(Cancer, "\n")
  FigFive(pandatList, Cancer, varcorcodes, 
          dataBasedTargetsList, dataBasedTargetsArray)  
}

#############################################################################
#############################################################################
### Apply TAS to all cancers using 31 (single-)cancer-derived targets only (Fig 6)

matrixWeights <- matrix(NA, length(cancerTypes), length(cancerTypes))
rownames(matrixWeights) <- colnames(matrixWeights) <- cancerTypes

for(j in 1:length(cancerTypes)){
  # print iteration
  cat("Cancer number", j, "of", length(cancerTypes), "\n", sep = " ")

	# Get data of cancer j
	pandatj <- pandat[PAN$Cancer_Type == cancerTypes[j],]

	# Use only data-derived targets
	arrayTargetsPAN <- dataBasedTargetsArray[,,-j]

	# Apply TAS
	tasPAN <- taShrink(t(pandatj), targets = arrayTargetsPAN, 
	                   without = 0, alpha = alpha, plots = FALSE)

	# Obtain posterior weights
	weightPAN <- targetWeights(tasPAN)
	matrixWeights[j, -j] <- weightPAN[-length(cancerTypes)]
	matrixWeights[j, j] <- weightPAN[length(cancerTypes)]
}

# create corrplot
corrplot(matrixWeights, is.corr = FALSE, mar = c(0.1, 2.1, 1.1, 2.1),
         tl.col='black', method="color", cl.lim=c(0,1), addgrid.col = "black",
         col=colorRampPalette(c("blue", "white", "black"))(200),
         order = 'FPC')
mtext("dataset", 2, adj=0.44, line=2, cex=1.5)
mtext("shrinkage targets", 3, adj=0.5, line=2, cex=1.5)