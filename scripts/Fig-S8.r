#############################################################################
#
# Title: script model-based simulation
#
# Authors: Harry Gray, Gwenael Leday and Catalina Vallejos
#
# Date: Sep 20, 2018
#
#############################################################################




#############################################################################
#############################################################################
### PRELIMINARIES

## Clean environment
rm(list=ls());gc()

## Working directory - set to base github folder
base.dir <-
setwd(base.dir)

## Required libraries and source code
library(TAS) # library(devtools);install_github("HGray384/TAS")
library(MCMCpack)
library(MASS)
library(Matrix)
library(corpcor)
library(ShrinkCovMat)
library(ggplot2)
library(data.table)
library(cgdsr)
library(org.Hs.eg.db)
library(BDgraph)
library(corrplot)
source('./scripts/scm_un.R')


## Data dimension values

# Number of features
p <- 100 
# Number of samples
n <- c(25, 50, 75)

## Auxiliary functions needed

# To generate matrix that stores constant correlation values 
# (diagonal elements set to zero)
corr.mat <- function(p = 10, rho = 0.9){
	mat <- diag(p) + rho
	diag(mat) <- 0
	return(mat)
}
# To generate matrix that stores decaying correlation values 
# (diagonal elements set to zero)
autocorr.mat <- function(p = 10, rho = 0.9){
	mat <- diag(p)
	mat <- rho^abs(row(mat)-col(mat))
	diag(mat) <- 0
	return(mat)
}

## Generating a set of covariance matrices
## One for each simulation scenario

# set seed
set.seed(25)

# param
mycor <- 0.3
myautocor <- -0.7
trueVars <- runif(p, 1, 5)

# Generate true covariance matrices
## Common variance + zero correlation
scenar1 <- 5*diag(p) + corr.mat(p = p, rho = 0)
## Unit variance + constant correlation
scenar2 <- diag(p) + corr.mat(p = p, rho = mycor)
## Unequal variance + decaying correlation
scenar3 <- diag(trueVars) + autocorr.mat(p = p, rho = myautocor)
## Unit variance + block correlation structure 
a <- as.matrix(bdiag(corr.mat(p = p/2, rho = 0.3), 
                     corr.mat(p = p/2, rho = 0.3)) + diag(p))
set.seed(2018)
scenar4 <- cov2cor(solve(BDgraph:::rwish(n=1, p=p, b=200, D = a)))

allScenar <- list(scenar1, scenar2, scenar3, scenar4)

# Check condition number of generated matrices
unlist(lapply(allScenar, kappa, exact=TRUE))
# Check determinant of generated matrices
unlist(lapply(allScenar, det, exact=TRUE))







#############################################################################
#############################################################################
### RUN SIMULATIONS

# Function to compute frobenius norm between every pairs of matrices in list
pairsFrob <- function(listMat){
	mat <- matrix(0, length(listMat), length(listMat))
	for(ii in 1:(length(listMat)-1)){
		for(jj in ii:length(listMat)){
			mat[ii,jj] <- Matrix::norm(listMat[[ii]]-listMat[[jj]], type = "F")^2
		}
	}
	
	mat <- mat + t(mat)

	return(mat)
}


# Function that return matrix of Frob norm 
simu <- function(scenario)
{
  print("1")
  # generate data
  X <- lapply(n, function(x){mvrnorm(n=x, mu=rep(0, p), Sigma=scenario)})
  X <- lapply(X, function(x){t(scale(x, scale=FALSE, center = TRUE))})

  # perform shrinkage and return frobenius loss
  listOfListTargets <- lapply(X, 
                  function(x){
                    varcorcodes <- which(matrix(TRUE, 3, 3), arr.ind=TRUE)
							listTargets <- sapply(1:nrow(varcorcodes), function(y){getTarget(x, var=varcorcodes[y, 1], cor=varcorcodes[y, 2])}, simplify=FALSE)
                    return(pairsFrob(c(listTargets, list(scenario))))
                  }
  )
  list(listOfListTargets)
}


nRep <- 100

# Scenario 1
set.seed(1)
res1 <- replicate(nRep, simu(scenario=scenar1))
mylabs <- c(paste0("T", 1:9), ":Sigma[1]")
list25 <- lapply(res1, function(xx){xx[[1]]})
avgFrob25 <- Reduce("+", list25)/length(list25)
colnames(avgFrob25) <- rownames(avgFrob25) <- mylabs
list50 <- lapply(res1, function(xx){xx[[2]]})
avgFrob50 <- Reduce("+", list50)/length(list50)
colnames(avgFrob50) <- rownames(avgFrob50) <- mylabs
list75 <- lapply(res1, function(xx){xx[[3]]})	
avgFrob75 <- Reduce("+", list75)/length(list75)
colnames(avgFrob75) <- rownames(avgFrob75) <- mylabs

corrplot(round(avgFrob25), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob50), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob75), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

# Scenario 2
set.seed(2)
res2 <- replicate(nRep, simu(scenario=scenar2))
mylabs <- c(paste0("T", 1:9), ":Sigma[2]")
list25 <- lapply(res2, function(xx){xx[[1]]})
avgFrob25 <- Reduce("+", list25)/length(list25)
colnames(avgFrob25) <- rownames(avgFrob25) <- mylabs
list50 <- lapply(res2, function(xx){xx[[2]]})
avgFrob50 <- Reduce("+", list50)/length(list50)
colnames(avgFrob50) <- rownames(avgFrob50) <- mylabs
list75 <- lapply(res2, function(xx){xx[[3]]})
avgFrob75 <- Reduce("+", list75)/length(list75)
colnames(avgFrob75) <- rownames(avgFrob75) <- mylabs

corrplot(round(avgFrob25), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob50), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob75), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")



# Scenario 3
set.seed(3)
res3 <- replicate(nRep, simu(scenario=scenar3))
mylabs <- c(paste0("T", 1:9), ":Sigma[3]")
list25 <- lapply(res3, function(xx){xx[[1]]})
avgFrob25 <- Reduce("+", list25)/length(list25)
colnames(avgFrob25) <- rownames(avgFrob25) <- mylabs
list50 <- lapply(res3, function(xx){xx[[2]]})
avgFrob50 <- Reduce("+", list50)/length(list50)
colnames(avgFrob50) <- rownames(avgFrob50) <- mylabs
list75 <- lapply(res3, function(xx){xx[[3]]})
avgFrob75 <- Reduce("+", list75)/length(list75)
colnames(avgFrob75) <- rownames(avgFrob75) <- mylabs

corrplot(round(avgFrob25), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob50), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob75), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")



# Scenario 4
set.seed(4)
res4 <- replicate(nRep, simu(scenario=scenar4))
mylabs <- c(paste0("T", 1:9), ":Sigma[4]")
list25 <- lapply(res4, function(xx){xx[[1]]})
avgFrob25 <- Reduce("+", list25)/length(list25)
colnames(avgFrob25) <- rownames(avgFrob25) <- mylabs
list50 <- lapply(res4, function(xx){xx[[2]]})
avgFrob50 <- Reduce("+", list50)/length(list50)
colnames(avgFrob50) <- rownames(avgFrob50) <- mylabs
list75 <- lapply(res4, function(xx){xx[[3]]})
avgFrob75 <- Reduce("+", list75)/length(list75)
colnames(avgFrob75) <- rownames(avgFrob75) <- mylabs


corrplot(round(avgFrob25), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob50), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")

corrplot(round(avgFrob75), is.corr = FALSE, mar = c(0.1, 0.1, 0.1, 1.1), method="color", addgrid.col = "black", tl.col='black', tl.cex=c(rep(1,9),1.5), col=colorRampPalette(c("blue", "white", "black"))(200), type="upper")




