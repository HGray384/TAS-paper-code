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

## Working directory - set to base directory of github folder
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

# To generate a pxp matrix that stores constant correlation values rho
# (diagonal elements set to zero)
corr.mat <- function(p = 10, rho = 0.9){
	mat <- diag(p) + rho
	diag(mat) <- 0
	return(mat)
}
# To generate a pxp matrix that stores decaying correlation values rho^|i-j|
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

# Function that generates data, compute estimates and return frobenius losses
simu <- function(scenario)
{
  print("1")
  # generate data
  X <- lapply(n, function(x){mvrnorm(n=x, mu=rep(0, p), Sigma=scenario)})
  X <- lapply(X, function(x){t(scale(x, scale=FALSE, center = TRUE))})
  
  # set shrinkage parameter
  alpha <- seq(0.01, 0.99, 0.01)
  
  # perform shrinkage and return frobenius loss
  frobs <- lapply(X, 
                  function(x){
                    # Single target linear shrinkage estimates 
                    # Coding:
                    # var = 1 (unit), 2 (constant), 3 (unequal)
                    # cor = 1 (r_ij = 0), 2 (r_ij = r), 3 (r_ij = r^abs(i-j))
                    gc11 <- gcShrink(x, alpha = alpha, var=1, cor=1, plots = FALSE, weighted = TRUE)$sigmahat
                    gc12 <- gcShrink(x, alpha = alpha, var=2, cor=1, plots = FALSE, weighted = TRUE)$sigmahat
                    gc13 <- gcShrink(x, alpha = alpha, var=3, cor=1, plots = FALSE, weighted = TRUE)$sigmahat
                    gc21 <- gcShrink(x, alpha = alpha, var=1, cor=2, plots = FALSE, weighted = TRUE)$sigmahat
                    gc22 <- gcShrink(x, alpha = alpha, var=2, cor=2, plots = FALSE, weighted = TRUE)$sigmahat
                    gc23 <- gcShrink(x, alpha = alpha, var=3, cor=2, plots = FALSE, weighted = TRUE)$sigmahat
                    gc31 <- gcShrink(x, alpha = alpha, var=1, cor=3, plots = FALSE, weighted = TRUE)$sigmahat
                    gc32 <- gcShrink(x, alpha = alpha, var=2, cor=3, plots = FALSE, weighted = TRUE)$sigmahat
                    gc33 <- gcShrink(x, alpha = alpha, var=3, cor=3, plots = FALSE, weighted = TRUE)$sigmahat
                    # Target average estimate
                    AUX <- taShrink(x, targets = "default", without = 0, alpha = alpha, plots = FALSE)
                    ta <- AUX$sigmahat
                    ta.weights <- apply(AUX$weights, 1, function(x) sum(x * alpha))
                    # Touloumis estimates
                    AT1 <- shrinkcovmat.identity(x, centered = TRUE)$Sigmahat
                    AT2 <- shrinkcovmat.equal(x, centered = TRUE)$Sigmahat
                    AT3 <- scm.un(x, centered = TRUE)$Sigmahat
                    # corpcor estimate
                    CPC <- cov.shrink(t(x), verbose = FALSE)
                    # Maximum likelihood estimate
                    S <- tcrossprod(x)/ncol(x);
                    ests <- list(gc11, gc12, gc13, gc21, gc22, gc23, gc31, gc32, gc33, ta, AT1, AT2, AT3, CPC, S)
                    # Calculating the Frobenious loss
                    frobloss <- lapply(ests, function(y){Matrix::norm(scenario-y, type = "F")^2})
                    
                    # Returns Frobenious loss and posterior weights (TAS only)
                    list("frobloss" = frobloss, "ta.weights" = ta.weights)
                  }
  )
  
  
  list("frobloss" = simplify2array(lapply(frobs, function(x) x$frobloss)),
       "ta.weights" = simplify2array(lapply(frobs, function(x) x$ta.weights)))
}

# repeat this function 100 times
nRep <- 100
set.seed(1)
res1 <- replicate(nRep, simu(scenario=scenar1))
set.seed(2)
res2 <- replicate(nRep, simu(scenario=scenar2))
set.seed(3)
res3 <- replicate(nRep, simu(scenario=scenar3))
set.seed(4)
res4 <- replicate(nRep, simu(scenario=scenar4))

# res objects are 100-dimensional arrays, where each dimension contains a list with:
# 1) Frobenious loss (w.r.t. S) for 14 estimators, 3 sample sizes
# 2) Target-specific weights (9 targets), 3 sample sizes (TAS only)

#############################################################################
#############################################################################
### PLOT RESULTS

# Labels
labelsMethods <- c(paste0("ST", 1:9), "TAS", paste0("AT", 1:3), "CPC", "S")

# Check missing values
unlist(lapply(list(res1, res2, res3, res4), function(yy){any(is.na(yy))}))

# Function that computes PRIAL
get.prial <- function(res, nindx=1){
  tp <- res[,nindx,]
  tp <- sapply(1:nrow(tp), function(x){unlist(tp[x,])}, simplify = FALSE)
  term1 <- rep(unlist(lapply(tp, sum))[length(labelsMethods)], length(labelsMethods)-1)
  term2 <- unlist(lapply(tp, sum, na.rm=FALSE))[-length(labelsMethods)]
  prial <- ((term1-term2)/term1)*100
  prial
}

res <- list(simplify2array(res1[1,]), simplify2array(res2[1,]),
            simplify2array(res3[1,]), simplify2array(res4[1,]))
            
# Get all prials for n=25
allPrials25 <- sapply(res, get.prial, nindx=1)

# Get all prials for n=50
allPrials50 <- sapply(res, get.prial, nindx=2)

# Get all prials for n=75
allPrials75 <- sapply(res, get.prial, nindx=3)

# Labels
A <- c("scenar1", "scenar2", "scenar3", "scenar4")
colnames(allPrials25) <- colnames(allPrials50) <- colnames(allPrials75) <- A
B <- labelsMethods[-length(labelsMethods)]
rownames(allPrials25) <- rownames(allPrials50) <- rownames(allPrials75) <- B

# Function used to personalise barplot 
# labs - the estimation method (e.g. TAS)
# vals - the prial results achieved (e.g. scenario 1 n=25)
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

# Plot results - PRIAL
for(ii in 1:4){
	assign(paste0("barplotScenar",ii),mybarplot(labs=labelsMethods[-length(labelsMethods)], vals=allPrials25[,ii]))
	plot(get(paste0("barplotScenar",ii)))
}
for(ii in 1:4){
  assign(paste0("barplotScenar",ii),mybarplot(labs=labelsMethods[-length(labelsMethods)], vals=allPrials50[,ii]))
  plot(get(paste0("barplotScenar",ii)))
}
for(ii in 1:4){
  assign(paste0("barplotScenar",ii),mybarplot(labs=labelsMethods[-length(labelsMethods)], vals=allPrials75[,ii]))
  plot(get(paste0("barplotScenar",ii)))
}

# Plot results - posterior weights (TAS only)
res <- list(simplify2array(res1[2,]), simplify2array(res2[2,]),
            simplify2array(res3[2,]), simplify2array(res4[2,]))

# Function used to personalise barplot 
myboxplot <- function(vals, size){
  tp <- data.table(t(vals[,size,]))
  names(tp) <- paste0("T", seq_len(9))
  tp$S <- 1 - rowSums(tp)
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

# Plot results - Weights
for(ii in 1:4){
  assign(paste0("weightsScenar",ii),myboxplot(vals = res[[ii]], size = 1))
  plot(get(paste0("weightsScenar",ii)))
}
for(ii in 1:4){
  assign(paste0("weightsScenar",ii),myboxplot(vals = res[[ii]], size = 2))
  plot(get(paste0("weightsScenar",ii)))
}
for(ii in 1:4){
  assign(paste0("weightsScenar",ii),myboxplot(vals = res[[ii]], size = 3))
  plot(get(paste0("weightsScenar",ii)))
}
