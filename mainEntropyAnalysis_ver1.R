
setwd("~/Documents/Anchor/Initiatives/Keyfitz' Error/Keyfitz_Aug22/Keyfitz_August2022")

# Clear Memory --------------------------------------------------------------------
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

# Libraries --------------------------------------------------------------------
library(scales)
library(Rcompadre)
library(Rage)
library(tidyverse)
library(popbio)
library(parallel)
# library(dplyr) 

# Functions --------------------------------------------------------------------
## ** Isolating Age-from-Stage - Bookkeeping functions for COMADRE/COMPADRE ----------

# Identify age-inclusive matrices that have a leslie format (no retrogression + stasis)
leslieElements <- function(A_matrix){
  holdMatrix <- A_matrix
  
  zeroVector <- rep(0, dim(A_matrix)[[1]])
  holdMatrix[1,] <- zeroVector
  
  for(i in 1:dim(A_matrix)[[1]]){
    if(i < dim(A_matrix)[[1]]){
      holdMatrix[i+1, i] <- zeroVector[i]
    }else if (i == dim(A_matrix)[[1]]){
      holdMatrix[i, i] <- zeroVector[i]
    }
  }
  nonZeroElements <- sum(holdMatrix)
  return(nonZeroElements==0)
}

# Evaluate age-inclusie matrices 
leslieApply <- function(database){
  leslieElemList <- list()
  for(i in 1:length(database$mat)){
    leslieElemList[[i]] <- leslieElements(database$mat[[i]]@matA)
  }
  return(which(unlist(leslieElemList)))
}


## ** Entropy-related functions ----------

# Formating adjustment for Rage - maps columns onto rows
prepMatrix <- function(matrix){
  holdmat <- try(as.matrix(matrix))
  try(rownames(holdmat) <- colnames(holdmat))
  try(colnames(holdmat) <- rownames(holdmat))
  return(holdmat)
}

# New Keyfitz Metric
NewKeyfitz <- function(U_matrix){
  maxage <- dim(U_matrix)[1]
  e0 <- rep(0,maxage)   #e0=zeros(maxage,1);
  e0[1] <- 1            #e0(1)=1;
  l <- rep(0,2*maxage) # double the duration
  temp <- rep(1,maxage)-colSums(U_matrix)  # temp=ones(1,maxage)-sum(U{i});
  M <- diag(temp)                         # M=diag(temp);
  N <- solve(diag(maxage)-U_matrix)        # N=(eye(maxage)-U{i})\eye(maxage); 
  eta1 <- rep(1,maxage)%*%N               # eta1=ones(1,maxage)*N
  B <- M%*%N                              # B=M*N;
  etadagger <- eta1%*%B                   # etadagger=eta1*B;
  eta0 <- eta1%*%e0                       # eta0=eta1*e0;
  H_new <- etadagger[1]/eta0              # H_new(i)=etadagger(1)/eta0
  return(H_new)
}


# Original Keyfitz' entropy metric from Finch's Hypothesis paper
OriginalKeyfitz <- function(U_matrix){
  try(lx <- Rage::mpm_to_lx(prepMatrix(U_matrix)))
  return(try(as.numeric(-t(lx[is.na(lx)==FALSE])%*%log(lx[is.na(lx)==FALSE])/sum(lx[is.na(lx)==FALSE]))))
}


# Life Expectancy from Rage
LifeExpectancy <- function(U_matrix){
  return(
    Rage::life_expect_mean(
      prepMatrix(U_matrix)
    )
  )
}


## ** Age-to-Stage Functions



## ** Age-to-Stage functions ---------------------------------------------------

buildMatrix <- function(survCurv){
  
  baseZeros <- data.frame(matrix(0,
                                 nrow=length(survCurv), 
                                 ncol=length(survCurv)))
  
  for(i in 1:(length(survCurv)-1)){
    baseZeros[i+1,i] <- survCurv[i]
  }
  for(i in length(survCurv)){
    baseZeros[i, i] <- survCurv[i]
  }
  baseZeros <- as.matrix(baseZeros)
  colnames(baseZeros) <- rownames(baseZeros)
  return(baseZeros)
}

applyAlpha <- function(lifeTable, alpha){
  return(1:seq(1:length(lifeTable[,2]))[lifeTable[,2]<alpha][1])
}

matToLifeTable <- function(database, index){
  return(Rage::mpm_to_table(prepMatrix(database@data$mat[[index]]@matU),
                            prepMatrix(database@data$mat[[index]]@matF)))
}

StageToAge <- function(database, index, alpha){
  px_N_alpha <- matToLifeTable(database, index)$px
  mx_N_alpha <- matToLifeTable(database, index)$mx
  newMatrix <- buildMatrix(px_N_alpha[applyAlpha(matToLifeTable(database, index), alpha)])
  newMatrix[1,] <- mx_N_alpha[applyAlpha(matToLifeTable(database, index), alpha)]
  return(newMatrix)
}


# Load Data ------------------------------------------------------------------------

# Load Data File from Finch List -- Similar characteristics to the spp list filtered by strong 
# criteria on duration, matrix dimesions, etc.
finchDataComadre <- read.csv("FinchComadreSpp.csv") # Import COMPADRE species list
finchDataCompadre <- read.csv("FinchCompadreSpp.csv") # Import COMPADRE species list

# length(finchDataComadre$Species)
# length(finchDataCompadre$Species)

# Load COMPADRE & COMADRE
# compadre <- cdb_fetch("compadre")
# comadre <- cdb_fetch("comadre")

comadre <- readRDS("comadre_Aug22.rds")
compadre <- readRDS("compadre_Aug22.rds")

# Isolate Age-From-Stage ------------------------------------------------------------------------

# Extract prospective age-only matrices (includes at least one age-based st/age)
comadre_age <- comadre[which(comadre$MatrixCriteriaAge=="Yes")] # Age-based matrices
compadre_age <- compadre[which(compadre$MatrixCriteriaAge=="Yes")] # Age-based matrices

comadre_leslie <- leslieApply(comadre_age)
compadre_leslie <- leslieApply(compadre_age)

comadre_leslie <- comadre[comadre_leslie]
compadre_leslie <- compadre[compadre_leslie]

# Isolate Matrices that Maximise Dimensional ---------------------------------------------

# Select matrix record that maximises  matrix dimension within each species
comadreMaxDim <- comadre_leslie %>%
  group_by(SpeciesAccepted) %>%
  slice(which.max(MatrixDimension))

compadreMaxDim <- compadre_leslie %>%
  group_by(SpeciesAccepted) %>%
  slice(which.max(MatrixDimension))

# Isolate Finch Papers ----------------------------------------------------------

comadreFinchIndex <- match(finchDataComadre$Species, comadreMaxDim$SpeciesAccepted)
compadreFinchIndex <- match(finchDataCompadre$Species, compadreMaxDim$SpeciesAccepted) # b/c age almost no overlap

# intersect(compadreMaxDim$SpeciesAccepted, finchDataCompadre$Species) # note only 2 overlap b/c age

comadreFinchIndex_noNa <- comadreFinchIndex[!is.na(comadreMaxDim$SpeciesAccepted[comadreFinchIndex])]
compadreFinchIndex_noNa <- compadreFinchIndex[!is.na(compadreMaxDim$SpeciesAccepted[compadreFinchIndex])]

comadre_finch <- comadreMaxDim[comadreFinchIndex_noNa]
compadre_finch <- compadreMaxDim[compadreFinchIndex_noNa]

length(comadreMaxDim$mat)
length(comadre_finch$mat)


# Calculating Entropy  ---------------------------------------------

calculateEntropy <- function(database){
  new.entropy <- list()
  original.entropy <- list()
  life.expectancy <- list()
  
  for(i in 1:length(database$MatrixID)){
    new.entropy[[match(database$MatrixID[[i]], database$MatrixID)]] <- NA
    original.entropy[[match(database$MatrixID[[i]], database$MatrixID)]] <- NA
    life.expectancy[[match(database$MatrixID[[i]], database$MatrixID)]] <- NA
    
    try(new.entropy[[match(database$MatrixID[[i]], database$MatrixID)]] <- 
      NewKeyfitz(database$mat[[i]]@matU))
    
    try(original.entropy[[match(database$MatrixID[[i]], database$MatrixID)]] <- 
      OriginalKeyfitz(database$mat[[i]]@matU))
    
    try(life.expectancy[[match(database$MatrixID[[i]], database$MatrixID)]] <- 
      LifeExpectancy(database$mat[[i]]@matU))
  }
  
  new.entropy <- as.numeric(unlist(new.entropy))
  original.entropy <- as.numeric(unlist(original.entropy))
  life.expectancy <- as.numeric(unlist(life.expectancy))
  
  outputs <- list(new.entropy, original.entropy, life.expectancy)
  names(outputs) <- c("NewEntropy", "OriginalEntropy", "lifeExpOut")
  return(outputs)
}


comadreEntropy <- calculateEntropy(comadre_finch)
compadreEntropy <- calculateEntropy(compadre_finch)


