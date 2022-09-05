

# Records Not Captured in the Age-Based Calculations ------------------------------

comadreDiff <- setdiff(finchDataComadre$Species,
                       comadre_finch$SpeciesAccepted)

compadreDiff <- setdiff(finchDataCompadre$Species,
                        compadre_finch$SpeciesAccepted)

nonAgeComadre <- match(comadreDiff,
                       comadre$SpeciesAccepted)

nonAgeCompadre <- match(compadreDiff,
                        compadre$SpeciesAccepted)


StageFinchComadre <- comadre[nonAgeComadre[!is.na(nonAgeComadre)]]
StageFinchCompadre <- compadre[nonAgeCompadre[!is.na(nonAgeCompadre)]]



# Stage Status Review -------------------------------------------------------

stageStat <- function(database, i){
  returnObject <- NA
  try(lambdaInitial <- lambda(database@data$mat[[i]]@matA))
  try(genTimeInitial <- generation.time(database@data$mat[[i]]@matA))
  try(netRepRateInitial <- net.reproductive.rate(database@data$mat[[i]]@matA))

  try(stageModel <- StageToAge(database, i, 0.01))
  try(lambdaStage <- lambda(stageModel))
  try(genTimeStage <- generation.time(stageModel))
  try(netRepRateStage <- net.reproductive.rate(stageModel))

  try(lambdaStatus <- abs(lambdaInitial-lambdaStage)<0.05)
  try(genTimeStatus <- abs(genTimeInitial-genTimeStage)<0.05)
  try(netRepRateStatus <- abs(netRepRateInitial-netRepRateStage)<0.05)

  try(if(sum(lambdaStatus, genTimeStatus, netRepRateStatus)==3){
    returnObject <- TRUE}else if(!sum(lambdaStatus, genTimeStatus, netRepRateStatus)==3){
      returnObject <- FALSE
    })

  if(is.na(returnObject)){returnObject <- FALSE}

  return(returnObject)
}

stageStatComadre <- mcmapply(stageStat,
                             i = 1:length(StageFinchComadre@data$mat),
                             MoreArgs = list(database = StageFinchComadre),
                             mc.silent = TRUE,
                             mc.cores = 6)

stageStatCompadre <- mcmapply(stageStat,
                             i = 1:length(StageFinchCompadre@data$mat),
                             MoreArgs = list(database = StageFinchCompadre),
                             mc.silent = TRUE,
                             mc.cores = 6)

stageGoodComadre <- StageFinchComadre[stageStatComadre]
stageGoodCompadre <- StageFinchCompadre[stageStatCompadre]

# Stage Status Review -------------------------------------------------------


matToLifeTable
StageToAge(stageGoodComadre, i, 0.1)


newAgeComadre <- mcmapply(StageToAge,
                          i = 1:length(stageGoodComadre@data$mat),
                          MoreArgs = list(database = stageGoodComadre,
                                          alpha = 0.1),
                          mc.silent = TRUE,
                          mc.cores = 6)

newAgeCompadre <- mcmapply(StageToAge,
                          i = 1:length(stageGoodCompadre@data$mat),
                          MoreArgs = list(database = stageGoodCompadre,
                                          alpha = 0.1),
                          mc.silent = TRUE,
                          mc.cores = 6)


# Stage of Keyfitz -------------------------------------------------------------

dim(newAgeCompadre[[1]])[1]

newAgeCompadre[[1]][1,]

calculateEntropyFromMat <- function(database){
  new.entropy <- list()
  original.entropy <- list()
  life.expectancy <- list()
  
  for(i in 1:length(database)){
    holdMatrix <- database[[i]]
    holdMatrix[1,] <- rep(0, dim(holdMatrix)[1])
    
    try(new.entropy[[i]] <- NewKeyfitz(holdMatrix))
    try(original.entropy[[i]] <- OriginalKeyfitz(holdMatrix))
    try(life.expectancy[[i]] <- LifeExpectancy(holdMatrix))
  }
  
  new.entropy <- as.numeric(unlist(new.entropy))
  original.entropy <- as.numeric(unlist(original.entropy))
  life.expectancy <- as.numeric(unlist(life.expectancy))
  
  outputs <- list(new.entropy, original.entropy, life.expectancy)
  names(outputs) <- c("NewEntropy", "OriginalEntropy", "lifeExpOut")
  return(outputs)
}


stageEntropyComadre <- calculateEntropyFromMat(newAgeComadre)
stageEntropyCompadre <- calculateEntropyFromMat(newAgeCompadre)



stageEntropyComadre <- calculateEntropyFromMat(newAgeComadre)








































