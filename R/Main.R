setwd("/Users/surgery/Project/HOME/github/MBSIT/R")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")
options(stringsAsFactors = F)
source("/Users/surgery/Project/HOME/myScript/zxli_lib.R")

# cellCountPerGene <- rowSums(counts(tmp_group)>=5)
# geneCountPerCell <- colSums(counts(tmp_group)>=5)
# totalReadCount <- colSums(counts(tmp_group))
# cvPerGene <- apply(counts(tmp_group), 1, sd)/apply(counts(tmp_group), 1, mean)
# ### cvPerCell <- apply(counts(tmp_group), 2, sd)/apply(counts(tmp_group), 2, mean)

epsDetection <- function(pcMatrix=mypcMatrix, start=1, end=100, fold=5, primary=T, minPts=10, percent=0.05, plot=T){
  # DBSCAN
  library(dbscan)
  if (primary==T){ by = 1 } else {by=(end-start)/100}
  epsLevels <- (seq(start, end, by=by)-1)*fold
  allnoiseCount <- c()
  alleps <- c()
  for ( eps in epsLevels){
    #eps <- 5*(i-1)
    alleps <- c(alleps, eps)
    # minPts <- 10
    res <- dbscan(pcMatrix, eps = eps, minPts = minPts) 
    noiseCount <- table(res$cluster)["0"]
    allnoiseCount <- c(allnoiseCount, noiseCount)
  }
  names(allnoiseCount) <- alleps
  allnoiseCount <- allnoiseCount[!is.na(allnoiseCount)]
  allnoiseCount <- allnoiseCount[allnoiseCount!=dim(pcMatrix)[1]]
  if (plot) {plot(allnoiseCount)}
  if (primary==T)
  {return(c(as.integer(names(allnoiseCount)[1]), as.integer(names(allnoiseCount)[length(allnoiseCount)])))}
  else {return(as.double(names(sort(allnoiseCount[allnoiseCount>=percent*dim(pcMatrix)[1]])[1])))}
  # {return(as.double(names(sort(allnoiseCount %% as.integer(percent*dim(pcMatrix)[1]))[1])))}
}

epsDetectionCor <- function(disMatrix=dis_matrix_dist, start=0, end=1, step=10, minPts=5, plot=F){
  # DBSCAN
  library(dbscan)
  epsLevels <- seq(start, end, by=(end-start)/step)
  allClusterCount <- c()
  alleps <- c()
  for ( eps in epsLevels){
    #eps <- 5*(i-1)
    alleps <- c(alleps, eps)
    # minPts <- 10
    res <- dbscan(disMatrix, eps = eps, minPts = minPts) 
    clusterCount <- length(table(res$cluster)) - 1
    print(paste("eps: ", eps, ", gene cluster number: ", clusterCount))
    allClusterCount <- c(allClusterCount, clusterCount)
  }
  names(allClusterCount) <- alleps
  maxStart <- which(allClusterCount==max(allClusterCount))[1]
  maxEnd <- which(allClusterCount==max(allClusterCount))[length(which(allClusterCount==max(allClusterCount)))]
  if (length(which(allClusterCount==max(allClusterCount)))>1) {return(as.double(names(allClusterCount[maxStart])))}
  return(c(as.double(names(allClusterCount[maxStart-1])), as.double(names(allClusterCount[maxStart+1]))))
}

centerDetection <- function(corM=dis_matrix, genelist=genelist){
  if (length(genelist)==1) {
    # print("WARN: only one gene in genelist!")
    return(genelist[1])}
  corMsub <- corM[genelist, genelist]
  corRowSum <- rowSums(corMsub)
  return(names(sort(corRowSum)[1]))
}

singleModuleDetectionCenter <- function(corM=dis_matrix, overlapM=geneOverlap,center=center, topCount=10, corThresd=0.5, overlapThresd=3, centerBlackList=c()){
  genelist <- c(center)
  while (length(genelist) < topCount) {
    centerBlackList=c(centerBlackList, center)
    nearestGeneV <- sort(corM[center,!colnames(corM)%in%genelist])[1]
    nearestGene <- names(nearestGeneV)
    if (nearestGeneV > corThresd) {
      # print("nearestGeneV > corThresd")
      return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
    } else if (overlapM[center, nearestGene] < overlapThresd) {
      # print("overlapM[center, nearestGene] < overlapThresd")
      return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
    } 
    genelist <- c(genelist, nearestGene)
    center <- centerDetection(corM, genelist)
    # print(paste("nearestGene is: ", nearestGene, ", new center is: ", center, sep=" "))
  }
  # genelist <- moduleShift(corM, genelist)
  # pheatmap(logcounts(tmp_group)[genelist,], show_colnames = F, cluster_rows = F)
  return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
}

singleModuleDetectionSum <- function(corM=dis_matrix, center=center, topCount=100, corThresd=0.5, cutThresd=0.6, centerBlackList=c()){
  # consider the densityDf
  genelist <- c(center)
  centerBlackList=c(centerBlackList, center)
  while (length(genelist) < topCount) {
    if (length(genelist)==1){
      nearestGeneV <- sort(corM[center,!colnames(corM)%in%genelist])[1]
      nearestGene <- names(nearestGeneV)
      if (nearestGeneV > corThresd) {return(genelist)} 
    } else if (length(genelist)>1) {
      nearestGeneV <- sort(colSums(corM[genelist,!colnames(corM)%in%genelist]))[1]
      nearestGene <- names(nearestGeneV)
      corDistribution <- corM[genelist,nearestGene] > corThresd
      # cut, avoid merge
      if (sum(corDistribution==T)/length(corDistribution) > cutThresd ) {return(genelist)}
    }
    # centerBlackList=c(centerBlackList, nearestGene)
    genelist <- c(genelist, nearestGene)
    if (nearestGene %in% centerBlackList) {
      print("nearestGene was in centerBlackList")
      return(genelist)
      next
    }
    # center <- centerDetection(corM, genelist)
    # print(paste("nearestGene is: ", nearestGene, ", new center is: ", center, sep=" "))
  }
  # genelist <- moduleShift(corM, genelist)
  # pheatmap(logcounts(tmp_group)[genelist,], show_colnames = F, cluster_rows = F)
  # return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
  return(genelist)
}

moduleShiftCenter <- function(corM=dis_matrix, genelist=genelist, centerBlackList=c()) {
  if (length(genelist)==1) {
    return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
  }
  while (T){
    center <- centerDetection(corM, genelist)
    centerBlackList=c(centerBlackList, center)
    nearestGene <- names(sort(corM[center,!colnames(corM)%in%genelist])[1])
    disSumOutlier <- sort(rowSums(corM[genelist, genelist]), decreasing = T)[1]
    genelist1 <- c(genelist, nearestGene)
    genelist1 <- genelist1[!genelist1%in%names(disSumOutlier)]
    nearestGeneDis <- sum(corM[nearestGene, genelist1])
    if (nearestGeneDis < disSumOutlier) {
      genelist <- genelist1
      # print(paste(nearestGene, "is closer than", names(disSumOutlier), ". replace it already!"))
    }
    else {
      # print("Gene set has local convergence")
      break}
  }
  return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
}

moduleShiftSum <- function(corM=dis_matrix, genelist=genelist, centerBlackList=c()) {
  if (length(genelist)==1) {
    return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
  }
  while (T){
    center <- centerDetection(corM, genelist)
    centerBlackList=c(centerBlackList, center)
    nearestGene <- names(sort(corM[center,!colnames(corM)%in%genelist])[1])
    disSumOutlier <- sort(rowSums(corM[genelist, genelist]), decreasing = T)[1]
    genelist1 <- c(genelist, nearestGene)
    genelist1 <- genelist1[!genelist1%in%names(disSumOutlier)]
    nearestGeneDis <- sum(corM[nearestGene, genelist1])
    if (nearestGeneDis < disSumOutlier) {
      genelist <- genelist1
      # print(paste(nearestGene, "is closer than", names(disSumOutlier), ". replace it already!"))
    }
    else {
      # print("Gene set has local convergence")
      break}
  }
  return(list(genelist=genelist, centerBlackList=unique(centerBlackList)))
}

fullModuleDetectionCenter <- function(corM=dis_matrix, marker_result=marker_result){
  moduleResult <- list()
  count <- 1
  #count2 <- 0
  j <- 0
  #totalCount <- dim(corM)[1]
  centerBlackList <- c()
  totalLength <- length(table(marker_result$cluster))
  for (i in 1:totalLength) {
    #for (center in rownames(corM)) {
    geneCluster <- marker_result[marker_result$cluster==i,]$gene
    center <- centerDetection(corM, geneCluster)
    #count2 <- count2 + 1
    #if (count2%%as.integer(totalCount/10)==0){
    j <- j+1
    print(paste(j, "of", totalLength, "was finished..."))
    # center <- i
    # print(center)
    result1 <- singleModuleDetection(corM=corM, overlapM=geneOverlap, center=center, centerBlackList=centerBlackList, topCount=100)
    genelist <- result1[["genelist"]]
    centerBlackList1 <- result1[["centerBlackList"]]
    result2 <- moduleShift(corM, genelist, centerBlackList=centerBlackList1)
    genelist <- result2[["genelist"]]
    centerBlackList2 <- result2[["centerBlackList"]]
    tempCenter <- centerDetection(corM, genelist)
    if (tempCenter%in%centerBlackList) {
      # print(paste(tempCenter, "is duplicated center!", sep=" "))
      next}
    centerBlackList <- c(centerBlackList, centerBlackList2)
    # print(genelist)
    if (length(genelist) >= 5){
      moduleResult[[count]] <- genelist
      count <- count + 1
    }
  }
  return(moduleResult)
}

# densityDf <- highDensityCenterDetection()
fullModuleDetectionSum <- function(corM=dis_matrix, localCenters=rownames(densityDf), minModuleGene=8, cutThresd=disThresd){
  moduleResult <- list()
  moduleResultbak <- list()
  count <- 1
  count2 <- 1
  #count2 <- 0
  j <- 0
  #totalCount <- dim(corM)[1]
  centerBlackList <- c()
  totalLength <- length(localCenters)
  # for (i in 1:totalLength) {
  for (center in localCenters) {
    if (center%in%centerBlackList) {
      print(paste("center: ", center, "is in centerBlackList!", sep=" "))
      next}
    j <- j+1
    genelist <- singleModuleDetectionSum(corM=corM, center=center, cutThresd=cutThresd, centerBlackList=centerBlackList, topCount=100)
    # genelist has duplicate with exsiting moduleResult
    if (length(intersect(genelist, centerBlackList)) > 0) {
      moduleResultbak[[count2]] <- genelist
      count2 <- count2 + 1
      next
    }
    centerBlackList <- unique(c(centerBlackList, genelist))
    if (length(genelist) >= 5){
      moduleResult[[count]] <- genelist
      count <- count + 1
    }
  }
  # merge moduleResultbak to moduleResult
  for (i in 1:length(moduleResultbak)) {
    for (j in 1:length(moduleResult)){
      if (length(intersect(moduleResult[[j]], moduleResultbak[[i]]))){
        moduleResult[[j]] <- unique(c(moduleResult[[j]], moduleResultbak[[i]]))
      }
    }
  }
  # moduleResultMerge <- list()
  count <- 1
  overlapPair <- list()
  for (i in 1:length(moduleResult)) {
    for (j in 1:length(moduleResult)){
      if (i >= j) {next}
      if ((length(intersect(moduleResult[[j]], moduleResult[[i]]))/(min(length(moduleResult[[j]]), length(moduleResult[[i]]))))>0.5) {
        #moduleResultMerge[[count]] <- unique(c(moduleResult[[j]], moduleResult[[i]]))
        overlapPair[[count]] <- c(i, j)
        count <- count + 1
      } 
    }
  }
  if (length(overlapPair)>0) {
    for (i in 1:length(overlapPair)) {
    start <- overlapPair[[i]][1]
    end <- overlapPair[[i]][2]
    moduleResult[[start]] <- unique(c(moduleResult[[start]], moduleResult[[end]]))
    }
    for (i in 1:length(overlapPair)) {
      end <- overlapPair[[i]][2]
      moduleResult[[end]] <- NULL
    }
  }
  
  moduleResultDf <- data.frame()
  for (i in 1:length(moduleResult)){
    if (length(moduleResult[[i]]) < minModuleGene) {next}
    for (j in moduleResult[[i]]) { 
      moduleResultDf <- rbind(moduleResultDf, c(j, i))
    }
  }
  colnames(moduleResultDf) <- c("gene", "module")
  moduleResultDf <- moduleResultDf[!duplicated(moduleResultDf$gene),]
  rownames(moduleResultDf) <- moduleResultDf$gene

  # moduleResultDf["EBF1",]
  # gene <- a
  # moduleResultDf[gene[gene%in%moduleResultDf$gene],]
  return(moduleResultDf)
}

localCentersDetection <- function(corM=dis_matrix, densityDf=densityDf, topGene=10, minModuleGene=8, cutThresd=disThresd){
  #corMsub <- corM[localCenters,localCenters]
  localCenters=rownames(densityDf)
  moduleResult <- list()
  moduleResultbak <- list()
  count <- 1
  count2 <- 1
  #count2 <- 0
  j <- 0
  #totalCount <- dim(corM)[1]
  centerBlackList <- c()
  totalLength <- length(localCenters)
  # for (i in 1:totalLength) {
  for (center in localCenters) {
    if (center%in%centerBlackList) {
      print(paste("new", center, "is duplicated center!", sep=" "))
      next}
    j <- j+1
    genelist <- singleModuleDetectionSum(corM=corM, center=center, cutThresd=cutThresd, centerBlackList=centerBlackList, topCount=topGene)
    # move all other genes to noise gene
    # if (length(genelist)==1) {
    genelist2 <- colnames(corM)[(corM[center,] <= cutThresd)]
      # } else {
      #   genelist2 <- colnames(corM)[colSums(corM[genelist,] <= cutThresd) > 0]
      # }
    centerBlackList <- unique(c(centerBlackList, genelist, genelist2))
    # creat module
    if (length(genelist) >= minModuleGene){
      moduleResult[[count]] <- genelist
      count <- count + 1
    }
  }
  return(moduleResult)
}

transferListToDf <- function(moduleResult=moduleResult) {
  # transfer to dataframe
  moduleResultDf <- data.frame()
  for (i in 1:length(moduleResult)){
    #if (length(moduleResult[[i]]) < minModuleGene) {next}
    for (j in moduleResult[[i]]) { 
      moduleResultDf <- rbind(moduleResultDf, c(j, i))
    }
  }
  colnames(moduleResultDf) <- c("gene", "module")
  moduleResultDf <- moduleResultDf[!duplicated(moduleResultDf$gene),]
  rownames(moduleResultDf) <- moduleResultDf$gene
  return(moduleResultDf)
}

fullModuleDetection <- function(corM=corM, moduleResult=moduleResult, cutThresd=cutThresd) {
  moduleResultFull <- moduleResult
  count <- 1
  centerBlackList=c()
  for (i in 1:length(moduleResultFull)) { 
    genelist <- moduleResultFull[[count]]
    fullGenelist <- colnames(corM)[colSums(corM[genelist,] < cutThresd) > length(genelist)/4]
    fullGenelist <- fullGenelist[!fullGenelist%in%centerBlackList]
    moduleResultFull[[count]] <- fullGenelist
    centerBlackList <- c(centerBlackList, moduleResultFull)
    count <- count + 1
  }
  return(moduleResultFull)
}

correctLocalCenters <- function(moduleResultFull=moduleResultFull, densityDf=densityDf) {
  count <- 1
  newModuleResult <- list()
  for (i in 1:length(moduleResultFull)) { 
    genelist <- moduleResultFull[[count]]
    densityDfSub <- densityDf[genelist[genelist%in%rownames(densityDf)],]
    densityDfSub <- densityDfSub[order(densityDfSub$centrality, decreasing=T),]
    newModuleResult[[count]] <- rownames(densityDfSub)[1:10]
    count <- count + 1
  }
  return(newModuleResult)
}

fullModuleDetectionAllHIDE <- function(corM=dis_matrix){
  moduleResult <- list()
  count <- 1
  count2 <- 0
  j <- 0
  totalCount <- dim(corM)[1]
  centerBlackList <- c()
  #for (i in 1:length(table(marker_result$cluster))) {
  for (center in rownames(corM)) {
    #geneCluster <- marker_result[marker_result$cluster==i,]$gene
    #center <- centerDetection(corM, geneCluster)
    count2 <- count2 + 1
    if (count2%%as.integer(totalCount/10)==0){
      j <- j+1
      print(paste(j*10, "percent was finished..."))}
    # center <- i
    # print(center)
    result1 <- singleModuleDetection(corM=corM, overlapM=geneOverlap, center=center, centerBlackList=centerBlackList, topCount=100)
    genelist <- result1[["genelist"]]
    centerBlackList1 <- result1[["centerBlackList"]]
    result2 <- moduleShift(corM, genelist, centerBlackList=centerBlackList1)
    genelist <- result2[["genelist"]]
    centerBlackList2 <- result2[["centerBlackList"]]
    tempCenter <- centerDetection(corM, genelist)
    if (tempCenter%in%centerBlackList) {
      # print(paste(tempCenter, "is duplicated center!", sep=" "))
      next}
    centerBlackList <- c(centerBlackList, centerBlackList2)
    # print(genelist)
    if (length(genelist) >= 5){
      moduleResult[[count]] <- genelist
      count <- count + 1
    }
  }
  return(moduleResult)
}

localCenterDetection <- function(corM=dis_matrix, disThresd=disThresd){
  start <- disThresd/10
  by <- disThresd/10
  densityDf <- data.frame(row.names=rownames(corM))
  for (i in seq(start,disThresd,by=by)){
    densityDf <- cbind(densityDf, rowSums(corM<i))
  }
  colnames(densityDf) <- as.character(seq(start,disThresd,by=by))
  # densityDf <- densityDf[densityDf$`0.7` > 5,]
  densityDf <- densityDf[apply(densityDf, 1, max) > 2,]
  densityDf <- as.data.frame(densityDf)
  densityDf$centrality <- apply(densityDf, 1, function(x) {return(sum(log10(1+x) * (max(seq(start,disThresd+by,by=by))-seq(start,disThresd,by=by))))})
  densityDf <- densityDf[order(densityDf$centrality, decreasing=T),]
  # how to merge the duplicated
  # a <- data.frame(x=names(densityDf[1,]), y=as.integer(densityDf[1,]))
  # a <- a[1:13,]
  # a$x <- as.double(a$x)
  # a$x <- 1-a$x
  # a$y <- log2(1+a$y)
  # plot(a, type="b", col="blue", xlab="1-k", ylab="log2(1+count(i))", main="Centrality of CLDN6")
  # write.csv(densityDf, file = "densityDf.csv")
  return(rownames(densityDf))
}

coreMarkerDetection <- function(corM=dis_matrix, exprM=expr_log3, densityDf=densityDf, overlap=0.8) {
  #exprM2 <- apply(exprM[,rownames(densityDf)]>=1,2,function(x) {storage.mode(x) <- 'integer'; x})
  #geneOverlap <- t(exprM2) %*% exprM2 
  densityDf$gene <- rownames(densityDf)
  blackList <- c()
  for (i in 1:dim(densityDf)[1]) {
    tmp_gene_list <- c()
    geneA <- densityDf[i,]$gene
    for (j in 1:dim(densityDf)[1]) {
      if (i >= j) {next}
      geneB <- densityDf[j,]$gene
      overlapPercentRaw <- table(colSums(corM[c(geneA, geneB),] < 0.5))
      if (is.na(overlapPercentRaw["2"])) {next}
      overlapPercent <- overlapPercentRaw["2"]/(overlapPercentRaw["1"]+overlapPercentRaw["2"])
      if (overlapPercent >= overlap) { 
        tmp_gene_list <- c(tmp_gene_list, geneB)
        blackList <- c(blackList, geneB)
        }
    }
    if (geneA %in% blackList) {next}
    blackList <- c(blackList, geneA)
    print(tmp_gene_list)
  }

}

sortmoduleResult <- function(corM, moduleResultDfTop=moduleResultDfTop){
  meanCor <- c()
  for (i in unique(moduleResultDfTop$module)){
    tmp <- moduleResultDfTop[moduleResultDfTop$module==i,]$gene
    tmpCorM <- corM[tmp,tmp]
    tmpMeanCor <- mean(tmpCorM[row(tmpCorM) < col(tmpCorM)])
    meanCor <- c(meanCor, tmpMeanCor)
  }
  names(meanCor) <- unique(moduleResultDfTop$module)
  meanCor <- sort(meanCor, decreasing = F)
  moduleResultDfTop$module <- factor(moduleResultDfTop$module, levels = names(meanCor))
  moduleResultDfTop$module <- as.integer(moduleResultDfTop$module)
  moduleResultDfTop <- moduleResultDfTop[order(moduleResultDfTop$module),]
  return(moduleResultDfTop)
}

selectTopMarker <- function(expr=expr, moduleResultDf=moduleResultDf, densityDf=densityDf, top=10) {
  topMarkers <- c()
  for (i in unique(moduleResultDf$module)) {
    tmp <- moduleResultDf[moduleResultDf$module==i,]$gene
    tmpCount <- rowSums(expr[tmp,] > 0 )
    tmp <- names(tmpCount)[tmpCount < dim(expr)[2]*0.95]
    tmpTop <- rownames(densityDf[tmp,][1:top,])
    if (length(tmpTop)<10) {next}
    topMarkers <- c(topMarkers, tmpTop)
  }
  return(moduleResultDf[moduleResultDf$gene%in%topMarkers,])
}

cellOrderInference <- function(expr=exprZscore, module=moduleResultDf[moduleResultDf$module=="1",]$gene){
  # module <- moduleResult[[1]]
  # exprM <- logcounts(tmp_group)[module,]
  exprM <- expr[module,]
  # sign <- cor_matrix_spearman[module,module[1]] < 0 
  # exprM[sign,] <- exprM[sign,]*(-1)
  cellOrder <- sort(apply(exprM, 2, mean))
  # cellOrderdf <- data.frame(as.vector(cellOrder))
  # library(breakpoint)
  # breakpoint.loc <- CE.Normal.MeanVar(cellOrderdf, Nmax=1)$BP.Loc-1
  plot(cellOrder, xlab="sorted cells", ylab="mean expression of the module")
  # abline(v=breakpoint.loc, col="red")
  # pheatmap(logcounts(tmp_group)[moduleResult[[1]],names(cellOrder)], show_colnames = F, cluster_rows = T, cluster_cols = F)
  # return(list(cellOrder=names(cellOrder), breakpoint=breakpoint.loc))
  return(names(cellOrder))
}

classifyAndInference <- function(expr=expr, moduleResultDfTop=moduleResultDfTop, thresd=0.7, minClusterNum=5) {
  moduleType <- data.frame()
  inferCluster <- data.frame(cell=colnames(expr), row.names = colnames(expr))
  tmp_colname <- c("cell")
  for (i in unique(moduleResultDfTop$module)) {
    module <- moduleResultDfTop[moduleResultDfTop$module==i,]$gene
    if (length(module)<10) (next)
    exprM <- expr[module,]
    if (sum(colSums(exprM>0) < 0.8*dim(exprM)[1]) >= minClusterNum) {
      moduleType <- rbind(moduleType, data.frame(module=i, type="discrete"))
      tmpCluster <- colSums(exprM>0) >= thresd*dim(exprM)[1]
      inferCluster <- cbind(inferCluster, i=tmpCluster[inferCluster$cell])
      tmp_colname <- c(tmp_colname, i)
    } else {
      moduleType <- rbind(moduleType, data.frame(module=i, type="continuous"))
    }
  }
  colnames(inferCluster) <- tmp_colname
  inferCluster <- inferCluster[,-1]
  # pheatmap(expr[moduleResultDfTop$gene,], show_colnames = F, cluster_rows = F, cluster_cols = T, show_rownames = T)
  return(inferCluster)
}

selectIntermediateMarkers <- function(expr=expr, moduleResultDfTop=moduleResultDfTop, inferCluster=inferCluster, overlapThresd=5) {
  intermediateModule <- c()
  normalModule <- c()
  for (i in 1:dim(inferCluster)[2]) {
    tmp1df <- inferCluster[,i]
    tmp2df <- inferCluster[,-i]
    tmp2dfMerge <- as.vector(rowSums(tmp2df)>0)
    # new above 5, or more than half is new
    if (sum((tmp1df - tmp2dfMerge) == 1) >= overlapThresd | sum((tmp1df - tmp2dfMerge)==1 & tmp1df==1) / sum(tmp1df==1) > 0.3) {
      normalModule <- c(normalModule, i)
    } else {
      intermediateModule <- c(intermediateModule, i)
    }
  }
  # expr[moduleResultDfTop$gene,]
  clusterDf <- data.frame(row.names=colnames(expr))
  for (j in normalModule) {
    tmpModule <- colSums(expr[moduleResultDfTop[moduleResultDfTop$module==j,]$gene,] > 0)
    clusterDf <- cbind(clusterDf, as.data.frame(tmpModule))
  }
  colnames(clusterDf) <- normalModule
  # take care of the rare one
  # clusterDf <- clusterDf[,names(sort(colSums(clusterDf), decreasing = F))] # no use
  # final clusters
  finalCluster <- as.data.frame(apply(clusterDf, 1, function(x) names(sort(x, decreasing = T)[1])))
  finalCluster$cell <- rownames(clusterDf)
  colnames(finalCluster) <- c("Cluster", "cell")
  if (sum(rowSums(clusterDf)==0)) {finalCluster[rowSums(clusterDf)==0,]$Cluster <- "unknow"}
  finalCluster2 <- data.frame(Cluster=finalCluster$Cluster, row.names=finalCluster$cell)
  # add intermediate module
  finalCluster4 <- finalCluster2
  order1 <- unique(finalCluster4$Cluster)
  for (j in intermediateModule) {
    # judge which two clusters it belong to
    tmp1 <- inferCluster[,j]
    tmp2 <- inferCluster[,-j]
    origin2Cluster <- names(sort(apply(tmp2, 2, function(x) sum(x&tmp1)), decreasing = T)[1:2])
    origin2Cluster <- sort(origin2Cluster)
    #
    finalCluster3 <- finalCluster
    tmpInter <- tmp1 & (inferCluster[,origin2Cluster[1]] | inferCluster[,origin2Cluster[2]])
    finalCluster3[tmpInter,]$Cluster <- as.character(j)
    order2 <- c(order1[1:which(order1==origin2Cluster[1])], j, order1[which(order1==origin2Cluster[2]):length(order1)])
    finalCluster4 <- cbind(finalCluster4, Cluster=finalCluster3$Cluster)
  }
  colnames(finalCluster4) <- c("MainClusters", paste("intermediate", j, sep=""))
  finalCluster4$intermediate4 <- factor(finalCluster4$intermediate4, levels=order2)
  # cellOrder <- rownames(finalCluster4)[order(finalCluster4$intermediate4)]
  finalCluster4 <- finalCluster4[order(finalCluster4$intermediate4),]
  cellOrder <- localSort(finalCluster4, local=as.character(j))
}

localSort <- function(finalCluster4, local=as.character(j)) {
  order3 <- levels(finalCluster4$intermediate4)
  former <- order3[1:(which(order3==local)-1)]
  latter <- order3[(which(order3==local)+1):length(order3)]
  formerCells <- rownames(finalCluster4)[finalCluster4$intermediate4 %in% former]
  latterCells <- rownames(finalCluster4)[finalCluster4$intermediate4 %in% latter]
  interCells <- finalCluster4[finalCluster4$intermediate4==local,]
  interCells2 <-  rownames(interCells)[order(interCells$MainClusters)]
  cellOrder <- c(formerCells, interCells2, latterCells)
  return(cellOrder)
}

getFullMarkerModule <- function(dis_matrix2, moduleResultDfTop) {
  # classify all the genes
  fullGeneClusters <- data.frame(row.names = rownames(dis_matrix2))
  for (i in unique(moduleResultDfTop$module)) {
    tmpGene <- moduleResultDfTop[moduleResultDfTop$module==i,]$gene
    fullGeneClusters <- cbind(fullGeneClusters, as.data.frame(rowSums(dis_matrix2[, tmpGene])))
  }
  colnames(fullGeneClusters) <- unique(moduleResultDfTop$module)
  as.data.frame(apply(fullGeneClusters, 1, function(x) names(sort(x, decreasing = T)[1])))
  fullGeneClusters2 <- as.data.frame(apply(fullGeneClusters, 1, function(x) names(sort(x, decreasing = F)[1])))
  colnames(fullGeneClusters2) <- "module"
  # select the top 100 genes
  fullGeneClusters3 <- list()
  for (j in unique(moduleResultDfTop$module)) {
    tmpGene <- rownames(fullGeneClusters2)[fullGeneClusters2$module==j]
    tmpGene2 <- moduleResultDfTop[moduleResultDfTop$module==j,]$gene
    tmpTop100 <- names(sort(rowSums(dis_matrix2[tmpGene, tmpGene2]), decreasing=F)[1:100])
    fullGeneClusters3[[j]] <- tmpTop100
    save(fullGeneClusters3, file="fullGeneClusters3.Rdata")
  }
}

# pdf(sprintf('results/%s_maturation_trajectory.pdf', result.bn), width = 7, height = 5)

expr <- logcounts(tmp_group)
# moduleResultDf <- transferListToDf(moduleResult=newModuleResult)
moduleResultDfTop <- selectTopMarker(expr=expr, moduleResultDf=moduleResultDf, densityDf=densityDf, top=10)
moduleResultDfTop <- moduleResultDfTop[moduleResultDfTop$module%in% names(table(moduleResultDfTop$module))[table(moduleResultDfTop$module)==10],]
moduleResultDfTop <- sortmoduleResult(corM=corM, moduleResultDfTop=moduleResultDfTop)
# setwd('D:\\2.Code\\github\\MBSIT\\R')
# expr <- expr_matrix2
exprZscore <- (expr-apply(expr, 1, mean))/apply(expr, 1, sd)
# for (i in 1:length(moduleResult)){
pdf('moduleResult.pdf')
library(pheatmap)
for (i in unique(moduleResultDfTop$module)) {
  # pheatmap(expr[moduleResult[[i]],], show_colnames = F, cluster_rows = T)
  module <- moduleResultDfTop[moduleResultDfTop$module==i,]$gene
  cellOrder <- cellOrderInference(expr=exprZscore, module=module)
  #breakpoint.loc <- cellOrderInference(exprZscore, module)$breakpoint
  pheatmap(expr[module,cellOrder], show_colnames = F, cluster_rows = T, cluster_cols = F, show_rownames = T)
}
dev.off()
pheatmap(expr[moduleResultDfTop$gene,], show_colnames = F, cluster_rows = F, cluster_cols = T, show_rownames = T)
pheatmap(expr[moduleResultDfTop$gene,cellOrder], show_colnames = F, cluster_rows = F, cluster_cols = F, show_rownames = T, annotation_col = finalCluster4)

#######################
# main code
result <- outlierDetection(M=logcounts(tmp_group))
table(tmp_group[,result$cluster==0]$rename_sc3_4)
tmp_group <- tmp_group[,result$cluster!=0]

# 
expr_log3 <- t(logcounts(tmp_group))
# expr_log3 <- t(expr_matrix)

# overlap
# expr_log4 <- apply(expr_log3>=1,2,function(x) {storage.mode(x) <- 'integer'; x})
# geneOverlap <- t(expr_log4) %*% expr_log4 
# use.gene <- rownames(geneOverlap)[apply(geneOverlap, 2, max) > 3]
expr_log3 <- expr_log3[,colSums(expr_log3>0)>3]

# library(WGCNA)
# standard deviation can't be zero
# cor_matrix_pearson <- WGCNA::cor(x = as.matrix((expr_log3)), method = "pearson")
cor_matrix_spearman <- cor(x = as.matrix((expr_log3)), method = "spearman")

cor_matrix_spearman[cor_matrix_spearman<0] <- 0
dis_matrix <- 1 - abs(cor_matrix_spearman)

dis_matrix[is.na(dis_matrix)] <- 1
dis_matrix[row(dis_matrix)==col(dis_matrix)] <- 1
# dis_matrix <- dis_matrix[use.gene, use.gene]
minPerRow <- apply(dis_matrix, 2, min)
# consider how to evaluate this value
disThresd <- quantile(minPerRow, probs = seq(0, 1, 0.1))[2]
#         0%        25%        50%        75%       100% 
# 0.06675682 0.65125326 0.72883328 0.76764952 0.82801885
dis_matrix2 <- dis_matrix
dis_matrix <- dis_matrix[names(minPerRow[minPerRow < disThresd]),names(minPerRow[minPerRow < disThresd])]
# dis_matrix_dist <- as.dist(dis_matrix)


## DBSCAN
## find most local centers
eps <- 0
start <- 0
end <- 1
count <- 0
while (T) {  
  count <- count + 1
  print(paste("Round:", count))
  epsPair <- epsDetectionCor(dis_matrix_dist, start=start, end=end)
  if (length(epsPair)==2){
    start <-  epsPair[1]
    end <-  epsPair[2]
  } else if (length(epsPair)==1){
    eps <- epsPair
    break
  }
}
minPts <- 5
res <- dbscan(dis_matrix_dist, eps = eps, minPts = minPts) 
marker_result <- data.frame(gene=rownames(dis_matrix),cluster=res$cluster,row.names = rownames(dis_matrix))
marker_result <- marker_result[marker_result$cluster!=0,]
marker_result$cluster <- as.integer(as.factor(marker_result$cluster))

## Hierarchical Clustering
# http://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
hc <- fastcluster::hclust(dis_matrix_dist, method="average")
# plot(hc)
clusterNum <- c()
for (i in 0:100){
  groups<-cutree(hc, h=(i/100))
  tempNum <- sum(table(groups)>=5)
  names(tempNum) <- i
  clusterNum <- c(clusterNum, tempNum)
}
plot(clusterNum)
# return(names(sort(clusterNum, decreasing = T)[1]))
i <- 92
groups<-cutree(hc, h=(i/100))

marker_result <- data.frame(gene=names(groups),cluster=as.vector(groups),row.names = names(groups))
filter_cluster <- names(table(marker_result$cluster))[table(marker_result$cluster)>=5]
marker_result <- marker_result[marker_result$cluster%in%filter_cluster,]
marker_result$cluster <- as.integer(as.factor(marker_result$cluster))

genelist <- marker_result[marker_result$cluster==1,]$gene


test <- function(){
  mydata <- t(logcounts(tmp_group))
  mydata <- mydata[,colSums(mydata)>0]

  #principal component analysis
  prin_comp <- prcomp(mydata, scale. = T, center = T)
  names(prin_comp) 
  #outputs the mean of variables
  prin_comp$center 
  #outputs the standard deviation of variables
  prin_comp$scale 
  #Each column of rotation matrix contains the principal component loading vector
  prin_comp$rotation 
  #resultant principal components
  # biplot(prin_comp, scale = 0)
  #compute standard deviation of each principal component
  std_dev <- prin_comp$sdev
  #compute variance
  pr_var <- std_dev^2
  #check variance of first 10 components
  pr_var[1:10]
  #proportion of variance explained
  prop_varex <- pr_var/sum(pr_var)
  prop_varex[1:20]
  #scree plot
  plot(prop_varex, xlab = "Principal Component",
         ylab = "Proportion of Variance Explained",
         type = "b")
  #cumulative scree plot
  plot(cumsum(prop_varex), xlab = "Principal Component",
         ylab = "Cumulative Proportion of Variance Explained",
         type = "b")
  plot(prin_comp$x[,1:2])

  library(pcaMethods)
  winePCAmethods <- pca(mydata, scale = "uv", center = T, nPcs = 100, method = "svd")
  # slplot(winePCAmethods, scoresLoadings = c(T,T), scol = wineClasses)
  prop_varex2 <- winePCAmethods@sDev^2/sum(winePCAmethods@sDev^2)
  plot(cumsum(prop_varex2), xlab = "Principal Component",
       ylab = "Cumulative Proportion of Variance Explained",
       type = "b")
  plot(winePCAmethods@scores[,1:2])

  ## scater
  prin_comp2 <- scater::plotPCA(tmp_group, exprs_values="logcounts", colour_by="rename_sc3_4")
    
}

