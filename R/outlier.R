setwd("/Users/surgery/Project/HOME/github/MBSIT/R")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")
options(stringsAsFactors = F)

cellCountPerGene <- rowSums(counts(tmp_group)>=5)
geneCountPerCell <- colSums(counts(tmp_group)>=5)
totalReadCount <- colSums(counts(tmp_group))
cvPerGene <- apply(counts(tmp_group), 1, sd)/apply(counts(tmp_group), 1, mean)
# cvPerCell <- apply(counts(tmp_group), 2, sd)/apply(counts(tmp_group), 2, mean)

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
fullModuleDetectionSum <- function(corM=dis_matrix, localCenters=densityDf, minModuleGene=5){
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
      print(paste(center, "is duplicated center!", sep=" "))
      next}
    j <- j+1
    genelist <- singleModuleDetectionSum(corM=corM, center=center, cutThresd=0.6, centerBlackList=centerBlackList, topCount=100)
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
  gene <- a
  moduleResultDf[gene[gene%in%moduleResultDf$gene],]
  return(moduleResultDf)
}

fullModuleDetectionAll <- function(corM=dis_matrix){
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

highDensityCenterDetection <- function(corM=dis_matrix){
  densityDf <- data.frame(row.names=rownames(corM))
  for (i in seq(0.05,0.7,by=0.05)){
    densityDf <- cbind(densityDf, rowSums(corM<i))
  }
  colnames(densityDf) <- as.character(seq(0.05,0.7,by=0.05))
  densityDf <- densityDf[densityDf$`0.7` > 5,]
  densityDf <- densityDf[names(sort(rowSums(densityDf), decreasing = T)),]
  densityDf <- densityDf[names(sort(rowSums(densityDf > 0), decreasing = T)),]
  return(rownames(densityDf))
}

sortmoduleResult <- function(corM, moduleResult=moduleResult){
  meanCor <- c()
  for (i in 1:length(moduleResult)){
    meanCor <- c(meanCor, mean(corM[moduleResult[[i]], moduleResult[[i]]]))
  }
  names(meanCor) <- 1:length(moduleResult)
  meanCor <- sort(meanCor)
}

cellOrderInference <- function(expr=exprZscore, module=moduleResultDf[moduleResultDf$module=="1",]$gene){
  # module <- moduleResult[[1]]
  # exprM <- logcounts(tmp_group)[module,]
  exprM <- expr[module,]
  # sign <- cor_matrix_spearman[module,module[1]] < 0 
  # exprM[sign,] <- exprM[sign,]*(-1)
  cellOrder <- sort(colSums(exprM))
  plot(cellOrder)
  # pheatmap(logcounts(tmp_group)[moduleResult[[1]],names(cellOrder)], show_colnames = F, cluster_rows = T, cluster_cols = F)
  return(names(cellOrder))
}

# pdf(sprintf('results/%s_maturation_trajectory.pdf', result.bn), width = 7, height = 5)

expr <- logcounts(tmp_group)
exprZscore <- (expr-apply(expr, 1, mean))/apply(expr, 1, sd)
# for (i in 1:length(moduleResult)){
pdf('moduleResult.pdf')
library(pheatmap)
for (i in unique(moduleResultDf$module)) {
  # pheatmap(expr[moduleResult[[i]],], show_colnames = F, cluster_rows = T)
  module <- moduleResultDf[moduleResultDf$module==i,]$gene
  cellOrder <- cellOrderInference(exprZscore, module)
  pheatmap(expr[module,cellOrder], show_colnames = F, cluster_rows = T, cluster_cols = F)
}
dev.off()



result <- outlierDetection(M=logcounts(tmp_group))
table(tmp_group[,result$cluster==0]$rename_sc3_4)
tmp_group <- tmp_group[,result$cluster!=0]
expr_log3 <- t(logcounts(tmp_group))
# overlap
expr_log4 <- apply(expr_log3>=1,2,function(x) {storage.mode(x) <- 'integer'; x})
geneOverlap <- t(expr_log4) %*% expr_log4 
use.gene <- rownames(geneOverlap)[apply(geneOverlap, 2, max) > 3]

library(WGCNA)
# standard deviation can't be zero
# cor_matrix_pearson <- WGCNA::cor(x = as.matrix((expr_log3)), method = "pearson")
cor_matrix_spearman <- WGCNA::cor(x = as.matrix((expr_log3)), method = "spearman")

cor_matrix_spearman[cor_matrix_spearman<0] <- 0
dis_matrix <- 1 - abs(cor_matrix_spearman)

dis_matrix[is.na(dis_matrix)] <- 1
dis_matrix[row(dis_matrix)==col(dis_matrix)] <- 1
dis_matrix <- dis_matrix[use.gene, use.gene]
minPerRow <- apply(dis_matrix, 2, min)
disThresd <- quantile(minPerRow, probs = seq(0, 1, 0.25))[2]
dis_matrix <- dis_matrix[names(minPerRow[minPerRow < disThresd]),names(minPerRow[minPerRow < disThresd])]
dis_matrix_dist <- as.dist(dis_matrix)

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

