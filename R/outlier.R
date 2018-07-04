setwd("/Users/surgery/Project/HOME/github/MBSIT/R")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")

logcounts(tmp_group)

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

bestEpsDetection <- function(corMatrix=dis_matrix_pearson2, start=1, end=100, fold=5, primary=T, minPts=6, percent=0.05, plot=T){
  # DBSCAN
  library(dbscan)
  corMatrix[is.na(corMatrix)] <- 1
  corMatrixDis <- as.dist(corMatrix)
  if (primary==T){ by = 1 } else {by=(end-start)/100}
  epsLevels <- (seq(start, end, by=by)-1)*fold
  allclusterCount <- c()
  alleps <- c()
  for ( eps in epsLevels){
    #eps <- 5*(i-1)
    alleps <- c(alleps, eps)
    # minPts <- 10
    res <- dbscan(corMatrixDis, eps = eps, minPts = minPts) 
    clusterCount <- length(table(res$cluster))
    allclusterCount <- c(allclusterCount, clusterCount)
  }
  names(allclusterCount) <- alleps
  allclusterCount <- allclusterCount[!is.na(allclusterCount)]
  allclusterCount <- allclusterCount[allclusterCount!=dim(corMatrixDis)[1]]
  if (plot) {plot(allnoiseCount)}
  if (primary==T)
  {return(c(as.integer(names(allnoiseCount)[1]), as.integer(names(allnoiseCount)[length(allnoiseCount)])))}
  else {return(as.double(names(sort(allnoiseCount[allnoiseCount>=percent*dim(corMatrixDis)[1]])[1])))}
  # {return(as.double(names(sort(allnoiseCount %% as.integer(percent*dim(corMatrix)[1]))[1])))}
}

outlierDetection <- function(M=logcounts(tmp_group), pcNum=100, method="prcomp", minPts=10, plot=T){
  exprMatrix <- t(M)
  exprMatrix <- exprMatrix[,colSums(exprMatrix)>0]
  pcMatrix <- 0
  start_time <- Sys.time()
  if (method=="prcomp"){
    prin_comp <- prcomp(exprMatrix, scale. = T, center = T)
    pcMatrix <- prin_comp$x[,1:pcNum]
    } else if (method=="pcaMethods"){
      library(pcaMethods)
      prin_comp <- pca(exprMatrix, scale = "uv", center = T, nPcs = pcNum, method = "svd")
      pcMatrix <- prin_comp@scores[,1:pcNum]
    }
    end_time <- Sys.time()
    print("matrix size: ")
    print(dim(M)) 
    print(paste("pcNum: ", pcNum, ",method: ", method))
    print(end_time - start_time)
    # return(pcMatrix)
    eps <- epsDetection(pcMatrix)
    eps2 <- epsDetection(pcMatrix, start=eps[1], end = eps[2], primary = F, fold = 1)
    res <- dbscan(pcMatrix, eps = eps2, minPts = minPts) 
    # rownames(winePCAmethods@scores)[res$cluster==0]
    pca <- as.data.frame(pcMatrix)
    pca$cluster <- as.character(res$cluster)
    if (plot) {
      library(ggplot2)
      print(ggplot(pca, aes(x=PC1, y=PC2, color=cluster)) + geom_point())}
    return(pca)
}

result <- outlierDetection(M=logcounts(tmp_group))
tmp_group <- tmp_group[,result$cluster!=0]
expr_log3 <- t(logcounts(tmp_group))
library(WGCNA)
# standard deviation can't be zero
cor_matrix_pearson <- WGCNA::cor(x = as.matrix((expr_log3)), method = "pearson")
cor_matrix_spearman <- WGCNA::cor(x = as.matrix((expr_log3)), method = "spearman")

dis_matrix_pearson2 <- 1 - abs(cor_matrix_pearson)
dis_matrix_pearson2[is.na(dis_matrix_pearson2)] <- 1
# dis_matrix_pearson2 <- round(dis_matrix_pearson2,3)
# dis_all <- c()
# for (i in 1:dim(dis_matrix_pearson2)[1]){
#   for (j in 1:dim(dis_matrix_pearson2)[2])
#     if (i > j) {
#       dis_all <- c(dis_all, dis_matrix_pearson2[i,j])
#     } else if (i == j) {
#       dis_matrix_pearson2[i,j] <- 1
#     }
# }
## Hierarchical Clustering
# http://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
dis_matrix_pearson3 <- as.dist(dis_matrix_pearson2)
hc <- fastcluster::hclust(dis_matrix_pearson3, method="average")
plot(hc)
clusterNum <- c()
for (i in 0:100){
  groups<-cutree(hc, h=(i/100))
  clusterNum <- c(clusterNum, sum(table(groups)>=5))
  plot(clusterNum)
  #marker_result <- data.frame(gene=names(groups),cluster=as.vector(groups),row.names = names(groups))
  #filter_cluster <- names(table(marker_result$cluster))[table(marker_result$cluster)>=5]
  #marker_result <- marker_result[marker_result$cluster%in%filter_cluster,]
  #plot_markers <- marker_result[marker_result$cluster==4,]
  #plot_markers <- plot_markers[plot_markers$gene%in%rownames(tmp_group),]
}



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

