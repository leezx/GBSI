setwd("/Users/surgery/Project/HOME/github/MBSIT/R")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")

logcounts(tmp_group)

cellCountPerGene <- rowSums(counts(tmp_group)>=5)
geneCountPerCell <- colSums(counts(tmp_group)>=5)
totalReadCount <- colSums(counts(tmp_group))
cvPerGene <- apply(counts(tmp_group), 1, sd)/apply(counts(tmp_group), 1, mean)
# cvPerCell <- apply(counts(tmp_group), 2, sd)/apply(counts(tmp_group), 2, mean)
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

# DBSCAN
library(dbscan)
allnoiseCount <- c()
alleps <- c()
for ( i in 1:100){
  eps <- 5*(i-1)
  alleps <- c(alleps, eps)
  minPts <- 10
  res <- dbscan(winePCAmethods@scores, eps = eps, minPts = minPts) 
  noiseCount <- table(res$cluster)["0"]
  allnoiseCount <- c(allnoiseCount, noiseCount)
}
names(allnoiseCount) <- alleps

mylevels <- seq(40, 165, by=(165-40)/100)
allnoiseCount <- c()
alleps <- c()
for ( eps in mylevels){
  #eps <- 5*(i-1)
  alleps <- c(alleps, eps)
  minPts <- 10
  res <- dbscan(winePCAmethods@scores, eps = eps, minPts = minPts) 
  noiseCount <- table(res$cluster)["0"]
  allnoiseCount <- c(allnoiseCount, noiseCount)
}
names(allnoiseCount) <- alleps
plot(allnoiseCount)

res <- dbscan(winePCAmethods@scores, eps = 126.25, minPts = 10) 
# rownames(prin_comp$x)[res$cluster==0]
rownames(winePCAmethods@scores)[res$cluster==0]
pca <- as.data.frame(winePCAmethods@scores)
pca$cluster <- res$cluster
ggplot(pca, aes(x=PC1, y=PC2, color=cluster)) + geom_point()

