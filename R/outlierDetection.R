setwd("/Users/surgery/Project/HOME/github/MBSIT/R")
# load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")
options(stringsAsFactors = F)
source("/Users/surgery/Project/HOME/myScript/zxli_lib.R") # my R lib

# download real dataset
# https://media.nature.com/original/nature-assets/nmeth/journal/v14/n4/extref/nmeth.4207-S2.zip
# git clone https://github.com/hemberg-lab/scRNA.seq.course.git

load("/Users/surgery/Downloads/Supplementary_Software/R/data/Test_3_Pollen.RData")
# Test_3_Pollen$results$y$cluster)
# Test_3_Pollen$true_labs
# Test_3_Pollen$in_X

deng <- readRDS("/Users/surgery/Project/HOME/github/scRNA.seq.course/deng/deng-reads.rds")
scater::plotPCA(deng, exprs_values = "logcounts", colour_by="cell_type1")

#expr <- logcounts(tmp_group)
#expr <- Test_3_Pollen$in_X
expr <- logcounts(deng)
expr <- expr[rowSums(expr>0)>3,]
# expr <- expr[rowSums(expr)>0,]

cellCountPerGene <- rowSums(Test_3_Pollen$in_X>0)
geneCountPerCell <- colSums(Test_3_Pollen$in_X>0)
totalReadCount <- colSums(Test_3_Pollen$in_X)
cvPerGene <- apply(Test_3_Pollen$in_X, 1, sd)/apply(Test_3_Pollen$in_X, 1, mean)

# get markers
library(SC3) # conflict with WGCNA
all_markers_raw <- get_marker_genes(logcounts(deng), deng$cell_type1)
all_markers <- all_markers_raw[all_markers_raw$auroc>0.95 & !is.na(all_markers_raw$auroc) & all_markers_raw$pvalue < 0.001,]
sc3_marker <- data.frame(name=rownames(deng)[as.numeric(rownames(all_markers))], cluster=all_markers$clusts)
rownames(sc3_marker) <- sc3_marker$name
annotation_col <- data.frame(Cluster = factor(deng$cell_type1), row.names = colnames(deng))
library(pheatmap)
markers <- sc3_marker[sc3_marker$cluster=="zygote",]$name
marker2 <- sort(colSums(abs(corMatrix[markers, markers])), decreasing = T)[1:100]
notmakers <- rownames(corMatrix)[!rownames(corMatrix)%in%dengMarkers$Gene][1:2000]
notmakers2 <- sort(colSums(abs(corMatrix[notmakers, notmakers])), decreasing = F)[1:100]
pheatmap(logcounts(deng)[names(marker2), order(deng$cell_type1)], cluster_rows = T, show_rownames = F, show_colnames = F, cluster_cols = F, annotation_col = annotation_col)
# remove
a <- logcounts(deng)[sc3_marker[sc3_marker$cluster=="2cell",]$name, order(deng$cell_type1)]
a <- a[rowSums(a>0)<100,]
b <- corMatrix[rownames(a), rownames(a)]
b <- b[apply(b, 2, min)>0.50,]

# outlier detection
pca <- outlierDetection(M=expr, minPts = 5, percent = 0.02)

pca$trueCluster <- as.character(Test_3_Pollen$true_labs$V1)
pca$trueCluster <- as.character(deng$cell_type1)
ggplot(pca, aes(x=PC1, y=PC2, color=cluster)) + geom_point()
# SIMLRcor <- as.data.frame(Test_3_Pollen$results$ydata)
# SIMLRcor$outlier <- pca$cluster
# SIMLRcor$trueCluster <- as.character(Test_3_Pollen$true_labs$V1)
# ggplot(SIMLRcor, aes(x=V1, y=V2, color=outlier)) + geom_point()

# simulate outliers
for (i in seq(1,15,by=2)) {
  print(i)
  deng1 <- deng
  cells <- colnames(deng)[deng$cell_type1=="2cell" | deng$cell_type1=="zygote"]
  count <- i
  # as.integer(length(cells)*i/100)
  randomCells <- sample(cells, count)
  randomGenes <- sample(rownames(deng1), 1000)
  logcounts(deng1)[randomGenes,randomCells] <- 100 # count%
  #
  expr <- logcounts(deng1)
  #expr <- expr[rowSums(expr>0)>3,]
  # print(scater::plotPCA(deng1, exprs_values = "logcounts", colour_by="cell_type1"))
  corMatrix2 <- buildGeneNetworkCor(expr=expr, method = "pearson")
  cor1 <- ((corMatrix2[names(notmakers2), names(notmakers2)]))
  cor1 <- cor1[row(cor1)>col(cor1)]
  cor2 <- ((corMatrix2[names(marker2), names(notmakers2)]))
  cor2 <- cor2[row(cor2)>col(cor2)]
  cor3 <- ((corMatrix2[names(marker2), names(marker2)]))
  cor3 <- cor3[row(cor3)>col(cor3)]
  # total_pearson <- data.frame(cor1, cor2, cor3)
  total_pearson <- cbind(total_pearson, cor1, cor2, cor3)
}
save(total_pearson, file = "total_pearson.Rdata")
barplot(apply(total_pearson, 2, mean))

detach("package:SC3", unload=TRUE)
library(WGCNA)
corMatrix <- buildGeneNetworkCor(expr=expr, method = "pearson")
# before
cor1 <- (sample(abs(corMatrix[!rownames(corMatrix)%in%dengMarkers$Gene, !colnames(corMatrix)%in%dengMarkers$Gene]), 1000))
cor2 <- (sample(abs(corMatrix[rownames(corMatrix)%in%dengMarkers.sig1$Gene, !colnames(corMatrix)%in%dengMarkers$Gene]), 1000))
cor3 <- (sample(abs(corMatrix[rownames(corMatrix)%in%dengMarkers.sig1$Gene, colnames(corMatrix)%in%dengMarkers.sig1$Gene]), 1000))

total_pearson <- data.frame(ctrl.cor1=cor1, ctrl.cor2=cor2, ctrl.cor3=cor3)
total_pearson <- cbind(percent1=total_pearson[,1:3], percent3=data.frame(cor1, cor2, cor3))
apply(total_pearson, 2, mean)

# my dataset
pca <- outlierDetection(M=logcounts(tmp_group), minPts = 5, percent = 0.01, pcNum = 2)
corMatrix <- buildGeneNetworkCor(expr=logcounts(tmp_group), method = "pearson")
#corMatrix2 <- buildGeneNetworkCor(expr=logcounts(tmp_group)[,rownames(pca)[pca$cluster=="1"]], , method = "pearson")
cor4 <- as.vector(corMatrix[a,a])
cor5 <- as.vector(corMatrix2[a,a])

# check deng known markers
dengMarkers <- read.csv("/Users/surgery/Project/HOME/github/MBSIT/SC3/nmeth.4236-S3.csv", header = T)
rownames(dengMarkers) <- dengMarkers$Gene
# dengMarkers.sig <- dengMarkers[dengMarkers$AUC>0.9 & dengMarkers$p.value<0.01,]
# dengMarkers.sig1 <- dengMarkers.sig[dengMarkers.sig$clusts=="7",]

buildGeneNetworkCor <- function(expr=Test_3_Pollen$in_X, thred=0, method="spearman", leastOverlapCell=3){
  # exprInteger <- apply(t(expr)>thred,2,function(x) {storage.mode(x) <- 'integer'; x})
  # geneOverlap <- t(exprInteger) %*% exprInteger 
  # use.gene <- rownames(geneOverlap)[apply(geneOverlap, 2, max) > 3]
  # expr <- expr[use.gene,]
  library(WGCNA)
  # standard deviation can't be zero
  corMatrix <- WGCNA::cor(x = as.matrix(t(expr)), method = method)
  # return(list(corMatrix=corMatrix, use.gene=use.gene))
  return(corMatrix)
}

# expr2 <- expr[,!colnames(expr)%in%c(c4n, c8n, c11n)]
save(expr2, file = "expr2.Rdata")
truel <- Test_3_Pollen$true_labs$V1[!colnames(expr)%in%c(c4n, c8n, c11n)]
simlrl <- Test_3_Pollen$results$y$cluster[!colnames(expr)%in%c(c4n, c8n, c11n)]

outlierDetection <- function(M=expr, pcNum=2, method="prcomp", minPts=4, percent=0.05, threads=1, plot=T){
  exprMatrix <- t(M)
  exprMatrix <- exprMatrix[,colSums(exprMatrix)>0]
  pcMatrix <- 0
  start_time <- Sys.time()
  if (method=="prcomp"){
    prin_comp <- prcomp(exprMatrix, scale. = T, center = T, rank. = pcNum)
    pcMatrix <- prin_comp$x
    #if (pcNum < dim(prin_comp$x)[2]) {
    #  pcMatrix <- prin_comp$x[,1:pcNum]
    #} else if (pcNum > dim(prin_comp$x)[2]) {
    #  print("ERROR: please choose a smaller pcNum!!!")
    #} else { pcMatrix <- prin_comp$x }
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
  # parallelDist
  library(parallelDist)
  # dis_matrix <- parDist(x = as.matrix(pcMatrix), method = "euclidean", threads=threads)
  dis_matrix <- parDist(x = as.matrix(pcMatrix), method = "mahalanobis", threads=threads)
  dis_matrix <- as.matrix(dis_matrix)
  dis_matrix[is.na(dis_matrix)] <- 0
  rownames(dis_matrix) <- rownames(pcMatrix)
  colnames(dis_matrix) <- rownames(dis_matrix)
  disOrder <- sort(apply(dis_matrix, 2, function(x) {return(mean(sort(x)[(minPts+1)]))}), decreasing=T)
  # plot(disOrder)
  # hist(disOrder, breaks = 20)
  outlier <- disOrder[1:(length(disOrder)*percent)]
  pca <- as.data.frame(pcMatrix)
  pca$outlier <- "NO"
  pca[names(outlier),]$outlier <- "YES"
  pca$outlier <- as.character(pca$outlier)
  # for simulate
  pca$groundTruth <- as.character(truel)
  pca$mvoutlier <- "NO"
  pca[!as.logical(outliers$wfinal01),]$mvoutlier <- "YES"
  # pca$SIMLR <- as.character(simlrl)
  # get x smallest 4: 1
  # c4 <- pca[pca$true=="4",]
  # c4n <- rownames(c4[order(c4$PC1, decreasing = F),][2:26,])
  # get y largest 8: 2
  # c8 <- pca[pca$true=="8",]
  # c8n <- rownames(c8[order(c8$PC2, decreasing = T),][3:42,])
  # get x largest 11: 5
  # c11 <- pca[pca$true=="11",]
  # c11n <- rownames(c11[order(c11$PC1, decreasing = T),][5:24,])
  # simulate data, remove some cells
  # pca2 <- pca[!rownames(pca)%in%c(c4n, c8n, c11n),]
  # DBSCAN
  # eps <- epsDetection(pcMatrix)
  #eps2 <- epsDetection(pcMatrix, start=eps[1], end = eps[2], primary = F, fold = 1)
  #res <- dbscan(pcMatrix, eps = eps2, minPts = minPts) 
  # rownames(winePCAmethods@scores)[res$cluster==0]
  #pca <- as.data.frame(pcMatrix)
  if (plot) {
    library(ggplot2)
    library(RColorBrewer)
    print(ggplot(pca, aes(x=PC1, y=PC2, color=outlier)) + 
            geom_point(size=2.5, alpha=1) + 
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
            scale_color_manual(values=brewer.pal(11,"Paired")[1:11]) )
    } # for cluster, c("#B2DF8A","#1F78B4"); #
  return(pca)
}

expr3 <- expr2[rowSums(expr2>0)>3600,]

expr3 <- t(Test_3_Pollen$in_X)
#expr3 <- t(logcounts(deng))
# expr3 <- expr3[,colSums(expr3>0) > 7]
#x.mad = apply(expr3, 2, mad)
#expr3 <- expr3[,(x.mad!=0)]
outliers <- mvoutlier::pcout(expr3, makeplot = FALSE,
                             explvar = 0.5, crit.M1 = 0.9,
                             crit.c1 = 5, crit.M2 = 0.9,
                             crit.c2 = 0.99, cs = 0.25,
                             outbound = 0.05)
!as.logical(outliers$wfinal01)

pca <- prcomp(t(Test_3_Pollen$in_X), scale. = F, center = F)
total.var <- sum(pca$sdev ^ 2)
percentVar <- pca$sdev ^ 2 / total.var
x <- sum(cumsum(percentVar) < 0.9)
#x1 <- cumsum(percentVar)[x+1]
y <- 0.9
# Cumulative Proportion of 
plot(cumsum(percentVar), xlab = "Principal Component",
              ylab = "Cumulative Variance Explained",
              type = "b")
abline(h=0.9, col="blue")
abline(v=x, col="blue")

pcax <- as.data.frame(pca$x)
pcax$trueCluster <- as.character(Test_3_Pollen$true_labs$V1)
pcax$mvoutlier <- "NO"
pcax[!as.logical(outliers$wfinal01),]$mvoutlier <- "YES"
library(ggplot2)
library(RColorBrewer)
pcax$trueCluster <- factor(pcax$trueCluster, levels=as.character(1:11))
print(ggplot(pcax, aes(x=PC1, y=PC2, color=trueCluster)) + 
        geom_point(size=2.5, alpha=1) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
        scale_color_manual(values=brewer.pal(11,"Paired")[1:11]) )


a <- t(Test_3_Pollen$in_X) %*% pca$rotation
library(MASS)
b <- pca$x %*% ginv(pca$rotation)
## add outlier to PCs
# outlierdf <- data.frame(PC1=c(-40,-40,-40,-200,-210,-220, rep(-250, 6), rep(10,6), rep(-100, 6)), PC2=c(25,0,5,-50,-40,-50, seq(30, 80, by=10), -seq(50, 100, by=10), seq(50, 100, by=10)))
outlierdf <- matrix(rep(0, 6982)*24, nrow = 6982, ncol = 24)
for (i in 1:24) {
  randomRow <- sample(1:6982, 400)
  outlierdf[randomRow, i] <- seq(1,1000, by=2.5)
}
colnames(outlierdf) <- 250:273
rownames(outlierdf) <- rownames(pollen)
d <- cbind(logcounts(pollen), outlierdf)
# pcax <- as.data.frame(pca$x)
c <- cbind(outlierdf, matrix(rep(0,247)*24, nrow = 24, ncol = 249))
colnames(c) <- colnames(pca$x)
d <- rbind(pca$x,c)
pcax2 <- as.data.frame(d)
pcax2$trueCluster <- "simulated outlier"
pcax2[rownames(pcax),]$trueCluster <- pcax$trueCluster

pcax2 <- pcax2[pcax2$trueCluster!="11",]
pcax2 <- pcax2[pcax2$trueCluster!="10",]
print(ggplot(pcax2, aes(x=PC1, y=PC2, color=trueCluster)) + 
        geom_point(size=2.5, alpha=1) + 
        theme_bw() + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
        scale_color_manual(values=brewer.pal(11,"Paired")[1:10]) )

recoverExpr <- as.matrix(d) %*% ginv(pca$rotation)
recoverExpr[recoverExpr<0] <- 0
save(recoverExpr, file="recoverExpr.Rdata")

# SC3 get marker
colData <- data.frame(sampleNames=rownames(pcax), cell_group=pcax$trueCluster, mvoutlier=pcax$mvoutlier, row.names=rownames(pcax))
rowData <- data.frame(feature_symbol=rownames(Test_3_Pollen$in_X), gene_short_name=rownames(Test_3_Pollen$in_X), row.names = rownames(Test_3_Pollen$in_X))
data <- Test_3_Pollen$in_X
#rownames(data) <- rownames(rowData)
#colnames(data) <- rownames(colData)
pollen <- SingleCellExperiment(assays = list(counts = data), colData = colData, rowData = rowData)
logcounts(pollen) <- log2(counts(pollen)+1)
pollen <- calculateQCMetrics(pollen)
library(SC3)
all_markers_raw <- get_marker_genes(logcounts(pollen), pollen$cell_group)
all_markers <- all_markers_raw[all_markers_raw$auroc>0.8 & !is.na(all_markers_raw$auroc) & all_markers_raw$pvalue < 0.05,]
sc3_marker <- data.frame(name=rownames(pollen)[as.numeric(rownames(all_markers))], cluster=all_markers$clusts)
rownames(sc3_marker) <- sc3_marker$name
annotation_col <- data.frame(Cluster = factor(pollen$cell_group), row.names = colnames(pollen))
library(pheatmap)
markers <- sc3_marker[sc3_marker$cluster=="8",]$name
corMatrix <- cor(x = as.matrix(t(logcounts(pollen))), method = "pearson")
marker2 <- sort(colSums(abs(corMatrix[markers, markers])), decreasing = T)[1:100]
notmakers <- rownames(corMatrix)[!rownames(corMatrix)%in%sc3_marker$name][1:2000]
notmakers2 <- sort(colSums(abs(corMatrix[notmakers, notmakers])), decreasing = F)[1:100]
pheatmap(logcounts(pollen)[names(marker2), order(pollen$cell_group)], cluster_rows = T, show_rownames = F, show_colnames = F, cluster_cols = F, annotation_col = annotation_col)

cor1 <- ((corMatrix[names(notmakers2), names(notmakers2)]))
cor1 <- cor1[row(cor1)>col(cor1)]
cor2 <- ((corMatrix[names(marker2), names(notmakers2)]))
cor2 <- cor2[row(cor2)>col(cor2)]
cor3 <- ((corMatrix[names(marker2), names(marker2)]))
cor3 <- cor3[row(cor3)>col(cor3)]
total_pearson <- data.frame(cor1, cor2, cor3)

save(d, file="simulated_outlier.Rdata")
for (i in seq(1,15,by=2)) {
  corMatrix2 <- cor(x = as.matrix(t(d[,1:(249+i)])), method = "pearson")
  corMatrix2 <- as.data.frame(corMatrix2)
  colnames(corMatrix2) <- 1:(dim(corMatrix2)[2])
  cor1 <- ((corMatrix2[names(notmakers2), names(notmakers2)]))
  cor1 <- cor1[row(cor1)>col(cor1)]
  cor2 <- ((corMatrix2[names(marker2), names(notmakers2)]))
  cor2 <- cor2[row(cor2)>col(cor2)]
  cor3 <- ((corMatrix2[names(marker2), names(marker2)]))
  cor3 <- cor3[row(cor3)>col(cor3)]
  temp <- data.frame(cor1=cor1, cor2=cor2, cor3=cor3)
  # cbind(temp, level=i)
  total_pearson <- rbind(total_pearson, cbind(temp, level=i))
}
save(notmakers2, marker2, total_pearson, file="makeOutlier.Rdata")

colnames(total_pearson) <- c("nonM-nonM", "M-nonM", "M-M", "level")

# total_pearson2 <- total_pearson[,c(3,6,9)]
# colnames(total_pearson2) <- c("Ctrl", "Outlier", "Outlier2")
library(reshape2)
total_pearson3 <- melt(total_pearson, id.vars = "level")
total_pearson3$level <- paste("outlier=", total_pearson3$level, sep="")
total_pearson3$level <- factor(total_pearson3$level, levels=paste("outlier=", c(0,1,3,5,7,9,11,13,15), sep=""))
colnames(total_pearson3) <- c("level", "group", "value")
ggplot(total_pearson3, aes(x=group, y=value)) + 
  geom_boxplot(aes(fill=group)) +
  facet_wrap(~ level,ncol=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_blank(),
    axis.text.x  = element_text(face="plain", angle=30, size = 12, color = "black", vjust=0.5),
    axis.text.y  = element_text(face="plain", size = 12, color = "black")) +
  labs(x = "Paired group", y = "Correlation\n", title = "Effect of outlier on correlation") +
  scale_fill_manual(values=brewer.pal(3,"Set3"))



