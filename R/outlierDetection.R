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

cellCountPerGene <- rowSums(Test_3_Pollen$in_X>0)
geneCountPerCell <- colSums(Test_3_Pollen$in_X>0)
totalReadCount <- colSums(Test_3_Pollen$in_X)
cvPerGene <- apply(Test_3_Pollen$in_X, 1, sd)/apply(Test_3_Pollen$in_X, 1, mean)

expr <- logcounts(tmp_group)
expr <- Test_3_Pollen$in_X
expr <- logcounts(deng)
expr <- expr[rowSums(expr>0)>3,]
# expr <- expr[rowSums(expr)>0,]

# get markers
all_markers_raw <- get_marker_genes(logcounts(deng), deng$cell_type1)
all_markers <- all_markers_raw[all_markers_raw$auroc>0.95 & !is.na(all_markers_raw$auroc) & all_markers_raw$pvalue < 0.001,]
sc3_marker <- data.frame(name=rownames(deng)[as.numeric(rownames(all_markers))], cluster=all_markers$clusts)
rownames(sc3_marker) <- sc3_marker$name
annotation_col <- data.frame(Cluster = factor(deng$cell_type1), row.names = colnames(deng))
pheatmap(logcounts(deng)[sc3_marker[sc3_marker$cluster=="2cell",]$name, order(deng$cell_type1)], cluster_rows = T, show_rownames = F, show_colnames = F, cluster_cols = F, annotation_col = annotation_col)
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

corMatrix <- buildGeneNetworkCor(expr=expr, method = "pearson")
# library("Matrix")
# corMatrix <- as(as.matrix(corMatrix), "dgCMatrix")
corMatrix2 <- buildGeneNetworkCor(expr=expr[,rownames(pca)[pca$cluster=="1"]], , method = "pearson")

# my dataset
pca <- outlierDetection(M=logcounts(tmp_group), minPts = 5, percent = 0.01, pcNum = 2)
corMatrix <- buildGeneNetworkCor(expr=logcounts(tmp_group), method = "pearson")
corMatrix2 <- buildGeneNetworkCor(expr=logcounts(tmp_group)[,rownames(pca)[pca$cluster=="1"]], , method = "pearson")
cor4 <- as.vector(corMatrix[a,a])
cor5 <- as.vector(corMatrix2[a,a])

# check deng known markers
dengMarkers <- read.csv("/Users/surgery/Project/HOME/github/MBSIT/SC3/nmeth.4236-S3.csv", header = T)
rownames(dengMarkers) <- dengMarkers$Gene
dengMarkers.sig <- dengMarkers[dengMarkers$AUC>0.9 & dengMarkers$p.value<0.01,]
dengMarkers.sig1 <- dengMarkers.sig[dengMarkers.sig$clusts=="7",]
# before
non.non.markersCor <- corMatrix[!rownames(corMatrix)%in%dengMarkers$Gene, !colnames(corMatrix)%in%dengMarkers$Gene]
non.yes.markersCor <- corMatrix[rownames(corMatrix)%in%rownames(b), !colnames(corMatrix)%in%dengMarkers$Gene]
yes.yes.markersCor <- corMatrix[rownames(corMatrix)%in%rownames(b), colnames(corMatrix)%in%rownames(b)]
# hist before
cor1 <- (sample(abs(non.non.markersCor), 1000))
cor2 <- (sample(abs(non.yes.markersCor), 1000))
cor3 <- (sample(abs(yes.yes.markersCor), 1000))
rm(non.non.markersCor, non.yes.markersCor, yes.yes.markersCor)
#
# after
# non-markers
non.non.markersCor2 <- corMatrix2[!rownames(corMatrix2)%in%dengMarkers$Gene, !colnames(corMatrix2)%in%dengMarkers$Gene]
non.yes.markersCor2 <- corMatrix2[rownames(corMatrix2)%in%rownames(b), !colnames(corMatrix2)%in%dengMarkers$Gene]
yes.yes.markersCor2 <- corMatrix2[rownames(corMatrix2)%in%rownames(b), colnames(corMatrix2)%in%rownames(b)]
# hist before
cor4 <- (sample(abs(non.non.markersCor2), 1000))
cor5 <- (sample(abs(non.yes.markersCor2), 1000))
cor6 <- (sample(abs(yes.yes.markersCor2), 1000))
rm(non.non.markersCor2, non.yes.markersCor2, yes.yes.markersCor2)


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
outlierDetection <- function(M=expr2, pcNum=2, method="prcomp", minPts=4, percent=0.05, threads=1, plot=T){
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
outliers <- mvoutlier::pcout(pcMatrix, makeplot = FALSE,
                             explvar = 0.5, crit.M1 = 0.9,
                             crit.c1 = 5, crit.M2 = 0.9,
                             crit.c2 = 0.99, cs = 0.25,
                             outbound = 0.05)
!as.logical(outliers$wfinal01)