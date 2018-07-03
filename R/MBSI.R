load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_20c7.Rdata")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/IMR90_ENCC.Rdata")
tmp_group <- readRDS("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eight_groups_2702/eight_groups_2702_cells.rds")
# expr_log <- t(logcounts(tmp_group))
expr_log <- t(logcounts(tmp_group[rowData(tmp_group)[,"gene_type"] == "protein_coding",]))
# expr_log2 <- expr_log[,plot_markers]
expr_log2 <- expr_log[,colSums(expr_log > 0) >= 3] # remove gene express less than 3 cells
expr_log2 <- expr_log2[,colSums(expr_log2 > 0) <= 0.9*dim(expr_log2)[1]] # remove gene express more than 0.9*total cells
# expr_log2 <- t(expr_log2)
expr_log3 <- scale(expr_log2)
# expr_log3[is.na(expr_log3)]
# judge if a gene express in a cell
expr_log4 <- apply(expr_log2>=1,2,function(x) {storage.mode(x) <- 'integer'; x})
overlap <- t(expr_log4) %*% expr_log4 

# simulate data
expr_log3 <- t(simulate_data())
# expr_log3 <- t(simulate_data(10,100,50,show_name = F))
# add outlier
expr_log3 <- rbind(expr_log3, full=rep(c(19,20),30), empty=rep(c(0,1),30))
pheatmap(t(expr_log3), cluster_rows = F, cluster_cols = F)

library(parallelDist)
dis_matrix <- parDist(x = as.matrix(t(expr_log3)), method = "euclidean", threads=3)
dis_matrix <- as.matrix(dis_matrix)
dis_matrix[is.na(dis_matrix)] <- 0
rownames(dis_matrix) <- colnames(expr_log3)
colnames(dis_matrix) <- rownames(dis_matrix)

pheatmap(dis_matrix, cluster_rows = F, cluster_cols = F)
# saveRDS(dis_matrix, file = "euclidean_dis_matrix_5c3.rds")
# saveRDS(dis_matrix, file = "euclidean_dis_matrix_20c7.rds")
# saveRDS(dis_matrix, file = "euclidean_dis_matrix_eight.rds")
# saveRDS(dis_matrix, file = "cells_euclidean_dis_matrix_eight.rds")
# dis_matrix <- apply(dis_matrix,2,function(x) {storage.mode(x) <- 'integer'; x})

library(WGCNA)
# standard deviation can't be zero
cor_matrix_pearson <- WGCNA::cor(x = as.matrix((expr_log3)), method = "pearson")
cor_matrix_spearman <- WGCNA::cor(x = as.matrix((expr_log3)), method = "spearman")
# cor_matrix_kendall <- WGCNA::cor(x = as.matrix((expr_log3)), method = "kendall") # too slow
pheatmap(cor_matrix_spearman, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T)

a <- c("SEMA3D","SMAGP","KRT7","CLDN4","CLDN7","CDH1","EFNA1","ACSM3","WNT6")
b <- c("PTPRZ1","MOXD1","SCRG1","EDNRB","INSC","SPP1","PLP1")
c <- c("MSX1","NELL2","PLCH1","MMRN1","TRPM3","WNT2B","COLEC12","LIX1","ZIC1","RBFOX1","OTX1","TMEM88","ZIC2")
d <- c("ISL1","NEUROD1","EYA2","SCG3","SIX1","SYT4","EBF1","NHLH1","ZMAT4","STMN2","ENO3","TAGLN3")
e <- c("FAP","PRRX1","TGFBI","TWIST1","LUM","CDH11")

simulate_df <- simulate_data(10,100,50,show_name = F)
res <- dbscan(simulate_df, eps = .5, minPts = 3) # col is feature, row is sample

# get gene cor matrix
# col is sample/cell, row is gene
get_cor_matrix_of_gene <- function(M, method="pearson", show_matrix=F, show_name=F){
  library(WGCNA)
  cor_matrix <- WGCNA::cor(x = as.matrix(t(M)), method = method)
  if (show_matrix){
    pheatmap(cor_matrix, cluster_rows = F, cluster_cols = F, show_rownames = show_name, show_colnames = show_name)
  }
  return(cor_matrix)
}

# get gene hub from a gene list
get_gene_hub_from_list <- function(dis_matrix, genelist){
  dis_matrix <- dis_matrix[,genelist]
  mean_distance <- apply(dis_matrix, 1, mean) # for each coordinate, get the mean
  distance1 <- apply(dis_matrix, 2, function(y){as.vector(y)-as.vector(mean_distance)}) # for each coordinate, get the distance to the mean
  distance2 <- apply(distance1,2,function(x){sum(x^2)}) # get square
  center_gene <- sort(distance2, decreasing=F)[1] # order
  return(names(center_gene))
}

# get gene hub of in a matrix
get_gene_hub_from_dismatrix <- function(M, cluster){
  dis_matrix <- t(M) # cell is the coordinate
  genehub <- c()
  for (i in names(table(cluster))) {
    genelist <- names(cluster[cluster==as.integer(i)])
    center_gene <- get_gene_hub_from_list(dis_matrix, genelist)
    genehub <- c(center_gene, genehub)
  }
  return(genehub)
}

# get gene hub by DBSCAN
get_gene_hub_by_DBSCAN <- function(M, eps=0.5, minPts=3, minGeneNum=3){
  library(dbscan)
  res <- dbscan(M, eps = eps, minPts = minPts) # col is feature, row is sample
  cluster <- res$cluster
  names(cluster) <- rownames(M)
  valid_cluster <- names(table(cluster))[table(cluster) >= minGeneNum]
  cluster <- cluster[cluster%in%valid_cluster]
  cluster <- cluster[cluster!=0]
  genehub <- get_gene_hub_from_dismatrix(M, cluster)
  return(genehub)
}

# get most correlated gene
get_most_correlated_gene <- function(gene, cor_matrix, store_gene=c(), pos=T) {
  print(paste("get nearest gene:", gene, sep = " "))
  store_gene <- unique(c(gene, store_gene))
  nearest_gene <- sort(cor_matrix[gene,!(colnames(cor_matrix) %in% store_gene)], decreasing = pos)[1]
  return(names(nearest_gene))
}

# get full gene group under a threshould
get_full_gene_group <- function(genehub, cor_matrix, cor_threshold = 0.5, minOverlap = 3) {
  
}

# get nearest point
get_nearest_point <- function(gene="SEMA3D", gene_list=c()){
  print(paste("get nearest gene:", gene, sep = " "))
  gene_list <- unique(c(gene, gene_list))
  nearest_gene <- sort(dis_matrix[gene,!(colnames(dis_matrix) %in% gene_list)], decreasing = F)[1]
  if (nearest_gene==0) {
    print(c("ERROR: ", gene, names(nearest_gene), "distance == 0"))
    quit(status=1)
  }else{
    return(names(nearest_gene))
  }}

# get center point
get_center_point <- function(gene_list=c("SEMA3D","SMAGP","KRT7","CLDN4")){
  if (length(gene_list)<=2){return(gene_list[1])}
  dis_matrix_local <- expr_log3[,gene_list]
  mean_distance <- apply(dis_matrix_local, 1, mean)
  distance1 <- apply(dis_matrix_local, 2, function(y){as.vector(y)-as.vector(mean_distance)})
  distance2 <- apply(distance1,2,function(x){sum(x^2)})
  center_gene <- sort(distance2, decreasing=F)[1]
  return(names(center_gene))
}

# get overlap cells
get_overlap_cells <- function(expr_log4, gene_list=c("SEMA3D","SMAGP","KRT7","CLDN4")){
  if (length(gene_list)==1){
    all_cells <- expr_log4[,gene_list]
    return(names(all_cells[all_cells==T]))
  } else {
    all_cells <- rowSums(expr_log4[,gene_list]) >= 0.8*length(gene_list)
    return(names(all_cells[all_cells==T]))
  }}

# overlap: 5 cells, 80% overall cells?
# distance: depends on overlap number, abs value, max value 110%
judge_condition_main <- function(gene_list=c("SEMA3D","SMAGP","KRT7","CLDN4"), nearest_point="CLDN7"){
  all_cells <- expr_log4[,nearest_point] == 1
  current_cells <- names(all_cells[all_cells==T])
  if (length(intersect(current_cells, get_overlap_cells(gene_list)))/length(get_overlap_cells(gene_list)) < 0.8){return(F)}
  else if (length(intersect(current_cells, get_overlap_cells(gene_list)))<10) {return(F)}
  else if (dis_matrix[get_center_point(gene_list),nearest_point]>15) {return(F)}
  else {return(T)}
}

judge_condition_rare <- function(gene_list=c("SEMA3D","SMAGP","KRT7","CLDN4"), nearest_point="CLDN7"){
  all_cells <- expr_log4[,nearest_point] == 1
  current_cells <- names(all_cells[all_cells==T])
  if (length(intersect(current_cells, get_overlap_cells(gene_list)))/length(get_overlap_cells(gene_list)) < 0.5){return(F)}
  if (length(intersect(current_cells, get_overlap_cells(gene_list)))<5) {return(F)}
  else if (dis_matrix[get_center_point(gene_list),nearest_point]>20) {return(F)}
  else {return(T)} }

library(kmed)
result <- fastkmed(dis_matrix, ncluster = 100, iterate = 50)
clusters <- names(table(result$cluster))
# clusters <- as.integer(clusters)
clusters_genes <- list()
for (i in clusters){
  i <- as.integer(i)
  # print(clusters[i])
  clusters_genes[[i]] <-  names(result$cluster[result$cluster==clusters[i]])
}

library(ClusterR)
# result2 <- kmeans(dis_matrix, centers = 30)
result2 <- KMeans_arma(dis_matrix, clusters = 30)
pr = predict_KMeans(dis_matrix, result2)
names(pr) <- rownames(dis_matrix)
clusters <- names(table(pr))
clusters_genes <- list()
for (i in clusters){
  i <- as.integer(i)
  # print(clusters[i])
  clusters_genes[[i]] <-  names(pr[pr==clusters[i]])
}

# not suggest
hc <- hclust(as.dist(dis_matrix), method="average")
# groups<-cutree(hc, h=0.7)
groups <- cutree(hc, k=30)
clusters <- names(table(groups))
clusters_genes <- list()
for (i in clusters){
  i <- as.integer(i)
  # print(clusters[i])
  clusters_genes[[i]] <-  names(groups[groups==clusters[i]])
}

start_gene_list <- c()
for (i in clusters){
  # print(clusters[i])
  i <- as.integer(i)
  start_gene_list <- c(start_gene_list, get_center_point(clusters_genes[[i]]))
}

get_common_markers <- function(start_gene="PLP1", exist_list=c()){
  #start_gene <- start_gene_list[14]
  gene_list <- unique(c(start_gene))
  exist_list <- unique(c(gene_list, exist_list))
  while (T) {
    nearest_point <- get_nearest_point(start_gene, exist_list)
    if (judge_condition_main(gene_list, nearest_point)){
      gene_list <- c(gene_list, nearest_point)
      exist_list <- unique(c(gene_list, exist_list))
    } else {break}
    if (length(gene_list)>= 20) {break}
    start_gene <- get_center_point(gene_list) }
  return(gene_list) 
  }

common_marker_group <- list()
exist_list <- c()
k <- 1
for (i in clusters){
  # print(clusters[i])
  i <- as.integer(i)
  common_markers <- get_common_markers(start_gene_list[i], exist_list)
  if (length(common_markers)<3){next}
  else { 
    common_marker_group[[k]] <-  common_markers 
    exist_list <- unique(c(exist_list, common_markers))
    k <- k+1
    }}

pheatmap(logcounts(tmp_group)[common_marker_group[[1]],], show_colnames = F)

simulate_data <- function(logExpr=10, ngene=10, ncell=10, show_matrix=T, show_name=T){
  ideal_matrix <- as.data.frame(matrix(rep(logExpr,ngene*ncell), nrow=ngene, ncol=ncell))
  rownames(ideal_matrix) <- ((ngene*0+1):(ngene*1))
  colnames(ideal_matrix) <- ((ncell*0+1):(ncell*1))
  ideal_matrix$name <- rownames(ideal_matrix)
  level_matrix <- as.data.frame(matrix(rep(logExpr,ngene*ncell), nrow=ngene, ncol=ncell)/sqrt(seq(1:ngene)))
  rownames(level_matrix) <- ((ngene*1+1):(ngene*2))
  colnames(level_matrix) <- ((ncell*1+1):(ncell*2))
  level_matrix$name <- rownames(level_matrix)
  level_matrix2 <- as.data.frame(t(t(as.data.frame(matrix(rep(logExpr,ngene*ncell), nrow=ngene, ncol=ncell)))/sqrt(seq(1:ncell))))
  rownames(level_matrix2) <- ((ngene*2+1):(ngene*3))
  colnames(level_matrix2) <- ((ncell*2+1):(ncell*3))
  level_matrix2$name <- rownames(level_matrix2)
  noise_matrix <- as.data.frame(matrix(rep(logExpr,ngene*ncell), nrow=ngene, ncol=ncell))
  rownames(noise_matrix) <- ((ngene*3+1):(ngene*4))
  colnames(noise_matrix) <- ((ncell*3+1):(ncell*4))
  noise_matrix$name <- rownames(noise_matrix)
  ideal_matrix2 <- as.data.frame(matrix(rep(0,ngene*ncell), nrow=ngene, ncol=ncell))
  rownames(ideal_matrix2) <- ((ngene*4+1):(ngene*5))
  colnames(ideal_matrix2) <- ((ncell*4+1):(ncell*5))
  ideal_matrix2$name <- rownames(ideal_matrix2)
  ideal_matrix3 <- as.data.frame(matrix(rep(logExpr,ngene*ncell), nrow=ngene, ncol=ncell))
  rownames(ideal_matrix3) <- ((ngene*5+1):(ngene*6))
  colnames(ideal_matrix3) <- ((ncell*5+1):(ncell*6))
  ideal_matrix3$name <- rownames(ideal_matrix3)
  total_matrix <- merge(merge(ideal_matrix, level_matrix,  by="name", all = T), level_matrix2, by="name", all=T)
  total_matrix <- merge(total_matrix, noise_matrix, by="name", all=T)
  total_matrix <- merge(total_matrix, ideal_matrix2, by="name", all=T)
  total_matrix <- merge(total_matrix, ideal_matrix3, by="name", all=T)
  total_matrix[is.na(total_matrix)] <- 0
  total_matrix <- total_matrix[order(as.numeric(total_matrix$name)),]
  rownames(total_matrix) <- paste("gene",total_matrix$name, sep="_")
  colnames(total_matrix) <- paste("cell",colnames(total_matrix), sep="_")
  total_matrix2 <- total_matrix
  total_matrix2$cell_name <- NULL
  # count <- ncell
  for (onegene in ((ngene*3+1):(ngene*4))){
    cell.no <- as.integer(runif(1, 0, ncell*3/3))
    #inside.no <- as.integer(runif(1, 0, cell.no))
    #outcell <- sample(((ncell*0+1):(ncell*2)), ncell-count)
    outcell <- sample(1:(6*ncell), cell.no)
    #incell <- sample(((ncell*2+1):(ncell*3)), ncell-count)
    #total_matrix2[onegene, incell] <- runif(1, 0, logExpr)
    total_matrix2[onegene, outcell] <- runif(1, 0, logExpr)
    #count <- count - 1
    #if (count <= 0) {count <- ncell}
  }
  total_matrix2[((ngene*4+1):(ngene*5)),(ncell*4-as.integer((ncell)/2)):((ncell*4-as.integer((ncell)/2))+ncell-1)] <- logExpr
  total_matrix2[((ngene*5+1):(ngene*6)),(ncell*4+1):(ncell*6)] <- logExpr
  if (show_matrix){
    library(pheatmap)
    pheatmap(total_matrix2, cluster_rows = F, cluster_cols = F, show_rownames = show_name, show_colnames = show_name)
  }
  return(total_matrix2)
}

get_distance_quantile <- function(M, unit=0.05, percent=1){
  M_half <- M[row(M) > col(M)]
  return(quantile(M_half, seq(0, 1, unit))[percent])
}

res <- dbscan(x, eps = .1, minPts = 3) # col is feature, row is sample


