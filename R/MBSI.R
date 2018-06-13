load("eachGroup/HSCR_5c3.Rdata")
load("eachGroup/HSCR_20c7.Rdata")
load("eachGroup/IMR90_ENCC.Rdata")
tmp_group <- readRDS("eight_groups_2702/eight_groups_2702_cells.rds")
expr_log <- t(logcounts(tmp_group))
expr_log <- t(logcounts(tmp_group[rowData(tmp_group)[,"gene_type"] == "protein_coding",]))
# expr_log2 <- expr_log[,plot_markers]
expr_log2 <- expr_log[,colSums(expr_log > 0) >= 3]
expr_log2 <- expr_log2[,colSums(expr_log2 > 0) <= 0.9*dim(expr_log2)[1]]
# expr_log2 <- t(expr_log2)
expr_log3 <- scale(expr_log2)
# expr_log3[is.na(expr_log3)]
# judge if a gene express in a cell
expr_log4 <- apply(expr_log2>=1,2,function(x) {storage.mode(x) <- 'integer'; x})
overlap <- t(expr_log4) %*% expr_log4 

library(parallelDist)
dis_matrix <- parDist(x = as.matrix(t(expr_log3)), method = "euclidean", threads=3)
dis_matrix <- as.matrix(dis_matrix)
dis_matrix[is.na(dis_matrix)] <- 0
rownames(dis_matrix) <- colnames(expr_log3)
colnames(dis_matrix) <- rownames(dis_matrix)
# saveRDS(dis_matrix, file = "euclidean_dis_matrix_5c3.rds")
# saveRDS(dis_matrix, file = "euclidean_dis_matrix_20c7.rds")
# saveRDS(dis_matrix, file = "euclidean_dis_matrix_eight.rds")
# saveRDS(dis_matrix, file = "cells_euclidean_dis_matrix_eight.rds")
# dis_matrix <- apply(dis_matrix,2,function(x) {storage.mode(x) <- 'integer'; x})

a <- c("SEMA3D","SMAGP","KRT7","CLDN4","CLDN7","CDH1","EFNA1","ACSM3","WNT6")
b <- c("PTPRZ1","MOXD1","SCRG1","EDNRB","INSC","SPP1","PLP1")
c <- c("MSX1","NELL2","PLCH1","MMRN1","TRPM3","WNT2B","COLEC12","LIX1","ZIC1","RBFOX1","OTX1","TMEM88","ZIC2")
d <- c("ISL1","NEUROD1","EYA2","SCG3","SIX1","SYT4","EBF1","NHLH1","ZMAT4","STMN2","ENO3","TAGLN3")
e <- c("FAP","PRRX1","TGFBI","TWIST1","LUM","CDH11")

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
get_overlap_cells <- function(gene_list=c("SEMA3D","SMAGP","KRT7","CLDN4")){
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
