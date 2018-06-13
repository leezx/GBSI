setwd("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/")
# source('R/lib.R') # Nature-2018 R lib
source("/Users/surgery/Project/HOME/myScript/zxli_lib.R") # my R lib
options(stringsAsFactors = F)
options(mc.cores = 3)
library(SingleCellExperiment)
library(SC3)
library(scater)
library(pheatmap)
library(RColorBrewer)
set.seed(1234567)

# library(parallel)
# Calculate the number of cores
# no_cores <- detectCores() - 1
# Initiate cluster
# cl <- makeCluster(no_cores)

load("eachGroup/HSCR_5c3.Rdata")
expr_log <- t(logcounts(tmp_group))
# expr_log2 <- data.frame(empty=rep(0,dim(expr_log)[1]), full=rep(1,dim(expr_log)[1]),expr_log)
#expr_log <- as.data.frame(t(logcounts(tmp_group)))
## check one gene
# hist(expr_log["CLDN4",], breaks = 50)
# plot(base::sort(expr_log["CLDN4",], decreasing=T))

# expr_log[expr_log>2] <- 1
q <- length(rownames(expr_log)[expr_log[,"CLDN4"]>1 & expr_log[,"CLDN7"]>1])
m <- sum(expr_log[,"CLDN4"]>1)
n <- dim(expr_log)[1]
k <- sum(expr_log[,"CLDN4"]>1)
# https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper
phyper(q, m, n, k, lower.tail=F)

my_sum <- function(a){
  return(list(a=sum(a),b=sum(a)+1))
}
#expr_log3 <- expr_log2[,colSums(expr_log2>1)>0]
#expr_log2 <- expr_log3
expr_log2[expr_log2>=1] <- 1
expr_log2[expr_log2<1] <- 0
overlap_1 <- as.matrix(expr_log2) %*% as.matrix(t(expr_log2))
overlap_2 <- as.matrix(t(expr_log2)) %*% as.matrix(expr_log2) 
use.gene2 <- apply(overlap_2, 1, max) >= 5
# rbind(expr_log2,c("empty",rep(0,dim(expr_log2)[2])))
test_log <- expr_log2[1:5,1:8]

get_pVal <- function(a,arg1){
  expr_log <- arg1
  apply(expr_log, 2, 
        function(x) phyper(sum(a>1 & x>1), sum(a>1), length(a), sum(a>1), lower.tail=F))}
get_cor <- function(a,arg1){
  expr_log <- arg1
  apply(expr_log, 2, function(x) WGCNA::cor(a>1, x>1))
}

## to get a p-value matrix
# apply(expr_log, 1, my_sum)
cor_matrix <- WGCNA::cor(as.matrix(expr_log2), method = "pearson")
pVal_matrix <- do.call(cbind, lapply(expr_log, get_cor, arg1=expr_log))
pVal_matrix <- do.call(cbind, parallel::parLapply(cl, expr_log, get_pVal, arg1=expr_log))
# do.call(cbind, lapply(iris[, 1:4], my_sum))

# to get distance
library(philentropy)
dis_matrix <- stats::dist(t(expr_log2), method = "euclidean")

# test which distance is the best!!!
expr_log3 <- scale(expr_log)
plot_markers <- c("SEMA3D","SMAGP","KRT7","CLDN4","CLDN7","CDH1","EFNA1","ACSM3","WNT6",
                  "PTPRZ1","MOXD1","SCRG1","EDNRB","INSC","SPP1","PLP1",
                  "MSX1","NELL2","PLCH1","MMRN1","TRPM3","WNT2B","COLEC12","LIX1","ZIC1","RBFOX1","OTX1","TMEM88","ZIC2",
                  "ISL1","NEUROD1","EYA2","SCG3","SIX1","SYT4","EBF1","NHLH1","ZMAT4","STMN2","ENO3","TAGLN3",
                  "FAP","PRRX1","TGFBI","TWIST1","LUM","CDH11")

library(parallelDist)
# dis_matrix <- parDist(x = as.matrix(t(test_log)), method = "binary")
dis_matrix <- parDist(x = as.matrix(t(expr_log2)), method = "binary", threads=3)
dis_matrix <- as.matrix(dis_matrix)
dis_matrix[is.na(dis_matrix)] <- 0
rownames(dis_matrix) <- colnames(expr_log2)
colnames(dis_matrix) <- rownames(dis_matrix)
saveRDS(dis_matrix, file = "dis_matrix_5c3.rds")
# dis_matrix <- readRDS("dis_matrix_5c3.rds")
empty_gene <- colnames(dis_matrix)[dis_matrix["empty",] < 0.5] # remove genes 
full_gene <- colnames(dis_matrix)[dis_matrix["full",] < 0.3]
use.gene <- rownames(dis_matrix)[!rownames(dis_matrix)%in%c(empty_gene,full_gene)]
dis_matrix2 <- dis_matrix[use.gene,use.gene]
use.gene3 <- rownames(dis_matrix2)[rownames(dis_matrix2)%in%names(use.gene2)[use.gene2==T]]
dis_matrix3 <- dis_matrix2[use.gene3,use.gene3]
#for (i in 1:dim(dis_matrix3)[1]){
#  dis_matrix3[i,i] <- 1 }
#unuse.gene <- apply(dis_matrix3, 1, min) > 0.8
#dis_matrix2 <- dis_matrix2[!unuse.gene, !unuse.gene]
# [plot_markers,plot_markers]
use.gene4 <- rowSums(dis_matrix3<0.7) >= 3
dis_matrix3 <- dis_matrix3[use.gene4, use.gene4]
# stopCluster(cl)

# K-medoids Distance-Based clustering
library(kmed)
result <- fastkmed(dis_matrix3, ncluster = 10, iterate = 50)
# genes <- rownames(dis_matrix[dis_matrix["empty",] < 0.5,])
# genes <- rownames(dis_matrix3[apply(dis_matrix3, 1, min) > 0.8,])

genes <- names(result$cluster)[result$cluster==3]
genes <- genes[genes%in%rownames(tmp_group)]
# pheatmap(logcounts(tmp_group)[genes,])
# hist(rowSums(logcounts(tmp_group)[genes,] > 0), breaks = 10)
b <- apply(overlap_2[genes,genes], 1, max) >= 5
pheatmap(logcounts(tmp_group)[genes[b],], show_colnames = F)

result2 <- fastkmed(dis_matrix2[genes[b],genes[b]], ncluster = 2, iterate = 50)
genes2 <- names(result2$cluster)[result2$cluster==2]
pheatmap(logcounts(tmp_group)[genes2,], show_colnames = F)

## Hierarchical Clustering
# http://stat.ethz.ch/R-manual/R-patched/library/stats/html/hclust.html
dis_matrix4 <- as.dist(dis_matrix3)
hc <- hclust(dis_matrix4, method="average")
plot(hc)
groups<-cutree(hc, h=0.7)
marker_result <- data.frame(gene=names(groups),cluster=as.vector(groups),row.names = names(groups))
filter_cluster <- names(table(marker_result$cluster))[table(marker_result$cluster)>=3]
marker_result <- marker_result[marker_result$cluster%in%filter_cluster,]
plot_markers <- marker_result[marker_result$cluster==4,]
plot_markers <- plot_markers[plot_markers$gene%in%rownames(tmp_group),]
pheatmap(logcounts(tmp_group)[plot_markers$gene,])

## MLE get parameter of a distritution
a <- data.frame(gene=genes ,count=rowSums(logcounts(tmp_group)[genes,] >= 1), row.names = genes)
# m1 <- zeroinfl(count~count, data = a, dist = "negbin", EM = TRUE)
# lambda_para <- (((mean(a$count))^2+var(a$count))/mean(a$count)) - 1
# pi_para <- (var(a$count)-mean(a$count))/(var(a$count)+(mean(a$count))^2-mean(a$count))
library(stats4)
library(VGAM)
NLLzinb = function(size, prob = NULL, munb = NULL, pstr0 = 0, log = FALSE) {
-sum(VGAM::dzinegbin(a$count, size, prob, munb, pstr0)) }
m = mean(a$count)
mzinbinom = mle(minuslogl=NLLzinb,start=list(size=0.5, munb=m, pstr0=0.1),lower=0.002,method="L-BFGS-B")
mzinbinom

dzinbinom = function(x,mu,size,zprob,log=FALSE) {
  v = ifelse(x==0,
             zprob+(1-zprob)*dnbinom(0,mu=mu,size=size),
             (1-zprob)*dnbinom(x,mu=mu,size=size))
  if (log) return(log(v)) else v}
NLLzinb = function(mu,k,zprob) {
  -sum(dzinbinom(a$count,mu=mu,size=k,zprob=zprob,log=TRUE))}
m = mean(a$count)
mzinbinom = mle(minuslogl=NLLzinb,start=list(mu=m,k=0.5,zprob=0.5), method = "L-BFGS-B", lower = 0.002)
mzinbinom




