library('Matrix')
library('parallel')
library('MASS')
library('diffusionMap')
library('FNN')
library('igraph')
library('princurve')
library('ggplot2')
library('inline')
library('gplots')

# for maturation trajectory
# fit maturation trajectory
maturation.trajectory <- function(md, expr, cm=c(), pricu.f=1/3) {
  cat('Fitting maturation trajectory\n')
  # genes <- apply(cm[rownames(expr), ] > 0, 1, mean) >= 0.02 & apply(cm[rownames(expr), ] > 0, 1, sum) >= 3
  expr <- logcounts(sce[rowData(sce)$gene_type=="protein_coding",sce$cellGroup %in% c("SAG4_ENCC","SAG10_ENCC","ENCC-derived_neurons")])
  rd <- dim.red(expr, max.dim=50, ev.red.th=0.04)
  # for a consisten look use Nes expression to orient each axis
  # for (i in 1:ncol(rd)) {
  #   if (cor(expr['Nes', ], rd[, i]) > 0) {
  #     rd[, i] <- -rd[, i]
  #   }
  # }
  
  # md <- md[, !grepl('^DMC', colnames(md))]
  # md <- cbind(md, rd)
  library(princurve)
  load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/1-batch/1-formal/SAG/SAG_PCA.Rdata")
  md <- SAG_PCA$data
  #md2 <- md[!(md$PC1 > 10 & md$PC2 < -4),]
  md1 <- (md[md$shape_by %in% c("Ctrl_ENCC","ENCC-derived_neurons"),])
  #md2 <- (md[md$shape_by %in% c("SAG4_ENCC","SAG10_ENCC","ENCC-derived_neurons"),])
  md2 <- (md[md$shape_by %in% c("SAG4_ENCC","SAG10_ENCC"),])
  #md2 <- md2[!(md2$PC1 > -5 & md2$PC2 < 2),]
  #md2 <- md2[!(md2$PC1 < -5 & md2$PC2 > 0),]
  rd1 <- as.matrix(md1[,1:2])
  # md2 <- md2[!rownames(md2)%in%sample(rownames(md2[md2$colour_by == "Cluster6",]), 450),]
  rd2 <- as.matrix(md2[,1:2])
  
  pricu1 <- principal.curve(rd1, smoother='lowess', trace=TRUE, f=0.6, stretch=333)
  # two DMCs
  pc.line1 <- as.data.frame(pricu1$s[order(pricu1$lambda), ])
  # lambda, for each point, its arc-length from the beginning of the curve. The curve is parametrized approximately by arc-length, and hence is unit-speed.
  maturation.score1 <- pricu1$lambda/max(pricu1$lambda)
  
  maturation.score <- my_project_to_curve(M=rd2, pc.line1=pc.line1, maturation.score1=maturation.score1)
  # start2 <- order(md2$colour_by)
  # md22 <- md2[md2$colour_by == "Cluster1",]
  # md22 <- md22[(md22$PC1 < 9 & md22$PC1 > -11 & md22$PC2 < 0 & md22$PC2 > -4),]
  # rd22 <- rd2[rownames(md22),]
  # pricu2start <- principal.curve(rd22, smoother='lowess', trace=TRUE, f=pricu.f, stretch=333)
  # pc.line22 <- as.data.frame(pricu2start$s[order(pricu2start$lambda), ])
  # start=rd22,
  # pricu2 <- principal.curve(rd2, smoother='lowess', trace=TRUE, f=0.6, stretch=333)
  # # two DMCs
  # pc.line2 <- as.data.frame(pricu2$s[order(pricu2$lambda), ])
  # # lambda, for each point, its arc-length from the beginning of the curve. The curve is parametrized approximately by arc-length, and hence is unit-speed.
  # maturation.score2 <- pricu2$lambda/max(pricu2$lambda)
  #
  # test <- project_to_curve(rd2, pricu1)
  # md3 <- md[!(md$PC1 > -5 & md$PC2 < 2),]
  # rd3 <- as.matrix(md3[,1:2])
  # pricu3 <- principal.curve(rd3, smoother='lowess', trace=TRUE, f=0.6, stretch=333)
  # # two DMCs
  # pc.line3 <- as.data.frame(pricu3$s[order(pricu3$lambda), ])

  # merge two score
  md$maturation.score <- 0
  md[rownames(md1),]$maturation.score <- maturation.score1
  md[rownames(md2),]$maturation.score <- maturation.score
  # orient maturation score using Nes expression
  # if (cor(md$maturation.score, expr['Nes', ]) > 0) {
  #   md$maturation.score <- -(md$maturation.score - max(md$maturation.score))
  # }
  
  # use 1% of neighbor cells to smooth maturation score
  # md$maturation.score.smooth <- nn.smooth(md$maturation.score, rd[, 1:2], round(ncol(expr)*0.01, 0))
  
  # pick maturation score cutoff to separate mitotic from post-mitotic cells
  # md$in.cc.phase <- md$cc.phase != 0
  # fit <- loess(as.numeric(md$in.cc.phase) ~ md$maturation.score.smooth, span=0.5, degree=2)
  # md$cc.phase.fit <- fit$fitted
  # pick MT threshold based on drop in cc.phase cells
  # ignore edges of MT because of potential outliers
  # mt.th <- max(subset(md, cc.phase.fit > mean(md$in.cc.phase)/2 & maturation.score.smooth >= 0.2 & maturation.score.smooth <= 0.8)$maturation.score.smooth)
  
  # md$postmitotic <- md$maturation.score.smooth > mt.th
  return(list(md=md, pricu=pricu, pc.line=pc.line, mt.th=mt.th))
}

my_project_to_curve <- function(M=rd2, pc.line1=pc.line1, maturation.score1=maturation.score1){
  nearestPoint <- c()
  for (i in 1:dim(M)[1]) {
    x <- M[i,1]
    y <- M[i,2]
    disM <- c()
    for (j in 1:dim(pc.line1)[1]) {
      distance <- (pc.line1[j,1]-x)^2 + (pc.line1[j,2]-y)^2
      disM <- c(disM, distance)
    }
    nearestPoint <- c(nearestPoint,rownames(pc.line1)[order(disM, decreasing = F)][1])
  }
  maturation.score <- maturation.score1[nearestPoint]
  return(maturation.score)
}

# visualization
ggplot(md, aes(PC1, PC2)) + geom_point(aes(color=maturation.score), size=1, shape=16) + 
  scale_color_gradientn(colours=my.cols.RYG, name='Pseudotime') +
  # stat_density2d(n=111, na.rm=TRUE, color='black', size=0.33, alpha=0.5) +
  geom_line(data=pc.line1, color='deeppink', size=0.77) +
  # geom_line(data=pc.line2, color='blue', size=0.77) +
  # geom_line(data=pc.line3, color='blue', size=0.77) +
  # geom_line(data=pc.line22, color='blue', size=0.77) +
  theme_grey(base_size=12) + labs(x='PC1', y='PC2')

# for smoothing maturation score
nn.smooth <- function(y, coords, k) {
  knn.out <- FNN::get.knn(coords, k)
  w <- 1 / (knn.out$nn.dist+.Machine$double.eps)
  w <- w / apply(w, 1, sum)
  v <- apply(knn.out$nn.index, 2, function(i) y[i])
  return(apply(v*w, 1, sum))
}

# for clustering
# ev.red.th: relative eigenvalue cutoff of 2%
dim.red <- function(expr, max.dim, ev.red.th, plot.title=NA, do.scale.result=FALSE) {
  cat('Dimensionality reduction via diffusion maps using', nrow(expr), 'genes and', ncol(expr), 'cells\n')
  if (sum(is.na(expr)) > 0) {
    dmat <- 1 - cor(expr, use = 'pairwise.complete.obs')
  } else {
    # distance 0 <=> correlation 1
    dmat <- 1 - cor(expr)
  }
  
  max.dim <- min(max.dim, nrow(dmat)/2)
  dmap <- diffuse(dmat, neigen=max.dim, maxdim=max.dim)
  ev <- dmap$eigenvals
  # relative eigenvalue cutoff of 2%, something like PCA
  ev.red <- ev/sum(ev)
  evdim <- rev(which(ev.red > ev.red.th))[1]
  
  if (is.character(plot.title)) {
    # Eigenvalues, We observe a substantial eigenvalue drop-off after the initial components, demonstrating that the majority of the variance is captured in the first few dimensions.
    plot(ev, ylim=c(0, max(ev)), main = plot.title)
    abline(v=evdim + 0.5, col='blue')
  }
  
  evdim <- max(2, evdim, na.rm=TRUE)
  cat('Using', evdim, 'significant DM coordinates\n')
  
  colnames(dmap$X) <- paste0('DMC', 1:ncol(dmap$X))
  res <- dmap$X[, 1:evdim]
  if (do.scale.result) {
    res <- scale(dmap$X[, 1:evdim])
  } 
  return(res)
}