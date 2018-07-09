# MSIC

# module-based subtype identification and clustering of single-cell RNA-seq data 
title of paper: Clustering and trajectory inference of single-cell RNA-seq data by marker-based subgroup identification method

# Step1: kinds of marker group for clustering

## marker1: ideal situation
unique expressed genes with the same expression level.
simulate a 20x20 expression matrix.

## marker2: noise of distribution
some cells miss the markers, other cells express few markers.
simulate a series of pattern: missing, adding, combinatin.

## marker3: difference of expression level
markers in the marker group have different expression level, thus the distance is not zero, but we should make it to be zero.

## marker4: combination

## marker5: continuous markers which express in all cells

## marker6: intermediate markers

note: find all this theoretical markers in real (published) dataset will make your paper excellent.

## large scale simulation data
our tool can be used in a specific situation. marker is not unique expressed gene, but the identity of a subgroup.

# Algorithm
## 1. calculate the distance/correlation of each gene
## 2. find local center by local centrality clustering
## 3. iterate until get all the markers of a subgroup
## 4. clustering
## 5. apply in pseudotime

# transfer to statistic question

# test dataset
Mainstream platform:
smart-seq
10x

# key feature

1. No need to choose a K;
2. More intuitive and easy for downstream analysis;
3. Avoid the effect of condounder for clustering;
4. Intermediate subgroup;
5. Identify rare subgroup;
6. find negative marker;
7. Pseudotime analysis?

# methods comparison

1. different between DBSCAN and other methods;

2. use gene X cell matrix or gene correlation matrix for DBSCAN?

3. how to find the best eps? just like cross validation.

4. How can I choose eps and minPts for DBSCAN?

5. Why cor is better than dis?

6. how to compare pearson and spearman?

7. How to evaluate the significance of a module? MVN and variance analysis

8. remove outlier (infer what lead to the outlier)

     
# limits of current clustering methods

1. working with a Gene X Cell matrix, features and cluster k are the key determining factor. while the feature selection is often controversial and the choosing of cluster k  is often too subjective. Cell and cell similarity.
2. overfitting problem
3. features such as cell cycle genes and apoptosis genes are often confounding features which can largely affect the clustering result but have less biological meanings. whether to remove it or not is quite controversial. If not removed, any current clustering method can not avoid the confounds.
4. features often have 10k~20k, cells often have 100~5000. so current clustering method often need to  consume a lot of computing resources and time.
5. hard to identify rare subgroups.
6. Can give the negative markers of a subtype
7. the final purpose of clustering is to find markers, key role of marker. how many markers are enough to define a subtype? Marker is very necessary.




# outlier

1. why removing outlier is necessary? using simulated data to explain this point
2. why current  outlier detection method is not suitable? use ground truth to demonstrate this point, compare different tools (removing the rare populations)
3. why using dimension reduction (cite paper)? top 100 or 95% explained variance PCs
4. try different dimension reduction method (PCA and DM), compare the results
5. outlier detection algorithm, minNum

# Hypothesis

1. Different marker types (MT1-MTn), explain the biological meaning of this markers. using simulated data to demonstrate this point.

2. why using correlation? why Spearman is better than Pearson? why not distance? figure out the meaning of the formula. extremum and outlier will seriously impact the Pearson correlation. Pearson is linear.  (The Spearman correlation is less sensitive than the Pearson correlation to strong outliers that are in the tails of both samples. That is because Spearman's rho limits the outlier to the value of its rank.)

3. key modules locate in the local center of the gene network. using simulated data to demonstrate this point.

4. why using overlap?

5. how to choose the threshold? cutting threshold, percentage, top... why?

   

# nearest neighbors searching algorithm



# cluster inference

10 markers are enough to make the inference

you'd better make this method unsupervised

# full module identification

how to set threshold?

key set & full set



# summary analysis

1. overlap, remove genes have less than k (3) overlap cells with any other genes.
2. quantile, get the information of correlation, in order to set the threshold.
3. slices method for network analysis for local center identification.



# co-expression analysis

something can be done by the way



# Visualisation of the results

Scatter plot

Heatmap



1. Low expression of most genes;
2. Contamination;
3. less than 3 cells far away from the main population;
4. How to prove my method is better for scRNA-seq data?
5. the factor cause the outlier
6. the easiest way: PCA and diffusion map (nature paper)



# Detecting High Dimensional Outliers    

Interpretability of outliers 

Which subspaces manifest the outliers or an assessment regarding the “outlying-ness” of the objects 

Data sparsity: data in high-D spaces are often sparse 

The distance between objects becomes heavily dominated by noise as the dimensionality increases 

Data subspaces 

Local behavior and patterns of data 

Scalability with respect to dimensionality 

The number of subspaces increases exponentially    



# build R bioconductor package

# Backgroud

# Reference
[SIMLR: Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning](https://www.nature.com/articles/nmeth.4207)

[SC3: consensus clustering of single-cell RNA-seq data](https://www.nature.com/articles/nmeth.4236)

[GiniClust: detecting rare cell types from single-cell gene expression data with Gini index](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1010-4)

[Classification of low quality cells from single-cell RNA-seq data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/)

[基于局部中心量度的聚类算法研究](http://blog.sciencenet.cn/blog-3273400-1097494.html)

[一种基于K-Means局部最优性的高效聚类算法](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.484.7290&rep=rep1&type=pdf)

[基于密度的聚类中心自动确定的混合属性数据聚类算法研究](https://www.cn-ki.net/doc_detail?dbcode=CJFQ&filename=MOTO201510011)

[dbscan: Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms](https://cran.r-project.org/web/packages/dbscan/index.html)

A Dynamic Method for Discovering Density Varied Clusters 

[Outliers and the correlation coefficient](http://andreykostenko.com/outliers-and-the-correlation-coefficient/)

[MVN: Multivariate Normality Tests](https://cran.r-project.org/web/packages/MVN/vignettes/MVN.pdf)

[MVN: An R Package for Assessing Multivariate Normality](https://journal.r-project.org/archive/2014-2/korkmaz-goksuluk-zararsiz.pdf)

```R
install.packages('knitr', repos = c('http://rforge.net', 'http://cran.rstudio.org'),
             type = 'source')
```

```R
test1 <- data.frame(a=rnorm(50, 0, 1), b=rnorm(50, 3, 2))
```



[Multivariate Negative Binomial Models for Insurance Claim Counts](https://www.actuaries.org/mexico2012/papers/Shi_Valdez.pdf)

[Practical Guide to Principal Component Analysis (PCA) in R & Python](https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/)

而loadings是用主成分反表示X变量的系数，也可叫载重，在matlab中是score。我的理解就是新的主成分变量占原来变量的比重。


大部分的聚类的出发点是考虑细胞之间的相似性，聚类，再找marker。
但是这必然面临者一个问题，怎么选择聚类的k，以及如何解读每个cluster，所以必然会寻找marker。不管什么聚类，都非常依赖与特征选择，这个太主观，而且风险很大。
为什么不换个角度，直接寻找能够定义一个subgroup的marker？然后再来定义一个cluster？这样才是从生物学的角度出发的方法。优势：直观、其他的可能。（急需一个数据集来证明）

一个subtype偶然出现的概率，blast偶然比对上的概率。e-value
我的工具必须被数学化、统计化，并添加一个实用的模型。
https://en.wikipedia.org/wiki/Distance_correlation
有哪几种可能？最简单的，有zero的，梯度的，都表达，但是有一群高表达的marker，考虑各种能够成为marker group的可能性。
https://www.cnblogs.com/arachis/p/Similarity.html

最可靠的应该是binary化，再计算距离，这样有什么劣势吗？这样对一些典型的marker能够做非常有效的鉴定，但是不能鉴定高低的差异。因为binary依赖阈值得设置。

但如果考虑表达值得化，就复杂了，必须normalize，用哪种距离度量比较好？
energy package for R
我的方法是否能够通过各种测试，解决其他工具解决不了的问题！！！模拟数据搭配真实数据来解释。

与现有的聚类工具对比，时间、资源消耗、准确度。nature method的套路。

现在我的做法是先找一些local的中心，然后依次为中心来搜索，直至迭代结束。怎么找local的中心也是个问题。
（科研就是把浅显的东西给搞深奥，好装逼忽悠）

找到subgroup后就可以开始merge一些共同的subgroup了。

最好结合一下单细胞的数据特征，不然就不够热门了。

现有的很多方法都可以做，有些方法有很多明显的缺陷，如何评估呢？

多看点nature method的单细胞文章，找找灵感。
各个平台的都测试一下，做得user friendly些，争取更多的引用。争取引爆引用。

# Firstly try NBT

[Boosting the power of single-cell analysis](https://www.nature.com/articles/nbt.4131)
[Spatial reconstruction of single-cell gene expression data](https://www.nature.com/articles/nbt.3192)
[Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain](https://www.nature.com/articles/nbt.4038)
[Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors](https://www.nature.com/articles/nbt.4091)
[Single-cell analysis at the threshold](https://www.nature.com/articles/nbt.3721)


# Secondly try Nature method

