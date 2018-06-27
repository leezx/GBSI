# MBSIT

# title of paper: Clustering and trajectory inference of single-cell RNA-seq data by marker-based subgroup identification method

example1: Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning

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

## large scale simulation data
our tool can be used in a specific situation. marker is not unique expressed gene, but the identity of a subgroup.

# Step2: kinds of marker group for pseudotime?

# Algorithem
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

     

Note: learn from the SIMLR.

# build R bioconductor package

# Backgroud

# Reference
[SIMLR: Visualization and analysis of single-cell RNA-seq data by kernel-based similarity learning](https://www.nature.com/articles/nmeth.4207)

[SC3: consensus clustering of single-cell RNA-seq data](https://www.nature.com/articles/nmeth.4236)

[基于局部中心量度的聚类算法研究](http://blog.sciencenet.cn/blog-3273400-1097494.html)

[一种基于K-Means局部最优性的高效聚类算法](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.484.7290&rep=rep1&type=pdf)

[基于密度的聚类中心自动确定的混合属性数据聚类算法研究](https://www.cn-ki.net/doc_detail?dbcode=CJFQ&filename=MOTO201510011)

[dbscan: Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms](https://cran.r-project.org/web/packages/dbscan/index.html)

A Dynamic Method for Discovering Density Varied Clusters 


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



