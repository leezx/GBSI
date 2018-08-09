

## goodSamplesGenes

Iterative filtering of samples and genes with too many missing entries

Filter outlier cells and noisy genes before network construction

## pickSoftThreshold

Analysis of scale free topology for soft-thresholding

Analysis of scale free topology for multiple soft thresholding powers. The aim is to help the user pick an appropriate soft-thresholding power for network construction.

## blockwiseModules

Automatic network construction and module detection

## adjacency

Calculates (correlation or distance) network adjacency from given expression data or from a similarity.

WGCNA constructs two matrices, first it defines a correlation matrix up to a power beta so the degree distribution will fit a small-word network.

This matrix gives only information about the expression correlation between genes.

## TOMsimilarity

Turn adjacency into topological overlap

Calculation of the topological overlap matrix, and the corresponding dissimilarity, from a given adjacency matrix.

WGCNA thinks that co-expression is not enough and the similarity between genes should be reflected at the expression and the network topology level.

 This is why it defines the TOM matrix which uses the co-expression Adjacency matrix and build another adjacency matrix that considers topological similarity.

Normally the TOM matrix is the final result of WGCNA.



I understand that the adjacency matrix is calculated using Pearson correlation as the metric and it gives us correlation between each pair of genes from the input normalized   expression matrix. But, if you look at the TOM matrix, we do not see any gene names. Can you tell me exactly how TOM is calculated. e.g. the first matrix uses correlation as the metric and the values in the first adjacency matrix are actually Pearson correlation values

what is the metric used for TOM. What do the values in TOM matrix actually mean, and what do the first row and first column of TOM matrix represent ?

https://www.researchgate.net/post/What_do_adjacency_matrix_and_Topology_Overlap_Matrix_from_WGCNA_package_tell_about_the_data



## moduleEigengenes

Calculates module eigengenes (1st principal component) of modules in a given single dataset.

## mergeCloseModules

Merges modules in gene expression networks that are too close as measured by the correlation of their eigengenes.

## orderMEs

Put close eigenvectors next to each other

## corPvalueStudent

Calculates Student asymptotic p-value for given correlations.

