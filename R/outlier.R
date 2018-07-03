setwd("/Users/surgery/Project/HOME/github/MBSIT/R")
load("/Users/surgery/Project/HOME/1-projects/1.scRNA-seq/2-smart-seq/eachGroup/HSCR_5c3.Rdata")

logcounts(tmp_group)

cellCountPerGene <- rowSums(counts(tmp_group)>=5)
geneCountPerCell <- colSums(counts(tmp_group)>=5)
totalReadCount <- colSums(counts(tmp_group))
cvPerGene <- apply(counts(tmp_group), 1, sd)/apply(counts(tmp_group), 1, mean)
# cvPerCell <- apply(counts(tmp_group), 2, sd)/apply(counts(tmp_group), 2, mean)

