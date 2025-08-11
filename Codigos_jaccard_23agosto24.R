
library(gplots)


heatmap.2(as.matrix(vegdist(x = t(t.my.asvs[,colSums(apply(t.my.asvs, MARGIN = 2, FUN = function(x){ifelse(x>0,1,0)}))>0]), method = "jaccard", binary = T))[1:10, 1:10])


heatmap.2(as.matrix(vegdist(x = t(t.my.asvs[,colSums(apply(t.my.asvs, MARGIN = 2, FUN = function(x){ifelse(x>0,1,0)}))>0]), method = "jaccard", binary = T)))


heatmap.2(as.matrix(vegdist(x = t(t.my.asvs[,colSums(apply(t.my.asvs, MARGIN = 2, FUN = function(x){ifelse(x>0,1,0)}))>0]), method = "jaccard", binary = T)), trace = "none", col = viridis(n = 10, direction = -1), density.info = "none")




cuf.off <- 13
my.jac <- as.matrix(vegdist(x = t(t.my.asvs[,colSums(apply(t.my.asvs, MARGIN = 2, FUN = function(x){ifelse(x>0,1,0)}))>cuf.off]), method = "jaccard", binary = T))
                                  
                                  
library(pheatmap)                              

