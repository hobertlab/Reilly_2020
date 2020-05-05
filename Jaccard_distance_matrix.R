Jaccard.Distance <- function (binary.data, binary.start) {
  #Shared.genes matrix
  matrix(data = NA, nrow = nrow(binary.data), ncol = nrow(binary.data), byrow = FALSE) -> Shared.genes
  rownames(binary.data) -> rownames(Shared.genes)
  rownames(binary.data) -> colnames(Shared.genes)
  #shared genes
  for(x in binary.start:nrow(binary.data)) {
    for (i in binary.start:nrow(binary.data)) {
      which(binary.data[x,] == 1) -> narray
      which(binary.data[i,] == 1) -> narray2
      sum(narray %in% narray2) -> Shared.genes[x,i]
    }
  } 
  #unshared.genes matrix
  matrix(data = NA, nrow = nrow(binary.data), ncol = nrow(binary.data), byrow = FALSE) -> unshared.genes
  rownames(binary.data) -> rownames(unshared.genes)
  rownames(binary.data) -> colnames(unshared.genes)
  #not shared genes
  for(x in binary.start:nrow(binary.data)) {
    for (i in binary.start:nrow(binary.data)) {
      which(binary.data[x,] == 1) -> narray
      which(binary.data[i,] == 1) -> narray2
      sum(!(narray %in% narray2),!(narray2 %in% narray)) -> unshared.genes[x,i]
    }
  } 
  #Jaccard similarity matrix
  Shared.genes/(Shared.genes+unshared.genes) -> Jaccard.Similarity
  #Jaccard distance matrix
  as.dist(1 - Jaccard.Similarity) ->> Jaccard.dist.matrix
  
  print(head(Jaccard.dist.matrix))
}

hclust(Jaccard.dist.matrix) -> h
plot(as.dendrogram(h), cex=20)
