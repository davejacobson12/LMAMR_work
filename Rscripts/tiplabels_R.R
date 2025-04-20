length(row.names(inFile)) -> totalTax
length(colnames(inFile)) -> totalSamp

finalMat <- matrix(ncol = 2, nrow = totalTax)
for (j in 1:totalTax){
row.names(inFile)[j] -> taxon

holder <- numeric(length = 1000)
for (i in 1:totalSamp){
  if(grepl("Had", colnames(inFile)[i]) == "TRUE"){
    inFile[j,i] -> holder[i]
  }
}
sum(holder) -> Had

holder <- numeric(length = 1000)
for (i in 1:totalSamp){
  if(grepl("bftm", colnames(inFile)[i]) == "TRUE"){
    inFile[j,i] -> holder[i]
  }
}
sum(holder) -> BFTM

holder <- numeric(length = 1000)
for (i in 1:totalSamp){
  if(grepl("Ita", colnames(inFile)[i]) == "TRUE"){
    inFile[j,i] -> holder[i]
  }
}
sum(holder) -> Ita

holder <- numeric(length = 1000)
for (i in 1:totalSamp){
    if(grepl("Ind", colnames(inFile)[i]) == "TRUE"){
        inFile[j,i] -> holder[i]
    }
}
sum(holder) -> Euro

holder <- numeric(length = 1000)
for (i in 1:totalSamp){
    if(grepl("SRR399", colnames(inFile)[i]) == "TRUE"){
        inFile[j,i] -> holder[i]
    }
}
sum(holder) -> Mongolia

max(c(BFTM, Had, Ita, Euro, Mongolia)) -> Max
matrix(ncol=2, nrow=1) -> tipMatrix
tipMatrix[1,1] <- taxon

if(Max == BFTM){
  tipMatrix[1,2] <- "bftm"
}

if (Max == Had){
  tipMatrix[1,2] <- "had"
}

if (Max == Ita){
  tipMatrix[1,2] <- "ita"
}

if (Max == Euro){
    tipMatrix[1,2] <- "euro"
}

if (Max == Mongolia){
    tipMatrix[1,2] <- "mongolia"
}

tipMatrix -> finalMat[j,]
}
colnames(finalMat) <- c("Taxon", "Pop")
as.data.frame(finalMat) -> tipDF


library(ggtree)
lirary(ggplot2)


