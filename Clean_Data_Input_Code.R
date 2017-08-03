setwd('~/Data-Science/Project-work/Lychee-Project/')
library(WGCNA)
options(stringsAsFactors = FALSE)


smlseed <- read.csv('Small_seeded_HSP.csv')
lrgseed <- read.csv('Large_seeded_HSP.csv')

names(smlseed) <- c("Contig", "Day6", "Day14", "Day0")
names(lrgseed) <- c("Contig", "Day6", "Day14", "Day0")

smlseed2 <- t(smlseed)
lrgseed2 <- t(t(lrgseed))


## For a special list
#gene_list <- read.table('Gene_list.txt', header = FALSE, sep = '\n')
#gene_list <- as.vector(gene_list, mode = 'character')
#gene_index <- match(gene_list[, 1], smlseed$Contig)

#smlseed <- smlseed[gene_index, ]
#lrgseed <- lrgseed[gene_index, ]



## Making multiset Data
nsets <- 2
seedlabels <- c('small_seed', 'large_seed')
multiExp <- vector(mode = "list", length=nsets)
multiExp[[1]] <- list(data = as.data.frame(t(smlseed[-1])))
names(multiExp[[1]]$data) <- smlseed$Contig
rownames(multiExp[[1]]$data) <- names(smlseed)[-1]
multiExp[[2]] <- list(data = as.data.frame(t(lrgseed[-1])))
names(multiExp[[2]]$data) <- lrgseed$Contig
rownames(multiExp[[2]]$data) <- names(lrgseed)[-1]

## Taking log of expression data in order to normalize it
library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

multiExp[[1]]$data <- log2(multiExp[[1]]$data + 1)
multiExp[[2]]$data <- log2(multiExp[[2]]$data + 1)



 #Calculation of variance of Genes across samples
var.small <- vector(mode = "numeric", length = dim(multiExp[[1]]$data)[2])
var.large <- vector(mode = "numeric", length = dim(multiExp[[2]]$data)[2])
for(i in 1:dim(multiExp[[1]]$data)[2]) {
  var.small[i] <- var(multiExp[[1]]$data[, i])
  var.large[i]  <- var(multiExp[[2]]$data[, i])
  
}
var.small.sorted <- sort(var.small, decreasing = TRUE)
var.large.sorted <- sort(var.large, decreasing = TRUE)

stopCluster(cl)

contig_to_keep <- match(var.small.sorted[1:8000], var.small)

multiExp[[1]]$data <- multiExp[[1]]$data[, contig_to_keep]
multiExp[[2]]$data <- multiExp[[2]]$data[, contig_to_keep]

#threshold.variance <- 15.0
#high.varGenes <- (var.small > threshold.variance) | (var.large > threshold.variance)

## Retaining only those genes who have high variance
#multiExp[[1]]$data <- multiExp[[1]]$data[, high.varGenes]
#multiExp[[2]]$data <- multiExp[[2]]$data[, high.varGenes]




# ### Checking Genes whose counts are less than 5 across all samples
# noiseGenesml <- vector(mode = "logical", length = dim(smlseed2)[1])
# noiseGenelrg <- vector(mode = "logical", length = dim(lrgseed2)[1])
# for(i in 1:dim(smlseed2)[1]) {
#   noiseGenesml[i] <- max(smlseed2[i, c(-1)]) < 5
#   noiseGenelrg[i] <- max(lrgseed2[i, c(-1)]) < 5
#   
# }
# table(noiseGenesml == noiseGenelrg)
# NoiseGenes.Consensus <- (noiseGenesml & noiseGenelrg)
# 
# for(i in 1:nsets) {
#   multiExp[[i]]$data <- multiExp[[i]]$data[, !NoiseGenes.Consensus]
# }
# 
# 
# 
# 
# var.vector.small <- vector(mode = "numeric", length = dim(multiExp[[1]]$data)[2])
# var.vector.large <- vector(mode = "numeric", length = dim(multiExp[[2]]$data)[2])
# for(i in 1:dim(multiExp[[1]]$data)[2]) {
#   var.vector.small[i] <- var(multiExp[[1]]$data[, i])
#   var.vector.large[i] <- var(multiExp[[2]]$data[, i])
#   
#   }

## Checking if multiExp has correct multiset structure
exprSize <- checkSets(multiExp)
exprSize

#gsg = goodSamplesGenesMS(multiExpr, verbose = 5)
#gsg$allOK


 ## Gene Clustering using average linkage hierarchical clustering
seedsampletrees <- list()
for (set in 1:exprSize$nSets) {
  seedsampletrees[[set]] <- hclust(dist(multiExp[[set]]$data), method = "average") 
}
##Plotting the results
pdf(file = "Seedsampleclustering.pdf", height = 12, width = 12)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nsets)
  plot(seedsampletrees[[set]], main = paste("Sample clustering across all genes in",
       seedlabels[set]), xlab="", sub="", cex = 0.7)
dev.off()


nGenes <- exprSize$nGenes
nSamples <- exprSize$nSamples

save(multiExp, nGenes, nSamples, nsets, seedlabels, exprSize,
     file = "Lychee_Consensus_Input_Data.Rdata")

