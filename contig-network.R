setwd('~/Data-Science/Project-work/Lychee-Project/')
library(WGCNA)
options(stringsAsFactors = FALSE)

lnames <- load(file = "Lychee_Consensus_Input_Data.Rdata")
#lnames2 <- load(file = "Network-data.Rdata")



library(doParallel)
#cl <- makeCluster(3)
#registerDoParallel(cl)
                             
# powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
# powertables <- vector(mode = "list", length = 
#)
# for(set in 1:nsets) {
#   powertables[[set]] <- list(data = pickSoftThreshold(multiExp[[set]]$data, 
#                                                       powerVector = powers,
#                                                       verbose = 4 )[[2]])
# }
# collectGarbage()
# 
# ## plotting of results
# colors = c("blue", "red")
# # Will plot these columns of the returned scale free analysis tables
# plotCols = c(2,5,6,7)
# colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
#              "Max connectivity")
# # Get the minima and maxima of the plotted points
# ylim = matrix(NA, nrow = 2, ncol = 4);
# for (set in 1:nSets)
# {
#   for (col in 1:length(plotCols))
#   {
#     if(col == 1) {
#       ylim[1, col] = min(ylim[1, col], -sign(powertables[[set]]$data[,3])*powertables[[set]]$data[, 2], na.rm = TRUE)
#       ylim[2, col] = max(ylim[2, col], -sign(powertables[[set]]$data[,3])*powertables[[set]]$data[, 2], na.rm = TRUE)
#       
#     }
#     else {
#     ylim[1, col] = min(ylim[1, col], powertables[[set]]$data[, plotCols[col]], na.rm = TRUE)
#     ylim[2, col] = max(ylim[2, col], powertables[[set]]$data[, plotCols[col]], na.rm = TRUE)
#     }
#     
# }
# }
# 
# pdf(file = "powerplots.pdf", width = 12, height = 12)
# par(mfcol = c(2,2))
# par(mar = c(4.2, 4.2 , 2.2, 0.5))
# cex1 = 0.7;
# for (col in 1:length(plotCols)) for (set in 1:nSets)
# {
#   if(set==1) {
#     plot(powertables[[set]]$data[,1], -sign(powertables[[set]]$data[,3])*powertables[[set]]$data[,2],
#          xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
#          main = colNames[col])
#     addGrid()
#   }
#     if (col==1)
#   {
#     text(powertables[[set]]$data[,1], -sign(powertables[[set]]$data[,3])*powertables[[set]]$data[,2],
#          labels=powers,cex=cex1,col=colors[set])
#   } else 
#     text(powertables[[set]]$data[,1], powertables[[set]]$data[,plotCols[col]],
#          labels=powers,cex=cex1,col=colors[set])
#   
#  
#     if (col==1)
#   {
#     legend("bottomright", legend = seedlabels, col = colors, pch = 2) 
#   } else
#     legend("topright", legend = seedlabels, col = colors, pch = 2) 
# }
# dev.off()
# 
#  ##Here I'm going to apply WGCNA on the obtain non-scale free network

#  softpower = 15;
#  similarity.matrix <- array(dim = c(nsets, nGenes, nGenes))
#  adjacency.matrix <-  array(dim = c(nsets, nGenes, nGenes))
#  for (set in 1:nsets) {
#    similarity.matrix[set, , ] = abs(cor(multiExp[[set]]$data, use = "pairwise.complete.obs"))
#    adjacency.matrix[set, , ] = similarity.matrix[set, , ]^softpower
#  }
# # 
#  TOM.distmatrix <- array(0, dim = c(nsets, nGenes, nGenes))
#  for(set in 1:nsets) {
#    TOM.distmatrix[set, , ] <- TOMdist(adjacency.matrix[set, , ], verbose = 4)
#  }
#  
# save(similarity.matrix, adjacency.matrix, TOM.distmatrix, 
#       file = "network_data.Rdata")

# lnames = load(file = "network_data.Rdata")
# lnames2 = load(file = "Lychee_Consensus_Input_Data.Rdata")
# smallseed_corr <- similarity.matrix[1, ,]
# largeseed_corr <- similarity.matrix[2, ,]
# rm(adjacency.matrix,  similarity.matrix)
# smlseed.connectivity <- softConnectivity.fromSimilarity(smallseed_corr, 
#                                       type = "unsigned", power = 15,
#                                       verbose = 4)
# lrgseed.connectivity <- softConnectivity.fromSimilarity(largeseed_corr, 
#                                                        type = "unsigned", power = 15,
#                                                        verbose = 4)
# connectivities <- data.frame(Small_Seed = smlseed.connectivity,
#                              Large_Seed = lrgseed.connectivity)
# write.csv(connectivities, file = 'contig_Connectivities.csv')
# 
# rm(smlseed.connectivity, lrgseed.connectivity, connectivities, smallseed_corr, largeseed_corr )

## For small seeded data
#dissTOM <- TOM.distmatrix[1, , ]

dataExp <- multiExp[[1]]$data
similarity.matrix <- cor(dataExp, use = 'p')
corr.matrix <- as.data.frame(similarity.matrix)
names(corr.matrix) <- names(dataExp)
rownames(corr.matrix) <- names(dataExp)
write.csv(corr.matrix, file = 'Small_Seeded_Pearson_Correlation_Matrix.csv')
rm(corr.matrix)

smlseed.connectivity <- vector(mode = "numeric", length = dim(dataExp)[2])
for(i in 1:length(smlseed.connectivity)) {
  smlseed.connectivity[i] = sum(abs(similarity.matrix[i, ])^15)
}
names(smlseed.connectivity) <- names(dataExp)

write.csv(data.frame(SmallSeed_Connectivity = smlseed.connectivity),
          file = "Small_Seeded_Contig_Connectivities.csv")
rm(smlseed.connectivity)

dissTOM <- TOMdist(abs(similarity.matrix)^15, verbose = 5)

dissTOM_data <- as.data.frame(dissTOM)
names(dissTOM_data) <- names(dataExp)
rownames(dissTOM_data) <- names(dataExp)
write.csv(dissTOM_data, file = 'Small_Seeded_TOM.csv')
rm(dissTOM_data)

geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

 ## Merging of modules whose expression profiles are very similar
 # Calculate eigengenes
MEList = moduleEigengenes(dataExp, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(dataExp, dynamicColors,
                          cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
#table(mergedColors)
mergedMEs = merge$newMEs
#sizeGrWindow(12, 9)
pdf(file = "./Small_Seeded_Modules.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors_small = mergedColors
write.csv(data.frame(Contigs = names(dataExp), Modules = moduleColors_small),
          file = 'Small_Seeded_Contig_Module_Assignment.csv')
dev.off()
rm(dataExp, similarity.matrix, dissTOM, geneTree)
 
# 
 ## For large seeded data

lrgseed <- read.csv('Large_seeded_HSP.csv')
names(lrgseed) <- c("Contig", "Day6", "Day14", "Day0")

## For  special list
gene_list <- read.table('Gene_list.txt', header = FALSE, sep = '\n')
#gene_list <- as.vector(gene_list, mode = 'character')
gene_index <- match(gene_list[, 1], lrgseed$Contig)
lrgseed <- lrgseed[gene_index, ]

dataExp <- as.data.frame(t(lrgseed[-1]))
names(dataExp) <- lrgseed$Contig
rownames(dataExp) <- names(lrgseed)[-1]
dataExp <- log2(dataExp + 1)
var.large <- vector(mode = "numeric", length = dim(dataExp)[2])
for(i in 1:dim(dataExp)[2]) {
  var.large[i]  <- var(dataExp[, i])
}
# 
var.large.sorted <- sort(var.large, decreasing = TRUE)
contig_to_keep <- match(var.large.sorted[1:8000], var.large)
dataExp <- dataExp[, contig_to_keep]

#dataExp = multiExp[[2]]$data
similarity.matrix <- cor(dataExp, use = 'p')
corr.matrix <- as.data.frame(similarity.matrix)
names(corr.matrix) <- names(dataExp)
rownames(corr.matrix) <- names(dataExp)
write.csv(corr.matrix, file = 'Large_Seeded_Pearson_Correlation_Matrix.csv')
rm(corr.matrix)

Lrgseed.connectivity <- vector(mode = "numeric", length = dim(dataExp)[2])
for(i in 1:length(Lrgseed.connectivity)) {
  Lrgseed.connectivity[i] = sum(abs(similarity.matrix[i, ])^15)
}

names(Lrgseed.connectivity) <- names(dataExp)

write.csv(data.frame(LrgSeed_Connectivity = Lrgseed.connectivity),
          file = "Large_Seeded_Contig_Connectivities.csv")
rm(Lrgseed.connectivity)
dissTOM <- TOMdist(abs(similarity.matrix)^15, verbose = 5)

dissTOM_data <- as.data.frame(dissTOM)
names(dissTOM_data) <- names(dataExp)
rownames(dissTOM_data) <- names(dataExp)
write.csv(dissTOM_data, file = 'Large_Seeded_TOM.csv')
rm(dissTOM_data)

geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
 ## Merging of modules whose expression profiles are very similar
 # Calculate eigengenes
MEList = moduleEigengenes(dataExp, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "p")
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(dataExp, dynamicColors,
                          cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
write.csv(data.frame(Contigs = names(dataExp), Modules = mergedColors),
          file = 'Large_Seeded_Contig_Module_Assignment.csv')
#connectivities <- data.frame(Small_Seed = smlseed.connectivity,
#                              Large_Seed = lrgseed.connectivity)
#write.csv(connectivities, file = 'contig_Connectivities.csv')

#rm(dataExp, similarity.matrix, dissTOM, geneTree)
table(mergedColors)
mergedMEs = merge$newMEs
# #sizeGrWindow(12, 9)
pdf(file = "./Large_Seeded_Modules.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                     c("Dynamic Tree Cut", "Merged dynamic"),
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)
dev.off()
 moduleColors = mergedColors
# # 
#  Consensus.TOMdist1 <- TOM.distmatrix[1, , ] 
#    #pmin(TOM.distmatrix[1, , ], TOM.distmatrix[1, , ]) 
#  constree <- hclust(as.dist(Consensus.TOMdist), method = "average")
#  minModulesize <- 100
#  unmergedlabels <- cutreeDynamic(dendro = constree, distM = Consensus.TOMdist,
#                                  deepSplit = 2, cutHeight = 0.995,
#                                  minClusterSize = minModulesize,
#                                  pamRespectsDendro = FALSE )
#  
#  unmergedcolors = labels2colors(unmergedlabels)
#  
#  unmergedMEs = multiSetMEs(multiExp, colors = NULL, universalColors = unmergedcolors)
#  consMEDiss = consensusMEDissimilarity(unmergedMEs)
#  consMETree = hclust(as.dist(consMEDiss), method = "average")
#  
# # 
#  sizeGrWindow(7,6)
#  par(mfrow = c(1,1))
#  plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
#       xlab = "", sub = "")
#  abline(h=0.25, col = "red")
#  merge = mergeCloseModules(multiExp, unmergedlabels,
#                            cutHeight = 0.25, verbose = 3)
# # 
#  modulelabels = merge$colors;
#  table(modulelabels)
#  moduleColors = labels2colors(modulelabels)
#  consMEs = merge$newMEs;
# # 
#  pdf(file = "SmallSeed_Contig_Modules.pdf", width = 12, height = 12)
#  par(mfcol = c(2,2))
#  par(mar = c(4.2, 4.2 , 2.2, 0.5))
#  cex1 = 0.7
#  plotDendroAndColors(constree, cbind(unmergedcolors, moduleColors),
#                      c("Unmerged", "Merged"),
#                      dendroLabels = FALSE, hang = 0.03,
#                      addGuide = TRUE, guideHang = 0.05)
#  
#  dev.off()
# #
 
 
#stopCluster(cl) 
# save(similarity.matrix, adjacency.matrix, TOM.distmatrix, powertables, 
#      file = "Network-data.Rdata")
# save(consMEs, moduleColors, modulelabels, constree, file = "Lychee_Consensus_Modules.Rdata")
