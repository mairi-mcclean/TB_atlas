### Load required modules

library(NMF)
library(dplyr)
library(igraph)
library(Matrix)
library(ggplot2)
library(CellChat) 
library(patchwork)
library(ggalluvial)
library(reticulate)
options(stringsAsFactors = FALSE)

### Read in data

cellchat <- readRDS('/Users/carlos.lopez/INBOX/tb_circuits/cellchat_objects/CaiY_activeTB-PBMC_cellchat.rds')
cellchat

### Calculate aggregated cell-cell communication

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)
options(repr.plot.width = 40, repr.plot.height = 40)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

### Visualise all pathways in a heatmap

options(repr.plot.width = 5, repr.plot.height = 5)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 14, height = 18, color.heatmap = "YlGnBu")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 14, height = 18, color.heatmap = "YlGnBu")
ht1 + ht2

options(repr.plot.width = 8, repr.plot.height = 10)
pathways.show <- c("LAMININ")
netAnalysis_contribution(cellchat, signaling = pathways.show)
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

options(repr.plot.width = 10, repr.plot.height = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("LAMININ", "MHC-I"))
gg1 + gg2


### Identify global communication patterns

options(repr.plot.width = 15, repr.plot.height = 15)
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,  width = 10, height = 10)

### Identify global communication patterns

options(repr.plot.width = 15, repr.plot.height = 15)
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,  width = 10, height = 10)

### Sankey plot for outgoing

options(repr.plot.width = 40, repr.plot.height = 22.5)
netAnalysis_river(cellchat, pattern = "outgoing", font.size = 2.5, font.size.title = 20)

netAnalysis_dot(cellchat, pattern = "outgoing")

pairLR.pathway <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.pathway[1,] # show one ligand-receptor pair

### Sankey plot for incoming

options(repr.plot.width = 40, repr.plot.height = 22.5)
netAnalysis_river(cellchat, pattern = "incoming", font.size = 2.5, font.size.title = 20)

netAnalysis_dot(cellchat, pattern = "incoming")

pairLR.pathway <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.pathway[1,] # show one ligand-receptor pair

# Hierarchy plot

vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

plotGeneExpression(cellchat, signaling = "LAMININ")

