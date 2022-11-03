## Script to compare CCI between PBMCs from different TB stages

### Load required packages

library(igraph)
library(ggplot2)                  
library(CellChat)
library(patchwork)

### Read in datasets

TBpbmc.healthy <- readRDS('/Users/carlos.lopez/INBOX/tb_circuits/cellchat_objects/CaiY_Healthy-PBMC_cellchat.rds')
TBpbmc.healthy <- updateCellChat(TBpbmc.healthy)
TBpbmc.healthy

TBpbmc.active <- readRDS('/Users/carlos.lopez/INBOX/tb_circuits/cellchat_objects/CaiY_activeTB-PBMC_cellchat.rds')
TBpbmc.active <- updateCellChat(TBpbmc.active)
TBpbmc.active

TBpbmc.latent <- readRDS('/Users/carlos.lopez/INBOX/tb_circuits/cellchat_objects/CaiY_latentTB-PBMC_cellchat.log.rds')
TBpbmc.latent <- updateCellChat(TBpbmc.latent)
TBpbmc.latent

TBpfmc.active <- readRDS('/Users/carlos.lopez/INBOX/tb_circuits/cellchat_objects/CaiY_activeTB-PFMC_cellchat.rds')
TBpfmc.active <- updateCellChat(TBpfmc.active)
TBpfmc.active

### Lift up CellChat object and merge together

group.new = levels(TBpbmc.healthy@idents)
TBpbmc.healthy <- liftCellChat(TBpbmc.healthy, group.new)

object.list <- list(PBMCh = TBpbmc.healthy, PBMCa = TBpbmc.active, PBMCl = TBpbmc.latent, PFMCa = TBpfmc.active)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

### Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

### Differential number of interactions or interaction strength among different cell populations

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(2,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2


### Visualize the inferred signaling network using the lifted object

pathways.show <- c("CXCL", 'CLEC')

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
vertex.receiver = seq(1,10) # Left portion of hierarchy plot the shows signaling to dermal cells and right portion shows signaling to epidermal cells
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, vertex.receiver = vertex.receiver, edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}