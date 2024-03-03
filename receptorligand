library(Seurat)
library(dplyr)
library(harmony)
library(symphony)
library(ggplot2)
library(RColorBrewer)
library(scRepertoire)
library(patchwork)
library(cowplot)
library(pheatmap)
library(ggpubr)
library(rstatix)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
library(UpSetR)
source("/Rao_Lab/AMP_Phase_2/manuscript/AMP_manuscript_functions.R")


setwd("/Desktop/Rao_Lab/AMP_Phase_2/manuscript/analysis/")


clust_30PC_res03_seurat <- readRDS("/Rao_Lab/AMP_Phase_2/data/mainPipeline/02clust_30PC_res03_seurat.rds")
TBclustercolors <- c(colorRampPalette(brewer.pal(8, "Blues"))(10)[3:10], colorRampPalette(brewer.pal(8, "Purples"))(10)[8], colorRampPalette(brewer.pal(8, "Reds"))(10)[4:9], colorRampPalette(brewer.pal(8, "Greens"))(6)[3:4])

Idents(clust_30PC_res03_seurat) <- clust_30PC_res03_seurat$tissue_source
current.cluster.ids <- c("PBL", "Syn")
new.cluster.ids <- c("PBL", "SYN")
clust_30PC_res03_seurat@active.ident <- plyr::mapvalues(x = clust_30PC_res03_seurat@active.ident, from = current.cluster.ids, to = new.cluster.ids)
clust_30PC_res03_seurat$tissue_source <- clust_30PC_res03_seurat@active.ident

clust_30PC_res03_seurat$tissue_source <- factor(clust_30PC_res03_seurat$tissue_source, levels = c("SYN", "PBL"))

Idents(clust_30PC_res03_seurat) <- clust_30PC_res03_seurat$ident
current.cluster.ids <- levels(clust_30PC_res03_seurat$ident)
new.cluster.ids <- c("T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "T cells", "Proliferating", "B cells", "B cells", "B cells", "B cells", "B cells", "B cells", "Fibroblasts", "Monocytes")
clust_30PC_res03_seurat@active.ident <- plyr::mapvalues(x = clust_30PC_res03_seurat@active.ident, from = current.cluster.ids, to = new.cluster.ids)
clust_30PC_res03_seurat$broadident <- clust_30PC_res03_seurat@active.ident

Idents(clust_30PC_res03_seurat) <- clust_30PC_res03_seurat$sample_ID
clust_30PC_res03_seurat <- subset(clust_30PC_res03_seurat, idents = "300_0415", invert = TRUE)

current.cluster.ids <- unique(clust_30PC_res03_seurat$sample_ID)
new.cluster.ids <- c("RA01", "RA02", "RA03", "RA04", "RA05", "RA06", "RA07", "RA08", "RA09", "RA10", "RA11", "RA12")
clust_30PC_res03_seurat@active.ident <- plyr::mapvalues(x = clust_30PC_res03_seurat@active.ident, from = current.cluster.ids, to = new.cluster.ids)
clust_30PC_res03_seurat$sample_ID_Pub <- clust_30PC_res03_seurat@active.ident

clust_30PC_res03_seurat$sampleTissue <- paste0(clust_30PC_res03_seurat$sample_ID, "_", clust_30PC_res03_seurat$tissue_source)
clust_30PC_res03_seurat$sampleTissue_Pub <- paste0(clust_30PC_res03_seurat$sample_ID_Pub, "_", clust_30PC_res03_seurat$tissue_source)

patient_metadata <- read.csv("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/metadata/amp_ra_patient_metadata_perCell.csv", row.names = 1)
clust_30PC_res03_seurat <- AddMetaData(clust_30PC_res03_seurat, patient_metadata)

bcellsmeta <- read.csv("~/Downloads/Bcell_ident_noSCT.csv")
bcellsmeta <- bcellsmeta %>% select(seurat_key, ident)
bcellsmeta$broadident <- "B cells"
bcellsmeta$lineage <- "B"
bcellsmeta <- bcellsmeta %>% filter(!is.na(ident))
old.b.names <- unique(bcellsmeta$ident)
new.b.names <- c("B: Naive(IgD-high)", "B: MT-high", "B: Naive(IgD-low)", "B: Memory", "B: ABC", "Plasma cells", "Plasmablasts", "B: Naive(HSPA1B+)", "B: Activated", "B: Clonal(LILRA4+)")
bcellsmeta$ident <- plyr::mapvalues(x = bcellsmeta$ident, from = old.b.names, to = new.b.names)

tcellsmeta <- read.csv("~/Downloads/amp2_finalCD8CD4idents.csv")
tcellsmeta$seurat_key <- tcellsmeta$X
tcellsmeta$ident <- tcellsmeta$clusters_names
tcellsmeta$broadident <- "T cells"
tcellsmeta <- tcellsmeta %>% select(seurat_key, ident, broadident, lineage)

#find intersecting barcodes
bt_intersection <- intersect(bcellsmeta$seurat_key, tcellsmeta$seurat_key)

#one intersecting cell, will remove from both though
bcellsmeta <- bcellsmeta %>% filter(!seurat_key %in% bt_intersection)
tcellsmeta <- tcellsmeta %>% filter(!seurat_key %in% bt_intersection)

btcellsmeta_bind <- rbind(bcellsmeta, tcellsmeta)
rownames(btcellsmeta_bind) <- btcellsmeta_bind$seurat_key
btcellsmeta_bind$seurat_key <- NULL

clust_30PC_res03_seurat <- AddMetaData(clust_30PC_res03_seurat, metadata = btcellsmeta_bind)
clust_30PC_res03_seurat <- subset(clust_30PC_res03_seurat, subset = broadident %in% c("B cells", "T cells"))




#CellChat between B, CD4, CD8 broad identities
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

Idents(clust_30PC_res03_seurat) <- clust_30PC_res03_seurat$lineage
data.input <- GetAssayData(clust_30PC_res03_seurat, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_broad_all <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_broad_all@DB <- CellChatDB
cellchat_broad_all <- subsetData(cellchat_broad_all) # This step is necessary even if using the whole database
cellchat_broad_all <- identifyOverExpressedGenes(cellchat_broad_all)
cellchat_broad_all <- identifyOverExpressedInteractions(cellchat_broad_all)

cellchat_broad_all <- computeCommunProb(cellchat_broad_all)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_broad_all <- filterCommunication(cellchat_broad_all, min.cells = 10)
cellchat_broad_all <- computeCommunProbPathway(cellchat_broad_all)

cellchat_broad_all
cellchat_broad_all <- aggregateNet(cellchat_broad_all)
cellchat_broad_all <- netAnalysis_computeCentrality(cellchat_broad_all, slot.name = "netP")

groupSize <- as.numeric(table(cellchat_broad_all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_broad_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_broad_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf("R-L_broad_all_circle.pdf", width = 4, height = 4.5)
netVisual_circle(cellchat_broad_all@net$count,  weight.scale = T, label.edge= T, title.name = "Number of Interactions\nAbove Significance Cutoff", arrow.width = 3, )
dev.off()

pdf("R-L_broad_all_heatmap.pdf", width = 3.8, height = 3)
netVisual_heatmap(cellchat_broad_all, color.heatmap = "Reds")
dev.off()

netVisual_bubble(cellchat_broad_all, remove.isolate = FALSE, angle.x = 45) + 
  theme(axis.text = element_text(size = 7)) + coord_flip()
ggsave("R-L_broad_all_bubbleAllInteractions.eps", width = 7.5, height = 3)

netVisual_bubble(cellchat_broad_all, remove.isolate = FALSE, signaling = c("BTLA", "CD40", "IFN-II"))  + 
  theme(axis.text = element_text(size = 7)) 
ggsave("R-L_broad_all_bubbleSelectedInteractions.eps", width = 5, height = 2)





#CellChat between B broad and CD4 narrow identities
clust_30PC_res03_seurat_BCD4 <- subset(clust_30PC_res03_seurat, subset = lineage == "CD8", invert = T)
clust_30PC_res03_seurat_BCD4$ident_Bbroad <- ifelse(clust_30PC_res03_seurat_BCD4$lineage == "B", "B", clust_30PC_res03_seurat_BCD4$ident)

Idents(clust_30PC_res03_seurat_BCD4) <- clust_30PC_res03_seurat_BCD4$ident_Bbroad
data.input <- GetAssayData(clust_30PC_res03_seurat_BCD4, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BCD4)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_broadBnarrowCD4_all <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_broadBnarrowCD4_all@DB <- CellChatDB
cellchat_broadBnarrowCD4_all <- subsetData(cellchat_broadBnarrowCD4_all) # This step is necessary even if using the whole database
cellchat_broadBnarrowCD4_all <- identifyOverExpressedGenes(cellchat_broadBnarrowCD4_all)
cellchat_broadBnarrowCD4_all <- identifyOverExpressedInteractions(cellchat_broadBnarrowCD4_all)

cellchat_broadBnarrowCD4_all <- computeCommunProb(cellchat_broadBnarrowCD4_all)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_broadBnarrowCD4_all <- filterCommunication(cellchat_broadBnarrowCD4_all, min.cells = 10)
cellchat_broadBnarrowCD4_all <- computeCommunProbPathway(cellchat_broadBnarrowCD4_all)

cellchat_broadBnarrowCD4_all
cellchat_broadBnarrowCD4_all <- aggregateNet(cellchat_broadBnarrowCD4_all)
cellchat_broadBnarrowCD4_all <- netAnalysis_computeCentrality(cellchat_broadBnarrowCD4_all, slot.name = "netP")

groupSize <- as.numeric(table(cellchat_broadBnarrowCD4_all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_broadBnarrowCD4_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_broadBnarrowCD4_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf("R-L_broadBnarrowCD4_all_circle.pdf", width = 7, height =7)
netVisual_circle(cellchat_broadBnarrowCD4_all@net$count, weight.scale = T, 
                 label.edge= T, title.name = "Number of Interactions\nAbove Significance Cutoff", targets.use = "B", vertex.label.cex = 0.8)
dev.off()

pdf("R-L_broadBnarrowCD4_all_circle_cd4cols.pdf", width = 7, height =7)
netVisual_circle(cellchat_broadBnarrowCD4_all@net$count, weight.scale = T, color.use = c("gray50", cd4cols[c(1,14,3,4,5,11,2,6,10,9,4,7,8,12)]),
                 label.edge= T, title.name = "Number of Interactions\nAbove Significance Cutoff", targets.use = "B", vertex.label.cex = 0.8, arrow.width = 1, arrow.size = 0.05)
dev.off()

pdf("R-L_broadBnarrowCD4_all_heatmap.pdf", width = 5, height = 4.5)
netVisual_heatmap(cellchat_broadBnarrowCD4_all, color.heatmap = "Reds", 
                  color.use = setNames(rep("black", 15), unique(clust_30PC_res03_seurat_BCD4$ident_Bbroad)), font.size = 6)
dev.off()

pdf("R-L_broadBnarrowCD4_all_heatmap_Btarget.pdf", width = 5, height = 4.5)
netVisual_heatmap(cellchat_broadBnarrowCD4_all, color.heatmap = "Reds", targets.use = "B",
                  color.use = setNames(rep("black", 15), unique(clust_30PC_res03_seurat_BCD4$ident_Bbroad)), font.size = 6)
dev.off()

netVisual_heatmap(cellchat_broadBnarrowCD4_all, color.heatmap = "Reds", sources.use = "B")


netVisual_bubble(cellchat_broadBnarrowCD4_all, remove.isolate = FALSE, signaling = c("CXCL", "IFN-II", "BTLA", "CD40"), targets.use = "B")

netVisual_bubble(cellchat_broadBnarrowCD4_all, remove.isolate = FALSE, signaling = c("CXCL", "IFN-II"), targets.use = "B", angle.x = 45) + 
  theme(axis.text = element_text(size = 7))
ggsave("R-L_broadBnarrowCD4_all_bubbleSelectedInteractions.eps", width = 5, height = 2.7)




#CellChat between B narrow and CD4 narrow identities
Idents(clust_30PC_res03_seurat_BCD4) <- clust_30PC_res03_seurat_BCD4$ident
data.input <- GetAssayData(clust_30PC_res03_seurat_BCD4, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BCD4)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_narrowBnarrowCD4_all <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_narrowBnarrowCD4_all@DB <- CellChatDB
cellchat_narrowBnarrowCD4_all <- subsetData(cellchat_narrowBnarrowCD4_all) # This step is necessary even if using the whole database
cellchat_narrowBnarrowCD4_all <- identifyOverExpressedGenes(cellchat_narrowBnarrowCD4_all)
cellchat_narrowBnarrowCD4_all <- identifyOverExpressedInteractions(cellchat_narrowBnarrowCD4_all)

cellchat_narrowBnarrowCD4_all <- computeCommunProb(cellchat_narrowBnarrowCD4_all)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_narrowBnarrowCD4_all <- filterCommunication(cellchat_narrowBnarrowCD4_all, min.cells = 10)
cellchat_narrowBnarrowCD4_all <- computeCommunProbPathway(cellchat_narrowBnarrowCD4_all)

cellchat_narrowBnarrowCD4_all
cellchat_narrowBnarrowCD4_all <- aggregateNet(cellchat_narrowBnarrowCD4_all)
cellchat_narrowBnarrowCD4_all <- netAnalysis_computeCentrality(cellchat_narrowBnarrowCD4_all, slot.name = "netP")

groupSize <- as.numeric(table(cellchat_narrowBnarrowCD4_all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_narrowBnarrowCD4_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_narrowBnarrowCD4_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds")

pdf("R-L_narrowBnarrowCD4_all_heatmap_Btarget.pdf", width = 8, height = 8)
netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds", targets.use = c("B: ABC", "B: Activated", "B: Memory"), 
                  color.use = setNames(rep("black", 24), unique(clust_30PC_res03_seurat_BCD4$ident)), font.size = 6)
dev.off()


pdf("R-L_narrowBnarrowCD4_all_heatmap_BtargetTsource.pdf", width = 6, height = 5.5)
netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds", targets.use = unique(bcellsmeta$ident), sources.use = unique(tcellsmeta$ident),  
                  color.use = setNames(rep("black", 24), unique(clust_30PC_res03_seurat_BCD4$ident)), font.size = 7)
dev.off()


netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds", targets.use = "B: ABC")
netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds", sources.use = "B: ABC")
netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds", targets.use = "CD4: Tph")
netVisual_heatmap(cellchat_narrowBnarrowCD4_all, color.heatmap = "Reds", sources.use = "CD4: Tph")


netVisual_bubble(cellchat_narrowBnarrowCD4_all, remove.isolate = FALSE, signaling = c("CXCL", "CD40", "IFN-II"))
netVisual_bubble(cellchat_narrowBnarrowCD4_all, remove.isolate = FALSE, signaling = c("CXCL", "IFN-II", "CD39"), targets.use = unique(bcellsmeta$ident))

cellchat_narrowBnarrowCD4_all <- netAnalysis_computeCentrality(cellchat_narrowBnarrowCD4_all, slot.name = "netP")
netAnalysis_signalingRole_scatter(cellchat_narrowBnarrowCD4_all, dot.size = c(1,6), label.size = 2.5)
ggsave("R-L_narrowBnarrowCD4_all_clusterStrengthScatter.eps", width = 7, height = 5)




#CellChat between B broad and CD4 narrow identities, split between blood and tissue
clust_30PC_res03_seurat_BCD4_SYN <- subset(clust_30PC_res03_seurat_BCD4, subset = tissue_source == "SYN")
Idents(clust_30PC_res03_seurat_BCD4_SYN) <- clust_30PC_res03_seurat_BCD4_SYN$ident_Bbroad
data.input <- GetAssayData(clust_30PC_res03_seurat_BCD4_SYN, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BCD4_SYN)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_broadBnarrowCD4_SYN <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_broadBnarrowCD4_SYN@DB <- CellChatDB
cellchat_broadBnarrowCD4_SYN <- subsetData(cellchat_broadBnarrowCD4_SYN) # This step is necessary even if using the whole database
cellchat_broadBnarrowCD4_SYN <- identifyOverExpressedGenes(cellchat_broadBnarrowCD4_SYN)
cellchat_broadBnarrowCD4_SYN <- identifyOverExpressedInteractions(cellchat_broadBnarrowCD4_SYN)

cellchat_broadBnarrowCD4_SYN <- computeCommunProb(cellchat_broadBnarrowCD4_SYN)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_broadBnarrowCD4_SYN <- filterCommunication(cellchat_broadBnarrowCD4_SYN, min.cells = 10)
cellchat_broadBnarrowCD4_SYN <- computeCommunProbPathway(cellchat_broadBnarrowCD4_SYN)

cellchat_broadBnarrowCD4_SYN
cellchat_broadBnarrowCD4_SYN <- aggregateNet(cellchat_broadBnarrowCD4_SYN)
cellchat_broadBnarrowCD4_SYN <- netAnalysis_computeCentrality(cellchat_broadBnarrowCD4_SYN, slot.name = "netP")


clust_30PC_res03_seurat_BCD4_PBL <- subset(clust_30PC_res03_seurat_BCD4, subset = tissue_source == "PBL")
Idents(clust_30PC_res03_seurat_BCD4_PBL) <- clust_30PC_res03_seurat_BCD4_PBL$ident_Bbroad
data.input <- GetAssayData(clust_30PC_res03_seurat_BCD4_PBL, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BCD4_PBL)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_broadBnarrowCD4_PBL <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_broadBnarrowCD4_PBL@DB <- CellChatDB
cellchat_broadBnarrowCD4_PBL <- subsetData(cellchat_broadBnarrowCD4_PBL) # This step is necessary even if using the whole database
cellchat_broadBnarrowCD4_PBL <- identifyOverExpressedGenes(cellchat_broadBnarrowCD4_PBL)
cellchat_broadBnarrowCD4_PBL <- identifyOverExpressedInteractions(cellchat_broadBnarrowCD4_PBL)

cellchat_broadBnarrowCD4_PBL <- computeCommunProb(cellchat_broadBnarrowCD4_PBL)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_broadBnarrowCD4_PBL <- filterCommunication(cellchat_broadBnarrowCD4_PBL, min.cells = 10)
cellchat_broadBnarrowCD4_PBL <- computeCommunProbPathway(cellchat_broadBnarrowCD4_PBL)

cellchat_broadBnarrowCD4_PBL
cellchat_broadBnarrowCD4_PBL <- aggregateNet(cellchat_broadBnarrowCD4_PBL)
cellchat_broadBnarrowCD4_PBL <- netAnalysis_computeCentrality(cellchat_broadBnarrowCD4_PBL, slot.name = "netP")

object.list <- list(SYN = cellchat_broadBnarrowCD4_SYN, PBL = cellchat_broadBnarrowCD4_PBL)
cellchat_broadBnarrowCD4_SYNPBLsplit <- mergeCellChat(object.list, add.names = names(object.list))

compareInteractions(cellchat_broadBnarrowCD4_SYNPBLsplit, show.legend = F, group = c(1,2))

netVisual_bubble(cellchat_broadBnarrowCD4_SYNPBLsplit, sources.use = "B",  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat_broadBnarrowCD4_SYNPBLsplit, targets.use = "B",  comparison = c(1, 2), angle.x = 45)




#CellChat between B narrow and CD4 narrow identities, split between blood and tissue
clust_30PC_res03_seurat_BCD4_SYN <- subset(clust_30PC_res03_seurat_BCD4, subset = tissue_source == "SYN")
Idents(clust_30PC_res03_seurat_BCD4_SYN) <- clust_30PC_res03_seurat_BCD4_SYN$ident
data.input <- GetAssayData(clust_30PC_res03_seurat_BCD4_SYN, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BCD4_SYN)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_narrowBnarrowCD4_SYN <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_narrowBnarrowCD4_SYN@DB <- CellChatDB
cellchat_narrowBnarrowCD4_SYN <- subsetData(cellchat_narrowBnarrowCD4_SYN) # This step is necessary even if using the whole database
cellchat_narrowBnarrowCD4_SYN <- identifyOverExpressedGenes(cellchat_narrowBnarrowCD4_SYN)
cellchat_narrowBnarrowCD4_SYN <- identifyOverExpressedInteractions(cellchat_narrowBnarrowCD4_SYN)

cellchat_narrowBnarrowCD4_SYN <- computeCommunProb(cellchat_narrowBnarrowCD4_SYN)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_narrowBnarrowCD4_SYN <- filterCommunication(cellchat_narrowBnarrowCD4_SYN, min.cells = 10)
cellchat_narrowBnarrowCD4_SYN <- computeCommunProbPathway(cellchat_narrowBnarrowCD4_SYN)

cellchat_narrowBnarrowCD4_SYN
cellchat_narrowBnarrowCD4_SYN <- aggregateNet(cellchat_narrowBnarrowCD4_SYN)
cellchat_narrowBnarrowCD4_SYN <- netAnalysis_computeCentrality(cellchat_narrowBnarrowCD4_SYN, slot.name = "netP")

clust_30PC_res03_seurat_BCD4_PBL <- subset(clust_30PC_res03_seurat_BCD4, subset = tissue_source == "PBL")
clust_30PC_res03_seurat_BCD4_PBL <- subset(clust_30PC_res03_seurat_BCD4_PBL, subset = ident == "B: Clonal(LILRA4+)", invert = T)
Idents(clust_30PC_res03_seurat_BCD4_PBL) <- clust_30PC_res03_seurat_BCD4_PBL$ident
data.input <- GetAssayData(clust_30PC_res03_seurat_BCD4_PBL, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BCD4_PBL)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_narrowBnarrowCD4_PBL <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_narrowBnarrowCD4_PBL@DB <- CellChatDB
cellchat_narrowBnarrowCD4_PBL <- subsetData(cellchat_narrowBnarrowCD4_PBL) # This step is necessary even if using the whole database
cellchat_narrowBnarrowCD4_PBL <- identifyOverExpressedGenes(cellchat_narrowBnarrowCD4_PBL)
cellchat_narrowBnarrowCD4_PBL <- identifyOverExpressedInteractions(cellchat_narrowBnarrowCD4_PBL)

cellchat_narrowBnarrowCD4_PBL <- computeCommunProb(cellchat_narrowBnarrowCD4_PBL)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_narrowBnarrowCD4_PBL <- filterCommunication(cellchat_narrowBnarrowCD4_PBL, min.cells = 10)
cellchat_narrowBnarrowCD4_PBL <- computeCommunProbPathway(cellchat_narrowBnarrowCD4_PBL)

cellchat_narrowBnarrowCD4_PBL
cellchat_narrowBnarrowCD4_PBL <- aggregateNet(cellchat_narrowBnarrowCD4_PBL)
cellchat_narrowBnarrowCD4_PBL <- netAnalysis_computeCentrality(cellchat_narrowBnarrowCD4_PBL, slot.name = "netP")

object.list <- list(SYN = cellchat_narrowBnarrowCD4_SYN, PBL = cellchat_narrowBnarrowCD4_PBL)
cellchat_narrowBnarrowCD4_SYNPBLsplit <- mergeCellChat(object.list, add.names = names(object.list))

compareInteractions(cellchat_narrowBnarrowCD4_SYNPBLsplit, show.legend = F, group = c(1,2)) +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = "Source")
ggsave("R-L_narrowBnarrowCD4_splitTissue_totalInteractions.eps", width = 3, height = 3)

netVisual_bubble(cellchat_narrowBnarrowCD4_SYNPBLsplit, targets.use = "B: ABC", 
                 sources.use = unique(tcellsmeta$ident),  comparison = c(1, 2), angle.x = 45) +
  theme(axis.text = element_text(size = 6))
ggsave("R-L_narrowBnarrowCD4_splitTissue_bubble_ABCtarget.eps", width = 8, height = 4)

netVisual_bubble(cellchat_narrowBnarrowCD4_SYNPBLsplit, sources.use = "B: ABC", 
                 targets.use = unique(tcellsmeta$ident),  comparison = c(1, 2), angle.x = 45) +
  theme(axis.text = element_text(size = 6))
ggsave("R-L_narrowBnarrowCD4_splitTissue_bubble_ABCsource.eps", width = 8, height = 6)

netVisual_bubble(cellchat_narrowBnarrowCD4_SYNPBLsplit, targets.use = c("B: ABC", "B: Activated", "B: Memory"), sources.use = c("CD4: Tph", "CD4: Tfh/Tph"),  comparison = c(1, 2), angle.x = 45)
ggsave("R-L_narrowBnarrowCD4_splitTissue_bubble_TphTfhsource.eps", width = 6, height = 4)

rankNet(cellchat_narrowBnarrowCD4_SYNPBLsplit, mode = "comparison", stacked = T, do.stat = F, title = "Differential Path Usage") + 
  scale_fill_manual(values = c( "#88a0c4", "#f2a21e"), name = "Source") +
  theme(axis.text.y = element_text(size = 8, colour = "black"),
        axis.title.x = element_blank())
ggsave("R-L_narrowBnarrowCD4_splitTissue_differentialUsage.eps", width = 3.5, height = 5)




#CellChat between B narrow and CD4 narrow identities, synovial, split between clonal and nonclonal CD4
clust_30PC_res03_seurat_BCD4_SYN <- subset(clust_30PC_res03_seurat_BCD4, subset = tissue_source == "SYN")

cd4_clonal_meta <- read.csv("~/Downloads/CD4_synovial_metadata_clonal.csv")
cd4_clonalidents <- cd4_clonal_meta %>% filter(clonal == "Yes")
cd4_nonclonalidents <- cd4_clonal_meta %>% filter(clonal == "No")

clust_30PC_res03_seurat_BCD4_SYN$clonalNonClonalB <- ifelse(clust_30PC_res03_seurat_BCD4_SYN$Barcode %in% cd4_clonalidents$Barcode, "Clonal", 
                                                            ifelse(clust_30PC_res03_seurat_BCD4_SYN$Barcode %in% cd4_nonclonalidents$Barcode, "Non-Clonal",
                                                            ifelse(colnames(clust_30PC_res03_seurat_BCD4_SYN) %in% bcellsmeta$seurat_key, "B", "Other")))

clust_30PC_res03_seurat_BallCD4clonal_SYN <- subset(clust_30PC_res03_seurat_BCD4_SYN, subset = clonalNonClonalB %in% c("B", "Clonal"))

Idents(clust_30PC_res03_seurat_BallCD4clonal_SYN) <- clust_30PC_res03_seurat_BallCD4clonal_SYN$ident
data.input <- GetAssayData(clust_30PC_res03_seurat_BallCD4clonal_SYN, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BallCD4clonal_SYN)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_narrowBnarrowclonalCD4_SYN <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_narrowBnarrowclonalCD4_SYN@DB <- CellChatDB
cellchat_narrowBnarrowclonalCD4_SYN <- subsetData(cellchat_narrowBnarrowclonalCD4_SYN) # This step is necessary even if using the whole database
cellchat_narrowBnarrowclonalCD4_SYN <- identifyOverExpressedGenes(cellchat_narrowBnarrowclonalCD4_SYN)
cellchat_narrowBnarrowclonalCD4_SYN <- identifyOverExpressedInteractions(cellchat_narrowBnarrowclonalCD4_SYN)

cellchat_narrowBnarrowclonalCD4_SYN <- computeCommunProb(cellchat_narrowBnarrowclonalCD4_SYN)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_narrowBnarrowclonalCD4_SYN <- filterCommunication(cellchat_narrowBnarrowclonalCD4_SYN, min.cells = 10)
cellchat_narrowBnarrowclonalCD4_SYN <- computeCommunProbPathway(cellchat_narrowBnarrowclonalCD4_SYN)

cellchat_narrowBnarrowclonalCD4_SYN
cellchat_narrowBnarrowclonalCD4_SYN <- aggregateNet(cellchat_narrowBnarrowclonalCD4_SYN)
cellchat_narrowBnarrowclonalCD4_SYN <- netAnalysis_computeCentrality(cellchat_narrowBnarrowclonalCD4_SYN, slot.name = "netP")


clust_30PC_res03_seurat_BallCD4nonclonal_SYN <- subset(clust_30PC_res03_seurat_BCD4_SYN, subset = clonalNonClonalB %in% c("B", "Non-Clonal"))

Idents(clust_30PC_res03_seurat_BallCD4nonclonal_SYN) <- clust_30PC_res03_seurat_BallCD4nonclonal_SYN$ident
data.input <- GetAssayData(clust_30PC_res03_seurat_BallCD4nonclonal_SYN, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_BallCD4nonclonal_SYN)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_narrowBnarrownonclonalCD4_SYN <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_narrowBnarrownonclonalCD4_SYN@DB <- CellChatDB
cellchat_narrowBnarrownonclonalCD4_SYN <- subsetData(cellchat_narrowBnarrownonclonalCD4_SYN) # This step is necessary even if using the whole database
cellchat_narrowBnarrownonclonalCD4_SYN <- identifyOverExpressedGenes(cellchat_narrowBnarrownonclonalCD4_SYN)
cellchat_narrowBnarrownonclonalCD4_SYN <- identifyOverExpressedInteractions(cellchat_narrowBnarrownonclonalCD4_SYN)

cellchat_narrowBnarrownonclonalCD4_SYN <- computeCommunProb(cellchat_narrowBnarrownonclonalCD4_SYN)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_narrowBnarrownonclonalCD4_SYN <- filterCommunication(cellchat_narrowBnarrownonclonalCD4_SYN, min.cells = 10)
cellchat_narrowBnarrownonclonalCD4_SYN <- computeCommunProbPathway(cellchat_narrowBnarrownonclonalCD4_SYN)

cellchat_narrowBnarrownonclonalCD4_SYN
cellchat_narrowBnarrownonclonalCD4_SYN <- aggregateNet(cellchat_narrowBnarrownonclonalCD4_SYN)
cellchat_narrowBnarrownonclonalCD4_SYN <- netAnalysis_computeCentrality(cellchat_narrowBnarrownonclonalCD4_SYN, slot.name = "netP")

object.list <- list(Clonal= cellchat_narrowBnarrowclonalCD4_SYN, Nonclonal = cellchat_narrowBnarrownonclonalCD4_SYN)
cellchat_narrowBnarrowCD4_clonalsplit <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)



#CellChat between CD8 broad and CD4 narrow identities
clust_30PC_res03_seurat_CD8CD4 <- subset(clust_30PC_res03_seurat, subset = lineage == "B", invert = T)
clust_30PC_res03_seurat_CD8CD4$ident_CD8broad <- ifelse(clust_30PC_res03_seurat_CD8CD4$lineage == "CD8", "CD8", clust_30PC_res03_seurat_CD8CD4$ident)

Idents(clust_30PC_res03_seurat_CD8CD4) <- clust_30PC_res03_seurat_CD8CD4$ident_CD8broad
data.input <- GetAssayData(clust_30PC_res03_seurat_CD8CD4, assay = "RNA", slot = "data") # normalized data matrix
data.input <- normalizeData(data.input, scale.factor = 10000, do.log = TRUE)
labels <- Idents(clust_30PC_res03_seurat_CD8CD4)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_broadCD8narrowCD4_all <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat_broadCD8narrowCD4_all@DB <- CellChatDB
cellchat_broadCD8narrowCD4_all <- subsetData(cellchat_broadCD8narrowCD4_all) # This step is necessary even if using the whole database
cellchat_broadCD8narrowCD4_all <- identifyOverExpressedGenes(cellchat_broadCD8narrowCD4_all)
cellchat_broadCD8narrowCD4_all <- identifyOverExpressedInteractions(cellchat_broadCD8narrowCD4_all)

cellchat_broadCD8narrowCD4_all <- computeCommunProb(cellchat_broadCD8narrowCD4_all)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_broadCD8narrowCD4_all <- filterCommunication(cellchat_broadCD8narrowCD4_all, min.cells = 10)
cellchat_broadCD8narrowCD4_all <- computeCommunProbPathway(cellchat_broadCD8narrowCD4_all)

cellchat_broadCD8narrowCD4_all
cellchat_broadCD8narrowCD4_all <- aggregateNet(cellchat_broadCD8narrowCD4_all)
cellchat_broadCD8narrowCD4_all <- netAnalysis_computeCentrality(cellchat_broadCD8narrowCD4_all, slot.name = "netP")

groupSize <- as.numeric(table(cellchat_broadCD8narrowCD4_all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_broadCD8narrowCD4_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_broadCD8narrowCD4_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

