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


setwd("/Rao_Lab/AMP_Phase_2/manuscript/analysis/")

load("20220319_AMP_initialclusterings.RData") #product of tcellprep



# Final subsetting innate T cells -----------------------------------------------------
rm(initialinnate)

Idents(clust_30PC_res03_seurat_Tcells) <- rownames(clust_30PC_res03_seurat_Tcells@meta.data)
finalinnate <- subset(clust_30PC_res03_seurat_Tcells, idents = c(innatecells_frominnate, innatecells_fromprolif))
Idents(finalinnate) <- finalinnate$seurat_clusters

finalinnate$tissue_source <- factor(finalinnate$tissue_source, levels = c("SYN", "PBL"))
finalinnate$sampleTissue_Pub <- factor(finalinnate$sampleTissue_Pub, levels = c("RA01_SYN", "RA01_PBL", "RA02_SYN", "RA02_PBL", "RA03_SYN", "RA03_PBL", "RA04_SYN", "RA04_PBL", "RA05_SYN", "RA05_PBL", "RA06_SYN", "RA06_PBL", "RA07_SYN", "RA07_PBL", "RA08_SYN", "RA09_SYN", "RA10_SYN", "RA10_PBL", "RA11_SYN", "RA11_PBL", "RA12_SYN", "RA12_PBL"))

finalinnate <- NormalizeData(finalinnate, normalization.method = "LogNormalize", scale.factor = 10000)
finalinnate <- FindVariableFeatures(finalinnate, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(finalinnate), 10)

all.genes <- rownames(finalinnate)
finalinnate <- ScaleData(finalinnate, features = all.genes)
finalinnate <- RunPCA(finalinnate, features = VariableFeatures(object = finalinnate))

dims_use = 1:10
finalinnate <- RunHarmony(finalinnate, group.by.vars=c("sample_ID_Pub", "tissue_source"))
finalinnate <- RunUMAP(object=finalinnate, reduction="harmony", dims=dims_use, verbose=FALSE)
finalinnate <- FindNeighbors(object=finalinnate, reduction="harmony", dims=dims_use, verbose=FALSE)
finalinnate <- FindClusters(object=finalinnate, resolution=0.4, verbose=FALSE)

innatecell_markers <- FindAllMarkers(object=finalinnate, only.pos = TRUE, logfc.threshold = 0.40, min.pct = 0.40)

finalinnate <- subset(finalinnate, idents = 7, invert = T)

finalinnate <- NormalizeData(finalinnate, normalization.method = "LogNormalize", scale.factor = 10000)
finalinnate <- FindVariableFeatures(finalinnate, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(finalinnate), 10)

all.genes <- rownames(finalinnate)
finalinnate <- ScaleData(finalinnate, features = all.genes)
finalinnate <- RunPCA(finalinnate, features = VariableFeatures(object = finalinnate))

dims_use = 1:10
finalinnate <- RunHarmony(finalinnate, group.by.vars=c("sample_ID_Pub", "tissue_source"))
finalinnate <- RunUMAP(object=finalinnate, reduction="harmony", dims=dims_use, verbose=FALSE)
finalinnate <- FindNeighbors(object=finalinnate, reduction="harmony", dims=dims_use, verbose=FALSE)
finalinnate <- FindClusters(object=finalinnate, resolution=0.4, verbose=FALSE)

innatecell_markers <- FindAllMarkers(object=finalinnate, only.pos = TRUE, logfc.threshold = 0.40, min.pct = 0.40)

innatetcellcount_add <- finalinnate@meta.data %>% group_by(sampleTissue_Pub) %>% add_tally(name = "sample_innatetcellcount")
finalinnate$innatetcellCount <- innatetcellcount_add$sample_innatetcellcount


Idents(finalinnate) <- finalinnate$seurat_clusters
current.cluster.ids <- c(0:6)
new.cluster.ids <- c("Vdelta1", "CD56dim NK", "Vdelta2", "MAIT", "ZNF683+ cells", "MT-high", "CD56bright NK")
finalinnate@active.ident <- plyr::mapvalues(x = finalinnate@active.ident, from = current.cluster.ids, to = new.cluster.ids)
finalinnate$clusters_names <- finalinnate@active.ident

finalinnate$clusters_names <- factor(finalinnate$clusters_names, levels = c("Vdelta1", "Vdelta2", "MAIT", "CD56dim NK", "CD56bright NK", "ZNF683+ cells", "MT-high"))
Idents(finalinnate) <- finalinnate$clusters_names

innateTcell_markers <- FindAllMarkers(object=finalinnate, only.pos = FALSE, logfc.threshold = 0.50, min.pct = 0.40)
write.csv(innateTcell_markers, "finalinnatetcellclustermarkers_res0.4.csv")

innatecols <- c(colorRampPalette(brewer.pal(8, "Paired"))(8))

DimPlot(finalinnate) + 
  scale_color_manual(values = innatecols, name = "Innate T Cell\nClusters") +
  ggtitle(NULL) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.text = element_text(size=10),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_text(size=10, face="bold"))
ggsave("innate_UMAP_namedclusters.eps", width = 5, height = 3)

StackedVlnPlot(finalinnate, features = c("CD3E", "CD8A", "CD4", "TRDV1", "TIGIT", 
                                         "TRDV2", "TRGV9", "SLC4A10", "AQP3", "ZBTB16", 
                                         "FCER1G", "KLRF1", "NKG7", "TYROBP", "PRF1", "GZMB",
                                         "GZMH", "XCL1",  "GZMK", "ZNF683", "MALAT1"), cols = innatecols)
ggsave("innate_stackedVln_markers_namedclusters.eps", width = 3.5, height = 5)

DotPlot(finalinnate, features = rev(c("CD3E", "CD8A", "CD4", "TRDV1", "TIGIT", 
                                   "TRDV2", "TRGV9", "SLC4A10", "AQP3", "ZBTB16", 
                                   "FCER1G", "KLRF1", "NKG7", "TYROBP", "PRF1", "GZMB",
                                   "GZMH", "XCL1",  "GZMK", "ZNF683", "MALAT1")), group.by = "clusters_names") + 
  coord_flip() + 
  RotatedAxis() + 
  xlab(NULL) + ylab(NULL) + 
  scale_color_distiller(palette = "Reds", direction = 1) +
  theme(legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        axis.text = element_text(size = 10))
ggsave("innate_dotplot_markers_namedclusters.eps", width = 4.5, height = 5.5)

FeaturePlot(finalinnate, features = c("TRDV1", "TRDV2", "SLC4A10", "FCER1G")) &
  scale_color_distiller(palette = "Reds", direction = 1) &
  theme(plot.title = element_text(size = 11),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("innate_TRDV1_TRDV2_SLC4A10_FCER1G.eps", width = 5, height = 4)


#loading in main pipeline Phase 2 T and NK cell data
tcell_RNA <- readRDS("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/data/mainPipeline/mRNA_counts_tcell_pass_QC.rds")
tcell_ADT <- readRDS("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/data/mainPipeline/protein_counts_tcell_pass_QC.rds")
tcell_meta <- readRDS("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/data/mainPipeline/meta_tcell_pass_QC.rds")
rownames(tcell_meta) <- tcell_meta$cell

AMP_phase2main_T <- CreateSeuratObject(tcell_RNA)
AMP_phase2main_T[["ADT"]] <- CreateAssayObject(tcell_ADT)

AMP_phase2main_T <- AddMetaData(AMP_phase2main_T, tcell_meta)
Idents(AMP_phase2main_T) <- AMP_phase2main_T$new_cluster_number

AMP_phase2main_T_innate <- subset(AMP_phase2main_T, idents = c(18,19,21,22,23))

nkcell_RNA <- readRDS("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/data/mainPipeline/mRNA_counts_nk_pass_QC.rds")
nkcell_ADT <- readRDS("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/data/mainPipeline/protein_counts_nk_pass_QC.rds")
nkcell_meta <- readRDS("/Users/gdunlap/Desktop/Rao_Lab/AMP_Phase_2/data/mainPipeline/meta_nk_pass_QC.rds")
rownames(nkcell_meta) <- nkcell_meta$cell

AMP_phase2main_NK <- CreateSeuratObject(nkcell_RNA)
AMP_phase2main_NK[["ADT"]] <- CreateAssayObject(nkcell_ADT)

AMP_phase2main_NK <- AddMetaData(AMP_phase2main_NK, nkcell_meta)
Idents(AMP_phase2main_NK) <- AMP_phase2main_NK$new_cluster_number

AMP_phase2main_T_NK <- merge(AMP_phase2main_T_innate, y = AMP_phase2main_NK)

AMP_phase2main_T_NK <- NormalizeData(AMP_phase2main_T_NK, normalization.method = "LogNormalize", scale.factor = 10000)

AMP_phase2main_T_NK$new_cluster_name <- factor(AMP_phase2main_T_NK$new_cluster_name, levels = c("T-22: Vdelta1", "T-23: Vdelta2", "T-21: Innate-like", "NK-0: CD56dim CD16+ IFNG-", "NK-1: CD56dim CD16+ IFNG+CD160+", "NK-2: CD56dim CD16+ IFNG+CD160-", "NK-3: CD56dim CD16+ GZMB-", "NK-6: CD56bright CD16- GNLY+", "NK-7: CD56bright CD16- GNLY+CD69+", "NK-4: CD56bright CD16- GZMA+CD160+", "NK-5: CD56bright CD16- GZMA+CD69+",  "NK-8: CD56bright CD16- IFN response", "NK-12: IL7R+ ILC", "NK-13: IL7R+CD161+ ILC", "NK-10: PCNA+ Proliferating", "NK-11: MKI67+ Proliferating", "T-18: Proliferating", "NK-9: MT-high", "T-19: MT-high (low quality)"))
Idents(AMP_phase2main_T_NK) <- AMP_phase2main_T_NK$new_cluster_name

#Symphony on repertoire innage T and NK cell data
symphT_reference <- symphony::buildReference(
  AMP_phase2main_T_NK@assays$RNA@data,
  AMP_phase2main_T_NK@meta.data,
  vars = c('sample'),         # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = FALSE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'new_cluster_name', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20                     # number of PCs
)

queryT <- mapQuery(finalinnate@assays$RNA@data,             # query gene expression (genes x cells)
                   finalinnate@meta.data,        # query metadata (cells x attributes)
                   symphT_reference,             # Symphony reference object
                   vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                   do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                   do_umap = FALSE)        # project query cells into reference UMAP

queryT <- knnPredict(queryT, symphT_reference, train_labels = symphT_reference$meta_data$new_cluster_name,
                     k = 10, save_as = 'cell_type_10_nn', confidence = TRUE)

finalinnate$symphony_narrow_type <- queryT$meta_data$cell_type_10_nn

## Cell type mappings heatmap
query_immune = queryT$meta_data
query_immune$symphony_narrow_type = droplevels(as.factor(query_immune$cell_type_10_nn))
res_immune = symphony:::evaluate(query_immune$clusters_names, query_immune$cell_type_10_nn)
Conf_immune = res_immune$Conf / rowSums(res_immune$Conf)

dev.off()
setEPS(width = 8.5, height = 7.5)
postscript("innate_symphony_mainPipeline.eps")
pheatmap::pheatmap(as.matrix(Conf_immune), cluster_rows = F, cluster_cols = F, angle_col = c("45"), cellwidth = 20, cellheight = 20, color=colorRampPalette(c("white", "#962f23", "#5e0e04"))(1000), border_color = "white")
dev.off()

#pseudobulk heatmap
col_cor = colorRamp2(c(1, 0), c("#5e0e04", "white"))
pseudobulk_byinnatecluster <- AverageExpression(finalinnate, group.by = "clusters_names")
cor.exp <- cor(pseudobulk_byinnatecluster$RNA)
Heatmap(as.matrix(cor.exp), col = col_cor)

#synovial vs blood skewing
#all clusters
innatemetadata <- finalinnate@meta.data

innateclusterbytissue <- innatemetadata %>% 
  group_by(clusters_names, sampleTissue_Pub) %>% 
  add_tally(name = "patientinnateclustersize") %>% 
  ungroup() %>%
  select(sample_ID_Pub, tissue_source, clusters_names, totalcellCount, tcellCount, innatetcellCount, patientinnateclustersize) %>%
  distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all = T) %>%
  mutate(freqinnate = patientinnateclustersize / innatetcellCount) %>%
  mutate(freqT = patientinnateclustersize / tcellCount)

#by innate T cells
stat.test_freqinnate <- innateclusterbytissue %>%
  group_by(clusters_names) %>%
  pairwise_t_test(freqinnate ~ tissue_source) %>%
  add_significance()

stat.test_freqinnate <- stat.test_freqinnate %>% add_xy_position(x = "tissue_source")

ggplot(innateclusterbytissue, aes(x = tissue_source, y = freqinnate)) +
  geom_line(aes(group = sample_ID_Pub)) +
  geom_point(aes(color = tissue_source)) + facet_wrap(~ clusters_names) + theme_classic() + ylim(0,1) +
  scale_color_manual(values = c("#f2a21e", "#88a0c4")) + theme(legend.position = "none") + xlab(NULL) +
  ylab("Frequency of subcluster in all innate cells, per patient") + stat_pvalue_manual(stat.test_freqinnate)
ggsave("innate_clusterfrequencyInnatedenom_byTissue.eps", width = 6, height = 5)

#by all T cells
stat.test_freqallT <- innateclusterbytissue %>%
  group_by(clusters_names) %>%
  pairwise_t_test(freqT ~ tissue_source) %>%
  add_significance()

stat.test_freqallT <- stat.test_freqallT %>% add_xy_position(x = "tissue_source")

ggplot(innateclusterbytissue, aes(x = tissue_source, y = freqT)) +
  geom_line(aes(group = sample_ID_Pub)) +
  geom_point(aes(color = tissue_source)) + facet_wrap(~clusters_names) + theme_classic() + ylim(0,0.20) +
  scale_color_manual(values = c("#f2a21e", "#88a0c4")) + theme(legend.position = "none") + xlab(NULL) +
  ylab("Frequency of subcluster in all T cells, per patient") + stat_pvalue_manual(stat.test_freqallT)
ggsave("innate_clusterfrequencyTcelldenom_byTissue.eps", width = 6, height = 5)

#select clusters
innateclusterbytissue <- innatemetadata %>% 
  filter(clusters_names %in% c("Vdelta1", "Vdelta2", "MAIT", "ZNF683+ cells")) %>%
  group_by(clusters_names, sampleTissue_Pub) %>% 
  add_tally(name = "patientinnateclustersize") %>% 
  ungroup() %>%
  select(sample_ID_Pub, tissue_source, clusters_names, totalcellCount, tcellCount, innatetcellCount, patientinnateclustersize) %>%
  distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all = T) %>%
  mutate(freqinnate = patientinnateclustersize / innatetcellCount) %>%
  mutate(freqT = patientinnateclustersize / tcellCount)

#by innate T cells
stat.test_freqinnate <- innateclusterbytissue %>%
  group_by(clusters_names) %>%
  pairwise_t_test(freqinnate ~ tissue_source) %>%
  add_significance()

stat.test_freqinnate <- stat.test_freqinnate %>% add_xy_position(x = "tissue_source")

ggplot(innateclusterbytissue, aes(x = tissue_source, y = freqinnate)) +
  geom_line(aes(group = sample_ID_Pub)) +
  geom_point(aes(color = tissue_source)) + facet_wrap(~ clusters_names) + theme_classic() + ylim(0,1) +
  scale_color_manual(values = c("#f2a21e", "#88a0c4")) + theme(legend.position = "none") + xlab(NULL) +
  ylab("Frequency of subcluster in\nall innate cells, per patient") + stat_pvalue_manual(stat.test_freqinnate)
ggsave("innate_clusterfrequencyInnatedenom_byTissue_selectClusters.eps", width = 4, height = 3)

#by all T cells
stat.test_freqallT <- innateclusterbytissue %>%
  group_by(clusters_names) %>%
  pairwise_t_test(freqT ~ tissue_source) %>%
  add_significance()

stat.test_freqallT <- stat.test_freqallT %>% add_xy_position(x = "tissue_source")

ggplot(innateclusterbytissue, aes(x = tissue_source, y = freqT)) +
  geom_line(aes(group = sample_ID_Pub)) +
  geom_point(aes(color = tissue_source)) + facet_wrap(~clusters_names) + theme_classic() + ylim(0,0.1) +
  scale_color_manual(values = c("#f2a21e", "#88a0c4")) + theme(legend.position = "none") + xlab(NULL) +
  ylab("Frequency of subcluster in\nall T cells, per patient") + stat_pvalue_manual(stat.test_freqallT, label = "p.adj")
ggsave("innate_clusterfrequencyTcelldenom_byTissue_selectClusters.eps", width = 4, height = 3)


#distribution by cluster and tissue source
DimPlot(finalinnate, group.by = "tissue_source", shuffle = T) +
  scale_color_manual(values = c("#f2a21e", "#88a0c4"), name = NULL, labels = c("SYN", "PBL")) +
  ggtitle(NULL) + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(legend.title = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())
ggsave("innate_UMAP_byTissue.eps", width = 4, height = 2.5)

finalinnate@meta.data %>%
  group_by(clusters_names, tissue_source) %>%
  dplyr::count() %>%
  group_by(clusters_names) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=tissue_source)) +
  geom_col() +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = "Source", labels = c("SYN", "PBL")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm')) 
ggsave("innate_barplot_tissueBreakdown_percent.eps", width = 4, height = 3)

#clonality by cluster

finalinnate$clonal <- replace_na(finalinnate$clonal, "No TCR")

finalinnate@meta.data %>%
  group_by(clusters_names, cloneType) %>%
  dplyr::count() %>%
  group_by(clusters_names) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=cloneType)) +
  geom_col() +
  scale_fill_brewer(palette = "Blues", direction = -1, name = "Clone Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm')) 
ggsave("innate_barplot_clonality_percent.eps", width = 5, height = 3)

DimPlot(finalinnate, group.by = "cloneType", shuffle = T) + 
  scale_color_brewer(palette = "Blues", direction = -1) +
  ggtitle(NULL) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.title = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.35, 'cm'))
ggsave("innate_UMAP_clonality.eps", width = 5.5, height = 3)

finalinnate@meta.data %>%
  group_by(clusters_names, cloneType) %>%
  dplyr::count() %>%
  group_by(cloneType) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=cloneType,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = innatecols) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=9),
        legend.title = element_blank(),
        legend.key.size = unit(0.35, 'cm'),
        plot.margin = margin(10, 10, 10, 50)) 
ggsave("innate_barplot_cloneType_clusterBreakdown.eps", width = 5, height = 4)

finalinnate@meta.data %>%
  group_by(clusters_names, clonal) %>%
  dplyr::count() %>%
  group_by(clonal) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clonal,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = innatecols) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(limits = c("Yes", "No", "No TCR"), labels = c("Clonal", "Not Clonal", "No TCR")) +
  ggtitle(NULL) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=9),
        legend.title = element_blank(),
        legend.key.size = unit(0.35, 'cm')) 
ggsave("innate_barplot_clonal_clusterBreakdown.eps", width = 3.5, height = 3)

#distribution by sample ID and tissue
finalinnate@meta.data %>%
  group_by(sampleTissue_Pub, clusters_names) %>%
  dplyr::count() %>%
  group_by(sampleTissue_Pub) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=sampleTissue_Pub,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = innatecols, name = "Cluster") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm')) 
ggsave("innate_barplot_patientSampleBreakdown_percent.eps", width = 6, height = 3)

finalinnate@meta.data %>%
  group_by(sampleTissue_Pub, clusters_names) %>%
  dplyr::count() %>%
  group_by(sampleTissue_Pub) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=sampleTissue_Pub,y=n, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = innatecols, name = "Cluster") +
  ggtitle(NULL) +
  ylab("Number of Cells") +
  xlab(NULL) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm')) 
ggsave("innate_barplot_patientSampleBreakdown_number.eps", width = 6.1, height = 3)

innate_p2 <- finalinnate@meta.data %>%
  group_by(clusters_names, tissue_source) %>%
  dplyr::count() %>%
  group_by(clusters_names) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=tissue_source)) +
  geom_col() +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = "Source", labels = c("SYN", "PBL")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm')) 

innate_p1 <- finalinnate@meta.data %>%
  group_by(clusters_names, sort_group) %>%
  dplyr::count() %>%
  group_by(clusters_names) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=n, fill = clusters_names)) +
  geom_col() +
  scale_fill_manual(name = NULL, values = innatecols) +
  scale_y_continuous(limits = c(0,1000)) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10)) 

innateplotalign <- align_plots(innate_p1, innate_p2, align="hv", axis="tblr")
innate_p1x <- ggdraw(innateplotalign[[1]])
innate_p2x <- ggdraw(innateplotalign[[2]])

save_plot("innate_barplot_numberCells_aligned.eps", innate_p1x, base_width = 4.8, base_height = 2)
save_plot("innate_barplot_TissueDist_aligned.eps", innate_p2x, base_width = 4.8, base_height = 2)


#examining invariant MAIT TCRs
finalinnate_chainsplit <- finalinnate

finalinnate_chainsplit@meta.data <- separate(data = finalinnate_chainsplit@meta.data, col = CTgene, into = c("TRAgene", "TRBgene"), sep = "_")
finalinnate_chainsplit@meta.data <- separate(data = finalinnate_chainsplit@meta.data, col = TRAgene, into = c("TRAV", "TRAJ", "TRAC"), sep = "\\.")
finalinnate_chainsplit@meta.data <- separate(data = finalinnate_chainsplit@meta.data, col = TRBgene, into = c("TRBV", "TRBJ", "TRBD", "TRBC"), sep = "\\.")

Idents(finalinnate_chainsplit) <- finalinnate_chainsplit$TRAV
DimPlot(finalinnate_chainsplit, cells.highlight = WhichCells(finalinnate_chainsplit, idents = "TRAV1-2"), sizes.highlight = 0.5) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("TRAV1-2", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("TRAV1-2") 
ggsave("innate_UMAP_TRAV1-2.eps", width = 3, height = 2.5)

Idents(finalinnate_chainsplit) <- finalinnate_chainsplit$TRAJ
DimPlot(finalinnate_chainsplit, cells.highlight = WhichCells(finalinnate_chainsplit, idents = "TRAJ33"), sizes.highlight = 0.5) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("TRAJ33", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("TRAJ33") 
ggsave("innate_UMAP_TRAJ33.eps", width = 3, height = 2.5)

DimPlot(finalinnate_chainsplit, cells.highlight = WhichCells(finalinnate_chainsplit, idents = "TRAJ12"), sizes.highlight = 0.5) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("TRAJ12", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("TRAJ12") 
ggsave("innate_UMAP_TRAJ12.eps", width = 3, height = 2.5)

DimPlot(finalinnate_chainsplit, cells.highlight = WhichCells(finalinnate_chainsplit, idents = "TRAJ20"), sizes.highlight = 0.5) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("TRAJ20", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("TRAJ20") 
ggsave("innate_UMAP_TRAJ20.eps", width = 3, height = 2.5)

finalinnate_chainsplit$TRAV_TRAJ <- paste0(finalinnate_chainsplit$TRAV, "_", finalinnate_chainsplit$TRAJ)
Idents(finalinnate_chainsplit) <- finalinnate_chainsplit$TRAV_TRAJ
DimPlot(finalinnate_chainsplit, cells.highlight = WhichCells(finalinnate_chainsplit, idents = c("TRAV1-2_TRAJ33","TRAV1-2_TRAJ12","TRAV1-2_TRAJ20")), sizes.highlight = 0.5) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("TRAV_TRAJ", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("TRAV1-2 - TRAJ33/20/12 TCR") 
ggsave("innate_UMAP_TRAV1-2-TRAJ33_20_12.eps", width = 3, height = 2.5)


finalinnate_chainsplit$MAIT_tra <- ifelse(finalinnate_chainsplit$TRAV == "TRAV1-2" & (finalinnate_chainsplit$TRAJ == "TRAJ33" | finalinnate_chainsplit$TRAJ == "TRAJ12" | finalinnate_chainsplit$TRAJ == "TRAJ20" ), "TRAV1-2 \u2013 TRAJ33/20/12", "Other")
finalinnate_chainsplit$MAIT_tra <- replace_na(finalinnate_chainsplit$MAIT_tra, "No TCR")

finalinnate_chainsplit@meta.data %>%
  group_by(clusters_names, MAIT_tra) %>%
  dplyr::count() %>%
  group_by(clusters_names) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=MAIT_tra)) +
  geom_col() +
  scale_fill_brewer(palette = "Blues", breaks = c("TRAV1-2 \u2013 TRAJ33/20/12", "Other", "No TCR"), labels = c("TRAV1-2 & TRAJ33/20/12", "Other TRA Pairing", "No TCR")) +
  ggtitle(NULL) +
  xlab(NULL) +
  labs(fill = NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#MAIT overlap
mait_syn <- finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(tissue_source == "SYN")
mait_pbl <- finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(tissue_source == "PBL")

Idents(finalinnate) <- finalinnate$CTaa
DimPlot(finalinnate, cells.highlight = WhichCells(finalinnate, idents = na.omit(intersect(mait_pbl$CTaa, mait_syn$CTaa))), sizes.highlight = 0.5) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Shared", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  ggtitle("MAIT Clones Shared SYN and PBL") 
ggsave("innate_UMAP_MAIToverlap.eps", width = 3, height = 2.5)
Idents(finalinnate) <- finalinnate$clusters_names

mait_all <- finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% na.omit(CTaa) %>% mutate(patientCTaa = paste0(sample_ID_Pub, CTaa))
mait_all_clones_tissue <- table(mait_all$patientCTaa, mait_all$tissue_source)
as.data.frame(mait_all_clones_tissue) %>% pivot_wider(names_from = Var2, values_from = Freq) %>%
  mutate(Patient = substr(Var1, 1, 4)) %>% 
  ggplot(aes(x = PBL, y = SYN, color = Patient)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.2), size = 2) + 
  theme_classic() + scale_color_brewer(palette = "Set3") +
  guides(color=guide_legend(ncol=2, override.aes = list(size=3))) + 
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6)) +
  xlab("Blood cells in clone") + ylab("Synovial cells in clone") +
  theme(legend.text = element_text(size=10),
        legend.key.size = unit(0.50, 'cm'),
        legend.title = element_text(size = 11))
ggsave("innate_mait_synpbl_share.eps", width = 4.7, height = 3)


Idents(finalinnate) <- finalinnate$clusters_names
finalinnate_mait <- subset(finalinnate, idents = "MAIT")
Idents(finalinnate_mait) <- finalinnate_mait$clonal
finalinnate_mait <- subset(finalinnate_mait, idents = "No TCR", invert = T)
finalinnate_mait$tissue_clonal <- paste0(finalinnate_mait$tissue_source, "_", finalinnate_mait$clonal)

finalinnate_mait_clonal <- subset(finalinnate_mait, idents = "Yes")
Idents(finalinnate_mait_clonal) <- finalinnate_mait_clonal$tissue_source

mait_clonal_tissueDEGs <- FindMarkers(finalinnate_mait_clonal, ident.1 = "SYN", ident.2 = "PBL", )
mait_clonal_tissueDEGs_cutoff <- mait_clonal_tissueDEGs %>% filter(p_val_adj < 0.01)

Idents(finalinnate) <- finalinnate$clusters_names

cd4signatures <- read.csv("~/Desktop/Rao_Lab/AMP_Phase_2/manuscript/cd4_signatures.csv")
finalinnate <- AddModuleScore(finalinnate, list(cd4signatures$Activation), name = "Activation_sig", ctrl = 1000)
maitstim <- read.csv("~/Downloads/mait_stim_genes.csv")
finalinnate <- AddModuleScore(finalinnate, list(maitstim$TCR), name = "MAIT_tcrstim1_sig", ctrl = 1000)
finalinnate <- AddModuleScore(finalinnate, list(maitstim$TCR2), name = "MAIT_tcrstim2_sig", ctrl = 1000)
finalinnate <- AddModuleScore(finalinnate, list(maitstim$Cytokine), name = "MAIT_cytostim1_sig", ctrl = 1000)
finalinnate <- AddModuleScore(finalinnate, list(maitstim$Cytokine2), name = "MAIT_cytostim2_sig", ctrl = 1000)

leng_tcr_noOver <- setdiff(maitstim$TCR, maitstim$Cytokine)
leng_cyto_noOver <- setdiff(maitstim$Cytokine, maitstim$TCR)

lam_tcr_noOver <- setdiff(maitstim$TCR2, maitstim$Cytokine2)
lam_cyto_noOver <- setdiff(maitstim$Cytokine2, maitstim$TCR2)

finalinnate <- AddModuleScore(finalinnate, list(leng_tcr_noOver), name = "MAIT_tcrstim1_noOverlap_sig", ctrl = 1000)
finalinnate <- AddModuleScore(finalinnate, list(lam_tcr_noOver), name = "MAIT_tcrstim2_noOverlap_sig", ctrl = 1000)
finalinnate <- AddModuleScore(finalinnate, list(leng_cyto_noOver), name = "MAIT_cytostim1_noOverlap_sig", ctrl = 1000)
finalinnate <- AddModuleScore(finalinnate, list(lam_cyto_noOver), name = "MAIT_cytostim2_noOverlap_sig", ctrl = 1000)

a <- finalinnate@meta.data %>% 
  filter(clusters_names %in% c("Vdelta1", "Vdelta2", "MAIT", "ZNF683+ cells")) %>% 
  group_by(sample_ID_Pub, tissue_source, clusters_names) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_source) %>% 
  pairwise_t_test(mean ~ tissue_source) %>% add_significance() %>% add_y_position() %>% mutate(y.position = 1*y.position)

finalinnate@meta.data %>% filter(clusters_names %in% c("Vdelta1", "Vdelta2", "MAIT", "ZNF683+ cells")) %>% 
  group_by(sample_ID_Pub, tissue_source, clusters_names) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_source, y = mean)) + 
  geom_boxplot(aes(fill = tissue_source)) + geom_point() + facet_wrap(~clusters_names) +
  stat_pvalue_manual(a, hide.ns = T) +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("Activation Signature") +
  theme_classic() + scale_y_continuous(limits = c(0.18, 0.32), breaks = c(0.2, 0.3)) +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") 
ggsave("innate_gdT_mait_znf_activationsig_SYNPBL.eps", width = 2.8, height = 3)

a <- finalinnate@meta.data %>% 
  filter(clusters_names %in% c("Vdelta1", "Vdelta2", "MAIT")) %>% 
  group_by(sample_ID_Pub, tissue_source, clusters_names) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_source) %>% 
  pairwise_t_test(mean ~ tissue_source) %>% add_significance() %>% add_y_position() %>% mutate(y.position = 1*y.position)

finalinnate@meta.data %>% filter(clusters_names %in% c("Vdelta1", "Vdelta2", "MAIT")) %>% 
  group_by(sample_ID_Pub, tissue_source, clusters_names) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_source, y = mean)) + 
  geom_boxplot(aes(fill = tissue_source)) + geom_point() + facet_wrap(~clusters_names) +
  stat_pvalue_manual(a, hide.ns = T) +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("Activation Signature") +
  theme_classic() + scale_y_continuous(limits = c(0.18, 0.32), breaks = c(0.2, 0.3)) +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") 
ggsave("innate_gdT_mait_activationsig_SYNPBL.eps", width = 4, height = 2)

a <- finalinnate@meta.data %>% 
  filter(clusters_names %in% c("Vdelta1", "Vdelta2", "ZNF683+ cells")) %>% 
  group_by(sample_ID_Pub, tissue_source, clusters_names) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_source) %>% 
  pairwise_t_test(mean ~ tissue_source) %>% add_significance() %>% add_y_position() %>% mutate(y.position = 1*y.position)

finalinnate@meta.data %>% filter(clusters_names %in% c("Vdelta1", "Vdelta2", "ZNF683+ cells")) %>% 
  group_by(sample_ID_Pub, tissue_source, clusters_names) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, clusters_names, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_source, y = mean)) + 
  geom_boxplot(aes(fill = tissue_source)) + geom_point() + facet_wrap(~clusters_names) +
  stat_pvalue_manual(a, hide.ns = T) +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("Activation Signature") +
  theme_classic() + scale_y_continuous(limits = c(0.18, 0.32), breaks = c(0.2, 0.3)) +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") 
ggsave("innate_gdT_znf_activationsig_SYNPBL.eps", width = 4, height = 2)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  group_by(sample_ID_Pub, tissue_source) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_source) %>% 
  pairwise_t_test(mean ~ tissue_source) %>% add_significance() %>% add_y_position() %>% mutate(y.position = 1.01*y.position)

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% 
  group_by(sample_ID_Pub, tissue_source) %>%
  mutate(mean = mean(Activation_sig1)) %>% distinct(sample_ID_Pub, tissue_source, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_source, y = mean)) + 
  geom_boxplot(aes(fill = tissue_source)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_fill_manual(values = c("#f2a21e", "#88a0c4"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("Activation Signature") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") 
ggsave("innate_mait_activationsig_SYNPBL.eps", width = 2.5, height = 2)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim1_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim1_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Dependent MAIT \nActivation (Leng et al, 2019)") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_tcractivationsig_SYNPBL_clonal_leng.eps", width = 3.3, height = 2.5)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim1_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim1_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Independent MAIT \nActivation (Leng et al, 2019)") +
  theme_classic() + scale_y_continuous(limits = c(0.1, 0.35), breaks = c(0.1, 0.2, 0.3)) +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_cytoactivationsig_SYNPBL_clonal_leng.eps", width = 3.3, height = 2.5)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim2_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim2_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Dependent MAIT \nActivation (Lamichhane et al, 2020)") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_tcractivationsig_SYNPBL_clonal_lamichhane.eps", width = 3.3, height = 2.7)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim2_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim2_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Independent MAIT \nActivation (Lamichhane et al, 2020)") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_cytoactivationsig_SYNPBL_clonal_lamichhane.eps", width = 3.3, height = 2.7)


a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim1_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim1_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Dependent MAIT \nActivation (Leng et al, 2019)") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_tcractivationsig_SYNPBL_clonal_leng_noOverlap.eps", width = 3.3, height = 2.5)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim1_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim1_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Independent MAIT \nActivation (Leng et al, 2019)") +
  theme_classic() + scale_y_continuous(limits = c(0.1, 0.4), breaks = c(0.1, 0.2, 0.3, 0.4)) +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_cytoactivationsig_SYNPBL_clonal_leng_noOverlap.eps", width = 3.3, height = 2.5)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim2_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_tcrstim2_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Dependent MAIT \nActivation (Lamichhane et al, 2020)") +
  theme_classic() +
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_tcractivationsig_SYNPBL_clonal_lamichhane_noOverlap.eps", width = 3.3, height = 2.7)

a <- finalinnate@meta.data %>% 
  filter(clusters_names == "MAIT") %>% 
  filter(clonal != "No TCR") %>% 
  mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim2_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ungroup(sample_ID_Pub, tissue_clonal) %>% 
  pairwise_wilcox_test(mean ~ tissue_clonal, p.adjust.method = "holm") %>% add_significance() %>% add_y_position() 

finalinnate@meta.data %>% filter(clusters_names == "MAIT") %>% filter(clonal != "No TCR") %>% mutate(tissue_clonal = paste0(tissue_source, "_", clonal)) %>%
  group_by(sample_ID_Pub, tissue_clonal) %>%
  mutate(mean = mean(MAIT_cytostim2_noOverlap_sig1)) %>% distinct(sample_ID_Pub, tissue_clonal, .keep_all=T) %>% select(mean) %>%
  ggplot(aes(x= tissue_clonal, y = mean)) + 
  geom_boxplot(aes(fill = tissue_clonal)) + geom_point() + 
  stat_pvalue_manual(a, hide.ns = T) +
  scale_x_discrete(limits = c("SYN_Yes", "SYN_No", "PBL_Yes", "PBL_No"), labels = c("SYN (>1)", "SYN (=1)", "PBL (>1)", "PBL (=1)")) +
  scale_fill_manual(values = c("#88a0c4", "#88a0c4", "#f2a21e", "#f2a21e"), name = NULL) +
  ggtitle(NULL) + xlab(NULL) + ylab("TCR-Independent MAIT \nActivation (Lamichhane et al, 2020)") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 10),
        legend.position = "none") + RotatedAxis()
ggsave("innate_mait_cytoactivationsig_SYNPBL_clonal_lamichhane_noOverlap.eps", width = 3.3, height = 2.7)



#TRDV2 syn vs pbl DEGs
finalinnate_trdv2 <- subset(finalinnate, idents = "Vdelta2")
Idents(finalinnate_trdv2) <- finalinnate_trdv2$tissue_source

Vdelta2_tissue_DEGs <- FindMarkers(finalinnate_trdv2, ident.1 = "SYN", ident.2 = "PBL")

EnhancedVolcano::EnhancedVolcano(Vdelta2_tissue_DEGs, lab = as.character(rownames(Vdelta2_tissue_DEGs)), x = 'avg_log2FC', y = 'p_val_adj', title = "Vdelta2 SYN/PBL DEGs", subtitle = NULL, legendPosition = 'none', legendLabels = NULL, caption = NULL, drawConnectors = TRUE, arrowheads = FALSE, labSize = 3) + theme_classic()
