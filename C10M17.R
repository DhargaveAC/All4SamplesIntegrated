options(future.global.maxSize = 2000000000000) # Increase object size for Seurat
memory.limit(size = 3000000000000) # Increase R global environment size efficiency


# Read the Seurat RDS files from DIEM Filtering
C10 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C10/C10')
#C19 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C19/C19')
M17 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M17/M17')
#M20 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M20/M20')


# Filter smaller libraries and low expressed genes
# For C10
genes <- apply(C10@assays$RNA@counts, 1,
               function(x){length(x[x>=1]) >= 100})
nuclei <- apply(C10@assays$RNA@counts, 2,
                function(x){length(x[x>=1]) >= 500})
C10seu <- C10[genes, nuclei]
dim(C10) - dim(C10seu)
# For C19
#genes <- apply(C19@assays$RNA@counts, 1,
#               function(x){length(x[x>=1]) >= 100})
#nuclei <- apply(C19@assays$RNA@counts, 2,
#                function(x){length(x[x>=1]) >= 500})
#C19seu <- C19[genes, nuclei]
#dim(C19) - dim(C19seu)
# For M17
genes <- apply(M17@assays$RNA@counts, 1,
               function(x){length(x[x>=1]) >= 100})
nuclei <- apply(M17@assays$RNA@counts, 2,
                function(x){length(x[x>=1]) >= 500})
M17seu <- M17[genes, nuclei]
dim(M17) - dim(M17seu)
# For M20
#genes <- apply(M20@assays$RNA@counts, 1,
#               function(x){length(x[x>=1]) >= 100})
#nuclei <- apply(M20@assays$RNA@counts, 2,
#                function(x){length(x[x>=1]) >= 500})
#M20seu <- M20[genes, nuclei]
#dim(M20) - dim(M20seu)


# Naming Batches
patient.names <- list.files("C:/Users/s4526789/Documents/snRNAseqData/Samples")
C10seu@meta.data$batch <- patient.names[1]
#C19seu@meta.data$batch <- patient.names[2]
M17seu@meta.data$batch <- patient.names[3]
#M20seu@meta.data$batch <- patient.names[4]


# Merging Objects Together
C10M17 <- merge(x = C10seu, y = M17seu, 
                add.cell.ids = c("HM_C10", "HM_M17"), project = "C10 + M17")
C10M17 <- FindVariableFeatures(C10M17, selection.method = "vst",
                                 nfeatures = nrow(C10M17))
species.list <- SplitObject(C10M17, split.by = "ident")
for (i in 1:length(species.list)) {
  species.list[[i]] <- FindVariableFeatures(species.list[[i]],
                                            selection.method = 'vst',
                                            nfeatures = nrow(species.list[[i]]))
  species.list[[i]] <- NormalizeData(species.list[[i]],
                                     scale.factor = nrow(species.list[[i]]))
}


# Integrating Objects Together
features <- SelectIntegrationFeatures(object.list = species.list,
                                      nfeatures = nrow(C10M17))
C10M17 <- IntegrateData(anchorset = FindIntegrationAnchors(
  object.list = species.list, anchor.features = features))
C10M17@project.name <- "C10 + M17"
levels(C10M17)
dim(C10M17)


# Normalisation and Clustering
C10M17 <- FindVariableFeatures(C10M17, selection.method = "vst",
                                 nfeatures = nrow(C10M17))
C10M17 <- ScaleData(C10M17, features = VariableFeatures(C10M17))
VlnPlot(C10M17,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2, group.by = "batch")
C10M17 <- RunPCA(C10M17, features = VariableFeatures(C10M17))
ElbowPlot(C10M17, ndims = 50, reduction = "pca")
C10M17 <- FindNeighbors(C10M17, dims = 1:30)
C10M17 <- FindClusters(C10M17, resolution = 0.25) # 7 clusters
DimPlot(C10M17, reduction = 'pca', label = TRUE) + NoLegend()
C10M17 <- RunUMAP(C10M17, dims = 1:30, min.dist = 0.1,
                    n.epochs = 500, reduction = 'pca')
DimPlot(C10M17, reduction = 'umap', group.by = 'seurat_clusters')
DimPlot(C10M17, reduction = 'umap', split.by = 'batch', ncol = 2)


# Finding Cluster Markers using MAST
p.markers <- FindAllMarkers(C10M17, min.pct = 0.25, 
                            logfc.threshold = 0.25, test.use = "MAST")
head(p.markers)
p.markers$p_val_adj <- p.adjust(p.markers$p_val,
                                method = "bonferroni", n = nrow(p.markers))
p.topmarkers <- p.markers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)
write.csv(p.topmarkers, 'C:/Users/s4526789/Documents/snRNAseqData/Analysis/C10M17_CluMarkers.csv') 

DoHeatmap(C10M17, features = p.topmarkers$gene) + NoLegend()


# Cluster Annotation using CIPR and BluePrint+ENCODE reference database
library(CIPR)
CIPR(p.markers, reference = c("blueprint", "hsrnaseq"),
     plot_ind = TRUE, global_results_obj = TRUE, global_plot_obj = TRUE)
DT::datatable(CIPR_top_results)
p.CIPRlabel <- data.frame(CIPR_top_results)
write.csv(p.markers, file = "C:/Users/s4526789/Documents/snRNAseqData/Analysis/CIPRlabels.csv")
p.CIPRlabel <- p.CIPRlabel %>% group_by(cluster) %>% 
  top_n(n = 1, wt = percent_pos_correlation)
C10M17[["CIPR.labels"]] <- p.CIPRlabel$long_name[
  match(C10M17[[]][["seurat_clusters"]], p.CIPRlabel$cluster)]
umap3 <- DimPlot(C10M17, reduction = "umap", group.by = "CIPR.labels")


# Cluster Annotations using SingleR with BluePrint+ENCODE database
library(SingleR)
bp.se <- BlueprintEncodeData()
bp.se <- as(bp.se, 'SingleCellExperiment')
C10M17.sce <- as.SingleCellExperiment(C10M17)
colData(C10M17.sce)$seurat_clusters <- as.character(Idents(C10M17))  # go from factor to character
head(colData(C10M17.sce))
annotation.cluster <- SingleR(test = C10M17.sce, method = 'cluster',
                              clusters = C10M17.sce$ident,
                              ref = bp.se, labels = bp.se$label.main, 
                              de.method='wilcox', genes = 'de')
table(annotation.cluster$labels)
library(pheatmap)
plotScoreHeatmap(annotation.cluster)
C10M17[["SingleR.labels"]] <- annotation.cluster$labels[
  match(C10M17[[]][["seurat_clusters"]], rownames(annotation.cluster))]
umap4 <- DimPlot(C10M17, reduction = "umap", group.by = "SingleR.labels")


# Cluster Annotations using scCATCH with CellMatch database
library(scCATCH)
clu_markers <- findmarkergenes(C10M17, species = 'Human', 
                               cluster = 'All', match_CellMatch = TRUE, cancer = NULL, 
                               cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05,
                               tissue = c('Muscle', 'Skeletal muscle', 'Spinal cord',
                                          'Adipose tissue', 'Bone', 'Cartilage',
                                          'Pluripotent stem cell', 'Epithelium',
                                          'Blood', 'Ligament', 'Plasma', 'Serum',
                                          'Dermis', 'Skin', 'Bone marrow'))
table(clu_markers$clu_markers$cluster)
clu_ann <- scCATCH(clu_markers$clu_markers, species = 'Human', cancer = NULL,
                   tissue = c('Muscle', 'Skeletal muscle', 'Spinal cord',
                              'Adipose tissue', 'Bone', 'Cartilage',
                              'Pluripotent stem cell', 'Epithelium',
                              'Blood', 'Ligament', 'Plasma', 'Serum',
                              'Dermis', 'Skin', 'Bone marrow'))
as.data.frame(clu_markers$clu_markers)
C10M17[["scCATCH.labels"]] <- clu_ann$cell_type[
  match(C10M17[[]][["seurat_clusters"]], clu_ann$cluster)]
umap5 <- DimPlot(C10M17, reduction = "umap", group.by = "scCATCH.labels")


multiplot(umap3, umap4, umap5, cols = 3)


# Cell Trajectory Analysis
library(slingshot)
all.sds <- slingshot(Embeddings(C10M17, "umap"), 
                     clusterLabels = C10M17$seurat_clusters) 
C10M17.sce$slingPseudotime_1 <- NULL
C10M17.sce <- slingshot(C10M17.sce, reducedDim = 'UMAP', 
                          clusterLabels = C10M17$seurat_clusters)
library(viridis) # colour palette
colors <- viridis(50, alpha = 1)
plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
lines(all.sds, lwd = 2, type = 'lineages', col = 'black')
plot(reducedDims(all.sds), 
     col = colors[cut(C10M17.sce$slingPseudotime_1, breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(C10M17.sce), lwd=2)
# Plotting Pseudotime Values for Each Lineage
pt <- slingPseudotime(all.sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(1, 2))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
  lines(all.sds, lwd = 2, col = 'black', type = 'lineages')
}


######
# GO Enrichment Analysis
# All Clusters analysed together to save the allRes file and pathway plot
allint.clumarkers <- p.markers
p.GOannotate <- list()
for(cluster in 0:(max(as.integer(C10M17$seurat_clusters))-1)){
  cluster.markers <- p.markers[p.markers$cluster==cluster,]
  gene.list <- cluster.markers$p_val_adj
  names(gene.list) <- rownames(cluster.markers)
  GOdata <- new("topGOdata", ontology = "BP", allGenes = gene.list,
                geneSelectionFun = function(x){x < 0.01},
                annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
  resultKS.w01 <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
  allRes <- GenTable(GOdata, 
                     weight = resultKS.w01,
                     orderBy = "weight", 
                     topNodes = 50)
  allRes <- allRes[allRes$weight <= 0.01,]
  data.frame(allRes)
  #library(Rgraphviz)
  showSigOfNodes(GOdata, score(resultKS.w01), firstSigNodes = 7)
  go.terms <- sapply(allRes$GO.ID, 
                     function(x){return(Term(GOTERM[[as.character(x)]]))})
  p.GOannotate[[paste0("Cluster", cluster)]] <- go.terms
}
p.GOannotate


# Up and Down Regulated Gene visualisation 
gene.up <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC > 0) 
            %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
gene.down <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC <= 0) 
              %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
up.features <- gene.up %>% filter(avg_logFC > 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
down.features <- gene.down %>% filter(avg_logFC <= 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
d1 <- DotPlot(C10M17, features = up.features$gene) + RotatedAxis() + 
  labs(x = "Up-Regulated Features", y = "Clusters")
d2 <- DotPlot(C10M17, features = down.features$gene) + RotatedAxis() + 
  labs(x = "Down-Regulated Features", y = "Clusters")

multiplot(d1, d2, cols = 2)


