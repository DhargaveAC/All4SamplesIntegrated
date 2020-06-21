# Read the Seurat RDS files from DIEM Filtering
C10 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C10/C10')
C19 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C19/C19')
M17 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M17/M17')
M20 <- readRDS('C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M20/M20')


# Filter smaller libraries and low expressed genes
# For C10
genes <- apply(C10@assays$RNA@counts, 1,
               function(x){length(x[x>=1]) >= 100})
nuclei <- apply(C10@assays$RNA@counts, 2,
                function(x){length(x[x>=1]) >= 500})
C10seu <- C10[genes, nuclei]
dim(C10) - dim(C10seu)
# For C19
genes <- apply(C19@assays$RNA@counts, 1,
               function(x){length(x[x>=1]) >= 100})
nuclei <- apply(C19@assays$RNA@counts, 2,
                function(x){length(x[x>=1]) >= 500})
C19seu <- C19[genes, nuclei]
dim(C19) - dim(C19seu)
# For M17
genes <- apply(M17@assays$RNA@counts, 1,
               function(x){length(x[x>=1]) >= 100})
nuclei <- apply(M17@assays$RNA@counts, 2,
                function(x){length(x[x>=1]) >= 500})
M17seu <- M17[genes, nuclei]
dim(M17) - dim(M17seu)
# For M20
genes <- apply(M20@assays$RNA@counts, 1,
               function(x){length(x[x>=1]) >= 100})
nuclei <- apply(M20@assays$RNA@counts, 2,
                function(x){length(x[x>=1]) >= 500})
M20seu <- M20[genes, nuclei]
dim(M20) - dim(M20seu)


# Naming Batches
patient.names <- list.files("C:/Users/s4526789/Documents/snRNAseqData/Samples")
C10seu@meta.data$batch <- patient.names[1]
C19seu@meta.data$batch <- patient.names[2]
M17seu@meta.data$batch <- patient.names[3]
M20seu@meta.data$batch <- patient.names[4]


# Normalisation and Dim Reduction
# For C10
C10seu <- NormalizeData(C10seu, scale.factor = nrow(C10seu))
C10seu <- FindVariableFeatures(C10seu, selection.method = 'vst',
                               nfeatures = nrow(C10seu))
C10seu <- ScaleData(C10seu, features = VariableFeatures(C10seu))
C10seu <- RunPCA(C10seu, features = VariableFeatures(C10seu))
ElbowPlot(C10seu, ndims = 50, reduction = "pca")
C10seu <- FindNeighbors(C10seu, dims = 1:30)
C10seu <- FindClusters(C10seu, resolution = 0.25) # 6 clusters
DimPlot(C10seu, reduction = 'pca', label = TRUE) + NoLegend()
C10seu <- RunUMAP(C10seu, dims = 1:30, min.dist = 0.1,
                    n.epochs = 500, reduction = 'pca')
DimPlot(C10seu, reduction = 'umap', group.by = 'seurat_clusters')
# For C19
C19seu <- NormalizeData(C19seu, scale.factor = nrow(C19seu))
C19seu <- FindVariableFeatures(C19seu, selection.method = 'vst',
                               nfeatures = nrow(C19seu))
C19seu <- ScaleData(C19seu, features = VariableFeatures(C19seu))
C19seu <- RunPCA(C19seu, features = VariableFeatures(C19seu))
ElbowPlot(C19seu, ndims = 50, reduction = "pca")
C19seu <- FindNeighbors(C19seu, dims = 1:30)
C19seu <- FindClusters(C19seu, resolution = 0.25) # 6 clusters
DimPlot(C19seu, reduction = 'pca', label = TRUE) + NoLegend()
C19seu <- RunUMAP(C19seu, dims = 1:30, min.dist = 0.1,
                  n.epochs = 500, reduction = 'pca')
DimPlot(C19seu, reduction = 'umap', group.by = 'seurat_clusters')
# For M17
M17seu <- NormalizeData(M17seu, scale.factor = nrow(M17seu))
M17seu <- FindVariableFeatures(M17seu, selection.method = 'vst',
                               nfeatures = nrow(M17seu))
M17seu <- ScaleData(M17seu, features = VariableFeatures(M17seu))
M17seu <- RunPCA(M17seu, features = VariableFeatures(M17seu))
ElbowPlot(M17seu, ndims = 50, reduction = "pca")
M17seu <- FindNeighbors(M17seu, dims = 1:30)
M17seu <- FindClusters(M17seu, resolution = 0.25) # 4 clusters
DimPlot(M17seu, reduction = 'pca', label = TRUE) + NoLegend()
M17seu <- RunUMAP(M17seu, dims = 1:30, min.dist = 0.1,
                  n.epochs = 500, reduction = 'pca')
DimPlot(M17seu, reduction = 'umap', group.by = 'seurat_clusters')
# For M20
M20seu <- NormalizeData(M20seu, scale.factor = nrow(M20seu))
M20seu <- FindVariableFeatures(M20seu, selection.method = 'vst',
                               nfeatures = nrow(M20seu))
M20seu <- ScaleData(M20seu, features = VariableFeatures(M20seu))
M20seu <- RunPCA(M20seu, features = VariableFeatures(M20seu))
ElbowPlot(M20seu, ndims = 50, reduction = "pca")
M20seu <- FindNeighbors(M20seu, dims = 1:30)
M20seu <- FindClusters(M20seu, resolution = 0.25) # 7 clusters
DimPlot(M20seu, reduction = 'pca', label = TRUE) + NoLegend()
M20seu <- RunUMAP(M20seu, dims = 1:30, min.dist = 0.1,
                  n.epochs = 500, reduction = 'pca')
DimPlot(M20seu, reduction = 'umap', group.by = 'seurat_clusters')


#############################################################################
########################## HM_C10 ###########################################
#############################################################################
# Finding Cluster Markers using MAST
p.markers <- FindAllMarkers(C10seu, min.pct = 0.25, 
                            logfc.threshold = 0.25, test.use = "MAST")
head(p.markers)
p.markers$p_val_adj <- p.adjust(p.markers$p_val,
                                method = "bonferroni", n = nrow(p.markers))
p.topmarkers <- p.markers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)
write.csv(p.topmarkers, 'C:/Users/s4526789/Documents/snRNAseqData/Analysis/C10_CluMarkers.csv') 
DoHeatmap(C10seu, features = p.topmarkers$gene) + NoLegend()
# Cluster Annotation using CIPR and BluePrint+ENCODE reference database
library(CIPR)
CIPR(p.markers, reference = c("blueprint", "hsrnaseq"),
     plot_ind = TRUE, global_results_obj = TRUE, global_plot_obj = TRUE)
DT::datatable(CIPR_top_results)
p.CIPRlabel <- data.frame(CIPR_top_results)
write.csv(p.markers, file = "C:/Users/s4526789/Documents/snRNAseqData/Analysis/C10_CIPRlabels.csv")
p.CIPRlabel <- p.CIPRlabel %>% group_by(cluster) %>% 
  top_n(n = 1, wt = identity_score)
C10seu[["CIPR.labels"]] <- p.CIPRlabel$long_name[
  match(C10seu[[]][["seurat_clusters"]], p.CIPRlabel$cluster)]
umap3 <- DimPlot(C10seu, reduction = "umap", group.by = "CIPR.labels")
# Cluster Annotations using SingleR with BluePrint+ENCODE database
library(SingleR)
bp.se <- BlueprintEncodeData()
bp.se <- as(bp.se, 'SingleCellExperiment')
C10seu.sce <- as.SingleCellExperiment(C10seu)
colData(C10seu.sce)$seurat_clusters <- as.character(Idents(C10seu))  # go from factor to character
head(colData(C10seu.sce))
annotation.cluster <- SingleR(test = C10seu.sce, method = 'cluster',
                              clusters = C10seu.sce$ident,
                              ref = bp.se, labels = bp.se$label.main, 
                              de.method='wilcox', genes = 'de')
table(annotation.cluster$labels)
library(pheatmap)
plotScoreHeatmap(annotation.cluster)
C10seu[["SingleR.labels"]] <- annotation.cluster$labels[
  match(C10seu[[]][["seurat_clusters"]], rownames(annotation.cluster))]
umap4 <- DimPlot(C10seu, reduction = "umap", group.by = "SingleR.labels")
# Cluster Annotations using scCATCH with CellMatch database
library(scCATCH)
clu_markers <- findmarkergenes(C10seu, species = 'Human', 
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
C10seu[["scCATCH.labels"]] <- clu_ann$cell_type[
  match(C10seu[[]][["seurat_clusters"]], clu_ann$cluster)]
umap5 <- DimPlot(C10seu, reduction = "umap", group.by = "scCATCH.labels")
# Save the MultiPlot Cluster Annotations
multiplot(umap3, umap4, umap5, cols = 3)
# Cell Trajectory Analysis
library(slingshot)
all.sds <- slingshot(Embeddings(C10seu, "umap"), 
                     clusterLabels = C10seu$seurat_clusters) 
C10seu.sce$slingPseudotime_1 <- NULL
C10seu.sce <- slingshot(C10seu.sce, reducedDim = 'UMAP', 
                          clusterLabels = C10seu$seurat_clusters)
library(viridis) # colour palette
colors <- viridis(50, alpha = 1)
plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
lines(all.sds, lwd = 2, type = 'lineages', col = 'black')
plot(reducedDims(all.sds), 
     col = colors[cut(C10seu.sce$slingPseudotime_2, breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(C10seu.sce), lwd=2)
# Plotting Pseudotime Values for Each Lineage
pt <- slingPseudotime(all.sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(1, 2))
# Save the Slingshot Lineage Curves
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
  lines(all.sds, lwd = 2, col = 'black', type = 'lineages')
}
# GO Enrichment Analysis
# All Clusters analysed together to save the allRes file and pathway plot
p.GOannotate <- list()
for(cluster in 0:(max(as.integer(C10seu$seurat_clusters))-1)){
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
p.GOannotate # Save the Cluster-wise GO Annotation
# Up and Down Regulated Gene visualisation 
gene.up <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC > 0) 
            %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
gene.down <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC <= 0) 
              %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
up.features <- gene.up %>% filter(avg_logFC > 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
down.features <- gene.down %>% filter(avg_logFC <= 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
d1 <- DotPlot(C10seu, features = up.features$gene) + RotatedAxis() + 
  labs(x = "Up-Regulated Features", y = "Clusters")
d2 <- DotPlot(C10seu, features = down.features$gene) + RotatedAxis() + 
  labs(x = "Down-Regulated Features", y = "Clusters")
# Save the MultiPlot DotPlot for Up- & Down-regulated Gene Markers
multiplot(d1, d2, cols = 2)


#############################################################################
########################## HM_C19 ###########################################
#############################################################################
# Finding Cluster Markers using MAST
p.markers <- FindAllMarkers(C19seu, min.pct = 0.25, 
                            logfc.threshold = 0.25, test.use = "MAST")
head(p.markers)
p.markers$p_val_adj <- p.adjust(p.markers$p_val,
                                method = "bonferroni", n = nrow(p.markers))
p.topmarkers <- p.markers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)
write.csv(p.topmarkers, 'C:/Users/s4526789/Documents/snRNAseqData/Analysis/C19_CluMarkers.csv') 
DoHeatmap(C19seu, features = p.topmarkers$gene) + NoLegend()
# Cluster Annotation using CIPR and BluePrint+ENCODE reference database
library(CIPR)
CIPR(p.markers, reference = c("blueprint", "hsrnaseq"),
     plot_ind = TRUE, global_results_obj = TRUE, global_plot_obj = TRUE)
DT::datatable(CIPR_top_results)
p.CIPRlabel <- data.frame(CIPR_top_results)
write.csv(p.markers, file = "C:/Users/s4526789/Documents/snRNAseqData/Analysis/C19_CIPRlabels.csv")
p.CIPRlabel <- p.CIPRlabel %>% group_by(cluster) %>% 
  top_n(n = 1, wt = identity_score)
C19seu[["CIPR.labels"]] <- p.CIPRlabel$long_name[
  match(C19seu[[]][["seurat_clusters"]], p.CIPRlabel$cluster)]
umap3 <- DimPlot(C19seu, reduction = "umap", group.by = "CIPR.labels")
# Cluster Annotations using SingleR with BluePrint+ENCODE database
library(SingleR)
bp.se <- BlueprintEncodeData()
bp.se <- as(bp.se, 'SingleCellExperiment')
C19seu.sce <- as.SingleCellExperiment(C19seu)
colData(C19seu.sce)$seurat_clusters <- as.character(Idents(C19seu))  # go from factor to character
head(colData(C19seu.sce))
annotation.cluster <- SingleR(test = C19seu.sce, method = 'cluster',
                              clusters = C19seu.sce$ident,
                              ref = bp.se, labels = bp.se$label.main, 
                              de.method='wilcox', genes = 'de')
table(annotation.cluster$labels)
library(pheatmap)
plotScoreHeatmap(annotation.cluster)
C19seu[["SingleR.labels"]] <- annotation.cluster$labels[
  match(C19seu[[]][["seurat_clusters"]], rownames(annotation.cluster))]
umap4 <- DimPlot(C19seu, reduction = "umap", group.by = "SingleR.labels")
# Cluster Annotations using scCATCH with CellMatch database
library(scCATCH)
clu_markers <- findmarkergenes(C19seu, species = 'Human', 
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
C19seu[["scCATCH.labels"]] <- clu_ann$cell_type[
  match(C19seu[[]][["seurat_clusters"]], clu_ann$cluster)]
umap5 <- DimPlot(C19seu, reduction = "umap", group.by = "scCATCH.labels")
# Save the MultiPlot Cluster Annotations
multiplot(umap3, umap4, umap5, cols = 3)
# Cell Trajectory Analysis
library(slingshot)
all.sds <- slingshot(Embeddings(C19seu, "umap"), 
                     clusterLabels = C19seu$seurat_clusters) 
C19seu.sce$slingPseudotime_1 <- NULL
C19seu.sce <- slingshot(C19seu.sce, reducedDim = 'UMAP', 
                        clusterLabels = C19seu$seurat_clusters)
library(viridis) # colour palette
colors <- viridis(50, alpha = 1)
plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
lines(all.sds, lwd = 2, type = 'lineages', col = 'black')
plot(reducedDims(all.sds), 
     col = colors[cut(C19seu.sce$slingPseudotime_2, breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(C19seu.sce), lwd=2)
# Plotting Pseudotime Values for Each Lineage
pt <- slingPseudotime(all.sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(1, 2))
# Save the Slingshot Lineage Curves
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
  lines(all.sds, lwd = 2, col = 'black', type = 'lineages')
}
# GO Enrichment Analysis
# All Clusters analysed together to save the allRes file and pathway plot
p.GOannotate <- list()
for(cluster in 0:(max(as.integer(C19seu$seurat_clusters))-1)){
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
p.GOannotate # Save the Cluster-wise GO Annotation
# Up and Down Regulated Gene visualisation 
gene.up <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC > 0) 
            %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
gene.down <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC <= 0) 
              %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
up.features <- gene.up %>% filter(avg_logFC > 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
down.features <- gene.down %>% filter(avg_logFC <= 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
d1 <- DotPlot(C19seu, features = up.features$gene) + RotatedAxis() + 
  labs(x = "Up-Regulated Features", y = "Clusters")
d2 <- DotPlot(C19seu, features = down.features$gene) + RotatedAxis() + 
  labs(x = "Down-Regulated Features", y = "Clusters")
# Save the MultiPlot DotPlot for Up- & Down-regulated Gene Markers
multiplot(d1, d2, cols = 2)


#############################################################################
########################## HM_M17 ###########################################
#############################################################################
# Finding Cluster Markers using MAST
p.markers <- FindAllMarkers(M17seu, min.pct = 0.25, 
                            logfc.threshold = 0.25, test.use = "MAST")
head(p.markers)
p.markers$p_val_adj <- p.adjust(p.markers$p_val,
                                method = "bonferroni", n = nrow(p.markers))
p.topmarkers <- p.markers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)
write.csv(p.topmarkers, 'C:/Users/s4526789/Documents/snRNAseqData/Analysis/M17_CluMarkers.csv') 
DoHeatmap(M17seu, features = p.topmarkers$gene) + NoLegend()
# Cluster Annotation using CIPR and BluePrint+ENCODE reference database
library(CIPR)
CIPR(p.markers, reference = c("blueprint", "hsrnaseq"),
     plot_ind = TRUE, global_results_obj = TRUE, global_plot_obj = TRUE)
DT::datatable(CIPR_top_results)
p.CIPRlabel <- data.frame(CIPR_top_results)
write.csv(p.markers, file = "C:/Users/s4526789/Documents/snRNAseqData/Analysis/M17_CIPRlabels.csv")
p.CIPRlabel <- p.CIPRlabel %>% group_by(cluster) %>% 
  top_n(n = 1, wt = identity_score)
M17seu[["CIPR.labels"]] <- p.CIPRlabel$long_name[
  match(M17seu[[]][["seurat_clusters"]], p.CIPRlabel$cluster)]
umap3 <- DimPlot(M17seu, reduction = "umap", group.by = "CIPR.labels")
# Cluster Annotations using SingleR with BluePrint+ENCODE database
library(SingleR)
bp.se <- BlueprintEncodeData()
bp.se <- as(bp.se, 'SingleCellExperiment')
M17seu.sce <- as.SingleCellExperiment(M17seu)
colData(M17seu.sce)$seurat_clusters <- as.character(Idents(M17seu))  # go from factor to character
head(colData(M17seu.sce))
annotation.cluster <- SingleR(test = M17seu.sce, method = 'cluster',
                              clusters = M17seu.sce$ident,
                              ref = bp.se, labels = bp.se$label.main, 
                              de.method='wilcox', genes = 'de')
table(annotation.cluster$labels)
library(pheatmap)
plotScoreHeatmap(annotation.cluster)
M17seu[["SingleR.labels"]] <- annotation.cluster$labels[
  match(M17seu[[]][["seurat_clusters"]], rownames(annotation.cluster))]
umap4 <- DimPlot(M17seu, reduction = "umap", group.by = "SingleR.labels")
# Cluster Annotations using scCATCH with CellMatch database
library(scCATCH)
clu_markers <- findmarkergenes(M17seu, species = 'Human', 
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
M17seu[["scCATCH.labels"]] <- clu_ann$cell_type[
  match(M17seu[[]][["seurat_clusters"]], clu_ann$cluster)]
umap5 <- DimPlot(M17seu, reduction = "umap", group.by = "scCATCH.labels")
# Save the MultiPlot Cluster Annotations
multiplot(umap3, umap4, umap5, cols = 3)
# Cell Trajectory Analysis
library(slingshot)
all.sds <- slingshot(Embeddings(M17seu, "umap"), 
                     clusterLabels = M17seu$seurat_clusters) 
M17seu.sce$slingPseudotime_1 <- NULL
M17seu.sce <- slingshot(M17seu.sce, reducedDim = 'UMAP', 
                        clusterLabels = M17seu$seurat_clusters)
library(viridis) # colour palette
colors <- viridis(50, alpha = 1)
plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
lines(all.sds, lwd = 2, type = 'lineages', col = 'black')
plot(reducedDims(all.sds), 
     col = colors[cut(M17seu.sce$slingPseudotime_2, breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(M17seu.sce), lwd=2)
# Plotting Pseudotime Values for Each Lineage
pt <- slingPseudotime(all.sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(1, 2))
# Save the Slingshot Lineage Curves
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
  lines(all.sds, lwd = 2, col = 'black', type = 'lineages')
}
# GO Enrichment Analysis
# All Clusters analysed together to save the allRes file and pathway plot
p.GOannotate <- list()
for(cluster in 0:(max(as.integer(M17seu$seurat_clusters))-1)){
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
p.GOannotate # Save the Cluster-wise GO Annotation
# Up and Down Regulated Gene visualisation 
gene.up <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC > 0) 
            %>% filter(p_val_adj < 0.05) %>% top_n(n = 50, wt = p_val_adj))
gene.down <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC <= 0) 
              %>% filter(p_val_adj < 0.05) %>% top_n(n = 50, wt = p_val_adj))
up.features <- gene.up %>% filter(avg_logFC > 0) %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC)
down.features <- gene.down %>% filter(avg_logFC <= 0) %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC)
d1 <- DotPlot(M17seu, features = up.features$gene) + RotatedAxis() + 
  labs(x = "Up-Regulated Features", y = "Clusters")
d2 <- DotPlot(M17seu, features = down.features$gene) + RotatedAxis() + 
  labs(x = "Down-Regulated Features", y = "Clusters")
# Save the MultiPlot DotPlot for Up- & Down-regulated Gene Markers
multiplot(d1, d2, cols = 2)


#############################################################################
########################## HM_M20 ###########################################
#############################################################################
# Finding Cluster Markers using MAST
p.markers <- FindAllMarkers(M20seu, min.pct = 0.25, 
                            logfc.threshold = 0.25, test.use = "MAST")
head(p.markers)
p.markers$p_val_adj <- p.adjust(p.markers$p_val,
                                method = "bonferroni", n = nrow(p.markers))
p.topmarkers <- p.markers %>% group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC)
write.csv(p.topmarkers, 'C:/Users/s4526789/Documents/snRNAseqData/Analysis/M20_CluMarkers.csv') 
DoHeatmap(M20seu, features = p.topmarkers$gene) + NoLegend()
# Cluster Annotation using CIPR and BluePrint+ENCODE reference database
library(CIPR)
CIPR(p.markers, reference = c("blueprint", "hsrnaseq"),
     plot_ind = TRUE, global_results_obj = TRUE, global_plot_obj = TRUE)
DT::datatable(CIPR_top_results)
p.CIPRlabel <- data.frame(CIPR_top_results)
write.csv(p.markers, file = "C:/Users/s4526789/Documents/snRNAseqData/Analysis/M20_CIPRlabels.csv")
p.CIPRlabel <- p.CIPRlabel %>% group_by(cluster) %>% 
  top_n(n = 1, wt = identity_score)
M20seu[["CIPR.labels"]] <- p.CIPRlabel$long_name[
  match(M20seu[[]][["seurat_clusters"]], p.CIPRlabel$cluster)]
umap3 <- DimPlot(M20seu, reduction = "umap", group.by = "CIPR.labels")
# Cluster Annotations using SingleR with BluePrint+ENCODE database
library(SingleR)
bp.se <- BlueprintEncodeData()
bp.se <- as(bp.se, 'SingleCellExperiment')
M20seu.sce <- as.SingleCellExperiment(M20seu)
colData(M20seu.sce)$seurat_clusters <- as.character(Idents(M20seu))  # go from factor to character
head(colData(M20seu.sce))
annotation.cluster <- SingleR(test = M20seu.sce, method = 'cluster',
                              clusters = M20seu.sce$ident,
                              ref = bp.se, labels = bp.se$label.main, 
                              de.method='wilcox', genes = 'de')
table(annotation.cluster$labels)
library(pheatmap)
plotScoreHeatmap(annotation.cluster)
M20seu[["SingleR.labels"]] <- annotation.cluster$labels[
  match(M20seu[[]][["seurat_clusters"]], rownames(annotation.cluster))]
umap4 <- DimPlot(M20seu, reduction = "umap", group.by = "SingleR.labels")
# Cluster Annotations using scCATCH with CellMatch database
library(scCATCH)
clu_markers <- findmarkergenes(M20seu, species = 'Human', 
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
M20seu[["scCATCH.labels"]] <- clu_ann$cell_type[
  match(M20seu[[]][["seurat_clusters"]], clu_ann$cluster)]
umap5 <- DimPlot(M20seu, reduction = "umap", group.by = "scCATCH.labels")
# Save the MultiPlot Cluster Annotations
multiplot(umap3, umap4, umap5, cols = 3)
# Cell Trajectory Analysis
library(slingshot)
all.sds <- slingshot(Embeddings(M20seu, "umap"), 
                     clusterLabels = M20seu$seurat_clusters) 
M20seu.sce$slingPseudotime_1 <- NULL
M20seu.sce <- slingshot(M20seu.sce, reducedDim = 'UMAP', 
                        clusterLabels = M20seu$seurat_clusters)
library(viridis) # colour palette
colors <- viridis(50, alpha = 1)
plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
lines(all.sds, lwd = 2, type = 'lineages', col = 'black')
plot(reducedDims(all.sds), 
     col = colors[cut(M20seu.sce$slingPseudotime_2, breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(M20seu.sce), lwd=2)
# Plotting Pseudotime Values for Each Lineage
pt <- slingPseudotime(all.sds)
nms <- colnames(pt)
pal <- viridis(100, end = 0.95)
par(mfrow = c(2, 2))
# Save the Slingshot Lineage Curves
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(all.sds), col = colors, pch = 16, cex = 0.5)
  lines(all.sds, lwd = 2, col = 'black', type = 'lineages')
}
# GO Enrichment Analysis
# All Clusters analysed together to save the allRes file and pathway plot
p.GOannotate <- list()
for(cluster in 0:(max(as.integer(M20seu$seurat_clusters))-1)){
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
p.GOannotate # Save the Cluster-wise GO Annotation
# Up and Down Regulated Gene visualisation 
gene.up <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC > 0) 
            %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
gene.down <- (p.markers %>% group_by(cluster) %>% filter(avg_logFC <= 0) 
              %>% filter(p_val_adj < 0.05) %>% top_n(n = 20, wt = p_val_adj))
up.features <- gene.up %>% filter(avg_logFC > 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
down.features <- gene.down %>% filter(avg_logFC <= 0) %>% group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
d1 <- DotPlot(M20seu, features = up.features$gene) + RotatedAxis() + 
  labs(x = "Up-Regulated Features", y = "Clusters")
d2 <- DotPlot(M20seu, features = down.features$gene) + RotatedAxis() + 
  labs(x = "Down-Regulated Features", y = "Clusters")
# Save the MultiPlot DotPlot for Up- & Down-regulated Gene Markers
multiplot(d1, d2, cols = 2)




             