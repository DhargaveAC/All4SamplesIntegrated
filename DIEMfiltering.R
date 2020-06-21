library(BiocManager)
library(Seurat)
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(scran)
library(ggplot2)
library(ggExtra)
library(scater)
library(diem)
library(magrittr)
library(topGO)
library(scMerge)
library(org.Hs.eg.db)
library(gridExtra)
library(MAST)





# ----------------------- Initialising Data ----------------------------

m17 <- Read10X_h5("C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M17/raw_feature_bc_matrix.h5")
m20 <- Read10X_h5("C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M20/raw_feature_bc_matrix.h5")
c10 <- Read10X_h5("C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C10/raw_feature_bc_matrix.h5")
c19 <- Read10X_h5("C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C19/raw_feature_bc_matrix.h5")

m17 <- create_SCE(m17)
m20 <- create_SCE(m20)
c10 <- create_SCE(c10)
c19 <- create_SCE(c19)




#########################################
#----------------HM_M17------------------
#########################################
mt_genes <- grep(pattern="^mt-", x=rownames(m17@gene_data), 
                 ignore.case=TRUE, value=TRUE)
m17 <- get_gene_pct(x = m17, genes=mt_genes, name="pct.mt")
malat <- grep(pattern="^malat1$", x=rownames(m17@gene_data), 
              ignore.case=TRUE, value=TRUE)
m17 <- get_gene_pct(x = m17, genes=malat, name="MALAT1")
dim(m17)
# Droplet Metadata
m17.drop_data <- droplet_data(m17)
head(m17.drop_data)
summary(m17.drop_data)
# Plotting QC Metrices of the HM_M17
datf <- m17.drop_data
par(mfrow=c(2,2))
plot(datf$total_counts, datf$n_genes, pch = 16,
     xlab="Total Counts", ylab="Number Genes")
plot(datf$n_genes, datf$pct.mt, pch = 16,
     xlab="Number Genes", ylab="MT%")
plot(datf$n_genes, datf$MALAT1, pch = 16,
     xlab="Number Genes", ylab="MALAT1%")
plot(datf$pct.mt, datf$MALAT1, pch = 16,
     xlab="MT%", ylab="MALAT1%")
# Specifying test set and debris set 
barcode_rank_plot(m17, title = "MND Patient HM_M17")
m17 <- set_debris_test_set(m17)
length(m17@test_set)
length(m17@bg_set)
m17 <- diem(m17)
m17 <- call_targets(m17)
clean_ids <- get_clean_ids(m17)
debris_ids <- get_removed_ids(m17)
m17.drop_data <- droplet_data(m17)
summary(m17.drop_data)
m17.dd_test <- droplet_data(m17, type = "test")
summary(m17.dd_test)
dd_clean <- droplet_data(m17, type = "clean")
dd_debris <- droplet_data(m17, type = "debris")
p3 <- ggplot(m17.dd_test, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  scale_x_log10() + 
  xlab("Total Number of Genes (HM_M17)") + 
  ylab("% of Mitochondrial Genes")
umi3 <- plot_umi_gene_call(m17, alpha=0.3, ret = T)
ggplot(m17.drop_data, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  xlab("Number of Genes (HM_M17)") + 
  ylab("% of Mitochondrial-Mapped UMIs")
counts <- raw_counts(m17)
dim(counts)
clean_ids <- get_clean_ids(m17)
counts.clean <- counts[, clean_ids]
dim(counts.clean)
#m17 <- convert_to_seurat(m17, project = "HM_M17")




#########################################
#----------------HM_M20------------------
#########################################
mt_genes <- grep(pattern="^mt-", x=rownames(m20@gene_data), 
                 ignore.case=TRUE, value=TRUE)
m20 <- get_gene_pct(x = m20, genes=mt_genes, name="pct.mt")
malat <- grep(pattern="^malat1$", x=rownames(m20@gene_data), 
              ignore.case=TRUE, value=TRUE)
m20 <- get_gene_pct(x = m20, genes=malat, name="MALAT1")
dim(m20)
# Droplet Metadata
m20.drop_data <- droplet_data(m20)
head(m20.drop_data)
summary(m20.drop_data)
# Plotting QC Metrices of the HM_M20
datf <- drop_data
par(mfrow=c(2,2))
plot(datf$total_counts, datf$n_genes, pch = 16,
     xlab="Total Counts", ylab="Number Genes")
plot(datf$n_genes, datf$pct.mt, pch = 16,
     xlab="Number Genes", ylab="MT%")
plot(datf$n_genes, datf$MALAT1, pch = 16,
     xlab="Number Genes", ylab="MALAT1%")
plot(datf$pct.mt, datf$MALAT1, pch = 16,
     xlab="MT%", ylab="MALAT1%")
# Specifying test set and debris set 
barcode_rank_plot(m20, title = "MND Patient HM_M20")
m20 <- set_debris_test_set(m20)
length(m20@test_set)
length(m20@bg_set)
m20 <- diem(m20)
m20 <- call_targets(m20)
clean_ids <- get_clean_ids(m20)
debris_ids <- get_removed_ids(m20)
m20.drop_data <- droplet_data(m20)
summary(m20.drop_data)
m20.dd_test <- droplet_data(m20, type = "test")
summary(m20.dd_test)
dd_clean <- droplet_data(m20, type = "clean")
dd_debris <- droplet_data(m20, type = "debris")
p4 <- ggplot(m20.dd_test, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  scale_x_log10() + 
  xlab("Total Number of Genes (HM_M20)") + 
  ylab("% of Mitochondrial Genes")
umi4 <- plot_umi_gene_call(m20, alpha=0.3, ret = T)
ggplot(m20.drop_data, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  xlab("Number of Genes (HM_M20)") + 
  ylab("% of Mitochondrial-Mapped UMIs")
counts <- raw_counts(m20)
dim(counts)
clean_ids <- get_clean_ids(m20)
counts.clean <- counts[, clean_ids]
dim(counts.clean)
#m20 <- convert_to_seurat(m20, project = "HM_M20")




#########################################
#----------------HM_C10------------------
#########################################
mt_genes <- grep(pattern="^mt-", x=rownames(c10@gene_data), 
                 ignore.case=TRUE, value=TRUE)
c10 <- get_gene_pct(x = c10, genes=mt_genes, name="pct.mt")
malat <- grep(pattern="^malat1$", x=rownames(c10@gene_data), 
              ignore.case=TRUE, value=TRUE)
c10 <- get_gene_pct(x = c10, genes=malat, name="MALAT1")
dim(c10)
# Droplet Metadata
c10.drop_data <- droplet_data(c10)
head(c10.drop_data)
summary(c10.drop_data)
# Plotting QC Metrices of the HM_C10
datf <- c10.drop_data
par(mfrow=c(2,2))
plot(datf$total_counts, datf$n_genes, pch = 16,
     xlab="Total Counts", ylab="Number Genes")
plot(datf$n_genes, datf$pct.mt, pch = 16,
     xlab="Number Genes", ylab="MT%")
plot(datf$n_genes, datf$MALAT1, pch = 16,
     xlab="Number Genes", ylab="MALAT1%")
plot(datf$pct.mt, datf$MALAT1, pch = 16,
     xlab="MT%", ylab="MALAT1%")
# Specifying test set and debris set 
barcode_rank_plot(c10, title = "MND Patient HM_C10")
c10 <- set_debris_test_set(c10)
length(c10@test_set)
length(c10@bg_set)
c10 <- diem(c10)
c10 <- call_targets(c10)
clean_ids <- get_clean_ids(c10)
debris_ids <- get_removed_ids(c10)
c10.drop_data <- droplet_data(c10)
summary(c10.drop_data)
c10.dd_test <- droplet_data(c10, type = "test")
summary(c10.dd_test)
dd_clean <- droplet_data(c10, type = "clean")
dd_debris <- droplet_data(c10, type = "debris")
p1 <- ggplot(c10.dd_test, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  scale_x_log10() + 
  xlab("Total Number of Genes (HM_C10)") + 
  ylab("% of Mitochondrial Genes")
umi1 <- plot_umi_gene_call(c10, alpha=0.3, ret = T)
ggplot(c10.drop_data, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  xlab("Number of Genes (HM_C10)") + 
  ylab("% of Mitochondrial-Mapped UMIs")
counts <- raw_counts(c10)
dim(counts)
clean_ids <- get_clean_ids(c10)
counts.clean <- counts[, clean_ids]
dim(counts.clean)
#c10 <- convert_to_seurat(c10, project = "HM_C10")





#########################################
#----------------HM_C19------------------
#########################################
mt_genes <- grep(pattern="^mt-", x=rownames(c19@gene_data), 
                 ignore.case=TRUE, value=TRUE)
c19 <- get_gene_pct(x = c19, genes=mt_genes, name="pct.mt")
malat <- grep(pattern="^malat1$", x=rownames(c19@gene_data), 
              ignore.case=TRUE, value=TRUE)
c19 <- get_gene_pct(x = c19, genes=malat, name="MALAT1")
dim(c19)
# Droplet Metadata
c19.drop_data <- droplet_data(c19)
head(c19.drop_data)
summary(c19.drop_data)
# Plotting QC Metrices of the HM_C19
datf <- drop_data
par(mfrow=c(2,2))
plot(datf$total_counts, datf$n_genes, pch = 16,
     xlab="Total Counts", ylab="Number Genes")
plot(datf$n_genes, datf$pct.mt, pch = 16,
     xlab="Number Genes", ylab="MT%")
plot(datf$n_genes, datf$MALAT1, pch = 16,
     xlab="Number Genes", ylab="MALAT1%")
plot(datf$pct.mt, datf$MALAT1, pch = 16,
     xlab="MT%", ylab="MALAT1%")
# Specifying test set and debris set 
barcode_rank_plot(c19, title = "MND Patient HM_C19")
c19 <- set_debris_test_set(c19)
length(c19@test_set)
length(c19@bg_set)
c19 <- diem(c19)
c19 <- call_targets(c19)
clean_ids <- get_clean_ids(c19)
debris_ids <- get_removed_ids(c19)
c19.drop_data <- droplet_data(c19)
summary(c19.drop_data)
c19.dd_test <- droplet_data(c19, type = "test")
summary(c19.dd_test)
dd_clean <- droplet_data(c19, type = "clean")
dd_debris <- droplet_data(c19, type = "debris")
p2 <- ggplot(c19.dd_test, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + 
  theme_minimal() + 
  scale_x_log10() + 
  xlab("Total Number of Genes (HM_C19)") + 
  ylab("% of Mitochondrial Genes")
umi2 <- plot_umi_gene_call(c19, alpha=0.3, ret = T)
ggplot(c19.drop_data, aes(x = n_genes, y = pct.mt, color = Call)) + 
  geom_point() + theme_minimal() + 
  xlab("Number of Genes (HM_C19)") + 
  ylab("% of Mitochondrial-Mapped UMIs")
counts <- raw_counts(c19)
dim(counts)
clean_ids <- get_clean_ids(c19)
counts.clean <- counts[,clean_ids]
dim(counts.clean)
#c19 <- convert_to_seurat(c19, project = "HM_C19")



multiplot(p1, p3, p2, p4, cols = 2)
multiplot(umi1, umi3, umi2, umi4, cols = 2)

#---------- Convert into Seurat Objects ---------- 
c10 <- convert_to_seurat(c10, project = "HM_C10")
c19 <- convert_to_seurat(c19, project = "HM_C19")
m17 <- convert_to_seurat(m17, project = "HM_M17")
m20 <- convert_to_seurat(m20, project = "HM_M20")

#---------------------- Creating Patient List ------------------------------
patient.names <- list.files("C:/Users/s4526789/Documents/snRNAseqData/Samples")

#---------------------- Save RDS Files ----------------------
saveRDS(c10, file = 'C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C10/C10')
saveRDS(c19, file = 'C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_C19/C19')
saveRDS(m17, file = 'C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M17/M17')
saveRDS(m20, file = 'C:/Users/s4526789/Documents/snRNAseqData/Samples/HM_M20/M20')
