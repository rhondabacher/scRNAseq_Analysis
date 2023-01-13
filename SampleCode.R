# BiocManager::install("scRNAseq")
# BiocManager::install("TENxPBMCData")

library(TENxPBMCData)
tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
sce <- tenx_pbmc4k

library(scater)
sce <- addPerCellQC(sce, subsets=list(Mito=grep("MT", rowData(sce)$Symbol_TENx)))
names(colData(sce))


library(gridExtra)
# Cell-level QC metrics
# We expect to see an increasing number of detected genes with increasing total count. Each point represents a cell.
p1 <- plotColData(sce, x = "sum", y = "subsets_Mito_percent") # Total Count vs Mito percent
p2 <- plotColData(sce, x = "sum", y = "detected") # Total Count vs Detected Features
grid.arrange(p1, p2, ncol=2)


p1 <- plotColData(sce, y="sum")  + ggtitle("Total count")
p2 <- plotColData(sce, y="detected")  + ggtitle("Detected features")
p3 <- plotColData(sce, y="subsets_Mito_percent") + ggtitle("Mito percent")
grid.arrange(p1, p2, p3, ncol = 3)


# Identifying low-quality cells
# # With fixed thresholds
qc.ntotal <- colData(sce)$sum > 2e5
qc.nexprs <- colData(sce)$detected > 4e3
qc.mito <- colData(sce)$subsets_Mito_percent > 10

outliers <- qc.ntotal | qc.nexprs | qc.mito

# discard <- qc.lib | qc.nexprs| qc.mito
# # Summarize the number of cells removed for each reason.
# DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
#           MitoProp=sum(qc.mito), Total=sum(discard))

# With adaptive thresholds
# Identifying outliers
qc.ntotal <- isOutlier(colData(sce)$sum, type="both")
qc.nexprs <- isOutlier(colData(sce)$detected, type="both")
qc.mito <- isOutlier(colData(sce)$subsets_Mito_percent, type="higher")

outliers <- qc.ntotal | qc.nexprs | qc.mito

DataFrame(LibSize=sum(qc.ntotal), NExprs=sum(qc.nexprs),
          MitoProp=sum(qc.mito), Total=sum(outliers))

sce.filt <- sce[, !outliers]


p1 <- plotColData(sce.filt, x = "sum", y = "subsets_Mito_percent") # Total Count vs Mito percent
p2 <- plotColData(sce.filt, x = "sum", y = "detected") # Total Count vs Detected Features
grid.arrange(p1, p2, ncol=2)

p1 <- plotColData(sce.filt, y="sum")  + ggtitle("Total count")
p2 <- plotColData(sce.filt, y="detected")  + ggtitle("Detected features")
p3 <- plotColData(sce.filt, y="subsets_Mito_percent") + ggtitle("Mito percent")
grid.arrange(p1, p2, p3, ncol = 3)

# setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
# save.image("QC.RDATA")



#######################################
############ Normalization ############
#######################################

setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
load("QC.RDATA")

# Normalization
# SCnorm
# https://www.biostat.wisc.edu/~kendzior/SCNORM/SCnorm_vignette.pdf
counts <- assay(sce.filt, "counts")
libsizes <- colSums(counts)


# scTransform
library(Seurat)
library(sctransform)

colnames(sce.filt) <- colData(sce.filt)$Barcode
rownames(sce.filt) <- rowData(sce.filt)$Symbol_TENx
counts_sparse <- as(counts(sce.filt), "dgCMatrix") # we need to convert the DelayedMatrix to a sparse format to apply sctransform.

pbmc <- CreateSeuratObject(counts_sparse, project = "pbmc4k") 
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
# DefaultAssay(object = pbmc)

# setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
# save(pbmc, file="Norm_seurat.RDATA")

# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage


# Dino
library(Dino)
colnames(sce.filt) <- colData(sce.filt)$Barcode
rownames(sce.filt) <- rowData(sce.filt)$Symbol_TENx
counts_sparse <- as(counts(sce.filt), "dgCMatrix") # we need to convert the DelayedMatrix to a sparse format to apply sctransform.

norm_dino <- Dino(counts_sparse)
vignette("Dino")


# scone
library(scone)
# https://bioconductor.org/packages/devel/bioc/vignettes/scone/inst/doc/sconeTutorial.html#running-and-scoring-normalization-workflows-with-scone

scaling=list(none=identity, 
             # sum = SUM_FN, tmm = TMM_FN, 
             # uq = UQ_FN, fq = FQT_FN, psi = PSINORM_FN,
             deseq = DESEQ_FN)

g.filt <- rownames(sce.filt)[Matrix::rowSums(counts(sce.filt)) > 3]

sce.filt = sce.filt[g.filt,]
counts(sce.filt) <- as.matrix(counts(sce.filt))

my_scone <- SconeExperiment(sce.filt)
res <- scone(my_scone, scaling=scaling, eval_kclust = 2,
             k_ruv=0, k_qc=0, return_norm = "in_memory")

res <- my_scone

get_scores(res)



# Scran
# vignette("scran")
# https://bioconductor.statistik.tu-dortmund.de/packages/3.7/bioc/vignettes/scran/inst/doc/scran.html
library(scran)
clusters <- quickCluster(sce.filt)
sce.filt <- computeSumFactors(sce.filt, clusters=clusters)
sce.filt <- logNormCounts(sce.filt)


# lib_size <- colSums(log1p(sce.filt@assays@data@listData[["counts"]]))
# norm_lib_size <- colSums(sce.filt@assays@data@listData[["logcounts"]])
# 
# plot(lib_size, norm_lib_size, pch=19, col=alpha("grey", 0.7), 
#      xlab="log sequencing depth", ylab="log-normalized sequencing depth")



# setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
# save(sce.filt, file="Norm_sce.RDATA")

#############################################
############ Dimension Deduction ############
#############################################

# Seurat:
library(Seurat)
setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
load("Norm_seurat.RDATA")

#Feature Selection
# SCTransform an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow

# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# genes.seurat <- VariableFeatures(pbmc)
# top10 <- head(genes.seurat, 10)

# Dimensionality reduction:
# PCA
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, verbose = FALSE)
ElbowPlot(pbmc, reduction = "pca", ndims = 50) # 1:30
DimPlot(pbmc, reduction = "pca")
# UMAP
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:30, verbose = FALSE) # Run UMAP map on first 30 PCs
DimPlot(pbmc, reduction = "umap")
# Tsne
pbmc <- RunTSNE(pbmc, reduction = "pca", dims = 1:30, verbose = FALSE)
DimPlot(pbmc, reduction = "tsne")

# Scran
setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
load("Norm_sce.RDATA")

library(scran)
gene.var <- modelGeneVar(sce.filt)
genes.scran <- getTopHVGs(gene.var, n=2000)


##########################################
############ Cluster Analysis ############
##########################################

# If we are able to annotate the clusters by markers, we can do:

library(Seurat)

pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.8)


library(gridExtra)
p1 <- DimPlot(pbmc, reduction = "umap")
p2 <- DimPlot(pbmc, reduction = "tsne")
grid.arrange(p1, p2, ncol=2)

library(dplyr)
Idents(object = pbmc) <- "seurat_clusters"
markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, features=pbmc@assays[["SCT"]]@var.features)
df_top <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC) %>%
  as.data.frame()

df_top$gene
VlnPlot(pbmc, features = df_top$gene)

# markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()
# 
# pbmc@meta.data$SCT_snn_new <- pbmc@meta.data$seurat_clusters
# pbmc@meta.data$SCT_snn_new[pbmc@meta.data$seurat_clusters == "7"] <- "6"
# pbmc@meta.data$SCT_snn_new[pbmc@meta.data$seurat_clusters == "11"] <- "6"
# pbmc@meta.data$SCT_snn_new[pbmc@meta.data$seurat_clusters == "10"] <- "9"
# pbmc@meta.data$SCT_snn_new <- droplevels(pbmc@meta.data$SCT_snn_new)
# 
# Idents(object = pbmc) <- "SCT_snn_new"
# markers_updates <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, features=pbmc@assays[["SCT"]]@var.features)
# 
# df_top <- markers_updates %>%
#   group_by(cluster) %>%
#   slice_max(n = 1, order_by = avg_log2FC) %>%
#   as.data.frame()
# 
# df_top$gene
# VlnPlot(pbmc, features = df_top$gene)
# 
# markers_updates %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()
# p1 <- DimPlot(pbmc, reduction = "umap")
# p2 <- DimPlot(pbmc, reduction = "tsne")
# grid.arrange(p1, p2, ncol=2)



# If we don't have reliable annotation information, we can do a supervised analysis guided by a reference data set:

# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
# remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
setwd("/Users/ruby/Desktop/Graduate/Research/BookChapter/")
# reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
# save(reference, file="pbmc_reference.RDATA")
load("pbmc_reference.RDATA")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

anchors <- FindTransferAnchors(reference = reference, query = pbmc,
  normalization.method = "SCT", reference.reduction = "spca", dims = 1:50)

pbmc <- MapQuery(anchorset = anchors, query = pbmc, reference = reference,
  refdata = list(celltype.l1 = "celltype.l1",celltype.l2 = "celltype.l2", celltype.l3 = "celltype.l3", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap")
p1 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p3 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l3", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2 + p3

Idents(object = pbmc) <- "predicted.celltype.l1"

# reference$id <- 'reference'
# pbmc$id <- 'query'
# refquery <- merge(reference, pbmc)
# refquery[["spca"]] <- merge(reference[["spca"]], pbmc[["ref.spca"]])
# refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
# DimPlot(refquery, group.by = 'id', shuffle = TRUE)



############################
############ DE ############
############################

# levels(pbmc)

de.markers <- FindMarkers(pbmc, ident.1 = "CD4 T", ident.2 = "CD8 T")
head(de.markers)

de <- rownames(de.markers)[de.markers$p_val_adj<0.05]
pbmc_sub <- subset(x = pbmc, subset = (predicted.celltype.l1 == "CD4 T" | predicted.celltype.l1 == "CD8 T"))
seurat_mat <- GetAssayData(pbmc_sub, slot = "scale.data")
heatmap_mat <- seurat_mat[rownames(seurat_mat)[rownames(seurat_mat) %in% de],]
dis <- dist(heatmap_mat)
clusters <- hclust(dis)
DoHeatmap(pbmc_sub, features = clusters$labels[clusters$order], size = 4, angle = 0)



##############################
############ GSEA ############
##############################

# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html

library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

original_gene_list <- de.markers$avg_log2FC

# name the vector
names(original_gene_list) <- rownames(de.markers)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "none")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


#############################################
############ Trajectory Analysis ############
#############################################
library(slingshot)
library(RColorBrewer)
pbmc_sub <- subset(x = pbmc, subset = (predicted.celltype.l1 == "CD4 T" ))
umap_coord <- Embeddings(pbmc_sub, reduction = "umap")

sds <- slingshot(umap_coord)
pse <- slingPseudotime(sds)
head(pse)


line.ord <- slingCurves(sds)[["Lineage1"]][["s"]]
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pse, breaks=100)]

plot(Embeddings(pbmc_sub, reduction = "umap"), pch = 16, cex = 0.5, col = alpha(plotcol,0.7))
lines(line.ord, lwd=2, col='black')
