#############################################################################
# COMPLETE SINGLE-CELL RNA-SEQ PIPELINE: RAW READS -> FINAL SEURAT ANALYSIS
#############################################################################
#
# This pipeline has two phases:
#   PHASE 0 (shell, not R): FASTQ -> count matrix, via Cell Ranger or STARsolo
#   PHASE 1 (R / Seurat):   count matrix -> QC -> clustering -> annotation -> DE
#
# Run PHASE 0 in a terminal (not R). Everything from "PHASE 1" onward is R.
#
#############################################################################


#############################################################################
# PHASE 0: RAW READS -> COUNT MATRIX (run in bash, NOT in R)
#############################################################################

# ---------------------------------------------------------------------------
# OPTION A: 10x Genomics data -> Cell Ranger (industry standard, easiest)
# ---------------------------------------------------------------------------
# Install: https://www.10xgenomics.com/support/software/cell-ranger
#
# cellranger count \
#     --id=sample1 \
#     --transcriptome=/path/to/refdata-gex-GRCh38-2024-A \
#     --fastqs=/path/to/fastqs/sample1 \
#     --sample=sample1 \
#     --create-bam=true \
#     --localcores=16 \
#     --localmem=64
#
# Output you need for R is here:
#   sample1/outs/filtered_feature_bc_matrix/
#       barcodes.tsv.gz
#       features.tsv.gz
#       matrix.mtx.gz
#
# For multiple samples, repeat cellranger count per sample, then aggregate
# (optional, Seurat can also just merge the individual matrices):
#
# cellranger aggr --id=all_samples --csv=aggregation.csv

# ---------------------------------------------------------------------------
# OPTION B: Non-10x data (any droplet/plate platform) -> STARsolo (free, fast)
# ---------------------------------------------------------------------------
# STAR --runMode genomeGenerate \
#      --genomeDir /path/to/star_index \
#      --genomeFastaFiles genome.fa \
#      --sjdbGTFfile annotation.gtf \
#      --runThreadN 16
#
# STAR --runMode alignReads \
#      --genomeDir /path/to/star_index \
#      --readFilesIn R2.fastq.gz R1.fastq.gz \
#      --readFilesCommand zcat \
#      --soloType CB_UMI_Simple \
#      --soloCBwhitelist 3M-february-2018.txt \
#      --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
#      --soloFeatures Gene \
#      --outSAMtype BAM SortedByCoordinate \
#      --runThreadN 16
#
# Output: Solo.out/Gene/filtered/{barcodes.tsv,features.tsv,matrix.mtx}

# ---------------------------------------------------------------------------
# OPTION C: Alignment-free -> salmon alevin / alevin-fry (fastest, low memory)
# ---------------------------------------------------------------------------
# salmon alevin -l ISR -1 R1.fastq.gz -2 R2.fastq.gz \
#     --chromiumV3 -i salmon_index -p 16 -o sample1_alevin \
#     --tgMap t2g.tsv --sketch
# alevin-fry generate-permit-list ... / collate ... / quant ...
# (see alevin-fry docs; output loads into R via tximport/fishpond)

#############################################################################
# PHASE 1: R / SEURAT ANALYSIS STARTS HERE
#############################################################################

## ---- 0. Setup ------------------------------------------------------------

# install.packages(c("Seurat","tidyverse","patchwork","Matrix"))
# install.packages("BiocManager")
# BiocManager::install(c("glmGamPoi","SingleR","celldex","clustree","DoubletFinder"))

library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB, raise if you parallelize

## ---- 1. Load data ---------------------------------------------------------

# --- Single sample from Cell Ranger output ---
counts <- Read10X(data.dir = "sample1/outs/filtered_feature_bc_matrix/")
seu <- CreateSeuratObject(
  counts = counts,
  project = "my_project",
  min.cells = 3,     # gene must appear in >=3 cells to be kept
  min.features = 200 # cell must express >=200 genes to be kept
)

# --- Multiple samples: load each, then merge, keeping sample IDs ---
# sample_dirs <- c(s1 = "cellranger/s1/outs/filtered_feature_bc_matrix",
#                   s2 = "cellranger/s2/outs/filtered_feature_bc_matrix",
#                   s3 = "cellranger/s3/outs/filtered_feature_bc_matrix")
# seu_list <- lapply(names(sample_dirs), function(nm) {
#   CreateSeuratObject(Read10X(sample_dirs[nm]), project = nm,
#                       min.cells = 3, min.features = 200)
# })
# seu <- merge(seu_list[[1]], y = seu_list[-1], add.cell.ids = names(sample_dirs))
# seu <- JoinLayers(seu)  # Seurat v5: merge layers before most operations

## ---- 2. QC metrics ---------------------------------------------------------

# Mitochondrial % (dying/stressed cells have high mito content)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")     # human
# seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")   # mouse

# Ribosomal % (optional, useful diagnostic)
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

# Visualize before filtering
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## ---- 3. Filtering ----------------------------------------------------------
# Thresholds are dataset-specific -- inspect the violin plots above and adjust.
# Typical starting point for 10x human tissue:

seu <- subset(
  seu,
  subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
           nCount_RNA > 500 &
           percent.mt < 15
)

## ---- 4. Doublet detection (recommended) ------------------------------------
# library(DoubletFinder)
# seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>%
#   RunPCA() %>% RunUMAP(dims = 1:20)
# sweep.res <- paramSweep(seu, PCs = 1:20, sct = FALSE)
# sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
# bcmvn <- find.pK(sweep.stats)
# optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
# nExp <- round(0.075 * ncol(seu))  # ~7.5% expected doublet rate for 10x
# seu <- doubletFinder(seu, PCs = 1:20, pN = 0.25, pK = optimal_pK, nExp = nExp)
# seu <- subset(seu, subset = DF.classifications_0.25_..._... == "Singlet")

## ---- 5. Normalization -------------------------------------------------------

# Option A (recommended): SCTransform -- regresses out technical noise,
# replaces NormalizeData + FindVariableFeatures + ScaleData in one step
seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)

# Option B (classic / lighter weight):
# seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
# seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
# seu <- ScaleData(seu, vars.to.regress = c("percent.mt", "nCount_RNA"))

## ---- 6. Dimensionality reduction --------------------------------------------

seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Pick number of PCs to carry forward
ElbowPlot(seu, ndims = 50)
DimHeatmap(seu, dims = 1:15, cells = 500, balanced = TRUE)
# Common choice: 20-30 PCs unless elbow plot clearly suggests otherwise

## ---- 7. Clustering -----------------------------------------------------------

seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.8, algorithm = 1)  # Louvain
# Try several resolutions and inspect stability with clustree if unsure:
# library(clustree)
# seu <- FindClusters(seu, resolution = seq(0.2, 1.6, 0.2))
# clustree(seu, prefix = "SCT_snn_res.")

## ---- 8. Non-linear embedding for visualization --------------------------------

seu <- RunUMAP(seu, dims = 1:30)
# seu <- RunTSNE(seu, dims = 1:30)  # alternative

DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()

## ---- 9. Integration across samples/batches (if multiple samples) --------------
# Use if UMAP shows samples separating by technical batch rather than biology.
#
# seu.list <- SplitObject(seu, split.by = "orig.ident")
# seu.list <- lapply(seu.list, SCTransform, vars.to.regress = "percent.mt")
# features <- SelectIntegrationFeatures(seu.list, nfeatures = 3000)
# seu.list <- PrepSCTIntegration(seu.list, anchor.features = features)
# anchors <- FindIntegrationAnchors(seu.list, normalization.method = "SCT",
#                                    anchor.features = features, dims = 1:30)
# seu <- IntegrateData(anchors, normalization.method = "SCT", dims = 1:30)
# DefaultAssay(seu) <- "integrated"
# # then repeat steps 6-8 (RunPCA, FindNeighbors, FindClusters, RunUMAP) on this object
#
# Or, faster alternative using Seurat v5's built-in method:
# seu <- IntegrateLayers(seu, method = HarmonyIntegration, orig.reduction = "pca",
#                         new.reduction = "harmony")
# seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30) %>% FindClusters()
# seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)

## ---- 10. Marker genes / cluster identity ---------------------------------------

DefaultAssay(seu) <- "SCT"  # or "RNA" if not using SCTransform
seu <- PrepSCTFindMarkers(seu)  # required once before FindAllMarkers with SCT

markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(markers, "all_cluster_markers.csv", row.names = FALSE)

# Visualize top markers
DoHeatmap(seu, features = top10$gene) + NoLegend()
DotPlot(seu, features = unique(top10$gene)) + RotatedAxis()

## ---- 11. Cell type annotation ----------------------------------------------

# Option A: Manual annotation from known canonical markers, e.g.:
# FeaturePlot(seu, features = c("CD3E","CD14","MS4A1","NKG7","PPBP","FCGR3A"))
# new.cluster.ids <- c("T cells","Monocytes","B cells","NK cells", ...)  # length = n clusters
# names(new.cluster.ids) <- levels(seu)
# seu <- RenameIdents(seu, new.cluster.ids)
# seu$cell_type <- Idents(seu)

# Option B: Automated reference-based annotation with SingleR
# library(SingleR); library(celldex)
# ref <- celldex::HumanPrimaryCellAtlasData()
# pred <- SingleR(test = GetAssayData(seu, assay = "SCT", slot = "data"),
#                  ref = ref, labels = ref$label.main)
# seu$SingleR_label <- pred$labels
# DimPlot(seu, group.by = "SingleR_label", label = TRUE, repel = TRUE)

## ---- 12. Differential expression between conditions -----------------------
# If you have multiple conditions (e.g. treated vs control) in seu$condition:

# Idents(seu) <- "cell_type"
# de_results <- FindMarkers(
#   seu,
#   ident.1 = "treated", ident.2 = "control",
#   group.by = "condition",
#   subset.ident = "T cells"   # run per cell type
# )
# write.csv(de_results, "DE_Tcells_treated_vs_control.csv")

## ---- 13. Additional standard visualizations --------------------------------

VlnPlot(seu, features = c("CD3E", "CD14", "MS4A1"), pt.size = 0)
FeaturePlot(seu, features = c("CD3E", "CD14", "MS4A1"), order = TRUE)
DimPlot(seu, group.by = "orig.ident")   # check batch effects
DimPlot(seu, split.by = "orig.ident")   # side-by-side samples

## ---- 14. Save results -------------------------------------------------------

saveRDS(seu, file = "seurat_final_object.rds")
# Reload any time with: seu <- readRDS("seurat_final_object.rds")

# Export a lightweight metadata + embeddings table for sharing/plotting elsewhere
export_df <- cbind(seu@meta.data, Embeddings(seu, "umap"))
write.csv(export_df, "cell_metadata_with_umap.csv")

#############################################################################
# END OF PIPELINE
#############################################################################
# Typical run order once thresholds are set: source() this file section by
# section rather than all at once -- QC filtering (step 3), PC count (step 6),
# and clustering resolution (step 7) all benefit from inspecting plots first.
#############################################################################
