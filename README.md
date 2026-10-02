#scseurat 

- my seurat analysis pipeline 
- it covers the entire analysis 

```
FASTQ → count matrix (shell commands, not R): Cell Ranger for 10x data, or STARsolo / salmon alevin as free alternatives.
Load & create Seurat object, single or multi-sample.
QC: mito %, feature/count distributions, filtering thresholds.
Doublet detection with DoubletFinder (optional but recommended).
Normalization via SCTransform (or classic LogNormalize/ScaleData).
PCA, elbow plot, PC selection.
Clustering (FindNeighbors/FindClusters) with resolution sweeping via clustree.
UMAP/t-SNE visualization.
Batch integration (Harmony or CCA-based) if you have multiple samples showing batch effects.
Marker genes per cluster, heatmaps/dot plots.
Cell type annotation, manual or automated (SingleR + celldex).
Differential expression between conditions within cell types.
Standard plots and saving the final object + exportable tables
```

Gaurav Sablok \
gsablok@proton.me
