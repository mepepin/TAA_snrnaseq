---
title: "snRNA-sequencing of Aortic Dissection"
author: "Mark E. Pepin, MD, PhD, MS"
date: "01/09/2023"
output:
  html_document:
    code_folding: hide
    keep_md: yes
    toc: yes
    toc_float: yes
header-includes:
- \usepackage{booktabs}
- \usepackage{longtable}
- \usepackage{array}
- \usepackage{multirow}
- \usepackage[table]{xcolor}
- \usepackage{wrapfig}
- \usepackage{float}
- \usepackage{colortbl}
- \usepackage{pdflscape}
- \usepackage{tabu}
- \usepackage{threeparttable}
mainfont: Times
fontsize: 10pt
always_allow_html: yes
editor_options: 
  markdown: 
    wrap: 72
---



**Code Authors**: Mark E. Pepin, MD, PhD, MS **Contact**:
[pepinme\@gmail.com](mailto:mepepin@bwh.harvard.edu){.email}\
**Institution**: Brigham and Women's Hospital | Broad Institute of Harvard and MIT\
**Location**: Boston, MA

# Data Pre-Processing


```r
library(Seurat)
# Upload data provided
hdat <- readRDS(file = "../1_Input/snRNA_Human_ControlDissection_integrated_seurat_v2.rds")
# Initialize the Seurat object with the raw (non-normalized data).
hdat[["percent.mt"]] <- PercentageFeatureSet(hdat, pattern = "^MT-")
hdat <- subset(hdat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
hdat <- NormalizeData(hdat, normalization.method = "LogNormalize", scale.factor = 10000)
hdat <- FindVariableFeatures(hdat, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(hdat)
hdat <- ScaleData(hdat, features = all.genes)
hdat[["Celltype_annotation_v2"]]$Celltype_annotation_v2 <- factor(hdat[["Celltype_annotation_v2"]]$Celltype_annotation_v2, c("VSMC1", "VSMC2", "Myofibroblast", "Fibroblast1", "Fibroblast2", "EC1", "EC2", "Macrophage", "NKT Cells", "Neuronal"))
##
TAA_list<-SplitObject(hdat, split.by = "Disease") # split dataset into a list of two seruat objects
# normalize and identify variable features for each dataset independently
TAA.list <- lapply(X = TAA_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = TAA.list)
# Preserve cell type markers in the anchoring dataset
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1")
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c( "NFASC",  "SAMD5","UACA","PRSS23") #  "TEX41",
Fibromyocyte_genes<-c("DGKG", "ADAMTS1", "RGS6", "TNC", "GRIP2", "ANGPT2")
EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") #, "NOTCH4"
EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
NKT_genes<-c("SKAP1", "RIPOR2", "RBPJ", "FYN", "ITGAL", "CD96")
Macrophage_genes<-c("MRC1", "LGMN", "F13A1", "RBM47")
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8")
features<-union(features, union(VSMC1_genes,VSMC2_genes))
#   VSMC1_genes, union(
#     VSMC2_genes, union(
#       Fibroblast1_genes, union(
#         Fibroblast2_genes, union(
#           Fibromyocyte_genes, union(
#             EC1_genes, union(
#               EC2_genes, union(
#                 NKT_genes, union(Macrophage_genes, Dendritic_genes)
#                 )
#               )
#             )
#           )
#         )
#       )
#   )
# )
# )
# Identify "anchors" (cells that are unchanged across datasets)
TAA.anchors <- FindIntegrationAnchors(object.list = TAA.list, anchor.features = features)
# this command creates an 'integrated' data assay
TAA.combined <- IntegrateData(anchorset = TAA.anchors)
# Specify the matrix used in the downstream analysis
DefaultAssay(TAA.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
TAA.combined <- ScaleData(TAA.combined, verbose = FALSE)
TAA.combined <- RunPCA(TAA.combined, npcs = 30, verbose = FALSE)
TAA.combined <- RunUMAP(TAA.combined, dims = 1:30)
TAA.combined <- RunTSNE(TAA.combined)
TAA.combined <- FindNeighbors(TAA.combined, dims = 1:30)
TAA.combined <- FindClusters(TAA.combined, resolution = 0.4)
```

```
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 9278
## Number of edges: 381125
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8951
## Number of communities: 11
## Elapsed time: 0 seconds
```

```r
saveRDS(TAA.combined, file = "TAA_FFPE_snRNA.rds")
#
plot1 <- VariableFeaturePlot(hdat)
top10 <- head(VariableFeatures(hdat), 10)
plot2 <- LabelPoints(plot = plot1, points = c(top10), repel = TRUE)
pdf(file = "../2_Output/Figure_1/FFPE_Variable_genes.pdf")
plot2
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Figure 1
## Clustering


```r
# Clustering
umap_disease <- DimPlot(hdat, reduction = "umap", group.by = "Disease")
umap_clusters <- DimPlot(TAA.combined, reduction = "umap", label = TRUE, repel = TRUE)
pca_disease <- DimPlot(TAA.combined, reduction = "pca", group.by = "Disease")
pca_clusters <- DimPlot(TAA.combined, reduction = "pca", label = TRUE, repel = TRUE)
tsne_disease <- DimPlot(TAA.combined, reduction = "tsne", group.by = "Disease")
tsne_clusters <- DimPlot(TAA.combined, reduction = "tsne", label = TRUE, repel = TRUE)
pdf("../2_Output/Supplemental_Figures/Clustering_Comparison.pdf", width = 11, height = 8)
umap_disease+pca_disease+tsne_disease+umap_clusters+pca_clusters+tsne_clusters
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Add tracings around the points
library(ggplot2)
library(ggtrace)
library(ggthemes)
library(ggrepel)
library(dplyr)
# UMAP
umap_data <- as.data.frame(TAA.combined@reductions$umap@cell.embeddings) # Create a data.frame from the UMAP (to be used in ggplot)
# Disease Clustering
umap_data$Disease <- TAA.combined$Disease
umap_data$Disease <- factor(umap_data$Disease, levels = c("Control", "Dissection"))
umap_disease <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = Disease)) +
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() +
  theme(text = element_text(size = 14))
umap_disease
```

![](snRNA_Seq_files/figure-html/cluster_combined-1.png)<!-- -->

```r
### UMAP Clustering
umap_data$Cluster <- TAA.combined$seurat_clusters # Create a variable based on UMAP clusters (Seurat)
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ Cluster, umap_data, mean) # Label Clusters
umap_clusters <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = Cluster)) +
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() + 
  geom_text_repel(data = cluster_centers, aes(label = Cluster), size = 6) +
  theme(legend.position = "none", text = element_text(size = 14))
umap_clusters
```

![](snRNA_Seq_files/figure-html/cluster_combined-2.png)<!-- -->

```r
pdf(file = "../2_Output/Supplemental_Figures/UMAP_overlay.pdf", width = 11.5, height = 5, bg = "transparent")
umap_disease + umap_clusters
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
#Figure 1C -  Proportional Graph
library(dittoSeq)
dittoBarPlot(
    object = hdat,
    var = "Celltype_annotation_v2",
    group.by = "Disease")
```

![](snRNA_Seq_files/figure-html/cluster_combined-3.png)<!-- -->

```r
dittoBarPlot(
    object = hdat,
    var = "Celltype_annotation_v2",
    group.by = "Sample_name",
    split.by = "Disease")
```

![](snRNA_Seq_files/figure-html/cluster_combined-4.png)<!-- -->

```r
######################################################
############## Dimensionality
#####################################################
# to control for any confounders, we can add the 'vars.to.regress' parameter, including any variable in the metadata <- ScaleData(pbmc, vars.to.regress = "percent.mt")
TAA.combined <- RunPCA(TAA.combined, features = VariableFeatures(object = TAA.combined))
# Visualize the features/genes
VizDimLoadings(TAA.combined, dims = 1:5, reduction = "pca")
```

![](snRNA_Seq_files/figure-html/cluster_combined-5.png)<!-- -->

```r
DimPlot(TAA.combined, reduction = "umap")
```

![](snRNA_Seq_files/figure-html/cluster_combined-6.png)<!-- -->

```r
pdf(file = "../2_Output/Supplemental_Figures/PC_Heatmaps.pdf")
DimHeatmap(TAA.combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Determine the number of dimensions represented by the data
TAA.combined <- JackStraw(TAA.combined, num.replicate = 100)
TAA.combined <- ScoreJackStraw(TAA.combined, dims = 1:20)
pdf(file = "../2_Output/Supplemental_Figures/JackStraw.pdf")
JackStrawPlot(TAA.combined, dims = 1:20)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
pdf(file = "../2_Output/Supplemental_Figures/ElbowPlot.pdf")
ElbowPlot(TAA.combined)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
######################################################
############## Cluster Annotation
#####################################################
library(ggplot2)
library(dplyr)
library(Seurat)
# TAA.combined<-readRDS(file = "TAA.ALL_snRNA.rds") # contains ACTA2, normal, and aortic aneurism (Cheu et. al)
# Find all gene markers
DEGs_Clusters<-FindAllMarkers(TAA.combined, assay = "RNA")
write.csv(DEGs_Clusters, "../2_Output/Supplemental_Figures/DEGs_Clusters.csv")
#Identify Clusters corresponding with known gene markers:
pdf(file = "../2_Output/Supplemental_Figures/UMAP_split.pdf", height = 4, width = 7)
DimPlot(TAA.combined,  label = T, split.by = "Disease") + NoLegend()
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
pdf(file = "../2_Output/Supplemental_Figures/UMAP_samples.pdf", height = 4, width = 11)
DimPlot(TAA.combined,  label = T, split.by = "Sample_name") + NoLegend()
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Create a dot-bplot of the top 5 markers for each 
library(scCustomize)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral2"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_Clusters, num_genes = 10, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/Supplemental_Figures/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = TAA.combined, features = top5_markers, k = 10, colors_use_exp = myColor)
```

```
## [[1]]
```

```
## 
## [[2]]
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
saveRDS(TAA.combined, file = "TAA_FFPE_snRNA.rds")
```

## Cell-type Annotation


```r
library(Seurat)
library(ggrepel)
TAA.combined<-readRDS(file = "TAA_FFPE_snRNA.rds")
#Identify Clusters corresponding with known gene markers:
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1")
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c("NFASC",  "SAMD5", "PRSS23") #"UACA","TEX41",
Fibromyocyte_genes<-c("ADAMTS1", "RGS6", "TNC") # , "ANGPT2", "DGKG", "GRIP2"
EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") # , "NOTCH4"
EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
NKT_genes<-c("SKAP1", "RIPOR2", "ITGAL", "CD96") #  "RBPJ", "FYN",
Macrophage_genes<-c("MRC1", "LGMN", "F13A1", "RBM47") #
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8") #
# Plot density function
library(Nebulosa)
VSMC1_density<-plot_density(TAA.combined, VSMC1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(TAA.combined, VSMC2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast1_density<-plot_density(TAA.combined, Fibroblast1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast2_density<-plot_density(TAA.combined, Fibroblast2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibromyocyte_density<-plot_density(TAA.combined, Fibromyocyte_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
EC1_density<-plot_density(TAA.combined, EC1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
EC2_density<-plot_density(TAA.combined, EC2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Macrophage_density<-plot_density(TAA.combined, Macrophage_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
NKT_density<-plot_density(TAA.combined, NKT_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
pdf(file = "../2_Output/Supplemental_Figures/Density_plots.pdf")
VSMC1_density
```

```
## $PKD1
```

```
## 
## $COL6A2
```

```
## 
## $PDLIM7
```

```
## 
## $FLNA
```

```
## 
## $SMTN
```

```
## 
## $`PKD1+ COL6A2+ PDLIM7+ FLNA+ SMTN+`
```

```r
VSMC2_density
```

```
## $MYH11
```

```
## 
## $ACTA2
```

```
## 
## $ITGA8
```

```
## 
## $PRKG1
```

```
## 
## $CRISPLD1
```

```
## 
## $`MYH11+ ACTA2+ ITGA8+ PRKG1+ CRISPLD1+`
```

```r
Fibroblast1_density
```

```
## $ABCA10
```

```
## 
## $C3
```

```
## 
## $ADGRD1
```

```
## 
## $FBLN1
```

```
## 
## $DCN
```

```
## 
## $`ABCA10+ C3+ ADGRD1+ FBLN1+ DCN+`
```

```r
Fibroblast2_density
```

```
## $NFASC
```

```
## 
## $SAMD5
```

```
## 
## $PRSS23
```

```
## 
## $`NFASC+ SAMD5+ PRSS23+`
```

```r
EC1_density
```

```
## $DIPK2B
```

```
## 
## $ARHGEF15
```

```
## 
## $STC1
```

```
## 
## $FLT1
```

```
## 
## $`DIPK2B+ ARHGEF15+ STC1+ FLT1+`
```

```r
EC2_density
```

```
## $VWF
```

```
## 
## $BMPER
```

```
## 
## $BMX
```

```
## 
## $NOS1
```

```
## 
## $`VWF+ BMPER+ BMX+ NOS1+`
```

```r
Macrophage_density
```

```
## $MRC1
```

```
## 
## $LGMN
```

```
## 
## $F13A1
```

```
## 
## $RBM47
```

```
## 
## $`MRC1+ LGMN+ F13A1+ RBM47+`
```

```r
NKT_density
```

```
## $SKAP1
```

```
## 
## $RIPOR2
```

```
## 
## $ITGAL
```

```
## 
## $CD96
```

```
## 
## $`SKAP1+ RIPOR2+ ITGAL+ CD96+`
```

```r
Fibromyocyte_density
```

```
## $ADAMTS1
```

```
## 
## $RGS6
```

```
## 
## $TNC
```

```
## 
## $`ADAMTS1+ RGS6+ TNC+`
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Overlay these gene markers onto the UMAP to identify clusters
pdf(file = "../2_Output/CellType_FeaturePlots.pdf")
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Fibroblast1_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Fibroblast2_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = EC1_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = EC2_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = NKT_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Macrophage_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = VSMC1_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = VSMC2_genes)
# FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Dendritic_genes)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
#Create a figure of cell-type specific markers overlying the UMAP
pdf(file = "../2_Output/Supplemental_Figures/CellType_Differentiation.pdf")
plot_density(TAA.combined, c("ACTA2", "FBLN1", "F13A1", "VWF", "STC1", "NFASC", "ITGAL", "ITGAX"), reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
```

```
## $ACTA2
```

```
## 
## $FBLN1
```

```
## 
## $F13A1
```

```
## 
## $VWF
```

```
## 
## $STC1
```

```
## 
## $NFASC
```

```
## 
## $ITGAL
```

```
## 
## $ITGAX
```

```
## 
## $`ACTA2+ FBLN1+ F13A1+ VWF+ STC1+ NFASC+ ITGAL+ ITGAX+`
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Change the identity of the clusters to cell types
Idents(TAA.combined) <- "Celltype_annotation_v2"
# Create UMAP with annotated cell-types
UMAP_CellTypes<-DimPlot(TAA.combined,  label = T) + NoLegend()
UMAP_CellTypes
```

![](snRNA_Seq_files/figure-html/merged.annotation-1.png)<!-- -->

```r
#######
umap_data <- as.data.frame(TAA.combined@reductions$umap@cell.embeddings) # Create a data.frame from the
umap_data$CellType<-TAA.combined@active.ident
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ CellType, umap_data, mean) # Label Clusters
umap_celltype <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = CellType)) +
  # geom_point() +
  # scale_fill_brewer(palette="Set2") + 
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() + 
  geom_text_repel(data = cluster_centers, aes(label = CellType), size = 5) +
  theme(legend.position = "none", text = element_text(size = 12))
## Fully Annotated Merged UMAP
pdf(file = "../2_Output/UMAP_Disease.Overlay.pdf", height = 5, width = 5)
umap_celltype
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
pdf(file = "../2_Output/Supplemental_Figures/UMAP_All.three_Labelled.pdf", height = 8, width = 16)
umap_disease+theme(legend.position = "topright", text = element_text(size = 14)) + 
umap_clusters+umap_celltype
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
## Use ridgeplots to identify bimodal gene marker distributions (enriched clusters)
pdf(file = "../2_Output/Celltype_RidgePlots.pdf", height = 10, width = 15)
RidgePlot(TAA.combined,features = Fibroblast1_genes, ncol = 2)
RidgePlot(TAA.combined, features = Fibroblast2_genes, ncol = 2)
RidgePlot(TAA.combined, features = EC1_genes, ncol = 2)
RidgePlot(TAA.combined, features = EC2_genes, ncol = 2)
RidgePlot(TAA.combined, features = NKT_genes, ncol = 2)
RidgePlot(TAA.combined, features = Macrophage_genes, ncol = 2)
RidgePlot(TAA.combined, features = VSMC1_genes, ncol = 2)
RidgePlot(TAA.combined, features = VSMC2_genes, ncol = 2)
RidgePlot(TAA.combined, features = Dendritic_genes, ncol = 2)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Differential Expression
# Export DEGs using cell-type clusters
DEGs_CellTypes<-FindAllMarkers(TAA.combined)
write.csv(DEGs_CellTypes, "../2_Output/DEGs_Clusters.csv")
# Create a dot-bplot of the top 5 markers for each 
library(scCustomize)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral2"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_CellTypes, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = TAA.combined, features = top5_markers, k = 10, colors_use_exp = myColor)
```

```
## [[1]]
```

```
## 
## [[2]]
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Heatmap of Clusters
DEGs_CellTypes %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
pdf(file = "../2_Output/Celltype_Heatmap.pdf", height = 10, width = 15)
DoHeatmap(TAA.combined, features = top10$gene, size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "coral2")) + labs(title = "Heatmap of Top-10 most variable genes within each cluster")
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Probability distribution of specific genes across all clusters/cell types
VlnPlot(TAA.combined, features = c("VWF", "HDAC9", "ACTA2"))
```

![](snRNA_Seq_files/figure-html/merged.annotation-2.png)<!-- -->

```r
#
DotPlot(TAA.combined, features = c(VSMC1_genes, VSMC2_genes, EC1_genes, EC2_genes, Fibroblast1_genes, Fibroblast2_genes, Fibromyocyte_genes, Macrophage_genes, NKT_genes)) + RotatedAxis()
```

![](snRNA_Seq_files/figure-html/merged.annotation-3.png)<!-- -->

```r
saveRDS(TAA.combined, file = "TAA_FFPE_snRNA.rds")
```


# Figure 2: VSMC-Specific Analysis

### VSMC Phenotypic Comparison


```r
library(Seurat)
library(dplyr)
library(ggtrace)
library(ggplot2)
library(ggrepel)
# Add tracings around the points
library(ggplot2)
library(ggtrace)
library(ggthemes)
library(ggrepel)
library(dplyr)
hdat_vsmc <- subset(TAA.combined, Celltype_annotation_v2==c("VSMC1", "VSMC2"))
UMAP_VSMC<-DimPlot(hdat_vsmc,  label = T) + NoLegend()
# Figure 1A - UMAP of VSMCs only
pdf("../2_Output/Figure_2/Fig2A_UMAP_VSMC.pdf", height = 3.5, width = 3.5)
UMAP_VSMC
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
library(Nebulosa)
VSMC1_density<-plot_density(hdat_vsmc, VSMC1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(hdat_vsmc, VSMC2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
pdf(file = "../2_Output/Figure_2/Fig2A_VSMC_Density_plots.pdf", height = 3.5, width = 3.5)
VSMC1_density
```

```
## $PKD1
```

```
## 
## $COL6A2
```

```
## 
## $PDLIM7
```

```
## 
## $FLNA
```

```
## 
## $SMTN
```

```
## 
## $`PKD1+ COL6A2+ PDLIM7+ FLNA+ SMTN+`
```

```r
VSMC2_density
```

```
## $MYH11
```

```
## 
## $ACTA2
```

```
## 
## $ITGA8
```

```
## 
## $PRKG1
```

```
## 
## $CRISPLD1
```

```
## 
## $`MYH11+ ACTA2+ ITGA8+ PRKG1+ CRISPLD1+`
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Visualization
p1 <- DimPlot(hdat_vsmc, reduction = "umap", group.by = "Disease")
p2 <- DimPlot(hdat_vsmc, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
```

![](snRNA_Seq_files/figure-html/subpopulation-1.png)<!-- -->

```r
umap_data <- as.data.frame(hdat_vsmc@reductions$umap@cell.embeddings) # Create a data.frame from the UMAP (to be used in ggplot)
# Disease Clustering
umap_data$Disease <- hdat_vsmc$Disease
umap_data$Disease <- factor(umap_data$Disease, levels = c("Control", "Dissection"))
p1 <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = Disease)) +
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() +
  theme(text = element_text(size = 14))
p1
```

![](snRNA_Seq_files/figure-html/subpopulation-2.png)<!-- -->

```r
umap_TAA<-umap_data %>% filter(Disease=="Control")
p_TAA <- ggplot(umap_TAA, aes(x = UMAP_1, y = UMAP_2, fill = Disease)) +
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() +
  theme(text = element_text(size = 14))
p_TAA
```

![](snRNA_Seq_files/figure-html/subpopulation-3.png)<!-- -->

```r
umap_ACTA<-umap_data %>% filter(Disease=="Dissection")
p_ACTA <- ggplot(umap_ACTA, aes(x = UMAP_1, y = UMAP_2, fill = Disease)) +
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() +
  theme(text = element_text(size = 14))
p1 + p_TAA + p_ACTA
```

![](snRNA_Seq_files/figure-html/subpopulation-4.png)<!-- -->

```r
# Create the UMAP of cell-type specific clusters
UMAP_CellTypes<-DimPlot(hdat_vsmc, label = T) + NoLegend()
UMAP_CellTypes
```

![](snRNA_Seq_files/figure-html/subpopulation-5.png)<!-- -->

```r
umap_data <- as.data.frame(hdat_vsmc@reductions$umap@cell.embeddings) # Create a data.frame from the
umap_data$Cluster<-hdat_vsmc@active.ident
cluster_centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ Cluster, umap_data, mean) # Label Clusters
p3 <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, fill = Cluster)) +
  # geom_point() +
  # scale_fill_brewer(palette="Set2") + 
  geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme_classic() + 
  geom_text_repel(data = cluster_centers, aes(label = Cluster), size = 6) +
  theme(legend.position = "none", text = element_text(size = 14))
p3
```

![](snRNA_Seq_files/figure-html/subpopulation-6.png)<!-- -->

```r
pdf(file = "../2_Output/Figure_2/VSMC_UMAP_CellTypes.pdf", width = 11.5, height = 5, bg = "transparent")
p1+theme(legend.position = "top", text = element_text(size = 14)) + p3
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Differential Expression between VSMCs
VSMC1.markers <- FindMarkers(hdat_vsmc, ident.1 = "VSMC1", min.pct = 0.25, logfc.threshold = 0.25)
VSMC2.markers <- FindMarkers(hdat_vsmc, ident.1 = "VSMC2", min.pct = 0.25, logfc.threshold = 0.25)

#Volcano Plot
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
# Read data from the web
options(ggrepel.max.overlaps = Inf)
results = mutate(VSMC2.markers, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.05 & log2FC>.5, 
                       "P < 0.05 and Fold-Change > 0.5", 
                       ifelse(p_val<0.05 & log2FC< 0-0.5,
                              "P < 0.05 and Fold-Change < -0.5", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.05 and Fold-Change < -0.5",
  "Not Sig",
  "P < 0.05 and Fold-Change > 0.5")
  )
max(results$minuslogpvalue, na.rm = TRUE)
```

```
## [1] 376.5639
```

```r
max(results$log2FC, na.rm = TRUE)
```

```
## [1] 2.104613
```

```r
min(results$log2FC, na.rm = TRUE)
```

```
## [1] -2.133374
```

```r
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank()
  ) +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=0.5, size=.5, linetype="dashed") +
  geom_vline(xintercept=-0.5, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.05), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  ylim(0, max(results$minuslogpvalue, na.rm = TRUE)) + geom_hline(yintercept = 0, size = 1) + 
  geom_vline(xintercept=0, size=1) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "coral2")) +
  scale_size_continuous(range = c(.1, 3))

  p+
  geom_text_repel(data=top_n(filter(results, log2FC< -0.5), 10, minuslogpvalue), aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>0.5), 10, minuslogpvalue), aes(label=gene_name)) +
  theme(text = element_text(size=20))
```

![](snRNA_Seq_files/figure-html/subpopulation-7.png)<!-- -->

```r
##
pdf(file = paste0("../2_Output/Figure_2/VSMC_1.vs.2_VolcanoPlot.pdf"), height = 5, width = 5)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< -0.5), 10, minuslogpvalue), aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>0.5), 10, minuslogpvalue), aes(label=gene_name)) +
    theme(text = element_text(size=14), legend.position="none")
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Create a file of differentially expressed genes between VSMC1 and VSMC2
openxlsx::write.xlsx(VSMC1.markers, file = "../2_Output/Figure_2/VSMC_1v2_DEGs.xlsx", rowNames=T)

VSMC_type_DEs <- results$gene_name[1:100]
##
plots_VSMC2 <- VlnPlot(hdat_vsmc, features = c("PALLD", "MAP2", "PRKG1"), group.by = "Celltype_annotation_v2",
    pt.size = 0) & theme(legend.position = 'none', axis.title.x = element_blank())
plots_VSMC1 <- VlnPlot(hdat_vsmc, features = c("B2M", "TNFAIP2", "STAB1"), group.by = "Celltype_annotation_v2",
    pt.size = 0) & theme(legend.position = 'none', axis.title.x = element_blank())

pdf(file = "../2_Output/Figure_2/VSMC2_Top.Violins.pdf", height = 3, width = 6)
plots_VSMC2
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
pdf(file = "../2_Output/Figure_2/VSMC1_Top.Violins.pdf", height = 3, width = 6)
plots_VSMC1
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
##################
DoHeatmap(hdat_vsmc, features = VSMC_type_DEs, group.by = "Celltype_annotation_v2",size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "coral2",midpoint = 0)
```

![](snRNA_Seq_files/figure-html/subpopulation-8.png)<!-- -->

```r
##################
pdf(file = "../2_Output/Figure_2/DotPlot_VSMC_1.v.2.pdf", height = 10, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, features = VSMC_type_DEs, group.by = "Celltype_annotation_v2", x_lab_rotate=F, k = 2)
```

```
## [[1]]
```

```
## 
## [[2]]
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

### VSMC Trajectory Analysis (Monocle 3)


```r
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
Merged_scaled<-hdat_vsmc
Merged_scaled@active.assay = "RNA"
cds<-SeuratWrappers::as.cell_data_set(Merged_scaled)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- estimate_size_factors(cds)
# Include gene names (not done by default by the seurat conversion)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
#  Run Monocle
cds <- cluster_cells(cds) # This step creates "partitions" that are used in the trajectory inference
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition") # this shows the partitions overlain on the UMAP
```

![](snRNA_Seq_files/figure-html/monocle3_Merged-1.png)<!-- -->

```r
cds <- learn_graph(cds, use_partition = TRUE) # creating trajectory inference within each partition
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |======================================================================| 100%
```

```r
# cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
root_cell <- "TACTCGCTCCCTGACT-1-1" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)

# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_branch_points = FALSE,
           label_leaves = FALSE)
```

![](snRNA_Seq_files/figure-html/monocle3_Merged-2.png)<!-- -->

```r
pdf(file = "../2_Output/Figure_2/VSMC_UMAP_Trajectory_Partition.pdf", height = 3, width = 3)
plot_cells(cds,
           color_cells_by = "partition",
           graph_label_size = 1,
           cell_size = .5,
           label_roots = F,
           label_branch_points = FALSE,
           label_leaves = FALSE)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
pdf(file = "../2_Output/Figure_2/VSMC_UMAP_Trajectory_Pseudotime.pdf", height = 3, width = 4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           graph_label_size = 1,
           cell_size = .7,
           label_branch_points = FALSE,
           label_leaves = F)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
```

```
## 
  |                                                              |   0%, ETA NA
  |                                                           |   0%, ETA 03:00
  |                                                           |   0%, ETA 00:33
  |                                                           |   0%, ETA 00:26
  |                                                           |   0%, ETA 00:24
  |                                                           |   0%, ETA 00:20
  |                                                           |   0%, ETA 00:17
  |                                                           |   0%, ETA 00:16
  |                                                           |   0%, ETA 00:14
  |                                                           |   0%, ETA 00:13
  |                                                           |   0%, ETA 00:10
  |                                                           |   0%, ETA 00:10
  |                                                           |   0%, ETA 00:10
  |                                                           |   1%, ETA 00:10
  |                                                           |   1%, ETA 00:10
  |                                                           |   1%, ETA 00:10
  |                                                           |   1%, ETA 00:09
  |                                                           |   1%, ETA 00:09
  |                                                           |   1%, ETA 00:09
  |                                                           |   1%, ETA 00:09
  |                                                           |   1%, ETA 00:09
  |                                                           |   1%, ETA 00:09
  |                                                           |   1%, ETA 00:09
  |=                                                          |   1%, ETA 00:09
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   1%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:07
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |=                                                          |   2%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:08
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   3%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |==                                                         |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   4%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:07
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   5%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |===                                                        |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   6%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   7%, ETA 00:06
  |====                                                       |   8%, ETA 00:06
  |====                                                       |   8%, ETA 00:06
  |====                                                       |   8%, ETA 00:06
  |====                                                       |   8%, ETA 00:06
  |====                                                       |   8%, ETA 00:06
  |====                                                       |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   8%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |=====                                                      |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |   9%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  10%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |======                                                     |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  11%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  12%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |=======                                                    |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  13%, ETA 00:06
  |========                                                   |  14%, ETA 00:06
  |========                                                   |  14%, ETA 00:06
  |========                                                   |  14%, ETA 00:06
  |========                                                   |  14%, ETA 00:06
  |========                                                   |  14%, ETA 00:06
  |========                                                   |  14%, ETA 00:06
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |========                                                   |  14%, ETA 00:05
  |=========                                                  |  14%, ETA 00:05
  |=========                                                  |  14%, ETA 00:05
  |=========                                                  |  14%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  15%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |=========                                                  |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  16%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  17%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |==========                                                 |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  18%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |===========                                                |  19%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  20%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |============                                               |  21%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:05
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  22%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |=============                                              |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  23%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  24%, ETA 00:04
  |==============                                             |  25%, ETA 00:04
  |==============                                             |  25%, ETA 00:04
  |==============                                             |  25%, ETA 00:04
  |==============                                             |  25%, ETA 00:04
  |==============                                             |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  25%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |===============                                            |  26%, ETA 00:04
  |================                                           |  26%, ETA 00:04
  |================                                           |  26%, ETA 00:04
  |================                                           |  26%, ETA 00:04
  |================                                           |  26%, ETA 00:04
  |================                                           |  26%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  27%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |================                                           |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  28%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  29%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |=================                                          |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  30%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |==================                                         |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  31%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  32%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |===================                                        |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  33%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  34%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |====================                                       |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:04
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  35%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |=====================                                      |  36%, ETA 00:03
  |======================                                     |  36%, ETA 00:03
  |======================                                     |  36%, ETA 00:03
  |======================                                     |  36%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  37%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |======================                                     |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  38%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  39%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |=======================                                    |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  40%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  41%, ETA 00:03
  |========================                                   |  42%, ETA 00:03
  |========================                                   |  42%, ETA 00:03
  |========================                                   |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  42%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |=========================                                  |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  43%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  44%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |==========================                                 |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  45%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  46%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |===========================                                |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  47%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |============================                               |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  48%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  49%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |=============================                              |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  50%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  51%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |==============================                             |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:03
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  52%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |===============================                            |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  53%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  54%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |================================                           |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  55%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  56%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |=================================                          |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  57%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |==================================                         |  58%, ETA 00:02
  |===================================                        |  58%, ETA 00:02
  |===================================                        |  58%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  59%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |===================================                        |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  60%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  61%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |====================================                       |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  62%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  63%, ETA 00:02
  |=====================================                      |  64%, ETA 00:02
  |=====================================                      |  64%, ETA 00:02
  |=====================================                      |  64%, ETA 00:02
  |=====================================                      |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  64%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |======================================                     |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  65%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  66%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |=======================================                    |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  67%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  68%, ETA 00:02
  |========================================                   |  69%, ETA 00:02
  |========================================                   |  69%, ETA 00:02
  |========================================                   |  69%, ETA 00:02
  |========================================                   |  69%, ETA 00:02
  |========================================                   |  69%, ETA 00:02
  |========================================                   |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  69%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |=========================================                  |  70%, ETA 00:02
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  70%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  71%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |==========================================                 |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  72%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  73%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |===========================================                |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  74%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |============================================               |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  75%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  76%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |=============================================              |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  77%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  78%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |==============================================             |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  79%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  80%, ETA 00:01
  |===============================================            |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  81%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |================================================           |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  82%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  83%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |=================================================          |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  84%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  85%, ETA 00:01
  |==================================================         |  86%, ETA 00:01
  |==================================================         |  86%, ETA 00:01
  |==================================================         |  86%, ETA 00:01
  |==================================================         |  86%, ETA 00:01
  |==================================================         |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  86%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |===================================================        |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  87%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  88%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |====================================================       |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  89%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:01
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  90%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |=====================================================      |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  91%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |======================================================     |  92%, ETA 00:00
  |=======================================================    |  92%, ETA 00:00
  |=======================================================    |  92%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  93%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |=======================================================    |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  94%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  95%, ETA 00:00
  |========================================================   |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  96%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |=========================================================  |  97%, ETA 00:00
  |========================================================== |  97%, ETA 00:00
  |========================================================== |  97%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  98%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |========================================================== |  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================|  99%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |===========================================================| 100%, ETA 00:00
  |=======================================================| 100%, Elapsed 00:05
```

```r
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes[modulated_genes$q_value < 0.05 & modulated_genes$status =="OK", ] # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05 & morans_I > 0.2))
openxlsx::write.xlsx(modulated_genes, "../2_Output/Figure_2/Pseudotime_DEGs.xlsx")

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y})) # Create a spline that smooths the pseudotime-based expression along the first 3 degrees of freedom.
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
pdf(file = "../2_Output/Figure_2/VSMC_Pseudotime.Heatmap.pdf", height = 4, width = 3)
print(hthc)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# Examine specific genes within the module(s) of interest (appears to emrich to metabolic genes)
library("viridis")
library(ggplot2)
GENES <- c("AEBP1", "LAMA2", "HDAC9", "ACTA2", "EXOC4", "PDE3A") #Genes selected based on visible trends
p1 <- plot_cells(cds, 
           genes = GENES[1:2],
           label_cell_groups = F,
           label_roots = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = F,
           min_expr = 1,
           alpha = 0.8,
           trajectory_graph_color = "grey28",
           ) +
  # geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme(
  axis.text = element_text(size = 6),    # Adjust the size as needed
  axis.title = element_text(size = 8),  # Adjust the size as needed
  legend.text = element_text(size = 4),  # Adjust the size as needed
  legend.title = element_text(size = 6),
  legend.key.size = unit(2, 'mm')) +
  scale_colour_viridis_c(option = "inferno")
p2 <- plot_cells(cds, 
           genes = GENES[3:4],
           label_cell_groups = F,
           label_roots = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = F,
           min_expr = 1,
           alpha = 0.8,
           trajectory_graph_color = "grey28",
           ) +
  # geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme(
  axis.text = element_text(size = 6),    # Adjust the size as needed
  axis.title = element_text(size = 8),  # Adjust the size as needed
  legend.text = element_text(size = 4),  # Adjust the size as needed
  legend.title = element_text(size = 6),
  legend.key.size = unit(2, 'mm')) +
  scale_colour_viridis_c(option = "inferno")
p3 <- plot_cells(cds, 
           genes = GENES[5:6],
           label_cell_groups = F,
           label_roots = F,
           label_branch_points = F,
           label_leaves = F,
           show_trajectory_graph = F,
           min_expr = 1,
           alpha = 0.8,
           trajectory_graph_color = "grey28",
           ) +
  # geom_point_trace(alpha = .25, stroke = .5, size = 1, color = "black")  +
  theme(
  axis.text = element_text(size = 6),    # Adjust the size as needed
  axis.title = element_text(size = 8),  # Adjust the size as needed
  legend.text = element_text(size = 4),  # Adjust the size as needed
  legend.title = element_text(size = 6),
  legend.key.size = unit(2, 'mm')) +
  scale_colour_viridis_c(option = "inferno")
pdf(file="../2_Output/Figure_2/UMAP_Pseudotime_Genes.pdf", height = 4, width = 3)
# ggarrange(p1, p2,p3, ncol=1, nrow=3, common.legend = TRUE, legend="right")
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
p1+p2+p3
```

![](snRNA_Seq_files/figure-html/monocle3_Merged-3.png)<!-- -->

```r
# Plot genes according to "pseudotime"
lineage_cds <- cds[rowData(cds)$gene_short_name %in% GENES, ] #colData(cds)$cell_type %in% c("VSMC")
lineage_cds<-order_cells(lineage_cds, root_cells = root_cell)
pdf(file="../2_Output/Figure_2/Pseudotime_curves.pdf", height = 7, width = 3)
monocle3::plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="ident",
                         min_expr=0)
dev.off()
```

```
## quartz_off_screen 
##                 2
```

######################################
### VSMC Differential Expression (Dissection vs. Control)
######################################


```r
library(ggpubr)
hdat_vsmc$VSMC_phenotype <- Idents(hdat_vsmc)
Idents(hdat_vsmc) <- "Disease"
de_VSMC <- FindMarkers(hdat_vsmc, ident.1 = "Dissection", ident.2 = "Control", verbose = FALSE, logfc.threshold = 0.25, min.pct = 0.25)
# Create a file of differentially expressed genes between VSMC1 and VSMC2
openxlsx::write.xlsx(de_VSMC, file = "../2_Output/Figure_2/VSMC_DEGs_Disease.xlsx", rowNames=T)
###
library(ggpubr)
plots_Disease_UP <- VlnPlot(hdat_vsmc, features = c("MACROD2", "CDK17", "FOS"), group.by = "Disease",
    pt.size = 0) & theme(axis.title.x = element_blank()) & scale_fill_manual(values = c("dodgerblue4", "goldenrod2"))
plots_Disease_DOWN <- VlnPlot(hdat_vsmc, features = c("PTPRQ", "HLA-B", "STAB1"), group.by = "Disease",
    pt.size = 0) & theme(axis.title.x = element_blank()) & scale_fill_manual(values = c("dodgerblue4", "goldenrod2"))
# ggarrange(plots_Disease_UP, plots_Disease_DOWN, ncol=1, nrow=2, common.legend = TRUE, legend="right")

pdf(file = "../2_Output/Figure_2/Dissection_up_Top.Violins.pdf", height = 3, width = 6)
plots_Disease_UP
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
pdf(file = "../2_Output/Figure_2/Dissection_down_Top.Violins.pdf", height = 3, width = 6)
plots_Disease_DOWN
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
##################
DoHeatmap(hdat_vsmc, features = rownames(de_VSMC), size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradient2(low = "dodgerblue4", mid = "white", high = "coral2",midpoint = 0)
```

![](snRNA_Seq_files/figure-html/DEGs_Disease-1.png)<!-- -->

```r
##################
pdf(file = "../2_Output/Figure_2/Dot.plot_top5markers.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = hdat_vsmc, features = rownames(de_VSMC), x_lab_rotate=F, k = 10)
```

```
## [[1]]
```

```
## 
## [[2]]
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
#Volcano Plot
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
# Read data from the web
options(ggrepel.max.overlaps = Inf)
results = mutate(de_VSMC, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.05 & log2FC>.5, 
                       "P < 0.05 and Fold-Change > 0.5", 
                       ifelse(p_val<0.05 & log2FC< 0-0.5,
                              "P < 0.05 and Fold-Change < -0.5", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.05 and Fold-Change < -0.5",
  "Not Sig",
  "P < 0.05 and Fold-Change > 0.5")
  )
max(results$minuslogpvalue, na.rm = TRUE)
```

```
## [1] 707.3396
```

```r
max(results$log2FC, na.rm = TRUE)
```

```
## [1] 1.164184
```

```r
min(results$log2FC, na.rm = TRUE)
```

```
## [1] -0.8923921
```

```r
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank()
  ) +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = 0.9) +
  geom_vline(xintercept=0.5, size=.5, linetype="dashed") +
  geom_vline(xintercept=-0.5, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.05), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  ylim(-0, max(results$minuslogpvalue, na.rm = TRUE)) + geom_hline(yintercept = 0, size = 1) + 
  geom_vline(xintercept=0, size=1) +
  scale_fill_manual(values=c("dodgerblue4", "darkgray", "goldenrod2")) +
  scale_size_continuous(range = c(.1, 3))

  p+
  geom_text_repel(data=top_n(filter(results, log2FC< -0.5), 20, minuslogpvalue), aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>0.5), 20, minuslogpvalue), aes(label=gene_name)) +
  theme(text = element_text(size=20))
```

![](snRNA_Seq_files/figure-html/DEGs_Disease-2.png)<!-- -->

```r
##
pdf(file = paste0("../2_Output/Figure_2/VSMC_Dissection.vs.Control_VolcanoPlot.pdf"), height = 5, width = 5)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< -0.5), 20, minuslogpvalue), aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>0.5), 20, minuslogpvalue), aes(label=gene_name)) +
    theme(text = element_text(size=14), legend.position="none")
dev.off()
```

```
## quartz_off_screen 
##                 2
```

```r
# 
# ##################
# DoHeatmap(hdat_vsmc, features = VSMC_type_DEs, group.by = "VSMC_phenotype",size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "coral2"))
##################
pdf(file = "../2_Output/Figure_2/DotPlot_VSMC_1.v.2.pdf", height = 10, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, features = VSMC_type_DEs, group.by = "VSMC_phenotype", x_lab_rotate=F, k = 2)
```

```
## [[1]]
```

```
## 
## [[2]]
```

```r
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Supplemental Table: R Session Information
All packages and setting are acquired using the following command:


```r
sinfo<-devtools::session_info()
sinfo$platform
```

```
##  setting  value
##  version  R version 4.2.2 (2022-10-31)
##  os       macOS Big Sur ... 10.16
##  system   x86_64, darwin17.0
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2023-01-14
##  pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
```

```r
sinfo$packages %>% kable( 
                         align="c", 
                         longtable=T, 
                         booktabs=T,
                         caption="Packages and Required Dependencies") %>% 
    kable_styling(latex_options=c("striped", "repeat_header", "condensed"))
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Packages and Required Dependencies</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:center;"> package </th>
   <th style="text-align:center;"> ondiskversion </th>
   <th style="text-align:center;"> loadedversion </th>
   <th style="text-align:center;"> path </th>
   <th style="text-align:center;"> loadedpath </th>
   <th style="text-align:center;"> attached </th>
   <th style="text-align:center;"> is_base </th>
   <th style="text-align:center;"> date </th>
   <th style="text-align:center;"> source </th>
   <th style="text-align:center;"> md5ok </th>
   <th style="text-align:center;"> library </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> abind </td>
   <td style="text-align:center;"> abind </td>
   <td style="text-align:center;"> 1.4.5 </td>
   <td style="text-align:center;"> 1.4-5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/abind </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/abind </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2016-07-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> assertthat </td>
   <td style="text-align:center;"> assertthat </td>
   <td style="text-align:center;"> 0.2.1 </td>
   <td style="text-align:center;"> 0.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/assertthat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/assertthat </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bslib </td>
   <td style="text-align:center;"> bslib </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/bslib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/bslib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cachem </td>
   <td style="text-align:center;"> cachem </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cachem </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cachem </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> callr </td>
   <td style="text-align:center;"> callr </td>
   <td style="text-align:center;"> 3.7.3 </td>
   <td style="text-align:center;"> 3.7.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/callr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/callr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-02 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cli </td>
   <td style="text-align:center;"> cli </td>
   <td style="text-align:center;"> 3.5.0 </td>
   <td style="text-align:center;"> 3.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cli </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cli </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-20 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cluster </td>
   <td style="text-align:center;"> cluster </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cluster </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cluster </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> codetools </td>
   <td style="text-align:center;"> codetools </td>
   <td style="text-align:center;"> 0.2.18 </td>
   <td style="text-align:center;"> 0.2-18 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/codetools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/codetools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-04 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> colorspace </td>
   <td style="text-align:center;"> colorspace </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/colorspace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/colorspace </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-02-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:center;"> cowplot </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cowplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/cowplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-12-30 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crayon </td>
   <td style="text-align:center;"> crayon </td>
   <td style="text-align:center;"> 1.5.2 </td>
   <td style="text-align:center;"> 1.5.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/crayon </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/crayon </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-29 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:center;"> data.table </td>
   <td style="text-align:center;"> 1.14.6 </td>
   <td style="text-align:center;"> 1.14.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/data.table </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/data.table </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DBI </td>
   <td style="text-align:center;"> DBI </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/DBI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/DBI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-06-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deldir </td>
   <td style="text-align:center;"> deldir </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/deldir </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/deldir </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-10-23 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:center;"> devtools </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/devtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/devtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> digest </td>
   <td style="text-align:center;"> digest </td>
   <td style="text-align:center;"> 0.6.31 </td>
   <td style="text-align:center;"> 0.6.31 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/digest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/digest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:center;"> dplyr </td>
   <td style="text-align:center;"> 1.0.10 </td>
   <td style="text-align:center;"> 1.0.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/dplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/dplyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-01 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ellipsis </td>
   <td style="text-align:center;"> ellipsis </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ellipsis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ellipsis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> evaluate </td>
   <td style="text-align:center;"> evaluate </td>
   <td style="text-align:center;"> 0.19 </td>
   <td style="text-align:center;"> 0.19 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/evaluate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/evaluate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-13 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fansi </td>
   <td style="text-align:center;"> fansi </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fansi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fansi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-24 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> farver </td>
   <td style="text-align:center;"> farver </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/farver </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/farver </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-06 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmap </td>
   <td style="text-align:center;"> fastmap </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fastmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fastmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-01-25 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fitdistrplus </td>
   <td style="text-align:center;"> fitdistrplus </td>
   <td style="text-align:center;"> 1.1.8 </td>
   <td style="text-align:center;"> 1.1-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-10 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fs </td>
   <td style="text-align:center;"> fs </td>
   <td style="text-align:center;"> 1.5.2 </td>
   <td style="text-align:center;"> 1.5.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/fs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future </td>
   <td style="text-align:center;"> future </td>
   <td style="text-align:center;"> 1.30.0 </td>
   <td style="text-align:center;"> 1.30.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/future </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/future </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future.apply </td>
   <td style="text-align:center;"> future.apply </td>
   <td style="text-align:center;"> 1.10.0 </td>
   <td style="text-align:center;"> 1.10.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/future.apply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/future.apply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-05 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> generics </td>
   <td style="text-align:center;"> generics </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/generics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/generics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-05 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:center;"> ggplot2 </td>
   <td style="text-align:center;"> 3.4.0 </td>
   <td style="text-align:center;"> 3.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-04 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrepel </td>
   <td style="text-align:center;"> ggrepel </td>
   <td style="text-align:center;"> 0.9.2 </td>
   <td style="text-align:center;"> 0.9.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ggrepel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ggrepel </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-06 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggridges </td>
   <td style="text-align:center;"> ggridges </td>
   <td style="text-align:center;"> 0.5.4 </td>
   <td style="text-align:center;"> 0.5.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ggridges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ggridges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-26 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> globals </td>
   <td style="text-align:center;"> globals </td>
   <td style="text-align:center;"> 0.16.2 </td>
   <td style="text-align:center;"> 0.16.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/globals </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/globals </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> glue </td>
   <td style="text-align:center;"> glue </td>
   <td style="text-align:center;"> 1.6.2 </td>
   <td style="text-align:center;"> 1.6.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/glue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/glue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-02-24 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> goftest </td>
   <td style="text-align:center;"> goftest </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/goftest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/goftest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-10-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridExtra </td>
   <td style="text-align:center;"> gridExtra </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/gridExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/gridExtra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtable </td>
   <td style="text-align:center;"> gtable </td>
   <td style="text-align:center;"> 0.3.1 </td>
   <td style="text-align:center;"> 0.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/gtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/gtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-01 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> highr </td>
   <td style="text-align:center;"> highr </td>
   <td style="text-align:center;"> 0.10 </td>
   <td style="text-align:center;"> 0.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/highr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/highr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmltools </td>
   <td style="text-align:center;"> htmltools </td>
   <td style="text-align:center;"> 0.5.4 </td>
   <td style="text-align:center;"> 0.5.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/htmltools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/htmltools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlwidgets </td>
   <td style="text-align:center;"> htmlwidgets </td>
   <td style="text-align:center;"> 1.6.1 </td>
   <td style="text-align:center;"> 1.6.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-01-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httpuv </td>
   <td style="text-align:center;"> httpuv </td>
   <td style="text-align:center;"> 1.6.7 </td>
   <td style="text-align:center;"> 1.6.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/httpuv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/httpuv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-14 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httr </td>
   <td style="text-align:center;"> httr </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/httr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/httr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ica </td>
   <td style="text-align:center;"> ica </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ica </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ica </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-08 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igraph </td>
   <td style="text-align:center;"> igraph </td>
   <td style="text-align:center;"> 1.3.5 </td>
   <td style="text-align:center;"> 1.3.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/igraph </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/igraph </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> irlba </td>
   <td style="text-align:center;"> irlba </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/irlba </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/irlba </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jquerylib </td>
   <td style="text-align:center;"> jquerylib </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/jquerylib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/jquerylib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jsonlite </td>
   <td style="text-align:center;"> jsonlite </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/jsonlite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/jsonlite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:center;"> kableExtra </td>
   <td style="text-align:center;"> 1.3.4 </td>
   <td style="text-align:center;"> 1.3.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/kableExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/kableExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-02-20 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KernSmooth </td>
   <td style="text-align:center;"> KernSmooth </td>
   <td style="text-align:center;"> 2.23.20 </td>
   <td style="text-align:center;"> 2.23-20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-05-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:center;"> knitr </td>
   <td style="text-align:center;"> 1.41 </td>
   <td style="text-align:center;"> 1.41 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/knitr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/knitr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> labeling </td>
   <td style="text-align:center;"> labeling </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/labeling </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/labeling </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-10-20 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> later </td>
   <td style="text-align:center;"> later </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/later </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/later </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-08-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lattice </td>
   <td style="text-align:center;"> lattice </td>
   <td style="text-align:center;"> 0.20.45 </td>
   <td style="text-align:center;"> 0.20-45 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lattice </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lattice </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lazyeval </td>
   <td style="text-align:center;"> lazyeval </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lazyeval </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lazyeval </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> leiden </td>
   <td style="text-align:center;"> leiden </td>
   <td style="text-align:center;"> 0.4.3 </td>
   <td style="text-align:center;"> 0.4.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/leiden </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/leiden </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-10 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lifecycle </td>
   <td style="text-align:center;"> lifecycle </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lifecycle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lifecycle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> limma </td>
   <td style="text-align:center;"> limma </td>
   <td style="text-align:center;"> 3.54.0 </td>
   <td style="text-align:center;"> 3.54.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/limma </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/limma </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-01 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> listenv </td>
   <td style="text-align:center;"> listenv </td>
   <td style="text-align:center;"> 0.9.0 </td>
   <td style="text-align:center;"> 0.9.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/listenv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/listenv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lmtest </td>
   <td style="text-align:center;"> lmtest </td>
   <td style="text-align:center;"> 0.9.40 </td>
   <td style="text-align:center;"> 0.9-40 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lmtest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/lmtest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> magrittr </td>
   <td style="text-align:center;"> magrittr </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/magrittr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/magrittr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-30 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS </td>
   <td style="text-align:center;"> MASS </td>
   <td style="text-align:center;"> 7.3.58.1 </td>
   <td style="text-align:center;"> 7.3-58.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/MASS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/MASS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Matrix </td>
   <td style="text-align:center;"> Matrix </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> 1.5-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Matrix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Matrix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> matrixStats </td>
   <td style="text-align:center;"> matrixStats </td>
   <td style="text-align:center;"> 0.63.0 </td>
   <td style="text-align:center;"> 0.63.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/matrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/matrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> memoise </td>
   <td style="text-align:center;"> memoise </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/memoise </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/memoise </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mime </td>
   <td style="text-align:center;"> mime </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/mime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/mime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-09-28 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miniUI </td>
   <td style="text-align:center;"> miniUI </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/miniUI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/miniUI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-05-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> munsell </td>
   <td style="text-align:center;"> munsell </td>
   <td style="text-align:center;"> 0.5.0 </td>
   <td style="text-align:center;"> 0.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/munsell </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/munsell </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-06-12 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nlme </td>
   <td style="text-align:center;"> nlme </td>
   <td style="text-align:center;"> 3.1.161 </td>
   <td style="text-align:center;"> 3.1-161 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/nlme </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/nlme </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-15 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> parallelly </td>
   <td style="text-align:center;"> parallelly </td>
   <td style="text-align:center;"> 1.33.0 </td>
   <td style="text-align:center;"> 1.33.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/parallelly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/parallelly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-14 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> patchwork </td>
   <td style="text-align:center;"> patchwork </td>
   <td style="text-align:center;"> 1.1.2 </td>
   <td style="text-align:center;"> 1.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/patchwork </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/patchwork </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pbapply </td>
   <td style="text-align:center;"> pbapply </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> 1.6-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pbapply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pbapply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pillar </td>
   <td style="text-align:center;"> pillar </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> 1.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pillar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pillar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgbuild </td>
   <td style="text-align:center;"> pkgbuild </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-27 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgconfig </td>
   <td style="text-align:center;"> pkgconfig </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgload </td>
   <td style="text-align:center;"> pkgload </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pkgload </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/pkgload </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plotly </td>
   <td style="text-align:center;"> plotly </td>
   <td style="text-align:center;"> 4.10.1 </td>
   <td style="text-align:center;"> 4.10.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/plotly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/plotly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plyr </td>
   <td style="text-align:center;"> plyr </td>
   <td style="text-align:center;"> 1.8.8 </td>
   <td style="text-align:center;"> 1.8.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/plyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/plyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> png </td>
   <td style="text-align:center;"> png </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/png </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/png </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-29 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polyclip </td>
   <td style="text-align:center;"> polyclip </td>
   <td style="text-align:center;"> 1.10.4 </td>
   <td style="text-align:center;"> 1.10-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/polyclip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/polyclip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-20 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> prettyunits </td>
   <td style="text-align:center;"> prettyunits </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/prettyunits </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/prettyunits </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> processx </td>
   <td style="text-align:center;"> processx </td>
   <td style="text-align:center;"> 3.8.0 </td>
   <td style="text-align:center;"> 3.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/processx </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/processx </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-26 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> profvis </td>
   <td style="text-align:center;"> profvis </td>
   <td style="text-align:center;"> 0.3.7 </td>
   <td style="text-align:center;"> 0.3.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/profvis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/profvis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-11-02 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> progressr </td>
   <td style="text-align:center;"> progressr </td>
   <td style="text-align:center;"> 0.12.0 </td>
   <td style="text-align:center;"> 0.12.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/progressr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/progressr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-13 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> promises </td>
   <td style="text-align:center;"> promises </td>
   <td style="text-align:center;"> 1.2.0.1 </td>
   <td style="text-align:center;"> 1.2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/promises </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/promises </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-02-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ps </td>
   <td style="text-align:center;"> ps </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ps </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ps </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-26 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:center;"> purrr </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/purrr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/purrr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-20 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R6 </td>
   <td style="text-align:center;"> R6 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/R6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/R6 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RANN </td>
   <td style="text-align:center;"> RANN </td>
   <td style="text-align:center;"> 2.6.1 </td>
   <td style="text-align:center;"> 2.6.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/RANN </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/RANN </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-01-08 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RColorBrewer </td>
   <td style="text-align:center;"> RColorBrewer </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-04-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rcpp </td>
   <td style="text-align:center;"> Rcpp </td>
   <td style="text-align:center;"> 1.0.9 </td>
   <td style="text-align:center;"> 1.0.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Rcpp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Rcpp </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-08 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppAnnoy </td>
   <td style="text-align:center;"> RcppAnnoy </td>
   <td style="text-align:center;"> 0.0.20 </td>
   <td style="text-align:center;"> 0.0.20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-27 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> remotes </td>
   <td style="text-align:center;"> remotes </td>
   <td style="text-align:center;"> 2.4.2 </td>
   <td style="text-align:center;"> 2.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/remotes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/remotes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape2 </td>
   <td style="text-align:center;"> reshape2 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/reshape2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/reshape2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reticulate </td>
   <td style="text-align:center;"> reticulate </td>
   <td style="text-align:center;"> 1.27 </td>
   <td style="text-align:center;"> 1.27 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/reticulate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/reticulate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-01-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rlang </td>
   <td style="text-align:center;"> rlang </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rlang </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rlang </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-24 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rmarkdown </td>
   <td style="text-align:center;"> rmarkdown </td>
   <td style="text-align:center;"> 2.19 </td>
   <td style="text-align:center;"> 2.19 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-15 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROCR </td>
   <td style="text-align:center;"> ROCR </td>
   <td style="text-align:center;"> 1.0.11 </td>
   <td style="text-align:center;"> 1.0-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ROCR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/ROCR </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstudioapi </td>
   <td style="text-align:center;"> rstudioapi </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> 0.14 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rtsne </td>
   <td style="text-align:center;"> Rtsne </td>
   <td style="text-align:center;"> 0.16 </td>
   <td style="text-align:center;"> 0.16 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Rtsne </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Rtsne </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-04-17 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rvest </td>
   <td style="text-align:center;"> rvest </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rvest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rvest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sass </td>
   <td style="text-align:center;"> sass </td>
   <td style="text-align:center;"> 0.4.4 </td>
   <td style="text-align:center;"> 0.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sass </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-24 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scales </td>
   <td style="text-align:center;"> scales </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/scales </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/scales </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-20 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scattermore </td>
   <td style="text-align:center;"> scattermore </td>
   <td style="text-align:center;"> 0.8 </td>
   <td style="text-align:center;"> 0.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/scattermore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/scattermore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-02-14 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sctransform </td>
   <td style="text-align:center;"> sctransform </td>
   <td style="text-align:center;"> 0.3.5 </td>
   <td style="text-align:center;"> 0.3.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sctransform </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sctransform </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sessioninfo </td>
   <td style="text-align:center;"> sessioninfo </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seurat </td>
   <td style="text-align:center;"> Seurat </td>
   <td style="text-align:center;"> 4.3.0 </td>
   <td style="text-align:center;"> 4.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Seurat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/Seurat </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeuratObject </td>
   <td style="text-align:center;"> SeuratObject </td>
   <td style="text-align:center;"> 4.1.3 </td>
   <td style="text-align:center;"> 4.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shiny </td>
   <td style="text-align:center;"> shiny </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/shiny </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/shiny </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-15 </td>
   <td style="text-align:center;"> CRAN (R 4.2.2) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sp </td>
   <td style="text-align:center;"> sp </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/sp </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.data </td>
   <td style="text-align:center;"> spatstat.data </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> 3.0-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.explore </td>
   <td style="text-align:center;"> spatstat.explore </td>
   <td style="text-align:center;"> 3.0.5 </td>
   <td style="text-align:center;"> 3.0-5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-10 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.geom </td>
   <td style="text-align:center;"> spatstat.geom </td>
   <td style="text-align:center;"> 3.0.3 </td>
   <td style="text-align:center;"> 3.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-25 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.random </td>
   <td style="text-align:center;"> spatstat.random </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> 3.0-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.sparse </td>
   <td style="text-align:center;"> spatstat.sparse </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> 3.0-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.utils </td>
   <td style="text-align:center;"> spatstat.utils </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> 3.0-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-19 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringi </td>
   <td style="text-align:center;"> stringi </td>
   <td style="text-align:center;"> 1.7.8 </td>
   <td style="text-align:center;"> 1.7.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/stringi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/stringi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:center;"> stringr </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> 1.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/stringr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/stringr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-02 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> survival </td>
   <td style="text-align:center;"> survival </td>
   <td style="text-align:center;"> 3.4.0 </td>
   <td style="text-align:center;"> 3.4-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/survival </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/survival </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-09 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> svglite </td>
   <td style="text-align:center;"> svglite </td>
   <td style="text-align:center;"> 2.1.0 </td>
   <td style="text-align:center;"> 2.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/svglite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/svglite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-02-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> systemfonts </td>
   <td style="text-align:center;"> systemfonts </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/systemfonts </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/systemfonts </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-02-11 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tensor </td>
   <td style="text-align:center;"> tensor </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tensor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tensor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2012-05-05 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:center;"> tibble </td>
   <td style="text-align:center;"> 3.1.8 </td>
   <td style="text-align:center;"> 3.1.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tibble </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tibble </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:center;"> tidyr </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tidyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tidyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-08 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyselect </td>
   <td style="text-align:center;"> tidyselect </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tidyselect </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/tidyselect </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-10 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> urlchecker </td>
   <td style="text-align:center;"> urlchecker </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/urlchecker </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/urlchecker </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:center;"> usethis </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/usethis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/usethis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-05-25 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> utf8 </td>
   <td style="text-align:center;"> utf8 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/utf8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/utf8 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-07-24 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> uwot </td>
   <td style="text-align:center;"> uwot </td>
   <td style="text-align:center;"> 0.1.14 </td>
   <td style="text-align:center;"> 0.1.14 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/uwot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/uwot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vctrs </td>
   <td style="text-align:center;"> vctrs </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/vctrs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/vctrs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-16 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridisLite </td>
   <td style="text-align:center;"> viridisLite </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/viridisLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/viridisLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-08-22 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> webshot </td>
   <td style="text-align:center;"> webshot </td>
   <td style="text-align:center;"> 0.5.4 </td>
   <td style="text-align:center;"> 0.5.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/webshot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/webshot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-26 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> withr </td>
   <td style="text-align:center;"> withr </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/withr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/withr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-03 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xfun </td>
   <td style="text-align:center;"> xfun </td>
   <td style="text-align:center;"> 0.36 </td>
   <td style="text-align:center;"> 0.36 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/xfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/xfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-12-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xml2 </td>
   <td style="text-align:center;"> xml2 </td>
   <td style="text-align:center;"> 1.3.3 </td>
   <td style="text-align:center;"> 1.3.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/xml2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/xml2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xtable </td>
   <td style="text-align:center;"> xtable </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/xtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/xtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yaml </td>
   <td style="text-align:center;"> yaml </td>
   <td style="text-align:center;"> 2.3.6 </td>
   <td style="text-align:center;"> 2.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/yaml </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/yaml </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-18 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zoo </td>
   <td style="text-align:center;"> zoo </td>
   <td style="text-align:center;"> 1.8.11 </td>
   <td style="text-align:center;"> 1.8-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/zoo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library/zoo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-09-17 </td>
   <td style="text-align:center;"> CRAN (R 4.2.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.2/Resources/library </td>
  </tr>
</tbody>
</table>
