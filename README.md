---
title: "snRNA-sequencing of Aortic Aneurysm"
author: "Mark E. Pepin, MD, PhD, MS"
date: "2024-10-16"
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



**Code Authors**: Mark E. Pepin, MD, PhD, MS **Contact**: [mepepin\@bwh.harvard.edu](mailto:mepepin@bwh.harvard.edu){.email}\
**Institution**: Brigham and Women's Hospital \| Broad Institute of Harvard and MIT\
**Location**: Boston, MA, USA

# Figure 2 - Integrated snRNA-Seq Analysis

## Integration of Chou et al. with Harmony

Data from Chou et al. 2022 were incorporated to provide a more generalizeable understanding of cellular heterogeneity at single-nuclear resolution. Once the individual snRNA-seq datasets were interrogated for quality and normalized, they were fully integrated using standard protocol within Seurat. To accomplish this, shared cell populations were matched across datasets ("anchors") were used to correct for technical differences between datasets (i.e. batch effects), as well as to perform comparative sn-RNA-seq analysis across experimental conditions (ACTA2-mut vs. TAA vs. CON).


``` r
start_time <- Sys.time()
options(future.globals.maxSize= 891289600000000000)
library(Seurat)
# Load current dataset
# library normalized the cells, log transformed the counts, and scaled the genes
taa_pepin <- readRDS(file = "../1_Input/snRNA_Human_ControlDissection_integrated_seurat_v2.rds")
taa_pepin$Disease<-dplyr::recode_factor(taa_pepin$Disease, Control="CABG Aortic Button", Dissection="aortic aneurysm")
taa_pepin$DataSet <- "Pepin.et.al"
taa_pepin$DonorID <- taa_pepin$Sample_name # rename to match the chou et al paper
taa_pepin$orig.ident <- dplyr::recode_factor(taa_pepin$orig.ident, SeuratProject="taa_Pepin.et.al")
taa_pepin$CellType <- taa_pepin$Celltype_annotation_v2
taa_pepin$Celltype_annotation_v2 <- NULL
taa_pepin$Sample_name <- NULL
taa_pepin$Sample <- NULL
taa_pepin$BroadCellTypes <- NULL
taa_pepin$louvain_05_clusters <- NULL
taa_pepin$leiden_05_c <- NULL
taa_pepin$louvain_05_c <- NULL
taa_pepin$umap_density_Disease <- NULL
taa_pepin$total_counts <- NULL
taa_pepin$n_genes <- NULL
taa_pepin$total_counts_mt <- NULL
taa_pepin$pct_counts_mt <- NULL
taa_pepin$n_genes_by_counts <- NULL
taa_pepin$scorect <- NULL
taa_pepin$n_counts_all <- NULL
# Import the Chou et al. dataset
taa_chou<-readRDS(file = "../1_Input/Chou.TAA_snRNA.rds")
taa_chou$Disease<-factor(taa_chou$Disease, levels = c("normal", "aortic aneurysm"))
taa_chou$DataSet <- "Chou.et.al"
taa_chou$CellType <- factor(taa_chou$CellType)
taa_chou$CellType<-dplyr::recode_factor(taa_chou$CellType,
                                       "08. Endothelial1"="EC1",
                                       "09. Endothelial2"="EC2",
                                       "06. Fibroblast"="Fibroblast1",
                                       "01. VSMC_B"="VSMC_B",
                                       "02. VSMC_A"="VSMC_A",
                                       "04. Pericyte"="Pericyte",
                                       "03. VSMC3"="VSMC3",
                                       "05. Mesothelial"="Mesothelial",
                                       "14. Macrophage"="Macrophage", 
                                       "13. Immune cells"="Immune Cell",
                                       "10. Lymphatic endothelial"="Lymphatic EC",
                                       "11. Neuronal"="Neuronal",
                                       "12. Lymphocyte"="Lymphocyte",
                                       "07. Adipocyte"="Adipocyte")
taa_chou$integrated_snn_res.0.3 <- NULL
taa_chou$seurat_clusters <- NULL
# Merge Seurat Objects
TAA_combined<-merge(taa_pepin, y = taa_chou, project = "TAA", merge.data = TRUE)
remove(taa_chou) # free memory by removing the prior datasets... we're going to need it ;)
remove(taa_pepin)
# Initialize the Seurat object with the raw (non-normalized data).
TAA_combined[["percent.mt"]] <- PercentageFeatureSet(TAA_combined, assay = "RNA", pattern = "^MT-")
TAA_combined <- subset(TAA_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
##################
TAA_list<-SplitObject(TAA_combined, split.by = "orig.ident") # split dataset into a list of two seruat objects
# normalize and identify variable features for each dataset independently
TAA.list <- lapply(X = TAA_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
    x <- ScaleData(x)
    x <- SCTransform(x)
})
# Find most variable features across samples to integrate
integ_features <- SelectIntegrationFeatures(object.list = TAA.list, nfeatures = 3000)
merged_seurat <- merge(x = TAA.list[[1]],
		       y = TAA.list[2:length(TAA.list)],
		       merge.data = TRUE)
DefaultAssay(merged_seurat) <- "SCT"
# Manually set variable features of merged Seurat object
VariableFeatures(merged_seurat) <- integ_features
# Calculate PCs using manually set variable features
merged_seurat <- RunPCA(merged_seurat, assay = "SCT", npcs = 50)
# Harmonize the integrated object
library(harmony)
harmonized_seurat <- RunHarmony(merged_seurat, 
				group.by.vars = c("DonorID", "orig.ident"), 
				reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
harmonized_seurat <- RunUMAP(harmonized_seurat, reduction = "harmony", assay = "SCT", dims = 1:40)
harmonized_seurat <- FindNeighbors(object = harmonized_seurat, reduction = "harmony")
harmonized_seurat <- FindClusters(harmonized_seurat, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
harmonized_seurat <- PrepSCTFindMarkers(harmonized_seurat, assay = "SCT", verbose = TRUE)
markers <- FindAllMarkers(
  object = harmonized_seurat,
  assay = "SCT",
  verbose = F
)
# saveRDS(harmonized_seurat, file = "../1_Input/TAA_Integration_snRNA.rds")
```

## Cell-Type Annotation

Cell types were identified from UMAP-based clustering using a combination of known gene markers and top differentially-expressed genes distiinct to each cellular community. Transcriptional distribution was visualized by overlaying expression with density plots.


``` r
library(Seurat)
library(scCustomize)
library(Nebulosa)
library(ggplot2)
library(ggtrace)
library(ggrepel)
TAA.combined<-readRDS(file = "../1_Input/TAA_Integration_snRNA.rds")
Idents(TAA.combined) <- "SCT_snn_res.0.4" # Change the identity of the clusters to cell types
UMAP_CellTypes<-DimPlot(TAA.combined,  label = T, split.by = "Disease") + NoLegend()
pdf(file = "../2_Output/Figure_2/UMAP_Clusters_snnres_0.2.pdf", height = 3, width = 7)
UMAP_CellTypes
dev.off()
# Unbiased cluster identification
DEGs_Clusters<-FindAllMarkers(TAA.combined, assay = "RNA")
write.csv(DEGs_Clusters, "../2_Output/Figure_2/DEGs_Clusters.csv")
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "lightsalmon"))(paletteLength)
top5_markers <- Extract_Top_Markers(marker_dataframe = DEGs_Clusters, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/Figure_2/DotPlot_Clusters_top5DEGs.pdf", height = 10, width = 7)
Clustered_DotPlot(seurat_object = TAA.combined, 
                  features = top5_markers, 
                  k = 10, 
                  colors_use_exp = myColor,
                  colors_use_idents = NA,
                  cluster_ident = F)
dev.off()
top30_markers <- Extract_Top_Markers(marker_dataframe = DEGs_Clusters, num_genes = 5, named_vector = FALSE,
    make_unique = TRUE)
pdf(file = "../2_Output/Figure_2/Heatmap_Clusters.pdf", height = 15, width = 10)
DoHeatmap(TAA.combined, features = top30_markers, size = 3, disp.min = -2, disp.max = 2) + scale_fill_gradientn(colors = c("dodgerblue4", "white", "lightsalmon"))
dev.off()
# Cluster Identification using gene markers
VSMC_TEST <- c("ITGA8", "SMTN")
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1", "VCAM1")
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c("NFASC",  "SMAD5", "PRSS23") #"UACA","TEX41",
Fibroblast_C_genes<-c("ADAMTS1", "RGS6", "TNC", "ANGPT2", "DGKG", "GRIP2", "LUM") # 
# EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") # , "NOTCH4"
# EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
EC1_genes <- c("VCAM1", "ELN", "CLU", "CYTL1", "BGN")
EC2_genes <- c("FLT1", "FABP5", "KDR", "SPARCL1")
##
NKT_genes<-c("SKAP1", "RIPOR2", "ITGAL", "CD96") #  "RBPJ", "FYN",
Macrophage_genes<-c("MRC1", "LGMN", "F13A1", "RBM47") #
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8") #
Neuronal_genes <- c("NRXN1", "CADM2", "SORCS1", "ARHGAP15")
Pericyte_genes <- c("VTN", "HIGD1B", "S1PR3", "MCAM", "IFITM1", "BAIAP3", "EHD3")
# Plot density function
VSMC_Test<-plot_density(TAA.combined, VSMC_TEST, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC1_density<-plot_density(TAA.combined, VSMC1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(TAA.combined, VSMC2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast1_density<-plot_density(TAA.combined, Fibroblast1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast2_density<-plot_density(TAA.combined, Fibroblast2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast_C_density<-plot_density(TAA.combined, Fibroblast_C_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
EC1_density<-plot_density(TAA.combined, EC1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
EC2_density<-plot_density(TAA.combined, EC2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Macrophage_density<-plot_density(TAA.combined, Macrophage_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
NKT_density<-plot_density(TAA.combined, NKT_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Dendritic_density<-plot_density(TAA.combined, Dendritic_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Neuronal_density<-plot_density(TAA.combined, Neuronal_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Pericyte_density<-plot_density(TAA.combined, Neuronal_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
pdf(file = "../2_Output/Figure_2/Density_plots.pdf", height = 3, width = 3)
VSMC_TEST
VSMC1_density
VSMC2_density
Fibroblast1_density
Fibroblast2_density
Fibroblast_C_density
EC1_density
EC2_density
Macrophage_density
NKT_density
Dendritic_density
dev.off()
# Overlay these gene markers onto the UMAP to identify clusters
pdf(file = "../2_Output/Figure_2/CellType_FeaturePlots.pdf")
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = VSMC1_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = VSMC2_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Fibroblast_C_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Fibroblast1_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Fibroblast2_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = EC1_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = EC2_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Macrophage_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = NKT_genes)
FeaturePlot(TAA.combined, reduction = "umap", label = T, features = Dendritic_genes)
dev.off()
#Create a figure of cell-type specific markers overlying the UMAP
pdf(file = "../2_Output/Figure_2/CellType_Differentiation.pdf")
plot_density(TAA.combined, c("ACTA2", "FBLN1", "F13A1", "VWF", "STC1", "NFASC", "ITGAL", "ITGAX"), reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
dev.off()
## Use ridgeplots to identify bimodal gene marker distributions (enriched clusters)
pdf(file = "../2_Output/Figure_2/Celltype_RidgePlots.pdf", height = 10, width = 15)
RidgePlot(TAA.combined,features = VSMC_TEST, ncol = 2)
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
# saveRDS(TAA.combined, file = "../1_Input/TAA_Clustering_snRNA.rds")
```

## Labelling Cells

Based on the cell identification above, clusters were renamed and saved for further sub-cluster analysis.


``` r
library(Seurat)
TAA.combined <- readRDS(file = "../1_Input/TAA_Clustering_snRNA.rds")
TAA.combined <- RenameIdents(TAA.combined,
             `0` = "VSMC_A",
             `1` = "VSMC_B",
             `2` = "Fibroblast_A",
             `3` = "VSMC_C",
             `4` = "EC",
             `5` = "Fibroblast_B",
             `6` = "Macrophage",
             `7` = "Fibroblast_C",
             `8` = "NKT",
             `9` = "Neuronal",
             `10` = "VSMC_B")
TAA.combined$CellTypev2 <- TAA.combined@active.ident
TAA.combined$CellTypev2 <- factor(TAA.combined$CellTypev2, levels = c("VSMC_A", 
                                                                      "VSMC_B", 
                                                                      "VSMC_C", 
                                                                      "Fibroblast_A", 
                                                                      "Fibroblast_B", 
                                                                      "Fibroblast_C",
                                                                      "EC", 
                                                                      "Macrophage", 
                                                                      "NKT", 
                                                                      "Neuronal"))
TAA.combined <- SetIdent(TAA.combined, value = "CellTypev2")
UMAP_CellTypes<-DimPlot(TAA.combined,  label = T, repel = T, label.size = 4,
                        cols = c("lightsalmon", 
                                 "firebrick4", 
                                 "coral2",
                                 "lightblue3", 
                                 "steelblue4",
                                 "steelblue",
                                 "goldenrod2",
                                 "azure4",
                                 "darkcyan",
                                 "seagreen3")) + NoLegend()
pdf(file = "../2_Output/Figure_2/Figure_2E_UMAP_Cell.Clusters.pdf", height = 5, width = 7)
UMAP_CellTypes
dev.off()
# Cell-type Specific Differential Expression
DEGs_CellTypes<-FindAllMarkers(TAA.combined, assay = "RNA")
openxlsx::write.xlsx(DEGs_CellTypes, "../2_Output/Figure_2/DEGs_CellTypes.xlsx")
DEGs_VSMCs <- FindMarkers(TAA.combined, assay = "RNA", ident.1 = "VSMC_B", ident.2 = "VSMC_A")
openxlsx::write.xlsx(DEGs_VSMCs, "../2_Output/Figure_3/DEGs_VSMC_BvA.xlsx")
# Save file
saveRDS(TAA.combined, file = "../1_Input/TAA_Labelling_snRNA.rds")

# Proportion of cells expressing specific genes
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
PrctCellExpringGene <- function(object, genes, group.by = "all"){
    if(group.by == "all"){
        prct = unlist(lapply(genes,calc_helper, object=object))
        result = data.frame(Markers = genes, Cell_proportion = prct)
        return(result)
    }

    else{        
        list = SplitObject(object, group.by)
        factors = names(list)

        results = lapply(list, PrctCellExpringGene, genes=genes)
        for(i in 1:length(factors)){
        results[[i]]$Feature = factors[i]
        }
        combined = do.call("rbind", results)
        return(combined)
    }
}

calc_helper <- function(object,genes){
    counts = object[['RNA']]@counts
    ncells = ncol(counts)
    if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
    }else{return(NA)}
}
library(Seurat)
PrctCellExpringGene(TAA.combined ,genes =c("ITGA8", "SMTN"), group.by = "CellTypev2")
TAA_VSM <- subset(TAA.combined, subset = CellTypev2 == c("VSMC_A", "VSMC_B", "VSMC_C", "Fibroblast_C"))
VlnPlot(TAA_VSM, features = c("ITGA8", "SMTN"), split.by = "CellTypev2", group.by = "Disease")
FeaturePlot(TAA_VSM, features = c("ITGA8", "SMTN"), split.by = "Disease",)
```

## UMAP and Cellular Proportions

Relative composition differences were identified from the single-nuclear RNA-sequencing data to better understand the likely culprit cell-type in the pathophysiology of aortic aneurysms and atherosclerosis.


``` r
library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined$CellTypev2 <- factor(TAA.combined$CellTypev2, levels = c("VSMC_A",
                                                                      "VSMC_B",
                                                                      "VSMC_C",
                                                                      "Fibroblast_A",
                                                                      "Fibroblast_B",
                                                                      "Fibroblast_C",
                                                                      "EC",
                                                                      "Macrophage",
                                                                      "NKT",
                                                                      "Neuronal"))
TAA.combined <- SetIdent(TAA.combined, value = "CellTypev2")
# Show Integration Success between cohorts/datasets
p1 <- DimPlot(TAA.combined, reduction = "umap", group.by = "Disease", split.by = "orig.ident") + ggtitle(NULL)# Chou vs. Pepin
p1
pdf(file = "../2_Output/Figure_2/Figure_2A.pdf", height = 4, width = 7)
p1
# Show distribution of clustering across diseases
dev.off()
p2 <- DimPlot(TAA.combined, reduction = "umap", group.by = "SCT_snn_res.0.4", split.by = "Disease", label = T, repel = TRUE,
              cols = c("lightsalmon", 
                       "firebrick4",
                       "steelblue",
                       "coral2",
                       "goldenrod2",
                       "steelblue4",
                       "azure4",
                       "lightblue3",
                       "darkcyan",
                       "seagreen3",
                       "lightsalmon")) + NoLegend() + ggtitle(NULL)
p2
pdf(file = "../2_Output/Figure_2/Figure_2B.pdf", height = 4, width = 7)
p2
dev.off()
p3 <- DimPlot(TAA.combined, reduction = "umap", group.by = "CellTypev2", split.by = "Disease", label = T, repel = TRUE,
              cols = c("lightsalmon", 
                       "firebrick4", 
                       "coral2",
                       "lightblue3", 
                       "steelblue",
                       "steelblue4",
                       "goldenrod2",
                       "azure4",
                       "darkcyan",
                       "seagreen3",
                       "lightsalmon")) + NoLegend() + ggtitle(NULL) # TAA vs. CABG vs. Control
p3

#Figure 2C -  Proportional Graph
library(dittoSeq)
pdf(file = "../2_Output/Figure_2/Fig2C_Proportional.Bar_Disease.pdf")
dittoBarPlot(
    object = TAA.combined,
    var = "CellTypev2",
    group.by = "Disease")+ ggtitle(NULL)
dev.off()
DittoPLOT <- dittoBarPlot(
    object = TAA.combined,
    var = "CellTypev2",
    group.by = "Disease",
    data.out = T)
write.csv(DittoPLOT$data, "../2_Output/Figure_2/Proportional.Cells.csv")
BarPlot <- dittoBarPlot(
    object = TAA.combined,
    var = "CellTypev2",
    group.by = "DonorID",
    split.by = "Disease", 
    data.out = T)
# Sample-specific proportional plot
pdf(file = "../2_Output/Figure_2/Fig2C_Proportional.Bar_Samples.pdf", height = 4, width = 6)
dittoBarPlot(
    object = TAA.combined,
    var = "CellTypev2",
    group.by = "DonorID",
    split.by = "Disease",
    retain.factor.levels = T,
    color.panel = c(VSMC_A = "lightsalmon",
                    VSMC_B = "firebrick4",
                     VSMC_C = "coral2",
                     Fibroblast_A = "steelblue",
                     Fibroblast_B = "steelblue4",
                     Fibroblast_C = "lightblue3",
                     EC = "goldenrod2",
                     Macrophage = "azure4",
                     NKT = "darkcyan",
                     Neuronal = "seagreen3")
    ) + 
  facet_wrap(~factor(Disease, levels = c("normal", "aortic aneurysm", "CABG Aortic Button")), 
             ncol = 3, 
             scales = "free_x") +
  ggtitle(NULL) + 
  labs(x = NULL)
dev.off()
# Extract Data for export
BarPlot_sampledata <- BarPlot$data %>% 
  select(-count, -label.count.total.per.facet) %>% 
  tidyr::pivot_wider(., names_from = "grouping", values_from = "percent")
openxlsx::write.xlsx(BarPlot_sampledata, "../2_Output/Figure_2/BarPlot_CellTypes.xlsx")
# Visualize the features/genes
VizDimLoadings(TAA.combined, dims = 1:5, reduction = "pca")
pdf(file = "../2_Output/Figure_2/PC_Heatmaps.pdf")
DimHeatmap(TAA.combined, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
```

# Figure 3: VSMC-Specific Analysis

From the observation that VSMC_B are differentially enriched in both aneurysmal and atherosclerotic aortic tissue, we sought to better understand the transcriptional differences between VSMC_B and the other VSMC sub-clusters.

### VSMC Phenotypic Comparison

The first step was to identify DEGs between VSMC_B and VSMC_A, which was accomplished as follows.


``` r
library(Seurat)
library(dplyr)
library(ggtrace)
library(ggplot2)
library(ggrepel)
# Add tracings
library(ggplot2)
library(ggtrace)
library(ggthemes)
library(ggrepel)
library(dplyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
# TAA.combined[["SCT"]]@SCTModel.list <- TAA.combined[["SCT"]]@SCTModel.list[1]
TAA_normal <- subset(TAA.combined, subset = Disease == "normal") # Use only the control samples for VSMC-specific analysis
# TAA_normal <- PrepSCTFindMarkers(TAA_normal, assay = "SCT", verbose = TRUE)
# markers <- FindAllMarkers(
#   object = TAA_normal,
#   assay = "SCT",
#   verbose = T
# )
# # Differential Expression between VSMCs
VSMC_A_DEGs <- FindMarkers(TAA_normal, ident.1 = "VSMC_B", ident.2 = "VSMC_A", assay = "SCT", min.pct = 0.1, logfc.threshold = 0)
VSMC_A_DEGs[VSMC_A_DEGs$p_val==0,"p_val"] <- 1e-300 # Change P-value of "0" to 1e-300
openxlsx::write.xlsx(VSMC_A_DEGs, "../2_Output/Figure_3/DEGs_VSMC_A.vs.B.xlsx", rowNames = T)
VSMC_A_DEGs <- openxlsx::read.xlsx("../2_Output/Figure_3/DEGs_VSMC_A.vs.B.xlsx", rowNames = T) %>% filter(pct.1>0.1, pct.2>0.1)
# Select only the central cluster
hdat_vsmc <- subset(TAA_normal, CellTypev2==c("VSMC_B", "VSMC_A", "VSMC_C")) #, "Fibroblast_C"
#Re-cluster
hdat_vsmc <- RunUMAP(hdat_vsmc, reduction = "harmony", assay = "SCT", dims = 1:40)
hdat_vsmc <- FindNeighbors(object = hdat_vsmc, reduction = "harmony")
hdat_vsmc <- FindClusters(hdat_vsmc, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
#
hdat_vsmc<-SetIdent(hdat_vsmc, value = "SCT_snn_res.0.2")
hdat_vsmc <- RenameIdents(hdat_vsmc, 
                          `0` = "VSMC_A",
                          `2` = "VSMC_B",
                          `1` = "VSMC_C")
UMAP_VSMC<-DimPlot(hdat_vsmc,  label = T,
                   cols = c("lightsalmon", 
                       "coral2", 
                       "firebrick4")) + NoLegend()
UMAP_VSMC
# Figure 1A - UMAP of VSMCs only
pdf("../2_Output/Figure_3/Fig3A_UMAP_VSMC.pdf", height = 3.5, width = 3.5)
UMAP_VSMC
dev.off()
# Upload Gene Markers
VSMC1_genes<-c("PKD1", "COL6A2", "PDLIM7", "FLNA", "SMTN")
VSMC2_genes<-c("MYH11", "ACTA2", "ITGA8", "PRKG1", "CRISPLD1")
VSMC_gene <- c("SERPINE1", "FN1")
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c("NFASC",  "SAMD5", "PRSS23") #"UACA","TEX41",
Fibroblast_C_genes<-c("ADAMTS1", "RGS6", "TNC", "ANGPT2", "DGKG", "GRIP2") # 
EC1_genes<-c("DIPK2B", "ARHGEF15", "STC1", "FLT1") # , "NOTCH4"
EC2_genes<-c("VWF", "BMPER", "BMX", "NOS1")
NKT_genes<-c("SKAP1", "RIPOR2", "ITGAL", "CD96") #  "RBPJ", "FYN",
Macrophage_genes<-c("MRC1", "LGMN", "F13A1", "RBM47") #
Dendritic_genes<-c("ITGAX", "S100A9", "CSF3R", "CXCL8") #
Neuronal_genes <- c("NRXN1", "CADM2", "SORCS1", "ARHGAP15")
Pericyte_genes <- c("VTN", "HIGD1B", "S1PR3", "MCAM", "IFITM1", "BAIAP3", "EHD3")
# Create VSMC Subset UMAP Density plot of Gene Markers
library(Nebulosa)
VSMC1_density<-plot_density(hdat_vsmc, VSMC1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC2_density<-plot_density(hdat_vsmc, VSMC2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC_density<-plot_density(hdat_vsmc, VSMC_gene, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
VSMC_density
pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Density_plots.pdf", height = 3.5, width = 3.5)
VSMC1_density
VSMC2_density
dev.off()
#Volcano Plot
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
# Read data from the web
options(ggrepel.max.overlaps = Inf)
results = mutate(VSMC_A_DEGs, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.0001 & log2FC>1, 
                       "P < 0.0001 and Log(Fold-Change) > 1", 
                       ifelse(p_val<0.0001 & log2FC< 0-1,
                              "P < 0.0001 and Log(Fold-Change) < -1", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -1",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > 1")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  # theme(axis.line = element_blank(),
  #       axis.ticks = element_blank()
  # ) +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=1, size=.5, linetype="dashed") +
  geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 10, -log2FC),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 10, log2FC), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
#
pdf(file = paste0("../2_Output/Figure_3/VSMC_B.vs.A_VolcanoPlot.pdf"), height = 4, width = 4)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< -1), 15,  -log2FC), aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 15, log2FC), aes(label=gene_name)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
#########
VSMC_type_DEs <- results %>% filter(p_val < 0.05, abs(log2FC)>1) %>% dplyr::select(gene_name)
VSMC_type_DEs <- VSMC_type_DEs$gene_name
##
plots_VSMC_A <- VlnPlot(hdat_vsmc, features = c("SERPINE1", "MALAT1", "PKD1"), group.by = "CellTypev2",
    pt.size = 0,
    cols = c("darkcyan", "lightsalmon", "black", "grey", "dodgerblue3")) & theme(legend.position = 'none', axis.title.x = element_blank())
plots_VSMC_B <- VlnPlot(hdat_vsmc, features = c("PDE4D", "PDE3A", "LPP"), group.by = "CellTypev2",
    pt.size = 0,
    cols = c("darkcyan", "lightsalmon", "black", "grey", "dodgerblue3")) & theme(legend.position = 'none', axis.title.x = element_blank())

pdf(file = "../2_Output/Figure_3/VSMC_A_Top.Violins.pdf", height = 3, width = 6)
plots_VSMC_A
dev.off()
pdf(file = "../2_Output/Figure_3/VSMC_B_Top.Violins.pdf", height = 3, width = 6)
plots_VSMC_B
dev.off()

########## Volcano Genes
library(scCustomize)
pdf(file = "../2_Output/Figure_3/Volcano_Genes.pdf", height = 5, width = 5)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod2"))(paletteLength)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = VSMC_type_DEs,
                  colors_use_idents = c("lightsalmon","coral2","firebrick4"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = T)
dev.off()
########## Contractility Genes
library(scCustomize)
contractility_genes <- c("ACTA2", "CNN1", "MYH11", "SMTN", "CARMN", "TGFB1I1", "SRF", "CALD1", "LMOD1", "KCNMB1", "SYNM", "ROCK2")
pdf(file = "../2_Output/Figure_3/Contractility_Genes.pdf", height = 5, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = contractility_genes,
                  colors_use_idents = c("lightsalmon","coral2","firebrick4"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = T)
dev.off()

library(scCustomize)
pdf(file = "../2_Output/Figure_3/DotPlot_VSMC_1.v.2.pdf", height = 10, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, features = VSMC_type_DEs,  x_lab_rotate=F, k = 3)
dev.off()
hdat_vsmc <- BuildClusterTree(object = hdat_vsmc)
PlotClusterTree(hdat_vsmc)
## Export 
saveRDS(hdat_vsmc, file = "../1_Input/TAA_VSMCs_snRNA.rds")
```

### VSMC GSEA (VSMC Subtypes)

We then used GO-term enrichment to identify the related pathways and gene-sets that are disproportionately represented among DEGs in aneurysmal and/or atherosclerotic aortas.


``` r
##Enrichr
# library(enrichR)
library(dplyr)
VSMC_A_DEGs <-openxlsx::read.xlsx("../2_Output/Figure_3/DEGs_VSMC_A.vs.B.xlsx", rowNames = T) %>% filter(p_val < 0.05)

## EnrichR Pathway Analysis
VSMC_A_UP <- VSMC_A_DEGs %>% filter(avg_log2FC > 0) %>% top_n(250, avg_log2FC)
VSMC_A_DOWN <- VSMC_A_DEGs %>% filter(avg_log2FC < -0) %>% top_n(250, -avg_log2FC)
##Enrichr
library(enrichR)
library("kableExtra")
library("knitr")
dbs <- c("WikiPathway_2023_Human")
enriched_UP <- enrichr(rownames(VSMC_A_UP), dbs)
enrich_UP<-enriched_UP[[dbs]] %>% filter(Adjusted.P.value < 0.05)
# Enrich DOWN
enriched_DOWN <- enrichr(rownames(VSMC_A_DOWN), dbs)
enrich_DOWN<-enriched_DOWN[[dbs]] %>% filter(Adjusted.P.value < 0.05)
write.csv(enrich_DOWN, "../2_Output/Figure_3/Fig3_PathwayDEGs_DOWN.csv")
write.csv(enrich_UP, "../2_Output/Figure_3/Fig3_PathwayDEGs_UP.csv")
DEGS_UP<-strsplit(enrich_UP$Genes[1], split = ";")[[1]]
DEGS_DOWN<-strsplit(enrich_DOWN$Genes[1], split = ";")[[1]]
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
hdat_vsmc <- subset(TAA.combined, CellTypev2==c("VSMC_B", "VSMC_A", "VSMC_C"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Pathway_DEGs.UP.pdf", height = 5, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = DEGS_UP,
                  colors_use_idents = c("lightsalmon","coral2", "firebrick4"),
                  k = 2,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Pathway_DEGs.DOWN.pdf", height = 5, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = DEGS_DOWN,
                  colors_use_idents = c("lightsalmon","coral2", "firebrick4"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
```

#### Upstream Regulator Identification


``` r
##Enrichr
library(dplyr)
VSMC_A_DEGs <-openxlsx::read.xlsx("../2_Output/Figure_3/DEGs_VSMC_A.vs.B.xlsx", rowNames = T) %>% filter(p_val < 0.05)

## EnrichR Pathway Analysis
VSMC_A_UP <- VSMC_A_DEGs %>% filter(avg_log2FC > 0) %>% top_n(1000, avg_log2FC)
VSMC_A_DOWN <- VSMC_A_DEGs %>% filter(avg_log2FC < 0) %>% top_n(1000, -avg_log2FC)

##Enrichr
library(enrichR)
dbs <- c("ChEA_2022")
enriched_UP <- enrichr(rownames(VSMC_A_UP), dbs)
enrich_UP<-enriched_UP[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_UP)
enriched_DOWN <- enrichr(rownames(VSMC_A_DOWN), dbs)
enrich_DOWN<-enriched_DOWN[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_DOWN)
write.csv(enrich_DOWN, "../2_Output/Figure_3/Fig3_Regulator_DEGs_DOWN.csv")
write.csv(enrich_UP, "../2_Output/Figure_3/Fig3_Regulator_DEGs_UP.csv")
DEGS_UP<-strsplit(enrich_UP$Genes[1], split = ";")[[1]]
DEGS_DOWN<-strsplit(enrich_DOWN$Genes[1], split = ";")[[1]]
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
hdat_vsmc <- subset(TAA.combined, CellTypev2==c("VSMC_B", "VSMC_A", "VSMC_C"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod2"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Regulator_DEGs.UP.pdf", height = 10, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = DEGS_UP,
                  colors_use_idents = c("lightsalmon","firebrick4","coral2","lightblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Regulator_DEGs.DOWN.pdf", height = 10, width = 5)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = DEGS_DOWN,
                  colors_use_idents = c("lightsalmon","firebrick4","coral2","lightblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
```

#### TCF21 Targets


``` r
library(dplyr)
library(openxlsx)
enrich_UP <- read.csv("../2_Output/Figure_3/Fig3_Regulator_DEGs_UP.csv")
DEGS_UP<-strsplit(enrich_UP$Genes[1], split = ";")[[1]]
VSMC_A_DEGs <-openxlsx::read.xlsx("../2_Output/Figure_3/DEGs_VSMC_A.vs.B.xlsx", rowNames = T) %>% filter(p_val_adj < 0.05)
VSMC_A_DEGs$GeneSymbol <- rownames(VSMC_A_DEGs)
VSMC_TCF21.Targets<-VSMC_A_DEGs %>% subset(., GeneSymbol %in% DEGS_UP)
######
library(biomaRt)
ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id_version", "external_gene_name", "chromosome_name", "start_position", "end_position"),  mart=ensembl)
write.csv(bm, "../1_Input/BiomaRt_Annotation.csv")
## Annotate DEGs with chromosomal location
VSMC_A.Sites <- merge(VSMC_A_DEGs, bm, by.x = "GeneSymbol", by.y = "external_gene_name") %>% rename(Chr=chromosome_name, Start = start_position, End = end_position)
VSMC_A.Sites$Chr<-paste0("chr",VSMC_A.Sites$Chr)
VSMC_A.Sites$Chr<-factor(VSMC_A.Sites$Chr, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chrX", "chrY"))
VSMC_A.Sites <- VSMC_A.Sites %>% na.omit()
## Annotate TCF21 Targets
Annotated_TCF21.Sites <- merge(VSMC_TCF21.Targets, bm, by.x = "GeneSymbol", by.y = "external_gene_name") %>% rename(Chr=chromosome_name, Start = start_position, End = end_position)
Annotated_TCF21.Sites$Chr<-paste0("chr",Annotated_TCF21.Sites$Chr)
Annotated_TCF21.Sites$Chr<-factor(Annotated_TCF21.Sites$Chr, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chrX", "chrY"))
Annotated_TCF21.Sites <- Annotated_TCF21.Sites %>% na.omit()
Gene_labels<-as.data.frame(Annotated_TCF21.Sites) %>% dplyr::select(chrom=Chr, chromStart=Start, chromEnd=End, GeneSymbol) %>% filter(GeneSymbol!="")
Gene_labels<-arrange(Gene_labels, chromStart)
Gene_labels$chrom<-factor(Gene_labels$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",  "chr13", "chr14", "chr15", "chr16",  "chr17", "chr18", "chr19", "chr20",  "chr21", "chr22", "chr23", "chrX", "chrY"))
Gene_labels<-Gene_labels[order(Gene_labels$chrom),]
Gene_labels<-Gene_labels[!duplicated(Gene_labels[,4]),] %>% tidyr::drop_na()
##Fold Change UP
Gene_FoldChange.UP<-dplyr::filter(VSMC_A.Sites, avg_log2FC>0) %>% dplyr::select(chrom=Chr, chromStart=Start, FoldChange_DEG=avg_log2FC) %>% dplyr::mutate(chromEnd=chromStart+1) %>% dplyr::select(chrom, chromStart, chromEnd, FoldChange_DEG) %>% distinct() %>% dplyr::arrange(chromStart)
Gene_FoldChange.UP$chrom<-factor(Gene_FoldChange.UP$chrom, levels=c("chr1", "chr2", "chr3", "chr4","chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chrX", "chrY"))
Gene_FoldChange.UP<-Gene_FoldChange.UP[order(Gene_FoldChange.UP$chrom),]
##Fold Change DOWN
Gene_FoldChange.DOWN<-dplyr::filter(VSMC_A.Sites, avg_log2FC<0) %>% dplyr::select(chrom=Chr, chromStart=Start, FoldChange_DEG=avg_log2FC) %>% dplyr::mutate(chromEnd=chromStart+1) %>% dplyr::select(chrom, chromStart, chromEnd, FoldChange_DEG) %>% distinct() %>% dplyr::arrange(chromStart)
Gene_FoldChange.DOWN$chrom<-factor(Gene_FoldChange.DOWN$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chr23", "chrX", "chrY"))
Gene_FoldChange.DOWN<-Gene_FoldChange.DOWN[order(Gene_FoldChange.DOWN$chrom),]
##Fold Change List
Gene_FoldChange_List<-list(Gene_FoldChange.UP, Gene_FoldChange.DOWN)
#Plot the Circos
library(circlize)
library(gtools)
library(dplyr)
circos.genomicDensity1 = function (data, ylim.force = FALSE, window.size = NULL, overlap = TRUE, col = ifelse(area, "grey", "black"), lwd = par("lwd"), lty = par("lty"), type = "l", area = TRUE, area.baseline = NULL, baseline = 0, border = NA, ...) { if (!is.null(area.baseline)) 
data = normalizeToDataFrame(data)
if (!is.dataFrameList(data)) {
data = list(data)
}
if (length(col) == 1) {
col = rep(col, length(data))
}
if (length(lwd) == 1) {
lwd = rep(lwd, length(data))
}
if (length(lty) == 1) {
lty = rep(lty, length(data))
}
if (length(type) == 1) {
type = rep(type, length(data))
}
if (length(area) == 1) {
area = rep(area, length(data))
}
if (length(baseline) == 1) {
    baseline = rep(baseline, length(data))
}
if (length(border) == 1) {
    border = rep(border, length(data))
}
s = sapply(get.all.sector.index(), function(si) get.cell.meta.data("xrange", 
    sector.index = si))
if (is.null(window.size)) {
    window.size = 10^nchar(sum(s))/1000
}
df = vector("list", length = length(data))
for (i in seq_along(data)) {
    all.chr = unique(data[[i]][[1]])
    for (chr in all.chr) {
        region = data[[i]][data[[i]][[1]] == chr, 2:3, drop = FALSE]
        dn = genomicDensity(region, window.size = window.size, 
            overlap = overlap)
        dn = cbind(rep(chr, nrow(dn)), dn)
        df[[i]] = rbind(df[[i]], dn)
    }
}
if (ylim.force) {
    ymax = 1
}
else {
    ymax = max(sapply(df, function(gr) max(gr[[4]])))
}
circos.genomicTrackPlotRegion(df, ylim = c(-ymax,0), panel.fun = function(region, 
    value, ...) {
    i = getI(...)

    circos.genomicLines(region, -value, col = col[i], lwd = lwd[i], 
        lty = lty[i], type = type[i], border = border[i], 
        area = area[i], baseline = baseline[i])
}, ...)
}

write.xlsx(Annotated_TCF21.Sites, "../2_Output/Figure_3/Circos_TCF21_Target.DEGs.xlsx",  overwrite = TRUE)

# Circlize
om = circos.par("track.margin")
oc = circos.par("cell.padding")
circos.par(track.margin = c(0, 0), cell.padding = c(0, 0, 0, 0))
circos.par(start.degree = -250)
pdf(file=paste0("../2_Output/Circos.pdf"))
circos.initializeWithIdeogram(track.height = 0.05, species = "hg19")
### Labels for inversely changing DMRs with DEG
# circos.genomicDensity(Methyl.UP, col = c("lightsalmon"), track.height = 0.1, baseline="bottom", bg.border ="white", track.margin = c(0, 0.0))
# # circos.genomicDensity1(Methyl.DOWN, col = c("darkcyan"), track.height = 0.1, baseline="top", bg.border ="white", track.margin = c(0, 0.0))

circos.genomicLabels(Gene_labels, labels.column=4, side='outside', cex=0.6)

##DEGs as putative downstream targets of given transcriptional regulator
circos.genomicTrackPlotRegion(Gene_FoldChange_List,
                              ylim = c(-3, 3), bg.border=NA,
                              panel.fun = function(region, value, ...) {
 col = ifelse(value[[1]] > 0, "darkgoldenrod1", "darkcyan")
 circos.genomicPoints(region, value, col = add_transparency(col, 0.2), cex = 0.8, pch = 16)
 cell.xlim = get.cell.meta.data("cell.xlim")
 for(h in c(-2, -1, 0, 1, 2)) {
   circos.lines(cell.xlim, c(h, h), col ="#00000040")
 }
}, track.height = 0.2)
# circos.par(track.margin=om, cell.padding=oc)
## Add link for all DEGs
Link_Anchor <- read.csv("../1_Input/Circos/Link_Anchor.csv")
# Link<-Gene_labels %>% dplyr::select(chrom=Chr, chromStart=Start, chromEnd=End) %>% distinct()
Link<-arrange(Gene_labels, chromStart)
Link$chrom<-factor(Link$chrom, levels=c("chr1", "chr2", "chr3", "chr4", 
                                                      "chr5", "chr6", "chr7", "chr8", 
                                                      "chr9", "chr10", "chr11", "chr12", 
                                                      "chr13", "chr14", "chr15", "chr16", 
                                                      "chr17", "chr18", "chr19", "chr20", 
                                                      "chr21", "chr22", "chr23", "chrX", 
                                                      "chrY"))
Link<-Link[order(Link$chrom),]
# Link<-Link[!duplicated(Link[,4]),]
Link_Anchor<-Link_Anchor %>% slice(rep(1, each = nrow(Link))) #create exact number of anchors as link ends.
Link_Anchor$chromStart<-as.numeric(Link_Anchor$chromStart)
Link_Anchor$chromEnd<-as.numeric(Link_Anchor$chromEnd)
Link$chromStart<-as.numeric(Link$chromStart)
Link$chromEnd<-as.numeric(Link$chromEnd)
circos.genomicLink(Link, Link_Anchor, col = add_transparency("#00000040", 0.8), lwd=0.5)
circos.clear()
dev.off()
```

### VSMC Trajectory Analysis (Monocle 3)

Trajectory analysis was employed to characterize the phenotypic shifts between VSMC_A and VSMCB among non-diseased aortic control tissue to better understand how VSMCs might transition from one functional state to another over time or in response to specific pathologic contexts.


``` r
options(future.globals.maxSize= 89128960000000)
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
TCF21_GENES <- c("ERRFI1", "SNED1", "PDXK", "COL16A1", "HSP90AB1", "C1R", "SERPINE1", "LTBP4", "LTBP2", "TFPI", "CABIN1", "ARHGAP22", "SLC22A15", "UBC", "SVEP1", "CCNL1", "DYNC2H1",  "PDGFRB", "PRPF38B", "COL27A1",  "TSC22D3", "STAT2", "CIRBP", "AP3D1",  "NAV1", "BICC1", "COL1A1", "SPSB1", "COL4A4",  "COL6A1", "ITGBL1", "ULK1", "PLA2R1", "RAD9A",  "PLEC")
TCF21_GENES <- c('A1CF','ADCY7','AFG3L2','AICDA','AKAP4','AKR7A3','ALDOA','ALG2','AQP5','ARL1','ARL9','BRSK1','BRWD1','C4B','CA2','CAMKK1','CBX3','CD3EAP','CDX2','CENPT','CLDN22','CMTM1','CNGA1','CRYGN','CTU2','DDAH2','DDX4','DHODH','DHRS7B','DNAJC30','EHMT2','EIF4A2','ELMO3','F7','FAM3B','FCMR','FERMT3','FETUB','GYS1','HIBCH','HNRNPA2B1','HNRNPA3','HRAS','HSPBP1','IL1RAP','IL1RAPL1','INS','ITIH4','ITPRIPL1','KLF6','KLK1','KRTAP13-2','LHX6','LIG1','LRRC56','MADCAM1','MEOX1','MRPL17','MTF1','NAT8','NID2','NPRL2','OTUD5','PASK','PBX3','PDE2A','PEAR1','PFKP','PHOX2A','PKDCC','PLA2G6','POLR2M','PPP1R7','PPP4C','PSMD8','RAB3C','RAB43','RECK','RNF166','RTF1','RUVBL2','S100A3','S100A4','SCAF8','SCX','SEC61B','SGSM1','SKP1','SLC22A7','SLC38A7','SLC51A','SPSB1','STK3','TAF6','TAS2R41','TECTB','TGFA','TGM2','THAP11','TIMM8B','TIPIN','TM9SF2','TMEM119','TMEM95','TMX4','TRAPPC6A','TSHR','TSPAN33','VAV2','VSIG8','WBSCR22')
Merged_scaled <- readRDS(file = "../1_Input/TAA_VSMCs_snRNA.rds")
# Merged_scaled <- subset(TAA.combined, subset = Disease == c("normal") & CellTypev2 == c("VSMC_A", "VSMC_B", "VSMC_C")) # Select only "normal" samples and VSMCs
# Merged_scaled<-TAA.combined
Merged_scaled@active.assay = "RNA"
cds<-SeuratWrappers::as.cell_data_set(Merged_scaled)
# cds <- preprocess_cds(cds, num_dim = 100) # creating errors!
cds <- estimate_size_factors(cds)
# Include gene names (not done by default by the seurat conversion)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
#  Run Monocle
cds <- cluster_cells(cds) #This step creates "partitions" that are used in the trajectory inference
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition") # this shows the partitions overlain on the UMAP
cds <- learn_graph(cds, use_partition = TRUE) # creating trajectory inference within each partition
# cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
root_cell <- "GGGATCCTCCGGACGT-1-4_2" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)
# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_branch_points = FALSE,
           label_leaves = FALSE)
pdf(file = "../2_Output/Figure_3/Fig3_VSMC_UMAP_Trajectory_Partition.pdf", height = 3, width = 4)
plot_cells(cds,
           color_cells_by = "partition",
           show_trajectory_graph = F,
           graph_label_size = 1,
           cell_size = .5,
           label_roots = F,
           label_branch_points = FALSE,
           label_leaves = FALSE)
dev.off()
pdf(file = "../2_Output/Figure_3/Fig3E_VSMC_UMAP_Trajectory_Pseudotime.pdf", height = 3, width = 4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           graph_label_size = 1,
           cell_size = .7,
           label_branch_points = FALSE,
           label_leaves = F)
dev.off()
# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes %>% filter(modulated_genes$q_value < 0.05 & modulated_genes$status =="OK") # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
modulated_genes <- top_n(modulated_genes, 500, -q_value)
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05))
openxlsx::write.xlsx(modulated_genes, "../2_Output/Figure_4/Pseudotime_DEGs.xlsx")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
pt.matrix <- as.data.frame(exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))])
cell_names <- colnames(pt.matrix)
Index<-as.data.frame(cds@colData) %>% dplyr::select(CellTypev2)
Index<-subset(Index, row.names(Index) %in% cell_names)
Index$CellTypev2 <- factor(Index$CellTypev2, levels = c("VSMC_A", "VSMC_B", "VSMC_C", "Fibroblast_C", "Fibroblast_A", "Fibroblast_B", "EC", "Macrophage", "NKT", "Neuronal"))
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=6)$y})) # Create a spline that smooths the pseudotime-based expression along 6 degrees of freedom.
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
colnames(pt.matrix) <- cell_names
###########
paletteLength <- 20
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod1"))(paletteLength)
ann_colors = list(Disease = c(`aortic aneurysm` = "black", 
                              normal="white", 
                              `CABG Aortic Button` = "grey25"),
                  CellTypev2 = c(VSMC_B="lightsalmon", 
                                 VSMC_A = "firebrick4",
                                 VSMC_C = "coral2",
                                 Fibroblast_A="steelblue", 
                                 Fibroblast_C = "lightblue3", 
                                 Neuronal = "darkcyan", 
                                 Fibroblast_B = "steelblue4", 
                                 EC="goldenrod2", 
                                 Macrophage="azure4", 
                                 NKT="seagreen3"))
genes_traj <- modulated_genes %>% top_n(., 25,  -p_value)
# ha = rowAnnotation(foo = anno_mark(at = which(rownames(modulated_genes)==rownames(genes_traj)), labels=rownames(modulated_genes),
#         labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
ha = rowAnnotation(foo = anno_mark(at = which(rownames(modulated_genes)==TCF21_GENES), 
                                   labels=TCF21_GENES,
        labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
heatmap_DMC<-pheatmap::pheatmap(pt.matrix, scale="row", 
                      cluster_cols = F, 
                      cluster_rows = TRUE,
                      cutree_rows = 3,
                      fontsize_col = 8,
                      color = myColor,
                      annotation_col = Index,
                      annotation_colors = ann_colors,
                      show_colnames = F,
                      show_rownames = F,
                      right_annotation = ha,
                      border_color = NA)
pdf("../2_Output/Figure_3/Fig3_Pseudotime_Pathway.Heatmap.pdf", height = 10, width = 7)
heatmap_DMC
dev.off()
################################
hc <-heatmap_DMC$tree_row
lbl <- cutree(hc, 4)
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)
cluster3<-which(lbl==3)
cluster4<-which(lbl==4)
Cluster1_data<-pt.matrix[cluster1,]
Cluster2_data<-pt.matrix[cluster2,]
Cluster3_data<-pt.matrix[cluster3,]
Cluster4_data<-pt.matrix[cluster4,]
Cluster1_GENES <- rownames(Cluster1_data)
Cluster2_GENES <- rownames(Cluster2_data)
Cluster3_GENES <- rownames(Cluster3_data)
Cluster4_GENES <- rownames(Cluster4_data)
write.csv(Cluster1_data, "Cluster1.csv")
write.csv(Cluster2_data, "Cluster2.csv")
write.csv(Cluster3_data, "Cluster3.csv")
write.csv(Cluster4_data, "Cluster4.csv")
##Enrichr
library(enrichR)
dbs <- c("BioPlanet_2019")
enriched_1 <- enrichr(Cluster1_GENES, dbs)
enrich_1<-enriched_1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_1)
enriched_2 <- enrichr(Cluster2_GENES, dbs)
enrich_2<-enriched_2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_2)
enriched_3 <- enrichr(Cluster3_GENES, dbs)
enrich_3<-enriched_3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_3)
enriched_4 <- enrichr(Cluster4_GENES, dbs)
enrich_4<-enriched_4[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_4)
Cluster1_Names<-unlist(strsplit(enrich_1$Genes[1], ";", fixed = T))
Cluster2_Names<-unlist(strsplit(enrich_2$Genes[1], ";", fixed = T))
Cluster3_Names<-unlist(strsplit(enrich_3$Genes[1], ";", fixed = T))
Cluster4_Names<-unlist(strsplit(enrich_4$Genes[1], ";", fixed = T))
############################################################################
GENES_HM<-c(Cluster1_Names,Cluster2_Names,Cluster3_Names, "PLA2G7")
Test<-as.data.frame(pt.matrix)
Test$gene_name<-rownames(pt.matrix)
ha = rowAnnotation(link = anno_mark(at = which(Test$gene_name %in% GENES_HM),
                   labels = as.character(Test[which(Test$gene_name %in% GENES_HM), "gene_name"]),
                   labels_gp = gpar(fontsize = 7),
                   padding = unit(1, "mm"))
                   )
heatmap_combination<-ComplexHeatmap::pheatmap(pt.matrix, scale="row",
                    cluster_cols = F,
                    cluster_rows = T,
                    cutree_rows = 4,
                    # cutree_cols = 3,
                     fontsize_col = 8,
                     color = myColor,
                    annotation_names_col = FALSE,
                    show_colnames = F,
                     show_rownames = F,
                     border_color = NA,
                    annotation_colors = ann_colors,
                    right_annotation = ha, 
                    annotation_col = Index,
                    border = TRUE)
pdf(file = "../2_Output/Figure_3/Trajectory_analysis_Heatmap.pdf", height = 6, width = 5)
heatmap_combination
dev.off()
# Examine specific genes within the module(s) of interest (appears to emrich to metabolic genes)
library("viridis")
library(ggplot2)
library(ggpubr)
contractility_genes <- c("ACTA2", "CNN1", "MYH11", "SMTN", "CARMN", "TGFB1I1", "SRF", "CALD1", "LMOD1", "KCNMB1", "SYNM", "ROCK2", "LAMA2", "")
Transitional_genes <- c("SETD5", "INSR", "OGT", "NOX4", "PLA2G7")
GENES <- Transitional_genes
# GENES <- contractility_genes
# GENES <- modulated_genes %>% top_n(., 8,  -p_value)
# GENES <- rownames(GENES)
p1 <- plot_cells(cds, 
           genes = GENES[1:4],
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
           genes = GENES[5:8],
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
pdf(file="../2_Output/Figure_3/Fig3G_UMAP_Pseudotime_Genes.pdf", height = 3, width = 3)
# ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="right")
p1+p2
dev.off()
# Plot genes according to "pseudotime"
lineage_cds <- cds[rowData(cds)$gene_short_name %in% GENES, ] #colData(cds)$cell_type %in% c("VSMC")
lineage_cds<-order_cells(lineage_cds, root_cells = root_cell)
pdf(file="../2_Output/Figure_3/Pseudotime_curves.pdf", height = 9, width = 3)
monocle3::plot_genes_in_pseudotime(lineage_cds,
                         color_cells_by="ident",
                         min_expr=0)
dev.off()
```

### Trajectory GSEA (VSMC Subtypes)

GO-term enrichment analysis of genes identified via trajectory analysis identified disproportionate representation of ECM-associated genes.


``` r
##Enrichr
# library(enrichR)
library(dplyr)
VSMC_A_DEGs <-openxlsx::read.xlsx("../2_Output/Figure_4/Pseudotime_DEGs.xlsx") %>% filter(p_value < 0.05, morans_I>0.05)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gene_list <- VSMC_A_DEGs$morans_I
names(gene_list) <- VSMC_A_DEGs$gene_name
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list <- names(gene_list)
#
library(genekitr)
# 2nd step: prepare gene set
gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database
# gs <- geneset::getGO(org = "human",ont = "mf")
# 3rd step: GSEA analysis
gse <- genORA(id = gene_list, geneset = gs)[1:5,] %>% mutate(geneID_symbol=geneID)
pdf(file = "../2_Output/Figure_3/Fig3I_Trajectory.Pathways_ORA.pdf", height = 3, width = 6)
plotEnrich(gse, 
           plot_type = "bar", 
           term_metric = "FoldEnrich", 
           stats_metric = "qvalue",
           wrap_length = 30)
dev.off()
Pathway_genes <- gse[1:5, "ID"]
Genes_pathways<-paste(gse[gse[, "ID"] %in% Pathway_genes,"geneID"], collapse = "/")
Top.Pathway_GeneSymbols<-unique(unlist(strsplit(Genes_pathways, "/")))
# Export sheet
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
hdat_vsmc <- subset(TAA.combined, CellTypev2==c("VSMC_B", "VSMC_A", "VSMC_C"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_3/Trajectory_Pathway_heatmap.pdf", height = 3, width = 10)
Clustered_DotPlot(seurat_object = hdat_vsmc, 
                  features = Top.Pathway_GeneSymbols,
                  colors_use_idents = c("#f3766d","#dc8e34","#48b04c"),
                  k = 3,
                  colors_use_exp = myColor,
                  # legend_title_size = 0.1,
                  cluster_ident = T,
                  flip = T)
dev.off()
######## Heatmap of Genes associated with pathways
pdf("../2_Output/Figure_3/Trajectory_gene.matrix.pdf", width = 6, height = 2)
plotEnrich(gse, 
           plot_type = "geneheat",
           wrap_length = 30)
dev.off()
# Export the sheet
expoSheet(data_list = gse,
                    data_name = names(gse),
                    filename = "KEGG_VSMC.xlsx",
                    dir = "../2_Output/Figure_3/")
```

# Figure 4

### VSMC Aneurysm vs. Normal

To identify disease-specific differentially-expressed genes within the VSMC population between aneurysmal and non-diseased control aortic tissues, the following was performed.


``` r
library(ggpubr)
library(Seurat)
library(dplyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = Disease == c("normal", "aortic aneurysm") & CellTypev2==c("VSMC_B"))
Idents(TAA.combined) <- "Disease"
de_VSMC <- FindMarkers(TAA.combined, ident.1 = "aortic aneurysm", ident.2 = "normal", verbose = FALSE, recorrect_umi = FALSE)
openxlsx::write.xlsx(de_VSMC, file = "../2_Output/Figure_4/VSMC_Aneurysm.v.CON_DEGs.xlsx", rowNames=T)
de_VSMC <- openxlsx::read.xlsx("../2_Output/Figure_4/VSMC_Aneurysm.v.CON_DEGs.xlsx", rowNames=T) %>% filter(pct.1>0.05 | pct.2>.05)
# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
results = mutate(de_VSMC, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.0001 & log2FC>1, 
                       "P < 0.0001 and Log(Fold-Change) > 1", 
                       ifelse(p_val<0.0001 & log2FC< 0-1,
                              "P < 0.0001 and Log(Fold-Change) < -1", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -1",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > 1")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=1, size=.5, linetype="dashed") +
  geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 15, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 15, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = "../2_Output/Figure_4/Fig4A_Volcano_Aneurysm.v.Con.pdf", height = 4, width = 5)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 10, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 10, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
## Figure 4B - GSEA
DEGs_UP <- results %>% filter(avg_log2FC > 0) %>% top_n(n = 100, -p_val_adj)
DEGs_DOWN <- results %>% filter(avg_log2FC < 0) %>% top_n(n = 100, -p_val_adj)
##Enrichr
library(enrichR)
dbs <- c("WikiPathway_2023_Human")
enriched_UP <- enrichr(rownames(DEGs_UP), dbs)
enrich_UP<-enriched_UP[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_UP)
enriched_DOWN <- enrichr(rownames(DEGs_DOWN), dbs)
enrich_DOWN<-enriched_DOWN[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_DOWN)
write.csv(enrich_DOWN, "../2_Output/Figure_4/Fig4_Aneurysm_PathwayDEGs_DOWN.csv")
write.csv(enrich_UP, "../2_Output/Figure_4/Fig4_Aneurysm_PathwayDEGs_UP.csv")
DEGS_UP<-unique(unlist(strsplit(enrich_UP$Genes[1:3], split = ";"))) # top 3 pathways
# DEGS_DOWN<-strsplit(enrich_DOWN$Genes[1], split = ";")[[1]]
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = CellTypev2==c("VSMC_B"))
TAA.combined$Disease<-factor(TAA.combined$Disease, levels = c("normal", "aortic aneurysm", "CABG Aortic Button"))
Idents(TAA.combined) <- "Disease"
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_4/Fig4_Aneurysm_Pathway_DEGs.UP.pdf", height = 5, width = 5)
Clustered_DotPlot(seurat_object = TAA.combined, 
                  features = DEGS_UP,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
# pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Pathway_DEGs.DOWN.pdf", height = 5, width = 5)
# Clustered_DotPlot(seurat_object = hdat_vsmc, 
#                   features = DEGS_DOWN,
#                   colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
#                   k = 3,
#                   colors_use_exp = myColor,
#                   cluster_ident = F)
# dev.off()
```

###################################### 

### VSMC CABG vs. Normal

###################################### 


``` r
library(ggpubr)
library(Seurat)
library(dplyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = Disease == c("normal", "CABG Aortic Button") & CellTypev2==c("VSMC_B"))
# TAA.combined <- PrepSCTFindMarkers(TAA.combined, assay = "SCT", verbose = TRUE)
# markers <- FindAllMarkers(
#   object = TAA.combined,
#   assay = "SCT",
#   verbose = T
# )
# hdat_vsmc <- readRDS("../1_Input/TAA_VSMCs_snRNA.rds")
# hdat_vsmc$VSMC_phenotype <- Idents(hdat_vsmc)
Idents(TAA.combined) <- "Disease"
de_VSMC <- FindMarkers(TAA.combined, ident.1 = "CABG Aortic Button", ident.2 = "normal", verbose = FALSE, recorrect_umi = FALSE)
openxlsx::write.xlsx(de_VSMC, file = "../2_Output/Figure_4/VSMC_CABG.v.CON_DEGs.xlsx", rowNames=T)
de_VSMC <- openxlsx::read.xlsx("../2_Output/Figure_4/VSMC_CABG.v.CON_DEGs.xlsx", rowNames=T)
# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
results = mutate(de_VSMC, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.0001 & log2FC>.5, 
                       "P < 0.0001 and Log(Fold-Change) > .5", 
                       ifelse(p_val<0.0001 & log2FC< 0-.5,
                              "P < 0.0001 and Log(Fold-Change) < -.5", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -.5",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > .5")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=1, size=.5, linetype="dashed") +
  geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 15, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 15, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = "../2_Output/Figure_4/Fig4D_Volcano_CABG.v.Con.pdf", height = 4, width = 5)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 10, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 10, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
## Figure 4 - GSEA
DEGs_UP <- results %>% filter(avg_log2FC > 0) %>% top_n(n = 100, -p_val_adj)
DEGs_DOWN <- results %>% filter(avg_log2FC < 0) %>% top_n(n = 100, -p_val_adj)
##Enrichr
library(enrichR)
dbs <- c("WikiPathway_2023_Human")
enriched_UP <- enrichr(rownames(DEGs_UP), dbs)
enrich_UP<-enriched_UP[[dbs]] %>% filter(P.value < 0.05)
head(enrich_UP)
enriched_DOWN <- enrichr(rownames(DEGs_DOWN), dbs)
enrich_DOWN<-enriched_DOWN[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_DOWN)
write.csv(enrich_DOWN, "../2_Output/Figure_4/Fig4_CABG_PathwayDEGs_DOWN.csv")
write.csv(enrich_UP, "../2_Output/Figure_4/Fig4_CABG_PathwayDEGs_UP.csv")
DEGS_UP<-unique(unlist(strsplit(enrich_UP$Genes[1:3], split = ";"))) # top 3 pathways
DEGS_DOWN<-unique(unlist(strsplit(enrich_DOWN$Genes[1:3], split = ";"))) # top 3 pathways
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = CellTypev2==c("VSMC_B"))
TAA.combined$Disease<-factor(TAA.combined$Disease, levels = c("normal", "aortic aneurysm", "CABG Aortic Button"))
Idents(TAA.combined) <- "Disease"
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_4/Fig4_CABG_Pathway_DEGs.UP.pdf", height = 5, width = 5)
Clustered_DotPlot(seurat_object = TAA.combined, 
                  features = DEGS_UP,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
# pdf(file = "../2_Output/Figure_3/Fig3_VSMC_Pathway_DEGs.DOWN.pdf", height = 5, width = 5)
# Clustered_DotPlot(seurat_object = hdat_vsmc, 
#                   features = DEGS_DOWN,
#                   colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
#                   k = 3,
#                   colors_use_exp = myColor,
#                   cluster_ident = F)
# dev.off()
```

### Venn Diagram: Aneurysm vs. CABG


``` r
library(dplyr)
library(pathview)
library(biomaRt)
library(openxlsx)
Aneurysm_DEGs<-read.xlsx("../2_Output/Figure_4/VSMC_Aneurysm.v.CON_DEGs.xlsx", rowNames = T) %>% filter(p_val < 0.05)
Aneurysm_DEGs$GeneName <- rownames(Aneurysm_DEGs)
CABG_DEGs<-read.xlsx("../2_Output/Figure_4/VSMC_CABG.v.CON_DEGs.xlsx", rowNames = T) %>% filter(p_val < 0.05)
CABG_DEGs$GeneName <- rownames(CABG_DEGs)
# Aneurysm Only DEGs
Aneurysm_UP<-dplyr::filter(Aneurysm_DEGs, avg_log2FC>0)
Aneurysm_DOWN<-filter(Aneurysm_DEGs, avg_log2FC<0)
Aneurysm_ONLY<-anti_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName")
Aneurysm_ONLY.UP<-Aneurysm_DEGs %>% filter(avg_log2FC>0)
Aneurysm_ONLY.DOWN<-Aneurysm_DEGs %>% filter(avg_log2FC<0)
write.xlsx(Aneurysm_ONLY, "../2_Output/Aneurysm_ONLY.xlsx", overwrite = TRUE)
# DKO Only DEGs
CABG_UP<-filter(CABG_DEGs, avg_log2FC>0)
CABG_DOWN<-filter(CABG_DEGs, avg_log2FC<0)
CABG_ONLY<-anti_join(CABG_DEGs, Aneurysm_DEGs, by = "GeneName")
CABG_ONLY.UP<-CABG_ONLY %>% filter(avg_log2FC>0)
CABG_ONLY.DOWN<-CABG_ONLY %>% filter(avg_log2FC<0)
write.xlsx(CABG_ONLY, "../2_Output/CABG.ONLY.xlsx", overwrite = TRUE)
# Overlapping DEGs
Conserved_DEGs<-inner_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName") 
rownames(Conserved_DEGs)<-make.unique(Conserved_DEGs$GeneName, sep = ".")
Conserved_DEGs <- Conserved_DEGs %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_CABG", "\\.x"="_Aneurysm")))
Conserved_Both.UP<-Conserved_DEGs %>% filter(avg_log2FC_Aneurysm > 0, avg_log2FC_CABG > 0)
Conserved_Both.DOWN<-Conserved_DEGs %>% filter(avg_log2FC_Aneurysm < 0, avg_log2FC_CABG < 0)
Conserved_Inverse<-Conserved_DEGs %>% filter((avg_log2FC_Aneurysm>0 & avg_log2FC_CABG<0) | (avg_log2FC_Aneurysm<0 & avg_log2FC_CABG>0))
write.xlsx(Conserved_DEGs, "../2_Output/Figure_4/Fig4G_Overlapping_DEGs.xlsx", overwrite = TRUE)
#Merge dataframe for IPA
Merged<-full_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName")
rownames(Merged)<-make.unique(Merged$GeneName, sep = ".")
Merged <- Merged %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_CABG", "\\.x"="_Aneurysm")))
write.xlsx(Merged,"../2_Output/Figure_4/Fig4H_Merged_DEGs.xlsx", overwrite = TRUE)
########### VENN DIAGRAM
library(ggVennDiagram)
x<-list(Aneurysm.vs.CON = Aneurysm_DEGs$GeneName, CABG.vs.CON = CABG_DEGs$GeneName)
pdf("Rplots.pdf", width = 5, height = 5)
ggVennDiagram(x, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white",high = "lightsalmon")
dev.off()
library(VennDiagram)
venn.diagram(x, fill = c("red", "grey"), alpha = c(0.75, 0.75), lty = 'blank', filename = "../2_Output/Figure_4/Fig4G_Venn_VSMC.DEGs.Overlap.svg", imagetype = "svg", na = "remove", disable.logging = T, width = 8, height = 8, units = "in")
#Write excel worksheet
wb_DESeq<-createWorkbook()
#Unfiltered
  addWorksheet(wb_DESeq, "CABG_ONLY_p05")
  writeData(wb_DESeq, "CABG_ONLY_p05", CABG_ONLY, startCol = 1)
#P-value Significant (0.05)
  addWorksheet(wb_DESeq, "Aneurysm_ONLY_p05")
  writeData(wb_DESeq, "Aneurysm_ONLY_p05", Aneurysm_ONLY, startCol = 1)
#Q-value Significant (0.05)
  addWorksheet(wb_DESeq, "Conserved_DEGs")
  writeData(wb_DESeq, "Conserved_DEGs", Conserved_DEGs, startCol = 1)
saveWorkbook(wb_DESeq, file = "../2_Output/Figure_4/Fig4A_Venn.Diagram_VSMCs.xlsx", overwrite = TRUE)
############################################
# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
ALL_DEGs <- full_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName")
rownames(ALL_DEGs) <- ALL_DEGs$GeneName
ALL_DEGs <- ALL_DEGs %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_CABG", "\\.x"="_Aneurysm")))
results <- ALL_DEGs %>% 
  mutate(., sig=ifelse(avg_log2FC_Aneurysm>0 & avg_log2FC_CABG>0, 
                       "avg_log2FC_Aneurysm > 0 avg_log2FC_CABG > 0", 
                       ifelse(avg_log2FC_Aneurysm<0 & avg_log2FC_CABG<0,
                              "avg_log2FC_Aneurysm < 0 and avg_log2FC_CABG < 0", 
                              "Not Sig")
                       )
         ) %>%
  mutate(., avg_pct = rowMeans(dplyr::select(.,pct.1_Aneurysm, pct.2_Aneurysm, pct.1_CABG, pct.2_CABG))) %>%
  mutate(., Abs_FC = abs(avg_log2FC_CABG)+abs(avg_log2FC_Aneurysm))
results$sig<-factor(results$sig, 
levels = c("avg_log2FC_Aneurysm < 0 and avg_log2FC_CABG < 0",
  "Not Sig",
  "avg_log2FC_Aneurysm > 0 avg_log2FC_CABG > 0")
  )
p = ggplot(results, aes(avg_log2FC_Aneurysm, avg_log2FC_CABG)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = avg_pct),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=0, size=.5, linetype="dashed") +
  # geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0, size=.5, linetype="dashed") +
  # geom_hline(yintercept=-1, size=.5, linetype="dashed") +
  labs(x="Aneurysm vs. CON (Log2FC)", y="CABG vs. CON (Log2FC)") + 
  xlim(min(results$avg_log2FC_Aneurysm, na.rm = TRUE),max(results$avg_log2FC_Aneurysm, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(min(results$avg_log2FC_CABG, na.rm = TRUE), max(results$avg_log2FC_CABG, na.rm = TRUE)), expand = c(0,0)) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG< -1 & avg_log2FC_Aneurysm < -1), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG>1 & avg_log2FC_Aneurysm > 1), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = "../2_Output/Figure_4/Fig4G_XY_Combined.pdf", height = 5, width = 5)
  p+
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG< -1 & avg_log2FC_Aneurysm < -1), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG>1 & avg_log2FC_Aneurysm > 1), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
###################
DEGs_UP <- results %>% filter(avg_log2FC_Aneurysm > 1, avg_log2FC_CABG > 1, p_val_adj_Aneurysm < 0.05, p_val_adj_CABG < 0.05)
DEGs_DOWN <- results %>% filter(avg_log2FC_Aneurysm < -1, avg_log2FC_CABG < -1, p_val_adj_Aneurysm < 0.05, p_val_adj_CABG < 0.05)
##Enrichr
library(enrichR)
dbs <- c("WikiPathway_2023_Human")
enriched_UP <- enrichr(rownames(DEGs_UP), dbs)
enrich_UP<-enriched_UP[[dbs]] %>% filter(P.value < 0.05)
head(enrich_UP)
enriched_DOWN <- enrichr(rownames(DEGs_DOWN), dbs)
enrich_DOWN<-enriched_DOWN[[dbs]] %>% filter(P.value < 0.05)
head(enrich_DOWN)
write.csv(enrich_DOWN, "../2_Output/Figure_4/Fig4_MERGE_PathwayDEGs_DOWN.csv")
write.csv(enrich_UP, "../2_Output/Figure_4/Fig4_MERGE_PathwayDEGs_UP.csv")
DEGS_UP<-unique(unlist(strsplit(enrich_UP$Genes[1:3], split = ";"))) # top 3 pathways
DEGS_DOWN<-unique(unlist(strsplit(enrich_DOWN$Genes[1:3], split = ";"))) # top 3 pathways
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = CellTypev2==c("VSMC_B"))
TAA.combined$Disease<-factor(TAA.combined$Disease, levels = c("normal", "aortic aneurysm", "CABG Aortic Button"))
Idents(TAA.combined) <- "Disease"
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_4/Fig4_MERGE_Pathway_DEGs.UP.pdf", height = 5, width = 5)
Clustered_DotPlot(seurat_object = TAA.combined, 
                  features = DEGS_UP,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F)
dev.off()
```

# Figure 5: Fibroblast-Specific Analysis

### Fibroblast Phenotypic Comparison


``` r
library(Seurat)
library(dplyr)
library(ggtrace)
library(ggplot2)
library(ggrepel)
# Add tracings
library(ggplot2)
library(ggtrace)
library(ggthemes)
library(ggrepel)
library(dplyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined[["SCT"]]@SCTModel.list <- TAA.combined[["SCT"]]@SCTModel.list[1]
TAA_normal <- subset(TAA.combined, subset = Disease == "normal") # Use only the control samples for Fibroblast-specific analysis
# TAA_normal <- PrepSCTFindMarkers(TAA_normal, assay = "SCT", verbose = TRUE)
# Differential Expression between Fibroblasts
Fibroblast_A_DEGs <- FindMarkers(TAA_normal, ident.1 = "Fibroblast_A", ident.2 = "Fibroblast_B", assay = "SCT", min.pct = 0.1, logfc.threshold = 0)
Fibroblast_A_DEGs[Fibroblast_A_DEGs$p_val==0,"p_val"] <- 1e-300 # Change P-value of "0" to 1e-300
openxlsx::write.xlsx(Fibroblast_A_DEGs, "../2_Output/Figure_5/DEGs_Fibroblast_A.vs.B.xlsx", rowNames = T)
Fibroblast_A_DEGs <- openxlsx::read.xlsx("../2_Output/Figure_5/DEGs_Fibroblast_A.vs.B.xlsx", rowNames = T) %>% filter(pct.1>0.1, pct.2>0.1)
# Select only the central cluster
hdat_Fibroblast <- subset(TAA_normal, CellTypev2==c("Fibroblast_B", "Fibroblast_A", "Fibroblast_C"))
#Re-cluster
hdat_Fibroblast <- RunUMAP(hdat_Fibroblast, reduction = "harmony", assay = "SCT", dims = 1:40)
hdat_Fibroblast <- FindNeighbors(object = hdat_Fibroblast, reduction = "harmony")
hdat_Fibroblast <- FindClusters(hdat_Fibroblast, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
#
hdat_Fibroblast<-SetIdent(hdat_Fibroblast, value = "SCT_snn_res.0.2")
hdat_Fibroblast <- RenameIdents(hdat_Fibroblast, 
                          `0` = "Fibroblast_A",
                          `1` = "Fibroblast_B",
                          `2` = "Fibroblast_C",
                          `3` = "Fibroblast_C")
UMAP_Fibroblast<-DimPlot(hdat_Fibroblast,  label = T, cols = c("steelblue", "steelblue4", "lightblue3")) + NoLegend()
UMAP_Fibroblast
# Figure 1A - UMAP of Fibroblasts only
pdf("../2_Output/Figure_5/Fig5A_UMAP_Fibroblast.pdf", height = 3.5, width = 3.5)
UMAP_Fibroblast
dev.off()
# Upload Gene Markers
Fibroblast1_genes<-c("ABCA10", "C3", "ADGRD1", "FBLN1", "DCN")
Fibroblast2_genes<-c("NFASC",  "SAMD5", "PRSS23") #"UACA","TEX41",
Fibroblast_C_genes<-c("ADAMTS1", "RGS6", "TNC", "ANGPT2", "DGKG", "GRIP2") # 
# Create Fibroblast Subset UMAP Density plot of Gene Markers
library(Nebulosa)
Fibroblast1_density<-plot_density(hdat_Fibroblast, Fibroblast1_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast2_density<-plot_density(hdat_Fibroblast, Fibroblast2_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
Fibroblast_C_density <- plot_density(hdat_Fibroblast, Fibroblast_C_genes, reduction = "umap", joint = TRUE, combine = FALSE, pal = "magma")
pdf(file = "../2_Output/Figure_5/Fig5_Fibroblast_Density_plots.pdf", height = 3.5, width = 3.5)
Fibroblast1_density
Fibroblast2_density
Fibroblast_C_density
dev.off()
#Volcano Plot
# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
# Read data from the web
options(ggrepel.max.overlaps = Inf)
results = mutate(Fibroblast_A_DEGs, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.0001 & log2FC>1, 
                       "P < 0.0001 and Log(Fold-Change) > 1", 
                       ifelse(p_val<0.0001 & log2FC< 0-1,
                              "P < 0.0001 and Log(Fold-Change) < -1", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -1",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > 1")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  # theme(axis.line = element_blank(),
  #       axis.ticks = element_blank()
  # ) +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=1, size=.5, linetype="dashed") +
  geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 10, -log2FC),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 10, log2FC), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
##
pdf(file = paste0("../2_Output/Figure_5/Fibroblast_B.vs.A_VolcanoPlot.pdf"), height = 4, width = 4)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< -1), 15,  -log2FC), aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 15, log2FC), aes(label=gene_name)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
#########
Fibroblast_type_DEs <- results %>% filter(p_val < 0.05, abs(log2FC)>.5) %>% dplyr::select(gene_name)
Fibroblast_type_DEs <- Fibroblast_type_DEs$gene_name
##
plots_Fibroblast_A <- VlnPlot(hdat_Fibroblast, features = c("SERPINE1", "MALAT1", "PKD1"), group.by = "CellTypev2",
    pt.size = 0,
    cols = c("darkcyan", "lightsalmon", "black", "grey", "dodgerblue3")) & theme(legend.position = 'none', axis.title.x = element_blank())
plots_Fibroblast_B <- VlnPlot(hdat_Fibroblast, features = c("PDE4D", "PDE3A", "LPP"), group.by = "CellTypev2",
    pt.size = 0,
    cols = c("darkcyan", "lightsalmon", "black", "grey", "dodgerblue3")) & theme(legend.position = 'none', axis.title.x = element_blank())

pdf(file = "../2_Output/Figure_5/Fibroblast_A_Top.Violins.pdf", height = 3, width = 6)
plots_Fibroblast_A
dev.off()
pdf(file = "../2_Output/Figure_5/Fibroblast_B_Top.Violins.pdf", height = 3, width = 6)
plots_Fibroblast_B
dev.off()
########## Volcano Genes
library(scCustomize)
pdf(file = "../2_Output/Figure_5/Volcano_Genes.pdf", height = 5, width = 5)
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "lightsalmon"))(paletteLength)
Clustered_DotPlot(seurat_object = hdat_Fibroblast, 
                  features = Fibroblast_type_DEs,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = T)
dev.off()
library(scCustomize)
pdf(file = "../2_Output/Figure_5/DotPlot_Fibroblast_1.v.2.pdf", height = 10, width = 5)
Clustered_DotPlot(seurat_object = hdat_Fibroblast, features = Fibroblast_type_DEs,  x_lab_rotate=F, k = 3)
dev.off()
## Export 
saveRDS(hdat_Fibroblast, file = "../1_Input/TAA_Fibroblasts_snRNA.rds")
```

### Fibroblast GSEA (Fibroblast Subtypes)


``` r
##Enrichr
# library(enrichR)
library(dplyr)
Fibroblast_A_DEGs <-openxlsx::read.xlsx("../2_Output/Figure_5/DEGs_Fibroblast_A.vs.B.xlsx", rowNames = T) %>% filter(p_val < 0.05)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gene_list <- Fibroblast_A_DEGs$avg_log2FC
names(gene_list) <- rownames(Fibroblast_A_DEGs)
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#
library(genekitr)
# 2nd step: prepare gene set
gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database

gs <- geneset::getGO(org = "human",ont = "mf")
# 3rd step: GSEA analysis
gse <- genGSEA(genelist = gene_list, geneset = gs)
plot_ridge <- plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         show_pathway = 25, 
         wrap_length = 1,
         colour = c("navyblue", "orange"))
plot_ridge
pdf(file = "../2_Output/Figure_5/Pathway_Fibroblast.pdf", height = 5, width = 6)
plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         show_pathway = 15, 
         wrap_length = 25,
         colour = c("navyblue", "orange"))
dev.off()

# Pathway_genes <- c("hsa03040", "hsa04261")
# Genes_pathways<-paste(gse$gsea_df[gse$gsea_df[, "ID"] %in% Pathway_genes,"geneID"], collapse = "/")
# Top.Pathway_GeneIDs<-unique(unlist(strsplit(Genes_pathways, "/")))
# Top.Pathway_GeneSymbols <- transId(Top.Pathway_GeneIDs, transTo = "sym")$symbol
# # pdf(file = "../2_Output/Figure_5/GSEA_Top5_Fibroblast.pdf", height = 5, width = 5)
# # plot_classic <- plotGSEA(gse, plot_type = "classic", show_pathway = gse$gsea_df$ID[1:5],label_by = "description")
# # dev.off()
# # Export sheet
# library(Seurat)
# library(scCustomize)
# TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
# hdat_Fibroblast <- subset(TAA.combined, CellTypev2==c("Fibroblast_B", "Fibroblast_A", "Fibroblast_C"))
# paletteLength <- 100
# myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
# ########## Pathway Genes
# pdf(file = "../2_Output/Figure_5/Heatmap_PathwayGenes.Top.pdf", height = 5, width = 5)
# Clustered_DotPlot(seurat_object = hdat_Fibroblast, 
#                   features = Top.Pathway_GeneSymbols,
#                   colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
#                   k = 3,
#                   colors_use_exp = myColor,
#                   cluster_ident = T)
# dev.off()
# expoSheet(data_list = gse,
#                     data_name = names(gse),
#                     filename = "KEGG_Fibroblast.xlsx",
#                     dir = "../2_Output/Figure_5/")
```

### Fibroblast Trajectory Analysis (Monocle 3)


``` r
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(ggtrace)
Merged_scaled <- readRDS(file = "../1_Input/TAA_Fibroblasts_snRNA.rds")
# Merged_scaled <- subset(TAA.combined, subset = Disease == c("normal") & CellTypev2 == c("Fibroblast_A", "Fibroblast_B")) # Select only "normal" samples and Fibroblasts
# Merged_scaled<-TAA.combined
Merged_scaled@active.assay = "RNA"
cds<-SeuratWrappers::as.cell_data_set(Merged_scaled)
# cds <- preprocess_cds(cds, num_dim = 100) # creating errors!
cds <- estimate_size_factors(cds)
# Include gene names (not done by default by the seurat conversion)
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
#  Run Monocle
cds <- cluster_cells(cds) #This step creates "partitions" that are used in the trajectory inference
plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition") # this shows the partitions overlain on the UMAP
cds <- learn_graph(cds, use_partition = TRUE) # creating trajectory inference within each partition
# cds <- order_cells(cds) # Used when setting the nodes; if already known, use the next line
root_cell <- "TCAAGTGGTGCTTATG-1-4_2" # names(which(cds@principal_graph_aux$UMAP$pseudotime==0))
cds<-order_cells(cds, root_cells = root_cell)

# Plot the pseudotime on UMAP
plot_cells(cds, 
           color_cells_by = "pseudotime",
           label_branch_points = FALSE,
           label_leaves = FALSE)
pdf(file = "../2_Output/Figure_5/Fig5_Fibroblast_UMAP_Trajectory_Partition.pdf", height = 3, width = 4)
plot_cells(cds,
           color_cells_by = "partition",
           show_trajectory_graph = F,
           graph_label_size = 1,
           cell_size = .5,
           label_roots = F,
           label_branch_points = FALSE,
           label_leaves = FALSE)
dev.off()
pdf(file = "../2_Output/Figure_5/Fig5E_Fibroblast_UMAP_Trajectory_Pseudotime.pdf", height = 3, width = 4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = F,
           graph_label_size = 1,
           cell_size = .7,
           label_branch_points = FALSE,
           label_leaves = F)
dev.off()

# Identify pseudotime
modulated_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8) # Identify differentially-expressed genes with pseudotime
modulated_genes <- na.omit(modulated_genes) # remove NA's
modulated_genes <- modulated_genes[modulated_genes$q_value < 0.05 & modulated_genes$status =="OK", ] # filter cds results down
modulated_genes <- modulated_genes[order(-modulated_genes$morans_test_statistic), ] # order by moran's test
#### Create a heatmap of genes with similar pseudotime kinetics
genes <- row.names(subset(modulated_genes, q_value < 0.05 & abs(morans_I)>0.1))
openxlsx::write.xlsx(modulated_genes, "../2_Output/Figure_5/Pseudotime_DEGs.xlsx")
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)
pt.matrix <- as.data.frame(exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))])
cell_names <- colnames(pt.matrix)
Index<-as.data.frame(cds@colData) %>% select(CellTypev2)
Index<-subset(Index, row.names(Index) %in% cell_names)
Index$CellTypev2 <- factor(Index$CellTypev2, levels = c("VSMC_A", "VSMC_B", "VSMC_C", "Fibroblast_C", "Fibroblast_A", "Fibroblast_B", "EC", "Macrophage", "NKT", "Neuronal"))
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=6)$y})) # Create a spline that smooths the pseudotime-based expression along 6 degrees of freedom.
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
colnames(pt.matrix) <- cell_names
###########
paletteLength <- 20
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod1"))(paletteLength)
ann_colors = list(Disease = c(`aortic aneurysm` = "black", 
                              normal="white", 
                              `CABG Aortic Button` = "grey25"),
                  CellTypev2 = c(VSMC_B="lightsalmon", 
                                 VSMC_A = "firebrick4",
                                 VSMC_C = "coral2",
                                 Fibroblast_A="steelblue", 
                                 Fibroblast_C = "lightblue3", 
                                 Neuronal = "darkcyan", 
                                 Fibroblast_B = "steelblue4", 
                                 EC="goldenrod2", 
                                 Macrophage="azure4", 
                                 NKT="seagreen3"))
genes_traj <- modulated_genes %>% top_n(., 25,  -p_value)
ha = rowAnnotation(foo = anno_mark(at = which(rownames(modulated_genes)==rownames(genes_traj)), labels=rownames(modulated_genes),
        labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))
heatmap_DMC<-pheatmap::pheatmap(pt.matrix, scale="row", 
                      cluster_cols = F, 
                      cluster_rows = TRUE,
                      cutree_rows = 3,
                      fontsize_col = 8,
                      color = myColor,
                      annotation_col = Index,
                      annotation_colors = ann_colors,
                      show_colnames = F,
                      show_rownames = F,
                      right_annotation = ha,
                      border_color = NA,
                      heatmap_legend_param = list(title = "Expression", 
                                                  at = c(-4,-2,0,2,4)))
pdf("../2_Output/Figure_5/Fig5_Pseudotime_Heatmap.pdf", height = 10, width = 7)
heatmap_DMC
dev.off()
################################
hc <-heatmap_DMC$tree_row
lbl <- cutree(hc, 3)
cluster1<-which(lbl==1)
cluster2<-which(lbl==2)
cluster3<-which(lbl==3)
#
Cluster1_data<-pt.matrix[cluster1,]
Cluster2_data<-pt.matrix[cluster2,]
Cluster3_data<-pt.matrix[cluster3,]
Cluster1_GENES <- rownames(Cluster1_data)
Cluster2_GENES <- rownames(Cluster2_data)
Cluster3_GENES <- rownames(Cluster3_data)
write.csv(Cluster1_data, "Cluster1.csv")
write.csv(Cluster2_data, "Cluster2.csv")
write.csv(Cluster3_data, "Cluster3.csv")
##Enrichr
library(enrichR)
dbs <- c("BioPlanet_2019")
enriched_1 <- enrichr(Cluster1_GENES, dbs)
enrich_1<-enriched_1[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_1)
enriched_2 <- enrichr(Cluster2_GENES, dbs)
enrich_2<-enriched_2[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_2)
enriched_3 <- enrichr(Cluster3_GENES, dbs)
enrich_3<-enriched_3[[dbs]] %>% filter(Adjusted.P.value < 0.05)
head(enrich_3)
Cluster1_Names<-unlist(strsplit(enrich_1$Genes[1], ";", fixed = T))
Cluster2_Names<-unlist(strsplit(enrich_2$Genes[1], ";", fixed = T))
Cluster3_Names<-unlist(strsplit(enrich_3$Genes[1], ";", fixed = T))
############################################################################
GENES_HM<-c(Cluster1_Names, Cluster2_Names,Cluster3_Names)
Test<-as.data.frame(pt.matrix)
Test$gene_name<-rownames(pt.matrix)
ha = rowAnnotation(link = anno_mark(at = which(Test$gene_name %in% GENES_HM),
                   labels = as.character(Test[which(Test$gene_name %in% GENES_HM), "gene_name"]),
                   labels_gp = gpar(fontsize = 7),
                   padding = unit(1, "mm"))
                   )
heatmap_combination<-ComplexHeatmap::pheatmap(pt.matrix, scale="row",
                    cluster_cols = F,
                    cluster_rows = T,
                    cutree_rows = 3,
                    # cutree_cols = 3,
                     fontsize_col = 8,
                     color = myColor,
                    annotation_names_col = FALSE,
                    show_colnames = F,
                     show_rownames = F,
                     border_color = NA,
                    annotation_colors = ann_colors,
                    right_annotation = ha, 
                    annotation_col = Index,
                    border = TRUE)
pdf(file = "../2_Output/Figure_5/Trajectory_analysis_Heatmap.pdf")
heatmap_combination
dev.off()

# Examine specific genes within the module(s) of interest (appears to emrich to metabolic genes)
library("viridis")
library(ggplot2)
library(ggpubr)
GENES <- modulated_genes %>% top_n(., 8,  -p_value)
GENES <- rownames(GENES)
# GENES <- contractility_genes
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
p4 <- plot_cells(cds, 
           genes = GENES[7:8],
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
p5 <- plot_cells(cds, 
           genes = GENES[9:10],
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
pdf(file="../2_Output/Figure_5/Fig5G_UMAP_Pseudotime_Genes.pdf", height = 4, width = 3)
ggarrange(p1, p2,p3, p4,p5, ncol=1, nrow=3, common.legend = TRUE, legend="right")
dev.off()
p1+p2+p3+p4+p5
# Plot genes according to "pseudotime"
lineage_cds <- cds[rowData(cds)$gene_short_name %in% GENES, ] #colData(cds)$cell_type %in% c("Fibroblast")
lineage_cds<-order_cells(lineage_cds, root_cells = root_cell)
pdf(file="../2_Output/Figure_5/Pseudotime_curves.pdf", height = 7, width = 3)
# monocle3::plot_genes_in_pseudotime(lineage_cds,
#                          color_cells_by="ident",
#                          min_expr=0)
# dev.off()
```

### Trajectory GSEA (Fibroblast Subtypes)


``` r
##Enrichr
# library(enrichR)
library(dplyr)
Fibroblast_A_DEGs <-openxlsx::read.xlsx("../2_Output/Figure_5/Pseudotime_DEGs.xlsx") %>% filter(p_value < 0.05, morans_I>0.05)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gene_list <- Fibroblast_A_DEGs$morans_I
names(gene_list) <- Fibroblast_A_DEGs$gene_name
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list <- names(gene_list)
#
library(genekitr)
# 2nd step: prepare gene set
gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database
# gs <- geneset::getGO(org = "human",ont = "mf")
# 3rd step: GSEA analysis
gse <- genORA(id = gene_list, geneset = gs)[1:5,] %>% mutate(geneID_symbol=geneID)
pdf(file = "../2_Output/Figure_5/Fig3I_Trajectory.Pathways_ORA.pdf", height = 3, width = 6)
plotEnrich(gse, 
           plot_type = "bar", 
           term_metric = "FoldEnrich", 
           stats_metric = "qvalue",
           wrap_length = 30)
dev.off()
Pathway_genes <- gse[1:5, "ID"]
Genes_pathways<-paste(gse[gse[, "ID"] %in% Pathway_genes,"geneID"], collapse = "/")
Top.Pathway_GeneSymbols<-unique(unlist(strsplit(Genes_pathways, "/")))
# Export sheet
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
hdat_Fibroblast <- subset(TAA.combined, CellTypev2==c("Fibroblast_B", "Fibroblast_A"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "coral4"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_5/Fibroblast_Trajectory_Pathway.pdf", height = 3, width = 10)
Clustered_DotPlot(seurat_object = hdat_Fibroblast, 
                  features = Top.Pathway_GeneSymbols,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  # legend_title_size = 0.1,
                  cluster_ident = T,
                  flip = T)
dev.off()
######## Heatmap of Genes associated with pathways
pdf("../2_Output/Figure_5/Trajectory_gene.matrix.pdf", width = 6, height = 2)
plotEnrich(gse, 
           plot_type = "geneheat",
           wrap_length = 30)
dev.off()
# Export the sheet
expoSheet(data_list = gse,
                    data_name = names(gse),
                    filename = "KEGG_Fibroblast.xlsx",
                    dir = "../2_Output/Figure_5/")
```

### Fibroblast Aneurysm vs. Normal


``` r
library(ggpubr)
library(Seurat)
library(dplyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = Disease == c("normal", "aortic aneurysm") & CellTypev2==c("Fibroblast_A")) # Use only the control samples for Fibroblast-specific analysis
Idents(TAA.combined) <- "Disease"
de_Fibroblast <- FindMarkers(TAA.combined, ident.1 = "aortic aneurysm", ident.2 = "normal", verbose = FALSE, recorrect_umi = FALSE)
openxlsx::write.xlsx(de_Fibroblast, file = "../2_Output/Figure_5/Fibroblast_Aneurysm.v.CON_DEGs.xlsx", rowNames=T)
de_Fibroblast <- openxlsx::read.xlsx("../2_Output/Figure_5/Fibroblast_Aneurysm.v.CON_DEGs.xlsx", rowNames=T) %>% filter(pct.1>0.05 | pct.2>.05)
# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
results = mutate(de_Fibroblast, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.0001 & log2FC>1, 
                       "P < 0.0001 and Log(Fold-Change) > 1", 
                       ifelse(p_val<0.0001 & log2FC< 0-1,
                              "P < 0.0001 and Log(Fold-Change) < -1", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -1",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > 1")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=1, size=.5, linetype="dashed") +
  geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 15, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 15, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = "../2_Output/Figure_5/Fig5A_Volcano_Aneurysm.v.Con.pdf", height = 5, width = 8)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 10, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 10, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
## Figure 4B - GSEA
library(dplyr)
Fibroblast_aneurysm <-de_Fibroblast %>% filter(p_val < 0.05, pct.1 > 0.05 | pct.2 > 0.05)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gene_list <- Fibroblast_aneurysm$avg_log2FC
names(gene_list) <- rownames(Fibroblast_aneurysm)
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#
library(genekitr)
# 2nd step: prepare gene set
# gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database
# gs <- geneset::getGO(org = "human",ont = "bp")
gs_sig <- geneset::getMsigdb(org = "human",category = "H")
# 3rd step: GSEA analysis
gse <- genGSEA(genelist = gene_list, geneset = gs_sig,
               p_cutoff = 0.05,
               q_cutoff = 0.15,
               padj_method = "bonferroni")
Pathways_aneurysm <- gse$gsea_df
plot_ridge <- plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         show_pathway = 10, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
plot_ridge
pdf(file = "../2_Output/Figure_5/Fig5B_Aneurysm.pdf", height = 5, width = 6)
plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         show_pathway = 5,
         stats_metric = "qvalue",
         wrap_length = 25,
         colour = c("darkcyan", "orange"))
dev.off()
## Figure 4C - Genes
Pathway_genes <- c("HALLMARK_HYPOXIA")
Genes_pathways<-paste(gse$gsea_df[gse$gsea_df[, "ID"] %in% Pathway_genes,"geneID"], collapse = "/")
Top.Pathway_GeneIDs<-unique(unlist(strsplit(Genes_pathways, "/")))
Top.Pathway_GeneSymbols <- transId(Top.Pathway_GeneIDs, transTo = "sym")$symbol
# pdf(file = "../2_Output/Figure_5/GSEA_Top5_Fibroblast.pdf", height = 5, width = 5)
# plot_classic <- plotGSEA(gse, plot_type = "classic", show_pathway = gse$gsea_df$ID[1:5],label_by = "description")
# dev.off()
# Export sheet
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined$Disease<-factor(TAA.combined$Disease, levels = c("normal","aortic aneurysm", "CABG Aortic Button"))
Idents(TAA.combined) <- "Disease"
hdat_Fibroblast <- subset(TAA.combined, CellTypev2==c("Fibroblast_B", "Fibroblast_A"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod2"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_5/Fig5C_Aneurysm_PathwayGenes.pdf", height = 6, width = 5)
Clustered_DotPlot(seurat_object = hdat_Fibroblast, 
                  features = Top.Pathway_GeneSymbols,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  row_label_size = 6,
                  cluster_ident = F,
                  flip = F)
dev.off()
expoSheet(data_list = gse,
                    data_name = names(gse),
                    filename = "KEGG_Fibroblast_Aneurysm.xlsx",
                    dir = "../2_Output/Figure_5/")
```

###################################### 

### Fibroblast CABG vs. Normal

###################################### 


``` r
library(ggpubr)
library(Seurat)
library(dplyr)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- subset(TAA.combined, subset = Disease == c("normal", "CABG Aortic Button") & CellTypev2==c("Fibroblast_A")) # Select only Fibroblasts for this analysis
Idents(TAA.combined) <- "Disease"
de_Fibroblast <- FindMarkers(TAA.combined, ident.1 = "CABG Aortic Button", ident.2 = "normal", verbose = FALSE, recorrect_umi = FALSE)
openxlsx::write.xlsx(de_Fibroblast, file = "../2_Output/Figure_5/Fibroblast_CABG.v.CON_DEGs.xlsx", rowNames=T)
de_Fibroblast <- openxlsx::read.xlsx("../2_Output/Figure_5/Fibroblast_CABG.v.CON_DEGs.xlsx", rowNames=T)
# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
results = mutate(de_Fibroblast, minuslogpvalue = -log(p_val), log2FC=avg_log2FC)
results<-results %>% filter(p_val!=0)
results$gene_name<-rownames(results)
results <- results %>% 
  mutate(., sig=ifelse(p_val<0.0001 & log2FC>1, 
                       "P < 0.0001 and Log(Fold-Change) > 1", 
                       ifelse(p_val<0.0001 & log2FC< 0-1,
                              "P < 0.0001 and Log(Fold-Change) < -1", 
                              "Not Sig")
                       )
         )
results$sig<-factor(results$sig, 
levels = c("P < 0.0001 and Log(Fold-Change) < -1",
  "Not Sig",
  "P < 0.0001 and Log(Fold-Change) > 1")
  )
max(results$minuslogpvalue, na.rm = TRUE)
max(results$log2FC, na.rm = TRUE)
min(results$log2FC, na.rm = TRUE)
p = ggplot(results, aes(log2FC, minuslogpvalue)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = minuslogpvalue),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=1, size=.5, linetype="dashed") +
  geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0-log(0.0001), size=.5, linetype="dashed") +
  labs(x=expression(Log[2](Fold-Change)), y=expression(-Log[10](P-value))) + 
  xlim(min(results$log2FC, na.rm = TRUE),max(results$log2FC, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(0, max(results$minuslogpvalue, na.rm = TRUE)), expand = c(0,0)) +
  # geom_hline(yintercept = 0, size = 1) + 
  # geom_vline(xintercept=0, size=0.5) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 15, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 15, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = "../2_Output/Figure_5/Fig5D_Volcano_CABG.v.Con.pdf", height = 5, width = 8)
  p+
  geom_text_repel(data=top_n(filter(results, log2FC< 0-1), 10, minuslogpvalue),
                  aes(label=gene_name)) +
  geom_text_repel(data=top_n(filter(results, log2FC>1), 10, minuslogpvalue), 
  aes(label=gene_name)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
## Figure 4B - GSEA
library(dplyr)
Fibroblast_CABG <-de_Fibroblast %>% filter(p_val < 0.05)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
gene_list <- Fibroblast_CABG$avg_log2FC
names(gene_list) <- rownames(Fibroblast_CABG)
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#
library(genekitr)
# 2nd step: prepare gene set
# gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database
# gs_enrich <- geneset::getEnrichrdb(org = "human", library = "Reactome_2016")
# gs <- geneset::getGO(org = "human", ont = "mf")
gs_sig <- geneset::getMsigdb(org = "human",category = "H")
# 3rd step: GSEA analysis
# go_ent <- genORA(gene_list, geneset = gs_enrich)
gse <- genGSEA(genelist = gene_list, geneset = gs_sig,
               p_cutoff = 1,
               q_cutoff = 1)
Pathways_CABG <- gse$gsea_df
plot_ridge <- plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         stats_metric = "pvalue",
         show_pathway = 5, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
pdf(file = "../2_Output/Figure_5/Fig5E_CABG.pdf", height = 3, width = 5)
plotGSEA(gse,
         plot_type = "ridge",
         label_by = "description",
         show_pathway = 5, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
dev.off()
## Figure 4C - Genes
Pathway_genes <- c("HALLMARK_UV_RESPONSE_UP")
Genes_pathways<-paste(gse$gsea_df[gse$gsea_df[, "ID"] %in% Pathway_genes,"geneID"], collapse = "/")
Top.Pathway_GeneIDs<-unique(unlist(strsplit(Genes_pathways, "/")))
Top.Pathway_GeneSymbols <- transId(Top.Pathway_GeneIDs, transTo = "sym")$symbol
# pdf(file = "../2_Output/Figure_5/GSEA_Top5_Fibroblast.pdf", height = 5, width = 5)
# plot_classic <- plotGSEA(gse, plot_type = "classic", show_pathway = gse$gsea_df$ID[1:5],label_by = "description")
# dev.off()
# Export sheet
library(Seurat)
library(scCustomize)
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined$Disease<-factor(TAA.combined$Disease, levels = c("normal","aortic aneurysm", "CABG Aortic Button"))
Idents(TAA.combined) <- "Disease"
hdat_Fibroblast <- subset(TAA.combined, CellTypev2==c("Fibroblast_B", "Fibroblast_A"))
paletteLength <- 100
myColor <- colorRampPalette(c("dodgerblue4", "white", "goldenrod2"))(paletteLength)
########## Pathway Genes
pdf(file = "../2_Output/Figure_5/Fig5F_CABG_PathwayGenes.pdf", height = 6, width = 5)
Clustered_DotPlot(seurat_object = hdat_Fibroblast, 
                  features = Top.Pathway_GeneSymbols,
                  colors_use_idents = c("lightsalmon", "darkcyan", "dodgerblue3"),
                  k = 3,
                  colors_use_exp = myColor,
                  cluster_ident = F,
                  flip = F)
dev.off()
expoSheet(data_list = gse,
                    data_name = names(gse),
                    filename = "KEGG_Fibroblast_CABG.xlsx",
                    dir = "../2_Output/Figure_5/")
```

### Venn Diagram: Aneurysm vs. CABG


``` r
library(dplyr)
library(pathview)
library(biomaRt)
library(openxlsx)
Aneurysm_DEGs<-read.xlsx("../2_Output/Figure_5/Fibroblast_Aneurysm.v.CON_DEGs.xlsx", rowNames = T) %>% filter(p_val < 0.05)
Aneurysm_DEGs$GeneName <- rownames(Aneurysm_DEGs)
CABG_DEGs<-read.xlsx("../2_Output/Figure_5/Fibroblast_CABG.v.CON_DEGs.xlsx", rowNames = T) %>% filter(p_val < 0.05)
CABG_DEGs$GeneName <- rownames(CABG_DEGs)
# Aneurysm Only DEGs
Aneurysm_UP<-dplyr::filter(Aneurysm_DEGs, avg_log2FC>0)
Aneurysm_DOWN<-filter(Aneurysm_DEGs, avg_log2FC<0)
Aneurysm_ONLY<-anti_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName")
Aneurysm_ONLY.UP<-Aneurysm_DEGs %>% filter(avg_log2FC>0)
Aneurysm_ONLY.DOWN<-Aneurysm_DEGs %>% filter(avg_log2FC<0)
write.xlsx(Aneurysm_ONLY, "../2_Output/Aneurysm_ONLY.xlsx", overwrite = TRUE)
# DKO Only DEGs
CABG_UP<-filter(CABG_DEGs, avg_log2FC>0)
CABG_DOWN<-filter(CABG_DEGs, avg_log2FC<0)
CABG_ONLY<-anti_join(CABG_DEGs, Aneurysm_DEGs, by = "GeneName")
CABG_ONLY.UP<-CABG_ONLY %>% filter(avg_log2FC>0)
CABG_ONLY.DOWN<-CABG_ONLY %>% filter(avg_log2FC<0)
write.xlsx(CABG_ONLY, "../2_Output/CABG.ONLY.xlsx", overwrite = TRUE)
# Overlapping DEGs
Conserved_DEGs<-inner_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName") 
rownames(Conserved_DEGs)<-make.unique(Conserved_DEGs$GeneName, sep = ".")
Conserved_DEGs <- Conserved_DEGs %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_CABG", "\\.x"="_Aneurysm")))
Conserved_Both.UP<-Conserved_DEGs %>% filter(avg_log2FC_Aneurysm > 0, avg_log2FC_CABG > 0)
Conserved_Both.DOWN<-Conserved_DEGs %>% filter(avg_log2FC_Aneurysm < 0, avg_log2FC_CABG < 0)
Conserved_Inverse<-Conserved_DEGs %>% filter((avg_log2FC_Aneurysm>0 & avg_log2FC_CABG<0) | (avg_log2FC_Aneurysm<0 & avg_log2FC_CABG>0))
write.xlsx(Conserved_DEGs, "../2_Output/Figure_5/Fig5G_Overlapping_DEGs.xlsx", overwrite = TRUE)
#Merge dataframe for IPA
Merged<-full_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName")
rownames(Merged)<-make.unique(Merged$GeneName, sep = ".")
Merged <- Merged %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_CABG", "\\.x"="_Aneurysm")))
write.xlsx(Merged,"../2_Output/Figure_5/Fig5H_Merged_DEGs.xlsx", overwrite = TRUE)
########### VENN DIAGRAM
library(ggVennDiagram)
x<-list(Aneurysm.vs.CON = Aneurysm_DEGs$GeneName, CABG.vs.CON = CABG_DEGs$GeneName)
pdf("Rplots.pdf", width = 5, height = 5)
ggVennDiagram(x, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white",high = "lightsalmon")
dev.off()
library(VennDiagram)
venn.diagram(x, fill = c("red", "grey"), alpha = c(0.75, 0.75), lty = 'blank', filename = "../2_Output/Figure_5/Fig5G_Venn_Fibroblast.DEGs.Overlap.svg", imagetype = "svg", na = "remove", disable.logging = T, width = 8, height = 8, units = "in")
#Write excel worksheet
wb_DESeq<-createWorkbook()
#Unfiltered
  addWorksheet(wb_DESeq, "CABG_ONLY_p05")
  writeData(wb_DESeq, "CABG_ONLY_p05", CABG_ONLY, startCol = 1)
#P-value Significant (0.05)
  addWorksheet(wb_DESeq, "Aneurysm_ONLY_p05")
  writeData(wb_DESeq, "Aneurysm_ONLY_p05", Aneurysm_ONLY, startCol = 1)
#Q-value Significant (0.05)
  addWorksheet(wb_DESeq, "Conserved_DEGs")
  writeData(wb_DESeq, "Conserved_DEGs", Conserved_DEGs, startCol = 1)
saveWorkbook(wb_DESeq, file = "../2_Output/Figure_5/Fig5A_Venn.Diagram_Fibroblasts.xlsx", overwrite = TRUE)
############################################

# Volcano Plot
library(dplyr)
library(ggplot2)
library(ggrepel)
library(openxlsx)
options(ggrepel.max.overlaps = Inf)
ALL_DEGs <- full_join(Aneurysm_DEGs, CABG_DEGs, by = "GeneName")
rownames(ALL_DEGs) <- ALL_DEGs$GeneName
ALL_DEGs <- ALL_DEGs %>% rename_all(~stringr::str_replace_all(.,c("\\.y"="_CABG", "\\.x"="_Aneurysm")))
results <- ALL_DEGs %>% 
  mutate(., sig=ifelse(avg_log2FC_Aneurysm>0 & avg_log2FC_CABG>0, 
                       "avg_log2FC_Aneurysm > 0 avg_log2FC_CABG > 0", 
                       ifelse(avg_log2FC_Aneurysm<0 & avg_log2FC_CABG<0,
                              "avg_log2FC_Aneurysm < 0 and avg_log2FC_CABG < 0", 
                              "Not Sig")
                       )
         ) %>%
  mutate(., avg_pct = rowMeans(dplyr::select(.,pct.1_Aneurysm, pct.2_Aneurysm, pct.1_CABG, pct.2_CABG))) %>%
  mutate(., Abs_FC = abs(avg_log2FC_CABG)+abs(avg_log2FC_Aneurysm))
results$sig<-factor(results$sig, 
levels = c("avg_log2FC_Aneurysm < 0 and avg_log2FC_CABG < 0",
  "Not Sig",
  "avg_log2FC_Aneurysm > 0 avg_log2FC_CABG > 0")
  )
p = ggplot(results, aes(avg_log2FC_Aneurysm, avg_log2FC_CABG)) + 
  theme_classic() +
  geom_point(aes(fill=sig, size = avg_pct),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_vline(xintercept=0, size=.5, linetype="dashed") +
  # geom_vline(xintercept=-1, size=0.5, linetype="dashed") +
  geom_hline(yintercept=0, size=.5, linetype="dashed") +
  # geom_hline(yintercept=-1, size=.5, linetype="dashed") +
  labs(x="Aneurysm vs. CON (Log2FC)", y="CABG vs. CON (Log2FC)") + 
  xlim(min(results$avg_log2FC_Aneurysm, na.rm = TRUE),max(results$avg_log2FC_Aneurysm, na.rm = TRUE)) + 
  scale_y_continuous(limits =c(min(results$avg_log2FC_CABG, na.rm = TRUE), max(results$avg_log2FC_CABG, na.rm = TRUE)), expand = c(0,0)) +
  scale_fill_manual(values=c("darkcyan", "darkgray", "darkgoldenrod1")) +
  scale_size_continuous(range = c(.1, 3))
  p+
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG< -1 & avg_log2FC_Aneurysm < -1), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG>1 & avg_log2FC_Aneurysm > 1), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
pdf(file = "../2_Output/Figure_5/Fig5G_XY_Combined.pdf", height = 5, width = 5)
  p+
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG< -1 & avg_log2FC_Aneurysm < -1), 15, Abs_FC),
                  aes(label=GeneName)) +
  geom_text_repel(data=top_n(filter(results, avg_log2FC_CABG>1 & avg_log2FC_Aneurysm > 1), 15, Abs_FC), 
  aes(label=GeneName)) +
  theme(text = element_text(size=20)) +
  theme(text = element_text(size=14), legend.position="none")
dev.off()
###################
## Figure 4I - GSEA Aneurysm
library(dplyr)
library("org.Hs.eg.db")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
DEGs_filtered <- results %>% filter(p_val_Aneurysm < 0.01) # %>% filter(abs(avg_log2FC_Aneurysm)>.5 & abs(avg_log2FC_CABG)>.5)
gene_list <- DEGs_filtered$avg_log2FC_Aneurysm
names(gene_list) <- DEGs_filtered$GeneName
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#
library(genekitr)
# 2nd step: prepare gene set
# gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database
gs <- geneset::getGO(org = "human",ont = "mf")
# 3rd step: GSEA analysis
gse <- genGSEA(genelist = gene_list, geneset = gs)
Pathways_CABG <- gse$gsea_df
plot_ridge <- plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         stats_metric = "pvalue",
         show_pathway = 25, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
pdf(file = "../2_Output/Figure_5/Fig5I_Merged_Aneurysm_GSEA.pdf", height = 3, width = 7)
plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         show_pathway = 25, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
dev.off()
### CABG
gene_list <- DEGs_filtered$avg_log2FC_CABG
names(gene_list) <- DEGs_filtered$GeneName
gene_list<-na.omit(gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#
library(genekitr)
# 2nd step: prepare gene set
# gs <- geneset::getKEGG(org = "human",category = 'pathway') # This is the pathway database
gs <- geneset::getGO(org = "human",ont = "bp")
# 3rd step: GSEA analysis
gse <- genGSEA(genelist = gene_list, geneset = gs)
Pathways_CABG <- gse$gsea_df
plot_ridge <- plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         stats_metric = "qvalue",
         show_pathway = 25, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
pdf(file = "../2_Output/Figure_5/Fig5I_Merged_CABG_GSEA.pdf", height = 4, width = 6)
plotGSEA(gse, 
         plot_type = "ridge",
         label_by = "description",
         stats_metric = "pvalue",
         show_pathway = 25, 
         wrap_length = 1,
         colour = c("darkcyan", "orange"))
dev.off()
```

# Figure 6: ECLIPSR


``` r
library(Seurat)
library(dplyr)
library(ggtrace)
library(ggplot2)
library(ggrepel)
# Add tracings
library(ggplot2)
library(ggtrace)
library(ggthemes)
library(ggrepel)
library(dplyr)
# PARAMETERS
NUMBER_DEGs <- 250
NUMBER_STAT <- 0.05
# # All cell types
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combied_control <- subset(TAA.combined, subset = Disease == "normal")
TAA.combined_aneurysm <- subset(TAA.combined, subset = Disease == "aortic aneurysm")
TAA.combined_CABG <- subset(TAA.combined, subset = Disease == "CABG Aortic Button")
DEGs_CON<-FindAllMarkers(TAA.combied_control, assay = "RNA", logfc.threshold = 0, min.pct = 0.1)
DEGs_TAA<-FindAllMarkers(TAA.combined_aneurysm, assay = "RNA", logfc.threshold = 0, min.pct = 0.1)
DEGs_CABG<-FindAllMarkers(TAA.combined_CABG, assay = "RNA", logfc.threshold = 0, min.pct = 0.1)
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "DEGs_CON")
  writeData(wb_DESeq, "DEGs_CON", DEGs_CON, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "DEGs_TAA")
  writeData(wb_DESeq, "DEGs_TAA", DEGs_TAA, startCol = 1, rowNames = T)
    addWorksheet(wb_DESeq, "DEGs_CABG")
  writeData(wb_DESeq, "DEGs_CABG", DEGs_CABG, startCol = 1, rowNames = T)
  saveWorkbook(wb_DESeq, file = "../2_Output/Figure_6/ECLIPSR_Disease.Specific.xlsx", overwrite = TRUE)
DEGs_CON <- openxlsx::read.xlsx("../2_Output/Figure_6/ECLIPSR_Disease.Specific.xlsx", sheet = "DEGs_CON")
DEGs_TAA <- openxlsx::read.xlsx("../2_Output/Figure_6/ECLIPSR_Disease.Specific.xlsx", sheet = "DEGs_TAA")  
DEGs_CABG <- openxlsx::read.xlsx("../2_Output/Figure_6/ECLIPSR_Disease.Specific.xlsx", sheet = "DEGs_CABG")
Merged_CON.TAA <- inner_join(DEGs_CON, DEGs_TAA, by = c("gene", "cluster"))

p = ggplot(Merged_CON.TAA, aes(avg_log2FC.x, avg_log2FC.y)) + 
  theme_classic() +
  geom_point(aes(fill=cluster),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_hline(yintercept=0, size=.5, linetype="dashed") +
  geom_vline(xintercept=0, size=0.5, linetype="dashed") +
  labs(x=expression(CON_Log[2](Fold-Change)), y=expression(TAA_Log[2](Fold-Change)))
  p

  
Merged_CON.CABG <- inner_join(DEGs_CON, DEGs_CABG, by = c("gene", "cluster"))

q = ggplot(Merged_CON.CABG, aes(avg_log2FC.x, avg_log2FC.y)) + 
  theme_classic() +
  geom_point(aes(fill=cluster),
             colour="black",
             shape=21,
             stroke = 0,
             alpha = .9) +
  geom_hline(yintercept=0, size=.5, linetype="dashed") +
  geom_vline(xintercept=0, size=0.5, linetype="dashed") +
  labs(x=expression(CON_Log[2](Fold-Change)), y=expression(CABG_Log[2](Fold-Change)))
  q
p + q
# Disease Specific
TAA.combined <- readRDS(file = "../1_Input/TAA_Labelling_snRNA.rds")
TAA.combined <- PrepSCTFindMarkers(TAA.combined, assay = "SCT", verbose = TRUE)
TAA.combined <- SetIdent(TAA.combined, value = "Disease")
# VSMC_A
TAA_VSMCA <- subset(TAA.combined, subset = CellTypev2 == "VSMC_A")
TAA_VSMCA <- PrepSCTFindMarkers(TAA_VSMCA, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_VSMCA.Aneurysm <- FindMarkers(
  object = TAA_VSMCA,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T ) #%>% filter(p_val < 0.05)
markers_VSMCA.Aneurysm <- markers_VSMCA.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_VSMCA.Aneurysm$Disease <- "TAA_v_CON"
markers_VSMCA.Aneurysm$CellType <- "VSMC_A"
## CABG
markers_VSMCA.CABG <- FindMarkers(
  object = TAA_VSMCA,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_VSMCA.CABG <- markers_VSMCA.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_VSMCA.CABG$Disease <- "CABG_v_CON"
markers_VSMCA.CABG$CellType <- "VSMC_A"
# VSMC_B
TAA_VSMCB <- subset(TAA.combined, subset = CellTypev2 == "VSMC_B")
TAA_VSMCB <- Seurat::PrepSCTFindMarkers(TAA_VSMCB, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_VSMCB.Aneurysm <- FindMarkers(
  object = TAA_VSMCB,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "RNA",
  test.use = "MAST",
  verbose = T )
markers_VSMCB.Aneurysm <- markers_VSMCB.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_VSMCB.Aneurysm$Disease <- "TAA_v_CON"
markers_VSMCB.Aneurysm$CellType <- "VSMC_B"
## CABG
markers_VSMCB.CABG <- FindMarkers(
  object = TAA_VSMCB,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "RNA",
  test.use = "MAST",
  verbose = T )
markers_VSMCB.CABG <- markers_VSMCB.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_VSMCB.CABG$Disease <- "CABG_v_CON"
markers_VSMCB.CABG$CellType <- "VSMC_B"
# VSMC_C
TAA_VSMCC <- subset(TAA.combined, subset = CellTypev2 == "VSMC_C")
TAA_VSMCC <- PrepSCTFindMarkers(TAA_VSMCC, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_VSMCC.Aneurysm <- FindMarkers(
  object = TAA_VSMCC,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_VSMCC.Aneurysm <- markers_VSMCC.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_VSMCC.Aneurysm$Disease <- "TAA_v_CON"
markers_VSMCC.Aneurysm$CellType <- "VSMC_C"
## CABG
markers_VSMCC.CABG <- FindMarkers(
  object = TAA_VSMCC,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_VSMCC.CABG <- markers_VSMCC.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_VSMCC.CABG$Disease <- "CABG_v_CON"
markers_VSMCC.CABG$CellType <- "VSMC_C"
# FibroA
TAA_FibroA <- subset(TAA.combined, subset = CellTypev2 == "Fibroblast_A")
TAA_FibroA <- PrepSCTFindMarkers(TAA_FibroA, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_FibroA.Aneurysm <- FindMarkers(
  object = TAA_FibroA,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_FibroA.Aneurysm <- markers_FibroA.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_FibroA.Aneurysm$Disease <- "TAA_v_CON"
markers_FibroA.Aneurysm$CellType <- "Fibroblast_A"
## CABG
markers_FibroA.CABG <- FindMarkers(
  object = TAA_FibroA,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_FibroA.CABG <- markers_FibroA.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_FibroA.CABG$Disease <- "CABG_v_CON"
markers_FibroA.CABG$CellType <- "Fibroblast_A"
# FibroB
TAA_FibroB <- subset(TAA.combined, subset = CellTypev2 == "Fibroblast_B")
TAA_FibroB <- PrepSCTFindMarkers(TAA_FibroB, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_FibroB.Aneurysm <- FindMarkers(
  object = TAA_FibroB,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_FibroB.Aneurysm <- markers_FibroB.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_FibroB.Aneurysm$Disease <- "TAA_v_CON"
markers_FibroB.Aneurysm$CellType <- "Fibroblast_B"
## CABG
markers_FibroB.CABG <- FindMarkers(
  object = TAA_FibroB,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_FibroB.CABG <- markers_FibroB.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_FibroB.CABG$Disease <- "CABG_v_CON"
markers_FibroB.CABG$CellType <- "Fibroblast_B"
# MyoFibro
TAA_MyoFibro <- subset(TAA.combined, subset = CellTypev2 == "Fibroblast_C")
TAA_MyoFibro <- PrepSCTFindMarkers(TAA_MyoFibro, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_MyoFibro.Aneurysm <- FindMarkers(
  object = TAA_MyoFibro,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_MyoFibro.Aneurysm <- markers_MyoFibro.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_MyoFibro.Aneurysm$Disease <- "TAA_v_CON"
markers_MyoFibro.Aneurysm$CellType <- "Fibroblast_C"
## CABG
markers_MyoFibro.CABG <- FindMarkers(
  object = TAA_MyoFibro,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_MyoFibro.CABG <- markers_MyoFibro.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_MyoFibro.CABG$Disease <- "CABG_v_CON"
markers_MyoFibro.CABG$CellType <- "Fibroblast_C"
# EC
TAA_EC <- subset(TAA.combined, subset = CellTypev2 == "EC")
TAA_EC <- PrepSCTFindMarkers(TAA_EC, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_EC.Aneurysm <- FindMarkers(
  object = TAA_EC,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_EC.Aneurysm <- markers_EC.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_EC.Aneurysm$Disease <- "TAA_v_CON"
markers_EC.Aneurysm$CellType <- "EC"
## CABG
markers_EC.CABG <- FindMarkers(
  object = TAA_EC,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_EC.CABG <- markers_EC.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_EC.CABG$Disease <- "CABG_v_CON"
markers_EC.CABG$CellType <- "EC"
# Macrophage
TAA_Macrophage <- subset(TAA.combined, subset = CellTypev2 == "Macrophage")
TAA_Macrophage <- PrepSCTFindMarkers(TAA_Macrophage, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_Macrophage.Aneurysm <- FindMarkers(
  object = TAA_Macrophage,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_Macrophage.Aneurysm <- markers_Macrophage.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_Macrophage.Aneurysm$Disease <- "TAA_v_CON"
markers_Macrophage.Aneurysm$CellType <- "Macrophage"
## CABG
markers_Macrophage.CABG <- FindMarkers(
  object = TAA_Macrophage,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_Macrophage.CABG <- markers_Macrophage.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_Macrophage.CABG$Disease <- "CABG_v_CON"
markers_Macrophage.CABG$CellType <- "Macrophage"
# NKT
TAA_NKT <- subset(TAA.combined, subset = CellTypev2 == "NKT")
TAA_NKT <- PrepSCTFindMarkers(TAA_NKT, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_NKT.Aneurysm <- FindMarkers(
  object = TAA_NKT,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_NKT.Aneurysm <- markers_NKT.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_NKT.Aneurysm$Disease <- "TAA_v_CON"
markers_NKT.Aneurysm$CellType <- "NKT"
## CABG
markers_NKT.CABG <- FindMarkers(
  object = TAA_NKT,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "SCT",
  test.use = "MAST",
  verbose = T )
markers_NKT.CABG <- markers_NKT.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_NKT.CABG$Disease <- "CABG_v_CON"
markers_NKT.CABG$CellType <- "NKT"
# Neuronal
TAA_Neuronal <- subset(TAA.combined, subset = CellTypev2 == "Neuronal")
TAA_Neuronal <- PrepSCTFindMarkers(TAA_Neuronal, assay = "SCT", verbose = TRUE)
## Aneurysm
markers_Neuronal.Aneurysm <- FindMarkers(
  object = TAA_Neuronal,
  ident.1 = "aortic aneurysm",
  ident.2 = "normal",
  assay = "RNA",
  test.use = "MAST",
  verbose = T )
markers_Neuronal.Aneurysm <- markers_Neuronal.Aneurysm %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_Neuronal.Aneurysm$Disease <- "TAA_v_CON"
markers_Neuronal.Aneurysm$CellType <- "Neuronal"
## CABG
markers_Neuronal.CABG <- FindMarkers(
  object = TAA_Neuronal,
  ident.1 = "CABG Aortic Button",
  ident.2 = "normal",
  assay = "RNA",
  test.use = "MAST",
  verbose = T )
markers_Neuronal.CABG <- markers_Neuronal.CABG %>% filter(p_val_adj < NUMBER_STAT) %>% top_n(., NUMBER_DEGs, abs(avg_log2FC))
markers_Neuronal.CABG$Disease <- "CABG_v_CON"
markers_Neuronal.CABG$CellType <- "Neuronal"
# Concatenate Results
markers_ALL <- rbind(markers_VSMCA.Aneurysm, markers_VSMCA.CABG, markers_VSMCB.Aneurysm, markers_VSMCB.CABG, markers_VSMCC.Aneurysm, markers_VSMCC.CABG, markers_FibroA.Aneurysm, markers_FibroA.CABG, markers_FibroB.Aneurysm, markers_FibroB.CABG, markers_MyoFibro.Aneurysm, markers_Macrophage.Aneurysm, markers_Macrophage.CABG, markers_NKT.Aneurysm, markers_NKT.CABG, markers_Neuronal.Aneurysm, markers_Neuronal.CABG)
#### Create File
library(openxlsx)
wb_DESeq<-createWorkbook()
  addWorksheet(wb_DESeq, "markers_ALL")
  writeData(wb_DESeq, "markers_ALL", markers_ALL, startCol = 1, rowNames = T)
  saveWorkbook(wb_DESeq, file = "../2_Output/Figure_6/ECLIPSR_List.xlsx", overwrite = TRUE)
# Pathway Analysis
library(enrichR)
dbs <- c("UK_Biobank_GWAS_v1") #UK_Biobank_GWAS_v1 OMIM_Expanded
#
enriched_VSMCA.Aneurysm <- enrichr(rownames(markers_VSMCA.Aneurysm), dbs)
enrich_VSMCA.Aneurysm<-enriched_VSMCA.Aneurysm[[dbs]] 
#
enriched_VSMCA.CABG <- enrichr(rownames(markers_VSMCA.CABG), dbs)
enrich_VSMCA.CABG<-enriched_VSMCA.CABG[[dbs]] 
#
enriched_VSMCB.Aneurysm <- enrichr(rownames(markers_VSMCB.Aneurysm), dbs)
enrich_VSMCB.Aneurysm<-enriched_VSMCB.Aneurysm[[dbs]] 
#
enriched_VSMCB.CABG <- enrichr(rownames(markers_VSMCB.CABG), dbs)
enrich_VSMCB.CABG<-enriched_VSMCB.CABG[[dbs]] 
#
enriched_VSMCC.Aneurysm <- enrichr(rownames(markers_VSMCC.Aneurysm), dbs)
enrich_VSMCC.Aneurysm<-enriched_VSMCC.Aneurysm[[dbs]] 
#
enriched_VSMCC.CABG <- enrichr(rownames(markers_VSMCC.CABG), dbs)
enrich_VSMCC.CABG<-enriched_VSMCC.CABG[[dbs]] 
#
enriched_FibroA.Aneurysm <- enrichr(rownames(markers_FibroA.Aneurysm), dbs)
enrich_FibroA.Aneurysm<-enriched_FibroA.Aneurysm[[dbs]] 
#
enriched_FibroA.CABG <- enrichr(rownames(markers_FibroA.CABG), dbs)
enrich_FibroA.CABG<-enriched_FibroA.CABG[[dbs]] 
#
enriched_FibroB.Aneurysm <- enrichr(rownames(markers_FibroB.Aneurysm), dbs)
enrich_FibroB.Aneurysm<-enriched_FibroB.Aneurysm[[dbs]] 
#
enriched_FibroB.CABG <- enrichr(rownames(markers_FibroB.CABG), dbs)
enrich_FibroB.CABG<-enriched_FibroB.CABG[[dbs]] 
#
enriched_MyoFibro.Aneurysm <- enrichr(rownames(markers_MyoFibro.Aneurysm), dbs)
enrich_MyoFibro.Aneurysm<-enriched_MyoFibro.Aneurysm[[dbs]] 
#
enriched_MyoFibro.CABG <- enrichr(rownames(markers_MyoFibro.CABG), dbs)
enrich_MyoFibro.CABG<-enriched_MyoFibro.CABG[[dbs]] 
#
enriched_EC.Aneurysm <- enrichr(rownames(markers_EC.Aneurysm), dbs)
enrich_EC.Aneurysm<-enriched_EC.Aneurysm[[dbs]] 
#
enriched_EC.CABG <- enrichr(rownames(markers_EC.CABG), dbs)
enrich_EC.CABG<-enriched_EC.CABG[[dbs]] 
#
enriched_Macrophage.Aneurysm <- enrichr(rownames(markers_Macrophage.Aneurysm), dbs)
enrich_Macrophage.Aneurysm<-enriched_Macrophage.Aneurysm[[dbs]] 
#
enriched_Macrophage.CABG <- enrichr(rownames(markers_Macrophage.CABG), dbs)
enrich_Macrophage.CABG<-enriched_Macrophage.CABG[[dbs]] 
#
enriched_NKT.Aneurysm <- enrichr(rownames(markers_NKT.Aneurysm), dbs)
enrich_NKT.Aneurysm<-enriched_NKT.Aneurysm[[dbs]] 
#
enriched_NKT.CABG <- enrichr(rownames(markers_NKT.CABG), dbs)
enrich_NKT.CABG<-enriched_NKT.CABG[[dbs]] 
#
enriched_Neuronal.Aneurysm <- enrichr(rownames(markers_Neuronal.Aneurysm), dbs)
enrich_Neuronal.Aneurysm<-enriched_Neuronal.Aneurysm[[dbs]]
```


# Supplemental Table: R Session Information

All packages and setting are acquired using the following command:


``` r
sinfo<-devtools::session_info()
sinfo$platform
```

```
##  setting  value
##  version  R version 4.4.1 (2024-06-14)
##  os       macOS 15.0.1
##  system   aarch64, darwin20
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/Los_Angeles
##  date     2024-10-16
##  pandoc   3.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64/ (via rmarkdown)
```

``` r
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
   <th style="text-align:left;">  </th>
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
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> 1.4-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/abind </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> beeswarm </td>
   <td style="text-align:center;"> beeswarm </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/beeswarm </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/beeswarm </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-06-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Biobase </td>
   <td style="text-align:center;"> Biobase </td>
   <td style="text-align:center;"> 2.64.0 </td>
   <td style="text-align:center;"> 2.64.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Biobase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Biobase </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> BiocGenerics </td>
   <td style="text-align:center;"> BiocGenerics </td>
   <td style="text-align:center;"> 0.50.0 </td>
   <td style="text-align:center;"> 0.50.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/BiocGenerics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bslib </td>
   <td style="text-align:center;"> bslib </td>
   <td style="text-align:center;"> 0.8.0 </td>
   <td style="text-align:center;"> 0.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/bslib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cachem </td>
   <td style="text-align:center;"> cachem </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cachem </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> circlize </td>
   <td style="text-align:center;"> circlize </td>
   <td style="text-align:center;"> 0.4.16 </td>
   <td style="text-align:center;"> 0.4.16 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/circlize </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/circlize </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cli </td>
   <td style="text-align:center;"> cli </td>
   <td style="text-align:center;"> 3.6.3 </td>
   <td style="text-align:center;"> 3.6.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cli </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cluster </td>
   <td style="text-align:center;"> cluster </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> 2.1.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cluster </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> codetools </td>
   <td style="text-align:center;"> codetools </td>
   <td style="text-align:center;"> 0.2.20 </td>
   <td style="text-align:center;"> 0.2-20 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/codetools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-31 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> colorspace </td>
   <td style="text-align:center;"> colorspace </td>
   <td style="text-align:center;"> 2.1.1 </td>
   <td style="text-align:center;"> 2.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/colorspace </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cowplot </td>
   <td style="text-align:center;"> cowplot </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/cowplot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> crayon </td>
   <td style="text-align:center;"> crayon </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> 1.5.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/crayon </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/crayon </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> curl </td>
   <td style="text-align:center;"> curl </td>
   <td style="text-align:center;"> 5.2.3 </td>
   <td style="text-align:center;"> 5.2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/curl </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/curl </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:center;"> data.table </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/data.table </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> DelayedArray </td>
   <td style="text-align:center;"> DelayedArray </td>
   <td style="text-align:center;"> 0.30.1 </td>
   <td style="text-align:center;"> 0.30.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/DelayedArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> deldir </td>
   <td style="text-align:center;"> deldir </td>
   <td style="text-align:center;"> 2.0.4 </td>
   <td style="text-align:center;"> 2.0-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/deldir </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> devtools </td>
   <td style="text-align:center;"> devtools </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> 2.4.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/devtools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> digest </td>
   <td style="text-align:center;"> digest </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> 0.6.37 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/digest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dittoSeq </td>
   <td style="text-align:center;"> dittoSeq </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dittoSeq </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dittoSeq </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dotCall64 </td>
   <td style="text-align:center;"> dotCall64 </td>
   <td style="text-align:center;"> 1.1.1 </td>
   <td style="text-align:center;"> 1.1-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dotCall64 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dplyr </td>
   <td style="text-align:center;"> dplyr </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/dplyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ellipsis </td>
   <td style="text-align:center;"> ellipsis </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> 0.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ellipsis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> enrichR </td>
   <td style="text-align:center;"> enrichR </td>
   <td style="text-align:center;"> 3.2 </td>
   <td style="text-align:center;"> 3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/enrichR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/enrichR </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-04-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> evaluate </td>
   <td style="text-align:center;"> evaluate </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/evaluate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fansi </td>
   <td style="text-align:center;"> fansi </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> 1.0.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fansi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fansi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> farver </td>
   <td style="text-align:center;"> farver </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/farver </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastDummies </td>
   <td style="text-align:center;"> fastDummies </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> 1.7.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastDummies </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fastmap </td>
   <td style="text-align:center;"> fastmap </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> 1.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fastmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fitdistrplus </td>
   <td style="text-align:center;"> fitdistrplus </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fitdistrplus </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> forcats </td>
   <td style="text-align:center;"> forcats </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/forcats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/forcats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fs </td>
   <td style="text-align:center;"> fs </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/fs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-25 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future </td>
   <td style="text-align:center;"> future </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> future.apply </td>
   <td style="text-align:center;"> future.apply </td>
   <td style="text-align:center;"> 1.11.2 </td>
   <td style="text-align:center;"> 1.11.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/future.apply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> generics </td>
   <td style="text-align:center;"> generics </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> 0.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/generics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDb </td>
   <td style="text-align:center;"> GenomeInfoDb </td>
   <td style="text-align:center;"> 1.40.1 </td>
   <td style="text-align:center;"> 1.40.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDb </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-16 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> GenomeInfoDbData </td>
   <td style="text-align:center;"> 1.2.12 </td>
   <td style="text-align:center;"> 1.2.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomeInfoDbData </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-02 </td>
   <td style="text-align:center;"> Bioconductor </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GenomicRanges </td>
   <td style="text-align:center;"> GenomicRanges </td>
   <td style="text-align:center;"> 1.56.1 </td>
   <td style="text-align:center;"> 1.56.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GenomicRanges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-16 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggbeeswarm </td>
   <td style="text-align:center;"> ggbeeswarm </td>
   <td style="text-align:center;"> 0.7.2 </td>
   <td style="text-align:center;"> 0.7.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggbeeswarm </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggbeeswarm </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggplot2 </td>
   <td style="text-align:center;"> ggplot2 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> 3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggplot2 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggprism </td>
   <td style="text-align:center;"> ggprism </td>
   <td style="text-align:center;"> 1.0.5 </td>
   <td style="text-align:center;"> 1.0.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggprism </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggprism </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrastr </td>
   <td style="text-align:center;"> ggrastr </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrastr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrastr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggrepel </td>
   <td style="text-align:center;"> ggrepel </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> 0.9.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggrepel </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggridges </td>
   <td style="text-align:center;"> ggridges </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> 0.5.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggridges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggthemes </td>
   <td style="text-align:center;"> ggthemes </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggthemes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggthemes </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ggtrace </td>
   <td style="text-align:center;"> ggtrace </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> 0.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggtrace </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ggtrace </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-06-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GlobalOptions </td>
   <td style="text-align:center;"> GlobalOptions </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> 0.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GlobalOptions </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/GlobalOptions </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-06-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> globals </td>
   <td style="text-align:center;"> globals </td>
   <td style="text-align:center;"> 0.16.3 </td>
   <td style="text-align:center;"> 0.16.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/globals </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> glue </td>
   <td style="text-align:center;"> glue </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> 1.8.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/glue </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> goftest </td>
   <td style="text-align:center;"> goftest </td>
   <td style="text-align:center;"> 1.2.3 </td>
   <td style="text-align:center;"> 1.2-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/goftest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-10-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gridExtra </td>
   <td style="text-align:center;"> gridExtra </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> 2.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gridExtra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2017-09-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> gtable </td>
   <td style="text-align:center;"> gtable </td>
   <td style="text-align:center;"> 0.3.5 </td>
   <td style="text-align:center;"> 0.3.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/gtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> harmony </td>
   <td style="text-align:center;"> harmony </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/harmony </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/harmony </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmltools </td>
   <td style="text-align:center;"> htmltools </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> 0.5.8.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmltools </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> htmlwidgets </td>
   <td style="text-align:center;"> htmlwidgets </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> 1.6.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/htmlwidgets </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httpuv </td>
   <td style="text-align:center;"> httpuv </td>
   <td style="text-align:center;"> 1.6.15 </td>
   <td style="text-align:center;"> 1.6.15 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httpuv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> httr </td>
   <td style="text-align:center;"> httr </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> 1.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/httr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ica </td>
   <td style="text-align:center;"> ica </td>
   <td style="text-align:center;"> 1.0.3 </td>
   <td style="text-align:center;"> 1.0-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ica </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-07-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igraph </td>
   <td style="text-align:center;"> igraph </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/igraph </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IRanges </td>
   <td style="text-align:center;"> IRanges </td>
   <td style="text-align:center;"> 2.38.1 </td>
   <td style="text-align:center;"> 2.38.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/IRanges </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/IRanges </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-03 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> irlba </td>
   <td style="text-align:center;"> irlba </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> 2.3.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/irlba </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-10-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> janitor </td>
   <td style="text-align:center;"> janitor </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> 2.2.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/janitor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/janitor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-02-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jquerylib </td>
   <td style="text-align:center;"> jquerylib </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> 0.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jquerylib </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> jsonlite </td>
   <td style="text-align:center;"> jsonlite </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/jsonlite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> kableExtra </td>
   <td style="text-align:center;"> kableExtra </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/kableExtra </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> KernSmooth </td>
   <td style="text-align:center;"> KernSmooth </td>
   <td style="text-align:center;"> 2.23.24 </td>
   <td style="text-align:center;"> 2.23-24 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/KernSmooth </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> knitr </td>
   <td style="text-align:center;"> knitr </td>
   <td style="text-align:center;"> 1.48 </td>
   <td style="text-align:center;"> 1.48 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/knitr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ks </td>
   <td style="text-align:center;"> ks </td>
   <td style="text-align:center;"> 1.14.3 </td>
   <td style="text-align:center;"> 1.14.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ks </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ks </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> later </td>
   <td style="text-align:center;"> later </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> 1.3.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/later </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lattice </td>
   <td style="text-align:center;"> lattice </td>
   <td style="text-align:center;"> 0.22.6 </td>
   <td style="text-align:center;"> 0.22-6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lattice </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lazyeval </td>
   <td style="text-align:center;"> lazyeval </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lazyeval </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> leiden </td>
   <td style="text-align:center;"> leiden </td>
   <td style="text-align:center;"> 0.4.3.1 </td>
   <td style="text-align:center;"> 0.4.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/leiden </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/leiden </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lifecycle </td>
   <td style="text-align:center;"> lifecycle </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> 1.0.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lifecycle </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> listenv </td>
   <td style="text-align:center;"> listenv </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> 0.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/listenv </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lmtest </td>
   <td style="text-align:center;"> lmtest </td>
   <td style="text-align:center;"> 0.9.40 </td>
   <td style="text-align:center;"> 0.9-40 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lmtest </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lubridate </td>
   <td style="text-align:center;"> lubridate </td>
   <td style="text-align:center;"> 1.9.3 </td>
   <td style="text-align:center;"> 1.9.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lubridate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/lubridate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-09-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> magrittr </td>
   <td style="text-align:center;"> magrittr </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/magrittr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-03-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MASS </td>
   <td style="text-align:center;"> MASS </td>
   <td style="text-align:center;"> 7.3.61 </td>
   <td style="text-align:center;"> 7.3-61 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MASS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Matrix </td>
   <td style="text-align:center;"> Matrix </td>
   <td style="text-align:center;"> 1.7.0 </td>
   <td style="text-align:center;"> 1.7-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Matrix </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MatrixGenerics </td>
   <td style="text-align:center;"> MatrixGenerics </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> 1.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/MatrixGenerics </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> matrixStats </td>
   <td style="text-align:center;"> matrixStats </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/matrixStats </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mclust </td>
   <td style="text-align:center;"> mclust </td>
   <td style="text-align:center;"> 6.1.1 </td>
   <td style="text-align:center;"> 6.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mclust </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mclust </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> memoise </td>
   <td style="text-align:center;"> memoise </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> 2.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/memoise </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mime </td>
   <td style="text-align:center;"> mime </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> 0.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mime </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-09-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> miniUI </td>
   <td style="text-align:center;"> miniUI </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> 0.1.1.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/miniUI </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2018-05-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> munsell </td>
   <td style="text-align:center;"> munsell </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> 0.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/munsell </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mvtnorm </td>
   <td style="text-align:center;"> mvtnorm </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mvtnorm </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/mvtnorm </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Nebulosa </td>
   <td style="text-align:center;"> Nebulosa </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> 1.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Nebulosa </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Nebulosa </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nlme </td>
   <td style="text-align:center;"> nlme </td>
   <td style="text-align:center;"> 3.1.166 </td>
   <td style="text-align:center;"> 3.1-166 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/nlme </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> openxlsx </td>
   <td style="text-align:center;"> openxlsx </td>
   <td style="text-align:center;"> 4.2.7.1 </td>
   <td style="text-align:center;"> 4.2.7.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/openxlsx </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/openxlsx </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> paletteer </td>
   <td style="text-align:center;"> paletteer </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> 1.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/paletteer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/paletteer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> parallelly </td>
   <td style="text-align:center;"> parallelly </td>
   <td style="text-align:center;"> 1.38.0 </td>
   <td style="text-align:center;"> 1.38.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/parallelly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> patchwork </td>
   <td style="text-align:center;"> patchwork </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/patchwork </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pbapply </td>
   <td style="text-align:center;"> pbapply </td>
   <td style="text-align:center;"> 1.7.2 </td>
   <td style="text-align:center;"> 1.7-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pbapply </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pheatmap </td>
   <td style="text-align:center;"> pheatmap </td>
   <td style="text-align:center;"> 1.0.12 </td>
   <td style="text-align:center;"> 1.0.12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pheatmap </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pheatmap </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-01-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pillar </td>
   <td style="text-align:center;"> pillar </td>
   <td style="text-align:center;"> 1.9.0 </td>
   <td style="text-align:center;"> 1.9.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pillar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgbuild </td>
   <td style="text-align:center;"> pkgbuild </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgbuild </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgconfig </td>
   <td style="text-align:center;"> pkgconfig </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> 2.0.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgconfig </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-09-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pkgload </td>
   <td style="text-align:center;"> pkgload </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> 1.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pkgload </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plotly </td>
   <td style="text-align:center;"> plotly </td>
   <td style="text-align:center;"> 4.10.4 </td>
   <td style="text-align:center;"> 4.10.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plotly </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> plyr </td>
   <td style="text-align:center;"> plyr </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> 1.8.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/plyr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> png </td>
   <td style="text-align:center;"> png </td>
   <td style="text-align:center;"> 0.1.8 </td>
   <td style="text-align:center;"> 0.1-8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/png </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/png </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-11-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> polyclip </td>
   <td style="text-align:center;"> polyclip </td>
   <td style="text-align:center;"> 1.10.7 </td>
   <td style="text-align:center;"> 1.10-7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/polyclip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pracma </td>
   <td style="text-align:center;"> pracma </td>
   <td style="text-align:center;"> 2.4.4 </td>
   <td style="text-align:center;"> 2.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pracma </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/pracma </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> profvis </td>
   <td style="text-align:center;"> profvis </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> 0.4.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/profvis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> progressr </td>
   <td style="text-align:center;"> progressr </td>
   <td style="text-align:center;"> 0.14.0 </td>
   <td style="text-align:center;"> 0.14.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/progressr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> promises </td>
   <td style="text-align:center;"> promises </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/promises </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> purrr </td>
   <td style="text-align:center;"> purrr </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> 1.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/purrr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> R6 </td>
   <td style="text-align:center;"> R6 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> 2.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/R6 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-08-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RANN </td>
   <td style="text-align:center;"> RANN </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> 2.6.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RANN </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-25 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RColorBrewer </td>
   <td style="text-align:center;"> RColorBrewer </td>
   <td style="text-align:center;"> 1.1.3 </td>
   <td style="text-align:center;"> 1.1-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RColorBrewer </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2022-04-03 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rcpp </td>
   <td style="text-align:center;"> Rcpp </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> 1.0.13 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rcpp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppAnnoy </td>
   <td style="text-align:center;"> RcppAnnoy </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> 0.0.22 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppAnnoy </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RcppHNSW </td>
   <td style="text-align:center;"> RcppHNSW </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> 0.6.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RcppHNSW </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rematch2 </td>
   <td style="text-align:center;"> rematch2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rematch2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rematch2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> remotes </td>
   <td style="text-align:center;"> remotes </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> 2.5.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/remotes </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reshape2 </td>
   <td style="text-align:center;"> reshape2 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> 1.4.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reshape2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-04-09 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> reticulate </td>
   <td style="text-align:center;"> reticulate </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> 1.39.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/reticulate </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rjson </td>
   <td style="text-align:center;"> rjson </td>
   <td style="text-align:center;"> 0.2.23 </td>
   <td style="text-align:center;"> 0.2.23 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rjson </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rjson </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-16 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rlang </td>
   <td style="text-align:center;"> rlang </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> 1.1.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rlang </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rmarkdown </td>
   <td style="text-align:center;"> rmarkdown </td>
   <td style="text-align:center;"> 2.28 </td>
   <td style="text-align:center;"> 2.28 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rmarkdown </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ROCR </td>
   <td style="text-align:center;"> ROCR </td>
   <td style="text-align:center;"> 1.0.11 </td>
   <td style="text-align:center;"> 1.0-11 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/ROCR </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2020-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RSpectra </td>
   <td style="text-align:center;"> RSpectra </td>
   <td style="text-align:center;"> 0.16.2 </td>
   <td style="text-align:center;"> 0.16-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/RSpectra </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstudioapi </td>
   <td style="text-align:center;"> rstudioapi </td>
   <td style="text-align:center;"> 0.16.0 </td>
   <td style="text-align:center;"> 0.16.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/rstudioapi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Rtsne </td>
   <td style="text-align:center;"> Rtsne </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> 0.17 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Rtsne </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-07 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Arrays </td>
   <td style="text-align:center;"> S4Arrays </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> 1.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Arrays </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Arrays </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S4Vectors </td>
   <td style="text-align:center;"> S4Vectors </td>
   <td style="text-align:center;"> 0.42.1 </td>
   <td style="text-align:center;"> 0.42.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/S4Vectors </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-03 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sass </td>
   <td style="text-align:center;"> sass </td>
   <td style="text-align:center;"> 0.4.9 </td>
   <td style="text-align:center;"> 0.4.9 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scales </td>
   <td style="text-align:center;"> scales </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> 1.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scales </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scattermore </td>
   <td style="text-align:center;"> scattermore </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> 1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scattermore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-06-12 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> scCustomize </td>
   <td style="text-align:center;"> scCustomize </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> 2.1.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scCustomize </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/scCustomize </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-28 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sctransform </td>
   <td style="text-align:center;"> sctransform </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> 0.4.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sctransform </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-19 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sessioninfo </td>
   <td style="text-align:center;"> sessioninfo </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> 1.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sessioninfo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-12-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Seurat </td>
   <td style="text-align:center;"> Seurat </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> 5.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/Seurat </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-10 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SeuratObject </td>
   <td style="text-align:center;"> SeuratObject </td>
   <td style="text-align:center;"> 5.0.2 </td>
   <td style="text-align:center;"> 5.0.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SeuratObject </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shape </td>
   <td style="text-align:center;"> shape </td>
   <td style="text-align:center;"> 1.4.6.1 </td>
   <td style="text-align:center;"> 1.4.6.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shape </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shape </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-02-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shiny </td>
   <td style="text-align:center;"> shiny </td>
   <td style="text-align:center;"> 1.9.1 </td>
   <td style="text-align:center;"> 1.9.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/shiny </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SingleCellExperiment </td>
   <td style="text-align:center;"> SingleCellExperiment </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> 1.26.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SingleCellExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SingleCellExperiment </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> snakecase </td>
   <td style="text-align:center;"> snakecase </td>
   <td style="text-align:center;"> 0.11.1 </td>
   <td style="text-align:center;"> 0.11.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/snakecase </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/snakecase </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-08-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sp </td>
   <td style="text-align:center;"> sp </td>
   <td style="text-align:center;"> 2.1.4 </td>
   <td style="text-align:center;"> 2.1-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/sp </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spam </td>
   <td style="text-align:center;"> spam </td>
   <td style="text-align:center;"> 2.10.0 </td>
   <td style="text-align:center;"> 2.10-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spam </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-23 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SparseArray </td>
   <td style="text-align:center;"> SparseArray </td>
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> 1.4.8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SparseArray </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SparseArray </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.data </td>
   <td style="text-align:center;"> spatstat.data </td>
   <td style="text-align:center;"> 3.1.2 </td>
   <td style="text-align:center;"> 3.1-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.data </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.explore </td>
   <td style="text-align:center;"> spatstat.explore </td>
   <td style="text-align:center;"> 3.3.2 </td>
   <td style="text-align:center;"> 3.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.explore </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.geom </td>
   <td style="text-align:center;"> spatstat.geom </td>
   <td style="text-align:center;"> 3.3.3 </td>
   <td style="text-align:center;"> 3.3-3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.geom </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.random </td>
   <td style="text-align:center;"> spatstat.random </td>
   <td style="text-align:center;"> 3.3.2 </td>
   <td style="text-align:center;"> 3.3-2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.random </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.sparse </td>
   <td style="text-align:center;"> spatstat.sparse </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.sparse </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.univar </td>
   <td style="text-align:center;"> spatstat.univar </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> 3.0-1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.univar </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-09-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.1) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spatstat.utils </td>
   <td style="text-align:center;"> spatstat.utils </td>
   <td style="text-align:center;"> 3.1.0 </td>
   <td style="text-align:center;"> 3.1-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/spatstat.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringi </td>
   <td style="text-align:center;"> stringi </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringi </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-06 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> stringr </td>
   <td style="text-align:center;"> stringr </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> 1.5.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/stringr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-11-14 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SummarizedExperiment </td>
   <td style="text-align:center;"> SummarizedExperiment </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> 1.34.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/SummarizedExperiment </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> survival </td>
   <td style="text-align:center;"> survival </td>
   <td style="text-align:center;"> 3.7.0 </td>
   <td style="text-align:center;"> 3.7-0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/survival </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-06-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> svglite </td>
   <td style="text-align:center;"> svglite </td>
   <td style="text-align:center;"> 2.1.3 </td>
   <td style="text-align:center;"> 2.1.3 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/svglite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-08 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> systemfonts </td>
   <td style="text-align:center;"> systemfonts </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> 1.1.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/systemfonts </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-05-15 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tensor </td>
   <td style="text-align:center;"> tensor </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> 1.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tensor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2012-05-05 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tibble </td>
   <td style="text-align:center;"> tibble </td>
   <td style="text-align:center;"> 3.2.1 </td>
   <td style="text-align:center;"> 3.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tibble </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-03-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyr </td>
   <td style="text-align:center;"> tidyr </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> 1.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyr </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-24 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tidyselect </td>
   <td style="text-align:center;"> tidyselect </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> 1.2.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/tidyselect </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-03-11 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> timechange </td>
   <td style="text-align:center;"> timechange </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> 0.3.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/timechange </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/timechange </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UCSC.utils </td>
   <td style="text-align:center;"> UCSC.utils </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> 1.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/UCSC.utils </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/UCSC.utils </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> urlchecker </td>
   <td style="text-align:center;"> urlchecker </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> 1.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/urlchecker </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2021-11-30 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> usethis </td>
   <td style="text-align:center;"> usethis </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> 3.0.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/usethis </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-29 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> utf8 </td>
   <td style="text-align:center;"> utf8 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> 1.2.4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/utf8 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/utf8 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-10-22 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> uwot </td>
   <td style="text-align:center;"> uwot </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> 0.2.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/uwot </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vctrs </td>
   <td style="text-align:center;"> vctrs </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> 0.6.5 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vctrs </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-01 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> vipor </td>
   <td style="text-align:center;"> vipor </td>
   <td style="text-align:center;"> 0.4.7 </td>
   <td style="text-align:center;"> 0.4.7 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vipor </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/vipor </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-18 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> viridisLite </td>
   <td style="text-align:center;"> viridisLite </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> 0.4.2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/viridisLite </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-05-02 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> withr </td>
   <td style="text-align:center;"> withr </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> 3.0.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/withr </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/withr </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-31 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> WriteXLS </td>
   <td style="text-align:center;"> WriteXLS </td>
   <td style="text-align:center;"> 6.7.0 </td>
   <td style="text-align:center;"> 6.7.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/WriteXLS </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/WriteXLS </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-20 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xfun </td>
   <td style="text-align:center;"> xfun </td>
   <td style="text-align:center;"> 0.47 </td>
   <td style="text-align:center;"> 0.47 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xfun </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-08-17 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xml2 </td>
   <td style="text-align:center;"> xml2 </td>
   <td style="text-align:center;"> 1.3.6 </td>
   <td style="text-align:center;"> 1.3.6 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xml2 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-12-04 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> xtable </td>
   <td style="text-align:center;"> xtable </td>
   <td style="text-align:center;"> 1.8.4 </td>
   <td style="text-align:center;"> 1.8-4 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/xtable </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2019-04-21 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> XVector </td>
   <td style="text-align:center;"> XVector </td>
   <td style="text-align:center;"> 0.44.0 </td>
   <td style="text-align:center;"> 0.44.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/XVector </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/XVector </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> yaml </td>
   <td style="text-align:center;"> yaml </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> 2.3.10 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/yaml </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-07-26 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zip </td>
   <td style="text-align:center;"> zip </td>
   <td style="text-align:center;"> 2.3.1 </td>
   <td style="text-align:center;"> 2.3.1 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zip </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zip </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-01-27 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zlibbioc </td>
   <td style="text-align:center;"> zlibbioc </td>
   <td style="text-align:center;"> 1.50.0 </td>
   <td style="text-align:center;"> 1.50.0 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zlibbioc </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zlibbioc </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2024-04-30 </td>
   <td style="text-align:center;"> Bioconductor 3.19 (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
  <tr>
   <td style="text-align:left;"> zoo </td>
   <td style="text-align:center;"> zoo </td>
   <td style="text-align:center;"> 1.8.12 </td>
   <td style="text-align:center;"> 1.8-12 </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/zoo </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2023-04-13 </td>
   <td style="text-align:center;"> CRAN (R 4.4.0) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library </td>
  </tr>
</tbody>
</table>
