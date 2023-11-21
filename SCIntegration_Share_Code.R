# Required Libraries
# -------------------
# Loading the necessary libraries for data manipulation, Seurat objects handling, and plotting
library(SeuratDisk)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(plyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(harmony)

# Data Loading Section
# --------------------
# Convert and load peri-gastruloids data from .h5ad format to Seurat object
Convert('../scRNA-seq data/Our peri-gastruloids data/hu_rna.h5ad', "h5seurat", overwrite = FALSE, assay = "RNA")
gastruloid <- LoadH5Seurat("../scRNA-seq data/Our peri-gastruloids data/hu_rna.h5seurat", misc = FALSE)
save(gastruloid, file = "Gastruloidllz.Rdata")

# Convert and load human reference data (week 3 post-conception) for comparative analysis
Convert('human_pcw3.h5ad', "h5seurat", overwrite = FALSE, assay = "RNA")
humanreference2 <- LoadH5Seurat("human_pcw3.h5seurat", misc = FALSE)
save(humanreference2, file = "humanreference2_pcw3.Rdata")

# Normalization Function (CPM)
# -----------------------------
# Function to normalize data using Counts Per Million (CPM) method
cpm2 <- function (expr_mat) {
  norm_factor <- colSums(expr_mat@data)
  expr_mat@data <- t(t(expr_mat@data) / norm_factor) * 10^6
  return(expr_mat)
}

# Processing Self Data
# --------------------
# Load processed data and add metadata to gastruloid object
load("../data/data1_self/Gastruloidllz.Rdata")
gastruloid@meta.data$datatype <- "Gastruloids"
gastruloid@meta.data$species <- "Human"
gastruloid@assays$RNA <- cpm2(gastruloid@assays$RNA)

# Annotate cell types based on clustering results
# Import cluster-to-cell-type mapping and assign cell types to gastruloid data
gastruloidcluster2celltype <- read.csv("cluster2celltype.csv", header = F)
gastruloidcellcluster <- read.table("cellcluster.txt", header = F)
colnames(gastruloidcellcluster) <- c("cellname", "cluster")
gastruloidcellcluster <- as.data.frame(gastruloidcellcluster)
gastruloidcellcluster$cluster <- as.factor(gastruloidcellcluster$cluster)
gastruloidcellcluster$celltype <- plyr::mapvalues(x = gastruloidcellcluster$cluster, 
                                                  from = gastruloidcluster2celltype$V1, 
                                                  to = gastruloidcluster2celltype$V2)

# Check for consistency in cell names between cluster data and Seurat object
if (identical(gastruloidcellcluster$cellname, colnames(gastruloid@assays$RNA))) {
  gastruloid@meta.data$celltype <- gastruloidcellcluster$celltype
}else{
  print("Cell Names are not same!")
}

# Human Reference Data Processing
# -------------------------------
# Load human reference data and process it for analysis
load("humanreference2_pcw3.Rdata")
humanreference2@meta.data$datatype<-"homo_PCW3"
humanreference2@meta.data$species<-"Human"
humanreference2@assays$RNA<-cpm2(humanreference2@assays$RNA)
hr2_w3meta<-read.csv("w3_meta.csv",header = T)
hrcellname<-colnames(humanreference2@assays$RNA)
hrcellname2<-unique(substring(hrcellname,19))
hr2_w3meta$barcode<-paste(hr2_w3meta$barcode,hrcellname2,sep = "")
index<-colnames(humanreference2@assays$RNA)%in%hr2_w3meta$barcode
humanreference2@meta.data$celltype<-NA
if (identical(colnames(humanreference2@assays$RNA)[index],hr2_w3meta$barcode)) {
  humanreference2@meta.data$celltype[index]<-hr2_w3meta$cell_type
}else{
  print("Cell Names are not same!")
}

# Monkey Ex Vivo Data Processing
# ------------------------------
# Load and process monkey ex vivo data for comparative analysis
monkey_ex<-readRDS("InvitroSeuratProject_74473Cells_V2.rds")
monkey_ex@meta.data$datatype<-"Monkey_exvivo"
monkey_ex@meta.data$species<-"Monkey"
monkey_ex@assays$RNA<-cpm2(monkey_ex@assays$RNA)
#Only use part of samples in this data set
#subset ME23,24,25. Male:Female=1:1
monkey_ex@meta.data$choose<-1:nrow(monkey_ex@meta.data)
male5500 = sample(monkey_ex@meta.data$choose[monkey_ex@meta.data$Day %in% c("ME23","ME24","ME25") & monkey_ex@meta.data$Gender=="Male"],5500)
female5500 = sample(monkey_ex@meta.data$choose[monkey_ex@meta.data$Day %in% c("ME23","ME24","ME25") & monkey_ex@meta.data$Gender=="Female"],5500)
monkeylist<-c(male5500,female5500)
monkey_ex@meta.data$choose[monkeylist]="TURE"
monkey_ex2<-subset(monkey_ex,choose== "TURE") 
# celltype_annotation
monkey_ex2@meta.data$celltype<-monkey_ex2@meta.data$CellType_new

# Monkey Vivo Data Processing
# ---------------------------
# Load and process monkey vivo data, including data cleaning and normalization
monkey_vivo_data<-Read10X(data.dir = "~/llz/data/data3_Macaca_zhai/MFE56636-processed/")
monkey_vivo_meta<-read.csv("~/llz/data/data3_Macaca_zhai/MFE56636-meta.csv",header = T,row.names = "X")   
#subset E20
monkey_vivo_meta_sub = monkey_vivo_meta[monkey_vivo_meta$stage=="E20",]
for (i in 1:nrow(monkey_vivo_meta_sub)) {
  rownames(monkey_vivo_meta_sub)[i]<-strsplit(rownames(monkey_vivo_meta_sub)[i],"_")[[1]][2]
}
monkey_vivo_data_sub<-monkey_vivo_data[,colnames(monkey_vivo_data)%in%rownames(monkey_vivo_meta_sub)]
monkey_vivo2 <- CreateSeuratObject(counts = monkey_vivo_data_sub, min.cells = 3, min.features = 300,
                                   meta.data=monkey_vivo_meta_sub)
monkey_vivo2@meta.data$datatype<-"Monkey_vivo"
monkey_vivo2@meta.data$species<-"Monkey"
monkey_vivo2[["percent.mt"]] <- PercentageFeatureSet(monkey_vivo2, pattern = "^MT-")   
monkey_vivo2 <- subset(x = monkey_vivo2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
monkey_vivo2@assays$RNA<-cpm2(monkey_vivo2@assays$RNA)
# celltype_annotation
monkey_vivo2@meta.data$celltype<-monkey_vivo2@meta.data$cell_type

# Data Integration and Analysis
# -----------------------------
# Merge all data sets (gastruloid, human reference, monkey ex vivo, monkey vivo) for integrated analysis
alldata_celltype <- merge(x = gastruloid, y = c(humanreference2, monkey_ex2, monkey_vivo2))
alldata <- merge(x = gastruloid, y = c(humanreference2, monkey_ex2, monkey_vivo2))

# Split, normalize, and identify variable features in the integrated dataset
alldata.list <- SplitObject(alldata, split.by = "datatype")
for (i in 1:length(x = alldata.list)) {
  alldata.list[[i]] <- NormalizeData(object = alldata.list[[i]], verbose = FALSE)
  alldata.list[[i]] <- FindVariableFeatures(object = alldata.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  print(i)
}

# Find integration anchors and integrate data across different datasets
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:60, k.anchor = 5, k.filter = 20)
alldata <- IntegrateData(anchorset = alldata.anchors, dims = 1:80)

# Scale the integrated data, perform PCA, and generate an elbow plot for dimensionality estimation
alldata <- ScaleData(alldata)
alldata <- RunPCA(object = alldata, npcs = 100)
pdf(file = "alldata_pca_use.pdf")
ElbowPlot(alldata, ndims = 100)
dev.off()

# UMAP Plotting and Clustering
# -----------------------------
# Setting parameters for UMAP analysis
pcause = 100  # Number of dimensions to consider before analysis
alldata <- FindNeighbors(object = alldata, dims = 1:pcause)
alldata <- FindClusters(object = alldata, resolution = 0.4)
alldata <- RunUMAP(object = alldata, dims = 1:pcause)

# Generating UMAP plots and saving them as PDF files
pdf(file = paste("alldata_tsne_cluster2", pcause, ".pdf", sep = ""), width = 8, height = 6)
DimPlot(object = alldata, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

pdf(file = "alldata_tsne_cluster2.pdf", width = 10, height = 8)
DimPlot(object = alldata, reduction = "umap", group.by = "datatype", pt.size = 0.1)
dev.off()

# Harmony Integration
# -------------------
# Running Harmony algorithm for data integration across different conditions
pdf(file = "alldata_harmony_use.pdf", width = 8, height = 6)
alldata_harmony <- RunHarmony(alldata, group.by.vars = "datatype", plot_convergence = TRUE)
dev.off()

# Plotting Harmony integrated data and saving the output
pdf(file = "alldata_harmony_cluster.pdf", width = 8, height = 6)
DimPlot(object = alldata_harmony, reduction = "harmony", group.by = "datatype", pt.size = 0.3)
dev.off()

# Further UMAP analysis on Harmony integrated data
alldata_harmony <- alldata_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:60) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  identity()

# Saving UMAP plots of Harmony integrated data
pdf(file = "alldata_harmony_cluster.pdf", width = 8, height = 6)
DimPlot(object = alldata_harmony, reduction = "umap", label = TRUE, pt.size = 0.05)
dev.off()

pdf(file = "alldata_harmony_cluster2.pdf", width = 8, height = 6)
DimPlot(object = alldata_harmony, reduction = "umap", group.by = "datatype", pt.size = 0.05)
dev.off()

# Adjusting metadata for Harmony data
alldata_harmony@meta.data$datatype <- factor(alldata_harmony@meta.data$datatype, levels = c("homo_PCW3", "Monkey_vivo", "Monkey_exvivo", "Gastruloids"))
toplist <- names(alldata_harmony$datatype[alldata_harmony$datatype == "Gastruloids"])

# Customized UMAP plots with specific color schemes
# Changing colors of the clusters and generating multiple UMAP plots
pdf(file = "alldata_harmony_cluster_changecolor_threecolor1.pdf", width = 8, height = 6)
DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#2979A9", "#A9A9A9", "#A9A9A9", "#C04844"), order = toplist, group.by = "datatype", pt.size = 0.01)
dev.off()

pdf(file = "alldata_harmony_cluster_changecolor_threecolor2.pdf", width = 8, height = 6)
DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#A9A9A9", "#2979A9", "#A9A9A9", "#C04844"), order = toplist, group.by = "datatype", pt.size = 0.01)
dev.off()

pdf(file = "alldata_harmony_cluster_changecolor_threecolor3.pdf", width = 8, height = 6)
DimPlot(object = alldata_harmony, reduction = "umap", cols = c("#A9A9A9", "#A9A9A9", "#2979A9", "#C04844"), order = toplist, group.by = "datatype", pt.size = 0.01)
dev.off()


# Celltype annotation for all single cells 
# ----------------------------------------
# loading all celltype annotation files
alldata_harmony@meta.data$celltype=NA
gd_celltype<-cbind(colnames(gastruloid@assays$RNA),as.matrix(gastruloid@meta.data$celltype))
gd_celltype[,1]=paste(gd_celltype[,1],"_1",sep = "")
for (i in 1:length(colnames(gastruloid@assays$RNA))) {
  if(sum(colnames(alldata_harmony@assays$RNA)==gd_celltype[i,1])){
    print(i)
    alldata_harmony@meta.data$celltype[colnames(alldata_harmony@assays$RNA)==gd_celltype[i,1]]=gd_celltype[i,2]
  }
}
hr2_celltype<-cbind(colnames(humanreference2@assays$RNA),as.matrix(humanreference2@meta.data$celltype))
hr2_celltype[,1]=paste(hr2_celltype[,1],"_2",sep = "")
for (i in 1:length(colnames(humanreference2@assays$RNA))) {
  if(sum(colnames(alldata_harmony@assays$RNA)==hr2_celltype[i,1])){
    print(i)
    alldata_harmony@meta.data$celltype[colnames(alldata_harmony@assays$RNA)==hr2_celltype[i,1]]=hr2_celltype[i,2]
  }
}
me_celltype<-cbind(colnames(monkey_ex2@assays$RNA),as.matrix(monkey_ex2@meta.data$celltype))
me_celltype[,1]=paste(me_celltype[,1],"_3",sep = "")
for (i in 1:length(colnames(monkey_ex2@assays$RNA))) {
  if(sum(colnames(alldata_harmony@assays$RNA)==me_celltype[i,1])){
    print(i)
    alldata_harmony@meta.data$celltype[colnames(alldata_harmony@assays$RNA)==me_celltype[i,1]]=me_celltype[i,2]
  }
}
mv_celltype<-cbind(colnames(monkey_vivo2@assays$RNA),as.matrix(monkey_vivo2@meta.data$celltype))
mv_celltype[,1]=paste(mv_celltype[,1],"_4",sep = "")
for (i in 1:length(colnames(monkey_vivo2@assays$RNA))) {
  if(sum(colnames(alldata_harmony@assays$RNA)==mv_celltype[i,1])){
    print(i)
    alldata_harmony@meta.data$celltype[colnames(alldata_harmony@assays$RNA)==mv_celltype[i,1]]=mv_celltype[i,2]
  }
}
write.csv(t(table(alldata_harmony$celltype)), file ="celltype_name_integration.txt")
alldata_harmony@meta.data$celltype<-as.factor(alldata_harmony@meta.data$celltype)

# Discuss Cell Type with Other Researchers
# -----------------------------------------
# Read in the new cell type classification data
celltypenew <- read.csv(file = "/export/home/lileijie/llz/data/data1_self/celltype_tonew.csv", header = TRUE)

# Assigning the active identities (clusters) and new cell type annotations to the Harmony object
alldata_harmony@meta.data$cluster <- alldata_harmony@active.ident
alldata_harmony@meta.data$celltypenew <- plyr::mapvalues(x = alldata_harmony@meta.data$celltype,
                                                         from = celltypenew$CellType,
                                                         to = celltypenew$Integration_CellType)

# Statistical Analysis of Cell Types
# ----------------------------------
# Computing the proportion of each cell type within each cluster
adct_statistic <- cbind(as.matrix(alldata_harmony@meta.data$celltypenew), as.matrix(alldata_harmony@meta.data$cluster))
for (i in levels(alldata_harmony@active.ident)) {
  i_statistic <- table(adct_statistic[adct_statistic[, 2] == i, 1]) / length(adct_statistic[adct_statistic[, 2] == i, 1])
  if (i == "0") {
    df.plot <- cbind(names(i_statistic), i_statistic, i)
  } else {
    df.plot <- rbind(df.plot, cbind(names(i_statistic), i_statistic, i))
  }
}
colnames(df.plot) <- c("CellType", "pct.ct", "cluster")
df.plot <- as.data.frame(df.plot)
df.plot$cluster <- factor(df.plot$cluster, levels = c("0", "18", "5", "6", "10", "8", "11", "16", "2", "13", "17", "1", "12", "3", "4", "7", "9", "15", "14"))

# Dot Plot Generation
# --------------------
# Creating a dot plot to visualize the proportion of each cell type within each cluster
p <- ggplot(df.plot, aes(x = CellType, y = as.numeric(cluster), size = as.numeric(pct.ct)/100), color = as.numeric(pct.ct)) +
  geom_point() + scale_size("% detected", range = c(0, 6)) +
  scale_color_gradientn(colours = c("white", "red")) + cowplot::theme_cowplot() +
  ylab("") + xlab("Markers") + theme_bw() + 
  scale_y_continuous(breaks = 1:length(levels(df.plot$cluster)), labels = levels(df.plot$cluster), sec.axis = dup_axis()) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Saving the dot plot as a PDF
pdf(file = "gastruloid2_celltype_dotplot3.pdf", width = 12, height = 6)
print(p)
dev.off()

setwd("../result")
mc2c <- read.csv(file = "merge_cluster2celltype.csv", header = TRUE)
alldata_harmony@meta.data$mergecelltype <- plyr::mapvalues(x = alldata_harmony@meta.data$cluster, from = mc2c$cluster, to = mc2c$cellname)

# UMAP Plot for Merged Cell Types
# -------------------------------
# Generating a UMAP plot based on merged cell types
pdf(file = "alldata_mergecelltype.pdf", width = 12, height = 8)
DimPlot(object = alldata_harmony, reduction = "umap", group.by = "mergecelltype", label = TRUE, pt.size = 0.3)
dev.off()

# Bar Plot for Cell Type Count across Data Sources
# ------------------------------------------------
# Preparing data for a bar plot by counting cell types across different data sources
merge_bar_all <- cbind(as.matrix(alldata_harmony$mergecelltype), as.matrix(alldata_harmony$datatype))
for (i in 1:length(unique(merge_bar_all[, 1]))) {
  sign <- table(merge_bar_all[merge_bar_all[, 1] == unique(merge_bar_all[, 1])[i], 2])
  singletype_table <- cbind(rep(unique(merge_bar_all[, 1])[i], length(sign)), names(sign), as.matrix(sign))
  rownames(singletype_table) <- NULL
  if (i == 1) {
    table_type <- singletype_table
  } else {
    table_type <- rbind(table_type, singletype_table)
  }
}
colnames(table_type) <- c("Celltype", "Datasource", "Cellnumber")

# Creating a dataframe and formatting it for plotting
table_sample_type <- as.data.frame(table_type)
table_sample_type$Datasource <- as.factor(table_sample_type$Datasource)
table_sample_type$Cellnumber <- as.numeric(as.matrix(table_sample_type$Cellnumber))

# Generating a bar plot showing cell type count across data sources
p_bar <- ggplot(table_sample_type, aes(x = Celltype, y = Cellnumber, fill = Datasource)) + geom_bar(stat = "identity") + theme_bw()
ggsave("allmerge_celltype_count.pdf", p_bar, width = 10, height = 12)

# Bar Plot for Cell Type Proportions
# ----------------------------------
# Computing the percentage weight of each cell type within each data source
ce <- ddply(table_sample_type, "Celltype", transform, percent_weight = Cellnumber / sum(Cellnumber) * 100)

# Setting the order of cell types based on certain criteria
celltype_level = c("Progenitors (Epi, PS)", ce$Celltype[ce$Datasource == "Gastruloids"][order(ce$percent_weight[ce$Datasource == "Gastruloids"], decreasing = FALSE)])
ce$Celltype <- factor(ce$Celltype, levels = celltype_level)
ce$Datasource <- factor(ce$Datasource, levels = c("Monkey_exvivo", "Monkey_vivo", "homo_PCW3", "Gastruloids"))

# Adding text labels to the bar plot
p_bar3 <- p_bar2 + geom_text(aes(label = round(percent_weight, 2), x = percent_weight), size = 5, color = "black")
ggsave("allmerge_celltype_percent.pdf", p_bar3, width = 10, height = 8)

# Finding DEGs for Cell Types
# ---------------------------
# Setting active identifiers to merged cell types for DEG analysis
alldata_harmony@active.ident <- alldata_harmony$mergecelltype

# Repeating the process of identifying and filtering DEGs for cell types
alldata.markers <- FindAllMarkers(object = alldata_harmony, only.pos = FALSE)
alldata.markers <- alldata.markers[alldata.markers$p_val_adj < 0.05,]
alldata.markers <- alldata.markers[abs(alldata.markers$avg_log2FC) > 1,]
alldata.markers %>% group_by(cluster) %>% top_n(25, avg_log2FC) -> top25_up
alldata.markers %>% group_by(cluster) %>% top_n(25, -avg_log2FC) -> top25_down
top25_up %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
top25_up <- top25_up[top25_up$avg_log2FC > 0,]
top25_down <- top25_down[top25_down$avg_log2FC < 0,]

# Writing the top DEGs for cell types to text files
write.table(x = top25_up, file = "alldata_celltype_top25_UP_DEGs.txt", sep = "\t", row.names = FALSE)
write.table(x = top25_down, file = "alldata_celltype_top25_Down_DEGs.txt", sep = "\t", row.names = FALSE)

# Creating a heatmap for top DEGs in cell types
pdf(file = "alldata_celltype_heatmap.pdf", width = 14, height = 12)
DoHeatmap(object = alldata_harmony, features = top10$gene, label = TRUE, group.bar = TRUE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()