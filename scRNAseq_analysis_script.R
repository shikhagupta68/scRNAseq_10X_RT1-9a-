{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww19400\viewh13300\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # Single-cell RNA sequencing analysis \
\
# library loading\
\
library(Seurat) \
library(Matrix) \
library(mclust) \
library(dplyr) \
library(reshape) \
library(rgl) \
library(pca3d) \
library(ggplot2)\
\
## creating a Seurat object\
\
amlseu_500_3 <- CreateSeuratObject(raw.data = amlseu_500_3.data, min.cells = 53, min.genes = 500, project = "AML1_ETO")\
\
## reading a Seurat object\
\
readRDS("~/Desktop/mm10_2/amlseu_500_3.rds")\
\
amlseu_500_3 <readRDS("~/Desktop/mm10_2/amlseu_500_3_clusterid.rds")\
\
## calculating percentage mitochondrial content\
\
mito.genes <- grep(pattern = "^MT-", x = rownames(x = amlseu_500_3@data), value = TRUE) \
\
percent.mito <- Matrix::colSums(amlseu_500_3@raw.data[mito.genes, ])/Matrix::colSums(amlseu_500_3@raw.data)\
\
## adding gene counts, UMI counts and percentage mitochondrial content as metadata\
\
amlseu_500_3 <- AddMetaData(object = amlseu_500_3, metadata = percent.mito, col.name = "percent.mito")\
\
VlnPlot(object = amlseu_500_3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)\
\
## plotting relationship between UMI counts and gene counts or mitochondrial content\
\
GenePlot(object = amlseu_500_3, gene1 = "nUMI", gene2 = "percent.mito") \
\
GenePlot(object = amlseu_500_3, gene1 = "nUMI", gene2 = "nGene")\
\
## Filtering cells with low and high threshold\
\
amlseu_500_3 <- FilterCells(object = amlseu_500_3, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))\
\
#Detection of variable genes\
\
amlseu_500_3 <- FindVariableGenes(object = amlseu_500_3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)\
\
length(x = amlseu_500_3@var.genes) ## getting the number of variable genes\
\
#Differential expression calculation using DESeq2\
\
amlseu_500_3 <- SetAllIdent(object = amlseu_500_3, id = 'ident')\
\
diff_56_46_deseq <- FindMarkers(object = amlseu_500_3, ident.1 = 'IN56', ident.2 = 'IN46', test.use = "DESeq2")\
\
amlseu_500_3 <- ScaleData(object = amlseu_500_3, vars.to.regress = c("nUMI", "percent.mito")). ## scaling the data\
\
#Linear Dimensionality reduction analysis\
\
amlseu_500_3 <- RunPCA(object = amlseu_500_3, pc.genes = amlseu_500_3@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)\
\
PCAPlot(object = amlseu_500_3, dim.1 = 1, dim.2 = 2) ## plotting PCA\
\
PCHeatmap(object = amlseu_500_3, pc.use = 10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)\
\
amlseu_500_3 <- JackStraw(object = amlseu_500_3, num.replicate = 100, display.progress = FALSE) ## Jackstraw plot for estimating the significant PCs\
\
JackStrawPlot(object = amlseu_500_3, PCs = 1:12) ## visualizing Jackstraw plot\
\
PCElbowPlot(object = amlseu_500_3) ## elbow plot to determine the significant PCs\
\
# Graph based clustering analysis\
\
amlseu_500_3 <- RunPCA(object = amlseu_500_3, pc.genes = amlseu_500_3@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)\
\
amlseu_500_3 <- FindClusters(object = amlseu_500_3, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)\
\
PrintFindClustersParams(object = amlseu_500_3)\
\
# Non-linear dimensionality reduction analysis\
\
amlseu_500_3 <- RunTSNE(object = amlseu_500_3, dims.use = 1:10, do.fast = TRUE) TSNEPlot(object = amlseu_500_3)\
\
# Creation of Cell Data Set object for pseudotime analysis using Monocle v3.0\
\
## reading barcodes from combined gene expression matrix\
\
amlseu_sample_sheet <read.table("~/Desktop/mm10_2/total.barcodes_monocle.txt", sep = "\\t")\
\
## reading gene IDs from combined gene expression matrix\
\
amlseu_gene_annotation_2 <read.table("~/Desktop/mm10_2/total.genes_monocle.txt", sep = "\\t")\
\
rownames(amlseu_sample_sheet) = amlseu_sample_sheet$V2\
\
colnames(amlseu_gene_annotation_2) <- c("V1", "gene_short_name")\
\
rownames(amlseu_gene_annotation_2) = amlseu_gene_annotation_2$gene_short_name\
\
pd_2 <- new("AnnotatedDataFrame", data = amlseu_sample_sheet) fd_2 <- new("AnnotatedDataFrame", data = amlseu_gene_annotation_2)\
\
amlseu_500_3_mono.data <- Read10X(data.dir = "~/Desktop/mm10_2/")\
\
## obtaining raw expression matrix using Seurat object\
\
amlseu_mono <- CreateSeuratObject(raw.data = amlseu_500_3_mono.data, project = "AMLETO_Monocle", min.cells = 53, min.genes = 500) expression_matrix <- amlseu_mono@data\
\
## creating a new cell data set object\
\
aml_m3_cds <- new_cell_data_set(expression_matrix, cell_metadata = amlseu_sample_sheet, gene_metadata = amlseu_gene_annotation_2)\
\
## performing pre-processing and normalization\
\
aml_m3_cds = preprocess_cds(aml_m3_cds, num_dim = 100) plot_pc_variance_explained(aml_m3_cds)\
\
## clustering the cells using k-means clustering\
\
aml_m3_cds <- cluster_cells(aml_m3_cds) plot_cells(aml_m3_cds)\
\
## colouring the cells based on sample ID\
\
colo <c(rep("IN42",515),rep("IN44",497),rep("IN46",374),rep("IN56",353))\
\
colData(aml_m3_cds)$sampleID <- colo\
\
plot_cells(aml_m3_cds, color_cells_by = 'sampleID')\
\
Performing dimensionality reduction analysis\
\
aml_m3_cds <- reduce_dimension(aml_m3_cds, reduction_method = "tSNE")\
\
## clustering the cells using k-means clustering and overlaying on tSNE plot\
\
aml_m3_cds <- cluster_cells(aml_m3_cds, reduction_method = "tSNE")\
\
Pseudotime trajectory analysis\
\
## learning the trajectory\
\
aml_m3_cds <- learn_graph(aml_m3_cds)\
\
plot_cells(aml_m3_cds, color_cells_by = 'sampleID', label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 1.5)\
\
## ordering the cells in pseudotime\
\
aml_m3_cds = order_cells(aml_m3_cds)\
\
plot_cells(aml_m3_cds, color_cells_by = "pseudotime", label_cell_groups = 'sampleID', label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)\
\
## obtaining the principal node\
\
get_earliest_principal_node <- function(aml_m3_cds, id = "IN46")\{ + cell_ids <- which(colData(aml_m3_cds)[,"sampleID"] == id)\
\
+ closest_vertex <aml_m3_cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_ vertex + closest_vertex <- as.matrix(closest_vertex[colnames(aml_m3_cds), ]) + root_pr_nodes <igraph::V(principal_graph(aml_m3_cds)[["UMAP"]])$name[as.numeric(nam es(which.max(table(closest_vertex[cell_ids, ]))))] + root_pr_nodes\}\
\
## ordering the cells relative to principal node\
\
aml_m3_cds = order_cells(aml_m3_cds, root_pr_nodes = get_earliest_principal_node(aml_m3_cds)) plot_cells(aml_m3_cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 1.5)\
\
# Pairwise distance calculation\
\
# Subsetting required matrix of cells\
\
Non_LeuP_amlseu <- SubsetData(object = amlseu_500_3, ident.use = c('B-cells', 'LMPPs', 'CMPs', 'Monocytes'))\
\
# Calculate mean, CV and DM\
\
means_Non_LeuP <- rowMeans(as.matrix(Non_LeuP_amlseu@data)) cv2_Non_LeuP <- apply(as.matrix(Non_LeuP_amlseu@data), 1, var)/means_Non_LeuP^2 dm_Non_LeuP <- DM(means_Non_LeuP, cv2_Non_LeuP)\
\
# Sort based on DM to identify top 500 highly variable genes\
\
dm_Non_LeuP_sorted <- sort(dm_Non_LeuP, decreasing = TRUE) top_500_Non_LeuP <- dm_Non_LeuP_sorted[1:500]\
\
# extracting matrix corresponding to highly variable genes \
\
Non_LeuP_subset_matrix <Non_LeuP_amlseu@data[names(top_500_Non_LeuP),] dim(Non_LeuP_subset_matrix)\
\
# calculating spearman correlation \
\
Non_LeuP_spearman <- cor(x = as.matrix(Non_LeuP_subset_matrix), y = NULL, method = "spearman")\
\
# calculating distance, d \
\
d_Non_LeuP <- sqrt((1-Non_LeuP_spearman)/2) dim(d_Non_LeuP)\
\
# for plotting merging every cell type \'91d\'92 in single vector \
\
d_Non_LeuP_vector <- as.vector(d_Non_LeuP) d_NonLeuP_w_LeuP <- c(d_Non_LeuP_vector, d_LeuP_vector) \
sampid_NonLeuP_w_LeuP <- c(rep("Non-Leukaemic Progenitors", length(d_Non_LeuP_vector)), rep("Leukaemic Progenitors", length(d_LeuP_vector)))\
\
# making a dataframe for \'91d\'92 and \'91sampleID\'92 \
total_NonLeuP_w_LeuP_dataframe <- data.frame(d_NonLeuP_w_LeuP, sampid_NonLeuP_w_LeuP) \
dim(total_NonLeuP_w_LeuP_dataframe)\
\
# assigning column names for defining aesthetics in ggplot \
colnames(total_NonLeuP_w_LeuP_dataframe) <- c("Pairwise_Distance", "Cell_Type") \
p_NonLeuP_w_LeuP_total <- ggplot(total_NonLeuP_w_LeuP_dataframe, aes(factor(Cell_Type), Pairwise_Distance)) \
p_NonLeuP_w_LeuP_total + geom_violin(aes(fill = factor(Cell_Type))) # colours based on cell type\
\
}