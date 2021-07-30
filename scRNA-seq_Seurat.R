##The codes were written following the tutorial provided by Seurat (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)

#Load all the necessary packages
library(dplyr)
library(patchwork)
library(cowplot)
library(Seurat)
library(tidyverse)
#Read Data
sc.data <- Read10X(data.dir = "E:/Single cell Files/filtered_feature_bc_matrix_scGIST-T1-4H-Imatinib")
#Create Seurat Object
sc<- CreateSeuratObject(counts = sc.data, min.cells = 3,min.features = 200, project = "sc4hT", assay = "RNA")

#QC
#Identify the percentage of mitochondrial genes

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

#Create plots with all the features
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)


#Based on the violin plots, set the limits for filtering
sc <- subset(sc, subset = percent.mt < 18 & nFeature_RNA > 1000)

#Normalization
sc <- NormalizeData(sc)

#Finding Variable Features
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

##Perform the following codes if you do not want to regress out the effect of cell cycle genes

#From here
all.genes <- row.names(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(object = sc, npcs = 30, verbose = FALSE)
DimPlot(sc)

#Till here

##If you want to regress out the effect of cell cycle genes, perform the following codes

#Regressing out cc genes
#Getting the cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Assign Cell Cycle Scores
sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#Scale data
sc <- ScaleData(sc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(sc))
#PCA with cc genes as features
sc <- RunPCA(sc, features = c(s.genes, g2m.genes))

#PCA with variable features
sc <- RunPCA(sc, features = VariableFeatures(sc), nfeatures.print = 10)
DimPlot(sc)

#Visualization of PCA
DimPlot(sc, reduction = "pca")

#Identification of significant Principle Components (PCs)

sc <- JackStraw(object = sc, reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)
sc <- ScoreJackStraw(object = sc, dims = 1:20, reduction = "pca")

#Visualization with p-value to check the significant PCs
JackStrawPlot(object = sc, dims = 1:20, reduction = "pca")

#In my case, all the 20 PCs have very low p values, therefore, all the 20 PCs are statistically signifcant

#Finding clusters
sc <- FindNeighbors(sc, reduction = "pca", dims = 1:20)
sc <- FindClusters(sc, resolution = 0.5, algorithm = 1)

##tSNE
sc <- RunTSNE(object = sc, dims.use = 1:20, do.fast = TRUE)
#Visualization of the tSNE plot
DimPlot(sc, reduction = "tsne")

##umap
sc <- RunUMAP(sc, reduction = "pca", dims=1:20)
#Visualization of the UMAP
DimPlot(sc, reduction = "umap")

#Visualize the expression of specific genes in the data. For my data, KIT and MYC were significant
FeaturePlot(sc, features = c("KIT","MYC"), cols = c("lightgrey", "red"))

##Identification and visualization of top 10 markers
sc_markers <- FindAllMarkers(sc, logfc.threshold = 0.25, 
                             only.pos = TRUE)
top10 <- sc_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(sc, features = top10$gene) +NoLegend()

#adding +NoLegend() will remove the legends from the heatmap, if it is relavant for the work

#Geneset Enrichment Analysis
library(presto)
library(msigdbr)
library(fgsea)
library(spatstat)

x_lab_name='Gene expression clusters'

#1. Output of all the expressed genes
sc.genes <- wilcoxauc(sc, group_by = "seurat_clusters", seurat_assay='RNA')
sc.genes

# 2. Hallmark geneset will be used for the gsea analysis

m_df<- msigdbr(species = "Homo sapiens", category = "H") %>%
  mutate(gs_name=gsub("HALLMARK_","",gs_name))
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)


# 3. fgsea analysis
fgsea_res_all <- data.frame(pathway=character(), pval=double(),
                            padj=double(), ES=double(),
                            NES=double(), nMoreExtreme=double(),
                            size=integer(), leadingEdge=list(),
                            time=character()
)

## for each group separately:
time_groups <- unique(sc.genes$group)
for(time_point in time_groups){
  cluster0.genes<- sc.genes %>%
    dplyr::filter(group == time_point) %>%
    arrange(desc(auc)) %>%
    dplyr::select(feature, auc)
  
  ranks<- deframe(cluster0.genes)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))%>%
    mutate(time=time_point)
  
  fgsea_res_all <- fgsea_res_all %>% bind_rows(fgseaResTidy)  
  
}

#plot
fgsea_res_all_plot <-  fgsea_res_all %>%
  filter(padj < 0.05) #%>%
#head(n= 200)


obj_plot<-ggplot(fgsea_res_all_plot, aes(time, pathway, fill= NES)) +
  geom_tile(aes(fill = NES),colour = "white")+
  scale_fill_gradient2(low = "dodgerblue4", high = "darkred", mid = "white",
                       midpoint = 0, space = "Lab",
                       name="NES") +
  coord_fixed(ratio=0.5)+
  ggtitle("4H_Treated") +
  xlab(x_lab_name) + ylab("Pathways")
obj_plot

#Annotation of cell types
library(celldex)
library(SingleR)

#In my case, I am using human primary cell atlas data as a reference
hpca.se <- celldex::HumanPrimaryCellAtlasData()

#treated
DefaultAssay(sc) <- "RNA"
counts <- GetAssayData(sc)
pred.hesc <- SingleR(test = counts, ref = hpca.se , labels = hpca.se$label.fine)
cell_types <- pred.hesc$first.labels

#Add cell type as metadata
sc <- AddMetaData(sc, cell_types, col.name = "cell_type")
#Visualization of the annotated data
p1 <- DimPlot(sc, reduction = "umap", group.by = "cell_type", label.size = 2)

#Add title of the plot and modify the appearance
p1 + ggtitle("Types of Cell") + 
  theme(plot.title = element_text(size=15),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
           legend.text = element_text(size=10))
