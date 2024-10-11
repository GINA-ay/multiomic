library(Seurat)
library(Signac) 
library(dplyr)
library(DoubletFinder) 
library(GenomicRanges) 
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)
library(ggplot2)
library(ggrepel)

#1.All clusters atlas
#a.UMAP
DefaultAssay(merged) <- "RNA"
RNAseq_umap <- DimPlot(merged,label=T,label.size=5,group.by = "rna_clusters")
RNAseq_samples_umap <- DimPlot(merged,label=T,label.size=5,group.by = "orig.ident")
#b.cell numeration
RNA_cells <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::select(rna_clusters, orig.ident) %>%
  group_by(rna_clusters, orig.ident) %>%
  tally() %>%
  ungroup()
RNA_cells_plot <- ggplot(RNA_cells, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
RNA_cellspro_plot <- ggplot(RNA_cells, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "Cell Cluster", y = "Proportion of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))
ggsave("proportion_cells_rna.png", plot = RNA_cellspro_plot, width = 10, height = 8, dpi = 300)
ggsave("num_of_cells_rna.png", plot = RNA_cells_plot, width = 10, height = 8, dpi = 300)
#c.examination
PD1_exp<-VlnPlot(object = merged, features = "Pdcd1", group.by = 'orig.ident')
TOX_exp<-VlnPlot(object = merged, features = "Tox", group.by = 'orig.ident')
TCF1_exp<-VlnPlot(object = merged, features = "Tcf7", group.by = 'orig.ident')
Il7r_exp<-VlnPlot(object = merged, features = "Il7r", group.by = 'orig.ident')
Klrg1_exp<-VlnPlot(object = merged, features = "Klrg1", group.by = 'orig.ident')

#ATAC
DefaultAssay(merged) <- "peaks"
#a.
ATACseq_umap <- DimPlot(merged, label=T,label.size=5,group.by = "atac_clusters")
ATACseq_samples_umap <- DimPlot(merged,label=T,label.size=5,group.by = "orig.ident")
#b.
ATAC_cells <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::select(atac_clusters, orig.ident) %>%
  group_by(atac_clusters, orig.ident) %>%
  tally() %>%
  ungroup()
ATAC_plot <- ggplot(ATAC_cells, aes(x = atac_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ATAC_proplot <- ggplot(ATAC_cells, aes(x = atac_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "Cell Cluster", y = "Proportion of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"))
ggsave("proportion_cells_atac.png", plot = ATAC_proplot, width = 10, height = 8, dpi = 300)
ggsave("num_of_cells_atac.png", plot = ATAC_plot, width = 10, height = 8, dpi = 300)
#c.Linkpeaks(ACR)
merged <- RegionStats(merged, genome = BSgenome.Mmusculus.UCSC.mm10)
merged <- LinkPeaks(
  object = merged,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Pdcd1", "Tox","Tcf7","Il7r","Klrg1")
)
Idents(merged)<-merged$orig.ident
merged <- SortIdents(merged)
PD1_cov <- CoveragePlot(
  object = merged,
  region = "Pdcd1",
  features = "Pdcd1",
  expression.assay = "RNA",
  extend.upstream = 2000,
  extend.downstream = 10000
)
TOX_cov <- CoveragePlot(
  object = merged,
  region = "Tox",
  features = "Tox",
  expression.assay = "RNA",
  extend.upstream = 2000,
  extend.downstream = 20000
)
TCF1_cov <- CoveragePlot(
  object = merged,
  region = "Tcf7",
  features = "Tcf7",
  expression.assay = "RNA",
  extend.upstream = 2000,
  extend.downstream = 10000
)
Il7r_cov <- CoveragePlot(
  object = merged,
  region = "Il7r",
  features = "Il7r",
  expression.assay = "RNA",
  extend.upstream = 2000,
  extend.downstream = 10000
)
Klrg1_cov <- CoveragePlot(
  object = merged,
  region = "Klrg1",
  features = "Klrg1",
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 20000
)

#Predict
DefaultAssay(merged) <- "peaks"
predict_colors <- c( "0"="#FF6347", "1"="#4682B4", "2"="#32CD32", "3"="#FFD700", "4"="#DA70D6","5"="#8A2BE2", "6"="#00CED1", "7"="#FF4500","8"="#7FFF00","9"="#FF69B4", "10"="#6495ED", "11"="#DC143C","12"="#00FA9A","13"="#BA55D3","14"="#1E90FF","15"="#FF8C00","False"="#BEBEBE")
merged$predicted <- ifelse(merged$prediction_correct==TRUE, as.character(merged$predicted.id),"False")
Activity_umap <- DimPlot(merged, label = TRUE, group.by = "predicted",cols = predict_colors)

#Annotate
DefaultAssay(merged)<-"RNA"
Idents(merged)<-"rna_clusters"
new.cluster.ids<-c("TE","Teff-infla","Tem","Teff","Texterm","Tscm","Teff-like"," "," ","Texprog"," ","Texint","Trans I"," "," "," ")
names(new.cluster.ids)<-levels(merged)
merged<-RenameIdents(merged,new.cluster.ids)
RNAseq_umap_annotated <- DimPlot(merged,label=T,label.size=3)
#
DefaultAssay(merged)<-"peaks"
Idents(merged)<-"atac_clusters"
atac_cluster_ids<-c("Tem","Teff","Teff-like","TE","Tscm","Texterm","Teff-infla","Teff","TE","Teff-infla","Teff-infla","Trans I","Tem","Texprog","Texint","Teff-like","Teff-like","Texterm")
names(atac_cluster_ids)<-levels(merged)
merged<-RenameIdents(merged, atac_cluster_ids)
ATACseq_umap_annotated <- DimPlot(merged,label=T,label.size=3)

#2.LCMV Arm 
#scRNAseq
DefaultAssay(merged) <- "RNA"
#
LCMV_Arm_RNA_clusters <- c("3","7","10","2","5","0") 
LCMV_Arm_samples <- c("LCMV_Arm_D7", "LCMV_Arm_D21")
#a.UMAP
merged$LCMV_Arm_RNAseq <- ifelse(merged$orig.ident %in% LCMV_Arm_samples & merged$rna_clusters %in% LCMV_Arm_RNA_clusters, as.character(merged$rna_clusters), "Other")
LCMV_Arm_RNA_colors <- c("3" = "#6BAED6", "5" = "#B0E0E6", "0" = "#E0FFFF", "2" = "#EAB8E4","7" = "#E6E6FA", "10"="#A0C4FF","Other" = "#BEBEBE")
LCMV_Arm_RNAsequmap <- DimPlot(merged, label = TRUE, group.by = "LCMV_Arm_RNAseq", cols = LCMV_Arm_colors)
#
merged$LCMV_Arm_cycle <- ifelse(merged$orig.ident %in% LCMV_Arm_samples & merged$rna_clusters %in% LCMV_Arm_RNA_clusters, as.character(merged$Phase), "Other")
LCMV_Arm_cycle_colors <- c("G1" = "#6BAED6", "S" = "#E6E6FA", "G2M" = "#E0FFFF","Other" = "#BEBEBE")
LCMV_Arm_cycle_umap <- DimPlot(merged, label = TRUE, group.by = "LCMV_Arm_cycle", cols = LCMV_Arm_cycle_colors)
#b.enumration   
LCMV_Arm_data <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(rna_clusters %in% LCMV_Arm_RNA_clusters & orig.ident %in% LCMV_Arm_samples) %>%
  dplyr::select(rna_clusters, orig.ident) %>%
  group_by(rna_clusters, orig.ident) %>%
  tally() %>%
  ungroup()
LCMV_Arm_plot <- ggplot(LCMV_Arm_data, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_Arm_num_of_cells.png", plot = LCMV_Arm_plot, width = 10, height = 8, dpi = 300)
#
LCMV_Arm_proplot <- ggplot(LCMV_Arm_data, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell Cluster", y = "Proportion of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_Arm_pro_of_cells.png", plot = LCMV_Arm_proplot, width = 10, height = 8, dpi = 300)
#cellcycle
LCMV_Arm_cycle_cells <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(rna_clusters %in% LCMV_Arm_RNA_clusters & orig.ident %in% LCMV_Arm_samples) %>%
  dplyr::select(rna_clusters, Phase) %>%
  group_by(rna_clusters, Phase) %>%
  tally() %>%
  ungroup()
LCMV_Arm_cycle_plot <- ggplot(LCMV_Arm_cycle_cells, aes(x = rna_clusters, y = n, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_Arm_num_of_cells_cycle.png", plot = LCMV_Arm_cycle_plot, width = 10, height = 8, dpi = 300)
#c.DEG+ACR
Idents(merged)<-"rna_clusters"
LCMV_Arm_Cx3cr1_exp<-VlnPlot(object = merged, features = "Cx3cr1", idents = LCMV_Arm_RNA_clusters)
LCMV_Arm_TCF1_exp<-VlnPlot(object = merged, features = "Tcf7", idents = LCMV_Arm_RNA_clusters)
LCMV_Arm_Klrg1_exp<-VlnPlot(object = merged, features = "Klrg1", idents = LCMV_Arm_RNA_clusters)
LCMV_Arm_Il7r_exp<-VlnPlot(object = merged, features = "Il7r", idents = LCMV_Arm_RNA_clusters)
LCMV_Arm_Cxcr3_exp<-VlnPlot(object = merged, features = "Cxcr3", idents = LCMV_Arm_RNA_clusters)
LCMV_Arm_CD27_exp<-VlnPlot(object = merged, features = "Cd27", idents = LCMV_Arm_RNA_clusters)
LCMV_Arm_S1pr5_exp<-VlnPlot(object = merged, features = "S1pr5", idents = LCMV_Arm_RNA_clusters)
#
DefaultAssay(merged) <- "peaks"
Idents(merged)<-merged$rna_clusters
merged <- SortIdents(merged)
merged <- RegionStats(merged, genome = BSgenome.Mmusculus.UCSC.mm10)
merged <- LinkPeaks(
  object = merged,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Cx3cr1", "Zeb2","Tcf7","Klrg1","Il7r"))
LCMV_Arm_Cx3cr1_cov <- CoveragePlot(
  object = merged,
  region = "Cx3cr1",
  features = "Cx3cr1",
  idents = LCMV_Arm_RNA_clusters,
  expression.assay = "RNA",
  extend.upstream = 3000,
  extend.downstream = 10000
)
LCMV_Arm_Zeb2_cov <- CoveragePlot(
  object = merged,
  region = "Zeb2",
  features = "Zeb2",
  idents = LCMV_Arm_RNA_clusters,
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
)
LCMV_Arm_TCF1_cov <- CoveragePlot(
  object = merged,
  region = "Tcf7",
  features = "Tcf7",
  idents = LCMV_Arm_RNA_clusters,
  expression.assay = "RNA",
  extend.upstream = 3000,
  extend.downstream = 10000
)
LCMV_Arm_Klrg1_cov <- CoveragePlot(
  object = merged,
  region = "Klrg1",
  features = "Klrg1",
  idents = LCMV_Arm_RNA_clusters,
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
)
LCMV_Arm_Il7r_cov <- CoveragePlot(
  object = merged,
  region = "Il7r",
  features = "Il7r",
  idents = LCMV_Arm_RNA_clusters,
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
)
#DEG
DefaultAssay(merged) <- "RNA"
#
cluster3_10.markers <- FindMarkers(merged, ident.1 = 3,ident.2 = 10)
cluster3_0.markers <- FindMarkers(merged, ident.1 = 3,ident.2 = 0)
cluster0_2.markers <- FindMarkers(merged, ident.1 = 0,ident.2 = 2)
cluster2_5.markers <- FindMarkers(merged, ident.1 = 2,ident.2 = 5)
cluster10_7.markers <- FindMarkers(merged, ident.1 = 10,ident.2 = 7)
cluster10_0.markers <- FindMarkers(merged, ident.1 = 10,ident.2 = 0)
cluster7_0.markers <- FindMarkers(merged, ident.1 = 7,ident.2 = 0)


#ATACseq
LCMV_Arm_samples <- c("LCMV_Arm_D7", "LCMV_Arm_D21")
LCMV_Arm_ATAC_clusters <- c("0","1","4","7","10","11","12")
#a.UMAP
DefaultAssay(merged) <- "peaks"
#
merged$LCMV_Arm_ATACseq <- ifelse(merged$orig.ident %in% LCMV_Arm_samples & merged$atac_clusters %in% LCMV_Arm_ATAC_clusters, as.character(merged$atac_clusters), "Other")
LCMV_Arm_ATAC_colors <- c("0" = "#6BAED6", "1" = "#B0E0E6", "4" = "#EAB8E4","7" = "#A0C4FF", "10" = "#E0FFFF", "12" = "#EE82EE","11" = "#E6E6FA", "Other" = "#BEBEBE")
LCMV_Arm_ATACseq_umap <- DimPlot(merged, label = TRUE, group.by = "LCMV_Arm_ATACseq", cols = LCMV_Arm_colors) 
#b.
LCMV_Arm_ATAC_data <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(atac_clusters %in% LCMV_Arm_ATAC_clusters & orig.ident %in% LCMV_Arm_samples) %>%
  dplyr::select(atac_clusters, orig.ident) %>%
  group_by(atac_clusters, orig.ident) %>%
  tally() %>%
  ungroup()
LCMV_Arm_ATAC_plot <- ggplot(LCMV_Arm_ATAC_data, aes(x = atac_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_Arm_num_of_cells_ATAC.png", plot = LCMV_Arm_ATAC_plot, width = 10, height = 8, dpi = 300)
#c.Gene activity
#d.DAG+ACR
DefaultAssay(merged) <- "peaks"
Idents(merged)<-"predicted.id"
LCMV_Arm_allACRs<-FindMarkers(merged,ident.1=c(2,5,0),ident.2=c(3,7,10)
#
DEG_data <- data.frame(
gene =rownames(LCMV_Arm_allACRs),
log2FoldChange = LCMV_Arm_allACRs$avg_log2FC,
pvalue = LCMV_Arm_allACRs$p_val,
pvalue_adj = LCMV_Arm_allACRs$p_val_adj)
DEG_data$pvalue_adj <- ifelse(DEG_data$pvalue_adj == 0, 1e-300, DEG_data$pvalue_adj)
DEG_data$pvalue <- ifelse(DEG_data$pvalue == 0, 1e-300, DEG_data$pvalue)
specific_data1 <- DEG_data %>% 
  dplyr::arrange(desc(abs(log2FoldChange))) %>% 
  head(5)
specific_data2 <- DEG_data %>% 
  dplyr::arrange(log2FoldChange) %>%
  head(5) 
open250<-specific_data$gene
open3710<-specific_data2$gene
closest_250 <- ClosestFeature(merged, open250)
closest_3710 <- ClosestFeature(merged, open3710)
#
BuildClusterTree(merged,reduction = "pca.lsi")

#3.LCMV Clone 13
library(ggrepel)
#LCMV C13
#
DefaultAssay(merged) <- "RNA"
LCMV_C13_rna_clusters <- c("6","8","12","13","4","9","11")
LCMV_C13_samples <- c("LCMV_C13_D7", "LCMV_C13_D21")
#a.UMAP
merged$LCMV_C13_RNAseq <- ifelse(merged$orig.ident %in% LCMV_C13_samples & merged$rna_clusters %in% LCMV_C13_rna_clusters, as.character(merged$rna_clusters), "Other")
LCMV_C13_rna_colors <- c("6" = "#E0FFFF", "8" = "#B0E0E6", "12" = "#6BAED6","13" = "#A0C4FF", "4" = "#EAB8E4", "9" = "#EE82EE","11" = "#E6E6FA", "Other" = "#BEBEBE")
LCMV_C13_RNAsequmap <- DimPlot(merged, label = TRUE, group.by = "LCMV_C13_RNAseq", cols = LCMV_C13_rna_colors)
#
merged$LCMV_C13_cycle <- ifelse(merged$orig.ident %in% LCMV_C13_samples & merged$rna_clusters %in% LCMV_C13_rna_clusters, as.character(merged$Phase), "Other")
LCMV_C13_cycle_colors <- c("G1" = "#6BAED6", "S" = "#E6E6FA", "G2M" = "#E0FFFF","Other" = "#BEBEBE")
LCMV_C13_cycle_umap <- DimPlot(merged, label = TRUE, group.by = "LCMV_C13_cycle", cols = LCMV_C13_cycle_colors)
#b.enumration   
LCMV_C13_data <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(rna_clusters %in% LCMV_C13_rna_clusters & orig.ident %in% LCMV_C13_samples) %>%
  dplyr::select(rna_clusters, orig.ident) %>%
  group_by(rna_clusters, orig.ident) %>%
  tally() %>%
  ungroup()
LCMV_C13_plot <- ggplot(LCMV_C13_data, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_C13_num_of_cells.png", plot = LCMV_C13_plot, width = 10, height = 8, dpi = 300)
#
LCMV_C13_plot_proplot <- ggplot(LCMV_C13_data, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell Cluster", y = "Proportion of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_C13_pro_of_cells.png", plot = LCMV_C13_plot_proplot, width = 10, height = 8, dpi = 300)
#
LCMV_C13_cycle_cells <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(rna_clusters %in% LCMV_C13_rna_clusters & orig.ident %in% LCMV_C13_samples) %>%
  dplyr::select(rna_clusters, Phase) %>%
  group_by(rna_clusters, Phase) %>%
  tally() %>%
  ungroup()
LCMV_C13_cycle_plot <- ggplot(LCMV_C13_cycle_cells, aes(x = rna_clusters, y = n, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("LCMV_C13_num_of_cells_cycle.png", plot = LCMV_C13_cycle_plot, width = 10, height = 8, dpi = 300)

#c.DEG+ACR
#
LCMV_C13_Ccr7_exp<-VlnPlot(object = merged, features = "Ccr7", idents = LCMV_C13_rna_clusters)
LCMV_C13_TCF1_exp<-VlnPlot(object = merged, features = "Tcf7", idents = LCMV_C13_rna_clusters)
LCMV_C13_Il7r_exp<-VlnPlot(object = merged, features = "Il7r", idents = LCMV_C13_rna_clusters)
LCMV_C13_Klrg1_exp<-VlnPlot(object = merged, features = "Klrg1", idents = LCMV_C13_rna_clusters)
LCMV_C13_PD1_exp<-VlnPlot(object = merged, features = "Pdcd1", idents = LCMV_C13_rna_clusters)
LCMV_C13_CD160_exp<-VlnPlot(object = merged, features = "Cd160", idents = LCMV_C13_rna_clusters)
LCMV_C13_ID3_exp<-VlnPlot(object = merged, features = "Id3", idents = LCMV_C13_rna_clusters)
LCMV_C13_Tox_exp<-VlnPlot(object = merged, features = "Tox", idents = LCMV_C13_rna_clusters)
LCMV_C13_Eomes_exp<-VlnPlot(object = merged, features = "Eomes", idents = LCMV_C13_rna_clusters)
LCMV_C13_Cx3cr1_exp<-VlnPlot(object = merged, features = "Cx3cr1", idents = LCMV_C13_rna_clusters)
LCMV_C13_Tbet_exp<-VlnPlot(object = merged, features = "Tbx21", idents = LCMV_C13_rna_clusters)
#DEG
cluster6_12.markers <- FindMarkers(merged, ident.1 = 6,ident.2 = 12)
cluster12_9.markers <- FindMarkers(merged, ident.1 = 12,ident.2 = 9)
cluster9_11.markers <- FindMarkers(merged, ident.1 = 9,ident.2 = 11)
cluster11_4.markers <- FindMarkers(merged, ident.1 = 11,ident.2 = 4)
cluster6_8.markers <- FindMarkers(merged, ident.1 = 6,ident.2 = 8)
cluster8_13.markers <- FindMarkers(merged, ident.1 = 8,ident.2 = 13)
cluster8_4.markers <- FindMarkers(merged, ident.1 = 8,ident.2 = 4)
cluster11_13.markers <- FindMarkers(merged, ident.1 = 11,ident.2 = 13)


#ATAC
DefaultAssay(merged) <- "peaks"
LCMV_C13_atac_clusters <- c("2","15","16","5","13","14","17")
#a.UMAP
merged$LCMV_C13_ATACseq <- ifelse(merged$orig.ident %in% LCMV_C13_samples & merged$atac_clusters %in% LCMV_C13_atac_clusters, as.character(merged$atac_clusters), "Other")
LCMV_C13_atac_colors <- c("2" = "#6BAED6", "15" = "#B0E0E6", "16" = "#EAB8E4","5" = "#A0C4FF", "13" = "#E0FFFF", "14" = "#EE82EE","17" = "#E6E6FA", "Other" = "#BEBEBE")
LCMV_C13_ATACsequmap <- DimPlot(merged, label = TRUE, group.by = "LCMV_C13_ATACseq", cols = LCMV_C13_atac_colors) 
#
DefaultAssay(merged) <- "peaks"
Idents(merged)<-"predicted.id"

#3.MCMV
library(ggrepel)
#scRNAseq
DefaultAssay(merged) <- "RNA"
#
MCMV_rna_clusters <- c("1","3","7","0","2","5") 
MCMV_samples <- c("MCMV_D7", "MCMV_D90")
#a.UMAP
merged$MCMV_RNAseq <- ifelse(merged$orig.ident %in% MCMV_samples & merged$rna_clusters %in% MCMV_rna_clusters, as.character(merged$rna_clusters), "Other")
MCMV_rna_colors <- c("1" = "#6BAED6", "3" = "#B0E0E6", "7" = "#E0FFFF", "0" = "#EAB8E4", "2"="#EE82EE","5"= "7EC8E3","Other" = "#BEBEBE")
MCMV_RNAsequmap <- DimPlot(merged, label = TRUE, group.by = "MCMV_RNAseq", cols = MCMV_rna_colors)
#
merged$MCMV_cycle <- ifelse(merged$orig.ident %in% MCMV_samples & merged$rna_clusters %in% MCMV_rna_clusters, as.character(merged$Phase), "Other")
MCMV_cycle_colors <- c("G1" = "#6BAED6", "S" = "#E6E6FA", "G2M" = "#E0FFFF","Other" = "#BEBEBE")
MCMV_cycle_umap <- DimPlot(merged, label = TRUE, group.by = "MCMV_cycle", cols = MCMV_cycle_colors)
#b.enumration     
MCMV_data <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(rna_clusters %in% MCMV_rna_clusters & orig.ident %in% MCMV_samples) %>%
  dplyr::select(rna_clusters, orig.ident) %>%
  group_by(rna_clusters, orig.ident) %>%
  tally() %>%
  ungroup()
MCMV_plot <- ggplot(MCMV_data, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("MCMV_num_of_cells.png", plot = MCMV_plot, width = 10, height = 8, dpi = 300)
#
MCMV_proplot <- ggplot(MCMV_data, aes(x = rna_clusters, y = n, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell Cluster", y = "Proportion of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("MCMV_pro_of_cells.png", plot = MCMV_proplot, width = 10, height = 8, dpi = 300)
#cellcycle
MCMV_cycle_cells <- merged@meta.data %>%
  as_tibble() %>%
  dplyr::filter(rna_clusters %in% MCMV_rna_clusters & orig.ident %in% MCMV_samples) %>%
  dplyr::select(rna_clusters, Phase) %>%
  group_by(rna_clusters, Phase) %>%
  tally() %>%
  ungroup()
MCMV_cycle_plot <- ggplot(MCMV_cycle_cells, aes(x = rna_clusters, y = n, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cell Cluster", y = "Number of Cells", fill = "Sample Type") +
  theme_minimal() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.background = element_rect(fill = "white"), panel.background = element_rect(fill = "white"))
ggsave("MCMV_num_of_cells_cycle.png", plot = MCMV_cycle_plot, width = 10, height = 8, dpi = 300)
#c.DEG+ACR
MCMV_Cx3cr1_exp<-VlnPlot(object = merged, features = "Cx3cr1", idents = MCMV_rna_clusters)
MCMV_TCF1_exp<-VlnPlot(object = merged, features = "Tcf7", idents = MCMV_rna_clusters)
MCMV_PD1_exp<-VlnPlot(object = merged, features = "Pdcd1", idents = MCMV_rna_clusters)
MCMV_Klrg1_exp<-VlnPlot(object = merged, features = "Klrg1", idents = MCMV_rna_clusters)
MCMV_Il7r_exp<-VlnPlot(object = merged, features = "Il7r", idents = MCMV_rna_clusters)
MCMV_Cxcr3_exp<-VlnPlot(object = merged, features = "Cxcr3", idents = MCMV_rna_clusters)
MCMV_Ccr7_exp<-VlnPlot(object = merged, features = "Ccr7", idents = MCMV_rna_clusters)
MCMV_Gzma_exp<-VlnPlot(object = merged, features = "Gzma", idents = MCMV_rna_clusters)
MCMV_CD27_exp<-VlnPlot(object = merged, features = "Cd27", idents = MCMV_rna_clusters)
MCMV_Tbet_exp<-VlnPlot(object = merged, features = "Tbx21", idents = MCMV_rna_clusters)
MCMV_Eomes_exp<-VlnPlot(object = merged, features = "Eomes", idents = MCMV_rna_clusters)
MCMV_TOX_exp<-VlnPlot(object = merged, features = "Tox", idents = MCMV_rna_clusters)
MCMV_S1pr5_exp<-VlnPlot(object = merged, features = "S1pr5", idents = MCMV_rna_clusters)

#
cluster1_3.markers <- FindMarkers(merged, ident.1 = 1,ident.2 = 3)
cluster3_7.markers <- FindMarkers(merged, ident.1 = 3,ident.2 = 7)
cluster1_0.markers <- FindMarkers(merged, ident.1 = 1,ident.2 = 0)
cluster0_2.markers <- FindMarkers(merged, ident.1 = 0,ident.2 = 2)
cluster2_5.markers <- FindMarkers(merged, ident.1 = 2,ident.2 = 5)
cluster1_5.markers <- FindMarkers(merged, ident.1 = 1,ident.2 = 5)
cluster0_3.markers <- FindMarkers(merged, ident.1 = 0,ident.2 = 3)
cluster7_5.markers <- FindMarkers(merged, ident.1 = 7,ident.2 = 5)
#ATAC
DefaultAssay(merged) <- "peaks"
Idents(merged)<-merged$rna_clusters
merged <- SortIdents(merged)
merged <- RegionStats(merged, genome = BSgenome.Mmusculus.UCSC.mm10)
merged <- LinkPeaks(
  object = merged,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("Cx3cr1", "Zeb2","Tcf7","Klrg1","Il7r"))
MCMV_Cx3cr1_cov <- CoveragePlot(
  object = merged,
  region = "Cx3cr1",
  features = "Cx3cr1",
  idents = MCMV_rna_clusters,
  expression.assay = "RNA",
  extend.upstream = 3000,
  extend.downstream = 10000
)
MCMV_Zeb2_cov <- CoveragePlot(
  object = merged,
  region = "Zeb2",
  features = "Zeb2",
  idents = MCMV_rna_clusters,
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
)
MCMV_TCF1_cov <- CoveragePlot(
  object = merged,
  region = "Tcf7",
  features = "Tcf7",
  idents = MCMV_rna_clusters,
  expression.assay = "RNA",
  extend.upstream = 3000,
  extend.downstream = 10000
)
MCMV_Klrg1_cov <- CoveragePlot(
  object = merged,
  region = "Klrg1",
  features = "Klrg1",
  idents = MCMV_rna_clusters,
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
)
MCMV_Il7r_cov <- CoveragePlot(
  object = merged,
  region = "Il7r",
  features = "Il7r",
  idents = MCMV_rna_clusters,
  expression.assay = "RNA",
  extend.upstream = 5000,
  extend.downstream = 5000
)
#DEG
DefaultAssay(merged) <- "RNA"
cluster0_4.markers <- FindMarkers(merged, ident.1 = 0,ident.2 = 4)
cluster2_3.markers <- FindMarkers(merged, ident.1 = 2,ident.2 = 3)
cluster0_2.markers <- FindMarkers(merged, ident.1 = 0,ident.2 = 2)
cluster3_7.markers <- FindMarkers(merged, ident.1 = 3,ident.2 = 7)
cluster14_8.markers <- FindMarkers(merged, ident.1 = 14,ident.2 = 8)
cluster14_3.markers <- FindMarkers(merged, ident.1 = 14,ident.2 = 3)


#ATACseq
MCMV_samples <- c("MCMV_D7", "MCMV_D90")
MCMV_atac_clusters <- c("6","9","10","3","8")
#a.UMAP
DefaultAssay(merged) <- "peaks"
#
merged$MCMV_ATACseq <- ifelse(merged$orig.ident %in% MCMV_samples & merged$atac_clusters %in% MCMV_atac_clusters, as.character(merged$atac_clusters), "Other")
MCMV_atac_colors <- c("6" = "#6BAED6", "9" = "#B0E0E6", "10" = "#EAB8E4","3" = "#A0C4FF", "8" = "#E0FFFF", "Other" = "#BEBEBE")
MCMV_ATACseq_umap <- DimPlot(merged, label = TRUE, group.by = "MCMV_ATACseq", cols = MCMV_atac_colors) 
#c.Gene activity
#d.ACR
DefaultAssay(merged) <- "peaks"
Idents(merged)<-"predicted.id"
MCMV_3_8_markers<-FindMarkers(merged,ident.1=3,ident.2=8)
BuildClusterTree(merged,reduction = "pca.lsi")

