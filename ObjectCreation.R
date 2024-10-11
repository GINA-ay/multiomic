#Load library
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

#Load matrix generating from cellranger 
MCMV_D7_raw<-Read10X('/home/twyang2/cellranger/counts/MCMV_D7/outs/filtered_feature_bc_matrix/') 
MCMV_D90_raw<-Read10X('/home/twyang2/cellranger/counts/MCMV_D90/outs/filtered_feature_bc_matrix/') 
LCMV_Arm_D7_raw<-Read10X('/home/twyang2/cellranger/counts/LCMV_Arm_D7/outs/filtered_feature_bc_matrix/') 
LCMV_Arm_D21_raw<-Read10X('/home/twyang2/cellranger/counts/LCMV_Arm_D21/outs/filtered_feature_bc_matrix/') 
LCMV_C13_D7_raw<-Read10X('/home/twyang2/cellranger/counts/LCMV_C13_D7/outs/filtered_feature_bc_matrix/') 
LCMV_C13_D21_raw<-Read10X('/home/twyang2/cellranger/counts/LCMV_C13_D21/outs/filtered_feature_bc_matrix/')

#Create seurat objects of gene expression data
MCMV_D7<-CreateSeuratObject(counts = MCMV_D7_raw$'Gene Expression',project="MCMV_D7")
MCMV_D90<-CreateSeuratObject(counts = MCMV_D90_raw$'Gene Expression',project='MCMV_D90') 
LCMV_Arm_D7<-CreateSeuratObject(counts = LCMV_Arm_D7_raw$'Gene Expression',project="LCMV_Arm_D7") 
LCMV_Arm_D21<-CreateSeuratObject(counts = LCMV_Arm_D21_raw$'Gene Expression',project="LCMV_Arm_D21") 
LCMV_C13_D7<-CreateSeuratObject(counts = LCMV_C13_D7_raw$'Gene Expression',project="LCMV_C13_D7") 
LCMV_C13_D21<-CreateSeuratObject(counts = LCMV_C13_D21_raw$'Gene Expression',project="LCMV_C13_D21")

#A list combining all six objects
seurat_list <- list(MCMV_D7 = MCMV_D7,
MCMV_D90 = MCMV_D90,
LCMV_Arm_D7 = LCMV_Arm_D7,
LCMV_Arm_D21 = LCMV_Arm_D21,
LCMV_C13_D7 = LCMV_C13_D7,
LCMV_C13_D21 = LCMV_C13_D21)


#peaks calling data generating from macs2
MCMV_D7_MACS2_peaks <-read.table("/home/twyang2/macs2/MCMV_D7/MCMV_D7_peaks.narrowPeak",header = FALSE) 
MCMV_D90_MACS2_peaks <-read.table("/home/twyang2/macs2/MCMV_D90/MCMV_D90_peaks.narrowPeak",header = FALSE) 
LCMV_Arm_D7_MACS2_peaks <-read.table("/home/twyang2/macs2/LCMV_Arm_D7/LCMV_Arm_D7_peaks.narrowPeak",header = FALSE) 
LCMV_Arm_D21_MACS2_peaks <-read.table("/home/twyang2/macs2/LCMV_Arm_D21/LCMV_Arm_D21_peaks.narrowPeak",header = FALSE) 
LCMV_C13_D7_MACS2_peaks <-read.table("/home/twyang2/macs2/LCMV_C13_D7/LCMV_C13_D7_peaks.narrowPeak",header = FALSE) 
LCMV_C13_D21_MACS2_peaks <-read.table("/home/twyang2/macs2/LCMV_C13_D21/LCMV_C13_D21_peaks.narrowPeak",header = FALSE)

#peaks annotation data
annotate <- GetGRangesFromEnsDb(EnsDb.Mmusculus.v79) 
seqlevels(annotate)<-paste0("chr",seqlevels(annotate)) 
genome(annotate)<"mm10"
MCMV_D7_gr <- GRanges(seqnames = MCMV_D7_MACS2_peaks$V1,ranges = IRanges(start = MCMV_D7_MACS2_peaks$V2,end = MCMV_D7_MACS2_peaks$V3)) 
MCMV_D90_gr <- GRanges(seqnames = MCMV_D90_MACS2_peaks$V1,ranges = IRanges(start = MCMV_D90_MACS2_peaks$V2,end = MCMV_D90_MACS2_peaks$V3)) 
LCMV_Arm_D7_gr <- GRanges(seqnames = LCMV_Arm_D7_MACS2_peaks$V1, ranges = IRanges(start = LCMV_Arm_D7_MACS2_peaks$V2,end = LCMV_Arm_D7_MACS2_peaks$V3)) 
LCMV_Arm_D21_gr <- GRanges(seqnames = LCMV_Arm_D21_MACS2_peaks$V1, ranges = IRanges(start = LCMV_Arm_D21_MACS2_peaks$V2, end = LCMV_Arm_D21_MACS2_peaks$V3)) 
LCMV_C13_D7_gr <- GRanges(seqnames = LCMV_C13_D7_MACS2_peaks$V1, ranges = IRanges(start = LCMV_C13_D7_MACS2_peaks$V2, end = LCMV_C13_D7_MACS2_peaks$V3)) 
LCMV_C13_D21_gr <- GRanges(seqnames = LCMV_C13_D21_MACS2_peaks$V1, ranges = IRanges(start = LCMV_C13_D21_MACS2_peaks$V2, end = LCMV_C13_D21_MACS2_peaks$V3)) 
combined_gr <- reduce(x = c(MCMV_D7_gr, MCMV_D90_gr,LCMV_Arm_D7_gr,
LCMV_Arm_D21_gr,LCMV_C13_D7_gr,LCMV_C13_D21_gr)) 
combined_gr <- combined_gr[width(combined_gr) < 10000 & width(combined_gr) > 20]

#peaks fragments
MCMV_D7_fragpath<-"/home/twyang2/cellranger/counts/MCMV_D7/outs/atac_fragments.tsv.gz" 
MCMV_D90_fragpath<-"/home/twyang2/cellranger/counts/MCMV_D90/outs/atac_fragments.tsv.gz" 
LCMV_Arm_D7_fragpath<-"/home/twyang2/cellranger/counts/LCMV_Arm_D7/outs/atac_fragments.tsv.gz" 
LCMV_Arm_D21_fragpath<-"/home/twyang2/cellranger/counts/LCMV_Arm_D21/outs/atac_fragments.tsv.gz" 
LCMV_C13_D7_fragpath<-"/home/twyang2/cellranger/counts/LCMV_C13_D7/outs/atac_fragments.tsv.gz" 
LCMV_C13_D21_fragpath<-"/home/twyang2/cellranger/counts/LCMV_C13_D21/outs/atac_fragments.tsv.gz" 

#frags path list 
fragpaths <- list(MCMV_D7=MCMV_D7_fragpath,
MCMV_D90=MCMV_D90_fragpath,
LCMV_Arm_D7=LCMV_Arm_D7_fragpath,
LCMV_Arm_D21=LCMV_Arm_D21_fragpath,
LCMV_C13_D7=LCMV_C13_D7_fragpath,
LCMV_C13_D21=LCMV_C13_D21_fragpath )

#Create ATAC assays in exsiting seurat object
for (name in names(fragpaths)) {
frags <- CreateFragmentObject( path = fragpaths[[name]], cells = colnames(seurat_list[[name]]))
counts <- FeatureMatrix(fragments = frags, features = combined_gr, sep = c(":", "-"), verbose = TRUE)
seurat_list[[name]][["peaks"]] <- CreateChromatinAssay(counts,fragments = frags,sep = c(":", "-"), annotation = annotate)}

#cell cycle genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- cc.genes$s.genes
s.genes <- str_to_title(tolower(s.genes))
g2m.genes <- str_to_title(tolower(g2m.genes))

#Preprocessing 
for (name in names(seurat_list)) {
df <- seurat_list[[name]]
DefaultAssay(df)<- "RNA"
##Quality Control and Filtering
mt_genes <- grep("mt-", rownames(df), value = TRUE)
df[["percent.mt"]] <- PercentageFeatureSet(df, features = mt_genes)
df <- subset(df,subset=nFeature_RNA<2000&nCount_RNA<5000&percent.mt<10&nFeature_peaks<20000&nCount_peaks<20000)
df<-NormalizeData(df)
#
existing_genes <- rownames(df)
s.genes_int <- intersect(s.genes, existing_genes)
g2m.genes_int <- intersect(g2m.genes, existing_genes)
df <- CellCycleScoring(df, s.features = s.genes_int, g2m.features = g2m.genes_int, set.ident = TRUE)
df<-SCTransform(df,vars.to.regress=c("S.Score","G2M.Score"))
df<-RunPCA(df,features=VariableFeatures(object = df))
df<-FindNeighbors(df,reduction="pca")
df<-FindClusters(df,cluster.name = "rna_clusters")
df <- RunUMAP(df,dims = 1:30, reduction = 'pca',reduction.name = "umap.pca")
sweep.res.list <- paramSweep(df, PCs = 1:30,sct=T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%as.character() %>%as.numeric()
DoubletRate <- 0.049 
homotypic.prop <- modelHomotypic(df$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(df))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
df <- doubletFinder(df,PCs = 1:30,pN = 0.25,pK = pK_bcmvn,nExp=nExp_poi.adj,reuse.pANN = F,sct = T)
last_meta_column_values <- df@meta.data[,ncol(df@meta.data)]
singlets <- df@meta.data[,ncol(df@meta.data)] == "Singlet"
df <- subset(df, cells = which(singlets))
seurat_list[[name]] <- df
}
#lsi normalisation
for (name in names(seurat_list)) {
df <- seurat_list[[name]]
DefaultAssay(df)<- "peaks"
df <- RunTFIDF(df)
df <- FindTopFeatures(df, min.cutoff = 20)
df <- RunSVD(df)
df <- RunUMAP(df, dims = 1:30, reduction = 'lsi',reduction.name = "umap.lsi")
df<-FindNeighbors(df,reduction="lsi")
df<-FindClusters(df,cluster.name = "atac_clusters")
seurat_list[[name]] <- df
}
#save objects from list
MCMV_D7 <- seurat_list$MCMV_D7
MCMV_D90 <- seurat_list$MCMV_D90
LCMV_Arm_D7 <- seurat_list$LCMV_Arm_D7
LCMV_Arm_D21 <- seurat_list$LCMV_Arm_D21
LCMV_C13_D7 <- seurat_list$LCMV_C13_D7
LCMV_C13_D21 <- seurat_list$LCMV_C13_D21

#merged 6 objects
merge_list<-list(MCMV_D90,LCMV_Arm_D7,LCMV_Arm_D21,LCMV_C13_D7,LCMV_C13_D21) 
merged<- merge(x=MCMV_D7,y=merge_list,project = "Merged")
merged$CC.Difference <- merged$S.Score - merged$G2M.Score
#RNA Assay Processing
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged,vars.to.regress="CC.Difference")
merged <- RunPCA(merged) 
merged <- FindNeighbors(merged,reduction="pca") 
merged <- FindClusters(merged,cluster.name = "rna_clusters",resolution=1)
merged <- RunUMAP(merged,dims = 1:30,reduction.name = "umap.rna")
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])
#ATAC Assay Processing
DefaultAssay(merged) <- "peaks"
merged <- RunTFIDF(merged)
merged <- FindTopFeatures(merged, min.cutoff = 20)
merged <- RunSVD(merged)
merged <- RunUMAP(merged, dims = 2:30, reduction = 'lsi',reduction.name = "umap.lsi")
merged<-FindNeighbors(merged,reduction="lsi")
merged<-FindClusters(merged,resolution=1,cluster.name = "atac_clusters")
#Activity Assay Processing
gene.activities <- GeneActivity(merged)
merged[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(merged) <- "ACTIVITY"
merged <- NormalizeData(
  object = merged,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(merged$nCount_ACTIVITY)
)
#Transfer labels
transfer_anchors <- FindTransferAnchors(
  reference = merged,
  query = merged,
  reference.assay = "RNA", 
  query.assay = "ACTIVITY",
  reduction = 'cca',
  dims = 1:30
)
clusters_predictions <- TransferData(anchorset = transfer_anchors, refdata = merged$rna_clusters, weight.reduction = merged[["lsi"]], dims = 2:30)
merged <- AddMetaData(merged, metadata = clusters_predictions)
merged$prediction_correct <- merged$predicted.id == merged$rna_clusters
atac_predict <- DimPlot(merged, group.by = "predicted.id", label = TRUE,reduction="umap.lsi") + NoLegend() + ggtitle("Predicted annotation")
atac_truth <- DimPlot(merged, group.by = "rna_clusters", label = TRUE,reduction="umap.lsi") + NoLegend() + ggtitle("Ground-truth annotation")

#Joint UMAP
merged <- FindMultiModalNeighbors(
  object = merged,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
merged <- RunUMAP(
  object = merged,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE,
  reduction.name = "umap.joint"
)
Joint_umap<-DimPlot(merged, label = TRUE, repel = TRUE, reduction = "umap.joint") + NoLegend()



