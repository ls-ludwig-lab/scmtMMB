## helper function: load10Xdata and get_QC
suppressMessages({
  #including 10X data loading and bridge integration
  obj.multi <- readRDS("~/Ludwig_lab/ref_data/obj.multi.rds")
  obj.rna.ext <- readRDS("~/Ludwig_lab/ref_data/obj.rna.ext.rds")
  obj.rna <- LoadH5Seurat("~/Ludwig_lab/ref_data/pbmc_multimodal.h5seurat")
})

# Get gene annotations
suppressWarnings(annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86))
# Change style to UCSC 
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

load10XData <- function(ID){
  counts <-  Read10X_h5(filename = paste0("../PBMC_large_data_files/input/", ID ,"/cellranger/filtered_peak_bc_matrix.h5"))
  fragpath <- paste0("../PBMC_large_data_files/input/", ID ,"/cellranger/fragments.tsv.gz")
  metadata <- read.csv(file= paste0("../PBMC_large_data_files/input/", ID ,"/cellranger/singlecell.csv"), header = TRUE, row.names = 1)
  
  # Create ChromatinAssay for ATAC data
  atac_pbmc_assay <- CreateChromatinAssay( 
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = fragpath,
    annotation = annotations)
  
  # Requantify query ATAC to have same features as multiome ATAC dataset
  requant_multiome_ATAC <- FeatureMatrix(
    fragments = Fragments(atac_pbmc_assay),
    features = granges(obj.multi[['ATAC']]),
    cells = Cells(atac_pbmc_assay))
  
  # Create assay with requantified ATAC data
  ATAC_assay <- CreateChromatinAssay(
    counts = requant_multiome_ATAC,
    fragments = fragpath,
    annotation = annotations)
  
  # Create Seurat sbject
  pbmc  <- CreateSeuratObject(counts = ATAC_assay,assay = 'ATAC', meta.data = metadata, project = ID)
  # pbmc[['peak.orig']] <- atac_pbmc_assay
  pbmc <- subset(pbmc, subset = nCount_ATAC < 3e4 & nCount_ATAC > 1000)
  
  # normalize query
  pbmc <- RunTFIDF(pbmc)
  
  bridge.anchor <- FindBridgeTransferAnchors(extended.reference = obj.rna.ext, 
                                             query = pbmc,
                                             reduction = "lsiproject",
                                             dims = 2:50)
  
  pbmc_mapped <- MapQuery(anchorset = bridge.anchor, 
                          reference = obj.rna, 
                          query = pbmc, 
                          refdata = list(
                            celltype.l1 = "celltype.l1",
                            celltype.l2 = "celltype.l2",
                            celltype.l3 = "celltype.l3"),
                          reduction.model = "wnn.umap")
  
  pbmc$predicted.celltype.l1 <- pbmc_mapped$predicted.celltype.l1
  pbmc$predicted.celltype.l2 <- pbmc_mapped$predicted.celltype.l2
  pbmc <- AddMetaData(pbmc, pbmc_mapped@reductions$ref.umap@cell.embeddings, 
                      col.name = c("refUMAP_1", "refUMAP_2"))
  
  return(pbmc)
}

# QC on original peaks
get_QC <- function(object, plotting = TRUE){
  
  SO <- object
  DefaultAssay(SO) <- 'ATAC'
  SO$pct_reads_in_peaks <- SO$peak_region_fragments / SO$passed_filters * 100
  SO <- TSSEnrichment(SO, fast = F)
  SO <- NucleosomeSignal(SO)
  return(SO)
}

```