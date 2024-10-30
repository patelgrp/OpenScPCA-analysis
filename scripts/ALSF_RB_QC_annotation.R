library(Seurat)
library(lisi)
library(infercnv)
library(dplyr)
library(SingleCellExperiment)
library(scCustomize)
options(Seurat.object.assay.version = "v5")

#define folder with unzipped datasets from ALSF scPCA retinoblastoma cohort (SCPCP000011)
dir <- "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/data/RB/SCPCP000011_SINGLE-CELL_SINGLE-CELL-EXPERIMENT_2024-10-16/"

#pull metadata table from ALSF scPCA file
sample.metadata <- read.table(
  file = paste0(dir, "single_cell_metadata.tsv"),
  header = T,
  sep = "\t"
)

#pull each dataset's SingleCellExperiment file (already filtered and auto-annotated by scPCA)
#create an empty list for sce and Seurat object
sample.seurat.object <- list()
cntmatrix.list <- list()

#use a for loop to cycle through the samples, convert ensembl ids to gene symbols
for (i in 1:length(rownames(sample.metadata)))
{
  sample.sce.object <- readRDS(
    file = paste0(
      dir,
      sample.metadata$scpca_sample_id[i],
      "/",
      sample.metadata$scpca_library_id[i],
      "_processed.rds"
    )
  )
  
  #convert the singlecellexperiment into a Seurat object - first by selecting genes that have approved
  #gene symbols from HGNC
  newgenes <- rowData(sample.sce.object)[!is.na(rowData(sample.sce.object)$gene_symbol), ]
  newgenes <- newgenes[!duplicated(newgenes$gene_symbol), ]
  sample.sce.object <- sample.sce.object[newgenes$gene_ids, ]
  rownames(sample.sce.object) <- rowData(sample.sce.object)$gene_symbol
  #add library id to avoid cell barcode overlap
  colnames(sample.sce.object) <- paste(sample.metadata$scpca_library_id[i],
                                       colnames(sample.sce.object),
                                       sep = "_")
  
  #convert sce to seurat object
  sample.seurat.object[[i]] <- as.Seurat(sample.sce.object)
  
  #assign metadata to the seurat objects (probably overkill, but hopefully will be useful later!)
  sample.seurat.object[[i]]$scpca_project_id <- sample.metadata$scpca_project_id[i]
  sample.seurat.object[[i]]$scpca_sample_id <- sample.metadata$scpca_sample_id[i]
  sample.seurat.object[[i]]$scpca_library_id <- sample.metadata$scpca_library_id[i]
  sample.seurat.object[[i]]$diagnosis <- sample.metadata$diagnosis[i]
  sample.seurat.object[[i]]$subdiagnosis <- sample.metadata$subdiagnosis[i]
  sample.seurat.object[[i]]$disease_timing <- sample.metadata$disease_timing[i]
  sample.seurat.object[[i]]$age <- sample.metadata$age[i]
  sample.seurat.object[[i]]$age_timing <- sample.metadata$age_timing[i]
  sample.seurat.object[[i]]$sex <- sample.metadata$sex[i]
  sample.seurat.object[[i]]$tissue_location <- sample.metadata$tissue_location[i]
  sample.seurat.object[[i]]$participant_id <- sample.metadata$participant_id[i]
  sample.seurat.object[[i]]$submitter_id <- sample.metadata$submitter_id[i]
  sample.seurat.object[[i]]$organism <- sample.metadata$organism[i]
  sample.seurat.object[[i]]$development_stage_ontology_term_id <- sample.metadata$development_stage_ontology_term_id[i]
  sample.seurat.object[[i]]$sex_ontology_term_id <- sample.metadata$sex_ontology_term_id[i]
  sample.seurat.object[[i]]$organism_ontology_id <- sample.metadata$organism_ontology_id[i]
  sample.seurat.object[[i]]$self_reported_ethnicity_ontology_term_id <- sample.metadata$self_reported_ethnicity_ontology_term_id[i]
  sample.seurat.object[[i]]$disease_ontology_term_id <- sample.metadata$disease_ontology_term_id[i]
  sample.seurat.object[[i]]$tissue_ontology_term_id <- sample.metadata$tissue_ontology_term_id[i]
  sample.seurat.object[[i]]$primary_site <- sample.metadata$primary_site[i]
  sample.seurat.object[[i]]$sample_type <- sample.metadata$sample_type[i]
  sample.seurat.object[[i]]$seq_unit <- sample.metadata$seq_unit[i]
  sample.seurat.object[[i]]$technology <- sample.metadata$technology[i]
  sample.seurat.object[[i]]$total_reads <- sample.metadata$total_reads[i]
  sample.seurat.object[[i]]$mapped_reads <- sample.metadata$mapped_reads[i]
  sample.seurat.object[[i]]$unfiltered_cells <- sample.metadata$unfiltered_cells[i]
  sample.seurat.object[[i]]$filtered_cell_count <- sample.metadata$filtered_cell_count[i]
  sample.seurat.object[[i]]$processed_cells <- sample.metadata$processed_cells[i]
  sample.seurat.object[[i]]$has_cellhash <- sample.metadata$has_cellhash[i]
  sample.seurat.object[[i]]$includes_anndata <- sample.metadata$includes_anndata[i]
  sample.seurat.object[[i]]$is_cell_line <- sample.metadata$is_cell_line[i]
  sample.seurat.object[[i]]$is_multiplexed <- sample.metadata$is_multiplexed[i]
  sample.seurat.object[[i]]$is_xenograft <- sample.metadata$is_xenograft[i]
  sample.seurat.object[[i]]$pi_name <- sample.metadata$pi_name[i]
  sample.seurat.object[[i]]$project_title <- sample.metadata$project_title[i]
  sample.seurat.object[[i]]$genome_assembly <- sample.metadata$genome_assembly[i]
  sample.seurat.object[[i]]$mapping_index <- sample.metadata$mapping_index[i]
  sample.seurat.object[[i]]$alevin_fry_version <- sample.metadata$alevin_fry_version[i]
  sample.seurat.object[[i]]$salmon_version <- sample.metadata$salmon_version[i]
  sample.seurat.object[[i]]$transcript_type <- sample.metadata$transcript_type[i]
  sample.seurat.object[[i]]$droplet_filtering_method <- sample.metadata$droplet_filtering_method[i]
  sample.seurat.object[[i]]$cell_filtering_method <- sample.metadata$cell_filtering_method[i]
  sample.seurat.object[[i]]$prob_compromised_cutoff <- sample.metadata$prob_compromised_cutoff[i]
  sample.seurat.object[[i]]$min_gene_cutoff <- sample.metadata$min_gene_cutoff[i]
  sample.seurat.object[[i]]$normalization_method <- sample.metadata$normalization_method[i]
  sample.seurat.object[[i]]$date_processed <- sample.metadata$date_processed[i]
  sample.seurat.object[[i]]$workflow <- sample.metadata$workflow[i]
  sample.seurat.object[[i]]$workflow_version <- sample.metadata$workflow_version[i]
  sample.seurat.object[[i]]$workflow_commit <- sample.metadata$workflow_commit[i]
  
  #create a list of the count matrices (useful for downstream inference of CNV)
  cntmatrix.list[[i]] <- sample.seurat.object[[i]]@assays$originalexp@counts
  
  remove(sample.sce.object, newgenes)
}

#Our overall strategy uses 3 pieces of information to annotate samples based on malignant status
#1. we do a simple merge (no batch correction or integration) to determine intermixing of cells and assign a
#LISI score (a la Korsunsky et al Nat Methods 2019). Non-malignant cells will intermix more than malignant cells
#2. we use auto-annotations (already computed by the scPCA pipeline using singleR and cellassign) to
#annotate immune and endothelial cells (which we can safely presume are not malignant)
#3. (most onerous) we use inferCNV using a normal tissue reference to identify cells/clusters with CNVs

#create a merge of all data - weird intermittent bug with indexing
# sample.merge <- merge(sample.seurat.object[[1]], sample.seurat.object[2:length(sample.seurat.object)], merge.dr = TRUE)
#
# #process merged data using standard Seurat pipeline
# sample.merge <- NormalizeData(sample.merge)
# sample.merge <- FindVariableFeatures(sample.merge)
# sample.merge <- ScaleData(sample.merge)
# sample.merge <- RunPCA(sample.merge, reduction.name = "merge.pca")
# sample.merge <- FindNeighbors(sample.merge, dims = 1:30, reduction = "merge.pca")
# sample.merge <- FindClusters(sample.merge, resolution = 0.2, cluster.name = "merge_cluster")
# sample.merge <- RunUMAP(
#   sample.merge,
#   dims = 1:30,
#   reduction = "merge.pca",
#   reduction.name = "merge.umap"
# )
#
# #compute LISI index for each cell/nucleus across sample mixing
# metadata_id <- sample.merge$scpca_library_id
# metadata_id <- data.frame(metadata_id)
# colnames(metadata_id) <- "scpca_library_id"
# merge.lisi <- compute_lisi(
#   X = Embeddings(sample.merge, reduction = "merge.umap"),
#   meta_data = metadata_id,
#   label_colnames = "scpca_library_id"
# )
# sample.merge <- AddMetaData(sample.merge, metadata = merge.lisi, col.name = "merge_lisi")
# saveRDS(sample.merge, file = "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/data/RB/SCPCP000011_SINGLE-CELL_SINGLE-CELL-EXPERIMENT_2024-10-16/merged RB Seurat objects.Rds")

sample.merge <- readRDS(file = "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/scripts/merged RB Seurat objects.Rds")
#start by running infercnv on each sample, using a fetal retina dataset as a reference
#import count matrix from developing retina dataset from Norrie et al Nat Comms 2021 (GSE116106)
# retina_ref <- Read10X(data.dir = "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/data/RB/GSE116106_retina_reference")
# retina_ref <- retina_ref[, sample(colnames(retina_ref), 5000)]
# saveRDS(retina_ref, file = "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/data/RB/retina_ref.RDS")

retina_ref <- readRDS(file = "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/data/RB/retina_ref.RDS")
dir.create(paste0(dir, "RB_infercnv/"))
remove(sample.seurat.object)
updated.seurat.list <- list()

for (i in 27:length(rownames(sample.metadata)))
{
  updated.seurat.list[[i]] <- subset(sample.merge,
                                     scpca_library_id == sample.metadata$scpca_library_id[i])
  combined_matrix <- cntmatrix.list[[i]]
  gene_subset <- intersect(rownames(combined_matrix), rownames(retina_ref))
  
  retina_ref_subset <- retina_ref[gene_subset, ]
  
  combined_matrix <- combined_matrix[gene_subset, ]
  combined_matrix <- cbind(combined_matrix, retina_ref_subset)
  write.table(
    combined_matrix,
    file = paste0(dir, sample.metadata$scpca_library_id[i], "_matrix.txt"),
     quote = F,
    col.names = colnames(combined_matrix)
  )
  
  reference_annotation <- data.frame(row.names = colnames(retina_ref_subset))
  reference_annotation$annotation = 'reference'
  tumor_annotation <- data.frame(row.names = Cells(updated.seurat.list[[i]]),
                                 annotation = updated.seurat.list[[i]]$cluster)
  annotation_table <- rbind(reference_annotation, tumor_annotation)
  write.table(
    annotation_table,
    file = paste0(dir, sample.metadata$scpca_library_id[i], "_annotation.txt"),
    col.names = F,
    quote = F,
    sep = "\t"
  )
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = paste0(dir, sample.metadata$scpca_library_id[i], "_matrix.txt"),
    gene_order_file = "/mnt/6TB/GITHUB_REPOS/OpenScPCA-analysis/data/hg38_gencode_v27.txt",
    annotations_file = paste0(dir, sample.metadata$scpca_library_id[i], "_annotation.txt"),
    ref_group_names = "reference",
  )
  
  infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = paste0(dir, "RB_infercnv/", sample.metadata$scpca_library_id[i], "/"),
    cluster_by_groups = TRUE,
    analysis_mode = "samples",
    denoise = TRUE,
    HMM = FALSE
  )
  
  pdf(
    paste0(
      dir,
      "RB_infercnv/",
      sample.metadata$scpca_library_id[i],
      "/umap_plots.pdf"
    )
  )
  print(DimPlot(
    updated.seurat.list[[i]],
    group.by = "cluster",
    reduction = "UMAP"
  ))
  print(
    DimPlot(
      updated.seurat.list[[i]],
      group.by = "singler_celltype_annotation",
      reduction = "UMAP",
      label = TRUE
    ) + NoLegend()
  )
  print(
    DimPlot(
      updated.seurat.list[[i]],
      group.by = "cellassign_celltype_annotation",
      reduction = "UMAP",
      label = TRUE
    ) + NoLegend()
  )
  print(DimPlot(
    updated.seurat.list[[i]],
    group.by = "is_xenograft",
    reduction = "UMAP"
  ))
  print(
    FeaturePlot(
      updated.seurat.list[[i]],
      features = c("VWF"),
      reduction = "UMAP",
      min.cutoff = "q5",
      max.cutoff = "q95",
      cols = c("lightgrey", "darkblue"),
      order = TRUE
    )
  )
  print(
    FeaturePlot(
      updated.seurat.list[[i]],
      features = c("PTPRC"),
      reduction = "UMAP",
      min.cutoff = "q5",
      max.cutoff = "q95",
      cols = c("lightgrey", "darkblue"),
      order = TRUE
    )
  )
  print(
    FeaturePlot(
      updated.seurat.list[[i]],
      features = c("COL1A1"),
      reduction = "UMAP",
      min.cutoff = "q5",
      max.cutoff = "q95",
      cols = c("lightgrey", "darkblue"),
      order = TRUE
    )
  )
  print(
    FeaturePlot(
      updated.seurat.list[[i]],
      features = c("merge_lisi"),
      reduction = "UMAP",
      min.cutoff = "q5",
      max.cutoff = "q95",
      cols = c("lightgrey", "darkblue"),
      order = TRUE
    )
  )
  dev.off()
  
  # setwd(paste0(dir, "RB_infercnv/"))
  # try(scevan.results <- pipelineCNA(
  #   count_mtx = updated.seurat.list[[i]]@assays$originalexp@counts,
  #   norm_cell = NULL,
  #   SUBCLONES = F,
  #   sample = paste0("SCEVAN_", sample.metadata$scpca_library_id[i]),
  #   par_cores = 4,
  # ))
  
  
  remove(
    infercnv_obj,
    combined_matrix,
    gene_subset,
    retina_ref_subset,
    reference_annotation,
    tumor_annotation,
    annotation_table
  )
}

reticulate::use_condaenv("scvi-env", required = TRUE)
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

temp_Seurat <- sample.merge
temp_Seurat[["RNA"]] <- as(object = temp_Seurat[["RNA"]], Class = "Assay")
merge_adata <- convertFormat(
  temp_Seurat,
  from = "seurat",
  to = "anndata",
  main_layer = "counts",
  drop_single_values = FALSE
)
remove(temp_Seurat)

#batches <- c('sample_id', 'Method')
scvi$model$SCVI$setup_anndata(
  merge_adata,
  batch_key = "sample_id",
  categorical_covariate_keys = list("Method")
)
model <- scvi$model$SCVI(adata = merge_adata)
model$train(accelerator = "gpu")#, max_epochs = 10L)
latent = model$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(sample.merge)

saveRDS(model, file = paste0("scVI_model.Rds"))
saveRDS(latent, file = paste0("/scVI_latent_reps.Rds"))

sample.merge[["integrated.scvi"]] <- CreateDimReducObject(embeddings = latent,
                                                          key = "scvi_",
                                                          assay = DefaultAssay(sample.merge))

sample.merge <- FindNeighbors(sample.merge, reduction = "integrated.scvi", dims = 1:10)
sample.merge <- FindClusters(sample.merge, resolution = 0.4, cluster.name = "scvi_clusters")
sample.merge <- RunUMAP(
  sample.merge,
  reduction = "integrated.scvi",
  dims = 1:10,
  reduction.name = "umap.scvi"
)
