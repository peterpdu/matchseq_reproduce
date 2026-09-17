library(Seurat)
library(dplyr)
library(Matrix)
library(qs2)
library(reticulate)

use_python("/home/groups/bassik/ppdu/miniconda3/envs/py3/bin/python")

np <- import('numpy')
sp <- import('scipy')
sc <- import('scanpy')

convertR2C <- function(x){
  as(as(x, 'CsparseMatrix'), 'dgCMatrix')
}

adata2seurat <- function(infile){
    print("Reading h5ad")
    data <- sc$read_h5ad(infile)
    
    print("Converting matrix")
    counts <- convertR2C(t(data$raw$X))
    cell_metadata <- data$obs
    gene_metadata <- data$var
    gene_metadata$gene_short_name <- row.names(gene_metadata)

    rownames(counts) <- row.names(gene_metadata)
    colnames(counts) <- row.names(cell_metadata)

    # store embeddings
    print("Storing embeddings")
    pca_emb <- data$obsm$get('X_pca')
    pca_loadings <- data$varm$get('PCs')
    row.names(pca_emb)<-row.names(cell_metadata)
    row.names(pca_loadings)<-row.names(gene_metadata)

    umap_emb <- data$obsm$get('X_umap')
    row.names(umap_emb)<-row.names(cell_metadata)
  
    # create object
    print("Creating object and normalizing")
    so <- CreateSeuratObject(counts, meta.data=cell_metadata)
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
    so <- ScaleData(so, features=row.names(so))
    so@meta.data <- cbind(so@meta.data, cell_metadata)
    all.genes <- rownames(so)
    so@reductions$pca <- CreateDimReducObject(embeddings=pca_emb, loadings=pca_loadings, key='PC_', assay='RNA')
    so@reductions$pca@stdev <- as.vector(sqrt(data$uns$get('pca')$variance))

    so@reductions$umap <- CreateDimReducObject(embeddings=umap_emb, key='UMAP_', assay='RNA')
    
    print("Writing to disk")
    # saveRDS(so, paste0(tools::file_path_sans_ext(infile), '_seurat.rds'))
    qs_save(so, paste0(tools::file_path_sans_ext(infile), '_seurat.qs2'))
    return(so)
}

args <- commandArgs(trailingOnly=TRUE)

adata2seurat(args[[1]])
