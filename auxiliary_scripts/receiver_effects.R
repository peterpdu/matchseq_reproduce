library(tidyverse)

library(Seurat)
library(tibble)
library(glmGamPoi)

library(qs2)

# usage: Rscript receiver_effects.R adata.qs2 label_matrix.txt column_key:celltype

args <- commandArgs(trailingOnly=TRUE)

so.file <- args[[1]]
df.file <- args[[2]]
key_celltype <- args[[3]]

stopifnot("arg3 not correctly formatted: column_key:celltype"= grepl(":", key_celltype))
key_celltype <- strsplit(key_celltype, ":")[[1]]
CELLTYPE_KEY <- key_celltype[1]
celltype <- key_celltype[2]

obj <- qs_read(so.file)
if (inherits(obj, "Seurat")) {
  counts <- obj[['RNA']]$counts
  meta <- obj@meta.data
} else {
  counts <- obj$counts
  meta <- obj$md
}

df <- read.table(df.file, row.names=1, header=T, stringsAsFactors=TRUE, sep='\t')
all_sgrna <- colnames(df)
df$size_factor <- meta[row.names(df), 'pan_size_factor']
df$celltype <- meta[row.names(df), CELLTYPE_KEY]
df$tumor <- meta[row.names(df), 'tumor']

df <- df[!is.na(df$tumor),]
df <- df[df$tumor != "",]

# subset to celltype
print(celltype)
cidx <- which(df$celltype == celltype)
stopifnot("fewer than 100 cells of cell type found, stopping"= length(cidx) > 100)

# subset to valid cells
counts <- counts[, row.names(df)]
counts <- counts[, cidx]

# remove low pct genes
nnz <- tabulate(counts@i + 1, nbins=nrow(counts))
pct.cells <- nnz / dim(counts)[2]
valid_genes <- (pct.cells > 0.1)
valid_genes <- names(valid_genes[valid_genes])
counts <- counts[valid_genes,]

print('Filtered count matrix dimensions')
print(dim(counts))

# set up sgRNAs
valid_sgrna <- all_sgrna[colSums(df[cidx,all_sgrna]) > 0]
valid_targets <- setdiff(unique(sapply(strsplit(valid_sgrna, '_'), "[", 1)), 'safe')
safe_sgrna <- valid_sgrna[grepl('safe_', valid_sgrna)]
df$total_safe <- rowSums(df[, safe_sgrna])

tumor_f <- factor(as.character(df$tumor[cidx]))
groups <- split(seq_along(tumor_f), tumor_f)
counts_dense <- as.matrix(counts)

# fit contrast models
print('Fitting contrast models vs safe sgRNA')
for(target in valid_targets){
  print(target)
  target_sgrna <- valid_sgrna[grepl(paste0(target, '_'), valid_sgrna)]
  pb.outpath <- paste0(tools::file_path_sans_ext(df.file),'--', CELLTYPE_KEY, '.',celltype,'.',target, '_vs_safe.glm_gp_contrast.txt')
  if (file.exists(pb.outpath)) next

  f <- as.formula(paste0('~ 1 + ', paste0(c(target_sgrna), collapse=' + ', " + total_safe")))
  c <- paste0('(', paste0(target_sgrna, collapse=' + '), ') /', length(target_sgrna), ' - total_safe')
  tryCatch(
    {
      mm <- model.matrix(f, df[cidx, ])

      pb_design <- do.call(rbind, lapply(groups, function(idx) DelayedMatrixStats::colMeans2(mm, rows=idx)))
      colnames(pb_design) <- colnames(mm)
      rownames(pb_design) <- names(groups)
      pb_counts <- do.call(cbind, lapply(groups, function(idx) DelayedMatrixStats::rowSums2(counts_dense, cols=idx)))
      rownames(pb_counts) <- rownames(counts_dense)
      colnames(pb_counts) <- names(groups)

      fit <- glm_gp(pb_counts, design=pb_design, size_factors="normed_sum", ridge_penalty=NULL, on_disk=FALSE, verbose=TRUE)

      pb.res <- test_de(fit, contrast=c, verbose=TRUE)
      pb.res$celltype <- celltype
      pb.res$target <- target
      write.table(pb.res, pb.outpath, sep='\t', quote=F, col.names=NA)
    },
    error = function(e) {
      print("could not fit model")
      print(e)
      NA
    }
  )
}
