library(deMULTIplex2)

# needs a lot of memory for EM steps (<100GB)

# arg1 = fastq directory
# arg2 = valid cell barcodes text file, no header (get after GEX processing); will remove "-1" suffix if found
# arg3 = barcode matrix csv file (bar.ref.new.csv)
args <- commandArgs(TRUE)

# set valid barcodes
valid.barcodes <- c('A1', 'A2', 'A3', 'A4', 'B1', 'B2', 'B3', 'B4', 'B5')


# directory contains R1, R2, I1, I2 fastqs
# sample name will be the name of this directory
fastq.dir <- normalizePath(args[1])
name <- basename(fastq.dir)
valid.cells.file <- args[2]
outdir <- file.path(fastq.dir, basename(tools::file_path_sans_ext(valid.cells.file)))
dir.create(outdir)

valid.cells <- read.table(valid.cells.file, header=F)$V1
suffix.found <- any(grepl("-1$", valid.cells))
if (suffix.found) {
	cat('10X "-1" cell barcode suffix found, removing\n')
	valid.cells <- sub("-1$", "", valid.cells)
}

tag.ref <- read.csv(args[3])
# make named vector
bar.ref <- setNames(as.character(tag.ref$Barcode_Sequence), tag.ref$Well_Position)
valid.bar.ref <- bar.ref[valid.barcodes]

read_table <- readTags(dir=fastq.dir,
                       name=name,
                       barcode.type="MULTIseq",
                       assay="RNA",
                       filter.cells=valid.cells)
saveRDS(read_table, file.path(outdir, 'read_table.rds'))

# use valid oligos
tag_mtx <- alignTags(read_table, valid.bar.ref)
saveRDS(tag_mtx, file.path(outdir, 'tag_mtx.rds'))

# demux
res <- demultiplexTags(tag_mtx,
                       plot.path=outdir,
                       plot.name=paste0(name, '_deMULTIplex_'),
                       plot.diagnostics=F)

saveRDS(res, file.path(outdir, 'demux.rds'))

out.plots <- paste0(outdir,'/',name,'_deMULTIPLEX_plots.pdf')
pdf(out.plots)

tagHist(tag_mtx = tag_mtx, 
        minUMI = 10,
        plotnUMI = T)


tagCallHeatmap(tag_mtx=tag_mtx,
               calls=res$final_assign)

dev.off()

