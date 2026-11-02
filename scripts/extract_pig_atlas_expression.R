#!/usr/bin/env Rscript

# Sum Seurat's raw RNA counts over the cells of each atlas (tissue, cell type). The
# Python step maps tissues, pools samples and scales each pool to a million.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/extract_pig_atlas_expression.R INPUT.rds OUTPUT.csv")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

atlas <- readRDS(args[[1]])
metadata <- atlas[[]]
counts <- GetAssayData(atlas, assay = "RNA", layer = "counts")
rm(atlas)

groups <- unique(metadata[, c("Tissue", "Celltype")])
group <- match(paste(metadata$Tissue, metadata$Celltype, sep = "\r"),
               paste(groups$Tissue, groups$Celltype, sep = "\r"))
membership <- sparseMatrix(i = seq_along(group), j = group, x = 1,
                           dims = c(length(group), nrow(groups)))
pooled <- as(counts %*% membership, "TsparseMatrix")

write.csv(data.frame(Gene = rownames(counts)[pooled@i + 1],
                     Tissue = groups$Tissue[pooled@j + 1],
                     `Cell type` = groups$Celltype[pooled@j + 1],
                     Counts = pooled@x, check.names = FALSE),
          args[[2]], row.names = FALSE)
