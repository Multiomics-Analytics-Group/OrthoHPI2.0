#!/usr/bin/env Rscript

# Aggregate Seurat's normalized sparse RNA matrix by its annotated tissue and
# cell type. This avoids materializing a genes-by-cells dense matrix.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript scripts/extract_pig_atlas_expression.R INPUT.rds OUTPUT.csv")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

atlas <- readRDS(args[[1]])
metadata <- atlas[[]][, c("Tissue", "Celltype"), drop = FALSE]
# SeuratObject 5 replaced the Seurat 4 `slot` argument with `layer`.
expression <- if (utils::compareVersion(as.character(utils::packageVersion("SeuratObject")),
                                        "5.0.0") >= 0) {
  GetAssayData(atlas, assay = "RNA", layer = "data")
} else {
  GetAssayData(atlas, assay = "RNA", slot = "data")
}
groups <- interaction(metadata$Tissue, metadata$Celltype, drop = TRUE, lex.order = TRUE)

output <- args[[2]]
write.csv(data.frame(Gene = character(), Tissue = character(), `Cell type` = character(),
                     nTPM = numeric(), check.names = FALSE),
          output, row.names = FALSE)

for (group in levels(groups)) {
  columns <- which(groups == group)
  if (!length(columns)) next

  parts <- strsplit(group, ".", fixed = TRUE)[[1]]
  means <- Matrix::rowMeans(expression[, columns, drop = FALSE])
  rows <- data.frame(Gene = rownames(expression), Tissue = parts[[1]],
                     `Cell type` = parts[[2]], nTPM = means, check.names = FALSE)
  rows <- rows[rows$nTPM > 0, ]
  write.table(rows, output, sep = ",", row.names = FALSE,
              col.names = FALSE, append = TRUE, quote = TRUE)
}
