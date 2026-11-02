# Pig cell atlas (`pig_atlas_20221014.rds`)

A Seurat v4.0.2 object (20,922 genes × 295,417 cells, 26 tissues, 71 cell types)
from <https://dreamapp.biomed.au.dk/pigatlas/pig_atlas_20221014.rds>,
downloaded on 2026-08-17.

The atlas's own cell-type labels are kept exactly as supplied. The
neonatal-ileum samples use unexplained abbreviations (`EPL`, `EPE`, `EP`, `EM`,
`EI`, `TA`, `TA-G1`, `TA-G2`, `EEC`, `Stem`, `Paneth`, `Tuft`, `Goblet`), `T cells`
coexists with the more specific CD4+/CD8+ subsets, and `Others` and
`Immune cells` are catch-alls. `Intestine` is 85% B cells while the neonatal
ileum is almost entirely epithelial.

## Using it for pig cell-type resolution

The app reads cell-type annotation from `data/tissues_cell_types.parquet`. The
pipeline's own input is `data/pig_atlas_cell_types.parquet`, which is committed;
the 2.1 GB `.rds` is needed only to regenerate it, and `pipeline.main` never
reads or downloads it. The atlas is a fixed 2022 snapshot, so a regeneration is
only called for when its tissue mapping changes or a newer atlas is adopted.

The `.rds` is not in the repository. Fetch it from the URL under
`urls.pig_atlas_url` in `config.yml` with a resumable transfer — the download is
long enough that an interrupted one is worth continuing rather than restarting:

```
curl -C - -o data/pig_atlas_20221014.rds \
    https://dreamapp.biomed.au.dk/pigatlas/pig_atlas_20221014.rds
```

Extracting from it needs R with Seurat and Matrix, which the `.venv` does not
provide. Build the pig contribution before rebuilding the pipeline artifact:

```
Rscript scripts/extract_pig_atlas_expression.R data/pig_atlas_20221014.rds \
    /tmp/pig_atlas_expression.csv
.venv/bin/python -m pipeline.build_pig_atlas_cell_types \
    --input /tmp/pig_atlas_expression.csv
.venv/bin/python -m pipeline.main
```

The R script reads Seurat's raw RNA counts, pools them over the cells of each
atlas tissue/cell-type pair and scales each pool to a million, as HPA does for
its single-cell nCPM, so the same `nCPM > 1` cutoff applies to both. The
Python step maps the atlas gene symbols and Ensembl IDs to pig STRING IDs,
normalizes atlas tissue names to the OrthoHPI vocabulary, and writes
`data/pig_atlas_cell_types.parquet`. `pipeline.main` then appends the rows to
the human HPA annotations before its left join with tissue-filtered host
proteins.

Only `Intestine`, each neonatal-ileum time point, and `PBMC` are normalized to
the pipeline's `intestine` and `blood` tissue labels. The ileum time points
remain separate during aggregation and are averaged only after they map to the
same `(STRING protein, tissue, cell type)` annotation; they are not pooled at
the cell level.
