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
    /tmp/pig_atlas_counts.csv
.venv/bin/python -m pipeline.build_pig_atlas_cell_types \
    --input /tmp/pig_atlas_counts.csv
.venv/bin/python -m pipeline.main
```

The R script sums Seurat's raw RNA counts over the cells of each atlas
tissue/cell-type pair. The Python step maps atlas tissues to the OrthoHPI
vocabulary, sums the counts of samples that share a tissue label and cell type,
and scales each pool to a million, as HPA does for its single-cell nCPM, so the
same `nCPM > 1` cutoff applies to both. It then maps the atlas gene symbols and
Ensembl IDs to pig STRING IDs and writes `data/pig_atlas_cell_types.parquet`.
`pipeline.main` appends the rows to the human HPA annotations before its left
join with tissue-filtered host proteins.

The ten brain samples (`Brain`, `Cerebellum`, the four lobes, `Hypothalamus`,
`OVoLT`, `AreaPostrema`, `SubfornicalOrgan`) map to `brain`, `Retina` to `eye`,
`PBMC` to `blood`, and `Intestine` and the six neonatal-ileum time points to
`intestine`; other tissues are kept when their lowercased name is already a
label. Samples are pooled at the count level, as for mouse and for the atlas's
own tissues that mix 10x and DNBelab cells, so a cell type's nCPM weights each
sample by its counts. The adult `Intestine` sample shares no cell-type labels
with the neonatal ileum, so the two are never pooled together.
