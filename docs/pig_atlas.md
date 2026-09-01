# Pig cell atlas (`pig_atlas_20221014.rds`)

A Seurat v4.0.2 object, project `"Pig"`, written by R 4.0.3. 14.1 GB
uncompressed — loading it whole is the real obstacle, not R.

Downloaded from <https://dreamapp.biomed.au.dk/pigatlas/pig_atlas_20221014.rds>
on 2026-08-17.

| | |
| --- | --- |
| Genes × cells | 20,922 × 295,417 |
| `counts` | `dgCMatrix`, 384,190,890 non-zeros |
| `data` | log-normalised, same sparsity |
| `scale.data` | dense 2,000 × 295,417 doubles ≈ 4.7 GB — the memory hog |
| `active.ident` | factor, 26 tissues (Liver, Heart, Kidney, Spleen, Adipose-S, …) |
| `meta.data` | 15 columns, incl. cell type and cell-cycle phase |
| Reductions | `pca` (20 PCs), `umap`, `tsne` |
| Provenance | `NormalizeData` → `FindVariableFeatures` (vst, 2000) → `ScaleData` → `RunPCA` → `RunUMAP` |

Gene names are mixed symbols and Ensembl IDs (`TBP`, `PSMB1`,
`ENSSSCG00000036983`); cell barcodes are 10x.

## Reading it without R

`scripts/read_rds.py` is a pure-Python RDS reader — no R, no dependencies. It
streams the gzip and skips bulk numeric payloads, so it peaks at a few hundred
MB rather than 14 GB.

```
.venv/bin/python scripts/read_rds.py pig_atlas_20221014.rds --depth 3
.venv/bin/python scripts/read_rds.py pig_atlas_20221014.rds --dump outdir
```

A full pass takes ~5 min. `--dump` produces:

| file | |
| --- | --- |
| `meta_data.csv` | 66 MB, 295,417 × 15 |
| `reduction_umap.csv` / `_tsne.csv` | 11 MB each |
| `genes.txt` | 20,922 genes |
| `active_ident.csv` | tissue per cell |

Two caveats:

- **PCA embeddings are not dumped.** At 295,417 × 20 they exceed the retention
  threshold, so the script skips them rather than write a truncated matrix.
  Re-run with `--keep 6000000` to get them.
- **The expression matrices are not extracted** — only annotations and
  embeddings. Pulling `counts` into a scipy sparse matrix is a further step;
  worth doing for per-tissue expression on the host-protein side of OrthoHPI.

The dumped files live in the scratchpad rather than `data/`, since `data/` is
tracked by git and these are ~100 MB.

## Metadata columns (15)

### Numeric (11)

| column | min | median | max |
| --- | --- | --- | --- |
| `nCount_RNA` | 226 | 2,367 | 52,233 |
| `nFeature_RNA` | 201 | 1,065 | 7,541 |
| `percent.mt` | 0 | 5.16 | 29.99 |
| `percent.rib` | 0 | 8.91 | 40.20 |
| `percent.exc` | 0 | 4.90 | 63.36 |
| `percent.cyt` | 7.46 | 44.36 | 81.05 |
| `percent.mem` | 0.69 | 26.09 | 60.17 |
| `percent.pcg` | 10.84 | 66.49 | 100.00 |
| `percent.tf` | 0.02 | 3.92 | 26.89 |
| `S.Score` | −9.44 | −0.03 | 6.61 |
| `G2M.Score` | −9.06 | −0.07 | 11.28 |

The `percent.*` columns are pre-computed gene-set fractions — mitochondrial,
ribosomal, excreted/secreted, cytoplasmic, membrane, protein-coding,
transcription factors. `percent.mem` and `percent.exc` are the ones relevant to
host-parasite work, since surface and secreted proteins are where interactions
happen.

### Categorical (4)

- `Phase` — G1 177,203 · S 74,181 · G2M 44,033
- `Celltype` — 71 values
- `Tissue` — 26 values (identical to `active.ident`)
- `Platform` — 10x Genomics 238,633 · DNBelab C4 56,784

## Cell types (71)

| # | type | cells | | # | type | cells |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | T cells | 26,624 | | 37 | Plasma cells | 2,040 |
| 2 | Oligodendrocytes | 21,262 | | 38 | Monocytes | 1,590 |
| 3 | Macrophages | 16,207 | | 39 | Retinal ganglion cells | 1,527 |
| 4 | Endothelial cells | 16,112 | | 40 | EEC | 1,521 |
| 5 | Microglia | 14,971 | | 41 | EI | 1,519 |
| 6 | B cells | 13,307 | | 42 | Neutrophils | 1,342 |
| 7 | Goblet | 12,843 | | 43 | Bipolar cells | 1,167 |
| 8 | Excitatory neurons | 11,072 | | 44 | Dendritic cells | 1,136 |
| 9 | Hepatocytes | 10,908 | | 45 | Cycling cells | 963 |
| 10 | EPL | 10,631 | | 46 | Others | 947 |
| 11 | Granule cells | 8,935 | | 47 | CD8+ Naive T cells | 651 |
| 12 | Mesenchymal cells | 8,290 | | 48 | Ciliated cells | 640 |
| 13 | Rod photoreceptor cells | 8,012 | | 49 | Naive B cells | 614 |
| 14 | EPE | 7,730 | | 50 | Memory B cells | 560 |
| 15 | EP | 6,167 | | 51 | Muller glia | 559 |
| 16 | TA | 6,153 | | 52 | Transient amplifying cells | 439 |
| 17 | EM | 6,084 | | 53 | Beige adipocytes | 390 |
| 18 | Kupffer cells | 5,859 | | 54 | Cone photoreceptor cells | 390 |
| 19 | Tuft | 5,779 | | 55 | Mast cells | 369 |
| 20 | Smooth muscle cells | 5,645 | | 56 | Purkinje cells | 322 |
| 21 | Stem | 5,113 | | 57 | Hepatic stellate cells | 257 |
| 22 | AT2 | 4,659 | | 58 | Loop of Henle's cells | 229 |
| 23 | Inhibitory neurons | 4,641 | | 59 | Immune cells | 225 |
| 24 | Astrocytes | 4,218 | | 60 | Ureteric epithelial cells | 182 |
| 25 | TA-G1 | 3,689 | | 61 | Distal convoluted tubule cells | 166 |
| 26 | TA-G2 | 3,605 | | 62 | AT1/AT2 | 138 |
| 27 | Oligodendrocyte progenitor cells | 3,575 | | 63 | Clara cells | 122 |
| 28 | NK cells | 3,481 | | 64 | Neural progenitor cells | 105 |
| 29 | Proximal tubule cells | 3,018 | | 65 | Megakaryocytes | 97 |
| 30 | Fibroblasts | 2,724 | | 66 | Basophils | 85 |
| 31 | CD8+ Cytotoxic T cells | 2,429 | | 67 | Collecting duct cells | 84 |
| 32 | Erythroid cells | 2,303 | | 68 | Podocytes | 75 |
| 33 | Cardiomyocytes | 2,300 | | 69 | Newly formed oligodendrocytes | 57 |
| 34 | AT1 | 2,294 | | 70 | Enterocytes | 27 |
| 35 | CD4+ Naive T cells | 2,170 | | 71 | Cholangiocytes | 15 |
| 36 | Paneth | 2,057 | | | | |

The neonatal-ileum samples use the labels `EPL`, `EPE`, `EP`, `EM`, `EI`,
`TA`, `TA-G1`, `TA-G2`, `EEC`, `Stem`, `Paneth`, `Tuft`, and `Goblet`. The RDS
metadata does not define the abbreviations or their biological relationships,
so they are retained exactly as supplied and must not be expanded or merged
without a cited source from the atlas authors.

Three labelling quirks to be aware of before joining on these:

- `T cells` coexists with the more specific CD4+/CD8+ subsets, so they are not
  mutually exclusive.
- `CD4+ Naive  T cells` has a double space.
- `Others` (947) and `Immune cells` (225) are catch-alls.
- `Transient amplifying cells` occurs in `Intestine`, while `TA` occurs in the
  neonatal-ileum samples. The metadata does not state whether they represent
  the same population.

## Cell types per tissue

| tissue | cells | dominant types |
| --- | --- | --- |
| Liver | 32,888 | Hepatocytes (10,908), T cells (6,559), Kupffer (5,859), Macrophages (4,148) |
| Lung | 25,725 | Macrophages (6,429), T cells (4,660), AT2 (4,659), AT1 (2,294) |
| Spleen | 21,791 | T cells (8,726), Macrophages (4,078), B cells (2,603), Plasma (1,411) |
| Adipose-V | 16,417 | Endothelial (7,020), Mesenchymal (3,324), Smooth muscle (2,402) |
| Brain | 15,794 | Microglia (9,086), Excitatory neurons (3,479), T cells (1,045) |
| Adipose-S | 12,814 | Mesenchymal (4,966), Endothelial (3,836), Smooth muscle (2,733) |
| Retina | 12,164 | Rod photoreceptors (8,012), Retinal ganglion (1,527), Bipolar (1,167) |
| Intestine | 10,518 | B cells (8,920), T cells (860), Transient amplifying (439) |
| Cerebellum | 10,225 | Granule cells (8,935), Astrocytes (443), Purkinje (322) |
| PBMC | 10,090 | CD8+ Cytotoxic T (2,429), CD4+ Naive T (2,170), T cells (1,920) |
| FrontalLobe | 8,812 | Oligodendrocytes (4,916), Microglia (1,887), Others (895) |
| OVoLT | 7,739 | Excitatory (2,094), Inhibitory (1,923), Astrocytes (1,635) |
| Hypothalamus | 7,302 | Oligodendrocytes (3,204), Microglia (1,558), OPC (1,223) |
| OccipitalLobe | 6,829 | Excitatory neurons (2,807), Oligodendrocytes (2,633) |
| ParietalLobe | 6,162 | Oligodendrocytes (3,373), Excitatory neurons (1,360) |
| AreaPostrema | 5,109 | Oligodendrocytes (3,597), Astrocytes (501), Inhibitory (434) |
| Kidney | 4,255 | Proximal tubule (3,018), Loop of Henle (229), Endothelial (212) |
| Heart | 3,220 | Cardiomyocytes (2,300), Fibroblasts (287), Macrophages (225) |
| TemporalLobe | 3,145 | Oligodendrocytes (1,632), Excitatory neurons (489) |
| SubfornicalOrgan | 1,527 | Inhibitory (512), Oligodendrocytes (364), Excitatory (281) |
| Ileum D0→D21 | 72,891 | shifts Tuft/Goblet → EPL → TA/Stem → EPE/TA-G1 |

Intestine is 85% B cells while the neonatal ileum samples are almost entirely
epithelial — those two sample different compartments despite both being gut, so
don't pool them.

## Using it for pig cell-type resolution

The app reads cell-type annotation from `data/tissues_cell_types.parquet`. Build
the pig contribution once before rebuilding that pipeline artifact:

```
Rscript scripts/extract_pig_atlas_expression.R data/pig_atlas_20221014.rds \
    /tmp/pig_atlas_expression.csv
.venv/bin/python -m pipeline.build_pig_atlas_cell_types \
    --input /tmp/pig_atlas_expression.csv
.venv/bin/python -m pipeline.main
```

The R script reads Seurat's normalized sparse RNA assay and writes the mean
normalized expression of every gene in each atlas tissue/cell-type pair. The
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
