# Tabula Muris Senis mouse cell atlas

Mouse cell-type annotation is derived from the Tabula Muris Senis raw droplet
dataset, downloaded from the atlas authors' Figshare record:
[Processed files (to use with Scanpy), v3](https://doi.org/10.6084/m9.figshare.8273102.v3).
The default inputs are the 4.06 GB
`tabula-muris-senis-droplet-official-raw-obj.h5ad` file and the 2.37 GB FACS
file from the same record.

Run the builder from the repository root:

```
.venv/bin/python -m pipeline.build_mouse_atlas_cell_types
.venv/bin/python -m pipeline.refresh_mouse_cell_type_annotations
```

The builder downloads the H5AD when it is not already present, selects only
`3m` mice, normalizes each cell to 10,000 counts, applies `log1p`, and averages
expression for each `(tissue, Cell Ontology cell type, gene)` combination. It
maps genes to mouse STRING IDs and writes `data/mouse_atlas_cell_types.parquet`;
the refresh command adds matching mouse rows to
`data/tissues_cell_types.parquet`, the annotation artifact read by the app. It
does not rerun the full prediction pipeline or download unrelated references.
The atlas download is streamed to a temporary file and is only used after its
HDF5 signature is validated, so an interrupted download will be discarded and
retried safely.

The droplet and Smart-seq2/FACS datasets use different count technologies, so
their values are never pooled. Droplet data is preferred for every tissue it
covers; FACS data is used only as a fallback for tissues missing from the
3-month droplet data, including brain and gut tissues. Pass `--droplet-only`
to omit that fallback, or `--input` / `--facs-input` to use explicitly
downloaded compatible H5AD files.

## Tissue mapping

The builder retains only atlas tissues with a direct, documented mapping to the
OrthoHPI tissue vocabulary: bladder, brain (myeloid and non-myeloid), heart,
kidney, large intestine, limb muscle, liver, lung, marrow, pancreas, skin, and
spleen. The droplet atlas's `tissue_free_annotation` is used where present, so
aorta cells in the combined `Heart_and_Aorta` source category are not labelled
as heart.

Tabula Muris Senis does not provide equivalent standalone data for every
mouse-relevant lifecycle tissue. In particular, it does not provide stomach,
mesenteric artery, nose, or mouth. Marrow remains `bone marrow`; it is not
silently treated as blood. The app shows cell-type resolution only where the
atlas has a retained, mapped tissue.
