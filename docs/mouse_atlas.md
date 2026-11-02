# Tabula Muris Senis mouse cell atlas

Mouse cell-type annotation is derived from the Tabula Muris Senis raw droplet
dataset, downloaded from the atlas authors' Figshare record:
[Processed files (to use with Scanpy), v3](https://doi.org/10.6084/m9.figshare.8273102.v3).
The default inputs are the 4.06 GB
`tabula-muris-senis-droplet-official-raw-obj.h5ad` file and the 2.37 GB FACS
file from the same record.

The H5ADs are not in the repository. Fetch them from the URLs under
`urls.tabula_muris_senis_droplet_url` and `urls.tabula_muris_senis_facs_url` in
`config.yml`, then run the builder from the repository root:

```
mkdir -p data/downloads/tabula_muris_senis && cd data/downloads/tabula_muris_senis
curl -L -C - -o tabula-muris-senis-droplet-official-raw-obj.h5ad \
    https://ndownloader.figshare.com/files/23938934
curl -L -C - -o tabula-muris-senis-facs-official-raw-obj.h5ad \
    https://ndownloader.figshare.com/files/23939711
cd -
.venv/bin/python -m pipeline.build_mouse_atlas_cell_types
.venv/bin/python -m pipeline.main
```

The builder selects only `3m` mice, pools the raw counts of each
`(tissue, Cell Ontology cell type)` and scales every pool to a million, as HPA
does for its single-cell nCPM, so the same `nCPM > 1` cutoff applies to both.
It maps genes to mouse STRING IDs and writes
`data/mouse_atlas_cell_types.parquet`; `pipeline.main` then appends the rows to
the human HPA annotations in `data/tissues_cell_types.parquet`.

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

`Large_Intestine` is written out under both `colon` and the coarser `intestine`,
the two labels TISSUES records the large intestine under and the only mouse cell
types either can have; its cells are largely colonic epithelium. `small intestine`
and `rectum` have no atlas tissue of their own and stay without cell types rather
than borrowing those rows.

Tabula Muris Senis does not provide equivalent standalone data for every
mouse-relevant lifecycle tissue. In particular, it does not provide blood,
stomach, mesenteric artery, nose, or mouth -- neither the droplet nor the FACS
dataset has a blood tissue at all. Marrow remains `bone marrow`; it is not
silently treated as blood. The app shows cell-type resolution only where the
atlas has a retained, mapped tissue.
