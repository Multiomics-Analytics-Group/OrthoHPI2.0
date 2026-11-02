# Prediction pipeline

The pipeline in `pipeline/` produces `data/predictions.parquet`, which the app
reads. The other files in `data/` are listed under [Data files](#data-files).

It is configured in `config.yml`:

- `urls`: data sources (STRING v12.0, EggNOG 6, ...)
- `hosts`: human (9606), rat (10116), mouse (10090) and pig (9823)
- `parasites`: 46 species, each with a STRING taxid, a label and colour, the BTO
  codes of the tissues it encounters during its life cycle, and the taxids of
  its hosts

## Stages

### 0. EggNOG members file (`pipeline/prepare_eggnog_members.py`)

Run once before `main.py`. The script streams `e6.og2seqs_and_species.tsv` (~10 GB,
all levels), keeps the Eukaryota (2759) groups and writes them to
`data/downloads/2759_members.tsv.gz` in the EggNOG 5 column layout.

### 1. Setup (`main.setup`)

Downloads every URL in `config.yml` except `eggNOG_members_url` (stage 0) and
the per-species STRING templates (`string_protein_url`, `string_ppi_url`,
`string_go_url`, `string_alias_url`, `string_sequences_url`), which contain a
`TAXID` placeholder and are fetched later. It also runs
`go.get_gene_ontology`, which downloads the GO ontology and STRING's GO
annotations and writes `gos.parquet` and `go_ontology.parquet`.

### 2. Proteins (`main.get_proteins`)

Downloads `protein.info.v12.0/{taxid}...txt.gz` from STRING for every host and
parasite and builds `{taxid: {protein_id: protein_name}}`. The filters below
narrow this set.

### 3. Parasite secretome (`filters.get_secretome_predictions`)

Only secreted or surface-exposed parasite proteins can reach the host, so the
rest are dropped before the homology transfer. The localizations come from
DeepLoc 2.1, run separately (see [DeepLoc secretome
prediction](#deeploc-secretome-prediction)), which leaves the proteins that
passed in `data/secretome/{taxid}.fasta`. This step keeps the parasite proteins
whose IDs are in that file.

### 4. Host tissues (`filters.apply_tissue_filter`, `get_tissues`)

For each host, downloads the TISSUES scores from its `tissues_url`
and keeps a host protein if its score is at least `tissue_cutoff` in one of the
BTO tissues listed for the host's parasites. The cutoff is 1.0 for every host.

Matches are collected in `{protein: [tissue, ...]}` across hosts and used again
in step 6.

The pool left after this step is per host. `filters.parasite_tissue_proteins`
splits it per parasite (`{parasite taxid: host proteins in its tissues}`), which
goes into the per-parasite targets (end of step 5). Proteins of a host without a
`tissues_url` pass for every parasite.

### 5. Host localization (`filters.apply_deeploc_filter`)

Keeps the host proteins a parasite can reach, using DeepLoc 2 (Accurate)
predictions from
`data/deeploc/output_accurate/deeploc_output_accurate/<taxid>/results_*.csv`.
The classes that count depend on the parasite's `niche` in `config.yml`
(`main.DEEPLOC_NICHE_CLASSES`):

- extracellular: `P(Extracellular) > 0.617` or `P(Cell membrane) > 0.565`
- intracellular: the same two, or `P(Cytoplasm) > 0.476` or
  `P(Nucleus) > 0.501`.

The niche is a property of the parasite but the pool is per host. The filter
narrows `valid_proteins[taxid]` to the union of the niches of the host's
parasites, and returns `{niche: host proteins}`. A parasite without a valid
`niche` is treated as extracellular.

`filters.parasite_targets` then combines the per-parasite parts of steps 4 and 5
into `{parasite taxid: host proteins it can be predicted to bind}`: proteins of
the hosts it infects, reached from its niche and expressed in a tissue it
infects. `homology.get_links` keeps a pair only if the host protein is in its
parasite's set.

`main.save_eligible_proteins` writes one `reachable_<niche>` column per niche
to `eligible_proteins.parquet`. The app uses it as the enrichment background,
so each parasite is tested against the host proteins it could reach.

After steps 3–5, `utils.merge_dict_of_dicts` flattens the per-taxid dicts into
one `{protein_id: name}` dict.

### 6. Tissue and cell-type annotation (`main.get_tissue_cell_type_annotation`)

Turns the tissue matches from step 4 into a `(Gene, Tissue)` table and
left-joins single-cell data onto it
(`cell_type_annotations.parse_cell_type_data`): the Human Protein Atlas for
human, and the Pig Cell Atlas and Tabula Muris Senis (mouse). Rat
proteins have no cell types.

The HPA data is the v25 single-cell cluster file, in nCPM, reduced to the
highest cluster per gene, tissue and cell type. HPA tissue names differ from
the labels in `config['tissues']`, so `cell_type_annotations.HPA_TISSUE_LABELS`
maps them first: `heart muscle` → `heart`, `skeletal muscle` → `muscle`,
`brain neurons` and `brain non-neurons` → `brain`, `retina` → `eye`. `colon`,
`rectum` and `small intestine` are also written under `intestine`, since
TISSUES uses both. When several HPA tissues or genes map to the same protein,
tissue and cell type, the highest nCPM is kept. Tissues HPA has no single-cell
data for (`mouth`, `nose`, `spinal cord`, `vagina`, `macrophage`, ...) get no
cell types.

The pig and mouse atlases are preprocessed outside the pipeline into
`data/pig_atlas_cell_types.parquet` and `data/mouse_atlas_cell_types.parquet`,
which `pipeline.main` reads but never rebuilds. Both are put on HPA's scale:
the raw counts of each tissue and cell type are pooled and scaled to a million,
so the same `nCPM > 1` cutoff applies to all three sources, and when several
atlas genes map to one protein the highest nCPM is kept.

- **Pig** ([pig_atlas.md](pig_atlas.md)): the Pig Cell Atlas Seurat object.
  `scripts/extract_pig_atlas_expression.R` sums the raw counts per atlas tissue
  and cell type, and `pipeline/build_pig_atlas_cell_types.py` maps the atlas
  tissues to `config['tissues']` labels, pools samples that share a label (the
  ten brain regions → `brain`, `Retina` → `eye`, `PBMC` → `blood`, the adult
  intestine and neonatal ileum → `intestine`) and scales each pool.
- **Mouse** ([mouse_atlas.md](mouse_atlas.md)): Tabula Muris Senis, 3-month
  mice, Cell Ontology cell types. `pipeline/build_mouse_atlas_cell_types.py`
  uses the droplet data, with FACS data only for tissues the droplet data
  lacks; the two technologies are never pooled. `Large_Intestine` is written
  under both `intestine` and `colon`.

The result is written to `data/tissues_cell_types.parquet`. It is used by the
app only and does not filter anything.

### 7. Homology transfer (`homology.get_eggnog_groups`, `homology.get_links`)

Interactions are transferred from STRING through orthology groups:

1. `get_eggnog_groups` reads `2759_members.tsv.gz` and keeps the groups that
   contain at least one protein left after step 5.
2. `get_links` reads STRING's `COG.links.detailed`, which scores evidence
   between orthology groups. A link is kept if both groups were kept in 1,
   `experimental_evidence >= 0.7` or `databases_evidence >= 0.7`, and it
   connects a parasite protein to a host protein in that parasite's targets
   (links within one species are skipped).

Every such protein pair goes into the predictions table, with the parasite
protein in the `source_*` columns and the host protein in `target_*`.
So if a parasite protein's group interacts with a host protein's group in
STRING, the two proteins are predicted to interact.

`COG.links.detailed` only names groups as `COG`, `KOG` or `NOG` ids, while most
level-2759 groups in EggNOG 6 have hash ids (e.g. `2QPHS`). In practice only
the ~4.8k KOG groups produce predictions.

### 8. UniProt accessions (`main.annotate_predictions`, `utils.annotate_alias_id`)

Adds `source_uniprot` and `target_uniprot` from STRING's alias files and writes
`data/predictions.parquet`.

The alias files differ between species, so `sources` is a preference list and
`parse_string_aliases` takes the first source that has an alias. Parasites use
`Uniprot`. Hosts use `Ensembl_HGNC_uniprot_ids` (one canonical accession per
protein, human only) and fall back to `UniProt_AC` for rat, mouse and pig,
which otherwise get no accession and so no structures in the app. Doesn't use
`Ensembl_UniProt` because it has gene names mixed in with the accessions.

## DeepLoc secretome prediction

DeepLoc produces the parasite secretome for step 3. It is run separately from
`pipeline/main.py` because it needs a GPU and, on the cluster, a different
Python environment. The scripts are in `deeploc/`, and they have to be run
first so that `data/secretome/{taxid}.fasta` exists.

There are two ways to run it. Both write `data/secretome/{taxid}.fasta`.

Locally, with `deeploc/run_deeploc.py`: downloads each parasite proteome from
STRING, runs `deeploc2 -f ... -d mps` (Apple GPU), filters the results and
writes the FASTA.

Three scripts:

1. `deeploc/prepare_deeploc_input.py` downloads the proteomes and writes
    FASTAs to `data/deeploc/input/`.
2. `deeploc/generate_hpc_jobs.py` writes one `bsub` job per parasite to `jobs/`,
   each running `deeploc2 -d cuda` on a GPU node. Submit them with
   `for f in jobs/*.sh; do bsub < $f; done`.
3. `deeploc/build_secretome_fastas.py` reads the `results_*.csv` files once
   they are copied back and writes the filtered FASTAs.

A protein is kept if it is extracellular, or in the cell membrane unless the
parasite has `multicellular: true` in `config.yml` (30 of them do, and keep
only extracellular proteins). `build_secretome_fastas.py` (`filter_by_score`)
reads this from the probabilities, `Extracellular > 0.617` and
`Cell membrane > 0.565`; `run_deeploc.py` (`filter_by_deeploc`) from DeepLoc's
predicted labels.

## Host orthologs (`scripts/build_host_orthologs.py`)

The app's cross-host page compares one parasite across its hosts. When an
interaction is predicted in one host and not another, either the other host has
no protein in that orthology group, or it has one that the filters removed
before step 7. `predictions.parquet` can't tell these apart.

This script reads `data/downloads/2759_members.tsv.gz` and writes, for each
host orthology group in `predictions.parquet`, the proteins each host has in
it:

```
$ python scripts/build_host_orthologs.py
```

The output is `data/host_orthologs.parquet` (`group`, `taxid`, `n_proteins`,
`proteins`), a few thousand rows covering the ~1,000 groups in the predictions.
Re-run it after `pipeline/main.py` when hosts or parasites change. If the file
is missing, the page leaves that section out.

## Data files

The parquet files in `data/`. `pipeline/main.py` writes the first four; the
others are built by their own script, to be re-run after `main.py` when hosts or
parasites change.

| File | Written by | One row per | Columns | Read by |
|---|---|---|---|---|
| `predictions.parquet` | `pipeline/main.py` (steps 7–8) | predicted host–parasite protein pair | `taxid1`/`source*` (parasite), `taxid2`/`target*` (host), `experimental_evidence_score`, `databases_evidence_score`, `weight`, `group1`, `group2`, `source_uniprot`, `target_uniprot` | the app and most scripts |
| `eligible_proteins.parquet` | `pipeline/main.py` (step 5) | host or parasite protein left after steps 3–5 | `protein`, `name`, `taxid`, `reachable_<niche>` | app enrichment background, filter funnel, AlphaFold job scripts |
| `tissues_cell_types.parquet` | `pipeline/main.py` (step 6) | host protein, tissue and cell type (one row per tissue if no cell types) | `Gene`, `Tissue`, `Gene name`, `Cluster`, `Cell type`, `Read count`, `nCPM` | app tissue and cell-type views |
| `gos.parquet` | `pipeline/go.py` (step 1) | protein and GO term, all species | `#string_protein_id`, `term`, `description`, `taxid` | app GO enrichment |
| `go_ontology.parquet` | `pipeline/go.py` (step 1) | GO parent–child relation | `parent`, `child` | app GO enrichment |
| `host_orthologs.parquet` | `scripts/build_host_orthologs.py` | host orthology group and host | `group`, `taxid`, `n_proteins`, `proteins` | app cross-host page |
| `host_tissues.parquet` | `scripts/build_host_tissues.py` | host protein and tissue, after step 4 but before step 5 | `protein`, `tissue` | app cross-host page, to tell "not expressed" from "out of reach" |
| `deeploc_localisations.parquet` | `pipeline/build_deeploc_localisations.py` | protein in the predictions | `protein`, `localizations`, `signals`, `membrane_types`, one probability column per class | app localisation of interactors |
| `protein_annotations.parquet` | `pipeline/build_protein_annotations.py` | protein in the predictions | `protein`, `description` | app network tooltips |
| `protein_sequences.parquet` | `pipeline/build_protein_sequences.py` | protein in the predictions | `protein`, `sequence` | app FASTA downloads |
| `figure_tissues.parquet` | `scripts/build_figure_tissues.py` | host protein and organ on the body figures | `Gene`, `Organ` | app body figure |
| `pig_atlas_cell_types.parquet` | `pipeline/build_pig_atlas_cell_types.py` | gene, tissue and cell type | `Gene`, `Tissue`, `Cell type`, `nCPM` | step 6 |
| `mouse_atlas_cell_types.parquet` | `pipeline/build_mouse_atlas_cell_types.py` | gene, tissue and cell type | `Gene`, `Tissue`, `Cell type`, `nCPM` | step 6 |
