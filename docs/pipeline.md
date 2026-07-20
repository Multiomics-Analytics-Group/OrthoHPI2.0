# Prediction pipeline

This document describes what happens when you run `python -m pipeline.main`
(the `main.py` invoked from the README as `python main.py`) — how the
host-parasite PPI predictions in `data/annotated_predictions.parquet` get
built. 

Everything is driven by `config.yml`: `urls` (data sources), `hosts` (only
human, 9606, currently), and `parasites` (35 species, each with a STRING
taxid, a display label/color, and the list of BTO tissue codes relevant to
that parasite's life cycle).

Data sources are STRING v12.0 and EggNOG 6.

## Stages

### 0. Build the EggNOG members file (`pipeline/prepare_eggnog_members.py`)

Run once, before `main.py`. EggNOG 6 no longer publishes per-tax-level
downloads, so the Eukaryota (2759) groups have to be extracted from
`e6.og2seqs_and_species.tsv` (~10 GB, every tax level). The script streams
and filters that file and writes `data/downloads/2759_members.tsv.gz` in
EggNOG 5's column layout, so the rest of the pipeline is unchanged. `main.py`
stops with an error if the file is missing.

### 1. Setup — download reference files (`main.setup`)

Downloads every URL in `config.yml` except the per-species templates
(`string_protein_url`, `string_ppi_url`, `string_go_url`,
`string_alias_url`, `string_sequences_url` — these contain a `TAXID`
placeholder and are fetched later, once per species) and
`eggNOG_members_url` (handled by stage 0). Also triggers
`go.get_gene_ontology`, which downloads the GO OBO ontology and STRING's
per-protein GO annotations and writes `gos.parquet` / `go_ontology.parquet`.

### 2. Get proteins (`main.get_proteins`)

For every host + parasite taxid, downloads
`protein.info.v12.0/{taxid}...txt.gz` from STRING and builds
`{taxid: {ensembl_protein_id: protein_name}}`. This is the set of
candidate proteins that every later filter narrows down.

### 3. Secretome filter — DeepLoc (`filters.get_secretome_predictions`)

Parasite proteins are only biologically able to interact with a host if
they're secreted or exposed on the cell surface — everything else (e.g.
purely cytoplasmic proteins) is discarded before homology transfer
starts. Subcellular localization is predicted by **DeepLoc 2.1**, run
*separately* from this pipeline (see [DeepLoc secretome
prediction](#deeploc-secretome-prediction) below) and expected to have
already produced `data/secretome/{taxid}.fasta` — a FASTA of just the
proteins that passed the localization filter.

`get_secretome_predictions` reads that FASTA per parasite and removes any
protein from `valid_proteins[parasite]` whose ID isn't in it. Hosts are
untouched at this stage.

### 4. Tissue filter (`filters.apply_tissue_filter` → `get_tissues`)

Downloads human tissue expression scores from `hosts_url` (jensenlab
tissues.jensenlab.org) and keeps a host protein only if:
- its expression score in some tissue is `>= cutoff` (2.5 in `main.py`), and
- that tissue is one of the BTO codes listed under the relevant parasite(s)
  in `config.yml` (i.e. a tissue the parasite actually encounters during its
  life cycle — gut, skin, blood, etc.).

Every match is recorded in a `tissues` dict (`{protein: [tissue, ...]}`)
which is returned and later joined with HPA cell-type data in
`main.get_tissue_cell_type_annotation`. Note only the last host processed
in the loop sets `tissues` in the current implementation — this only works
correctly because there is exactly one host (human) configured.

### 5. Compartment filter (`filters.apply_compartment_filter` → `get_compartments`)

Same shape as the tissue filter, but against jensenlab's subcellular
compartment scores (`compartments.jensenlab.org`) instead of tissues. A host
protein passes if its compartment score is `>= cutoff` and the compartment
GO term matches the parasite's configured `compartments` entry (defaulting
to `GO:0005886`, plasma membrane, if unset) — i.e. the protein needs to be
somewhere the parasite could physically reach it.

After steps 3–5, `proteins = utils.merge_dict_of_dicts(proteins)` flattens
the per-taxid dicts into one `{protein_id: name}` dict for the rest of the
pipeline.

### 6. Tissue / cell-type annotation (`main.get_tissue_cell_type_annotation`)

Reshapes the `tissues` dict from step 4 into a long `(Gene, Tissue)` table,
then left-joins it with Human Protein Atlas single-cell RNA data
(`hpa.parse_hpa`) to add cell-type resolution on top of tissue resolution.
Written to `data/tissues_cell_types.parquet` — this is annotation data
consumed by the Streamlit app, not a filter.

### 7. Homology transfer (`homology.get_eggnog_groups`, `homology.get_links`)

The interaction predictions are made via orthology-based transfer
of STRING interactions:

1. `get_eggnog_groups` scans the EggNOG `2759_members.tsv.gz` file (level
   2759 = Eukaryota, built in stage 0) and keeps only groups that contain at
   least one of the filtered proteins from step 5.
2. `get_links` scans STRING's `COG.links.detailed` file, which lists
   experimental/database evidence scores between orthology groups. Note this
   file only names groups as `COG####`/`KOG####`/`NOG####`, while most
   EggNOG groups at level 2759 have hash-style ids (e.g. `2QPHS`) that do not appear
   in it. Only the ~4.8k `KOG` groups can therefore contribute
   predictions.
   A link is kept if:
   - both groups are in `valid_groups`,
   - `experimental_evidence >= 0.7` OR `databases_evidence >= 0.7`,
   - one member is a host protein and the other a parasite protein (same
     taxid pairs are skipped — no intra-species links).

   For every protein pair in the two groups meeting these conditions, an
   inter-species edge is written to `data/predictions.parquet`, with the
   parasite protein always in `source_*` columns and the host protein in
   `target_*` columns.

   In short: if a parasite protein's orthology group is known to interact
   (in STRING) with a human protein's orthology group, that interaction is
   transferred to this specific parasite–host protein pair.

### 8. UniProt alias annotation (`main.py` `__main__`, `utils.annotate_alias_id`)

Adds `source_uniprot` / `target_uniprot` columns by mapping STRING protein
IDs to UniProt accessions via STRING's alias file, for parasites and hosts
respectively. Final output: `data/annotated_predictions.parquet`, which the
Streamlit app reads.

## DeepLoc secretome prediction

DeepLoc predicts subcellular localization from sequence and is what backs
the secretome filter in step 3. It runs as a separate pipeline from
`pipeline/main.py`, since it needs a GPU and, on the HPC cluster, a
different Python environment than the main pipeline. The scripts live in
`deeploc/` and must be run before `pipeline/main.py` so that
`data/secretome/{taxid}.fasta` exists when step 3 runs.

There are two ways to run it, both ending at the same place
(`data/secretome/{taxid}.fasta`):

**A — locally in one step (`deeploc/run_deeploc.py`)**
Downloads each parasite's full proteome from STRING, runs
`deeploc2 -f ... -d mps` (Apple GPU) directly via subprocess, filters the
results, and writes the filtered FASTA. Good for a laptop with an M-series
GPU; slow for 24 large proteomes.

**B — on an HPC / LSF cluster (three scripts, run in order)**
1. `deeploc/prepare_deeploc_input.py` — downloads full proteomes from
   STRING and writes uncompressed FASTAs to `data/deeploc/input/`, meant to
   be copied to the cluster.
2. `deeploc/generate_hpc_jobs.py` — writes one LSF (`bsub`) job script per
   parasite to `jobs/`, each running `deeploc2 -d cuda` on a GPU node.
   Submit with `for f in jobs/*.sh; do bsub < $f; done`.
3. `deeploc/build_secretome_fastas.py` — once job outputs
   (`results_*.csv`) are back locally, reads each CSV and filters the
   source FASTA down to the surviving protein IDs, writing
   `data/secretome/{taxid}.fasta`.

**Localization filter** (shared by both paths,
`filter_by_localization`/`filter_by_deeploc`): a protein is kept if DeepLoc
labels it `Extracellular`, or — only if the parasite is not flagged
`multicellular: true` in `config.yml` — `Cell membrane`. No parasite in the
current `config.yml` sets `multicellular`, so all 24 currently get both
categories; the flag exists to let a multicellular species (e.g. a
helminth) be restricted to `Extracellular` only, but isn't in use yet.
