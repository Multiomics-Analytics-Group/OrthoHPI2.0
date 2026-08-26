# Prediction pipeline

This document describes what happens when you run `python -m pipeline.main`
(the `main.py` invoked from the README as `python main.py`) — how the
host-parasite PPI predictions in `data/annotated_predictions.parquet` get
built. 

Everything is driven by `config.yml`: `urls` (data sources), `hosts` (human
9606, rat 10116, mouse 10090, pig 9823), and `parasites` (40 species, each with
a STRING taxid, a display label/color, the list of BTO tissue codes relevant to
that parasite's life cycle, and a `hosts` list of the host taxids it infects).

A parasite is only paired with the hosts named in its `hosts` list
(`homology.get_links`), so e.g. a pig-only parasite never generates
human predictions.

A host may carry an optional `group`, which only affects the app: hosts sharing
a group are offered as a single option in the host selector and their
predictions are pooled, while the predictions themselves keep their own species
label, taxid and colour. Rat and mouse share `group: Rodent` — `Parasites_db.xlsx`
only curates rat as a host, and we assume the same parasites infect mouse, so
every parasite with 10116 in its `hosts` list also lists 10090. Mouse is worth
having alongside rat because rat's jensenlab tissue data is sparse (no blood,
skin, macrophage, nose or mouth above the 2.5 cutoff) where mouse's is not.

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

Runs once per host: downloads that host's tissue expression scores from its
`tissues_url` (jensenlab tissues.jensenlab.org) and keeps a host protein only
if:
- its expression score in some tissue is `>= cutoff` (2.5 in `main.py`), and
- that tissue is one of the BTO codes listed under the relevant parasite(s)
  in `config.yml` (i.e. a tissue the parasite actually encounters during its
  life cycle — gut, skin, blood, etc.).

Every match is recorded in a `tissues` dict (`{protein: [tissue, ...]}`),
aggregated across all hosts, and later joined with HPA cell-type data in
`main.get_tissue_cell_type_annotation`. HPA is human-only, so non-human host
proteins get no cell-type resolution (left join → NaN). Note that not every
BTO tissue exists in every host's jensenlab data (e.g. rat lacks skin and
blood), so a parasite whose tissues are all absent for its host yields no
predictions.

### 5. DeepLoc host filter (`filters.apply_deeploc_filter`)

Keeps only surface-exposed host proteins, using DeepLoc 2 (Accurate)
subcellular-localization predictions read from
`data/deeploc/output_accurate/deeploc_output_accurate/<taxid>/results_*.csv`.
A host protein passes if `P(Extracellular) >= 0.617` or
`P(Cell membrane) >= 0.565` — i.e. the protein is at the cell surface or
secreted, somewhere the parasite could physically reach it. Hosts without a
DeepLoc run are left unfiltered (with a warning).

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

The alias `sources` are a preference list — `parse_string_aliases` takes the
alias of the first source that has one — because STRING's alias files are not
uniform across species. Parasites use `Uniprot`. Hosts use
`Ensembl_HGNC_uniprot_ids`, one canonical accession per protein but present
only in human's file, falling back to `UniProt_AC` for rat, mouse and pig,
which would otherwise get no accession at all (and so no structures on the
app's structure page). `Ensembl_UniProt` looks like a third option but must not
be used: it mixes gene names in among the accessions.

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

## Host orthologs for the cross-host page (`scripts/build_host_orthologs.py`)

The app's *One parasite, several hosts* page compares the same parasite in
two or more hosts, and has to separate two reasons an interaction predicted
in one host can be missing from another: the host has no protein in that
orthology group at all, or it has one and the pipeline's filters dropped it
before the transfer in step 7. `predictions.parquet` cannot tell them apart —
a group with no predicted interaction is simply absent from it either way.

This script streams `data/downloads/2759_members.tsv.gz` (the same members
file step 7 reads) and writes, for the host orthology groups that appear in
`predictions.parquet`, which proteins of each host belong to them:

```
$ python scripts/build_host_orthologs.py
```

Output: `data/host_orthologs.parquet` (`group`, `taxid`, `n_proteins`,
`proteins`) — a few hundred rows, since only the ~200 groups the predictions
reach are kept. Re-run it after `pipeline/main.py` whenever the hosts or the
parasites change. The page treats a missing file as "not known" and leaves
that section out rather than failing.
