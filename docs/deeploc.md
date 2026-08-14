# DeepLoc 2.1

DeepLoc predicts subcellular localization and is run on a GPU, outside
`pipeline/main.py`, in its own venv. Two consumers: the parasite secretome
filter and the host surface filter ([pipeline.md](pipeline.md) steps 3 and 5).
A third keeps the probabilities for the app — `build_deeploc_localisations.py`,
run after the main pipeline, which is where the shape of the dots in the
shared-interactors matrix of the "Parasites in a host" page comes from.

## Scripts

| Script | In | Out |
| --- | --- | --- |
| `deeploc/prepare_deeploc_input.py` | STRING proteomes | `data/deeploc/input/{taxid}.fasta` |
| `deeploc/generate_hpc_jobs.py` | `config.yml` | `jobs/{fast,accurate}/{taxid}_{label}.sh` |
| *(DeepLoc on the HPC)* | `WORK_DIR/data/{taxid}.fasta` | `WORK_DIR/data/deeploc_output_accurate/{taxid}/results_*.csv` |
| `deeploc/build_secretome_fastas.py` | results CSVs | `data/secretome/{taxid}.fasta` (parasites) |
| `pipeline/build_deeploc_localisations.py` | results CSVs + `predictions.parquet` | `data/deeploc_localisations.parquet` (app) |

Host results need no post-processing — copy the output dirs to
`data/deeploc/output_accurate/deeploc_output_accurate/{taxid}/` and
`apply_deeploc_filter` reads `results_*.csv` directly.

`deeploc/run_deeploc.py` does everything locally in one step (`-d mps`) — fine
for a few species on a laptop.

## Running

```
python deeploc/prepare_deeploc_input.py
python deeploc/generate_hpc_jobs.py --work-dir WORK_DIR --model Accurate [--taxids 9606,10116,9823]
# copy FASTAs to WORK_DIR/data/, then on the HPC:
mkdir -p WORK_DIR/logs/accurate          # LSF won't create it; logs are lost silently
for f in jobs/accurate/*.sh; do bsub < $f; done
```

`--taxids` accepts host taxids too. Fast (ESM-1b) and Accurate (ProtT5) runs
keep separate job, log and output directories so they never mix.

## HPC setup (once)

Request DeepLoc 2.1 [here](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=deeploc&version=2.1&packageversion=2.1&platform=All),
copy it over, then:

```
module load python3/3.10.14              # same module the job scripts load
python3 -m venv ~/deeploc_env
source ~/deeploc_env/bin/activate
pip install DeepLoc-2.1.0.tar.gz
pip install "transformers==4.45.2"       # after DeepLoc; required for ProtT5/Accurate

# GPU nodes have no internet, so cache ProtT5 on the login node first
export HF_HOME=WORK_DIR/.cache/huggingface && deeploc2 -f test.fasta -o /tmp/warmup -m Accurate
```

`HF_HOME` must be on shared scratch (visible from the compute nodes, and not
your `$HOME` quota); the generated job scripts export the same path. The cache
doesn't expire.

## Thresholds

DeepLoc's own per-class cutoffs for the Accurate model (`label_threshold` in
`DeepLoc2/deeploc2.py` — note the array is offset by one against the label
list). Host values live in `pipeline/main.py`:

```python
DEEPLOC_EXTRACELLULAR_CUTOFF = None       # 0.61728516 to also keep secreted
DEEPLOC_MEMBRANE_CUTOFF = 0.56464844
```

Parasite values are `build_secretome_fastas.py` flags (`--extracellular-cutoff`
/ `--membrane-cutoff`, defaults 0.617 / 0.524); membrane only counts for species
not flagged `multicellular: true` in `config.yml`.
