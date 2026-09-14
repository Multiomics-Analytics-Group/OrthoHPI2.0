'''Generate LSF job scripts for running DeepLoc 2.1 on HPC, one per parasite.'''

import argparse
import os
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils

RESOURCES = {
    "Fast":     {"mem": "8GB",  "walltime": "01:00"},
    "Accurate": {"mem": "40GB", "walltime": "04:00"},
}

OUTPUT_SUBDIR = {
    "Fast": "deeploc_output",
    "Accurate": "deeploc_output_accurate",
}

TEMPLATE = """\
#!/bin/sh
### General options
#BSUB -q gpua100
#BSUB -J {taxid}_{model}
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem={mem}]"
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W {walltime}
#BSUB -B
#BSUB -N
#BSUB -o {work_dir}/logs/{model_lower}/deeploc_{taxid}_%J.out
#BSUB -e {work_dir}/logs/{model_lower}/deeploc_{taxid}_%J.err

module load cuda/12.1
module load python3/3.10.14
source ~/deeploc_env/bin/activate

# Pin the HuggingFace cache to shared scratch (not $HOME, whose quota a 32GB
# ProtT5 checkpoint can exceed) so it matches wherever it was pre-downloaded on
# the login node -- see the warm-up step in the workflow notes.
export HF_HOME={work_dir}/.cache/huggingface

mkdir -p {work_dir}/data/{output_subdir}/{taxid}

echo "Start: $(date)"
time deeploc2 \\
    -f {work_dir}/data/{taxid}.fasta \\
    -o {work_dir}/data/{output_subdir}/{taxid} \\
    -m {model} \\
    -d cuda
echo "End: $(date)"
"""


def select_species(config_file, taxids):
    '''Return {taxid: info} to run for.'''
    hosts = utils.read_config(filepath=config_file, field="hosts") or {}
    parasites = utils.read_config(filepath=config_file, field="parasites") or {}
    if taxids is None:
        return parasites
    species = {**hosts, **parasites}
    return {t: species.get(t, {}) for t in taxids}


def main(config_file, work_dir, taxids=None, model="Fast"):
    selected = select_species(config_file, taxids)
    resources = RESOURCES[model]
    output_subdir = OUTPUT_SUBDIR[model]
    model_lower = model.lower()

    jobs_dir = os.path.join("jobs", model_lower)
    os.makedirs(jobs_dir, exist_ok=True)

    for taxid, info in selected.items():
        label = info.get("label", str(taxid)).replace(" ", "_")
        script = TEMPLATE.format(taxid=taxid, work_dir=work_dir, model=model,
                                  model_lower=model_lower, output_subdir=output_subdir,
                                  **resources)

        filename = os.path.join(jobs_dir, f"{taxid}_{label}.sh")
        with open(filename, "w") as f:
            f.write(script)
        print(f"Written {filename}")

    print(f"\nBefore submitting, the log directory must exist on the HPC side (LSF needs it "
          f"upfront for -o/-e):\n  mkdir -p {work_dir}/logs/{model_lower}")
    print(f"\nSubmit all jobs with:\n  for f in jobs/{model_lower}/*.sh; do bsub < $f; done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate LSF job scripts for DeepLoc on HPC")
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--work-dir", default="/work3/idamei/orthohpi",
                        help="Base directory on HPC")
    parser.add_argument("--taxids", help="comma-separated taxids to run for (default: config parasites); "
                                         "accepts host taxids too, e.g. 9606,10116,9823")
    parser.add_argument("--model", choices=["Fast", "Accurate"], default="Fast",
                        help="DeepLoc model to run: Fast (ESM1b) or Accurate (ProtT5). Default: Fast")
    args = parser.parse_args()
    taxids = [int(t) for t in args.taxids.split(",")] if args.taxids else None
    main(args.config, args.work_dir, taxids=taxids, model=args.model)
