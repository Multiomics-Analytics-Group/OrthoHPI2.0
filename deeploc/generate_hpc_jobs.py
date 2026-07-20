"""
Generate LSF job scripts for running DeepLoc 2.1 on HPC, one per parasite.

Usage:
    python generate_hpc_jobs.py --work-dir /work3/idamei/orthohpi

This assumes:
  - Protein sequence FASTAs are at {work_dir}/data/TAXID.protein.sequences.v12.0.fa
  - DeepLoc venv is at {work_dir}/deeploc_venv
  - Output goes to {work_dir}/data/deeploc_output/TAXID/

Job scripts are written to jobs/ and can be submitted with:
    bsub < jobs/5759_Entamoeba_histolytica.sh
or all at once:
    for f in jobs/*.sh; do bsub < $f; done
"""

import argparse
import os
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils

TEMPLATE = """\
#!/bin/sh
### General options
#BSUB -q gpua100
#BSUB -J {taxid}
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8GB]"
#BSUB -gpu "num=1:mode=exclusive_process"
#BSUB -W 04:00
#BSUB -B
#BSUB -N
#BSUB -o {work_dir}/logs/deeploc_{taxid}_%J.out
#BSUB -e {work_dir}/logs/deeploc_{taxid}_%J.err

module load cuda/12.1
module load python3/3.10.14
source ~/deeploc_env/bin/activate

mkdir -p {work_dir}/data/deeploc_output/{taxid}

echo "Start: $(date)"
time deeploc2 \\
    -f {work_dir}/data/{taxid}.fasta \\
    -o {work_dir}/data/deeploc_output/{taxid} \\
    -d cuda
echo "End: $(date)"
"""


def main(config_file, work_dir):
    parasites = utils.read_config(filepath=config_file, field="parasites")

    jobs_dir = "jobs"
    os.makedirs(jobs_dir, exist_ok=True)

    for taxid, info in parasites.items():
        label = info.get("label", str(taxid)).replace(" ", "_")
        script = TEMPLATE.format(taxid=taxid, work_dir=work_dir)

        filename = os.path.join(jobs_dir, f"{taxid}_{label}.sh")
        with open(filename, "w") as f:
            f.write(script)
        print(f"Written {filename}")

    print(f"\nSubmit all jobs with:\n  for f in jobs/*.sh; do bsub < $f; done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate LSF job scripts for DeepLoc on HPC")
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--work-dir", default="/work3/idamei/orthohpi",
                        help="Base directory on HPC")
    args = parser.parse_args()
    main(args.config, args.work_dir)
