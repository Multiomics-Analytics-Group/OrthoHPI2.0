"""
Prepare input FASTA files for DeepLoc 2 on HPC.
Downloads protein sequences from STRING for each parasite and writes the full
proteome to data/secretome_pred_input_data/input_data/{taxid}.fasta.

These files can then be copied to the HPC and used as input for DeepLoc 2.
After DeepLoc runs, use build_secretome_fastas.py to filter by localization.

Usage:
    python prepare_deeploc_input.py [--config config.yml] [--data-dir data]
"""

import argparse
import gzip
import os
import shutil

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


def main(config_file, data_dir):
    parasites = utils.read_config(filepath=config_file, field="parasites")
    urls = utils.read_config(filepath=config_file, field="urls")
    seq_url_template = urls["string_sequences_url"]
    out_dir = os.path.join(data_dir, "deeploc", "input")
    os.makedirs(out_dir, exist_ok=True)

    for taxid, info in parasites.items():
        label = info.get("label", str(taxid))
        print(f"[{taxid}] {label}")

        gz_path = utils.download_file(url=seq_url_template.replace("TAXID", str(taxid)), data_dir=os.path.join(data_dir, "downloads", "species", str(taxid)))

        out_path = os.path.join(out_dir, f"{taxid}.fasta")
        print(f"  Decompressing to {out_path}")
        with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare input FASTAs for DeepLoc 2 on HPC")
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--data-dir", default="data")
    args = parser.parse_args()
    main(args.config, args.data_dir)
