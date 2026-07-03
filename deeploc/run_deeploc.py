"""
Run DeepLoc 2.1 on parasite proteomes to predict secretome and membrane proteins,
producing FASTA files in data/secretome_pred_input_data/input_data/{taxid}.fasta
that are consumed by get_secretome_predictions() in filters.py.

Usage:
    python run_deeploc.py [--config config.yml] [--data-dir data]

Requires deeploc2 installed:
    pip install . (from the deeploc2_package directory)

For each parasite in config:
  - Downloads protein sequences from STRING (TAXID.protein.sequences.v11.5.fa.gz)
  - Runs DeepLoc 2.1 (Fast model) on the full proteome
  - Keeps Extracellular proteins for all parasites, plus Cell membrane for unicellular
  - Writes filtered sequences to data/secretome_pred_input_data/input_data/{taxid}.fasta
"""

import argparse
import glob
import gzip
import os
import shutil
import subprocess

import pandas as pd
from Bio import SeqIO

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils

def download_sequences(taxid, data_dir, url_template):
    url = url_template.replace("TAXID", str(taxid))
    gz_path = utils.download_file(url=url, data_dir=data_dir)
    fa_path = gz_path.replace(".gz", "")
    if not os.path.exists(fa_path):
        print(f"  Decompressing {gz_path}")
        with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return fa_path


def run_deeploc(fasta_path, output_dir):
    existing = sorted(glob.glob(os.path.join(output_dir, "results_*.csv")))
    if existing:
        print(f"  Using existing DeepLoc output: {existing[-1]}")
        return existing[-1]

    os.makedirs(output_dir, exist_ok=True)
    print(f"  Running DeepLoc 2.1 on {fasta_path}")
    subprocess.run(["deeploc2", "-f", fasta_path, "-o", output_dir, "-d", "mps"], check=True)

    csvs = sorted(glob.glob(os.path.join(output_dir, "results_*.csv")))
    if not csvs:
        raise FileNotFoundError(f"DeepLoc produced no output CSV in {output_dir}")
    return csvs[-1]


def filter_by_deeploc(csv_path, multicellular):
    df = pd.read_csv(csv_path)
    target_locs = {"Extracellular"}
    if not multicellular:
        target_locs.add("Cell membrane")

    valid = set()
    for _, row in df.iterrows():
        locs = set(row["Localizations"].split("|"))
        if locs & target_locs:
            valid.add(row["Protein_ID"])
    return valid


def write_filtered_fasta(fa_path, valid_ids, output_path):
    records = [r for r in SeqIO.parse(fa_path, "fasta") if r.id in valid_ids]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    SeqIO.write(records, output_path, "fasta")
    return len(records)


def main(config_file, data_dir):
    parasites = utils.read_config(filepath=config_file, field="parasites")
    urls = utils.read_config(filepath=config_file, field="urls")
    seq_url_template = urls["string_sequences_url"]
    out_dir = os.path.join(data_dir, "secretome")
    deeploc_dir = os.path.join(data_dir, "deeploc", "output")

    for taxid, info in parasites.items():
        label = info.get("label", taxid)
        multicellular = info.get("multicellular", False)
        print(f"\n[{taxid}] {label} ({'multicellular' if multicellular else 'unicellular'})")

        fa_path = download_sequences(taxid, os.path.join(data_dir, "downloads"), seq_url_template)
        csv_path = run_deeploc(fa_path, os.path.join(deeploc_dir, str(taxid)))

        valid_ids = filter_by_deeploc(csv_path, multicellular)
        print(f"  Proteins passing DeepLoc filter: {len(valid_ids)}")

        output_fasta = os.path.join(out_dir, f"{taxid}.fasta")
        n = write_filtered_fasta(fa_path, valid_ids, output_fasta)
        print(f"  Written {n} sequences to {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run DeepLoc 2.1 to predict secretome/membrane proteins")
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--data-dir", default="data")
    args = parser.parse_args()
    main(args.config, args.data_dir)
