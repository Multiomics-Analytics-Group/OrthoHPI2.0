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


def select_species(config_file, taxids):
    """Return {taxid: info} to run for. Defaults to config parasites; if taxids is
    given, pick those from hosts+parasites (so host taxids can be targeted too)."""
    hosts = utils.read_config(filepath=config_file, field="hosts") or {}
    parasites = utils.read_config(filepath=config_file, field="parasites") or {}
    if taxids is None:
        return parasites
    species = {**hosts, **parasites}
    return {t: species.get(t, {}) for t in taxids}


def main(config_file, data_dir, taxids=None):
    selected = select_species(config_file, taxids)
    urls = utils.read_config(filepath=config_file, field="urls")
    seq_url_template = urls["string_sequences_url"]
    out_dir = os.path.join(data_dir, "deeploc", "input")
    os.makedirs(out_dir, exist_ok=True)

    for taxid, info in selected.items():
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
    parser.add_argument("--taxids", help="comma-separated taxids to run for (default: config parasites); "
                                         "accepts host taxids too, e.g. 9606,10116,9823")
    args = parser.parse_args()
    taxids = [int(t) for t in args.taxids.split(",")] if args.taxids else None
    main(args.config, args.data_dir, taxids=taxids)
