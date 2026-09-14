'''Build secretome FASTA files from DeepLoc 2 output.'''

import argparse
import glob
import gzip
import os

import pandas as pd

# DeepLoc 2 Accurate-model thresholds from DeepLoc2/deeploc2.py label_threshold, which is
# read at i+1: Extracellular is labels[2], Cell membrane labels[3]. See docs/deeploc.md
EXTRACELLULAR_CUTOFF = 0.61728516
MEMBRANE_CUTOFF = 0.56464844
from Bio import SeqIO

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils


def load_deeploc_csv(deeploc_dir, taxid):
    pattern = os.path.join(deeploc_dir, str(taxid), "results_*.csv")
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No DeepLoc results found for taxid {taxid} in {deeploc_dir}")
    return matches[-1]


def filter_by_score(df, multicellular, extracellular_cutoff, membrane_cutoff):
    '''Keep proteins whose per-class probability clears DeepLoc's threshold.'''
    keep = df["Extracellular"] > extracellular_cutoff
    if not multicellular:
        keep = keep | (df["Cell membrane"] > membrane_cutoff)
    return set(df.loc[keep, "Protein_ID"])


def write_filtered_fasta(source_fasta, valid_ids, output_path):
    opener = gzip.open if source_fasta.endswith(".gz") else open
    with opener(source_fasta, "rt") as handle:
        records = [r for r in SeqIO.parse(handle, "fasta") if r.id in valid_ids]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    SeqIO.write(records, output_path, "fasta")
    return len(records)


def source_fasta_path(config_file, data_dir, taxid):
    '''Path to the STRING sequence FASTA for a species, named after string_sequences_url.'''
    urls = utils.read_config(filepath=config_file, field="urls")
    filename = urls["string_sequences_url"].split("/")[-1].replace("TAXID", str(taxid))
    species_dir = os.path.join(data_dir, "downloads", "species", str(taxid))
    uncompressed = os.path.join(species_dir, filename.removesuffix(".gz"))

    return uncompressed if os.path.exists(uncompressed) else os.path.join(species_dir, filename)


def main(config_file, data_dir, deeploc_dir, extracellular_cutoff, membrane_cutoff):
    parasites = utils.read_config(filepath=config_file, field="parasites")
    out_dir = os.path.join(data_dir, "secretome")

    for taxid, info in parasites.items():
        label = info.get("label", str(taxid))
        multicellular = info.get("multicellular", False)
        print(f"\n[{taxid}] {label} ({'multicellular' if multicellular else 'unicellular'})")

        try:
            csv_path = load_deeploc_csv(deeploc_dir, taxid)
        except FileNotFoundError as e:
            print(f"  SKIP: {e}")
            continue

        df = pd.read_csv(csv_path)
        print(f"  Total proteins in DeepLoc output: {len(df)}")
        valid_ids = filter_by_score(df, multicellular, extracellular_cutoff, membrane_cutoff)
        print(f"  Proteins passing score filter "
              f"(Extracellular>{extracellular_cutoff}"
              f"{'' if multicellular else f', Cell membrane>{membrane_cutoff}'}): {len(valid_ids)}")

        source_fasta = source_fasta_path(config_file, data_dir, taxid)
        if not os.path.exists(source_fasta):
            print(f"  SKIP: source FASTA not found: {source_fasta}")
            continue

        output_fasta = os.path.join(out_dir, f"{taxid}.fasta")
        n = write_filtered_fasta(source_fasta, valid_ids, output_fasta)
        print(f"  Written {n} sequences to {output_fasta}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build secretome FASTAs from DeepLoc 2 output")
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--deeploc-dir",
                        default="data/deeploc/output_accurate/deeploc_output_accurate")
    parser.add_argument("--extracellular-cutoff", type=float, default=EXTRACELLULAR_CUTOFF,
                        help=f"Minimum P(Extracellular), exclusive. Default: {EXTRACELLULAR_CUTOFF}")
    parser.add_argument("--membrane-cutoff", type=float, default=MEMBRANE_CUTOFF,
                        help=f"Minimum P(Cell membrane), exclusive, unicellular only. "
                             f"Default: {MEMBRANE_CUTOFF}")
    args = parser.parse_args()
    main(args.config, args.data_dir, args.deeploc_dir,
         args.extracellular_cutoff, args.membrane_cutoff)
