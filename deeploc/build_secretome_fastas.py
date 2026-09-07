"""
Build secretome FASTA files from DeepLoc 2 output.

Reads DeepLoc CSV results (Accurate model) and filters proteins on the per-class
probability, against DeepLoc's own thresholds for that model:
  - P(Extracellular) > 0.61728516 (all parasites)
  - P(Cell membrane) > 0.56464844 (unicellular parasites only)

The same thresholds the host proteins are filtered on (pipeline/main.py). See
docs/deeploc.md for where the numbers come from and why they are not the
assigned Localizations column.

Writes filtered FASTAs to data/secretome_pred_input_data/input_data/{taxid}.fasta

Usage:
    python build_secretome_fastas.py [--config config.yml] [--data-dir data]
        [--deeploc-dir data/deeploc/output_accurate/deeploc_output_accurate]
        [--extracellular-cutoff 0.61728516] [--membrane-cutoff 0.56464844]
"""

import argparse
import glob
import gzip
import os

import pandas as pd

# DeepLoc 2's own per-class thresholds for the Accurate (ProtT5) model, taken from
# DeepLoc2/deeploc2.py label_threshold. That array carries one entry more than there are
# classes and convert_label2string reads it at i+1, so the threshold of labels[i] is
# label_threshold[i+1]: Extracellular is labels[2] -> 0.61728516 and Cell membrane is
# labels[3] -> 0.56464844. Reading the array straight gives the wrong pair, and taking it
# from the Fast model gives another wrong pair (0.64638672 / 0.52368164).
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
    """Keep proteins whose per-class probability clears DeepLoc's threshold.

    Extracellular is kept for every parasite; cell membrane is additionally kept
    for unicellular parasites, a multicellular one reaching its host with
    secreted proteins alone.

    Strictly greater than, which is how DeepLoc itself calls a class
    (DeepLoc2/utils.py convert_label2string). Read the probability rather than
    the Localizations column: DeepLoc writes a class into Localizations even
    when nothing crosses a threshold, falling back to whichever class came
    closest, and those proteins are not surface-exposed in any useful sense.
    """
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
    """Path to the STRING sequence FASTA for a species, named after string_sequences_url."""
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
