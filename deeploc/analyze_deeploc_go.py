"""
Compare DeepLoc surface-localization scores against GO Cellular Component terms for a
species, to help choose a score threshold for host filtering.

For each protein it computes a "surface score" = max(Extracellular, Cell membrane)
DeepLoc probability, and flags the protein as "GO-surface" if it carries a Cellular
Component GO term matching a surface/secreted keyword (extracellular, plasma membrane,
cell surface, ...). It then reports, across a sweep of score thresholds:
  - how many proteins pass the DeepLoc threshold,
  - precision: fraction of passing proteins that are GO-surface,
  - recall:    fraction of all GO-surface proteins that pass.

GO Cellular Component terms are read from STRING's enrichment file (gos.parquet only
keeps Biological Process, which is not localization). With --out, also writes a
per-protein table for closer inspection.

Usage:
    python deeploc/analyze_deeploc_go.py --taxid 9606 \
        --deeploc-dir data/deeploc/output [--out human_deeploc_go.csv]
"""

import argparse
import glob
import gzip
import os

import pandas as pd

# DeepLoc localization classes treated as "reachable by a parasite" (secreted / surface)
SURFACE_CLASSES = ["Extracellular", "Cell membrane"]

# substrings marking a Cellular Component GO term as surface / secreted
SURFACE_GO_KEYWORDS = ["extracellular", "plasma membrane", "cell surface",
                       "external side", "cell periphery", "anchored component of plasma membrane",
                       "external encapsulating structure"]

STRING_CC_CATEGORY = "Cellular Component (Gene Ontology)"


def load_deeploc_results(deeploc_dir, taxid):
    """Load the DeepLoc results CSV for a taxid (newest if several)."""
    matches = sorted(glob.glob(os.path.join(deeploc_dir, str(taxid), "results_*.csv")))
    if not matches:
        raise FileNotFoundError(f"No DeepLoc results for taxid {taxid} in {deeploc_dir}")
    df = pd.read_csv(matches[-1])
    missing = [c for c in ["Protein_ID"] + SURFACE_CLASSES if c not in df.columns]
    if missing:
        raise ValueError(f"DeepLoc CSV missing columns {missing}; found {list(df.columns)}")
    df["surface_score"] = df[SURFACE_CLASSES].max(axis=1)
    return df[["Protein_ID", "Localizations", "surface_score"]] if "Localizations" in df.columns \
        else df[["Protein_ID", "surface_score"]]


def load_go_surface_flags(enrichment_file, taxid):
    """
    Return (surface_proteins, cc_terms_per_protein):
      surface_proteins: set of protein ids with a surface/secreted CC GO term
      cc_terms_per_protein: {protein: set(cc term descriptions)}
    """
    opener = gzip.open if enrichment_file.endswith(".gz") else open
    cc_terms = {}
    surface = set()
    with opener(enrichment_file, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        pid_i, cat_i, desc_i = idx["#string_protein_id"], idx["category"], idx["description"]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[cat_i] != STRING_CC_CATEGORY:
                continue
            pid, desc = parts[pid_i], parts[desc_i]
            cc_terms.setdefault(pid, set()).add(desc)
            if any(k in desc.lower() for k in SURFACE_GO_KEYWORDS):
                surface.add(pid)
    return surface, cc_terms


def enrichment_path(data_dir, taxid):
    return os.path.join(data_dir, "downloads", "species", str(taxid),
                        f"{taxid}.protein.enrichment.terms.v12.0.txt.gz")


def main(taxid, deeploc_dir, data_dir, enrichment_file, out):
    deeploc = load_deeploc_results(deeploc_dir, taxid)
    enrichment_file = enrichment_file or enrichment_path(data_dir, taxid)
    surface_proteins, cc_terms = load_go_surface_flags(enrichment_file, taxid)

    deeploc["go_surface"] = deeploc["Protein_ID"].isin(surface_proteins)
    total_go_surface = deeploc["go_surface"].sum()
    print(f"taxid {taxid}: {len(deeploc)} proteins scored by DeepLoc, "
          f"{total_go_surface} carry a surface/secreted CC GO term\n")

    print(f"{'cutoff':>7}{'n_pass':>9}{'precision':>11}{'recall':>9}")
    print("-" * 36)
    for cutoff in [round(0.1 * i, 1) for i in range(1, 10)]:
        passing = deeploc[deeploc["surface_score"] >= cutoff]
        n = len(passing)
        prec = passing["go_surface"].mean() if n else 0.0
        rec = passing["go_surface"].sum() / total_go_surface if total_go_surface else 0.0
        print(f"{cutoff:>7}{n:>9}{prec:>11.2f}{rec:>9.2f}")

    if out:
        deeploc["cc_go_terms"] = deeploc["Protein_ID"].map(
            lambda p: "; ".join(sorted(cc_terms.get(p, set()))))
        deeploc.sort_values("surface_score", ascending=False).to_csv(out, index=False)
        print(f"\nWrote per-protein table to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare DeepLoc surface scores with GO CC terms")
    parser.add_argument("--taxid", type=int, required=True)
    parser.add_argument("--deeploc-dir", default="data/deeploc/output")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--enrichment-file", help="override path to STRING enrichment.terms file")
    parser.add_argument("--out", help="optional path to write the per-protein table as CSV")
    args = parser.parse_args()
    main(args.taxid, args.deeploc_dir, args.data_dir, args.enrichment_file, args.out)
