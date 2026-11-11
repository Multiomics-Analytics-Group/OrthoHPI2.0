"""Overlap between the predictions and curated host-parasite interactions.

Queries IntAct for every host-parasite pair in config.yml, keeps the curated
inter-species interactions, and compares them with the predictions at three
levels: the exact protein pair, the pair of orthologous groups the prediction
was transferred from, and the host protein alone.

Both sides are resolved to STRING identifiers through the STRING alias files,
since IntAct names a protein by whatever the curating publication used -- a
UniProt accession of another isolate, a locus tag of another strain, a gene
symbol -- and none of those are the names the pipeline carries. A curated
protein that resolves to no STRING id is absent from the reference proteome and
could never have been predicted, which is reported apart from the proteins the
filters removed.

Writes paper/tables/known_ppis.csv (one row per curated interaction, with how
it matches) and paper/tables/known_ppis_summary.csv (one row per host-parasite
pair with curated data).

    .venv/bin/python scripts/build_known_ppis.py [--cache DIR]

The IntAct responses are cached, so a second run is offline.
"""

import argparse
import time
from glob import glob
from pathlib import Path
from urllib.parse import quote
from urllib.request import urlopen

import pandas as pd
import yaml

PSICQUIC = (
    "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/"
    "current/search/query/"
)
# PSI-MITAB 2.5 columns we use, by position.
COLS = {
    "a": 0,
    "b": 1,
    "alt_a": 2,
    "alt_b": 3,
    "alias_a": 4,
    "alias_b": 5,
    "method": 6,
    "author": 7,
    "pmid": 8,
    "taxid_a": 9,
    "taxid_b": 10,
    "type": 11,
}
ROOT = Path(__file__).resolve().parents[1]
SPECIES = ROOT / "data" / "downloads" / "species"


def fetch(query, cache_dir, pause=0.5):
    """The MITAB rows IntAct returns for one MIQL query, cached by query."""
    cache = cache_dir / (quote(query, safe="") + ".tab")
    if cache.exists():
        return cache.read_text().splitlines()
    url = PSICQUIC + quote(query, safe=":() ").replace(" ", "%20")
    with urlopen(url, timeout=120) as response:
        text = response.read().decode("utf-8", "replace")
    cache.write_text(text)
    time.sleep(pause)
    return text.splitlines()


def taxid(field):
    """The numeric taxid in a MITAB organism field, or None."""
    if not field or field == "-":
        return None
    first = field.split("|")[0]
    if not first.startswith("taxid:"):
        return None
    digits = first[len("taxid:") :].split("(")[0]
    return int(digits) if digits.lstrip("-").isdigit() else None


def accessions(*fields):
    """Every UniProt accession named in the given MITAB identifier fields."""
    found = set()
    for field in fields:
        if not field or field == "-":
            continue
        for part in field.split("|"):
            if part.startswith("uniprotkb:"):
                value = part[len("uniprotkb:") :].split("(")[0]
                # Drop isoform and chain suffixes: P12345-2, P12345-PRO_0000001.
                found.add(value.split("-")[0])
    return found


def names(*fields):
    """Every alias in the given MITAB alias fields, upper-cased."""
    found = set()
    for field in fields:
        if not field or field == "-":
            continue
        for part in field.split("|"):
            if "(" in part and ":" in part:
                value = part.split(":", 1)[1].split("(")[0]
                found.add(value.upper())
    return found


class Proteomes:
    """
    The STRING identifiers of each species, by every alias STRING lists for them.

    STRING's alias file carries the accessions, locus tags, gene symbols and
    database cross-references of a protein under one id, which is what makes a
    curated interaction comparable with a prediction at all.
    """

    def __init__(self, species_dir=SPECIES):
        self.species_dir = Path(species_dir)
        self._index = {}

    def index(self, taxid):
        taxid = str(taxid)
        if taxid not in self._index:
            self._index[taxid] = self._read(taxid)
        return self._index[taxid]

    def _read(self, taxid):
        files = sorted(
            glob(str(self.species_dir / taxid / f"{taxid}.protein.aliases.*.txt.gz"))
        )
        if not files:
            print(f"  ! no STRING alias file for taxid {taxid}")
            return {}
        rows = pd.read_csv(
            files[-1],
            sep="\t",
            header=0,
            names=["protein", "alias", "source"],
            dtype=str,
        )
        index = {}
        for protein, alias in zip(rows.protein, rows.alias):
            if isinstance(alias, str):
                index.setdefault(alias.upper(), set()).add(protein)
        # A STRING id resolves to itself, so a caller may pass either.
        for protein in rows.protein.dropna().unique():
            index.setdefault(protein.upper(), set()).add(protein)
        return index

    def resolve(self, taxid, *identifiers):
        """The STRING ids of one species that any of these identifiers name."""
        index = self.index(taxid)
        found = set()
        for group in identifiers:
            for identifier in group:
                found |= index.get(str(identifier).upper(), set())
        return found


def curated_pairs(config, cache_dir):
    """Curated inter-species interactions for every host-parasite pair."""
    hosts = {int(t): c["label"].split(" (")[0] for t, c in config["hosts"].items()}
    rows = []
    for parasite_taxid, parasite in config["parasites"].items():
        species = parasite["label"]
        for host_taxid in parasite["hosts"]:
            host = hosts[int(host_taxid)]
            query = f'species:"{species}" AND species:"{host}"'
            try:
                lines = fetch(query, cache_dir)
            except Exception as error:  # noqa: BLE001 - report and carry on
                print(f"  ! {species} / {host}: {error}")
                continue
            for line in lines:
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) <= COLS["type"]:
                    continue
                ta, tb = taxid(f[COLS["taxid_a"]]), taxid(f[COLS["taxid_b"]])
                if ta is None or tb is None or ta == tb:
                    continue
                # One side must be the host; the other is then taken as the
                # parasite, since the query already fixed the two species.
                if ta == int(host_taxid):
                    side = "b"
                elif tb == int(host_taxid):
                    side = "a"
                else:
                    continue
                other = "a" if side == "b" else "b"
                row = {
                    "parasite": species,
                    "parasite_taxid": int(parasite_taxid),
                    "host": host,
                    "host_taxid": int(host_taxid),
                    "parasite_acc": sorted(
                        accessions(f[COLS[side]], f[COLS["alt_" + side]])
                    ),
                    "host_acc": sorted(
                        accessions(f[COLS[other]], f[COLS["alt_" + other]])
                    ),
                    "parasite_names": sorted(names(f[COLS["alias_" + side]])),
                    "host_names": sorted(names(f[COLS["alias_" + other]])),
                    "method": f[COLS["method"]],
                    "type": f[COLS["type"]],
                    "pmid": f[COLS["pmid"]],
                }
                rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache",
        type=Path,
        default=ROOT / "data" / "tmp" / "intact",
        help="where to cache the IntAct responses",
    )
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)

    config = yaml.safe_load((ROOT / "config.yml").read_text())
    niches = {
        str(t): str(c.get("niche", "")).strip().lower()
        for t, c in config["parasites"].items()
    }
    proteomes = Proteomes()

    print("Querying IntAct for every host-parasite pair...")
    curated = curated_pairs(config, args.cache)
    if curated.empty:
        print("No curated inter-species interactions found.")
        return

    predictions = pd.read_parquet(ROOT / "data" / "predictions.parquet")
    # The proteins that reached the transfer step, to tell a protein the filters
    # removed from one they kept but could not transfer an interaction to. A
    # parasite protein carries both niche columns, so membership is enough for it;
    # a host protein counts only in the half of the pool its parasite's niche reaches.
    eligible = pd.read_parquet(ROOT / "data" / "eligible_proteins.parquet")
    eligible_ids = set(eligible.protein)
    reachable = {
        niche: set(eligible.protein[eligible[f"reachable_{niche}"]])
        for niche in ("intracellular", "extracellular")
        if f"reachable_{niche}" in eligible.columns
    }
    # The orthologous groups of each STRING id, for the group-level comparison.
    groups = {}
    for ids, column in (
        (predictions.source, predictions.group1),
        (predictions.target, predictions.group2),
    ):
        for protein, group in zip(ids, column):
            if isinstance(protein, str):
                groups.setdefault(protein, set()).add(group)

    predicted_pairs = set(zip(predictions.source, predictions.target))
    predicted_pairs |= {(b, a) for a, b in predicted_pairs}
    predicted_groups = set(zip(predictions.group1, predictions.group2))
    predicted_groups |= {(b, a) for a, b in predicted_groups}
    # The host proteins reached by any parasite, and by each parasite on its own.
    any_hosts = set(predictions.target.dropna())
    by_parasite = {
        species.split(" (")[0]: set(rows.target.dropna())
        for species, rows in predictions.groupby("taxid1_label")
    }

    def match(row):
        parasite_ids = proteomes.resolve(
            row.parasite_taxid, row.parasite_acc, row.parasite_names
        )
        host_ids = proteomes.resolve(row.host_taxid, row.host_acc, row.host_names)
        pgroups = {g for i in parasite_ids for g in groups.get(i, ())}
        hgroups = {g for i in host_ids for g in groups.get(i, ())}
        host_pool = reachable.get(niches.get(str(row.parasite_taxid)), eligible_ids)
        return pd.Series(
            {
                "parasite_string": ";".join(sorted(parasite_ids)),
                "host_string": ";".join(sorted(host_ids)),
                # Whether each side is in the reference proteome at all, which a
                # curated protein of another strain or isolate is not.
                "parasite_in_proteome": bool(parasite_ids),
                "host_in_proteome": bool(host_ids),
                # Whether each side survived the secretome, tissue and DeepLoc filters.
                "parasite_eligible": bool(parasite_ids & eligible_ids),
                "host_eligible": bool(host_ids & host_pool),
                "pair_predicted": any(
                    (p, h) in predicted_pairs for p in parasite_ids for h in host_ids
                ),
                "group_predicted": any(
                    (p, h) in predicted_groups for p in pgroups for h in hgroups
                ),
                # The same host protein reached by the same parasite, and by any.
                "host_predicted": bool(
                    host_ids & by_parasite.get(row.parasite, set())
                ),
                "host_predicted_any": bool(host_ids & any_hosts),
            }
        )

    curated = pd.concat([curated, curated.apply(match, axis=1)], axis=1)
    # IntAct reports one row per detection method and publication, so the same
    # interaction can appear several times. Count each distinct pair of
    # proteins once, keeping the methods and publications that support it.
    curated["key"] = [
        ";".join(sorted(p)) + "|" + ";".join(sorted(h))
        for p, h in zip(curated.parasite_acc, curated.host_acc)
    ]
    flags = [
        "pair_predicted",
        "group_predicted",
        "parasite_in_proteome",
        "host_in_proteome",
        "parasite_eligible",
        "host_eligible",
        "host_predicted",
        "host_predicted_any",
    ]
    curated = (
        curated.groupby(["parasite", "host", "key"], as_index=False)
        .agg(
            {
                "parasite_taxid": "first",
                "host_taxid": "first",
                "parasite_acc": "first",
                "host_acc": "first",
                "parasite_names": "first",
                "host_names": "first",
                "parasite_string": "first",
                "host_string": "first",
                "method": lambda s: ";".join(sorted(set(s))),
                "type": lambda s: ";".join(sorted(set(s))),
                "pmid": lambda s: ";".join(sorted(set(s))),
                **{flag: "max" for flag in flags},
            }
        )
        .drop(columns="key")
    )

    out = ROOT / "paper" / "tables"
    flat = curated.assign(
        parasite_acc=curated.parasite_acc.str.join(";"),
        host_acc=curated.host_acc.str.join(";"),
        parasite_names=curated.parasite_names.str.join(";"),
        host_names=curated.host_names.str.join(";"),
    )
    flat.to_csv(out / "known_ppis.csv", index=False)

    # The curated record repeats one parasite protein across many host preys, so
    # the number of distinct parasite proteins is counted beside the interactions.
    summary = (
        curated.groupby(["parasite", "host"])
        .agg(
            curated=("pmid", "size"),
            parasite_proteins=("parasite_string", "nunique"),
            publications=("pmid", lambda s: len({p for v in s for p in v.split(";")})),
            pair_predicted=("pair_predicted", "sum"),
            group_predicted=("group_predicted", "sum"),
            host_predicted=("host_predicted", "sum"),
            host_predicted_any=("host_predicted_any", "sum"),
            parasite_in_proteome=("parasite_in_proteome", "sum"),
            host_in_proteome=("host_in_proteome", "sum"),
            parasite_eligible=("parasite_eligible", "sum"),
            host_eligible=("host_eligible", "sum"),
        )
        .reset_index()
        .sort_values("curated", ascending=False)
    )
    summary.to_csv(out / "known_ppis_summary.csv", index=False)

    proteins = curated.drop_duplicates(["parasite", "parasite_string", "pmid"])
    print(f"\n{len(curated)} distinct curated inter-species interactions")
    print(f"  from {len(proteins)} distinct parasite proteins")
    print(f"  exact protein pair also predicted:   {curated.pair_predicted.sum()}")
    print(f"  orthologous-group pair predicted:    {curated.group_predicted.sum()}")
    print(f"  host protein reached, same parasite: {curated.host_predicted.sum()}")
    print(f"  host protein reached, any parasite:  {curated.host_predicted_any.sum()}")

    print("\nThe parasite proteins, and how far each got:")
    for _, row in proteins.iterrows():
        if not row.parasite_in_proteome:
            state = "not in the STRING proteome"
        elif not row.parasite_eligible:
            state = "in the proteome, removed by the filters"
        elif not row.pair_predicted:
            state = "eligible, but no transferable orthologous group"
        else:
            state = "predicted"
        label = ";".join(row.parasite_names[:2]) or ";".join(row.parasite_acc) or "?"
        print(f"  {row.parasite:<22} {label:<34} {state}")

    print(f"\n{summary.to_string(index=False)}")
    print(f"\nwrote {out / 'known_ppis.csv'} and {out / 'known_ppis_summary.csv'}")


if __name__ == "__main__":
    main()
