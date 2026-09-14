'''Builds the EggNOG 6 orthology-group members file used by the prediction pipeline.'''

import argparse
import gzip
import os

import requests

import utils

LEVEL = "2759"  # Eukaryota


def _emit(line, out):
    '''Rewrite one EggNOG 6 row into EggNOG 5 layout.'''
    data = line.rstrip("\n").split("\t")
    if len(data) < 6 or data[0] != LEVEL:
        return False

    group, species, proteins = data[1], data[4], data[5]
    out.write(
        "\t".join(
            [
                data[0],
                group,
                str(len(proteins.split(","))),
                str(len(species.split(","))),
                proteins,
                species,
            ]
        )
        + "\n"
    )
    return True


def build_members_file(url, output_filepath, source_filepath=None):
    '''Writes the level-2759 EggNOG groups to output_filepath (gzipped, EggNOG 5 layout).'''
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)
    kept = seen = 0

    # written to a temporary file and renamed, so a dropped connection leaves the existing
    # file untouched
    tmp_filepath = output_filepath + ".tmp"
    try:
        with gzip.open(tmp_filepath, "wt") as out:
            if source_filepath is not None:
                print(f"  Reading {source_filepath}...")
                with open(source_filepath) as handle:
                    for line in handle:
                        seen += 1
                        kept += _emit(line, out)
            else:
                print(f"  Streaming {url} (~10 GB, filtering as it downloads)...")
                with requests.get(url, stream=True, timeout=60) as r:
                    r.raise_for_status()
                    # the server sends no charset, so decode_unicode would yield bytes
                    r.encoding = r.encoding or "utf-8"
                    for raw in r.iter_lines(decode_unicode=True):
                        if not raw:
                            continue
                        seen += 1
                        kept += _emit(raw, out)
                        if seen % 2_000_000 == 0:
                            print(f"    {seen:,} groups scanned, {kept:,} at level {LEVEL}")
    except BaseException:
        if os.path.exists(tmp_filepath):
            os.remove(tmp_filepath)
        raise

    if kept == 0:
        os.remove(tmp_filepath)
        print(f"  {seen:,} groups scanned, none at level {LEVEL} — {output_filepath} left unchanged")
        return 0

    os.replace(tmp_filepath, output_filepath)
    print(f"  {seen:,} groups scanned, {kept:,} kept at level {LEVEL}")
    print(f"  Wrote {output_filepath}")

    return kept


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config.yml")
    parser.add_argument("--data-dir", default="data")
    parser.add_argument(
        "--source-file",
        help="path to an already-downloaded e6.og2seqs_and_species.tsv",
    )
    args = parser.parse_args()

    urls = utils.read_config(filepath=args.config, field="urls")
    output_filepath = os.path.join(args.data_dir, "downloads", "2759_members.tsv.gz")

    print("Building EggNOG 6 members file...")
    kept = build_members_file(
        url=urls["eggNOG_members_url"],
        output_filepath=output_filepath,
        source_filepath=args.source_file,
    )
    if kept == 0:
        raise SystemExit(f"No groups found at level {LEVEL} — check eggNOG_members_url")
