'''
Which orthology groups each host proteome has a protein in, for the host groups that carry a
predicted interaction.
'''

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils

# columns of the EggNOG 5-layout members file written by pipeline/prepare_eggnog_members.py
GROUP_COLUMN, PROTEINS_COLUMN, SPECIES_COLUMN = 1, 4, 5


def get_host_groups(predictions_file):
    '''
    The orthology groups of the host side of the predictions, and the hosts to look for
    in them.
    '''
    predictions = utils.read_parquet_file(input_file=predictions_file)

    return set(predictions['group2']), set(predictions['taxid2'].astype(str))


def count_members(members_file, groups, taxids):
    '''
    Streams the members file and counts, for each wanted group, how many proteins of
    each host it holds.
    '''
    counts = []
    with utils.read_gzipped_file(members_file) as members:
        for i, line in enumerate(members, 1):
            if i % 500_000 == 0:
                print(f"  {i:,} groups scanned, {len(counts)} host memberships found")
            data = line.decode("utf-8").rstrip("\n").split("\t")
            if len(data) <= SPECIES_COLUMN or data[GROUP_COLUMN] not in groups:
                continue
            present = taxids.intersection(data[SPECIES_COLUMN].split(","))
            if not present:
                continue
            proteins = data[PROTEINS_COLUMN].split(",")
            for taxid in sorted(present):
                # the proteins are STRING-style taxid.identifier
                prefix = f"{taxid}."
                members_of_host = [p for p in proteins if p.startswith(prefix)]
                # the ids are kept so the app can ask which are annotated to an infected
                # tissue
                counts.append([data[GROUP_COLUMN], taxid, len(members_of_host),
                               ",".join(members_of_host)])

    return counts


def main(data_dir):
    predictions_file = os.path.join(data_dir, "predictions.parquet")
    members_file = os.path.join(data_dir, "downloads", "2759_members.tsv.gz")
    output_file = os.path.join(data_dir, "host_orthologs.parquet")

    groups, taxids = get_host_groups(predictions_file)
    print(f"{len(groups)} host orthology groups to look up, {len(taxids)} hosts")

    counts = count_members(members_file, groups, taxids)
    orthologs = pd.DataFrame(counts, columns=["group", "taxid", "n_proteins", "proteins"])
    # what the app asks is whether the host has a protein of the group at all
    orthologs = orthologs[orthologs["n_proteins"] > 0]

    utils.save_to_parquet(orthologs, output_file)
    found = orthologs.groupby("taxid")["group"].nunique()
    for taxid, n in found.items():
        print(f"  {taxid}: proteins in {n} of the {len(groups)} groups")
    print(f"Written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data", help="directory holding predictions.parquet")
    main(parser.parse_args().data_dir)
