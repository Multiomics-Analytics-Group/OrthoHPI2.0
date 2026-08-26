"""
Which orthology groups each host proteome has a protein in, for the host groups that
carry a predicted interaction.

The cross-host page has to tell two different absences apart. A parasite predicted to
interact with a human protein but not with its pig counterpart is either a parasite
that cannot reach that protein in pig -- because pig has no protein in the group at
all -- or the same interaction filtered out of the pig predictions, because no pig
protein of the group is annotated to a tissue the parasite infects, or DeepLoc did not
call one surface-exposed, or STRING has no pig entry for it. Only the first is about
the host's biology; the second is about how deeply the host is annotated, and reading
it as the first is the mistake this file exists to prevent.

predictions.parquet cannot answer that: a group with no predicted interaction in pig is
absent from it whichever of the two happened. The EggNOG members file can, since it
lists every protein of every group.

Only the host groups that appear in predictions.parquet are kept (193 of them at the
time of writing), so the output is a few kilobytes rather than a scan of all of
Eukaryota.

    python scripts/build_host_orthologs.py

Re-run it whenever predictions.parquet is rebuilt with new hosts or new parasites; the
app treats a missing file as "not known" and simply leaves the section out.
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils

# columns of the EggNOG 5-layout members file written by pipeline/prepare_eggnog_members.py
GROUP_COLUMN, PROTEINS_COLUMN, SPECIES_COLUMN = 1, 4, 5


def get_host_groups(predictions_file):
    """
    The orthology groups of the host side of the predictions, and the hosts to look for
    in them. Both come from the predictions rather than from the config so that the
    output covers exactly what the app can ask about.

    :param str predictions_file: path to predictions.parquet
    :return: (set of group ids, set of host taxids as str)
    """
    predictions = utils.read_parquet_file(input_file=predictions_file)

    return set(predictions['group2']), set(predictions['taxid2'].astype(str))


def count_members(members_file, groups, taxids):
    """
    Streams the members file and counts, for each wanted group, how many proteins of each
    host it holds. The proteins of a group are one comma-separated field, so the species
    field is read first and the proteins are only split for a group that has a host in it.

    :param str members_file: path to the gzipped EggNOG members file
    :param set groups: group ids to keep
    :param set taxids: host taxids to count, as strings
    :return: list of [group, taxid, n_proteins, proteins] for every host present in a
             wanted group, proteins comma-separated as they are in the members file
    """
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
                # the proteins are STRING-style taxid.identifier, so the host's are the
                # ones under its own prefix
                prefix = f"{taxid}."
                members_of_host = [p for p in proteins if p.startswith(prefix)]
                # the ids are kept, and not only counted, so the app can ask the next
                # question of them: whether any of the host's proteins of the group is
                # annotated to a tissue the parasite infects
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
    # a group whose only member of a host is in EggNOG but not in STRING is still a
    # protein that host has, so the zero counts are kept out rather than kept as zeros:
    # what the app asks is whether the host has one at all
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
