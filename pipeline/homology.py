import pandas as pd
import utils

# STRING evidence scores are stored as integers 0-1000; divide to get a 0-1 score.
STRING_SCORE_SCALE = 1000
# keep an inter-group link only if its experimental or database evidence is at least this
EVIDENCE_CUTOFF = 0.7

# columns of the predicted-links table written by get_links
LINK_COLUMNS = ["taxid1", "taxid1_label", "source_color", "source_shape", "source", "source_name",
                "taxid2", "taxid2_label", "target_color", "target_shape", "target", "target_name",
                "experimental_evidence_score", "databases_evidence_score", "weight",
                "group1", "group2", "edge_type"]


def get_eggnog_groups(filepath, proteins):
    """
    Obtains all the EggNOG groups which contains a list of given proteins
    :param str filepath: path to the EggNOG groups file
    :param list proteins: list of Ensembl protein identifiers
    :return: dictionary with all the valid EggNOG groups. Key -> group, value -> list proteins in the group
    """
    print("  Scanning EggNOG groups...")
    valid_groups = {}
    protein_set = set(proteins)
    with utils.read_gzipped_file(filepath) as groups:
        for i, line in enumerate(groups, 1):
            if i % 100_000 == 0:
                print(f"    {i:,} groups scanned, {len(valid_groups)} matched")
            data = line.decode("utf-8").rstrip().split('\t')
            group, gproteins = data[1], data[4].split(',')
            matched = protein_set.intersection(gproteins)
            if matched:
                valid_groups[group] = list(matched)

    return valid_groups


def get_links(filepath, valid_groups, proteins, ouput_filepath, config_file):
    """
    Obtain the transferred interactions at the EggNOG group level from STRING.
    Writes the predicted links to output_filepath as parquet (columns: LINK_COLUMNS)
    and returns them as a DataFrame. Each link is a parasite-host protein pair whose
    orthology groups interact in STRING, restricted to hosts the parasite infects.

    :param str filepath: path to STRING file with the groups links
    :param dict valid_groups: dictionary with all the valid groups
    :param dict proteins: mapping from ENSP to protein name
    :param str ouput_filepath: path to output parquet file
    :param str config_file: path to the configuration file
    :return: DataFrame of predicted links (columns: LINK_COLUMNS)
    """
    links = []
    seen = set()
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')

    print("  Scanning COG links...")
    with utils.read_gzipped_file(filepath) as cog_links:
        next(cog_links, None)  # skip header
        for line in cog_links:
            data = line.decode("utf-8").rstrip().split(' ')
            group1, group2 = data[0], data[1]
            if group1 not in valid_groups or group2 not in valid_groups:
                continue
            experimental_evidence = round(int(data[6])/STRING_SCORE_SCALE, 3)
            databases_evidence = round(int(data[7])/STRING_SCORE_SCALE, 3)
            if experimental_evidence < EVIDENCE_CUTOFF and databases_evidence < EVIDENCE_CUTOFF:
                continue
            average_score = round((experimental_evidence + databases_evidence) / 2, 3)

            for protein1 in valid_groups[group1]:
                taxid1 = protein1.split('.')[0]
                for protein2 in valid_groups[group2]:
                    taxid2 = protein2.split('.')[0]
                    taxid1_is_host = int(taxid1) in hosts
                    taxid2_is_host = int(taxid2) in hosts
                    # keep only host-parasite pairs; skip host-host and parasite-parasite
                    if taxid1_is_host == taxid2_is_host or taxid1 == taxid2:
                        continue
                    if (protein1, protein2) in seen:
                        continue
                    if taxid1_is_host:
                        target_taxid, target_group, target_protein = taxid1, group1, protein1
                        source_taxid, source_group, source_protein = taxid2, group2, protein2
                    else:
                        target_taxid, target_group, target_protein = taxid2, group2, protein2
                        source_taxid, source_group, source_protein = taxid1, group1, protein1
                    # restrict to hosts this parasite actually infects
                    allowed_hosts = parasites[int(source_taxid)].get('hosts')
                    if allowed_hosts is not None and int(target_taxid) not in allowed_hosts:
                        continue
                    links.append([source_taxid, parasites[int(source_taxid)]['label'],
                                parasites[int(source_taxid)]['color'], 'diamond', source_protein,
                                proteins[source_protein],
                                target_taxid, hosts[int(target_taxid)]['label'],
                                hosts[int(target_taxid)]['color'], 'dot', target_protein,
                                proteins[target_protein],
                                str(experimental_evidence), str(databases_evidence), str(average_score),
                                source_group, target_group, "inter-species"])
                    seen.add((protein1, protein2))
                    seen.add((protein2, protein1))

    links_df = pd.DataFrame(links, columns=LINK_COLUMNS)
    print(f"  {len(links_df)} total interactions")
    if not links_df.empty:
        for taxid, count in links_df.groupby('taxid1')['source'].count().items():
            label = parasites[int(taxid)]['label']
            print(f"    {label} ({taxid}): {count} interactions")
    utils.save_to_parquet(links_df, ouput_filepath)

    return links_df