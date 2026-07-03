import time
import pandas as pd
import utils


def get_eggnog_groups(filepath, proteins):
    """
    Obtains all the EggNOG groups which contains a list of given proteins
    :param str filepath: path to the EggNOG groups file
    :param list proteins: list of Ensembl protein identifiers
    :return: dictionary with all the valid EggNOG groups. Key -> group, value -> list proteins in the group
    """
    print("  Scanning EggNOG groups...")
    sum_prots = 0
    valid_groups = {}
    protein_set = set(proteins)
    groups = utils.read_gzipped_file(filepath)
    first = True
    i = 0
    t0 = time.time()
    for line in groups:
        if first:
            first = False
            continue
        i += 1
        if i % 100_000 == 0:
            print(f"    {i:,} groups scanned, {len(valid_groups)} matched ({time.time()-t0:.0f}s)")
        data = line.decode("utf-8").rstrip().split('\t')
        group = data[1]
        gproteins = data[4].split(',')
        int_proteins = list(protein_set.intersection(gproteins))
        if len(int_proteins) > 0:
            valid_groups[group] = int_proteins
            sum_prots += len(int_proteins)

    return valid_groups


def get_links(filepath, valid_groups, proteins, ouput_filepath, config_file):
    """
    Obtain the transferred interactions at the EggNOG group level from STRING
    Writes into a file 'predictions.tsv' with the list of predicted links based on homology.
    Structure of the file: 
    ["taxid1", "taxid1_label", "source_color", "source_shape", "source", "source_name", \
                            "taxid2", "taxid2_label", "target_color", "target_shape", "target", "target_name", \
                            "experimental_evidence_score", "databases_evidence_score", "weight", \
                            "group1", "group2", "edge_type"]

    :param str filepath: path to STRING file with the groups links
    :param dict valid_groups: dictionary with all the valid groups
    :param dict proteins: mapping from ENSP to protein name
    :param str output_filepath: path to output file
    :param str config_file: path to the configuration file
    """
    links = []
    seen = set()
    cog_links = utils.read_gzipped_file(filepath)
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    first = True
    cols = ["taxid1", "taxid1_label", "source_color", "source_shape", "source", "source_name", \
                            "taxid2", "taxid2_label", "target_color", "target_shape", "target", "target_name", \
                            "experimental_evidence_score", "databases_evidence_score", "weight", \
                            "group1", "group2", "edge_type"]
    print("  Counting COG links...")
    total_lines = sum(1 for _ in cog_links) - 1
    cog_links.close()
    print(f"  {total_lines:,} total lines")
    cog_links = utils.read_gzipped_file(filepath)
    first = True
    i = 0
    t0 = time.time()
    print("  Scanning COG links...")
    for line in cog_links:
        if first:
            first = False
            continue
        i += 1
        if i % 1_000_000 == 0:
            print(f"    {i:,} / {total_lines:,} lines scanned, {len(links)} interactions found ({time.time()-t0:.0f}s)")
        data = line.decode("utf-8").rstrip().split(' ')
        group1 = data[0]
        group2 = data[1]
        experimental_evidence = round(int(data[6])/1000, 3)
        databases_evidence = round(int(data[7])/1000, 3)
        
        if group1 in valid_groups and group2 in valid_groups:
            if experimental_evidence >= 0.7 or databases_evidence >= 0.7:
                average_score = (experimental_evidence + databases_evidence) / 2
                average_score = round(average_score, 3)
                for protein1 in valid_groups[group1]:
                    taxid1 = protein1.split('.')[0]
                    for protein2 in valid_groups[group2]:
                        taxid2 = protein2.split('.')[0]
                        if int(taxid1) in hosts or int(taxid2) in hosts:
                            if taxid1 != taxid2:
                                if (protein1, protein2) not in seen:
                                    if int(taxid1) in hosts:
                                        target_taxid = taxid1
                                        target_group = group1
                                        target_protein = protein1
                                        source_taxid = taxid2
                                        source_group = group2
                                        source_protein = protein2
                                    else:
                                        target_taxid = taxid2
                                        target_group = group2
                                        target_protein = protein2
                                        source_taxid = taxid1
                                        source_group = group1
                                        source_protein = protein1
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
    links_df = pd.DataFrame(links, columns=cols)
    print(f"  {len(links_df)} total interactions")
    if not links_df.empty:
        for taxid, count in links_df.groupby('taxid1')['source'].count().items():
            label = parasites[int(taxid)]['label']
            print(f"    {label} ({taxid}): {count} interactions")
    utils.save_to_parquet(links_df, ouput_filepath)