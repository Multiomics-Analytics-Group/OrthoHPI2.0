import utils


def get_eggnog_groups(filepath, proteins):
    """
    Obtains all the EggNOG groups which contains a list of given proteins
    :param str filepath: path to the EggNOG groups file
    :param list proteins: list of Ensembl protein identifiers
    :return: dictionary with all the valid EggNOG groups. Key -> group, value -> list proteins in the group
    """
    valid_groups = {}
    groups = utils.read_gzipped_file(filepath)
    first = True
    for line in groups:
        if first:
            first = False
            continue
        data = line.rstrip().split('\t')
        group = data[1]
        gproteins = data[4].split(',')
        int_proteins = list(set(proteins).intersection(gproteins))
        if  len(int_proteins) > 0:
            valid_groups[group] = int_proteins
    
    return valid_groups
            
def get_links(filepath, valid_groups, proteins, ouput_filepath, config_file):
    """
    Obtain the transferred interactions at the EggNOG group level from STRING
    Writes into a file 'predictions.tsv' with the list of predicted links based on homology.
    Structure of the file: taxid1, protein1, taxid2, protein2, str(experimental_evidence), str(databases_evidence), group1, group2
    :param str filepath: path to STRING file with the groups links
    :param dict valid_groups: dictionary with all the valid groups
    :param str output_filepath: path to output file
    """
    links = []
    seen = set()
    cog_links = utils.read_gzipped_file(filepath)
    hosts = utils.read_config(filepath=config_file, field='hosts')
    first = True
    with open(ouput_filepath, 'w') as out:
        out.write("\t".join(["taxid1", "source", "source_name", "taxid2", "target", "target_name", "experimental_evidence_score", "databases_evidence_score", "average_score", "group1", "group2", "edge_type"])+"\n")
        for line in cog_links:
            if first:
                first = False
                continue
            data = line.rstrip().split(' ')
            group1 = data[0]
            group2 = data[1]
            experimental_evidence = int(data[6])/1000
            databases_evidence = int(data[7])/1000
            
            if group1 in valid_groups and group2 in valid_groups:
                if experimental_evidence >= 0.7 or databases_evidence >= 0.7:
                    average_score = (experimental_evidence + databases_evidence) / 2
                    for protein1 in valid_groups[group1]:
                        for protein2 in valid_groups[group2]:
                            taxid1 = protein1.split('.')[0]
                            taxid2 = protein2.split('.')[0]
                            if (int(taxid1) in hosts or int(taxid2) in hosts) and (taxid1 != taxid2):
                                if (protein1, protein2) not in seen:
                                    if int(taxid1) in hosts:
                                        aux_id = taxid1
                                        aux_prot = protein1
                                        aux_group = group1
                                        taxid1 = taxid2
                                        protein1 = protein2
                                        group1 = group2
                                        taxid2 = aux_id
                                        protein2 = aux_prot
                                        group2 = aux_group
                                    out.write("\t".join([taxid1, protein1, proteins[protein1], taxid2, protein2, proteins[protein2], str(experimental_evidence), str(databases_evidence), str(average_score), group1, group2, "inter-species"])+"\n")
                                    seen.add((protein1, protein2))
                                    seen.add((protein2, protein1))
                                        

    
