import os
import pandas as pd
import utils


def get_go_annotations(config_file, output_dir):
    """
    Retrieve gene ontology biological processes for all valid proteins
    :param str config_file: path to config file

    :return: dictionary with all biological processes annotated for all proteins for all species. Key -> tax id, value -> dictionary: key -> protein id, value -> go term
    """

    go_frames = []
    build_go_hierarchy(config_file=config_file, output_directory=output_dir)
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')
    if "string_go_url" in urls and hosts is not None and parasites is not None:
        string_url = urls['string_go_url']
        for taxid in list(hosts.keys()) + list(parasites.keys()):
            go_frames.append(get_species_go(string_url, taxid))

    gos = pd.concat(go_frames) if go_frames else pd.DataFrame(columns=['#string_protein_id', 'term', 'description', 'taxid'])

    utils.save_to_parquet(gos, os.path.join(output_dir, 'gos.parquet'))


def get_species_go(string_url, taxid):
    """
    Retrieve gos for a given species
    :param str string_url: url template to the STRING GO enrichment terms file (contains TAXID)
    :param int taxid: taxonomic id of the species of interest

    :return: dataframe with the biological-process GO terms per protein for this species,
             identified by GO id with the description kept for reading
    """
    if string_url is None:
        return pd.DataFrame()

    filename = utils.download_file(url=string_url.replace('TAXID', str(taxid)), data_dir=os.path.join('data/downloads/species', str(taxid)))
    data = pd.read_csv(filename, sep='\t', compression='gzip')
    data = data[data['category'] == 'Biological Process (Gene Ontology)']
    # the GO id is kept as the identity of a term: descriptions are what a person reads,
    # but they do not survive a round trip through the ontology intact
    data = data[['#string_protein_id', 'term', 'description']]
    data['taxid'] = taxid

    return data


def build_go_hierarchy(config_file, output_directory):
    """
    The parent-child relations of the ontology, as GO ids.

    Terms are related by id rather than by name: a name is not a stable key -- STRING and
    the ontology punctuate the same term differently -- and matching on it silently drops
    the relations of the terms that do not match. Both `is_a` and `part_of` are followed,
    since a biological process is as often a part of a larger process as a kind of one.
    """
    urls = utils.read_config(filepath=config_file, field='urls')
    if 'go_ontology_url' not in urls:
        return

    filename = utils.download_file(url=urls['go_ontology_url'], data_dir='data/downloads')
    graph = utils.convertOBOtoNet(filename)

    # parent -> child relations, over the two relations that make a term narrower than
    # the one above it
    relations = []
    for term, attr in graph.nodes(data=True):
        parents = list(attr.get('is_a', []))
        parents += [rel.split(' ', 1)[1] for rel in attr.get('relationship', [])
                    if rel.startswith('part_of ')]
        for parent in parents:
            relations.append([parent, term])

    relations = pd.DataFrame(relations, columns=['parent', 'child']).drop_duplicates()

    utils.save_to_parquet(relations, os.path.join(output_directory, 'go_ontology.parquet'))
