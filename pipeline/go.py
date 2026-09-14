import os
import pandas as pd
import utils


def get_go_annotations(config_file, output_dir):
    '''Retrieve gene ontology biological processes for all valid proteins'''

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
    '''Retrieve gos for a given species'''
    if string_url is None:
        return pd.DataFrame()

    filename = utils.download_file(url=string_url.replace('TAXID', str(taxid)), data_dir=os.path.join('data/downloads/species', str(taxid)))
    data = pd.read_csv(filename, sep='\t', compression='gzip')
    data = data[data['category'] == 'Biological Process (Gene Ontology)']
    # the GO id is the identity of a term; descriptions do not survive the ontology intact
    data = data[['#string_protein_id', 'term', 'description']]
    data['taxid'] = taxid

    return data


def build_go_hierarchy(config_file, output_directory):
    '''The parent-child relations of the ontology, as GO ids.'''
    urls = utils.read_config(filepath=config_file, field='urls')
    if 'go_ontology_url' not in urls:
        return

    filename = utils.download_file(url=urls['go_ontology_url'], data_dir='data/downloads')
    graph = utils.convertOBOtoNet(filename)

    # parent -> child relations, over the two relations that make a term narrower
    relations = []
    for term, attr in graph.nodes(data=True):
        parents = list(attr.get('is_a', []))
        parents += [rel.split(' ', 1)[1] for rel in attr.get('relationship', [])
                    if rel.startswith('part_of ')]
        for parent in parents:
            relations.append([parent, term])

    relations = pd.DataFrame(relations, columns=['parent', 'child']).drop_duplicates()

    utils.save_to_parquet(relations, os.path.join(output_directory, 'go_ontology.parquet'))
