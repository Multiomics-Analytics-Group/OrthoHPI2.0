import os
import pandas as pd
import utils


def get_gene_ontology(config_file, output_dir):
    """
    Retrieve gene ontology biological processes for all valid proteins
    :param str config_file: path to config file
    :param dict valid_proteins: dictionary of valid proteins for each specie. Key -> tax id, value -> dictionary: key -> protein id, value -> protein name
    
    :return: dictionary with all biological processes annotated for all proteins for all species. Key -> tax id, value -> dictionary: key -> protein id, value -> go term
    """
    
    gos = []
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')
    if "string_go_url" in urls:
        string_file = urls['string_go_url']
        if hosts is not None and parasites is not None:
            taxids = list(hosts.keys()) + list(parasites.keys())
            for taxid in taxids:
                gos.append(parse_gene_ontology(string_file, taxid))
    
    gos = pd.concat(gos)

    gos.to_csv(os.path.join(output_dir, 'gos.tsv'), sep='\t', header=True, index=False, doublequote=None)



def parse_gene_ontology(string_file, taxid):
    """
    Retrieve gos for a given specie
    :param str string_file: url to string PPI file
    :param int taxid: taxonomic id of the species of interest
    :param dict valid_proteins: dictionary of valid proteins for the species of interest. key -> protein id, value -> protein name
    
    :return: dictionary with all go terms (description). Key -> Ensembl protein id, value -> go term
    """
    data = pd.DataFrame()
    if string_file is not None:
        filename = utils.download_file(url=string_file.replace('TAXID', str(taxid)), data_dir='data')
        data = pd.read_csv(filename, sep='\t', compression='gzip')
        data = data[data['category'] == 'Biological Process (Gene Ontology)']
        data = data[['#string_protein_id', 'description']]
        data['taxid'] = taxid
        
    
    return data

def parse_ontology(config_file, output_directory):
    urls = utils.read_config(filepath=config_file, field='urls')

    terms = {}
    rels = []
    if 'go_ontology_url' in urls:
        filename = utils.download_file(url=urls['go_ontology_url'], data_dir='data')
        graph = utils.convertOBOtoNet(filename)
        for term, attr in graph.nodes(data=True):
            if "name" in attr:
                terms[term] = attr["name"].capitalize()
            if "is_a" in attr:
                for isa in attr["is_a"]:
                    if isa in terms:
                        isa = terms[isa]
                    rels.append([isa, terms[term]])
    
    rels = pd.DataFrame(rels, columns=['parent', 'child'])
    mapped_terms = []
    for i, row in rels.iterrows():
        term = row['parent']
        if term in terms:
            term = terms[term]
        mapped_terms.append(term)
    
    rels['parent'] = mapped_terms
    rels.to_csv(os.path.join(output_directory, 'go_ontology.tsv'), sep='\t', header=True, index=False, doublequote=None)

if __name__ == "__main__":
    config_file = 'config.yml'
    
    #get_gene_ontology(config_file, output_dir='data/')
    parse_ontology(config_file, output_directory='data/')