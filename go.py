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


if __name__ == "__main__":
    config_file = 'config_short.yml'
    
    get_gene_ontology(config_file, output_dir='data/')