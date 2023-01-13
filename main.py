import os
import homology
import utils
import filters
import hpa
import go
import pandas as pd
   

def get_proteins(config_file):
    """
    Retrieve all proteins for all species
    :param str config_file: path to config file
    
    :return: dictionary with all proteins for all species. Key -> tax id, value -> dictionary: key -> protein id, value -> protein name"""
    proteins = {}
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')
    if "string_protein_url" in urls:
        string_file = urls['string_protein_url']
        if hosts is not None and parasites is not None:
            taxids = list(hosts.keys()) + list(parasites.keys())
            for taxid in taxids:
                proteins[taxid] = parse_proteins(string_file, taxid)
    
    return proteins


def parse_proteins(string_file, taxid):
    """
    Retrieve proteins for a given specie
    :param str string_file: url to string PPI file
    :param int taxid: taxonomic id of the species of interest
    
    :return: dictionary with all proteins. Key -> Ensembl protein id, value -> protein name
    """
    proteins = {}
    if string_file is not None:
        filename = utils.download_file(url=string_file.replace('TAXID', str(taxid)), data_dir='data')
        sp = utils.read_gzipped_file(filename)
        first = True
        for line in sp:
            if first:
                first = False
                continue
            
            data = line.decode("utf-8").rstrip().split('\t')
            identifier = data[0]
            name = data[1]
            proteins[identifier] = name
            
    return proteins


def get_tissue_cell_type_annotation(tissues, output_file):
    tissues_df = pd.concat({k: pd.Series(v) for k, v in tissues.items()}).reset_index()
    tissues_df = tissues_df.iloc[:, [0, 2]]
    tissues_df.columns = ['Gene', 'Tissue']
    tissues_df = tissues_df[tissues_df['Gene'].isin(proteins.keys())]
    hpa_data = hpa.parse_hpa(config_file, valid_proteins=proteins.keys())
    tissues_df = pd.merge(tissues_df, hpa_data, on=['Gene', 'Tissue'], how='left')
    
    utils.save_to_parquet(tissues_df, output_file)


def setup(config_file, output_file_path):
    """
    Downloads all necessary files according to the urls specified in the configuration file
    except the protein files and host PPIs that are downloaded for each species. The go terms
    will also be downloaded and formatted.
    
    :param str config_file: path to the configuration file
    """
    urls = utils.read_config(filepath=config_file, field='urls')
    for url_name in urls:
        url = urls[url_name]
        if url_name != "string_protein_url" and url_name != "string_ppi_url" and url_name != "string_go_url":
            filename = utils.download_file(url=url, data_dir=output_file_path)
    
    go.get_gene_ontology(config_file, output_dir=output_file_path)


if __name__ == "__main__":
    data_dir = 'data'
    config_file = 'config.yml'
    
    setup(config_file=config_file, output_file_path=data_dir)

    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    
    #Get host and parasite proteins
    proteins = get_proteins(config_file)

    #Apply filters -- secretome, tissue, cellular compartment context
    proteins = filters.get_secretome_predictions(config_file=config_file, secretome_dir='data/secretome_pred_input_data/input_data', valid_proteins=proteins)
    tissues = filters.apply_tissue_filter(config_file, proteins, cutoff=2.5)
    compartments = filters.apply_compartment_filter(config_file, proteins, cutoff=2.5)
    proteins = utils.merge_dict_of_dicts(dict_of_dicts=proteins)
    
    #Annotate tissue and cell type expression
    get_tissue_cell_type_annotation(tissues, output_file=os.path.join(data_dir, 'tissues_cell_types.parquet'))
    
    #Get eggnog groups and transfer PPIs
    valid_groups = homology.get_eggnog_groups(filepath=os.path.join(data_dir, '2759_members.tsv.gz'), proteins=proteins.keys())
    homology.get_links(filepath=os.path.join(data_dir, 'COG.links.detailed.v11.5.txt.gz'), valid_groups=valid_groups, proteins=proteins,
              ouput_filepath=os.path.join(data_dir, 'predictions.parquet'), config_file=config_file)

    predictions = pd.read_parquet(os.path.join(data_dir, 'predictions.parquet'))
    predictions = utils.annotate_alias_id(predictions_df=predictions, 
                            taxids=list(parasites.keys()), config_file=config_file, 
                            sources=['BLAST_UniProt_AC'], new_col="source_uniprot", 
                            mapping_col="source")
    
    predictions = utils.annotate_alias_id(predictions_df=predictions, 
                            taxids=list(hosts.keys()), config_file=config_file, 
                            sources=['Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)'], 
                            new_col="target_uniprot", mapping_col="target")
    
    utils.save_to_parquet(df=predictions, output_file=os.path.join(data_dir, 'annotated_predictions.parquet'))
