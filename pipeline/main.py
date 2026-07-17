import os
import pandas as pd
import utils
from . import homology, filters, hpa, go
   

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
        filename = utils.download_file(url=string_file.replace('TAXID', str(taxid)), data_dir=os.path.join('data/downloads/species', str(taxid)))
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


PER_SPECIES_URLS = {"string_protein_url", "string_ppi_url", "string_go_url", "string_alias_url", "string_sequences_url"}


def setup(config_file, output_file_path):
    """
    Downloads all necessary files according to the urls specified in the configuration file
    except the ones templated per species (contain a TAXID placeholder), which are
    downloaded elsewhere once a taxid is known. The go terms will also be downloaded and formatted.

    :param str config_file: path to the configuration file
    """
    urls = utils.read_config(filepath=config_file, field='urls')
    for url_name in urls:
        url = urls[url_name]
        if url_name not in PER_SPECIES_URLS:
            filename = utils.download_file(url=url, data_dir=os.path.join(output_file_path, 'downloads'))
    
    go.get_gene_ontology(config_file, output_dir=output_file_path)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    args = parser.parse_args()
    data_dir = args.data_dir
    downloads_dir = os.path.join(data_dir, 'downloads')
    config_file = args.config

    print("Setup: downloading reference files...")
    setup(config_file=config_file, output_file_path=data_dir)

    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')

    print("Getting proteins...")
    proteins = get_proteins(config_file)
    total_proteins = sum(len(v) for v in proteins.values())
    print(f"  {total_proteins} proteins before filtering")

    print("Applying secretome/tissue/compartment filters...")
    proteins = filters.get_secretome_predictions(config_file=config_file, secretome_dir=os.path.join(data_dir, 'secretome'), valid_proteins=proteins)
    tissues = filters.apply_tissue_filter(config_file, proteins, cutoff=2.5)
    compartments = filters.apply_compartment_filter(config_file, proteins, cutoff=2.5)
    proteins = utils.merge_dict_of_dicts(dict_of_dicts=proteins)
    print(f"  {len(proteins)} proteins after filtering")

    print("Annotating tissue and cell type expression...")
    get_tissue_cell_type_annotation(tissues, output_file=os.path.join(data_dir, 'tissues_cell_types.parquet'))

    print("Getting EggNOG groups and transferring PPIs...")
    cog_filename = urls['string_COG_url'].split('/')[-1]
    valid_groups = homology.get_eggnog_groups(filepath=os.path.join(downloads_dir, '2759_members.tsv.gz'), proteins=proteins.keys())
    print(f"  {len(valid_groups)} valid EggNOG groups")
    from collections import Counter
    taxid_counts = Counter()
    for prots in valid_groups.values():
        for p in prots:
            taxid_counts[p.split('.')[0]] += 1
    for taxid, count in sorted(taxid_counts.items()):
        print(f"    taxid {taxid}: {count} proteins in EggNOG groups")
    homology.get_links(filepath=os.path.join(downloads_dir, cog_filename), valid_groups=valid_groups, proteins=proteins,
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
