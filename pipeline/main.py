import os
from collections import Counter

import pandas as pd

import utils
from . import cell_type_annotations, homology, filters, go

# jensenlab confidence score below which tissue evidence is ignored;
# hosts.<taxid>.tissue_cutoff overrides it
TISSUE_CUTOFF = 2.5

# DeepLoc 2 Accurate-model per-class thresholds (DeepLoc2/deeploc2.py label_threshold, read
# at i+1); see docs/deeploc.md
DEEPLOC_ACCURATE_DIR = os.path.join('deeploc', 'output_accurate', 'deeploc_output_accurate')
DEEPLOC_CUTOFFS = {
    'Extracellular': 0.61728516,
    'Cell membrane': 0.56464844,
    'Cytoplasm': 0.47612305,
    'Nucleus': 0.50136719,
}

# which classes each niche reaches. Both keep the host surface; intracellular adds cytosol
# and nucleus. The organelle classes are left out: adding them keeps 93% of the human
# proteome
DEEPLOC_NICHE_CLASSES = {
    'extracellular': ('Extracellular', 'Cell membrane'),
    'intracellular': ('Extracellular', 'Cell membrane', 'Cytoplasm', 'Nucleus'),
}
# the narrower niche, so a missing value cannot widen a host pool by accident
DEEPLOC_DEFAULT_NICHE = 'extracellular'
DEEPLOC_NICHE_CUTOFFS = {niche: {c: DEEPLOC_CUTOFFS[c] for c in classes}
                         for niche, classes in DEEPLOC_NICHE_CLASSES.items()}


def get_proteins(config_file):
    '''Retrieve all proteins for all species'''
    proteins = {}
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')
    if "string_protein_url" not in urls or hosts is None or parasites is None:
        return proteins

    string_url = urls['string_protein_url']
    for taxid in list(hosts.keys()) + list(parasites.keys()):
        proteins[taxid] = get_species_proteins(string_url, taxid)

    return proteins


def get_species_proteins(string_url, taxid):
    '''Download and parse the STRING protein.info file for a single species.'''
    proteins = {}
    if string_url is None:
        return proteins

    filename = utils.download_file(url=string_url.replace('TAXID', str(taxid)), data_dir=os.path.join('data/downloads/species', str(taxid)))
    with utils.read_gzipped_file(filename) as handle:
        next(handle, None)  # skip header
        for line in handle:
            identifier, name = line.decode("utf-8").rstrip().split('\t')[:2]
            proteins[identifier] = name

    return proteins


def filter_proteins(config_file, data_dir, proteins):
    '''
    Narrow every species to the proteins an interaction could be predicted between: the
    parasite proteins the secretome predictions call secreted or membrane-bound, and the
    host proteins expressed in a tissue some parasite of the config infects and put by
    DeepLoc in a localisation a parasite infecting that host can reach.
    '''
    # proteins stays {taxid: {protein: name}} through the filters, then is flattened for the
    # homology transfer
    proteins = filters.get_secretome_predictions(config_file=config_file, secretome_dir=os.path.join(data_dir, 'secretome'), valid_proteins=proteins)
    tissues = filters.apply_tissue_filter(config_file=config_file, valid_proteins=proteins, cutoff=TISSUE_CUTOFF)
    # the per-niche halves of the pool go on to get_links, which applies them parasite by
    # parasite
    reachable = filters.apply_deeploc_filter(config_file=config_file, valid_proteins=proteins,
                                             deeploc_dir=os.path.join(data_dir, DEEPLOC_ACCURATE_DIR),
                                             niche_cutoffs=DEEPLOC_NICHE_CUTOFFS,
                                             default_niche=DEEPLOC_DEFAULT_NICHE)
    infected = filters.parasite_tissue_proteins(config_file=config_file, tissues=tissues,
                                                valid_proteins=proteins)

    return utils.merge_dict_of_dicts(dict_of_dicts=proteins), tissues, reachable, infected


def save_eligible_proteins(proteins, output_file, reachable=None):
    '''Write the proteins the filters passed, over every species of the config.'''
    eligible = pd.DataFrame(sorted(proteins.items()), columns=['protein', 'name'])
    eligible['taxid'] = eligible['protein'].str.split('.').str[0]
    if reachable:
        # every host protein of the pool is in some niche's set, so what is in none is a
        # parasite protein
        parasite_side = ~eligible['protein'].isin(set().union(*reachable.values()))
        for niche, host_proteins in sorted(reachable.items()):
            eligible[f'reachable_{niche}'] = (eligible['protein'].isin(host_proteins)
                                              | parasite_side)

    utils.save_to_parquet(eligible, output_file)


def get_tissue_cell_type_annotation(tissues, proteins, config_file, output_file):
    '''Build the (Gene, Tissue, cell-type) annotation table and write it to parquet.'''
    tissues_df = pd.DataFrame(
        [(gene, tissue) for gene, ts in tissues.items() for tissue in ts],
        columns=['Gene', 'Tissue'],
    )
    tissues_df = tissues_df[tissues_df['Gene'].isin(proteins.keys())]
    # HPA annotates human; an optional preprocessed pig atlas adds pig cell types
    cell_type_data = cell_type_annotations.parse_cell_type_data(
        config_file, data_dir=os.path.dirname(output_file), valid_proteins=proteins.keys())
    tissues_df = pd.merge(tissues_df, cell_type_data, on=['Gene', 'Tissue'], how='left')

    utils.save_to_parquet(tissues_df, output_file)


PER_SPECIES_URLS = {"string_protein_url", "string_ppi_url", "string_go_url", "string_alias_url", "string_sequences_url"}

# the EggNOG 6 source file is ~10 GB; pipeline/prepare_eggnog_members.py filters it to level
# 2759
PREPROCESSED_URLS = {"eggNOG_members_url"}


def setup(config_file, output_file_path):
    '''
    Downloads all necessary files according to the urls specified in the configuration
    file except the ones templated per species (contain a TAXID placeholder), which are
    downloaded elsewhere once a taxid is known, and the ones prepared by a separate
    script.
    '''
    urls = utils.read_config(filepath=config_file, field='urls')
    for url_name in urls:
        url = urls[url_name]
        if url_name not in PER_SPECIES_URLS and url_name not in PREPROCESSED_URLS:
            utils.download_file(url=url, data_dir=os.path.join(output_file_path, 'downloads'))
    
    go.get_go_annotations(config_file, output_dir=output_file_path)


def print_group_counts(valid_groups):
    '''Print how many proteins each species contributes to the matched EggNOG groups.'''
    taxid_counts = Counter()
    for prots in valid_groups.values():
        for p in prots:
            taxid_counts[p.split('.')[0]] += 1
    for taxid, count in sorted(taxid_counts.items()):
        print(f"    taxid {taxid}: {count} proteins in EggNOG groups")


def annotate_predictions(predictions, hosts, parasites, config_file):
    '''Add source_uniprot / target_uniprot columns by mapping STRING ids to UniProt.'''
    predictions = utils.annotate_alias_id(predictions_df=predictions,
                            taxids=list(parasites.keys()), config_file=config_file,
                            sources=['Uniprot'], new_col="source_uniprot",
                            mapping_col="source")
    predictions = utils.annotate_alias_id(predictions_df=predictions,
                            taxids=list(hosts.keys()), config_file=config_file,
                            sources=['Ensembl_HGNC_uniprot_ids', 'UniProt_AC'],
                            new_col="target_uniprot", mapping_col="target")
    return predictions


def run(config_file, data_dir, verbose=False):
    '''Run the full prediction pipeline and write the parquet outputs into data_dir.'''
    downloads_dir = os.path.join(data_dir, 'downloads')

    print("Setup: downloading reference files...")
    setup(config_file=config_file, output_file_path=data_dir)

    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')

    print("Getting proteins...")
    proteins = get_proteins(config_file)
    total_proteins = sum(len(v) for v in proteins.values())
    print(f"  {total_proteins} proteins before filtering")

    print("Applying secretome/tissue/DeepLoc filters...")
    proteins, tissues, reachable, infected = filter_proteins(config_file=config_file,
                                                            data_dir=data_dir, proteins=proteins)
    print(f"  {len(proteins)} proteins after filtering")

    print("Writing the proteins the filters passed...")
    save_eligible_proteins(proteins=proteins, reachable=reachable,
                           output_file=os.path.join(data_dir, 'eligible_proteins.parquet'))

    print("Annotating tissue and cell type expression...")
    get_tissue_cell_type_annotation(tissues=tissues, proteins=proteins, config_file=config_file, output_file=os.path.join(data_dir, 'tissues_cell_types.parquet'))

    print("Getting EggNOG groups and transferring PPIs...")
    # setup() saved the COG links file under its URL basename
    cog_filename = urls['string_COG_url'].split('/')[-1]
    members_file = os.path.join(downloads_dir, '2759_members.tsv.gz')
    if not os.path.isfile(members_file):
        raise SystemExit(f"{members_file} not found — run 'python -m pipeline.prepare_eggnog_members' first")
    valid_groups = homology.get_eggnog_groups(filepath=members_file, proteins=proteins.keys())
    print(f"  {len(valid_groups)} valid EggNOG groups")
    if verbose:
        print_group_counts(valid_groups)
    predictions = homology.get_links(filepath=os.path.join(downloads_dir, cog_filename), valid_groups=valid_groups,
              proteins=proteins, config_file=config_file, reachable=reachable,
              default_niche=DEEPLOC_DEFAULT_NICHE, infected=infected)

    print("Annotating predictions with UniProt accessions...")
    predictions = annotate_predictions(predictions=predictions, hosts=hosts, parasites=parasites, config_file=config_file)
    utils.save_to_parquet(df=predictions, output_file=os.path.join(data_dir, 'predictions.parquet'))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--verbose', action='store_true', help='print per-taxid EggNOG protein counts')
    args = parser.parse_args()

    run(config_file=args.config, data_dir=args.data_dir, verbose=args.verbose)
