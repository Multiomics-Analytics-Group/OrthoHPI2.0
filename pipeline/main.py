import os
from collections import Counter

import pandas as pd

import utils
from . import homology, filters, hpa, go

# Default jensenlab confidence score below which a host protein's tissue evidence is
# ignored. A host can override it with hosts.<taxid>.tissue_cutoff in config.yml.
TISSUE_CUTOFF = 2.5

# DeepLoc 2 (Accurate) probability cut-offs for keeping a host protein as surface-exposed.
# Values are DeepLoc's own per-class thresholds for the Accurate (ProtT5) model
# (DeepLoc2/deeploc2.py label_threshold, offset by one: labels[i] -> threshold[i+1]).
# A host protein is kept if it is at the cell membrane or secreted; set
# DEEPLOC_EXTRACELLULAR_CUTOFF to None to keep only Cell membrane proteins.
DEEPLOC_ACCURATE_DIR = os.path.join('deeploc', 'output_accurate', 'deeploc_output_accurate')
DEEPLOC_EXTRACELLULAR_CUTOFF = 0.61728516  # None to disable
DEEPLOC_MEMBRANE_CUTOFF = 0.56464844


def get_proteins(config_file):
    """
    Retrieve all proteins for all species
    :param str config_file: path to config file
    
    :return: dictionary with all proteins for all species. Key -> tax id, value -> dictionary: key -> protein id, value -> protein name"""
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
    """
    Download and parse the STRING protein.info file for a single species.
    :param str string_url: url template to the STRING protein.info file (contains TAXID)
    :param int taxid: taxonomic id of the species of interest

    :return: dictionary with all proteins. Key -> Ensembl protein id, value -> protein name
    """
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
    """
    Narrow every species to the proteins an interaction could be predicted between: the
    parasite proteins the secretome predictions call secreted or membrane-bound, and the
    host proteins expressed in a tissue some parasite of the config infects and called
    surface-exposed by DeepLoc.

    Kept apart from run() so that the pool can be rebuilt for a data directory without
    repeating the orthology transfer, which is the expensive half of the pipeline.

    :param str config_file: path to the configuration file
    :param str data_dir: directory holding the secretome and DeepLoc inputs
    :param dict proteins: {taxid: {protein: name}} before filtering; filtered in place
    :return: ({protein: name} over every species, {protein: [tissue, ...]} for the hosts)
    """
    # proteins stays a {taxid: {protein: name}} dict through all three filters,
    # then is flattened to one {protein: name} dict for the homology transfer
    proteins = filters.get_secretome_predictions(config_file=config_file, secretome_dir=os.path.join(data_dir, 'secretome'), valid_proteins=proteins)
    tissues = filters.apply_tissue_filter(config_file=config_file, valid_proteins=proteins, cutoff=TISSUE_CUTOFF)
    # host proteins kept if DeepLoc calls them Extracellular or Cell membrane (surface-exposed)
    filters.apply_deeploc_filter(config_file=config_file, valid_proteins=proteins,
                                 deeploc_dir=os.path.join(data_dir, DEEPLOC_ACCURATE_DIR),
                                 extracellular_cutoff=DEEPLOC_EXTRACELLULAR_CUTOFF,
                                 membrane_cutoff=DEEPLOC_MEMBRANE_CUTOFF)

    return utils.merge_dict_of_dicts(dict_of_dicts=proteins), tissues


def save_eligible_proteins(proteins, output_file):
    """
    Write the proteins the filters passed, over every species of the config.

    This is what the app tests a network's enrichment against. A network holds the
    proteins it holds because these are the ones it could have been built from, so the
    whole proteome is the wrong background: tested against it, a network returns the
    filters that made it -- surface processes for the hosts, secretion for the parasites
    -- whichever parasite is being asked about. The parasites need this file, having no
    tissue table to be read out of.

    :param dict proteins: {protein: name} after the filters, over every species
    :param str output_file: parquet path to write
    """
    eligible = pd.DataFrame(sorted(proteins.items()), columns=['protein', 'name'])
    eligible['taxid'] = eligible['protein'].str.split('.').str[0]

    utils.save_to_parquet(eligible, output_file)


def get_tissue_cell_type_annotation(tissues, proteins, config_file, output_file):
    """
    Build the (Gene, Tissue, cell-type) annotation table and write it to parquet.
    :param dict tissues: {protein_id: [tissue, ...]} from the tissue filter
    :param dict proteins: valid {protein_id: name} after all filters
    :param str config_file: path to the configuration file
    :param str output_file: parquet path to write
    """
    tissues_df = pd.DataFrame(
        [(gene, tissue) for gene, ts in tissues.items() for tissue in ts],
        columns=['Gene', 'Tissue'],
    )
    tissues_df = tissues_df[tissues_df['Gene'].isin(proteins.keys())]
    # HPA annotates human; an optional preprocessed pig atlas adds pig cell types.
    hpa_data = hpa.parse_cell_type_data(config_file, data_dir=os.path.dirname(output_file),
                                        valid_proteins=proteins.keys())
    tissues_df = pd.merge(tissues_df, hpa_data, on=['Gene', 'Tissue'], how='left')

    utils.save_to_parquet(tissues_df, output_file)


PER_SPECIES_URLS = {"string_protein_url", "string_ppi_url", "string_go_url", "string_alias_url", "string_sequences_url"}

# The EggNOG 6 source file is ~10 GB and covers every tax level; it is streamed and
# filtered down to level 2759 by pipeline/prepare_eggnog_members.py instead.
PREPROCESSED_URLS = {"eggNOG_members_url"}


def setup(config_file, output_file_path):
    """
    Downloads all necessary files according to the urls specified in the configuration file
    except the ones templated per species (contain a TAXID placeholder), which are
    downloaded elsewhere once a taxid is known, and the ones prepared by a separate
    script. The go terms will also be downloaded and formatted.

    :param str config_file: path to the configuration file
    """
    urls = utils.read_config(filepath=config_file, field='urls')
    for url_name in urls:
        url = urls[url_name]
        if url_name not in PER_SPECIES_URLS and url_name not in PREPROCESSED_URLS:
            utils.download_file(url=url, data_dir=os.path.join(output_file_path, 'downloads'))
    
    go.get_go_annotations(config_file, output_dir=output_file_path)


def print_group_counts(valid_groups):
    """Print how many proteins each species contributes to the matched EggNOG groups."""
    taxid_counts = Counter()
    for prots in valid_groups.values():
        for p in prots:
            taxid_counts[p.split('.')[0]] += 1
    for taxid, count in sorted(taxid_counts.items()):
        print(f"    taxid {taxid}: {count} proteins in EggNOG groups")


def annotate_predictions(predictions, hosts, parasites, config_file):
    """
    Add source_uniprot / target_uniprot columns by mapping STRING ids to UniProt.

    Parasites and hosts use different STRING alias source names in v12 (the v11.5
    names BLAST_UniProt_AC / Ensembl_HGNC_UniProt_ID(supplied_by_UniProt) no longer
    exist).

    Host sources are a preference list, because the alias files are not uniform:
    Ensembl_HGNC_uniprot_ids holds one canonical accession per protein but exists
    only for human, so rat, mouse and pig fall back to UniProt_AC. Ensembl_UniProt
    is deliberately not used -- it mixes gene names in with the accessions.
    """
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
    """
    Run the full prediction pipeline and write the parquet outputs into data_dir.

    :param str config_file: path to the configuration file
    :param str data_dir: directory for downloads and output parquet files
    :param bool verbose: print the per-taxid EggNOG protein counts
    """
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
    proteins, tissues = filter_proteins(config_file=config_file, data_dir=data_dir,
                                        proteins=proteins)
    print(f"  {len(proteins)} proteins after filtering")

    print("Writing the proteins the filters passed...")
    save_eligible_proteins(proteins=proteins,
                           output_file=os.path.join(data_dir, 'eligible_proteins.parquet'))

    print("Annotating tissue and cell type expression...")
    get_tissue_cell_type_annotation(tissues=tissues, proteins=proteins, config_file=config_file, output_file=os.path.join(data_dir, 'tissues_cell_types.parquet'))

    print("Getting EggNOG groups and transferring PPIs...")
    # setup() saved the COG links file under its URL basename; rebuild that name to find it
    cog_filename = urls['string_COG_url'].split('/')[-1]
    members_file = os.path.join(downloads_dir, '2759_members.tsv.gz')
    if not os.path.isfile(members_file):
        raise SystemExit(f"{members_file} not found — run 'python -m pipeline.prepare_eggnog_members' first")
    valid_groups = homology.get_eggnog_groups(filepath=members_file, proteins=proteins.keys())
    print(f"  {len(valid_groups)} valid EggNOG groups")
    if verbose:
        print_group_counts(valid_groups)
    predictions = homology.get_links(filepath=os.path.join(downloads_dir, cog_filename), valid_groups=valid_groups,
              proteins=proteins, config_file=config_file)

    print("Annotating predictions with UniProt accessions...")
    predictions = annotate_predictions(predictions=predictions, hosts=hosts, parasites=parasites, config_file=config_file)
    # single output: predictions plus source_uniprot / target_uniprot columns
    utils.save_to_parquet(df=predictions, output_file=os.path.join(data_dir, 'predictions.parquet'))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--data-dir', default='data')
    parser.add_argument('--verbose', action='store_true', help='print per-taxid EggNOG protein counts')
    args = parser.parse_args()

    run(config_file=args.config, data_dir=args.data_dir, verbose=args.verbose)
