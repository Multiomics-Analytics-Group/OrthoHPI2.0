import os
import homology
import utils
import graph
   

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
            
            data = line.rstrip().split('\t')
            identifier = data[0]
            name = data[1]
            proteins[identifier] = name
            
    return proteins

def get_gene_ontology(config_file, valid_proteins):
    """
    Retrieve gene ontology biological processes for all valid proteins
    :param str config_file: path to config file
    :param dict valid_proteins: dictionary of valid proteins for each specie. Key -> tax id, value -> dictionary: key -> protein id, value -> protein name
    
    :return: dictionary with all biological processes annotated for all proteins for all species. Key -> tax id, value -> dictionary: key -> protein id, value -> go term
    """
    
    gos = {}
    hosts = utils.read_config(filepath=config_file, field='hosts')
    parasites = utils.read_config(filepath=config_file, field='parasites')
    urls = utils.read_config(filepath=config_file, field='urls')
    if "string_go_url" in urls:
        string_file = urls['string_go_url']
        if hosts is not None and parasites is not None:
            taxids = list(hosts.keys()) + list(parasites.keys())
            for taxid in taxids:
                gos[taxid] = parse_gene_ontology(string_file, taxid, valid_proteins[taxid])
    
    return gos


def parse_gene_ontology(string_file, taxid, valid_proteins):
    """
    Retrieve gos for a given specie
    :param str string_file: url to string PPI file
    :param int taxid: taxonomic id of the species of interest
    :param dict valid_proteins: dictionary of valid proteins for the species of interest. key -> protein id, value -> protein name
    
    :return: dictionary with all go terms (description). Key -> Ensembl protein id, value -> go term
    """
    gos = {}
    if string_file is not None:
        filename = utils.download_file(url=string_file.replace('TAXID', str(taxid)), data_dir='data')
        sp = utils.read_gzipped_file(filename)
        first = True
        for line in sp:
            if first:
                first = False
                continue
            
            data = line.rstrip().split('\t')
            identifier = data[0]
            category = data[1]
            term = data[2]
            description = data[3]
            if category == 'Biological Process (Gene Ontology)' and identifier in valid_proteins:
                gos[identifier] = description
            
    return gos

def get_host_ppi(config_file, valid_proteins, ouput_filepath):
    """
    Retrieve Host intra-species proteins
    :param str config_file: path to config file
    :param dict valid_proteins: dictionary with valid proteins
    :param str ouput_filepath: path to edges file
    protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score
    """
    seen = set()
    hosts = utils.read_config(filepath=config_file, field='hosts')
    urls = utils.read_config(filepath=config_file, field='urls')
    with open(ouput_filepath, 'a') as out:
        if "string_ppi_url" in urls:
            string_file = urls['string_ppi_url']
            if hosts is not None:
                for host in hosts:
                    if string_file is not None:
                        filename = utils.download_file(url=string_file.replace('HOST', str(host)), data_dir='data')
                        sp = utils.read_gzipped_file(filename)
                        first = True
                        for line in sp:
                            if first:
                                first = False
                                continue
                            data = line.rstrip().split(' ')
                            protein1 = data[0]
                            protein2 = data[1]
                            exp_score = float(data[6])/1000
                            db_score = float(data[7])/1000
                            avg_score = (exp_score + db_score)/2
                            if protein1 in valid_proteins and protein2 in valid_proteins:
                                if exp_score >= 0.7 or db_score >= 0.7:
                                    if (protein1, protein2) not in seen:
                                        out.write("\t".join([str(host), protein1, str(host), protein2, str(exp_score), str(db_score), str(avg_score), "", "", "intra-species"])+"\n")
                                        seen.add((protein1, protein2))
                                        seen.add((protein2, protein1))

        
def setup(config_file, output_file_path):
    """
    Downloads all necessary files according to the urls specified in configuration
    except the protein files and host PPIs
    
    :param str config_file: path to the configuration file
    """
    urls = utils.read_config(filepath=config_file, field='urls')
    for url_name in urls:
        url = urls[url_name]
        if url_name != "string_protein_url" and url_name != "string_ppi_url" and url_name != "string_go_url":
            filename = utils.download_file(url=url, data_dir=output_file_path)


def apply_tissue_filter(config_file, valid_proteins, cutoff):
    hosts = utils.read_config(filepath=config_file, field='hosts')
    for taxid in hosts:
        proteins = valid_proteins[taxid]
        #print("T before", len(proteins))
        if 'tissues_url' in hosts[taxid]:
            url = hosts[taxid]['tissues_url']
            filename = utils.download_file(url=url, data_dir='data')
            tissues, proteins = get_tissues(filename, proteins, cutoff)
    
        valid_proteins[taxid] = proteins
        #print("T after", len(valid_proteins[taxid]))

    return tissues

def get_tissues(tissues_file, valid_proteins, cutoff):
    """
    Get protein tissue expression for relevant tissues in the lifecycle of the
    studied parasites
    
    :param str tissues_file: path to file with tissue expression (tissues.jensenlab.org)
    :param list valid_proteins: all proteins studied
    :param float cutoff: minimum confidence score accepted (tissues.jensenlab.org)
    
    :return tissues: dictionary with protein tissue expression
    """
    tissues = {}
    filters = {}
    valid_tissues = set()
    parasites = utils.read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        t = parasites[parasite]['tissues']
        valid_tissues.update(t)
    
    first = True
    with open(tissues_file, 'r') as f:
        for line in f:
            if first:
                first = False
                continue
            
            data = line.rstrip().split('\t')
            protein  = "9606."+data[0]
            tissue = data[2]
            score = float(data[6])
            
            if protein in valid_proteins and score > cutoff and tissue in valid_tissues:
                if protein not in tissues:
                    tissues[protein] = []
                
                tissues[protein].append(tissue)
                filters[protein] = valid_proteins[protein]
    
    return tissues, filters


def apply_compartment_filter(config_file, valid_proteins, cutoff):
    hosts = utils.read_config(filepath=config_file, field='hosts')
    for taxid in hosts:
        proteins = valid_proteins[taxid]
        #print("C before", len(proteins))
        if 'compartments_url' in hosts[taxid]:
            url = hosts[taxid]['compartments_url']
            filename = utils.download_file(url=url, data_dir='data')
            compartments, proteins = get_compartments(filename, proteins, cutoff)
    
        valid_proteins[taxid] = proteins
        #print("C after", len(valid_proteins[taxid]))
    
    return compartments


def get_compartments(compartments_file, valid_proteins, cutoff):
    """
    Get protein cellular compartment expression relevant in the lifecycle of the
    studied parasites
    
    :param str compartments_file: path to file with cellular compartment expression (compartments.jensenlab.org)
    :param dict valid_proteins: dictionary with annotations in valid proteins
    :param float cutoff: minimum confidence score accepted (tissues.jensenlab.org)
    
    :return filtered_dict: dictionary with only proteins in relevant compartments
    """
    compartments = {}
    filters = {}
    valid_compartments = set()
    parasites = utils.read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        if "compartments" in parasites[parasite]:
            t = parasites[parasite]['compartments']
        else:
            t = 'GO:0005886'
        valid_compartments.add(t)
    
    first = True
    with open(compartments_file, 'r') as f:
        for line in f:
            if first:
                first = False
                continue
            
            data = line.rstrip().split('\t')
            protein  = "9606."+data[0]
            compartment = data[2]
            score = float(data[4])
            if protein in valid_proteins and score > cutoff and compartment in valid_compartments:
                if protein not in compartments:
                    compartments[protein] = []
                
                compartments[protein].append(compartment)
                filters[protein] = valid_proteins[protein]
  
    return compartments, filters

def get_secretome_predictions(secretome_dir, valid_proteins):
    """
    Filter out proteins that are not secreted or membrane from the list of parasite proteins
    
    :param str secretome_dir: path to the directory where the prediction files are
    :param dict valid_proteins: dictionary with annotations in valid proteins
        
    :return filtered_dict: dictionary with only secreted or membrane parasite proteins
    """
    parasites = utils.read_config(filepath=config_file, field='parasites')
    for parasite in parasites:
        filepath = os.path.join(secretome_dir, str(parasite)+'.fasta')
        sequences = utils.read_fasta(filepath)
        filter_out_ids = utils.filter_sequences(sequences, valid_proteins[parasite])
        for k in filter_out_ids:
            valid_proteins[parasite].pop(k, None)
    
    return valid_proteins


if __name__ == "__main__":
    data_dir = 'data'
    config_file = 'config.yml'
    
    setup(config_file=config_file, output_file_path=data_dir)
    proteins = get_proteins(config_file)
    proteins = get_secretome_predictions(secretome_dir='data/secretome_pred_input_data/input_data', valid_proteins=proteins)
    tissues = apply_tissue_filter(config_file, proteins, cutoff=3.0)
    compartments = apply_compartment_filter(config_file, proteins, cutoff=3.5)
    functions = get_gene_ontology(config_file, valid_proteins=proteins)
    proteins = utils.merge_dict_of_dicts(dict_of_dicts=proteins)
    valid_groups = homology.get_eggnog_groups(filepath=os.path.join(data_dir, '2759_members.tsv.gz'), proteins=proteins.keys())
    
    homology.get_links(filepath=os.path.join(data_dir, 'COG.links.detailed.v11.5.txt.gz'), valid_groups=valid_groups,
              ouput_filepath=os.path.join(data_dir, 'predictions.tsv'), config_file=config_file)
    get_host_ppi(config_file=config_file, valid_proteins=proteins, ouput_filepath=os.path.join(data_dir, 'predictions.tsv'))
    graph.generate_cytoscape_network(edges_file_path=os.path.join(data_dir, 'predictions.tsv'), proteins=proteins,
                                        tissues=tissues, config_file=config_file, output_dir_path='web/public', n=6)
    graph.generate_common_cytoscape_network(edges_file_path=os.path.join(data_dir, 'predictions.tsv'), proteins=proteins,
                                        tissues=tissues, config_file=config_file, output_dir_path='web/public', min_common=15)
    graph.generate_gos_cytoscape_network(edges_file_path=os.path.join(data_dir, 'predictions.tsv'), functions=functions, output_dir_path=data_dir)
    
    
    
