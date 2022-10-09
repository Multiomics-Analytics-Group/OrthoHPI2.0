import os
import utils


def apply_tissue_filter(config_file, valid_proteins, cutoff):
    hosts = utils.read_config(filepath=config_file, field='hosts')
    tissue_mapping = utils.read_config(filepath=config_file, field='tissues')
    for taxid in hosts:
        proteins = valid_proteins[taxid]
        if 'tissues_url' in hosts[taxid]:
            url = hosts[taxid]['tissues_url']
            filename = utils.download_file(url=url, data_dir='data')
            tissues, proteins = get_tissues(config_file, filename, proteins, cutoff, tissue_mapping)
    
        valid_proteins[taxid] = proteins

    return tissues

def get_tissues(config_file, tissues_file, valid_proteins, cutoff, mapping):
    """
    Get protein tissue expression for relevant tissues in the lifecycle of the
    studied parasites
    
    :param str config_file: path to the configuration file
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
            
            if protein in valid_proteins and score >= cutoff and tissue in valid_tissues:
                if protein not in tissues:
                    tissues[protein] = []
                
                tissues[protein].append(mapping[tissue])
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
            compartments, proteins = get_compartments(config_file, filename, proteins, cutoff)
    
        valid_proteins[taxid] = proteins
        #print("C after", len(valid_proteins[taxid]))
    
    return compartments


def get_compartments(config_file, compartments_file, valid_proteins, cutoff):
    """
    Get protein cellular compartment expression relevant in the lifecycle of the
    studied parasites
    
    :param str config_file: path to the configuration file
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
            if protein in valid_proteins and score >= cutoff and compartment in valid_compartments:
                if protein not in compartments:
                    compartments[protein] = []
                
                compartments[protein].append(compartment)
                filters[protein] = valid_proteins[protein]
  
    return compartments, filters

def get_secretome_predictions(config_file, secretome_dir, valid_proteins):
    """
    Filter out proteins that are not secreted or membrane from the list of parasite proteins
    
    :param str config_file: path to the configuration file
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
