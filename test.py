import homology


def test_config_file(config_file):
    '''Test wether the configuration file exists or not'''
    pass

def test_config_file_yaml_format(config_file):
    '''Test if the config file format is yaml'''
    pass

def test_config_field(config, field):
    '''Test whether the configuration field exists or not'''
    pass

def test_get_eggnog_homologs(url, proteins, taxids, group):
    '''Test if the get eggnog orthologs function returns the right homologs
        expected response from url: {'orthologs': [[['10090', 'Mus musculus', ['10090.ENSMUSP00000026572', '10090.ENSMUSP00000029445', '10090.ENSMUSP00000032399']],
        [['44689', 'Dictyostelium discoideum', ['44689.DDB0214827']], ['5786', 'Dictyostelium purpureum', ['5786.XP_003294405.1']]]]]}
    '''
    expected_paralogs = ['10090', 'Mus musculus', ['10090.ENSMUSP00000026572', '10090.ENSMUSP00000029445', '10090.ENSMUSP00000032399']]
    expected_orthologs = [['44689', 'Dictyostelium discoideum', ['44689.DDB0214827']], ['5786', 'Dictyostelium purpureum', ['5786.XP_003294405.1']]]
    url = "http://eggnogapi5.embl.de/pairwise_orthologs_by_member/json/PROTEINS/TAXIDS/GROUP"
    proteins = ['10090.ENSMUSP00000032399']
    taxids = ['5786', '44689']
    group = 'KOG0395'
    orthologs, paralogs = homology.get_eggnog_orthologs(url=url, proteins=proteins, taxids=taxids, group=group)
    
    assert(expected_paralogs == paralogs)
    assert(expected_orthologs == orthologs)
