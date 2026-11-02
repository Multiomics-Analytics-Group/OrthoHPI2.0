'''Convert pooled pig-atlas counts into the app's cell-type nCPM Parquet.'''
import argparse

import pandas as pd

import utils

PIG_TAXID = '9823'
COUNTS_PER_MILLION = 1_000_000
# atlas tissue -> config['tissues'] label; other atlas tissues are kept only when their
# lowercased name is already a label
BRAIN_REGIONS = ['Brain', 'Cerebellum', 'FrontalLobe', 'ParietalLobe', 'OccipitalLobe',
                 'TemporalLobe', 'Hypothalamus', 'OVoLT', 'AreaPostrema', 'SubfornicalOrgan']
NEONATAL_ILEUM = [f'Neonatal Ileum_D{day}' for day in ['0', '01', '03', '07', '14', '21']]
TISSUE_MAPPING = {
    **{region: 'brain' for region in BRAIN_REGIONS},
    **{sample: 'intestine' for sample in NEONATAL_ILEUM},
    'Intestine': 'intestine',
    'PBMC': 'blood',
    'Retina': 'eye',
}


def map_expression(input_file, config_file):
    '''Pool the counts of each (tissue, cell type) across samples and scale to a million.'''
    data = pd.read_csv(input_file)
    valid_tissues = {name.lower() for name in
                     utils.read_config(filepath=config_file, field='tissues').values()}
    data['Tissue'] = data['Tissue'].map(lambda t: TISSUE_MAPPING.get(t, t.lower()))
    data = data[data['Tissue'].isin(valid_tissues)]

    data = data.groupby(['Gene', 'Tissue', 'Cell type'], as_index=False)['Counts'].sum()
    totals = data.groupby(['Tissue', 'Cell type'])['Counts'].transform('sum')
    data['nCPM'] = data['Counts'] * COUNTS_PER_MILLION / totals

    aliases = utils.parse_string_aliases(
        config_file=config_file, sources=['Ensembl_gene', 'UniProt_GN_Name'], taxid=PIG_TAXID)
    data['Gene'] = data['Gene'].map(aliases)
    data = data[data['Gene'].notna() & (data['nCPM'] > 0)]
    return data.groupby(['Gene', 'Tissue', 'Cell type'], as_index=False)['nCPM'].mean()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', required=True,
                        help='CSV produced by scripts/extract_pig_atlas_expression.R')
    parser.add_argument('--config', default='config.yml')
    parser.add_argument('--output', default='data/pig_atlas_cell_types.parquet')
    args = parser.parse_args()
    data = map_expression(input_file=args.input, config_file=args.config)
    utils.save_to_parquet(data, args.output)
    print(f'Wrote {len(data):,} pig cell-type expression rows to {args.output}')
