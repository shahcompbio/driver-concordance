import os
import warnings

import click
import tqdm
import pandas as pd

import wgs_analysis
import gdan_hcmi_tools.copynumber
from gdan_hcmi_tools.process import make_anndata, make_obs_var

@click.command()
@click.argument('inputs', nargs=-1)
@click.option('-o', '--output', help="output h5ad path")
def make_remixt_anndata(inputs, output):
    wgs_analysis.refgenome.set_genome_version('hg38')
    warnings.simplefilter(action='ignore')
    seg_cols = [
        'major_readcount', 'minor_readcount', 'readcount', 'allele_ratio',
        'major_depth', 'minor_depth', 'total_depth', 'major_0', 'minor_0',
        'major_1', 'minor_1', 'major_2', 'minor_2', 'major_raw', 'minor_raw',
        'major_depth_e', 'minor_depth_e', 'total_depth_e', 'major_e', 'minor_e',
        'total_e', 'major_raw_e', 'minor_raw_e', 'major_diff', 'minor_diff',
        'prob_is_outlier_total', 'prob_is_outlier_allele',
        'total_likelihood_mask', 'allele_likelihood_mask',
    ]
    data = {}
    for path in tqdm.tqdm(inputs):
        _, fname = os.path.split(path)
        assert '.cn.csv' in fname, fname
        sample_id = fname.replace('.cn.csv', '')
        sample_data = pd.read_csv(path)
        sample_data.set_index(['chr', 'start', 'end'], inplace=True)
        data[sample_id] = sample_data
    data, obs, var = make_obs_var(data)
    adata = make_anndata(data, obs, var, seg_cols)
    adata.write_h5ad(output)

if __name__ == '__main__':
    make_remixt_anndata()
