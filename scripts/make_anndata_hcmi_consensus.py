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
def make_consensus_anndata(inputs, output):
    wgs_analysis.refgenome.set_genome_version('hg38')
    warnings.simplefilter(action='ignore')
    seg_cols = [
        "total_copy_ratio", "modal_total_cn", "expected_total_cn", "total_HZ", 
        "total_amp", "corrected_total_cn", "rescaled_total_cn", "bi.allelic", 
        "copy.ratio", "hscr.a1", "hscr.a2", "modal.a1", "modal.a2", "expected.a1", "expected.a2", 
        "subclonal.a1", "subclonal.a2", "cancer.cell.frac.a1", "ccf.ci95.low.a1", "ccf.ci95.high.a1", 
        "cancer.cell.frac.a2", "ccf.ci95.low.a2", "ccf.ci95.high.a2", "LOH", "HZ", "SC_HZ", 
        "amp.a1", "amp.a2", "rescaled.cn.a1", "rescaled.cn.a2", 
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
    make_consensus_anndata()
