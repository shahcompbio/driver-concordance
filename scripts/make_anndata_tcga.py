import os
import glob
import warnings

import click
import tqdm
import pandas as pd
import anndata as ad

import wgs_analysis
import gdan_hcmi_tools.copynumber

def parse_tcga_data(cn_dir, bin_size, seg_cols):
    tcga_cn_paths = glob.glob(f'{cn_dir}/T*.segments.tsv')

    data = {}

    for path in tqdm.tqdm(tcga_cn_paths):
        _, file_name = os.path.split(path)
        sample_id = file_name.split('.')[0]
        case_id = sample_id.rsplit('-', 1)[0]
        sample_data = pd.read_csv(path, sep='\t').assign(case=case_id, sample=sample_id)
        sample_data.columns = sample_data.columns.str.lower()

        sample_data = gdan_hcmi_tools.copynumber.rebin(sample_data, bin_size, seg_cols, cn_cols=['chromosome', 'start', 'end'])
        sample_data = sample_data.rename(columns={'chromosome': 'chr'})
        sample_data = sample_data.set_index(['chr', 'start', 'end'])
        
        data[sample_id] = sample_data
    return data

def make_obs_var(data):
    sample_ids = list(data.keys())
    data = pd.concat(data, axis=1)

    obs = pd.DataFrame(index=sample_ids)

    var = data.reset_index()[['chr', 'start', 'end']]
    var['bin'] = var['chr'].astype(str) + ':' + var['start'].astype(str) + '-' + var['end'].astype(str)
    var.columns = var.columns.get_level_values(0)
    var = var.set_index('bin')

    return data, obs, var

def make_anndata(data, obs, var, seg_cols):
    layers = {}
    for layer in tqdm.tqdm(seg_cols):
        layer_data = data.loc[:, (slice(None), layer)]
        layer_data.columns = layer_data.columns.droplevel(1)
        layer_data = layer_data.reset_index()
        layer_data['bin'] = layer_data['chr'].astype(str) + ':' + layer_data['start'].astype(str) + '-' + layer_data['end'].astype(str)
        layer_data = layer_data.set_index('bin')
        layer_data = layer_data.loc[var.index, obs.index].T
        layers[layer] = layer_data

    adata = ad.AnnData(var=var, obs=obs, layers=layers)
    return adata

@click.command()
@click.option('-c', '--cn', help='input copy number data', required=True)
@click.option('-b', '--bin_size', type=str, default='50k', help="bin size str")
@click.option('-o', '--output', help="output h5ad path", required=True)
def make_tcga_anndata(cn, bin_size, output):
    if type(bin_size) == str:
        bin_size = int(bin_size.replace('k', '0'*3).replace('M', '0'*6))
    wgs_analysis.refgenome.set_genome_version('hg38')
    warnings.simplefilter(action='ignore')
    seg_cols = ['copy_number', 'minor_copy_number', 'major_copy_number']
    data = parse_tcga_data(cn, bin_size, seg_cols)
    data, obs, var = make_obs_var(data)
    adata = make_anndata(data, obs, var, seg_cols)
    adata.write_h5ad(output)

if __name__ == '__main__':
    make_tcga_anndata()
