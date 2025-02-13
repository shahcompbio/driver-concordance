import pickle

import click
import tqdm
import numpy as np
import pandas as pd
import anndata as ad

import scgenome

import gdan_hcmi_tools.process


def compare_loh(data):
    data = data[['is_loh_tumor', 'is_loh_model']].dropna()
    loh_either = data[(data['is_loh_tumor'] == 1) | (data['is_loh_model'] == 1)]
    return {
        'loh_union': ((data['is_loh_tumor'] == 1) | (data['is_loh_model'] == 1)).mean(),
        'loh_match': (loh_either['is_loh_tumor'] == loh_either['is_loh_model']).mean(),
        'loh_shared': (data['is_loh_tumor'] == data['is_loh_model']).mean(),
        'loh_model': ((data['is_loh_tumor'] == 0) & (data['is_loh_model'] == 1)).mean(),
        'loh_tumor': ((data['is_loh_tumor'] == 1) & (data['is_loh_model'] == 0)).mean(),
    }


def compute_cn_distance_wgd(data, tumor_n_wgd, model_n_wgd):
    data = data[['total_cn_tumor', 'total_cn_model']].dropna().round()
    if tumor_n_wgd > model_n_wgd:
        data['total_cn_model'] *= 2 ** (tumor_n_wgd - model_n_wgd)
    else:
        data['total_cn_tumor'] *= 2 ** (model_n_wgd - tumor_n_wgd)
    return {
        'cn_equal_wgd': (data['total_cn_tumor'] == data['total_cn_model']).mean(),
    }


def compute_cn_distance(data):
    data = data[['total_cn_tumor', 'total_cn_model']].dropna().round()
    return {
        'cn_equal': (data['total_cn_tumor'] == data['total_cn_model']).mean(),
    }


def compute_cn_mean_sq_diff(data):
    data = data[['total_cn_tumor', 'total_cn_model']].dropna().round()
    return {
        'cn_mean_sq_diff': np.square(data['total_cn_tumor'] - data['total_cn_model']).mean(),
    }


def compute_cn_mean_sq_diff_wgd(data, tumor_n_wgd, model_n_wgd):
    data = data[['total_cn_tumor', 'total_cn_model']].dropna().round()
    if tumor_n_wgd > model_n_wgd:
        data['total_cn_model'] *= 2 ** (tumor_n_wgd - model_n_wgd)
    else:
        data['total_cn_tumor'] *= 2 ** (model_n_wgd - tumor_n_wgd)
    return {
        'cn_mean_sq_diff_wgd': np.square(data['total_cn_tumor'] - data['total_cn_model']).mean(),
    }


def compute_cn_correlation(data):
    correlation = data[['total_cn_tumor', 'total_cn_model']].dropna().round().corr().loc['total_cn_tumor', 'total_cn_model']
    return {
        'cn_correlation': correlation
    }


def get_pairwise_info(models, adata):

    adata.layers['is_loh'] = 1 * (adata.layers['minor_cn'] < 0.5)

    pairwise_info = []

    for idx, row in tqdm.tqdm(models.iterrows()):
        model_id = row['model']
        model_type = row['model_type']
        tumor_id = row['tumor']
        case = row['case']

        if model_type == 'model_expanded': # omit expanded models as pairs
            continue

        if tumor_id not in adata.obs.index:
            continue

        if model_id not in adata.obs.index:
            continue

        tumor_data = scgenome.tl.get_obs_data(adata, tumor_id, layer_names=['total_cn', 'major_cn', 'minor_cn', 'is_loh'])
        tumor_n_wgd = adata.obs.loc[tumor_id, 'n_wgd']
        tumor_fraction_diploid = ((tumor_data['minor_cn'] == 1) & (tumor_data['major_cn'] == 1)).mean()

        cancer_type = adata.obs.loc[tumor_id, 'cancer_type']
        
        model_data = scgenome.tl.get_obs_data(adata, model_id, layer_names=['total_cn', 'major_cn', 'minor_cn', 'is_loh'])
        model_n_wgd = adata.obs.loc[model_id, 'n_wgd']
        model_fraction_diploid = ((model_data['minor_cn'] == 1) & (model_data['major_cn'] == 1)).mean()

        combined_data = tumor_data.merge(model_data, left_index=True, right_index=True, suffixes=['_tumor', '_model'])

        stats = {
            'case': case,
            'tumor_id': tumor_id,
            'tumor_n_wgd': tumor_n_wgd,
            'tumor_fraction_diploid': tumor_fraction_diploid,
            'model_id': model_id,
            'model_type': model_type,
            'model_n_wgd': model_n_wgd,
            'model_fraction_diploid': model_fraction_diploid,
            'cancer_type': cancer_type,
        }

        stats.update(compare_loh(combined_data))
        stats.update(compute_cn_distance(combined_data))
        stats.update(compute_cn_distance_wgd(combined_data, tumor_n_wgd, model_n_wgd))
        stats.update(compute_cn_mean_sq_diff(combined_data))
        stats.update(compute_cn_mean_sq_diff_wgd(combined_data, tumor_n_wgd, model_n_wgd))
        stats.update(compute_cn_correlation(combined_data))

        pairwise_info.append(stats)

    pairwise_info = pd.DataFrame(pairwise_info)

    # Add tumor and model purity and ploidy
    pairwise_info['tumor_consensus_purity'] = adata.obs.loc[pairwise_info['tumor_id'].values, 'consensus_purity'].values
    pairwise_info['tumor_consensus_ploidy'] = adata.obs.loc[pairwise_info['tumor_id'].values, 'consensus_ploidy'].values
    pairwise_info['model_consensus_purity'] = adata.obs.loc[pairwise_info['model_id'].values, 'consensus_purity'].values
    pairwise_info['model_consensus_ploidy'] = adata.obs.loc[pairwise_info['model_id'].values, 'consensus_ploidy'].values

    # Add WGD status
    pairwise_info['is_wgd_tumor'] = pairwise_info['tumor_n_wgd'] > 0
    pairwise_info['is_wgd_model'] = pairwise_info['model_n_wgd'] > 0

    pairwise_info['wgd_status'] = 'Both diploid'
    pairwise_info.loc[(pairwise_info['is_wgd_model'] & ~pairwise_info['is_wgd_tumor']), 'wgd_status'] = 'Model WGD'
    pairwise_info.loc[(pairwise_info['is_wgd_tumor'] & ~pairwise_info['is_wgd_model']), 'wgd_status'] = 'Tumor WGD'
    pairwise_info.loc[(pairwise_info['is_wgd_tumor'] & pairwise_info['is_wgd_model']), 'wgd_status'] = 'Both WGD'

    return pairwise_info


@click.command()
@click.option('--input_tm_link', required=True)
@click.option('--input_cn_adata', required=True)
@click.option('--output_pairwise_csv', required=True)
def make_fig2C_data(input_tm_link, input_cn_adata, output_pairwise_csv):
    tm_link = gdan_hcmi_tools.process.read_tm_link(input_tm_link)
    models = gdan_hcmi_tools.process.generate_model_table(tm_link)

    adata = ad.read_h5ad(input_cn_adata)
    adata.var['chr'] = 'chr' + adata.var['chr'].astype(str)
    
    pairwise_info = get_pairwise_info(models, adata)
    pairwise_info.to_csv(output_pairwise_csv, index=False)


if __name__ == "__main__":
    make_fig2C_data()
