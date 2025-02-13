import pickle

import click
import anndata as ad
import pandas as pd

from gdan_hcmi_tools.process import get_samples_to_exclude

@click.command()
@click.option('--cohort', type=str, help="input cohort pickle path")
@click.option('--cn_anndata', type=str, help="input gene-level cn anndata")
@click.option('--mutation_anndata', type=str, help="output annotated anndata for mutations")
@click.option('--qc', type=str, help="input qc table for sample filtering")
@click.option('-o', '--output', type=str, help="output cohort pickle path")
def make_samples_list(cohort, cn_anndata, mutation_anndata, qc, output):
    with open(cohort, 'rb') as handle:
        cohort = pickle.load(handle)
    adata = ad.read_h5ad(mutation_anndata)
    gdata = ad.read_h5ad(cn_anndata)
    shared_ids = set(gdata.obs.index) & set(adata.obs.index)
    obs = gdata.obs.copy().reset_index()
    obs = obs[obs['sample_id'].isin(shared_ids)]
    obs['sample_type'] = obs['sample_type'].replace({'model_expanded':'model'})
    sample_type_counts = obs.groupby('case')['sample_type'].value_counts().unstack(fill_value=0)
    sample_type_counts = sample_type_counts[(sample_type_counts['model']>0) & (sample_type_counts['tumor']>0)]
    shared_cases = sample_type_counts.index.tolist()

    gdata = gdata[(gdata.obs['case'].isin(shared_cases)) & (gdata.obs.index.isin(shared_ids))]
    model_ids = gdata.obs[gdata.obs['sample_type'].str.startswith('model')].index.tolist()
    tumor_ids = gdata.obs[gdata.obs['sample_type']=='tumor'].index.tolist()
    samples_to_exclude = get_samples_to_exclude(qc)
    model_ids = [s for s in model_ids if s not in samples_to_exclude]
    tumor_ids = [s for s in tumor_ids if s not in samples_to_exclude]

    cohort.model_ids = model_ids # tumor_ids get automatic update as well d/t property decorator
    cohort.tumor_ids = tumor_ids # model ids get automatic update as well d/t property decorator
    print(f'cohort.sample_ids set: {len(cohort.sample_ids)}')
    with open(output, 'wb') as handle:
        pickle.dump(cohort, handle)

if __name__ == "__main__":
    make_samples_list()
