import os
import click
import glob
import pickle

import tqdm
import numpy as np
import pandas as pd
import anndata as ad

from gdan_hcmi_tools.process import VEP_meta


def init_layers_data(samples, driver_genes, vep_consequences=('T', 'I', 'D', 'M')):
    """
    Initialize a dictionary of DataFrames representing layers of mutation events for samples and driver genes.

    This function creates a dictionary of DataFrames, where each DataFrame corresponds to a specific type of variant effect predictor (VEP) consequence.
    Each DataFrame is initialized with zero values and has samples as its index and driver genes as columns.

    Parameters
    ----------
    samples : list
        A list of unique sample identifiers.
    driver_genes : list
        A list of driver gene names that serve as column headers for each DataFrame.
    vep_consequences : tuple of str, optional
        A tuple of VEP consequence types (e.g., 'T', 'I', 'D', 'M'), each representing a type of mutation effect.
        These are used as keys for the output dictionary. Default is ('T', 'I', 'D', 'M').

    Returns
    -------
    dict
        A dictionary where:
        - Each key is a VEP consequence (e.g., 'T', 'I', 'D', 'M').
        - Each value is a DataFrame with rows indexed by `samples` and columns named after `driver_genes`.
        - All DataFrames are initialized to zero values.
    """
    assert len(samples) == len(set(samples)), samples
    n_samples = len(samples)
    n_drivers = len(driver_genes)
    layers_data = {
        v: pd.DataFrame(index=samples, columns=driver_genes,
                        data=np.zeros((n_samples, n_drivers)))
        for v in vep_consequences
    }
    return layers_data

def make_layers_data(variant_data_dir, sample_ids, driver_genes, vep_config):
    """
    Create and populate mutation event layers data from variant files for a set of samples and driver genes.

    This function initializes data layers for specific VEP consequence types and fills them based on variant event data for each sample.
    Each layer represents a different consequence type, and within each layer, mutations are recorded per sample and gene.

    Parameters
    ----------
    variant_data_dir : str
        The directory path containing variant data files. Each file should be named with the sample ID and have a '.csv' extension.
    sample_ids : list
        A list of unique sample identifiers. Each sample ID corresponds to a file in `variant_data_dir`.
    driver_genes : list
        A list of driver gene names to be used as column headers in the layer DataFrames.
    vep_config : object
        A configuration object containing:
        - `layer_map` : dict
            A mapping of event types (e.g., 'T', 'I', 'D', 'M') to specific VEP consequence labels.
            These labels determine the types of mutations recorded in each layer.

    Returns
    -------
    dict
        A dictionary of DataFrames representing layers for different event types:
        - Each key is an event type, and each value is a DataFrame.
        - Each DataFrame has rows indexed by `sample_ids` and columns for `driver_genes`.
        - Cells contain binary values (1 for an event present in a sample for a gene, 0 if absent).
    """
    layers_data = init_layers_data(sample_ids, driver_genes, vep_consequences=vep_config.layer_map.values())
    for sample_id in tqdm.tqdm(sample_ids):
        var_data_path = f'{variant_data_dir}/{sample_id}.csv'
        # case_id = sample_id[:17]
        variant_data = pd.read_csv(var_data_path)
        for _, row in variant_data.iterrows():
            sample_id = row['sample']
            gene_symbol = row['gene']
            event = row['event']
            layers_data[event].loc[sample_id, gene_symbol] = 1 # e.g. if >1 missense in same gene, collapsed to 1
    return layers_data

def get_mutation_anndata(variant_data_dir, driver_genes):
    vep_config = VEP_meta()
    var_data_paths = glob.glob(f'{variant_data_dir}/*.csv')
    sample_ids = [os.path.split(p)[-1].replace('.csv', '') for p in var_data_paths]
    layers_data = make_layers_data(variant_data_dir, sample_ids, driver_genes, vep_config)
    obs = pd.DataFrame(index=sample_ids)
    var = pd.DataFrame(index=driver_genes)
    adata = ad.AnnData(obs=obs, var=var, layers=layers_data)
    return adata

@click.command()
@click.option('--cohort_data', type=str, help="input pickle path for the processed cohort data")
@click.option('--variant_data_dir', type=str, help="input directory with simple somatic mutation event tables")
@click.option('--drivers_mutation_path', type=str, help="input driver genes table for simple somatic mutations")
@click.option('--mutation_anndata', type=str, help="output annotated anndata for mutations")
def make_snvindel_anndata(cohort_data, variant_data_dir, drivers_mutation_path, mutation_anndata):
    with open(cohort_data, 'rb') as handle:
        cohort = pickle.load(handle)
    drivers_mutation_df = pd.read_table(drivers_mutation_path, comment='#')
    drivers_mutation = drivers_mutation_df['Gene'].unique().tolist()

    if os.path.exists(mutation_anndata):
        adata = ad.read_h5ad(mutation_anndata)
    else:
        adata = get_mutation_anndata(variant_data_dir, drivers_mutation)
        adata.obs.index.name = 'sample_id'
        adata.obs['case'] = adata.obs.index.str.slice(0, 17)
        adata.obs['sample_type'] = adata.obs.index.map(cohort.sample2sampletype)
        adata.obs['cancer_type'] = adata.obs['case'].map(cohort.case2cancertype)
        adata.obs['consensus_purity'] = adata.obs.index.map(cohort.sample2purity)
        adata.obs['consensus_ploidy'] = adata.obs.index.map(cohort.sample2ploidy)
        adata.write_h5ad(mutation_anndata)
    print('adata read')

if __name__ == "__main__":
    make_snvindel_anndata()