import os
import glob
import click
import pickle

import tqdm
import pandas as pd


def make_aggregated_table(variant_data_dir):
    """
    Aggregates variant data from multiple CSV files in a directory into a single DataFrame.

    Parameters
    ----------
    variant_data_dir : str or Path
        Path to the directory containing CSV files with variant data. Each file is expected 
        to represent a single sample, with the filename (minus '.csv') used as the sample ID.

    Returns
    -------
    DataFrame
        A concatenated DataFrame containing variant data from all files in the specified directory.
        The resulting DataFrame aggregates all sample data for downstream analysis.

    """
    var_data_paths = glob.glob(f'{variant_data_dir}/*.csv')
    sample_ids = [os.path.split(p)[-1].replace('.csv', '') for p in var_data_paths]
    agg = pd.DataFrame()
    for sample_id in tqdm.tqdm(sample_ids):
        var_data_path = f'{variant_data_dir}/{sample_id}.csv'
        variant_data = pd.read_csv(var_data_path)
        agg = pd.concat([agg, variant_data])
    return agg

def add_sample_data(agg, cohort):
    agg['sample_type'] = agg['sample'].map(cohort.sample2sampletype)
    agg['cancer_type'] = agg['case'].map(cohort.case2cancertype)
    return agg

@click.command()
@click.option('--cohort_data', type=str, help="input pickle path for the processed cohort data")
@click.option('--variant_data_dir', type=str, help="input directory with simple somatic mutation event tables")
@click.option('--aggregated_output', type=str, help="output aggregated mutation event table")
def aggregate_variant_tables(cohort_data, variant_data_dir, aggregated_output):
    with open(cohort_data, 'rb') as handle:
        cohort = pickle.load(handle)
    
    agg = make_aggregated_table(variant_data_dir)
    agg = add_sample_data(agg, cohort)
    agg.to_csv(aggregated_output, index=False)

if __name__ == "__main__":
    aggregate_variant_tables()