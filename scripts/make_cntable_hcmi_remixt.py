import os
import warnings

import click
import pandas as pd

import wgs_analysis
import gdan_hcmi_tools.copynumber

def parse_remixt_data(input, bin_size, seg_cols, sample_id):
    """
    Parses and processes consensus copy number data from multiple files in a directory, rebins it, 
    and structures it for downstream analysis.

    Parameters
    ----------
    cn_dir : str or Path
        Directory path containing consensus copy number segmentation files with a ".segtab.txt" extension.
    bin_size : int
        Bin size for rebinning the copy number data.
    seg_cols : list of str
        List of column names for segmentation data required by the `rebin` function. 
        This typically includes columns like 'Chromosome', 'Start.bp', 'End.bp'.

    Returns
    -------
    dict
        A dictionary where keys are `sample_id`s and values are DataFrames containing rebinned copy 
        number data for each sample. Each DataFrame is indexed by ['chr', 'start', 'end'].
    """
    case_id = sample_id.rsplit('-', 1)[0]
    sample_data = pd.read_table(input).assign(case=case_id, sample=sample_id)
    sample_data = gdan_hcmi_tools.copynumber.rebin(sample_data, bin_size, seg_cols, cn_cols=['chromosome', 'start', 'end'])
    sample_data = sample_data.rename(columns={'chromosome': 'chr'})
    sample_data = sample_data.set_index(['chr', 'start', 'end'])
    return sample_data


@click.command()
@click.option('-i', '--input', help='input paths file')
@click.option('-b', '--bin_size', type=str, default='50k', help="bin size str")
@click.option('-o', '--output', help="output sample cn data table")
@click.option('-s', '--sample_id', help="sample ID")
def make_cntable_hcmi_consensus(input, bin_size, output, sample_id):
    if type(bin_size) == str:
        bin_size = int(bin_size.replace('k', '0'*3).replace('M', '0'*6))
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
    cntable = parse_remixt_data(input, bin_size, seg_cols, sample_id)
    cntable.to_csv(output, index=True)

if __name__ == '__main__':
    make_cntable_hcmi_consensus()
