import os
import glob

import pandas as pd
import pyranges as pr
import anndata as ad
import numpy as np
import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import wgs_analysis.refgenome
from wgs_analysis.plots.cnv import plot_cnv_segments
import matplotlib


def add_fraction_genome_altered(adata):
    """
    Add the fraction of the genome altered (FGA) to the given AnnData object.

    Args:
        adata (AnnData): An AnnData object containing copy number data in `layers`
                         and genomic coordinates in `var`. The 'major_cn' and 'minor_cn'
                         layers should represent major and minor copy numbers respectively.
    
    Adds:
        fga (pandas.Series): A column 'fga' in `adata.obs` that contains the fraction 
                             of the genome altered for each observation.
    """
    haplo_ploidy = 1 # major and minor CN should be this
    assert 'major_cn' in adata.layers, adata.layers
    assert 'minor_cn' in adata.layers, adata.layers
    major_cn = adata.layers['major_cn'].round(0)
    minor_cn = adata.layers['minor_cn'].round(0)
    altered = ((major_cn != haplo_ploidy) | (minor_cn != haplo_ploidy)).astype(int)
    lengths = adata.var['end'] - adata.var['start']
    assert (lengths < 0).sum() == 0, lengths.value_counts()
    fga = np.dot(altered, lengths) / lengths.sum()
    adata.obs['fga'] = fga

def add_loh_fraction(adata):
    """
    Calculate and add the fraction of the genome with loss of heterozygosity (LOH) to the given AnnData object.

    Args:
        adata (AnnData): An AnnData object containing LOH data in the `layers` and genomic coordinates in `var`. 
                         The 'LOH' layer should represent the presence of LOH across genomic regions.

    Adds:
        loh (pandas.Series): A column 'loh' in `adata.obs` that contains the calculated fraction of the genome 
                             affected by LOH for each observation.
    """
    assert 'LOH' in adata.layers, adata.layers
    loh = adata.layers['LOH']
    lengths = adata.var['end'] - adata.var['start']
    assert (lengths < 0).sum() == 0, lengths.value_counts()
    segment_loh = (loh.round(0) == 1).astype(int)
    loh_fraction = np.dot(segment_loh, lengths) / lengths.sum()
    assert loh_fraction.max() <= 1 and loh_fraction.min() >= 0
    adata.obs['loh'] = loh_fraction

def add_calculated_ploidy(adata, max_seg_cn_for_ploidy=6):
    """
    Calculate and add the average ploidy to the given AnnData object.

    Args:
        adata (AnnData): An AnnData object containing copy number data in the `layers` and
                         genomic coordinates in `var`. The 'total_cn' layer should represent
                         the total copy number.

    Adds:
        ploidy (pandas.Series): A column 'ploidy' in `adata.obs` that contains the calculated
                                average ploidy for each observation.
    """
    default_cn = 2
    cn = np.nan_to_num(adata.layers['total_cn'], nan=default_cn, posinf=max_seg_cn_for_ploidy, neginf=0).round(0)
    cn = cn.clip(min=0, max=max_seg_cn_for_ploidy)
    lengths = np.abs(adata.var['end'] - adata.var['start'])
    assert (lengths < 0).sum() == 0, lengths.value_counts()
    assert cn.shape[1] == lengths.shape[0], f"shape mismatch: {cn.shape} == {lengths.shape}"
    ploidy = np.dot(cn, lengths) / lengths.sum()
    adata.obs['ploidy'] = ploidy
    return adata

def plot_cnv_genome(ax, cnv, maxcopies=4, minlength=1000, major_col='major_raw', minor_col='minor_raw', scatter=False, squashy=False, rasterized=False):
    """
    Plot major/minor copy number across the genome

    Args:
        ax (matplotlib.axes.Axes): plot axes
        cnv (pandas.DataFrame): `cnv_site` table
        maxcopies (int): maximum number of copies for setting y limits
        minlength (int): minimum length of segments to be drawn
        major_col (str): name of column to use as major copy number
        minor_col (str): name of column to use as minor copy number
        scatter (boolean): display segments as scatter points not segments
        squashy (boolean): squash the y axis to display all copy numbers

    """

    segment_color_major = plt.get_cmap('RdBu')(0.1)
    segment_color_minor = plt.get_cmap('RdBu')(0.9)

    cnv = cnv.copy()

    squash_coeff = 0.15
    squash_f = lambda a: np.tanh(squash_coeff * a)
    if squashy:
        cnv[major_col] = squash_f(cnv[major_col])
        cnv[minor_col] = squash_f(cnv[minor_col])

    if 'length' not in cnv:
        cnv['length'] = cnv['end'] - cnv['start']

    # cnv = cnv[['chromosome', 'start', 'end', 'length', major_col, minor_col]]

    cnv = cnv[cnv['length'] >= minlength]
    
    cnv = cnv[cnv['chromosome'].isin(wgs_analysis.refgenome.info.chromosomes)]

    cnv.set_index('chromosome', inplace=True)
    cnv['chromosome_start'] = wgs_analysis.refgenome.info.chromosome_start
    cnv.reset_index(inplace=True)

    cnv['start'] = cnv['start'] + cnv['chromosome_start']
    cnv['end'] = cnv['end'] + cnv['chromosome_start']

    if scatter:
        cnv['mid'] = 0.5 * (cnv['start'] + cnv['end'])
        for column, color in ((minor_col, segment_color_minor), (major_col, segment_color_major)):
            clipped_cnv = cnv[cnv[column] < maxcopies]
            amp_cnv = cnv[cnv[column] >= maxcopies]
            sns.scatterplot(ax=ax, x='mid', y=column, data=clipped_cnv, color=color, s=1, linewidth=0, alpha=1, rasterized=rasterized)
            # ax.scatter(clipped_cnv['mid'], clipped_cnv[column], color=color, s=10, alpha=0.1)
            # ax.scatter(amp_cnv['mid'], np.ones(amp_cnv.shape[0]) * maxcopies, color=color, s=30)

    else:
        plot_cnv_segments(ax, cnv, column=major_col, segment_color=segment_color_major)
        plot_cnv_segments(ax, cnv, column=minor_col, segment_color=segment_color_minor)

    #  ax.fill_between(x=cnv['start'], y1=8*cnv['is_subclonal'], where=cnv['is_subclonal']>0,
                    #  color='tab:orange', interpolate=False, step='post', alpha=0.3, rasterized=rasterized)

    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlim((-0.5, wgs_analysis.refgenome.info.chromosome_end.max()))
    ax.set_xlabel('chromosome')
    ax.set_xticks([0] + list(wgs_analysis.refgenome.info.chromosome_end.values))
    ax.set_xticklabels([])

    if squashy:
        yticks = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        ax.set_yticks(yticks_squashed)
        ax.set_yticklabels(ytick_labels)
        ax.set_ylim((-0.01, 1.01))
        ax.spines['left'].set_bounds(0, 1)
    else:
        ax.set_ylim((-0.05*maxcopies, maxcopies))
        ax.set_yticks(range(0, int(maxcopies) + 1))
        ax.spines['left'].set_bounds(0, maxcopies)

    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(wgs_analysis.refgenome.info.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter([a.replace('chr', '') for a in wgs_analysis.refgenome.info.chromosomes]))

    ax.xaxis.grid(True, which='major', linestyle=':')
    ax.yaxis.grid(True, which='major', linestyle=':')


def plot_cn_genome(ax, df, major_col='major_raw_e', minor_col='minor_raw_e', maxcopies=5, genome_version='grch38', scatter=False, squashy=False):
    """
    From given major and minor CN cols and a dataframe, plot genome-wide CN plots in an axis
    """
    refinfo = wgs_analysis.refgenome.RefGenomeInfo('grch38')
    chroms = [c.replace('chr', '') for c in refinfo.chromosomes]
    df = df.rename(columns={'chr':'chromosome'}).copy()
    df['chromosome'] = 'chr'+df['chromosome'].astype(str).str.replace('chr', '')
    wgs_analysis.refgenome.set_genome_version(genome_version)
    refinfo = wgs_analysis.refgenome.RefGenomeInfo(genome_version)
    plot_cnv_genome(ax, df, minor_col=minor_col, major_col=major_col, maxcopies=maxcopies, scatter=scatter, squashy=squashy, rasterized=True)
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refinfo.chromosome_mid))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(chroms))


def rebin(data:pd.DataFrame, bin_size:int, cols:list,  # TODO: add test
          ref_cols:list=['Chromosome', 'Start', 'End'], 
          cn_cols:list=['Chromosome', 'Start.bp', 'End.bp'],
          genome_version:str='grch38') -> pd.DataFrame:
    """ 
    Rebin the given data to the given bin size

    Args:
        data (pd.DataFrame): dataframe with columns to rebin
        bin_size (int): size of bins
        cols (list): columns to rebin
        ref_cols (list): columns for parsing reference table
        cn_cols (list): columns for parsing cn data input
        genome_version (str): genome version to use

    Returns:
        pd.DataFrame: rebinned data
    """
    wgs_analysis.refgenome.set_genome_version(genome_version)
    refinfo = wgs_analysis.refgenome.RefGenomeInfo(genome_version)
    chromsizes = pd.DataFrame(refinfo.chromosome_lengths).reset_index()
    chromsizes.columns = [ref_cols[0], ref_cols[2]]
    chromsizes[ref_cols[0]] = chromsizes[ref_cols[0]].astype(str).str.replace('chr', '')
    chromsizes[ref_cols[1]] = 0
    chromsizes = chromsizes[ref_cols]
    chromsizes = pr.PyRanges(chromsizes)

    bins = pr.gf.tile_genome(chromsizes, bin_size)
    bins = bins.insert(pd.Series(range(len(bins)), name='bin'))

    data[cn_cols[0]] = data[cn_cols[0]].astype(str).str.replace('chr', '')
    data2 = pr.PyRanges(data.rename(columns=dict(zip(cn_cols, ref_cols))))

    intersect_1 = data2.intersect(bins)
    intersect_2 = bins.intersect(data2)

    intersect = pd.merge(
        intersect_1.as_df(), intersect_2.as_df(),
        on=ref_cols, how='outer')

    intersect['length'] = intersect['End'] - intersect['Start'] + 1

    for col in cols:
        intersect[col] = intersect[col] * intersect['length'] / bin_size
    intersect = intersect.groupby('bin')[cols].sum().reset_index()
    intersect = intersect.merge(bins.as_df())
    intersect = intersect.rename(columns={
        ref_cols[0]: 'chromosome',
        ref_cols[1]: 'start',
        ref_cols[2]: 'end',
    })

    return intersect[['chromosome', 'start', 'end'] + cols]


def get_remixt_anndata(path):
    remixt = ad.read_h5ad(path)
    assert remixt.layers['major_raw'].clip(max=8).sum() > 78830607
    assert remixt.layers['minor_raw'].clip(max=8).sum() > 38556104
    assert remixt.layers['major_raw'].shape == (659, 58962)
    
    remixt.layers['LOH'] = remixt.layers['minor_raw_e'] < 0.5
    remixt.obs['LOH_fraction'] = np.array(remixt.layers['LOH'].mean(axis=1))

    return remixt


def fill_ploidies_if_empty(paired_data, ploidies, col, modes=['tumor', 'model']):
    """
    Calculate and fill in empty ploidy value if 0
    paired_data: dataframe with _tumor, _model suffices
    ploidies: list of [tumor_purity, model_purity]. if 0, calculate sample ploidy
    col: suffix filling in column names for major_{} and minor_{}
    modes: suffices for paired data
    """
    inf_map = {-np.inf: 0, np.inf: 8}
    for pix, ploidy in enumerate(ploidies):
        mode = modes[pix]
        major = paired_data[f'major_{col}_{mode}'].replace(inf_map)
        minor = paired_data[f'minor_{col}_{mode}'].replace(inf_map)
        total = major + minor
        if ploidy == 0:
            ploidy = (total * paired_data['length']).sum() / paired_data['length'].sum()
            ploidies[pix] = ploidy
    return ploidies

def fill_ploidies_if_empty_consensus(paired_data, ploidies, major_col='rescaled.cn.a2', minor_col='rescaled.cn.a1', modes=['tumor', 'model']):
    """
    Calculate and fill in empty ploidy value if 0
    paired_data: dataframe with _tumor, _model suffices
    ploidies: list of [tumor_purity, model_purity]. if 0, calculate sample ploidy
    col: suffix filling in column names for major_{} and minor_{}
    modes: suffices for paired data
    """
    inf_map = {-np.inf: 0, np.inf: 8}
    for pix, ploidy in enumerate(ploidies):
        mode = modes[pix]
        major = paired_data[f'{major_col}_{mode}'].replace(inf_map)
        minor = paired_data[f'{minor_col}_{mode}'].replace(inf_map)
        total = major + minor
        if ploidy == 0:
            ploidy = (total * paired_data['length']).sum() / paired_data['length'].sum()
            ploidies[pix] = ploidy
    return ploidies


def compute_cn_distance(paired_data, ploidies, col="raw_e", modes=['tumor', 'model']):
    """
    Compute distance  between the CN of two CN samples, considering WGD status
    paired_data: dataframe with _tumor, _model suffices
    ploidies: list of [tumor_purity, model_purity]. if 0, calculate sample ploidy
    col: suffix filling in column names for major_{} and minor_{}
    modes: suffices to use for paired data
    """
    paired_data = paired_data.copy()
    dist = pd.DataFrame()
    paired_data['length'] = paired_data['end'] - paired_data['start']
    dist['length'] = paired_data['length']
    # print(f'initial ploidies: {ploidies}')
    ploidy1, ploidy2 = fill_ploidies_if_empty(paired_data, ploidies, "raw_e", modes=modes)
    # print(f'fixed ploidies: {tumor_ploidy}, {model_ploidy}')

    dip, wgd = modes[::-1]#'model', 'tumor'
    if ploidy1 < ploidy2: # if model ploidy is bigger
        dip, wgd = modes #'tumor', 'model' # presume wgd sample is model
    best_dist = np.inf
    best_wix = -1
    for wix in [1, 2]: 
        for mix in ['major', 'minor']:
            dist[f'{mix}_dist_{wix}'] = np.abs((wix * paired_data[f'{mix}_{col}_{dip}']) - paired_data[f'{mix}_{col}_{wgd}'])
        dist_major = (dist[f'major_dist_{wix}'] * dist['length']).sum() / dist['length'].sum()
        dist_minor = (dist[f'minor_dist_{wix}'] * dist['length']).sum() / dist['length'].sum()
        dist_total = (dist_major + dist_minor) / 2
        if dist_total < best_dist:
            best_dist = dist_total
            best_wix = wix
    for mix in ['major', 'minor']:
        paired_data[f'{mix}_distance'] = dist[f'{mix}_dist_{best_wix}']
    logs = {'wgd':wgd, 'best_dist':best_dist, 'best_wix':best_wix}
    return paired_data, logs

def compute_cn_distance_consensus(paired_data, ploidies, major_col="rescaled.cn.a2", minor_col="rescaled.cn.a1", modes=['tumor', 'model']):
    """
    Compute distance  between the CN of two CN samples, considering WGD status
    paired_data: dataframe with _tumor, _model suffices
    ploidies: list of [tumor_purity, model_purity]. if 0, calculate sample ploidy
    col: suffix filling in column names for major_{} and minor_{}
    modes: suffices to use for paired data
    """
    paired_data = paired_data.copy()
    dist = pd.DataFrame()
    paired_data['length'] = paired_data['end'] - paired_data['start']
    dist['length'] = paired_data['length']
    # print(f'initial ploidies: {ploidies}')
    ploidy1, ploidy2 = fill_ploidies_if_empty_consensus(paired_data, ploidies, modes=modes)
    # print(f'fixed ploidies: {tumor_ploidy}, {model_ploidy}')

    dip, wgd = modes[::-1]#'model', 'tumor'
    if ploidy1 < ploidy2: # if model ploidy is bigger
        dip, wgd = modes #'tumor', 'model' # presume wgd sample is model
    best_dist = np.inf
    best_wix = -1
    for wix in [1, 2]: 
        dist[f'major_dist_{wix}'] = np.abs((wix * paired_data[f'{major_col}_{dip}']) - paired_data[f'{major_col}_{wgd}'])
        dist[f'minor_dist_{wix}'] = np.abs((wix * paired_data[f'{minor_col}_{dip}']) - paired_data[f'{minor_col}_{wgd}'])
        dist_major = (dist[f'major_dist_{wix}'] * dist['length']).sum() / dist['length'].sum()
        dist_minor = (dist[f'minor_dist_{wix}'] * dist['length']).sum() / dist['length'].sum()
        dist_total = (dist_major + dist_minor) / 2
        if dist_total < best_dist:
            best_dist = dist_total
            best_wix = wix
    for mix in ['major', 'minor']:
        paired_data[f'{mix}_distance'] = dist[f'{mix}_dist_{best_wix}']
    logs = {'wgd':wgd, 'best_dist':best_dist, 'best_wix':best_wix}
    return paired_data, logs

def compute_mean_distance(paired_data:pd.DataFrame, metric='major_distance'):
    """
    Compute the mean of input distance metric from a pandas DataFrame
    """
    assert 'length' in paired_data.columns, paired_data.head(2)
    assert paired_data[metric].isna().sum() == 0, paired_data[metric].isna().sum()
    assert np.isinf(paired_data[metric]).sum() == 0, paired_data[metric].sum()
    distance = (paired_data[metric] * paired_data['length']).sum() / paired_data['length'].sum()
    return distance


def calculate_cn_stats(adata):
    assert 'ploidy' in adata.obs.columns, adata.obs.columns
    lengths = (adata.var['end'] - adata.var['start'] + 1).values
    minor_cn = adata.layers['minor_cn']
    major_cn = adata.layers['major_cn']
    mean_allele_diff = (major_cn - minor_cn)
    cn_gt2 = (major_cn >= 2)
    cn_gt3 = (major_cn >= 3)
    cn_gt4 = (major_cn >= 4)
    is_loh = (minor_cn == 0)
    total_cn = adata.layers['total_cn']
    adata.obs['fraction_loh'] = np.array(np.nansum(is_loh * lengths, axis=1) / np.nansum(lengths))
    adata.obs['cn_gt2'] = np.array(np.nansum(cn_gt2 * lengths, axis=1) / np.nansum(lengths))
    adata.obs['cn_gt3'] = np.array(np.nansum(cn_gt3 * lengths, axis=1) / np.nansum(lengths))
    adata.obs['cn_gt4'] = np.array(np.nansum(cn_gt4 * lengths, axis=1) / np.nansum(lengths))
    adata.obs['mean_allele_diff'] = np.array(np.nansum(mean_allele_diff * lengths, axis=1) / np.nansum(lengths))
    
    return adata


def calculate_n_wgd(adata):
    """
    Calculate the number of whole genome duplications (WGD) based on copy number data.

    Parameters:
        adata (AnnData): Annotated data matrix containing copy number data.

    Returns:
        int: The number of whole genome duplications.
    """
    adata = calculate_cn_stats(adata)    
    adata.obs['n_wgd'] = 0
    adata.obs.loc[adata.obs['cn_gt2'] > 0.5, 'n_wgd'] = 1
    adata.obs.loc[adata.obs['cn_gt3'] > 0.5, 'n_wgd'] = 2
    
    return adata 

def gene_over_min_cn(gene_cn, gene_cn_remixt, min_gene_cn):
    over_cn_min_consensus = gene_cn >= min_gene_cn 
    over_cn_min_remixt = gene_cn_remixt >= min_gene_cn
    over_cn_min = (over_cn_min_consensus or over_cn_min_remixt)
    return over_cn_min

def has_hlamp(gene_cn, hlamp_cutoff, ploidy, call_cn_amp_cn=6.0, min_ploidy_cutoff=2.0, debug=False):
    """
    Determine if a gene has high-level amplification (hlamp).

    Args:
        gene_cn (float): Copy number of the gene.
        hlamp_cutoff (float): Cutoff value for high-level amplification.
        ploidy (float): Ploidy of the sample.
        call_cn_amp_cn (float): Absolute copy number threshold to call amplifications.
        min_ploidy_cutoff (float): Minimum ploidy cutoff to consider.
        debug (bool): If True, print debug information.

    Returns:
        bool: True if the gene has high-level amplification, False otherwise.
    """
    
    gene_cn_over_relative_cutoff = gene_cn >= round(hlamp_cutoff, 0)
    gene_cn_over_absolute_cutoff = gene_cn >= round(call_cn_amp_cn, 0)
    gene_cn_over_hlamp_cutoff = (gene_cn_over_relative_cutoff or gene_cn_over_absolute_cutoff)
    gene_cn_bigger_than_ploidy = (gene_cn >= ploidy + 1 and ploidy >= min_ploidy_cutoff)
    gene_has_hlamp = gene_cn_over_hlamp_cutoff and gene_cn_bigger_than_ploidy
    if debug:
        print(f'gene_cn_over_relative_cutoff:{gene_cn_over_relative_cutoff} = {gene_cn} >= {hlamp_cutoff}')
        print(f'gene_cn_over_absolute_cutoff:{gene_cn_over_absolute_cutoff} = {gene_cn} >= {call_cn_amp_cn}')
        print(f'gene_cn_over_hlamp_cutoff:{gene_cn_over_relative_cutoff} = {gene_cn_over_relative_cutoff} or {gene_cn_over_absolute_cutoff}')
        print(f'gene_cn_bigger_than_ploidy:{gene_cn_bigger_than_ploidy} = {gene_cn} >= {ploidy+1} and {ploidy} >= {min_ploidy_cutoff}')
        print(f'gene_has_hlamp:{gene_has_hlamp} = {gene_cn_over_hlamp_cutoff} and {gene_cn_bigger_than_ploidy}')
    return gene_has_hlamp

def eval_cn_event(gene_cn, gene_cn_remixt, ploidy_consensus, ploidy_remixt, hlamp_cutoff_consensus, hlamp_cutoff_remixt, 
                  min_gene_cn=2.0, call_cn_amp_cn=5.0, min_ploidy_cutoff=2, homdel_cutoff=0.5):
    """
    Evaluate copy number events for a given gene based on consensus and remixt data.

    Args:
        gene_cn (float): Copy number of the gene from consensus data.
        gene_cn_remixt (float): Copy number of the gene from remixt data.
        ploidy_consensus (float): Ploidy value from consensus data.
        ploidy_remixt (float): Ploidy value from remixt data.
        hlamp_cutoff_consensus (float): High-level amplification cutoff for consensus data.
        hlamp_cutoff_remixt (float): High-level amplification cutoff for remixt data.
        min_gene_cn (float, optional): Minimum gene copy number to consider. Default is 2.0.
        call_cn_amp_cn (float, optional): Absolute copy number threshold to call amplifications. Default is 5.0.
        min_ploidy_cutoff (float, optional): Minimum ploidy cutoff to consider. Default is 2.
        homdel_cutoff (float, optional): Cutoff value for homozygous deletions. Default is 0.5.

    Returns:
        str: 'hlamp' if the gene has high-level amplification, 'homdel' if the gene has homozygous deletion, 
             or None if no significant event is detected.
    """
    if gene_over_min_cn(gene_cn, gene_cn_remixt, min_gene_cn):
        if has_hlamp(gene_cn, hlamp_cutoff_consensus, ploidy_consensus, call_cn_amp_cn=call_cn_amp_cn, min_ploidy_cutoff=min_ploidy_cutoff):
            return 'hlamp'
        elif not np.isnan(gene_cn_remixt):
            if has_hlamp(gene_cn_remixt, hlamp_cutoff_remixt, ploidy_remixt, call_cn_amp_cn=call_cn_amp_cn, min_ploidy_cutoff=min_ploidy_cutoff):
                return 'hlamp'
                
    if gene_cn < homdel_cutoff:
        return 'homdel'
    elif not np.isnan(gene_cn_remixt):
        if gene_cn_remixt < homdel_cutoff:
            return 'homdel'
    return None

def get_cnv_events_hcmi(gdata, rgdata, gene_sets=None, 
                        homdel_cutoff=0.5, a_consensus=2, b_consensus=0.5, a_remixt=2, b_remixt=0.5, model_hlamp_cutoff_offset=-1,
                        call_cn_amp_cn=6.0, min_ploidy_cutoff=2, min_gene_cn=2):
    """
    Identify copy number variation (CNV) events such as high-level amplifications (hlamp) and homozygous deletions (homdel) 
    in the given gene data (gdata) and remixt gene data (rgdata).

    Args:
        gdata (AnnData): Annotated data matrix containing consensus copy number data.
        rgdata (AnnData): Annotated data matrix containing remixt copy number data.
        gene_sets (list, optional): List of gene sets to consider. Default is None.
        homdel_cutoff (float, optional): Cutoff value for homozygous deletions. Default is 0.5.
        a_consensus (float, optional): Coefficient for consensus high-level amplification cutoff. Default is 2.
        b_consensus (float, optional): Offset for consensus high-level amplification cutoff. Default is 0.5.
        a_remixt (float, optional): Coefficient for remixt high-level amplification cutoff. Default is 2.
        b_remixt (float, optional): Offset for remixt high-level amplification cutoff. Default is 0.5.
        model_hlamp_cutoff_offset (float, optional): Offset for model high-level amplification cutoff. Default is -1.
        call_cn_amp_cn (float, optional): Absolute copy number threshold to call amplifications. Default is 6.0.
        min_ploidy_cutoff (float, optional): Minimum ploidy cutoff to consider. Default is 2.
        min_gene_cn (float, optional): Minimum gene copy number to consider. Default is 2.

    Returns:
        dict: Dictionary with keys 'hlamp' and 'homdel', each containing a numpy array indicating the presence of the respective CNV events.
    """

    event_types = ['hlamp', 'homdel']
    total_cn = pd.DataFrame(gdata.layers['total_cn'], index=gdata.obs.index, columns=gdata.var.index)
    total_cn_remixt = pd.DataFrame(rgdata.layers['total_cn'], index=rgdata.obs.index, columns=rgdata.var.index)
    cn_events = {}
    for event_type in event_types:
        cn_events[event_type] = np.zeros(gdata.shape)
    
    for ix, sample in tqdm.tqdm(enumerate(gdata.obs.index), total=gdata.shape[0]): # per sample
        ploidy_consensus = gdata.obs.loc[sample, 'ploidy']
        hlamp_cutoff_consensus = (a_consensus * ploidy_consensus + b_consensus)
        hlamp_cutoff_remixt = hlamp_cutoff_consensus
        has_remixt_data = sample in rgdata.obs.index
        if has_remixt_data:
            ploidy_remixt = rgdata.obs.loc[sample, 'ploidy']
            hlamp_cutoff_remixt = (a_remixt * ploidy_remixt + b_remixt)
        
        sample_type = gdata.obs.loc[sample, 'sample_type']
        if sample_type.startswith('model'):
            hlamp_cutoff_consensus += model_hlamp_cutoff_offset
            hlamp_cutoff_remixt += model_hlamp_cutoff_offset
        
        for jx, gene in enumerate(gdata.var.index): # per gene
            gene_cn = total_cn.loc[sample, gene]
            gene_cn_remixt = np.nan
            if has_remixt_data:
                gene_cn_remixt = total_cn_remixt.loc[sample, gene]

            event = eval_cn_event(gene_cn, gene_cn_remixt, ploidy_consensus, ploidy_remixt,
                                  hlamp_cutoff_consensus, hlamp_cutoff_remixt,
                                  min_gene_cn=min_gene_cn, call_cn_amp_cn=call_cn_amp_cn, min_ploidy_cutoff=min_ploidy_cutoff, homdel_cutoff=homdel_cutoff)
            if event in event_types:
                cn_events[event][ix, jx] = 1
    return cn_events

def get_cnv_events_tcga(gdata, min_gene_cn_all=2.0, call_cn_amp_cn=6.0,
                  homdel_cutoff=0.5, a_cn=2, b_cn=0.5, min_ploidy_cutoff=2):
    """
    Get presenve of CNV events based on Consensus and ReMixT gene-level CN
    
    Arguments
    =========
    gdata: Consensus-based AnnData with obs as samples, var as genes
    tt2genes (optional): tumor type -> genes map

    Returns
    =======
    cn_events: dict with 'homdel', 'hlamp' keys, values being pd.DataFrame of event booleans
    """
    event_types = ['hlamp', 'homdel']
    total_cn = pd.DataFrame(gdata.layers['total_cn'], index=gdata.obs.index, columns=gdata.var.index)
    cn_events = {}
    for event_type in event_types:
        cn_events[event_type] = np.zeros(gdata.shape)
    
    for ix, sample in tqdm.tqdm(enumerate(gdata.obs.index), total=gdata.shape[0]):
        ploidy = gdata.obs.loc[sample, 'ploidy']
        hlamp_cutoff = (a_cn * ploidy + b_cn)
        
        for jx, gene in enumerate(gdata.var.index):
            flags = {e: False for e in event_types}
            gene_in_del_gene_set = True
            gene_cn = total_cn.loc[sample, gene]
            
            if gene_cn >= min_gene_cn_all:
                if (gene_cn >= hlamp_cutoff or gene_cn >= call_cn_amp_cn) and gene_cn >= ploidy + 1:
                    flags['hlamp'] = True
                if True: # else
                    if flags['hlamp'] and ploidy > min_ploidy_cutoff:
                        cn_events['hlamp'][ix, jx] = 1
    
            if gene_in_del_gene_set:
                if gene_cn < homdel_cutoff:
                    cn_events['homdel'][ix, jx] = 1
    return cn_events


def load_consensus_cn(sample_id, cn_dir = '../data/wgs/cn/consensus'):
    """
    Input: `sample_id` and `cn_dir`
    Return: pd.DataFrame of segtab.txt, or None if file nonexistent
    """
    cn_path = f'{cn_dir}/{sample_id}.segtab.txt'
    if not os.path.exists(cn_path):
        cn_paths = glob.glob(f'{cn_dir}/{sample_id}*.segtab.txt')
        assert len(cn_paths) == 1, cn_paths
        cn_path = cn_paths[0]
        
    if os.path.exists(cn_path):
        cn = pd.read_table(cn_path)
        cn.rename(columns={'Chromosome':'chromosome', 'Start.bp':'start', 'End.bp':'end'}, inplace=True)
        cn['chromosome'] = 'chr'+cn['chromosome'].astype(str).str.replace('chr', '')
        return cn
    return None
def load_remixt_cn(sample_id, remixt, rebin=False, bin_size=50000):
    if sample_id not in remixt.obs.index:
        raise MissingDataError()
    seg_cols = [
        'major_readcount', 'minor_readcount', 
        'major_0', 'minor_0',
        'major_1', 'minor_1', 'major_2', 'minor_2', 'major_raw', 'minor_raw',
        'major_depth_e', 'minor_depth_e', 'total_depth_e', 'major_e', 'minor_e',
        'total_e', 'major_raw_e', 'minor_raw_e', 'major_diff', 'minor_diff',
        'total_likelihood_mask', 'allele_likelihood_mask',
    ]
    case = sample_id[:17]
    data = remixt.var.reset_index(drop=True) 
    data['case'] = case
    data['sample'] = sample_id
    for layer in seg_cols:
        # print(layer)
        layer_vec = remixt[remixt.obs.index == sample_id, :].layers[layer].squeeze().copy()
        data[layer] = layer_vec
    data['segment_length'] = data['end'] - data['start'] + 1
    if rebin:
        assert type(bin_size) == int, bin_size
        data = data.rename(columns={'chromosome': 'chr'})
        data = rebin(data, bin_size, seg_cols)
        data = data.set_index(['chr', 'start', 'end'])
    data['major_subclonal'] = data['major_1'].round(0) != data['major_2'].round(0)
    data['minor_subclonal'] = data['minor_1'].round(0) != data['minor_2'].round(0)
    data['is_subclonal'] = data['major_subclonal'] | data['minor_subclonal']
    data['major_ccf'] = (data['major_raw_e'] - data['major_raw_e'].round()).abs()
    data['minor_ccf'] = (data['minor_raw_e'] - data['minor_raw_e'].round()).abs()
    data['total_raw'] = data['major_raw'] + data['minor_raw']
    
    return data

class MissingDataError(Exception):
    pass