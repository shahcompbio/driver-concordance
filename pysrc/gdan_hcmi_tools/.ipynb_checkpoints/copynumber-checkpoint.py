import pandas as pd
import pyranges as pr
import anndata as ad
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import wgs_analysis.refgenome
from wgs_analysis.plots.cnv import plot_cnv_segments
import gdan_hcmi_tools.copynumber
import matplotlib


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

def load_consensus_cn(sample_id, cn_dir = '../data/wgs/cn/consensus'):
    """
    Input: `sample_id` and `cn_dir`
    Return: pd.DataFrame of segtab.txt, or None if file nonexistent
    """
    cn_path = f'{cn_dir}/{sample_id}.segtab.txt'
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
    # Filtering
    data['segment_length'] = data['end'] - data['start'] + 1
    # data['length_ratio'] = data['length'] / data['segment_length']
    # data = data.query('length_ratio > 0.8').copy()
    # data = data[data['minor_readcount'] > 100]

    if rebin:
        assert type(bin_size) == int, bin_size
        data = data.rename(columns={'chromosome': 'chr'})
        data = gdan_hcmi_tools.copynumber.rebin(data, bin_size, seg_cols)
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


def rebin(data:pd.DataFrame, bin_size:int, cols:list, 
          ref_cols:list=['Chromosome', 'Start', 'End'], 
          cn_cols:list=['Chromosome', 'Start.bp', 'End.bp'],
          genome_version:str='grch38') -> pd.DataFrame:
    """ 
    Merge/split data bins to given bin_size
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


