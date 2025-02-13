import click
import pandas as pd
import numpy as np
import anndata as ad

import gdan_hcmi_tools.copynumber

pcawg_tissue_map = {
    'Liver-HCC': 'LIHC',
    'Prost-AdenoCa': 'PRAD',
    'Panc-AdenoCa': 'PAAD',
    'Breast-AdenoCa': 'BRCA',
    'CNS-Medullo': 'MB', #
    'Ovary-AdenoCa': 'OV', 
    'Kidney-CCRCC': 'KIRC', 
    'Lymph-BNHL': 'NHL', #
    'Skin-Melanoma': 'SKCM', 
    'Eso-AdenoCa': 'ESCA', 
    'Lymph-CLL': 'CLL',  #
    'CNS-PiloAstro': 'LGG',
    'Panc-Endocrine': 'PAEN',  #
    'Stomach-AdenoCa': 'STAD', 
    'ColoRect-AdenoCa': 'COAD', 
    'Head-SCC': 'HNSC',
    'Myeloid-MPN': 'MPN',  #
    'Uterus-AdenoCa': 'UCEC', 
    'Thy-AdenoCa': 'THCA', 
    'Lung-SCC': 'LUSC',
    'Kidney-ChRCC': 'KICH', 
    'CNS-GBM': 'GBM', 
    'Bone-Osteosarc': 'SARC', 
    'Lung-AdenoCa': 'LUAD',
    'Biliary-AdenoCa': 'CHOL', 
    'Kidney-PapRCC': 'KIRP', 
    'Bladder-TCC': 'BLCA',
    'SoftTissue-Liposarc': 'SARC', 
    'Cervix-SCC': 'CESC', 
    'CNS-Oligo': 'LGG', 
    'Bone-Benign': 'BB',  #
    'SoftTissue-Leiomyo': 'SARC', 
    'Breast-LobularCa': 'BRCA', 
    'Myeloid-AML': 'LAML', 
    'Bone-Epith': 'BE', #
    'Myeloid-MDS': 'MDS', #
    'Breast-DCIS': 'BCIS',  #
    'Cervix-AdenoCa': 'CESC',
}

def rename_layers(adata, layer_name_map, delete_source=False):
    for src_name, dst_name in layer_name_map.items():
        if src_name in adata.layers:
            adata.layers[dst_name] = adata.layers[src_name].round(0)
            if delete_source:
                del(adata.layers[src_name])

def rename_layers_hcmi_consensus(adata):
    layer_name_map = {
        'rescaled.cn.a1': 'minor_cn',
        'rescaled.cn.a2': 'major_cn',
        'rescaled_total_cn': 'total_cn',
    }
    rename_layers(adata, layer_name_map, delete_source=True)

def rename_layers_hcmi_remixt(adata):
    layer_name_map = {
        'minor_raw': 'minor_cn',
        'major_raw': 'major_cn',
    }
    rename_layers(adata, layer_name_map, delete_source=False)
    adata.layers['total_cn'] = adata.layers['minor_cn'] + adata.layers['major_cn']
    adata.layers['LOH'] = (adata.layers['minor_cn'] == 0).astype(int)

def rename_layers_tcga(adata):
    layer_name_map = {
        'minor_copy_number': 'minor_cn',
        'major_copy_number': 'major_cn',
        'copy_number': 'total_cn',
    }
    rename_layers(adata, layer_name_map, delete_source=True)
    adata.layers['LOH'] = (adata.layers['minor_cn'] == 0).astype(int)

def rename_layers_pcawg(adata):
    layer_name_map = {
        'minor_cn': 'minor_cn',
        'major_cn': 'major_cn',
        'total_cn': 'total_cn',
    }
    rename_layers(adata, layer_name_map, delete_source=True)
    adata.layers['LOH'] = (adata.layers['minor_cn'].round() == 0).astype(int)

def read_and_normalize_adata(cn, source):
    """
    Rename layer names and add LOH if the layer does not exist.
    """
    adata = ad.read_h5ad(cn)
    adata.obs.index.name = 'sample_id'
    if source == 'HCMI-consensus':
        rename_layers_hcmi_consensus(adata)
    elif source == 'HCMI-remixt':
        rename_layers_hcmi_remixt(adata)
    elif source == 'TCGA':
        rename_layers_tcga(adata)
    elif source == 'PCAWG':
        rename_layers_pcawg(adata)
    else:
        raise ValueError(f'Wrong value for source={source}')
    dst_names = {'major_cn', 'minor_cn', 'total_cn'} # assert presence of these layers
    assert set(adata.layers.keys()) & dst_names == dst_names, adata.layers
    return adata

def process_metadata_from_analysis_tracker(metadata):
    """
    Process pre-saved Analysis Tracker spreadsheet and
    map sample types as 'tumor', 'model', or 'model_expanded'
    """
    detailed_sample_type_map = {
        'next generation cancer model': 'model',
        'primary tumor': 'tumor',
        'metastatic': 'tumor',
        'expanded next generation cancer model': 'model_expanded',
        'post neo-adjuvant therapy': 'tumor',
        'slides': 'tumor',
        'recurrent tumor': 'tumor',
        'human tumor original cells': 'tumor',
        'ffpe scrolls': 'tumor',
        'neoplasms of uncertain and unknown behavior': 'tumor',
        'additional metastatic': 'tumor',
        'ffpe recurrent': 'tumor',
    }
    meta = pd.read_csv(metadata)
    meta = meta[['id3', 'sample_id', 'sample_type', 'detailed_sample_type']]
    meta.rename(columns={'id3':'case'}, inplace=True)
    meta['detailed_sample_type'] = meta['detailed_sample_type'].str.lower().str.replace('_', ' ')
    meta['sample_type'] = meta['detailed_sample_type'].replace(detailed_sample_type_map)
    assert set(meta['sample_type'].value_counts(dropna=False).index) == {'model', 'tumor', 'model_expanded'}
    return meta


def add_consensus_purity_ploidy(adata, pp):
    """
    Parse purity and ploidy spreadsheet and add the info to anndata obs
    """
    pp = pd.read_table(pp)
    assert pp['sample_id'].value_counts().max() == 1, pp['sample_id'].value_counts()
    pp.rename(columns={'id3':'case'}, inplace=True)
    pp.drop('sample_type', axis=1, inplace=True)
    pp.set_index('sample_id', inplace=True)
    sample2purity = dict(zip(pp.index, pp['consensus_purity']))
    sample2ploidy = dict(zip(pp.index, pp['consensus_ploidy']))
    # adata.uns['consensus_purity_ploidy'] = pp
    adata.obs['consensus_purity'] = adata.obs.index.map(sample2purity) # no change in values but more explicit
    adata.obs['consensus_ploidy'] = adata.obs.index.map(sample2ploidy)

# def add_wgd_status(adata, use_consensus_ploidy=True): # DEPRECATED
#     assert 'loh' in adata.obs, adata.obs
#     assert 'ploidy' in adata.obs or 'consensus_ploidy' in adata.obs, adata.obs
#     loh = adata.obs['loh']
#     if use_consensus_ploidy and 'consensus_ploidy' in adata.obs:
#         ploidy = adata.obs['consensus_ploidy']
#     else:
#         ploidy = adata.obs['ploidy']
#     adata.obs['n_wgd'] = (2.9 -2*loh <= ploidy).astype(int)

def extract_tcga_sample_ids(tcga_metadata):
    """
    Process TCGA metadata to extract unique tumor sample identifiers and add them to the DataFrame.

    This function iterates through the TCGA metadata DataFrame to identify and save unique
    tumor sample identifiers based on file and sample IDs. It filters out normal samples and
    ensures each tumor sample is uniquely represented, adding a 'sample' column with tumor IDs.

    Args:
        tcga_metadata (DataFrame): DataFrame containing metadata for TCGA samples, including
                                   file IDs, file names, project IDs, and sample IDs.

    Returns:
        DataFrame: The updated DataFrame with an additional 'sample' column, containing
                   identifiers for tumor samples only, with non-tumor samples removed.
    """
    df = tcga_metadata.copy()
    df['sample_id'] = ''
    saved = set()
    for rix, row in df.iterrows():
        sample_ids = row['Sample ID'].split(', ')
        cancer_type = row['Project ID'].replace('TCGA-', '')
        sample_types = row['Sample Type'].split(', ') # e.g. "Primary Tumor, Blood Derived Normal"
        assert len(sample_types) == 2, sample_types
        tumor_ixs = [sample_types.index(s) for s in sample_types if s.lower().count('normal') == 0]
        assert len(tumor_ixs) == 1, tumor_ixs
        tumor_ix = tumor_ixs[0]
        tumor_id = sample_ids[tumor_ix]
        case_id = '-'.join(tumor_id.split('-')[:3])
        if tumor_id in saved: 
            # print(f'{tumor_id} already saved; from {sample_ids}', file=sys.stderr)
            continue
        saved.add(tumor_id)
        df.loc[rix, 'sample_id'] = tumor_id
        df.loc[rix, 'case'] = case_id
        df.loc[rix, 'cancer_type'] = cancer_type
    df = df[df['sample_id'].str.len() > 0]
    return df.reset_index(drop=True)

def map_metadata_to_obs_hcmi(adata, ct):
    assert 'Diagnosis_TCGA_HCMI_Code' in ct.columns, ct.columns # HCMI: TCGA_Codes.tsv (Dina ElHarouni)
    adata.obs['case'] = adata.obs.index.str.slice(0, 17)
    ct = ct[['Case_ID', 'Diagnosis_TCGA_HCMI_Code']].drop_duplicates()
    ct.rename(columns={'Case_ID':'case', 'Diagnosis_TCGA_HCMI_Code':'cancer_type'}, inplace=True)
    ct_map = dict(zip(ct['case'], ct['cancer_type']))
    adata.obs['cancer_type'] = adata.obs['case'].map(ct_map) # np.nan should be removed downstream
    adata.obs['cohort'] = 'HCMI'

def map_metadata_to_obs_tcga(adata, ct):
    assert 'Project ID' in ct.columns, ct.columns # TCGA: sample_sheet.tsv
    ct = extract_tcga_sample_ids(ct)
    ct_map = dict(zip(ct['case'], ct['cancer_type']))
    adata.obs['case'] = adata.obs.index.str.rsplit('-', n=1).str[0]
    adata.obs['cancer_type'] = adata.obs['case'].map(ct_map)
    adata.obs['cohort'] = 'TCGA'
    cols_to_add_nan = ['sample_type', 'detailed_sample_type', 'consensus_purity', 'consensus_purity']
    adata.obs[cols_to_add_nan] = np.nan

def map_metadata_to_obs_pcawg(adata, ct):
    assert 'icgc_sample_id' in ct.columns, ct.columns
    assert 'icgc_donor_id' in ct.columns, ct.columns
    assert 'WGD' in ct.columns, ct.columns
    assert 'cancer_type' in ct.columns, ct.columns
    ct.rename(columns={'icgc_sample_id':'sample_id', 'icgc_donor_id':'case', 'WGD':'n_wgd'}, inplace=True)
    ct['cancer_type'] = ct['tissue'].map(pcawg_tissue_map)
    ct['n_wgd'] = ct['wgd'].astype(int)
    ct = ct[['uuid', 'sample_id', 'case', 'cancer_type', 'ploidy', 'purity', 'n_wgd']]
    ct_keys = ['case', 'cancer_type', 'ploidy', 'purity', 'n_wgd']
    for ct_key in ct_keys: # map to guarentee preservation of column order
        ct_map = dict(zip(ct['sample_id'], ct[ct_key]))
        adata.obs[ct_key] = adata.obs.index.map(ct_map)
    cols_to_add_nan = ['sample_type', 'detailed_sample_type', 'consensus_purity', 'consensus_purity']
    adata.obs[cols_to_add_nan] = np.nan

def add_cancer_types_to_obs(adata, ct, source):
    """Map 'cancer_type' to adata obs

    Args:
        adata: sample x (gene or segment) anndata for addition of `obs`
        ct: dataframe with `sample_id` and `cancer_type` columns
        source: data source, either HCMI-consensus, HCMI-remixt, TCGA, or PCAWG
    """
    assert adata.obs.index.unique().shape[0] == adata.shape[0], adata.obs.index.value_counts() # assert unique sample IDs
    if ct.endswith('.csv') or ct.endswith('.csv.gz'):
        ct = pd.read_csv(ct)
    else:
        ct = pd.read_table(ct)

    if source.startswith('HCMI'):
        map_metadata_to_obs_hcmi(adata, ct)
    elif source == 'TCGA':
        map_metadata_to_obs_tcga(adata, ct)
    elif source == 'PCAWG':
        map_metadata_to_obs_pcawg(adata, ct)
    else:
        raise ValueError(f'Wrong value for source={source}')        

def add_metadata_to_obs_hcmi(adata, cancer_types, pairs, pp, metadata):
    """
    Adds metadata from the Human Cancer Models Initiative (HCMI) to the AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix of shape (n_obs, n_vars) where `n_obs` represents the number of observations.
    cancer_types : DataFrame
        DataFrame containing mappings of observation indices to specific cancer types.
    pairs : str or Path
        Path to the file containing pairs of case data, with columns 'case', 'id3', 'model', 
        'model_expanded', 'tumor', 'normal', and 'freeze1'. This file maps sample pairs and 
        is expected in a tab-separated format.
    pp : str or Path
        Path to the purity and ploidy metadata file, containing consensus purity and ploidy data.
    metadata : DataFrame
        DataFrame containing metadata from the analysis tracker. This includes `sample_id`, 
        `sample_type`, and `detailed_sample_type` columns that will be mapped to `adata.obs`.
    """
    assert pairs, f'pairs={pairs}'
    assert pp, f'pp (purity and purity metadata path)={pp}'
    assert metadata, f'metadata={metadata}'

    meta = process_metadata_from_analysis_tracker(metadata)
    sample_type_map = dict(zip(meta['sample_id'], meta['sample_type']))
    detailed_sample_type_map = dict(zip(meta['sample_id'], meta['detailed_sample_type']))

    pairs = pd.read_table(pairs, dtype=str)
    pairs = pairs.iloc[:, :7]
    pairs.columns = ['case', 'id3', 'model', 'model_expanded', 'tumor', 'normal', 'freeze1']
    pairs.fillna('', inplace=True)

    adata.uns['pairs'] = pairs
    adata.obs['sample_type'] = adata.obs.index.map(sample_type_map)
    adata.obs['detailed_sample_type'] = adata.obs.index.map(detailed_sample_type_map)
    add_cancer_types_to_obs(adata, cancer_types, 'HCMI')

    add_consensus_purity_ploidy(adata, pp)

@click.command()
@click.option('-i', '--cn', help='input copy number anndata', required=True)
@click.option('-p', '--pairs', help='input tumor model pairs table')
@click.option('-pp', '--pp', help='input consensus purity/ploidy table')
@click.option('-m', '--metadata', help="input metadata csv path (from Analysis Tracker: dna_analysis_updated)")
@click.option('-c', '--cancer_types', help='input cancer types table', required=True)
@click.option('-s', '--source', help='input data source label (e.g. TCGA)', required=True)
@click.option('-o', '--output', help="output h5ad path", required=True)
def add_obs(cn, pairs, pp, metadata, cancer_types, source, output):
    adata = read_and_normalize_adata(cn, source)
    assert source in {'HCMI-consensus', 'HCMI-remixt', 'TCGA', 'PCAWG'}, f'Wrong input for source={source}.'
    if source.startswith('HCMI'):
        add_metadata_to_obs_hcmi(adata, cancer_types, pairs, pp, metadata)
    else:
        add_cancer_types_to_obs(adata, cancer_types, source)

    gdan_hcmi_tools.copynumber.add_fraction_genome_altered(adata)
    gdan_hcmi_tools.copynumber.add_loh_fraction(adata)

    adata = gdan_hcmi_tools.copynumber.add_calculated_ploidy(adata)
    adata = gdan_hcmi_tools.copynumber.calculate_n_wgd(adata)

    adata.write_h5ad(output)

if __name__ == "__main__":
    add_obs()
