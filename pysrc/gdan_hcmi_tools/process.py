import re
from operator import xor
from collections import defaultdict

import tqdm
import pandas as pd
import numpy as np
import anndata as ad

class Cohort:
    event_map = {
        'missense': 1,
        'truncating': 8,
        'hlamp': 16,
        'homdel': 32,
        'deletion': 4,
        'insertion': 2,
    }
    
    def __init__(self, links):
        self._links = links
        self._tumor_ids = []
        self._model_ids = []
        self._tumor2models_raw = {}
        self._model2tumor_raw = {}
        self._tumor2models = {}
        self._model2tumor = {}
        self._sample2sampletype = {}

        # Initialize data
        self.init_mapping()

    def init_mapping(self):
        tumor2models = {}
        model2tumor = {}

        for _, row in self._links.iterrows():
            tumor_id = row['tumor']
            model_id = row['model']
            model_expanded_id = row['model_expanded']

            if tumor_id == '': # no tumor -> skip
                print(f"Skipping {model_id} {model_expanded_id} d/t missing tumor ({tumor_id})")
                continue

            if tumor_id not in tumor2models:
                tumor2models[tumor_id] = []

            for idx, model in enumerate([model_id, model_expanded_id]):
                is_expanded = idx == 1
                tumor_type = 'tumor'
                model_type = 'model_expanded' if is_expanded else 'model'

                if model != '': # only pair existing model with tumor
                    tumor2models[tumor_id].append(model)
                    model2tumor[model] = tumor_id
                    self._sample2sampletype[model] = model_type
                    self._sample2sampletype[tumor_id] = tumor_type

        self._tumor2models = tumor2models
        self._tumor2models_raw = tumor2models.copy()
        self._model2tumor = model2tumor
        self._model2tumor_raw = model2tumor.copy()

        self._tumor_ids = list(tumor2models.keys()) # init tumor list
        self._model_ids = list(model2tumor.keys()) # init model list

    def _update_mappings(self): # update mapping dicts
        self._model2tumor = {
            model_id: tumor_id
            for tumor_id in self._tumor2models_raw
            for model_id in self._tumor2models_raw[tumor_id]
            if model_id in self._model_ids
        }
        self._tumor2models = {
            tumor_id: [
                model_id
                for model_id in self._tumor2models_raw[tumor_id]
                if model_id in self._model_ids
            ]
            for tumor_id in self._tumor_ids
        }

    def _update_ids(self, tumor_ids=None, model_ids=None): # update logic
        if tumor_ids is not None: # update tumors, re-enumerate models
            self._tumor_ids = list(set(tumor_ids))
            self._model_ids = [
                model_id
                for tumor_id in self._tumor_ids
                for model_id in self._tumor2models_raw.get(tumor_id, [])
            ]
            self._update_mappings()

        if model_ids is not None: # update models, re-enumerate tumors
            self._model_ids = list(set(model_ids))
            self._tumor_ids = list({
                self._model2tumor_raw[model_id]
                for model_id in self._model_ids
                if model_id in self._model2tumor_raw
            })
            self._update_mappings()

    @property
    def tumor_ids(self):
        return self._tumor_ids
    @tumor_ids.setter
    def tumor_ids(self, value):
        self._update_ids(tumor_ids=value)

    @property
    def model_ids(self):
        return self._model_ids
    @model_ids.setter
    def model_ids(self, value):
        self._update_ids(model_ids=value)

    @property
    def tumor2models(self):
        return self._tumor2models

    @property
    def model2tumor(self):
        return self._model2tumor
    
    @property
    def sample2sampletype(self):
        return self._sample2sampletype

    @property
    def sample_ids(self):
        assert isinstance(self._tumor_ids, list), self._tumor_ids
        assert isinstance(self._model_ids, list), self._model_ids
        return self._tumor_ids + self._model_ids
        
    def add_purity_ploidy(self, pp_path):
        pp = pd.read_table(pp_path)
        pp = pp[['id3', 'sample_id', 'consensus_purity', 'consensus_ploidy']]
        pp.rename(columns={'id3':'case_id'}, inplace=True)
        self.sample2purity = dict(zip(pp['sample_id'], pp['consensus_purity']))
        self.sample2ploidy = dict(zip(pp['sample_id'], pp['consensus_ploidy']))

    def add_cancer_type(self, ct_path, others_name='Others'):
        ct = pd.read_table(ct_path)
        ct = ct[['Case_ID', 'Diagnosis_TCGA_HCMI_Code']]
        ct.columns = ['case', 'cancer_type']
        ct.loc[ct['cancer_type'].isna(), 'cancer_type'] = others_name
        self.case2cancertype = dict(zip(ct['case'], ct['cancer_type']))
        self.sample2cancertype = self.assign_cancer_type_to_samples()

    def assign_cancer_type_to_samples(self, others_name='Others'):
        assert isinstance(self.case2cancertype, dict), self.case2cancertype
        sample2cancertype = {}
        for sample_id in self.sample_ids:
            case_id = sample_id.rsplit('-', 1)[0]
            if case_id in self.case2cancertype:
                sample2cancertype[sample_id] = self.case2cancertype[case_id]
            else:
                sample2cancertype[sample_id] = others_name
        return sample2cancertype


class VEP_meta:
    vep_fmts = ('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|'
        'HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|'
        'Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID'
    ).split('|')
    fmt_size = len(vep_fmts)
    consequences_map = { # from https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
        "transcript_ablation": 'T',
        "splice_acceptor_variant": 'T',
        "splice_donor_variant": 'T',
        "stop_gained": 'T',
        "frameshift_variant": 'T',
        "stop_lost": 'T',
        "start_lost": 'T',
        "transcript_amplification": 'T',
        "feature_elongation": 'T',
        "feature_truncation": 'T',
        "inframe_insertion": 'I',
        "inframe_deletion": 'D',
        "missense_variant": 'M',
        'upstream_gene_variant': 'P',
    }
    layer_map = {'T':'truncating', 'I':'insertion', 'D':'deletion', 'M':'missense', 'P':'promoter'}
    impact_selected = {'MODIFIER', 'MODERATE', 'HIGH'}
    event_rank = {
        'homdel': 7,
        'hlamp': 6,
        'truncating': 5,
        'promoter': 4,
        'deletion': 3,
        'insertion': 2,
        'missense': 1,
    }

class VEP_germline:
    vep_fmts = ('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|'
                'HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|'
                'Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|'
                'TSL|APPRIS|GIVEN_REF|USED_REF|SOURCE|GENE_PHENO|NEAREST|SIFT|PolyPhen|DOMAINS|'
                'HGVS_OFFSET|AF_1000G|AFR_AF_1000G|AMR_AF_1000G|EAS_AF_1000G|EUR_AF_1000G|SAS_AF_1000G|'
                'CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|'
                'CADD_PHRED|CADD_RAW|FATHMM|FATHMM_SOMATIC|'
                'CosmicCoding|CosmicCoding_CNT|CosmicCoding_CDS|CosmicCoding_AA|CosmicNonCoding|'
                'CLN|CLN_VARIATION_ID|CLN_MOLECULAR_CONSEQUENCE|CLN_CLINICAL_SIGNIFICANCE|CLN_CONFLICTED|'
                'CLN_REVIEW_STATUS|CLN_TRAITS|CLN_PMIDS|CLN_XREFS|CLN_ORIGIN|'
                'GnomadExomes|GnomadExomes_AF|GnomadExomes_AN|GnomadExomes_Hom|'
                'GnomadExomes_AF_AFR|GnomadExomes_AF_AMR|GnomadExomes_AF_ASJ|GnomadExomes_AF_EAS|'
                'GnomadExomes_AF_FIN|GnomadExomes_AF_NFE|GnomadExomes_AF_OTH|'
                'GnomadGenomes|GnomadGenomes_AF|GnomadGenomes_AN|GnomadGenomes_Hom|GnomadGenomes_AF_AFR|'
                'GnomadGenomes_AF_AMR|GnomadGenomes_AF_ASJ|GnomadGenomes_AF_EAS|GnomadGenomes_AF_FIN|'
                'GnomadGenomes_AF_NFE|GnomadGenomes_AF_OTH|mtfl|mtfl_FUNC_LOC|'
                'mitimpact|mitimpact_OXPHOS_complex|mitomap_disease|mitomap_disease_AC|'
                'mitomap_disease_AF|mitomap_disease_homoplasmy|mitomap_disease_heteroplasmy|'
                'mitomap_disease_PubmedIDs|mitomap_disease_Disease|mitomap_disease_DiseaseStatus|'
                'mitomap_poly|mitomap_poly_AC|mitomap_poly_AF|ACMG59|ACMG59_GENE|ACMG59_DISEASE|'
                'CCRS|CCRS_ccr_pct|CLN|CLN_DiseaseName|phastcons100|phyloP100|AR|AR_AR_GENE|'
                'PGx|PGx_pgx_rsid|REVEL|REVEL_REVEL_SCORE').split('|')
    fmt_size = len(vep_fmts)
    consequences_map = { # from https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
        "transcript_ablation": 'T',
        "splice_acceptor_variant": 'T',
        "splice_donor_variant": 'T',
        "splice_region_variant": 'T',
        "stop_gained": 'T',
        "frameshift_variant": 'T',
        "stop_lost": 'T',
        "start_lost": 'T',
        "transcript_amplification": 'T',
        "feature_elongation": 'T',
        "feature_truncation": 'T',
        "inframe_insertion": 'I',
        "inframe_deletion": 'D',
        "missense_variant": 'M',
        "synonymous_variant": 'S',
    }
    flag_map = {
        0: 'not found',
        1: 'oncogene in database',
        2: 'truncating variant',
        3: 'oncogene in database, truncating variant',
        4: 'known variant',
        5: 'oncogene in database, known variant',
        6: 'truncating variant, known variant',
        7: 'oncogene in database, truncating variant, known variant',
        8: 'tumor suppressor in database',
        10: 'truncating variant, tumor suppressor in database',
        12: 'known variant, tumor suppressor in database',
        14: 'truncating variant, known variant, tumor suppressor in database',
    }
    flags_to_keep = { # 4, 5, 6, 7, 10, 12, 14
        'known variant',
        'oncogene in database, known variant',
        'truncating variant, known variant',
        'oncogene in database, truncating variant, known variant',
        'truncating variant, tumor suppressor in database',
        'known variant, tumor suppressor in database',
        'truncating variant, known variant, tumor suppressor in database',
    }


def read_tm_link(tm_link_filename):
    tm_link = pd.read_table(tm_link_filename, sep='\t', low_memory=False)
    tm_link = tm_link.rename(columns={
        'Case ID3': 'case',
        'Model ID3': 'unused',
        'Cancer Model Aliquot': 'model',
        'Matched Tumor Aliquot': 'tumor',
        'Matched Normal Aliquot': 'normal',
        'Matched Expanded Model Aliquot': 'model_expanded',
    })
    tm_link = tm_link.fillna('') # avoid using isinstance for str

    return tm_link


def generate_model_table(tm_link):
    """ Create model centric table

    Parameters
    ----------
    tm_link : DataFrame
        tumor model linkage table with columns case, tumor, model, model_expanded

    Returns
    -------
    DataFrame : table with columns case, model, tumor, model_type
    """

    model_table = []

    for _, row in tm_link.iterrows():
        tumor_id = row['tumor']
        case = row['case']
        if not tumor_id:
            continue
        for model_type in ('model', 'model_expanded'):
            model_id = row[model_type]
            if not model_id:
                continue
            model_table.append({
                'case': case,
                'model': model_id,
                'tumor': tumor_id,
                'model_type': model_type,
            })

    model_table = pd.DataFrame(model_table)

    return model_table

def interp_nan(vector):
    base_value = 0.
    vector = vector.copy()
    mask = np.isnan(vector)
    if np.all(mask):
        vector[mask] = base_value
    else:
        vector[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), vector[~mask])
    return vector

def interp_nan_per_chrom(data, chrom_mask):
    data = data.copy()
    if np.all(np.isnan(data[:, chrom_mask])):
        data[:, chrom_mask] = 0.
    else:
        data[:, chrom_mask] = np.apply_along_axis(interp_nan, 1, data[:, chrom_mask])
    return data

def get_gene_annotation_from_refflat(cnv_genes, refflat):
    canonical_chroms_short = [str(c) for c in range(1, 23)] + ['X', 'Y']
    canonical_chroms = ['chr'+c for c in canonical_chroms_short]
    data = []
    genes_saved = set()
    for tt, tt_events in cnv_genes.items():
        for event, genes in tt_events.items():
            for gene in genes:
                if gene in genes_saved: 
                    continue
                else:
                    genes_saved.add(gene)
                assert gene in set(refflat['gene_symbol'].tolist()), (tt, event, gene)
                refflat_genes = refflat[(refflat['gene_symbol']==gene) & (refflat['chrom'].isin(canonical_chroms))]
                chroms = refflat_genes['chrom'].unique()
                assert chroms.shape[0] == 1, chroms
                chrom = chroms[0]
                start = refflat_genes['cdna_start'].min()
                end = refflat_genes['cdna_end'].max()
                field = [gene, chrom, start, end]
                data.append(field)
    gene_var = pd.DataFrame(data, columns=['gene', 'chrom', 'start', 'end'])
    gene_var['chrom'] = pd.Categorical(gene_var['chrom'], categories=canonical_chroms, ordered=True)
    gene_var.sort_values(['chrom', 'start', 'end'], inplace=True)
    return gene_var

def make_curated_gene_anndata(gene_var, adata):
    var = adata.var
    gene_cols = ['gene', 'chr', 'start', 'end']
    gene_table = pd.DataFrame(columns=gene_cols)
    data_dict = defaultdict(list)
    for riw, row in gene_var.iterrows():
        gene = row['gene']
        chrom = str(row['chrom']).replace('chr', '')
        start = int(row['start'])
        end = int(row['end'])
        gene_field = [gene, chrom, start, end]
        df = var[
            (var['chr']==chrom) &
            (var['start'] < end) &
            (var['end'] >= start)
        ]
        if df.shape[0] == 0:
            continue
        gene_table.loc[gene_table.shape[0]] = gene_field
        for layer_name in adata.layers:
            data_vec = adata[:, df.index].layers[layer_name].mean(axis=1)
            data_dict[layer_name].append(data_vec)

    for layer_name in data_dict:
         data_dict[layer_name] = np.array(data_dict[layer_name])
    gene_table.set_index('gene', inplace=True)
    gdata = ad.AnnData(var=adata.obs, obs=gene_table, layers=data_dict).T
    return gdata


def remove_events_not_curated(gdata, cnv_genes):
    obs = gdata.obs.copy()
    event_cols = [c for c in obs.columns if (c.endswith('_amp') or c.endswith('_del'))]
    for tt, tt_obs in obs.groupby('cancer_type'):
        sample_ids = tt_obs.index
        if tt in cnv_genes:
            amp_genes = cnv_genes[tt]['hlamp']
            del_genes = cnv_genes[tt]['homdel']
            amp_cols = set([f'{g}_amp' for g in amp_genes])
            del_cols = set([f'{g}_del' for g in del_genes])
            for col in event_cols:
                if col in amp_cols or col in del_cols:
                    continue
                else:
                    obs.loc[sample_ids, col] = 0
    gdata.obs = obs


def get_gene_level_anndata(adata, gene_var):
    layer_name = 'total_cn'
    for chrom in adata.var['chr'].unique():
        chrom_mask = (adata.var['chr']==chrom)
        adata.layers[layer_name] = interp_nan_per_chrom(adata.layers[layer_name], chrom_mask)

    gdata = make_curated_gene_anndata(gene_var, adata)
    gdata = gdata[gdata.obs.dropna().index] # remove samples with nulls in annotation
    return gdata


def extract_sample_ids(pair_id):
    assert '--' in pair_id
    sample1, sample2 = pair_id.split('--')
    tags1 = sample1.split('-')
    tags2 = sample2.split('-')
    short1 = '-'.join(tags1[:5])
    short2 = '-'.join(tags2[:5])
    is_normal1 = tags1[4][0] == '1'
    is_normal2 = tags2[4][0] == '1'
    assert xor(is_normal1, is_normal2)
    if is_normal1:
        return (sample2, sample1, short2, short1)
    else:
        return (sample1, sample2, short1, short2)

def get_samples_to_exclude(qc):
    df = pd.read_csv(qc)
    assert 'model_id' in df.columns, df
    # samples_to_exclude = df['model_id'].tolist()
    samples_to_exclude = [] # don't remove samples for genomics
    return samples_to_exclude

def read_vcf(vcf_path, short_id_length=21):
    header = ['chrom', 'pos', 'ID', 'ref', 'alt', 'qual', 'FILTER', 'info', 'FORMAT']
    header_cnt = 0
    for line in open(vcf_path):
        if line.startswith('#CHROM'): # occurs only once
            header_cnt += 1
            field = line.rstrip().split()
            sample_cols = field[9:] # VCF sample columns
            sample_cols = [s[:short_id_length] for s in sample_cols] # shortened ID
            header += sample_cols
    assert header_cnt == 1, header_cnt
    df = pd.read_table(vcf_path, comment='#', names=header)
    return df

def extract_vaf(row, sample_id):
    fmt = row['FORMAT'].split(':')
    gt = row[sample_id].split(':')
    gt = dict(zip(fmt, gt))
    assert 'AD' in gt, f'ERROR: AD key not in genotype of {sample_id} for:\n{row}'
    ads = [int(d) if d != '.' else 0 for d in gt['AD'].split(',')]
    depth = sum(ads)
    var_count = depth - ads[0] # subtract ref count
    vaf = var_count / depth if depth else 0.0
    return vaf


def make_obs_var(data):
    """
    Constructs observation and variable tables from a dictionary of data, representing genomic copy number 
    or segmentation information across multiple samples.

    Parameters
    ----------
    data : dict
        Dictionary where keys are sample IDs and values are DataFrames with copy number or genomic 
        segmentation data. Each DataFrame should have 'chr', 'start', and 'end' columns to define 
        genomic regions.

    Returns
    -------
    tuple
        A tuple containing:
        - `data` (DataFrame): Concatenated DataFrame of all sample data, with each sample as a column.
        - `obs` (DataFrame): Observation DataFrame with an index of sample IDs.
        - `var` (DataFrame): Variable DataFrame with genomic region bins as the index, structured 
          with 'chr', 'start', and 'end' columns to represent each genomic segment.
    """
    sample_ids = list(data.keys())
    data = pd.concat(data, axis=1)

    obs = pd.DataFrame(index=sample_ids)

    var = data.reset_index()[['chr', 'start', 'end']]
    var['bin'] = var['chr'].astype(str) + ':' + var['start'].astype(str) + '-' + var['end'].astype(str)
    var.columns = var.columns.get_level_values(0)
    var.set_index('bin', inplace=True)

    return data, obs, var

def make_anndata(data, obs, var, seg_cols):
    """
    Creates an AnnData object from genomic data layers, observations, and variables.

    Parameters
    ----------
    data : DataFrame
        A DataFrame containing genomic data for multiple samples, with multi-level columns where 
        the second level contains layer names defined in `seg_cols`.
    obs : DataFrame
        DataFrame representing observations (samples) with sample IDs as the index.
    var : DataFrame
        DataFrame representing variables (genomic regions) with region bins as the index and columns 
        'chr', 'start', and 'end' to define genomic positions.
    seg_cols : list of str
        List of column names in `data` that represent different layers of genomic data to be added 
        to the AnnData object.

    Returns
    -------
    AnnData
        An AnnData object with `obs` as observation data, `var` as variable data, and `layers` 
        representing different layers of genomic segmentation information for each sample.
    """
    layers = {}
    for layer in tqdm.tqdm(seg_cols):
        layer_data = data.loc[:, (slice(None), layer)]
        layer_data.columns = layer_data.columns.droplevel(1)
        layer_data = layer_data.reset_index()
        layer_data['bin'] = layer_data['chr'].astype(str) + ':' + layer_data['start'].astype(str) + '-' + layer_data['end'].astype(str)
        layer_data = layer_data.set_index('bin')
        layer_data = layer_data.loc[var.index, obs.index].T
        layers[layer] = layer_data # index: sample, column: region

    var['chr'] = var.index.str.split(':').str[0]
    chrom_categories = [str(c) for c in range(1, 23)] + ['X', 'Y']
    var['chr'] = pd.Categorical(var['chr'], categories=chrom_categories, ordered=True)

    # Check indices match
    for layer_name in layers:
        assert np.all(layers[layer_name].index == obs.index)
        assert np.all(layers[layer_name].columns == var.index)

    adata = ad.AnnData(var=var, obs=obs, layers=layers)

    # Sort anndata, be sure to do this after the anndata creation
    sorted_index = adata.var.sort_values(['chr', 'start', 'end']).index
    adata = adata[:, sorted_index]

    return adata

def transform_aa_change(aa_change):
    three2one = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Ter": "*",
    }
    pattern = r'[A-Z][a-z]{2}' # e.g. Ala
    result = aa_change
    for search in re.finditer(pattern, aa_change):
        src = search.group()
        if src in three2one:
            dst = three2one[src]
            result = result.replace(src, dst, 1)
    return result