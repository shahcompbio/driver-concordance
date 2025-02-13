import pickle

import yaml
import anndata as ad
import pandas as pd
import click

from gdan_hcmi_tools.process import (get_gene_level_anndata,
    get_gene_annotation_from_refflat, remove_events_not_curated)
from gdan_hcmi_tools.copynumber import get_cnv_events_hcmi

def annotate_cnv_events(gdata, rgdata, cnv_genes):
    """
    Annotates copy number variation (CNV) events in two AnnData objects, filtering based on curated CNV genes 
    and adding high-level amplification and homozygous deletion events.

    Parameters
    ----------
    gdata : AnnData
        AnnData object containing genomic or gene expression data for CNV annotation. CNV events are 
        filtered and joined to the observations in this object.
    rgdata : AnnData
        ReMixT AnnData object used in conjunction with `gdata` for CNV annotation, similarly filtered 
        and annotated based on curated CNV genes.
    cnv_genes : list of str
        List of curated genes to retain in the CNV events annotation process, removing any non-curated 
        CNV events.

    Returns
    -------
    None
        Modifies `gdata` and `rgdata` in place by:
        - Filtering and annotating CNV events based on `cnv_genes`.
        - Adding CNV event tables for high-level amplifications ('hlamp') and homozygous deletions 
          ('homdel') to `gdata` and `rgdata`.
    """
    remove_events_not_curated(gdata, cnv_genes)
    remove_events_not_curated(rgdata, cnv_genes)
    cn_events = get_cnv_events_hcmi(gdata, rgdata, gene_sets=None,
        homdel_cutoff=0.5,
        a_consensus=2, b_consensus=0, 
        a_remixt=2, b_remixt=0, 
        model_hlamp_cutoff_offset = 0)
    join_cnv_events_tables_to_obs(gdata, cn_events)
    remove_events_not_curated(gdata, cnv_genes)
    remove_events_not_curated(rgdata, cnv_genes)

def join_cnv_events_tables_to_obs(gdata, cn_events):
    """
    Joins copy number variation (CNV) event tables to the `obs` attribute of an AnnData object, 
    allowing for integration of amplification and deletion events into the main observation data.

    Parameters
    ----------
    gdata : AnnData
        An AnnData object containing gene expression or other genomic data, where CNV event layers 
        will be added and observed events will be joined with `obs`.
    cn_events : dict
        A dictionary containing CNV event matrices for different event types. The keys should 
        be 'hlamp' (high-level amplification) and 'homdel' (homozygous deletion), and the values 
        should be matrices of CNV event data.
    """
    event_types = ['hlamp', 'homdel']
    event_short_map = {'hlamp':'amp', 'homdel':'del'}
    
    for event_type in event_types:
        gdata.layers[event_type] = cn_events[event_type]
    
    for event_type in event_types:
        event_short = event_short_map[event_type]
        event_cols = [f'{g}_{event_short}' for g in gdata.var.index]
        if len(set(gdata.obs.columns) & set(event_cols)) == len(set(event_cols)):
            gdata.obs = gdata.obs.drop(event_cols, axis=1)
        cnv_event_table = pd.DataFrame(gdata.layers[event_type], index=gdata.obs.index, columns=event_cols)
        gdata.obs = gdata.obs.join(cnv_event_table)

def make_portable_event_table(gdata):
    """
    Creates a sharable event table from genomic data, extracting gene amplification and deletion events 
    with additional sample details.

    Parameters
    ----------
    gdata : AnnData
        Annotated data matrix (`AnnData` object) containing observation data with event columns. 
        Event columns in `gdata.obs` should be named in the format 'gene_amp' or 'gene_del'.

    Returns
    -------
    DataFrame
        A DataFrame with columns:
        - 'sample': Sample ID for each event.
        - 'sample_type': Type of sample (e.g., tumor, normal).
        - 'gene': Gene associated with the event.
        - 'event': Type of event, either 'homdel' for deletion or 'hlamp' for amplification.
        - 'case': Case ID associated with the sample.
        - 'cancer_type': Cancer type for the sample.
    """
    data = []
    event_map = {'del':'homdel', 'amp':'hlamp'}
    event_cols = [c for c in gdata.obs.columns if (c.endswith('_amp') or c.endswith('_del'))]
    for sample_id, row in gdata.obs.iterrows():
        sample_type = row['sample_type']
        case_id = row['case']
        cancer_type = row['cancer_type']
        if row[event_cols].sum():
            sample_events = {
                c.split('_')[0]:c.split('_')[1] for c in event_cols
                if row[c] > 0
            }
        else: # skip if no events
            continue
        for gene, event in sample_events.items():
            event = event_map[event]
            field = [sample_id, sample_type, gene, event, case_id, cancer_type]
            data.append(field)
    events_df = pd.DataFrame(data, columns=['sample', 'sample_type', 'gene', 'event', 'case', 'cancer_type'])
    return events_df

@click.command()
@click.option('-ic', '--input', help='input consensus h5ad with obs', required=True)
@click.option('-ir', '--input_remixt', help='input remixt h5ad with obs', required=True)
@click.option('-d', '--driver_genes', help='input CNV driver genes', required=True)
@click.option('-r', '--refflat', help='input refflat txt path', required=True)
@click.option('-c', '--cohort', help='input cohort data path', required=True)
@click.option('-og', '--output', help="output h5ad path", required=True)
@click.option('-or', '--output_remixt', help="output remixt h5ad path", required=True)
@click.option('-oe', '--output_events', help="output events table path", required=True)
def make_genelevel_anndata(input, input_remixt, driver_genes, refflat, cohort, output, output_remixt, output_events):
    adata = ad.read_h5ad(input)
    rdata = ad.read_h5ad(input_remixt)
    with open(cohort, 'rb') as handle:
        cohort = pickle.load(handle)
    # cohort.sample_ids = ['HCM-BROD-0002-C71-01A', 'HCM-BROD-0002-C71-86A', 'HCM-BROD-0002-C71-85A'] ##@##
    adata = adata[adata.obs.index.isin(cohort.sample_ids)] # remove samples with no tumor-model pairing
    rdata = rdata[rdata.obs.index.isin(cohort.sample_ids)] # remove samples with no tumor-model pairing
    cnv_genes = yaml.safe_load(open(driver_genes).read())
    refflat_cols = ['gene_symbol', 'gene_id', 'chrom', 'strand', 'cdna_start', 'cdna_end', 'orf_start', 'orf_end', 'exons', 'exon_starts', 'exon_ends']
    refflat = pd.read_table(refflat, names=refflat_cols)
    gene_var = get_gene_annotation_from_refflat(cnv_genes, refflat)
    gdata = get_gene_level_anndata(adata, gene_var)
    rgdata = get_gene_level_anndata(rdata, gene_var)
    annotate_cnv_events(gdata, rgdata, cnv_genes)
    events_df = make_portable_event_table(gdata)
    events_df.to_csv(output_events, index=False)
    rgdata.write_h5ad(output_remixt)
    gdata.write_h5ad(output)

if __name__ == "__main__":
    make_genelevel_anndata()
