import yaml
import anndata as ad
import pandas as pd
import click

from gdan_hcmi_tools.process import (
    interp_nan_per_chrom, 
    make_curated_gene_anndata, 
    get_gene_annotation_from_refflat, 
    remove_events_not_curated
)
from gdan_hcmi_tools.copynumber import get_cnv_events_tcga


def join_cnv_events_tables_to_obs(gdata, cn_events):
    """
    Integrates copy number variation (CNV) event tables into the `obs` attribute of an AnnData object,
    adding amplification and deletion event columns to the main observation data.

    Parameters
    ----------
    gdata : AnnData
        An AnnData object containing genomic or gene expression data. This function will add CNV 
        event layers for amplification and deletion events and integrate these events into `gdata.obs`.
    cn_events : dict
        Dictionary containing CNV event matrices. Keys should be event types ('hlamp' for high-level
        amplifications and 'homdel' for homozygous deletions), with each value being a matrix of CNV 
        event data to be added to the respective layer in `gdata`.
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

def make_melted_events_from_layers(gdata):
    """
    Creates a melted table of CNV events from specific layers in an AnnData object, summarizing
    samples, genes, event types, and cancer types in a long-format DataFrame.

    Parameters
    ----------
    gdata : AnnData
        An AnnData object containing layers with CNV event data. The layers should include 'homdel' 
        (homozygous deletion) and 'hlamp' (high-level amplification) event data, indexed by sample ID.

    Returns
    -------
    DataFrame
        A melted DataFrame with columns:
        - 'sample_id': Sample ID for each event.
        - 'gene': Gene associated with the CNV event.
        - 'event': Type of event ('homdel' for deletions and 'hlamp' for amplifications).
        - 'cancer_type': Cancer type for each sample, sourced from `gdata.obs`.
    """
    data = []
    layer_names = ['homdel', 'hlamp']
    events_cols = ['sample_id', 'event']
    events = pd.DataFrame(columns=events_cols)
    for layer_name in layer_names:
        df = pd.DataFrame(gdata.layers[layer_name], index=gdata.obs.index, columns=gdata.var.index)
        for gene in df.columns:
            samples_with_event = df[gene][df[gene] > 0].index.tolist()
            for sample_id in samples_with_event:
                cancer_type = gdata.obs.loc[sample_id, 'cancer_type']
                field = [sample_id, gene, layer_name, cancer_type]
                data.append(field)
    events = pd.DataFrame(data, columns=['sample_id', 'gene', 'event', 'cancer_type'])
    return events

@click.command()
@click.option('-i', '--input', help='input consensus h5ad with obs', required=True)
@click.option('-d', '--driver_genes', help='input CNV driver genes', required=True)
@click.option('-r', '--refflat', help='input refflat txt path', required=True)
@click.option('-og', '--output', help="output h5ad path", required=True)
@click.option('-oe', '--output_events', help="output events table path", required=True)
def make_genelevel_anndata(input, driver_genes, refflat, output, output_events):
    adata = ad.read_h5ad(input)
    layer_name = 'total_cn'
    for chrom in adata.var['chr'].unique():
        chrom_mask = (adata.var['chr']==chrom)
        adata.layers[layer_name] = interp_nan_per_chrom(adata.layers[layer_name], chrom_mask)
    cnv_genes = yaml.safe_load(open(driver_genes).read())
    refflat_cols = ['gene_symbol', 'gene_id', 'chrom', 'strand', 'cdna_start', 'cdna_end', 
                    'orf_start', 'orf_end', 'exons', 'exon_starts', 'exon_ends']
    refflat = pd.read_table(refflat, names=refflat_cols)
    gene_var = get_gene_annotation_from_refflat(cnv_genes, refflat)
    gdata = make_curated_gene_anndata(gene_var, adata)
    cn_events = get_cnv_events_tcga(
        gdata, min_gene_cn_all=3.0, call_cn_amp_cn=7.0,
        homdel_cutoff=0.3,
        a_cn=2, b_cn=0., 
        min_ploidy_cutoff=1.5
    )
    join_cnv_events_tables_to_obs(gdata, cn_events)
    remove_events_not_curated(gdata, cnv_genes)
    events = make_melted_events_from_layers(gdata)
    gdata.write_h5ad(output)
    events.to_csv(output_events, index=False)

if __name__ == "__main__":
    make_genelevel_anndata()