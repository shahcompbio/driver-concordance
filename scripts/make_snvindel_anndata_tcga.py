import click

import tqdm
import pandas as pd
import anndata as ad

class MAF_meta:
    variant_classification_map = {
        "Frame_Shift_Del": "T",
        "Frame_Shift_Ins": "T",
        "In_Frame_Del": "D",
        "In_Frame_Ins": "I",
        "Missense_Mutation": "M",
        "Nonsense_Mutation": "T",
        "Splice_Site": "T",
        "Translation_Start_Site": "M",
        "Nonstop_Mutation": "M",
    }
    layer_map = {'T':'truncating', 'I':'insertion', 'D':'deletion', 'M':'missense'}

def make_anndata(samples, genes, metadata, maf_config):
    """
    Make anndata from samples, genes, and maf paths metadata
    """
    gene_set = set(genes)
    layers_data = {
        t: pd.DataFrame(0, index=samples, columns=genes)
        for t in maf_config.layer_map
    }
    obs = pd.DataFrame(index=samples, columns=['case', 'cancer_type'])
    var = pd.DataFrame(index=genes)
    for sample_id, df in tqdm.tqdm(metadata.groupby('sample_id')):
        cases = df['case'].unique()
        assert len(cases) == 1, cases
        case = cases[0]
        cancer_types = df['cancer_type'].unique()
        assert len(cancer_types) == 1, cancer_types
        cancer_type = cancer_types[0]
        paths = df['path']
        gene2effect = {}
        obs.loc[sample_id, 'case'] = case
        obs.loc[sample_id, 'cancer_type'] = cancer_type
        for path in paths:
            maf = pd.read_table(path, comment='#')
            maf = maf[maf['Variant_Classification'].isin(maf_config.variant_classification_map)]
            maf['effect'] = maf['Variant_Classification'].map(maf_config.variant_classification_map)
            _gene2effect = dict(zip(maf['Hugo_Symbol'], maf['effect']))
            gene2effect.update(_gene2effect)
        for gene, effect in gene2effect.items(): # effect in {'T', 'I', 'D', 'M'}
            if gene not in gene_set: 
                continue
            layers_data[effect].loc[sample_id, gene] = 1
    adata = ad.AnnData(obs=obs, var=var, layers=layers_data)
    return adata

def make_melted_events_from_layers(gdata):
    """
    Generate a melted DataFrame of events based on non-zero values in specified layers of input data.

    This function processes specific layers ('T', 'D', 'I', 'M') of a given AnnData object (`gdata`),
    creating a DataFrame of samples that contain events (non-zero values) for each gene in each layer.
    For each detected event, it records the sample ID, gene name, event type (layer name), and cancer type.

    Parameters
    ----------
    gdata : AnnData
        An AnnData object with multiple layers and metadata, where:
        - `gdata.layers` contains data matrices indexed by layers (e.g., 'T', 'D', 'I', 'M').
        - `gdata.obs.index` provides sample identifiers.
        - `gdata.var.index` provides gene names.
        - `gdata.obs` includes metadata, specifically a 'cancer_type' column.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns:
        - 'sample_id': identifier of the sample containing an event.
        - 'gene': the gene in which an event was detected.
        - 'event': the layer name representing the type of event.
        - 'cancer_type': the cancer type associated with the sample.
    """
    data = []
    layer_names = ['T', 'D', 'I', 'M']
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
@click.option('--maf_path_table', type=str, help="input table of TCGA maf files")
@click.option('--drivers_mutation_path', type=str, help="input driver genes table for simple somatic mutations")
@click.option('--mutation_anndata', type=str, help="output annotated anndata for mutations")
@click.option('--mutation_events', type=str, help="output event table for mutations")
def make_snvindel_anndata(maf_path_table, drivers_mutation_path, mutation_anndata, mutation_events):
    meta = pd.read_csv(maf_path_table)
    drivers_mutation_df = pd.read_table(drivers_mutation_path, comment='#')
    driver_genes = drivers_mutation_df['Gene'].unique().tolist()
    maf_config = MAF_meta()
    sample_ids = meta['sample_id'].unique().tolist()
    gene_names = driver_genes[:]
    adata = make_anndata(sample_ids, gene_names, meta, maf_config)
    adata.write_h5ad(mutation_anndata)
    events = make_melted_events_from_layers(adata)
    events.to_csv(mutation_events, index=False)

if __name__ == "__main__":
    make_snvindel_anndata()